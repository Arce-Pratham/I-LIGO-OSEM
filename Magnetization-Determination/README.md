# Magnet Magnetization Fit (Amperian Loop Model)

Estimates the axial magnetization $M$ of two identical, coaxial cylindrical
permanent magnets from lab-measured axial repulsive force vs. separation
data, using the exact **Amperian-loop** (equivalent bulk current) model —
no small-separation or point-dipole approximation.

The force law is $F(s) = M^2 \, I(s)$, where $I(s)$ is a purely geometric
factor (elliptic-integral double integral, no dependence on $M$). Rather
than inverting this pointwise and averaging the results, $M$ is estimated
by a single **weighted least-squares regression** across all data points —
see [`mag_theory.pdf`](./mag_theory.pdf) for the full physics derivation
and the statistical justification for why this is the correct approach
(and why the naive point-by-point method is biased).

## Contents

| File | Description |
|---|---|
| `fit_magnetization.py` | Main pipeline: physics model + regression fit + diagnostic plot |
| `params.json` | Magnet geometry (diameter, length) |
| `obs.csv` | Lab observations (separation, force[, force uncertainty]) |
| `mag_theory.pdf` | Full physics + statistics derivation |
| `requirements.md` | Dependency list and install instructions |

## Physics & statistics summary

Each magnet (radius $R$, length $L$) is modeled as a stack of circular
Amperian current loops. The mutual axial force between two coaxial loops
of equal radius, separated by $z_\text{sep}$, has an exact closed form in
terms of the complete elliptic integrals $K,E$:

$$
F(s) = M^2 \int_0^{L}\!\!\int_{L+s}^{2L+s} f(z_2-z_1,\,R)\; dz_2\,dz_1
\;\equiv\; M^2\, I(s)
$$

Since $F$ is **linear** in $\theta \equiv M^2$ with known "$x$-value"
$I(s)$, $\theta$ (and hence $M=\sqrt\theta$) is estimated by weighted
least squares across *all* observations at once — using every data point
(including any with $F\le 0$) and taking a single square root at the very
end. This avoids two biases present in a naive per-point
"invert-then-average" approach: a downward bias from Jensen's inequality,
and an upward bias from silently discarding negative force readings. Full
derivation and worked examples: [`mag_theory.pdf`](./mag_theory.pdf).

## Installation

See [`requirements_magnetization.md`](./requirements_magnetization.md). In short:

```bash
pip install numpy pandas scipy plotly
```

## Usage

Run with the default file names (`params.json`, `obs.csv` in the current directory):

```bash
python fit_magnetization.py
```

Point at specific files:

```bash
python fit_magnetization.py --params my_params.json --obs my_obs.csv
```

### Command-line options

| Flag | Effect | Default |
|---|---|---|
| `--params` | Path to geometry JSON | `params.json` |
| `--obs` | Path to observations CSV | `obs.csv` |
| `--txt-out` | Output results text file | `amperian_fit_results.txt` |
| `--html-out` | Output interactive plot file | `amperian_fit_plot.html` |
| `--show` | Open the plot in a browser after running | off |

### Input files

**`params.json`** — magnet geometry, SI units (meters):

```json
{
  "D": 0.01,
  "L": 0.01
}
```

| Field | Meaning |
|---|---|
| `D` | Magnet diameter (radius $R = D/2$ is used internally) |
| `L` | Magnet length |

**`obs.csv`** — one row per measurement:

| Column | Meaning | Required? |
|---|---|---|
| `z` | Face-to-face separation $s$ (meters) | yes |
| `F` | Measured axial force (Newtons) | yes |
| `sigma_F` | Per-point force uncertainty (Newtons) | optional |

If `sigma_F` is present, it's used as a proper inverse-variance regression
weight ($w_i = 1/\sigma_i^2$) and the standard error on $M$ is computed
directly from it. If absent, all points are weighted equally and the noise
level is instead estimated from the fit's own residual scatter — both
cases are handled automatically; see `mag_theory.pdf` §5.2 for why the
formula differs between them.

### Outputs

- **`amperian_fit_results.txt`** — geometry, every observation alongside
  its computed $I(s)$, fitted force, and residual, followed by a summary
  block: $\hat\theta=\hat M^2$, $\mathrm{SE}(\hat\theta)$, $\hat M$,
  $\mathrm{SE}(\hat M)$, residual scale (if applicable), and $R^2$.
- **`amperian_fit_plot.html`** — interactive 3-panel Plotly diagnostic:
  1. Force vs. separation, with the fitted curve $\hat M^2 I(s)$ overlaid
  2. The actual linear regression view — force vs. $I(s)$, with the
     through-origin fit line (this is the plot the estimate is computed
     from)
  3. Residuals vs. separation, to visually check the fit quality (a
     systematic pattern here would suggest a model mismatch — e.g.
     non-identical magnets, misalignment, or an unmodeled force
     contribution)

  The plot title states $\hat M \pm \mathrm{SE}(\hat M)$ and $R^2$
  directly.

## A note on the statistics

$R^2$ here is measured **about a zero baseline** (predicting $F=0$
everywhere), not about the mean of $F$, since this model has no intercept
by construction — the force law is exactly zero at infinite separation,
so "no effect" is the natural baseline for a through-origin fit like this
one, rather than the sample mean used in an ordinary (intercept-fitting)
$R^2$.
