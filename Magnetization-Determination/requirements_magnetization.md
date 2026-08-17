# Requirements

## Python

Python 3.8 or newer (the code uses `dataclasses`, f-strings, and type hints).

## Dependencies

| Package | Minimum version | Used for |
|---|---|---|
| [`numpy`](https://numpy.org/) | 1.20 | Array math, weighted least-squares arithmetic |
| [`pandas`](https://pandas.pydata.org/) | 1.3 | Reading/validating `obs.csv`, tabular results |
| [`scipy`](https://scipy.org/) | 1.7 | `scipy.constants.mu_0`, complete elliptic integrals (`scipy.special.ellipk`, `ellipe`), 2D numerical integration (`scipy.integrate.dblquad`) |
| [`plotly`](https://plotly.com/python/) | 5.0 | Interactive 3-panel HTML diagnostic plot |

Everything else used (`argparse`, `json`, `sys`, `dataclasses`, `typing`) is part of the Python standard library.

## Install

```bash
pip install numpy pandas scipy plotly
```

Or, saving the block below as `requirements.txt`:

```
numpy>=1.20
pandas>=1.3
scipy>=1.7
plotly>=5.0
```

```bash
pip install -r requirements.txt
```

## Notes

- No GPU or compiled extensions are required.
- Runtime is dominated by `scipy.integrate.dblquad` calls (one per observation, plus ~60 more to draw the smooth fitted curve on the plot). Each call is fast (elliptic integrals are cheap to evaluate), so even a few hundred observations typically finish in well under a minute.
- The interactive plot is written as a self-contained `.html` file that loads the Plotly.js library from a CDN (`include_plotlyjs='cdn'`), keeping the output file small; an internet connection is needed the first time it's opened in a browser (swap to `include_plotlyjs=True` in the code to embed the library for fully offline viewing).
- `scipy.special.ellipk`/`ellipe` take the parameter `m = k²` directly (not the modulus `k`) — this matches the convention already used throughout the code, so no conversion is needed if you modify it.
