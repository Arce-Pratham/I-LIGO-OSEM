"""
fit_magnetization.py
=====================

Estimate the axial magnetization M of two identical, coaxial cylindrical
magnets from lab-measured axial repulsive force vs. separation data, using
the Amperian-loop (equivalent bulk current) model.

Physics (unchanged from the original derivation)
--------------------------------------------------
Each magnet of radius R, length L is modeled as a stack of circular
Amperian current loops. The axial force between the two magnets, separated
face-to-face by distance s, is

    F(s) = M^2 * I(s)

where I(s) is a purely geometric factor (no dependence on M) obtained by
double-integrating the exact elliptic-integral mutual-force law for two
coaxial loops over both magnets' lengths. See `theory/mag_theory.pdf` for
the full derivation.

Statistics (the actual change from the original script)
----------------------------------------------------------
F(s) = M^2 * I(s) is linear in theta = M^2, with I(s) as the known
"x-value" and F as the noisy "y-value". Rather than inverting each point
individually (M_i = sqrt(F_i / I(s_i))) and averaging the results -- which
is provably biased, both by Jensen's inequality and by silently discarding
any point with F_i <= 0 -- this script fits theta by (weighted) least
squares across *all* points simultaneously, and takes a single square root
at the end:

    theta_hat = sum(w_i * I(s_i) * F_i) / sum(w_i * I(s_i)^2)
    M_hat     = sqrt(theta_hat)

Uncertainty on theta_hat comes from ordinary least-squares theory (using
per-point weights if available, otherwise the residual scatter), and is
propagated to M_hat via the delta method. See `theory/mag_theory.pdf`,
Part II, for the full statistical justification.

Usage
-----
    python fit_magnetization.py --params params.json --obs obs.csv
"""

import argparse
import json
import sys
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scipy.constants as const
from scipy.special import ellipk, ellipe
from scipy.integrate import dblquad
import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ============================================================
# Configuration
# ============================================================

@dataclass
class MagnetGeometry:
    """Identical-magnet geometry, SI units (meters)."""
    R: float   # radius
    L: float   # length


def load_geometry(path: str) -> MagnetGeometry:
    """Load magnet radius/length from params.json (D = diameter, L = length, meters)."""
    with open(path, 'r') as f:
        params = json.load(f)
    try:
        L = params['L']
        D = params['D']
    except KeyError as e:
        raise KeyError(f"'{path}' is missing required key {e}; it must contain 'D' and 'L'.")
    return MagnetGeometry(R=D / 2.0, L=L)


def load_observations(path: str) -> pd.DataFrame:
    """
    Load (z, F) observations, coercing to numeric and dropping invalid rows.
    An optional 'sigma_F' column (per-point force uncertainty) is used as
    regression weights if present.
    """
    df = pd.read_csv(path)
    if 'z' not in df.columns or 'F' not in df.columns:
        raise ValueError(f"'{path}' must contain 'z' and 'F' columns.")

    df['z'] = pd.to_numeric(df['z'], errors='coerce')
    df['F'] = pd.to_numeric(df['F'], errors='coerce')
    if 'sigma_F' in df.columns:
        df['sigma_F'] = pd.to_numeric(df['sigma_F'], errors='coerce')

    n_before = len(df)
    df = df.dropna(subset=['z', 'F']).reset_index(drop=True)
    if len(df) < n_before:
        print(f"[warning] dropped {n_before - len(df)} row(s) with invalid z/F.")

    return df


# ============================================================
# Physics: exact Amperian-loop force model
# ============================================================

def loop_pair_force_kernel(z_sep: float, R: float) -> float:
    """
    Axial force per unit I1*I2 between two coaxial circular current loops
    of equal radius R, separated axially by z_sep (SI units), via the
    exact elliptic-integral formula (no small-separation or point-dipole
    approximation).
    """
    if np.abs(z_sep) < 1e-10:
        return 0.0

    k_sq = (4.0 * R ** 2) / ((2.0 * R) ** 2 + z_sep ** 2)
    Kk = ellipk(k_sq)
    Ek = ellipe(k_sq)

    prefactor = (const.mu_0 * z_sep) / np.sqrt((2.0 * R) ** 2 + z_sep ** 2)
    main_term = ((2.0 * R ** 2 + z_sep ** 2) * Ek / z_sep ** 2) - Kk

    return prefactor * main_term


def geometry_factor(s: float, R: float, L: float) -> float:
    """
    I(s): the purely geometric double integral of the loop-pair kernel over
    both magnets' cross-sections (stacks of loops), independent of M.

        source magnet: z1 in [0, L]
        target magnet: z2 in [L+s, 2L+s]

    F(s) = M^2 * I(s).
    """
    def integrand(z2, z1):
        return loop_pair_force_kernel(z2 - z1, R)

    I_val, _ = dblquad(integrand, 0.0, L, lambda z1: L + s, lambda z1: 2.0 * L + s)
    return I_val


# ============================================================
# Statistics: weighted regression estimator for M
# ============================================================

@dataclass
class FitResult:
    theta_hat: float       # M^2
    se_theta: float
    M_hat: float
    se_M: float
    sigma_hat: float       # residual scale (only meaningful when weights are uniform)
    r_squared: float        # about a zero baseline (model has no intercept)
    weighted: bool


def fit_magnetization(I_vals: np.ndarray, F_vals: np.ndarray,
                       sigma_vals: Optional[np.ndarray] = None) -> FitResult:
    """
    Weighted least-squares fit of F_i = theta * I_i (theta = M^2), then a
    single square root to recover M. See mag_theory.pdf Part II for the
    full derivation of every formula used here.
    """
    n = len(I_vals)
    weighted = sigma_vals is not None

    if weighted:
        w = 1.0 / sigma_vals ** 2
    else:
        w = np.ones(n)

    sum_wI2 = np.sum(w * I_vals ** 2)
    theta_hat = np.sum(w * I_vals * F_vals) / sum_wI2

    residuals = F_vals - theta_hat * I_vals

    if weighted:
        # Per-point noise already known -> standard error follows directly.
        se_theta = np.sqrt(1.0 / sum_wI2)
        sigma_hat = np.nan  # not estimated; per-point sigmas were supplied instead
    else:
        # Noise level unknown -> estimate it from how far the data falls
        # from the fitted line (the residuals), then propagate.
        sigma_hat = np.sqrt(np.sum(residuals ** 2) / (n - 1))
        se_theta = sigma_hat / np.sqrt(np.sum(I_vals ** 2))

    M_hat = np.sqrt(theta_hat)
    se_M = se_theta / (2.0 * M_hat)  # delta method, d(sqrt theta)/d theta = 1/(2 sqrt theta)

    # R^2 measured about a zero baseline, since this model has no intercept
    # (predicting F=0 everywhere is the natural "no effect" baseline here).
    ss_res = np.sum(w * residuals ** 2)
    ss_tot = np.sum(w * F_vals ** 2)
    r_squared = 1.0 - ss_res / ss_tot

    return FitResult(theta_hat=theta_hat, se_theta=se_theta, M_hat=M_hat,
                      se_M=se_M, sigma_hat=sigma_hat, r_squared=r_squared,
                      weighted=weighted)


# ============================================================
# Plotting
# ============================================================

def build_diagnostic_plot(df: pd.DataFrame, geom: MagnetGeometry, fit: FitResult) -> go.Figure:
    """
    Three-panel interactive diagnostic:
      1. F vs separation s, with the fitted curve F_hat(s) = M_hat^2 * I(s)
      2. The linearized regression view: F vs I(s), with the fitted line
         through the origin -- this is the plot the fit is actually done on
      3. Residuals vs s, to visually check the fit quality
    """
    s_vals = df['z'].to_numpy()
    F_vals = df['F'].to_numpy()
    I_vals = df['I_s'].to_numpy()
    resid = df['residual'].to_numpy()

    s_smooth = np.linspace(s_vals.min(), s_vals.max(), 60)
    I_smooth = np.array([geometry_factor(s, geom.R, geom.L) for s in s_smooth])
    F_smooth = fit.theta_hat * I_smooth

    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=(
            "Force vs. separation, Using fitted model",
            "Regression view: force vs. geometry factor I(s)",
            "Residuals: F_obs \u2212 fitted F",
        ),
        vertical_spacing=0.09,
    )

    # --- Panel 1: F vs s ---
    fig.add_trace(go.Scatter(
        x=s_vals, y=F_vals, mode='markers', name='observed F',
        marker=dict(color='#BC8F8F', size=8),
        hovertemplate='s=%{x:.4g} m<br>F=%{y:.4g} N<extra></extra>',
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=s_smooth, y=F_smooth, mode='lines', name='fitted F(s)=M\u0302\u00b2 I(s)',
        line=dict(color='black', width=2.2),
        hovertemplate='s=%{x:.4g} m<br>fitted F=%{y:.4g} N<extra></extra>',
    ), row=1, col=1)

    # --- Panel 2: F vs I(s), the linear regression ---
    fig.add_trace(go.Scatter(
        x=I_vals, y=F_vals, mode='markers', name='(I(s_i), F_i)',
        marker=dict(color='#E2725B', size=8), showlegend=False,
        hovertemplate='I(s)=%{x:.4g}<br>F=%{y:.4g} N<extra></extra>',
    ), row=2, col=1)
    I_line = np.linspace(0, I_vals.max() * 1.05, 2)
    fig.add_trace(go.Scatter(
        x=I_line, y=fit.theta_hat * I_line, mode='lines',
        name=f'fit: slope = M\u0302\u00b2 = {fit.theta_hat:.4e}',
        line=dict(color='black', width=2.2),
    ), row=2, col=1)

    # --- Panel 3: residuals ---
    fig.add_trace(go.Scatter(
        x=s_vals, y=resid, mode='markers', name='residual', showlegend=False,
        marker=dict(color='#E1AD01', size=8, line=dict(color='black', width=0.6)),
        hovertemplate='s=%{x:.4g} m<br>residual=%{y:.4g} N<extra></extra>',
    ), row=3, col=1)
    fig.add_hline(y=0, line=dict(color='black', width=1, dash='dash'), row=3, col=1)

    fig.update_xaxes(title_text='separation s (m)', row=1, col=1)
    fig.update_yaxes(title_text='force F (N)', row=1, col=1)
    fig.update_xaxes(title_text='geometry factor I(s)', row=2, col=1)
    fig.update_yaxes(title_text='force F (N)', row=2, col=1)
    fig.update_xaxes(title_text='separation s (m)', row=3, col=1)
    fig.update_yaxes(title_text='residual (N)', row=3, col=1)

    weight_note = "per-point sigma_F weights" if fit.weighted else "equal weights (\u03c3 estimated from residuals)"
    fig.update_layout(
        height=1000,
        template='plotly_white',
        title=dict(
            text=(f"Amperian-loop fit  \u2014  "
                  f"M\u0302 = {fit.M_hat:.4e} \u00b1 {fit.se_M:.2e} A/m  "
                  f"(R\u00b2={fit.r_squared:.4f}, {weight_note})"),
            x=0.5,
        ),
        hovermode='closest',
    )
    return fig


# ============================================================
# Main pipeline
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description="Fit magnet axial magnetization from force-vs-separation data "
                    "via the exact Amperian-loop model and weighted regression.")
    parser.add_argument('--params', default='params.json', help='Path to geometry JSON (D, L in meters)')
    parser.add_argument('--obs', default='obs.csv', help='Path to observations CSV (columns: z, F[, sigma_F])')
    parser.add_argument('--txt-out', default='amperian_fit_results.txt', help='Output results text file')
    parser.add_argument('--html-out', default='amperian_fit_plot.html', help='Output interactive plot file')
    parser.add_argument('--show', action='store_true', help='Open the plot in a browser after running')
    args = parser.parse_args()

    try:
        geom = load_geometry(args.params)
        df = load_observations(args.obs)
    except (FileNotFoundError, KeyError, ValueError) as e:
        print(f"Error: {e}")
        sys.exit(1)

    print(f"Loaded geometry: R={geom.R:.6g} m, L={geom.L:.6g} m")
    print(f"Loaded {len(df)} observation(s) from {args.obs}")

    print("Computing geometry factor I(s) for each observation "
          "(elliptic-integral double integral)...")
    df['I_s'] = [geometry_factor(s, geom.R, geom.L) for s in df['z']]

    sigma_vals = df['sigma_F'].to_numpy() if 'sigma_F' in df.columns else None
    fit = fit_magnetization(df['I_s'].to_numpy(), df['F'].to_numpy(), sigma_vals)
    df['F_fit'] = fit.theta_hat * df['I_s']
    df['residual'] = df['F'] - df['F_fit']

    print("\n--- Fit summary ---")
    print(f"theta_hat (= M^2)        : {fit.theta_hat:.6e}")
    print(f"SE(theta_hat)            : {fit.se_theta:.6e}")
    print(f"M_hat                    : {fit.M_hat:.6e} A/m")
    print(f"SE(M_hat)                : {fit.se_M:.6e} A/m")
    if not fit.weighted:
        print(f"sigma_hat (residual scale): {fit.sigma_hat:.6e} N")
    print(f"R^2 (about zero baseline) : {fit.r_squared:.6f}")

    # --- Save results text file ---
    with open(args.txt_out, 'w') as f:
        f.write("--- Amperian-loop magnetization fit ---\n")
        f.write(f"geometry: R={geom.R:.6e} m, L={geom.L:.6e} m\n")
        f.write(f"n_observations: {len(df)}\n")
        f.write(f"weighted: {fit.weighted}\n\n")
        f.write(df[['z', 'F', 'I_s', 'F_fit', 'residual']].to_string(index=False))
        f.write("\n\n--- Summary ---\n")
        f.write(f"theta_hat (M^2) = {fit.theta_hat:.6e}\n")
        f.write(f"SE(theta_hat)   = {fit.se_theta:.6e}\n")
        f.write(f"M_hat           = {fit.M_hat:.6e} A/m\n")
        f.write(f"SE(M_hat)       = {fit.se_M:.6e} A/m\n")
        if not fit.weighted:
            f.write(f"sigma_hat       = {fit.sigma_hat:.6e} N\n")
        f.write(f"R^2 (about zero)= {fit.r_squared:.6f}\n")
    print(f"\nSaved results to '{args.txt_out}'")

    # --- Save interactive plot ---
    fig = build_diagnostic_plot(df, geom, fit)
    fig.write_html(args.html_out, include_plotlyjs='cdn')
    print(f"Saved plot to '{args.html_out}'")
    if args.show:
        fig.show()

    print("\nDone.")


if __name__ == "__main__":
    main()
