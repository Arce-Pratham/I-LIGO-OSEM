"""
osem_forces.py
===============

Axial magnetic force between a solenoid coil and a cylindrical permanent
magnet, as used in a LIGO-style OSEM (Optical Sensor and ElectroMagnetic)
actuator.
"""

import argparse
import json
import os
import time
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import plotly.graph_objects as go
from numpy.polynomial.legendre import leggauss

# Permeability of free space (H/m)
MU0 = 4e-7 * np.pi


#
# Configuration data structures
#

##Solenoid coil geometry and drive current (SI units, meters/amps).
@dataclass
class Coil:
    L: float                       # axial length
    R: float                       # inner radius
    Ra: float                      # outer radius
    turns_total: int               # total number of turns
    current_A: float               # drive current per turn
    sigma_override: Optional[float] = None  # if set, bypass N*I/volume calc


##Cylindrical permanent magnet geometry and magnetization (SI units).
@dataclass
class Magnet:
    a: float                       # radius
    h: float                       # height
    Mz: float                      # axial magnetization (A/m)


##Gauss-Legendre quadrature orders and the z-scan range for the sweep.
@dataclass
class IntegrationCfg:
    n_theta: int                   # azimuthal angle nodes
    n_r1: int                      # coil radial nodes
    n_z1: int                      # coil axial nodes
    n_r2: int                      # magnet radial nodes
    n_z2: int                      # magnet axial nodes
    z_scan_min: float              # sweep start, magnet center position (m)
    z_scan_max: float              # sweep end, magnet center position (m)
    z_scan_num: int                # number of sweep points


##Output plot settings. save_path should end in .html (interactive plotly plot).
@dataclass
class PlotCfg:
    save_path: str


#
# JSON config I/O
#

###Reasonable default OSEM.json contents, used to bootstrap a fresh config.
def default_config_dict() -> dict:
    return {
        "coil": {
            "length_mm": 10.0,
            "r_inner_mm": 4.0,
            "r_outer_mm": 6.0,
            "turns_total": 750,
            "current_A": 0.02,
            "sigma_override_A_per_m2": None
        },
        "magnet": {
            "radius_mm": 3.175,
            "height_mm": 3.175,
            "Mz_A_per_m": 8.6e5
        },
        "integration": {
            "n_theta": 24,
            "n_r1": 8,
            "n_z1": 8,
            "n_r2": 8,
            "n_z2": 8,
            "z_scan_mm": {
                "min": -10.0,
                "max": 10.0,
                "num": 41
            }
        },
        "plot": {
            "save_path": "Fz_vs_z.html"
        }
    }


##Write a default config file at `path` if one does not already exist.
def ensure_json(path: str) -> None:
    if not os.path.exists(path):
        cfg = default_config_dict()
        with open(path, 'w') as f:
            json.dump(cfg, f, indent=2)
        print(f"[info] '{path}' not found. Wrote default config.")


##Load OSEM.json and convert all lengths from mm to meters.
def load_config(path: str) -> Tuple[Coil, Magnet, IntegrationCfg, PlotCfg]:
    with open(path, 'r') as f:
        cfg = json.load(f)

    coil = Coil(
        L=cfg['coil']['length_mm'] * 1e-3,
        R=cfg['coil']['r_inner_mm'] * 1e-3,
        Ra=cfg['coil']['r_outer_mm'] * 1e-3,
        turns_total=cfg['coil']['turns_total'],
        current_A=cfg['coil']['current_A'],
        sigma_override=cfg['coil']['sigma_override_A_per_m2'],
    )

    magnet = Magnet(
        a=cfg['magnet']['radius_mm'] * 1e-3,
        h=cfg['magnet']['height_mm'] * 1e-3,
        Mz=cfg['magnet']['Mz_A_per_m'],
    )

    integ = IntegrationCfg(
        n_theta=cfg['integration']['n_theta'],
        n_r1=cfg['integration']['n_r1'],
        n_z1=cfg['integration']['n_z1'],
        n_r2=cfg['integration']['n_r2'],
        n_z2=cfg['integration']['n_z2'],
        z_scan_min=cfg['integration']['z_scan_mm']['min'] * 1e-3,
        z_scan_max=cfg['integration']['z_scan_mm']['max'] * 1e-3,
        z_scan_num=cfg['integration']['z_scan_mm']['num'],
    )

    plotcfg = PlotCfg(save_path=cfg['plot']['save_path'])
    return coil, magnet, integ, plotcfg


#
# Physics & numerics
#

def compute_sigma(coil: Coil) -> float:
    """
    Effective volumetric current density (A/m^2) of the coil.

    Equal to (total ampere-turns) / (coil cross-sectional volume), unless
    `sigma_override` is set in the config, in which case that value is used
    directly.
    """
    if coil.sigma_override is not None:
        return coil.sigma_override
    volume = np.pi * (coil.Ra**2 - coil.R**2) * coil.L
    return (coil.turns_total * coil.current_A) / volume


##Gauss-Legendre quadrature nodes/weights mapped from [-1, 1] to [a, b].
def gl_nodes_weights(n: int, a: float, b: float) -> Tuple[np.ndarray, np.ndarray]:
    xi, wi = leggauss(n)
    x = 0.5 * (b - a) * xi + 0.5 * (a + b)
    w = 0.5 * (b - a) * wi
    return x, w


def force_at_z(zc: float, coil: Coil, magnet: Magnet,
                integ: IntegrationCfg, sigma: float) -> float:
    """
    Axial force (N) on the magnet when its center is at axial position `zc`,
    via the 5D Gauss-Legendre quadrature:

        Fz(zc) = -3*mu0*sigma*Mz*pi *
                 ∫∫∫∫∫  [r1^2 * r2 * (r1 - r2*cos(theta)) * (z2 - z1)] / D^5  d(theta) dr1 dz1 dr2 dz2

    where:
        D^2 = r1^2 + r2^2 - 2*r1*r2*cos(theta) + (z2 - z1)^2

        (r1, z1) : coil source point, r1 in [R, Ra], z1 in [-L/2, L/2]
        (r2, z2) : magnet field point, r2 in [0, a],  z2 in [zc - h/2, zc + h/2]
        theta    : azimuthal separation between source and field filaments,
                   integrated over [0, 2*pi]

    The nested-loop structure below mirrors the mathematical integral term
    for term (innermost: theta: then r1, z1; then r2, z2) and is left
    algebraically identical to the original implementation.
    """
    # Magnet integrates over its own height, centered on zc
    z2_min, z2_max = zc - magnet.h / 2.0, zc + magnet.h / 2.0

    theta_nodes, theta_w = gl_nodes_weights(integ.n_theta, 0.0, 2.0 * np.pi)
    r1_nodes, r1_w = gl_nodes_weights(integ.n_r1, coil.R, coil.Ra)
    z1_nodes, z1_w = gl_nodes_weights(integ.n_z1, -coil.L / 2.0, coil.L / 2.0)
    r2_nodes, r2_w = gl_nodes_weights(integ.n_r2, 0.0, magnet.a)
    z2_nodes, z2_w = gl_nodes_weights(integ.n_z2, z2_min, z2_max)

    cos_theta = np.cos(theta_nodes)

    prefactor = -3.0 * MU0 * sigma * magnet.Mz * np.pi
    Fz = 0.0

    # Outer two sums: magnet cross-section (field point)
    for z2, wz2 in zip(z2_nodes, z2_w):
        for r2, wr2 in zip(r2_nodes, r2_w):
            magnet_point_accum = 0.0

            # Inner two sums: coil cross-section (source point)
            for z1, wz1 in zip(z1_nodes, z1_w):
                dz = z2 - z1
                for r1, wr1 in zip(r1_nodes, r1_w):
                    # Squared distance between the two filaments, integrated
                    # over the azimuthal angle theta between them
                    D2 = r1*r1 + r2*r2 - 2.0*r1*r2*cos_theta + dz*dz
                    D5 = D2**2.5

                    numerator = r1*r1 * r2 * (r1 - r2 * cos_theta) * dz
                    integrand_theta = numerator / D5

                    theta_integral = np.sum(theta_w * integrand_theta)
                    magnet_point_accum += wz1 * wr1 * theta_integral

            Fz += wz2 * wr2 * magnet_point_accum

    return prefactor * Fz


#
# Simulation pipeline
#

##Apply optional command-line overrides on top of the loaded JSON config.
def apply_cli_overrides(args: argparse.Namespace, integ: IntegrationCfg) -> None:
    if args.orders is not None:
        integ.n_theta = integ.n_r1 = integ.n_z1 = integ.n_r2 = integ.n_z2 = int(args.orders)
    if args.zmin is not None:
        integ.z_scan_min = args.zmin * 1e-3
    if args.zmax is not None:
        integ.z_scan_max = args.zmax * 1e-3
    if args.znum is not None:
        integ.z_scan_num = int(args.znum)


## Sweep the magnet center position over the configured z range and
## evaluate the axial force at each point.
def run_z_scan(coil: Coil, magnet: Magnet, integ: IntegrationCfg,
                sigma: float, log) -> Tuple[np.ndarray, np.ndarray]:
    z_grid = np.linspace(integ.z_scan_min, integ.z_scan_max, integ.z_scan_num)
    F = np.zeros_like(z_grid)

    # Write a zeroed placeholder up front, matching the original behavior
    save_force_csv(z_grid, F)

    log("Computing F(z)...")
    progress_step = max(1, integ.z_scan_num // 20)
    for i, zc in enumerate(z_grid):
        F[i] = force_at_z(zc, coil, magnet, integ, sigma)
        if (i + 1) % progress_step == 0:
            log(f"  {i + 1}/{integ.z_scan_num} done")

    return z_grid, F


## Save the z (mm) vs. force (N) sweep results to CSV.
def save_force_csv(z_grid: np.ndarray, F: np.ndarray,
                    path: str = "Fz_vs_z.csv") -> None:
    np.savetxt(
        path,
        np.column_stack([z_grid * 1e3, F]),
        delimiter=",",
        header="z_mm,F_N",
        comments="",
    )


""" Build an interactive HTML plot of force vs. magnet position (via plotly)
    and save it to `save_path`. Each computed point is individually
    hoverable (shows exact z and Fz), and the force-maximizing ('sweet
    spot') position is marked with a vertical line and a highlighted marker.
"""
def save_force_plot(z_grid: np.ndarray, F: np.ndarray, z_max_mm: float,
                     save_path: str, show: bool) -> None:

    z_mm = z_grid * 1e3
    idx_max = int(np.argmax(np.abs(F)))

    fig = go.Figure()

    # Main curve, with every quadrature/scan point individually hoverable
    fig.add_trace(go.Scatter(
        x=z_mm,
        y=F,
        mode='lines+markers',
        name='Fz(z)',
        line=dict(width=2),
        marker=dict(size=6),
        hovertemplate='z = %{x:.3f} mm<br>Fz = %{y:.6e} N<extra></extra>',
    ))

    # Highlight the sweet spot as its own trace so it stands out on hover
    fig.add_trace(go.Scatter(
        x=[z_mm[idx_max]],
        y=[F[idx_max]],
        mode='markers',
        name=f'Sweet spot (z={z_max_mm:.2f} mm)',
        marker=dict(size=12, symbol='star', color='black'),
        hovertemplate='Sweet spot<br>z = %{x:.3f} mm<br>Fz = %{y:.6e} N<extra></extra>',
    ))

    fig.add_vline(
        x=z_max_mm, line_dash='dash', line_color='black', opacity=0.4,
        annotation_text=f'z_max={z_max_mm:.2f} mm', annotation_position='top',
    )

    fig.update_layout(
        title='Axial force vs z (full 5D integral)',
        xaxis_title='Magnet center z (mm)',
        yaxis_title='Axial force Fz (N)',
        template='plotly_white',
        hovermode='closest',
    )

    fig.write_html(save_path, include_plotlyjs='cdn')
    if show:
        fig.show()


#
# CLI / entry point
#

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="OSEM axial force via 5D integral")
    parser.add_argument('--json', default='OSEM.json',
                         help='Path to JSON config (auto-created if missing)')
    parser.add_argument('--orders', type=int, default=None,
                         help='Set same Gauss-Legendre order for all 5 dims')
    parser.add_argument('--zmin', type=float, default=None,
                         help='Override z-scan min (mm)')
    parser.add_argument('--zmax', type=float, default=None,
                         help='Override z-scan max (mm)')
    parser.add_argument('--znum', type=int, default=None,
                         help='Override z-scan num points')
    parser.add_argument('--show', action='store_true',
                         help='Show plot interactively')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ensure_json(args.json)
    coil, magnet, integ, plotcfg = load_config(args.json)
    apply_cli_overrides(args, integ)

    log_lines = []

    def log(msg: str) -> None:
        print(msg)
        log_lines.append(msg)

    log("=== OSEM axial force (fresh) ===")
    log(f"json: {args.json}")
    log(f"orders: theta={integ.n_theta}, r1={integ.n_r1}, z1={integ.n_z1}, "
        f"r2={integ.n_r2}, z2={integ.n_z2}")

    sigma = compute_sigma(coil)
    log(f"sigma (A/m^2): {sigma:.6e}")

    t0 = time.time()
    z_grid, F = run_z_scan(coil, magnet, integ, sigma, log)
    t1 = time.time()

    save_force_csv(z_grid, F)

    idx_max = np.argmax(np.abs(F))
    z_max_mm = z_grid[idx_max] * 1e3
    log(f"Sweet spot: z = {z_max_mm:.3f} mm")
    log(f"Runtime: {t1 - t0:.2f} s")
    log("Saved plot data")

    save_force_plot(z_grid, F, z_max_mm, plotcfg.save_path, show=args.show)
    log(f"Saved plot to {plotcfg.save_path}")

    with open('osem_log.txt', 'w') as f:
        f.write("\n".join(log_lines))


if __name__ == '__main__':
    main()
