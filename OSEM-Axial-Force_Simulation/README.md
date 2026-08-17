<!-- GitHub Actions Workflow Status Badge -->
[![Run OSEM Forces Simulation Tests](https://github.com/Arce-Pratham/I-LIGO-OSEM/.github/workflows/tests.yml/badge.svg)](https://github.com/Arce-Pratham/I-LIGO-OSEM/.github/workflows/tests.yml)

## Overview
This repository contains the simulation code for analyzing axial forces. 

### Automated Testing
The status badge above reflects the live operational state of the automated test suite running against `OSEM-Axial-Force_simulation/osem_forces.py`. 
- **Passing (Green):** The script executes successfully without syntax errors or runtime exceptions.
- **Failing (Red):** An error was introduced or a dependency is missing.
- **Running (Yellow):** GitHub Actions is currently executing the simulation pipeline.


# OSEM Axial Force Simulation

Numerical simulation of the axial magnetic force between the coil and
permanent magnet in a LIGO-style **OSEM** (Optical Sensor and
ElectroMagnetic) actuator.

The coil is modeled as a bulk azimuthal current distribution and the
magnet as a continuum of axial magnetic dipoles. The mutual force is
computed as a full **five-dimensional integral** via Gauss-Legendre
quadrature — no elliptic-integral or point-dipole approximations are used.
See [`osem_theory.pdf`](./osem_theory.pdf) for the complete derivation.

## Contents

| File | Description |
|---|---|
| `osem_forces.py` | Main simulation script |
| `OSEM.json` | Geometry, material, and quadrature configuration (auto-generated with defaults on first run if missing) |
| `osem_theory.pdf`| Physics derivation of the force integral the code evaluates |
| `requirements.md` | Dependency list and install instructions |

## Physics summary

For a coil of inner/outer radius $R,R_a$, length $L$, carrying bulk current
density $\sigma$, and a coaxial magnet of radius $a$, height $h$,
magnetization $M_z$, centered at axial position $z_c$, the axial force is

$$
F_z(z_c) = -3\pi\mu_0\,\sigma\, M_z\int_{z2}
\int_{r2}\int_{z1}\int_{r1}\int_{\theta}
\frac{r_1^2\, r_2\,(r_1 - r_2\cos\theta)\,(z_2-z_1)}{D^5}
\,d\theta\,dr_1\,dz_1\,dr_2\,dz_2
$$

integrated over the coil cross-section ($r_1,z_1$), the magnet
cross-section ($r_2,z_2$), and the relative azimuthal angle $\theta$
between a coil filament and a magnet filament, with
$D^2 = r_1^2+r_2^2-2r_1r_2\cos\theta+(z_2-z_1)^2$. Full derivation,
figures, and the numerical evaluation scheme are in `osem_theory.pdf`.

## Installation

See [`requirements.md`](./requirements.md). In short:

```bash
pip install numpy plotly
```

## Usage

Run with defaults (auto-creates `OSEM.json` on first run):

```bash
python osem_forces.py
```

Point at a specific config file:

```bash
python osem_forces.py --json my_config.json
```

### Command-line overrides

These override whatever is in the JSON config for a single run, without
editing the file:

| Flag | Effect |
|---|---|
| `--orders N` | Set all five Gauss-Legendre quadrature orders to `N` |
| `--zmin MM` | Override z-scan start (mm) |
| `--zmax MM` | Override z-scan end (mm) |
| `--znum N` | Override number of z-scan points |
| `--show` | Open the interactive plot in a browser after running |

Example — quick, low-resolution scan for testing:

```bash
python osem_forces.py --orders 6 --znum 15 --show
```

### Outputs

Each run produces:

- **`Fz_vs_z.csv`** — raw sweep data (`z_mm, F_N`)
- **`Fz_vs_z.html`** (or whatever `plot.save_path` is set to) — interactive
  Plotly plot; hover over any point to see its exact `z` and `Fz`, with the
  force-maximizing ("sweet spot") position marked
- **`osem_log.txt`** — a copy of everything printed to the console
  (quadrature orders, computed `sigma`, sweet-spot position, runtime)

## Configuration (`OSEM.json`)

All lengths are specified in **millimeters** in the file and converted to
meters internally; other quantities are SI.

```json
{
  "coil": {
    "length_mm": 10.0,
    "r_inner_mm": 4.0,
    "r_outer_mm": 6.0,
    "turns_total": 750,
    "current_A": 0.02,
    "sigma_override_A_per_m2": null
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
    "z_scan_mm": { "min": -10.0, "max": 10.0, "num": 41 }
  },
  "plot": {
    "save_path": "Fz_vs_z.html"
  }
}
```

| Field | Meaning |
|---|---|
| `coil.length_mm`, `r_inner_mm`, `r_outer_mm` | Coil geometry |
| `coil.turns_total`, `current_A` | Used to compute the bulk current density $\sigma = NI/V_\text{coil}$ |
| `coil.sigma_override_A_per_m2` | If set (not `null`), used directly as $\sigma$ instead of computing it from turns/current — useful for matching a known/measured value |
| `magnet.radius_mm`, `height_mm` | Magnet geometry |
| `magnet.Mz_A_per_m` | Magnet's axial magnetization |
| `integration.n_theta`, `n_r1`, `n_z1`, `n_r2`, `n_z2` | Gauss-Legendre quadrature order for each of the five integration variables — higher = more accurate, slower |
| `integration.z_scan_mm` | Range and point count of the magnet-center sweep |
| `plot.save_path` | Output path for the interactive HTML plot |

Increasing the quadrature orders or the z-scan point count trades runtime
for accuracy/resolution; `--orders`, `--znum` etc. let you experiment
without editing the file.

## A note on accuracy

This simulation performs the **exact extended-body force integral** — it
does not assume the magnet is small compared to the length scale over
which the coil's field varies. Point-dipole
approximation is only accurate in that small-magnet limit, which a compact OSEM geometry does not always satisfy. If
a point-dipole-based reference shows a
discrepancy, this is the likely reason — see `osem_theory.pdf` for details.

## This repository has been created and published by Pratham Patil
