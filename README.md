# OSEM Voice Coil Actuator Characterization

Python-based computational tools for the theoretical and experimental characterization of the **OSEM/BOSEM voice coil actuator** used in LIGO suspension systems.

The repository contains two complementary components:

## Repository Structure

```text
.
├── OSEM-Axial-Force-Simulation/
│   ├── OSEM.json
│   ├── README.md
│   ├── osem_forces.py
│   ├── osem_theory.pdf
│   ├── requirements.md
│   └── requirements.txt
│
├── Magnetization-Determination/
│   ├── Results/
│   │   ├── amperian_fit_plot.html
│   │   ├── amperian_fit_results.txt
│   ├── README.md
│   ├── fit_magnetization.py
│   ├── obs.csv
│   ├── params.json
│   ├── requirements.txt
│   └── requirements_magnetization.md
│
└── .github/
    └── workflows/
        └── ...
```

### `OSEM-Axial-Force-Simulation`

Contains `osem_forces.py`, a Python implementation for calculating the **axial force produced by the OSEM voice coil actuator** as a function of magnet position and actuator parameters.

The simulation follows the electromagnetic and cylindrical-geometry formulation developed for the project, providing a computational route for analysing the coil–magnet interaction and force profile.

### `Magnetization-Determination`

Contains `fit_magnetization.py`, which analyses experimentally measured **axial repulsive force between cylindrical permanent magnets**.

The code applies regression and statistical analysis to the force–separation observations and determines the magnetization (M) using the **Amperian surface-current model** for uniformly magnetized cylindrical magnets.

## Automated Testing

The `.github/workflows/` directory contains GitHub Actions workflows that automatically test the project code.

The workflows provide execution feedback and update the corresponding directory-level `README.md` files with the current test status, making the repository's computational health visible alongside each component.

## Project Context

The computational work accompanies an experimental study of OSEM actuator force characteristics, with particular emphasis on determining the magnetization of the permanent magnets and improving agreement between theoretical predictions and experimental observations.

For the underlying theory, experimental methodology, and detailed results, refer to the project report.

---

**Language:** Python
**Domain:** Electromagnetism · Magnetic Actuation · Numerical Simulation · Experimental Data Analysis
**Application:** LIGO OSEM / BOSEM Characterization
