# Requirements

## Python

Python 3.8 or newer (the code uses `dataclasses` and f-strings).

## Dependencies

| Package | Minimum version | Used for |
|---|---|---|
| [`numpy`](https://numpy.org/) | 1.20 | Gauss-Legendre quadrature nodes/weights (`numpy.polynomial.legendre`), array math |
| [`plotly`](https://plotly.com/python/) | 5.0 | Interactive HTML plot of `Fz` vs. `z` |

Everything else used (`argparse`, `json`, `os`, `time`, `dataclasses`, `typing`) is part of the Python standard library.

## Install

```bash
pip install numpy plotly
```

Or, saving the block below as `requirements.txt`:

```
numpy>=1.20
plotly>=5.0
```

```bash
pip install -r requirements.txt
```

## Notes

- No GPU, compiled extensions, or external solvers are required — this is pure NumPy/Plotly.
- Runtime scales with the product of all five Gauss-Legendre orders (`n_theta * n_r1 * n_z1 * n_r2 * n_z2`) and the number of z-scan points, so increasing quadrature order or scan resolution in `OSEM.json` trades accuracy for runtime.
- The interactive plot is written as a self-contained `.html` file that loads the Plotly.js library from a CDN (`include_plotlyjs='cdn'`), keeping the output file small; an internet connection is needed the first time it's opened in a browser (or `include_plotlyjs=True` can be set in the code to embed the library offline).
