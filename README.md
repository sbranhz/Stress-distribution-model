# Stress Distribution Model — Trapezoidal Beam

Numerical model built in Python that calculates and visualizes the full stress distribution across a trapezoidal cross-section beam under combined loading conditions.

## What it does

Given a trapezoidal beam with configurable geometry and loading, the model computes:

- **σyy** — Normal stress along the vertical axis
- **σxx** — Normal stress along the horizontal axis
- **τxy** — Shear stress
- **σ_max / σ_min** — Principal stresses
- **τ_max** — Maximum shear stress

Results are displayed as 2D contour plots across the beam's cross-section.

## Configuration

All parameters are defined at the top of the script:

| Parameter | Description | Default |
|-----------|-------------|---------|
| `ac` | Top width (mm) | 50 |
| `ap` | Bottom width (mm) | 70 |
| `h` | Beam height (mm) | 100 |
| `b` | Thickness (mm) | 1 |
| `f` | Applied distributed load (N/mm) | 1.04e6 |
| `gc` | Specific weight of material (N/mm³) | 2.3 × g |
| `r` | Mesh resolution (points) | 10 |
| `body` | Include self-weight | `True` |
| `uni` | Uniform load (`True`) or triangular load (`False`) | `False` |

## How to run

```bash
pip install numpy matplotlib
python Modelo_Mecanica_Solidos.py
```

The script will display all 6 stress plots and save them as `stress_distribution.png`.

## Tech stack

- **Python 3**
- **NumPy** — numerical computation and matrix operations
- **Matplotlib** — contour plot visualization

## Background

Developed as part of a Solid Mechanics course at Universidad EAFIT. Applies beam theory and stress transformation equations to model how internal forces distribute across non-uniform cross-sections — relevant to structural analysis in civil and mechanical engineering.
