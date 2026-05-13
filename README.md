# Quadratic Volatility Numerics

Lightweight numerical scaffold for quadratic normal volatility pricing and spectral diagnostics.

This repository is intended as a clean handoff point for continuing the numerical work. It focuses on:

- discriminant regimes for $\sigma(x)=\alpha x^2+\beta x+\gamma$;
- self-adjoint operator spectra and effective potential landscapes;
- finite-difference reference pricing;
- cached spectral online pricing diagnostics;
- runtime tables for computational agility.

## Run

```bash
pip install -r requirements.txt
python scripts/run_experiments.py
```

Outputs are written to:

- `results/figures/`
- `results/tables/`

Core theory notes are in `theory.md`.

A short numerical handoff report is in `NUMERICS_REPORT.md`.

Current figure outputs:

- `results/figures/regime_potential_landscapes.pdf`
- `results/figures/spectrum_by_regime.pdf`
- `results/figures/smile_by_regime.pdf`
- `results/figures/quadratic_surface_space.pdf`

Current benchmark table:

- `results/tables/speed_benchmark.csv`
