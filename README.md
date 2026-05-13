# Quadratic Volatility Numerics

Quick, lightweight replication of the quadratic normal volatility pricing/spectral logic.

The goal is not a full paper scaffold. This directory focuses on:

- discriminant regimes for $\sigma(x)=\alpha x^2+\beta x+\gamma$;
- self-adjoint operator spectra and effective potential landscapes;
- finite-difference reference pricing;
- cached spectral online pricing diagnostics;
- runtime tables for computational agility.

## Run

```bash
cd /home/joelasaucedo/Development/jetbundle/quadratic-vol
python scripts/run_experiments.py
```

Outputs are written to:

- `results/figures/`
- `results/tables/`

Core theory notes are in `theory.md`.

The quick numerical report is in `NUMERICS_REPORT.md`.

Current figure outputs:

- `results/figures/regime_potential_landscapes.pdf`
- `results/figures/spectrum_by_regime.pdf`
- `results/figures/smile_by_regime.pdf`
- `results/figures/quadratic_surface_space.pdf`

Current benchmark table:

- `results/tables/speed_benchmark.csv`
