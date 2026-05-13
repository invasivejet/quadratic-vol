# Numerics Report

## Scope

This is a quick standalone replication of the quadratic normal volatility numerics:

$$
\sigma(x)=\alpha x^2+\beta x+\gamma,\qquad D(x)=\frac12\sigma(x)^2.
$$

It tests three finite-domain regimes:

- elliptic: $\Delta<0$;
- parabolic: $\Delta=0$ with the double root outside the main pricing window;
- hyperbolic: $\Delta>0$ with both roots outside the main pricing window.

The positivity-on-grid check matters because $D(x)>0$ is required for the gauge transform and self-adjoint operator. Hyperbolic discriminant alone is not sufficient for stable pricing if a real root enters the computational domain.

## Generated Figures

- `results/figures/regime_potential_landscapes.pdf`
- `results/figures/spectrum_by_regime.pdf`
- `results/figures/smile_by_regime.pdf`
- `results/figures/quadratic_surface_space.pdf`

## Generated Tables

- `results/tables/regime_summary.csv`
- `results/tables/regime_potential_landscapes.csv`
- `results/tables/spectrum_by_regime.csv`
- `results/tables/smile_by_regime.csv`
- `results/tables/quadratic_surface_space.csv`
- `results/tables/speed_benchmark.csv`

## Current Benchmark

For 63 contracts in the hyperbolic finite-domain regime:

- Crank--Nicolson reference loop: about `1.33s`.
- Spectral precompute: about `0.20s`.
- Cached spectral online pricing: about `0.009s`.
- Online speedup: about `149x`.
- First-run speedup including spectrum precompute: about `6.4x`.
- Max absolute price gap over the benchmark grid: about `2.9e-4`.

This is a lightweight computational-agility diagnostic, not a final production accuracy claim.

## Quadratic Surface Space

The random coefficient scan sampled 3000 quadratic surfaces. In this run:

- 511 had $\Delta>0$.
- 22 were both hyperbolic and positive on the pricing grid.

This shows why hyperbolic market states must be filtered by both discriminant and domain positivity before using the gauge/spectral machinery.

## Caveats

- The spectral solver is a finite-interval homogeneous-boundary approximation in transformed variables.
- Crank--Nicolson with financial boundary conditions is the reference pricing path.
- The Pöschl--Teller connection is represented through effective-potential and spectrum diagnostics here; analytic associated-Legendre eigenfunctions are not yet implemented in this quick scaffold.
