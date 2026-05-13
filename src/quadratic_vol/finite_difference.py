from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import diags
from scipy.sparse.linalg import factorized

from .model import QuadraticVolParams, diffusion


FloatArray = NDArray[np.float64]


def call_right_boundary(x_max: float, strike_log_moneyness: float, tau: float, r: float) -> float:
    return float(np.exp(x_max) - np.exp(strike_log_moneyness) * np.exp(-r * tau))


def crank_nicolson_call_grid(
    x_grid: FloatArray,
    params: QuadraticVolParams,
    tau: float,
    strike_log_moneyness: float,
    n_steps: int = 300,
) -> FloatArray:
    x = np.asarray(x_grid, dtype=np.float64)
    dx = float(x[1] - x[0])
    dt = tau / n_steps
    strike = float(np.exp(strike_log_moneyness))
    values = np.maximum(np.exp(x) - strike, 0.0)

    xi = x[1:-1]
    d = diffusion(xi, params)
    drift = -(d - params.r)
    lower = d / dx**2 - drift / (2.0 * dx)
    diag = -2.0 * d / dx**2 - params.r
    upper = d / dx**2 + drift / (2.0 * dx)
    n = xi.size

    operator = diags(
        [lower[1:], diag, upper[:-1]],
        offsets=[-1, 0, 1],
        shape=(n, n),
        format="csr",
    )
    identity = diags([np.ones(n)], [0], format="csr")
    lhs = identity - 0.5 * dt * operator
    rhs_matrix = identity + 0.5 * dt * operator
    solve = factorized(lhs.tocsc())

    def boundary_vec(t: float) -> FloatArray:
        vec = np.zeros(n, dtype=np.float64)
        vec[-1] = upper[-1] * call_right_boundary(float(x[-1]), strike_log_moneyness, t, params.r)
        return vec

    for step in range(n_steps):
        t0 = step * dt
        t1 = (step + 1) * dt
        rhs = rhs_matrix @ values[1:-1]
        rhs += 0.5 * dt * (boundary_vec(t0) + boundary_vec(t1))
        values[1:-1] = solve(rhs)
        values[0] = 0.0
        values[-1] = call_right_boundary(float(x[-1]), strike_log_moneyness, t1, params.r)
    return values
