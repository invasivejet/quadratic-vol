from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class QuadraticVolParams:
    alpha: float
    beta: float
    gamma: float
    r: float = 0.05
    name: str = "regime"

    @property
    def discriminant(self) -> float:
        return self.beta * self.beta - 4.0 * self.alpha * self.gamma

    @property
    def regime(self) -> str:
        delta = self.discriminant
        if delta > 1.0e-12:
            return "hyperbolic"
        if delta < -1.0e-12:
            return "elliptic"
        return "parabolic"

    def roots(self) -> tuple[float, float] | None:
        delta = self.discriminant
        if delta < 0.0:
            return None
        root_delta = np.sqrt(delta)
        return (
            float((-self.beta - root_delta) / (2.0 * self.alpha)),
            float((-self.beta + root_delta) / (2.0 * self.alpha)),
        )


def as_array(x: ArrayLike) -> FloatArray:
    return np.asarray(x, dtype=np.float64)


def sigma(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    x_arr = as_array(x)
    return params.alpha * x_arr * x_arr + params.beta * x_arr + params.gamma


def sigma_prime(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    x_arr = as_array(x)
    return 2.0 * params.alpha * x_arr + params.beta


def sigma_second(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    x_arr = as_array(x)
    return np.full_like(x_arr, 2.0 * params.alpha, dtype=np.float64)


def diffusion(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    sig = sigma(x, params)
    return 0.5 * sig * sig


def diffusion_prime(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    sig = sigma(x, params)
    sig_p = sigma_prime(x, params)
    return sig * sig_p


def diffusion_second(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    sig = sigma(x, params)
    sig_p = sigma_prime(x, params)
    sig_pp = sigma_second(x, params)
    return sig_p * sig_p + sig * sig_pp


def effective_potential(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    d = diffusion(x, params)
    dp = diffusion_prime(x, params)
    ddp = diffusion_second(x, params)
    r = params.r
    return (
        0.25 * d
        + 0.5 * r
        + (r * r) / (4.0 * d)
        - (r * dp) / (2.0 * d)
        + (dp * dp) / (4.0 * d)
        - 0.5 * ddp
    )


def geometric_potential(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    d = diffusion(x, params)
    dp = diffusion_prime(x, params)
    ddp = diffusion_second(x, params)
    return (dp * dp) / (4.0 * d) - 0.5 * ddp


def financial_potential(x: ArrayLike, params: QuadraticVolParams) -> FloatArray:
    d = diffusion(x, params)
    dp = diffusion_prime(x, params)
    r = params.r
    return 0.25 * d + 0.5 * r + (r * r) / (4.0 * d) - (r * dp) / (2.0 * d)


def gauge_phi_numeric(x: FloatArray, params: QuadraticVolParams) -> FloatArray:
    """Numerical gauge with Phi(0)=0 on a sorted grid."""

    d = diffusion(x, params)
    phi_prime = 0.5 - params.r / (2.0 * d)
    phi = np.zeros_like(x)
    increments = 0.5 * (phi_prime[:-1] + phi_prime[1:]) * np.diff(x)
    phi[1:] = np.cumsum(increments)
    return phi - np.interp(0.0, x, phi)


def normalized_call_payoff(x: ArrayLike, strike_log_moneyness: float) -> FloatArray:
    x_arr = as_array(x)
    return np.maximum(np.exp(x_arr) - np.exp(strike_log_moneyness), 0.0)


def default_regimes() -> list[QuadraticVolParams]:
    return [
        QuadraticVolParams(0.0200, 0.0000, 0.2000, name="elliptic_delta_negative"),
        QuadraticVolParams(0.0100, 0.1000, 0.2500, name="parabolic_delta_zero"),
        QuadraticVolParams(0.0100, 0.1200, 0.3500, name="hyperbolic_roots_outside"),
    ]
