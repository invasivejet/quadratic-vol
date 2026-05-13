from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh_tridiagonal
from scipy.integrate import trapezoid

from .model import (
    QuadraticVolParams,
    diffusion,
    effective_potential,
    gauge_phi_numeric,
    normalized_call_payoff,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class Spectrum:
    eigenvalues: FloatArray
    eigenvectors: FloatArray
    interior_vectors: FloatArray


class SpectralPricer:
    def __init__(self, x_grid: FloatArray, params: QuadraticVolParams):
        self.x = np.asarray(x_grid, dtype=np.float64)
        self.params = params
        self._cache: dict[int, Spectrum] = {}

    def _diagonals(self) -> tuple[FloatArray, FloatArray]:
        dx = float(self.x[1] - self.x[0])
        d = diffusion(self.x, self.params)
        u = effective_potential(self.x, self.params)
        d_half = 0.5 * (d[:-1] + d[1:])
        offdiag = -d_half[1:-1] / dx**2
        diag = (d_half[:-1] + d_half[1:]) / dx**2 + u[1:-1]
        return diag, offdiag

    def spectrum(self, n_modes: int) -> Spectrum:
        if n_modes in self._cache:
            return self._cache[n_modes]
        diag, offdiag = self._diagonals()
        eigenvalues, interior = eigh_tridiagonal(
            diag, offdiag, select="i", select_range=(0, n_modes - 1)
        )
        full = np.zeros((n_modes, self.x.size), dtype=np.float64)
        full[:, 1:-1] = interior.T
        norms = np.sqrt(trapezoid(full * full, self.x, axis=1))
        full /= norms[:, None]
        spec = Spectrum(eigenvalues=eigenvalues, eigenvectors=full, interior_vectors=interior)
        self._cache[n_modes] = spec
        return spec

    def price_grid(self, tau: float, strike_log_moneyness: float, n_modes: int) -> FloatArray:
        spec = self.spectrum(n_modes)
        x = self.x
        d = diffusion(x, self.params)
        phi = gauge_phi_numeric(x, self.params)
        payoff = normalized_call_payoff(x, strike_log_moneyness)
        initial = np.exp(-phi) * payoff / np.sqrt(d)
        q = spec.interior_vectors
        coeffs = q.T @ initial[1:-1]
        weights = coeffs * np.exp(-spec.eigenvalues * tau)
        evolved = np.zeros_like(x)
        evolved[1:-1] = q @ weights
        return np.exp(phi) * np.sqrt(d) * evolved

    def price_at(self, x0: float, tau: float, strike_log_moneyness: float, n_modes: int) -> float:
        return float(np.interp(x0, self.x, self.price_grid(tau, strike_log_moneyness, n_modes)))
