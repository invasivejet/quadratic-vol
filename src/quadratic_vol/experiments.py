from __future__ import annotations

from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.optimize import brentq
from scipy.stats import norm

from .finite_difference import crank_nicolson_call_grid
from .model import (
    QuadraticVolParams,
    default_regimes,
    diffusion,
    effective_potential,
    financial_potential,
    geometric_potential,
    sigma,
)
from .spectral import SpectralPricer


def configure_plots() -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
            "figure.figsize": (5.2, 3.7),
            "lines.linewidth": 1.55,
            "axes.linewidth": 0.8,
            "savefig.dpi": 300,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.family": "serif",
            "mathtext.fontset": "cm",
        }
    )


def save_csv(path: Path, header: str, data: np.ndarray) -> None:
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def black_scholes_call(k_log: float, tau: float, r: float, vol: float) -> float:
    strike = float(np.exp(k_log))
    st = vol * np.sqrt(tau)
    d1 = (np.log(1.0 / strike) + (r + 0.5 * vol * vol) * tau) / st
    d2 = d1 - st
    return norm.cdf(d1) - strike * np.exp(-r * tau) * norm.cdf(d2)


def implied_vol(price: float, k_log: float, tau: float, r: float) -> float:
    strike = float(np.exp(k_log))
    lower = max(1.0 - strike * np.exp(-r * tau), 0.0)
    if not np.isfinite(price) or price <= lower + 1.0e-12 or price >= 1.0 - 1.0e-12:
        return np.nan

    def objective(vol: float) -> float:
        return black_scholes_call(k_log, tau, r, vol) - price

    try:
        return float(brentq(objective, 1.0e-6, 5.0))
    except ValueError:
        return np.nan


def fd_price_at(x_grid: np.ndarray, params: QuadraticVolParams, tau: float, k: float, steps: int) -> float:
    values = crank_nicolson_call_grid(x_grid, params, tau, k, n_steps=steps)
    return float(np.interp(0.0, x_grid, values))


def run_regime_summary(tables_dir: Path) -> None:
    rows = []
    x = np.linspace(-3.0, 3.0, 601)
    for p in default_regimes():
        sig = sigma(x, p)
        rows.append(
            [
                p.alpha,
                p.beta,
                p.gamma,
                p.discriminant,
                p.regime == "hyperbolic",
                p.regime == "parabolic",
                p.regime == "elliptic",
                np.min(sig),
                np.max(sig),
                np.min(diffusion(x, p)),
                np.max(diffusion(x, p)),
            ]
        )
    save_csv(
        tables_dir / "regime_summary.csv",
        "alpha,beta,gamma,discriminant,is_hyperbolic,is_parabolic,is_elliptic,min_sigma,max_sigma,min_D,max_D",
        np.asarray(rows, dtype=np.float64),
    )


def run_regime_potentials(figures_dir: Path, tables_dir: Path) -> None:
    x = np.linspace(-3.0, 3.0, 801)
    rows = []
    fig, axes = plt.subplots(2, 1, figsize=(5.8, 5.1), sharex=True)
    colors = plt.cm.magma(np.linspace(0.18, 0.88, 3))
    for color, p in zip(colors, default_regimes()):
        sig = sigma(x, p)
        u = effective_potential(x, p)
        geom = geometric_potential(x, p)
        fin = financial_potential(x, p)
        label = f"{p.regime}, $\\Delta={p.discriminant:.3g}$"
        axes[0].plot(x, sig, color=color, label=label)
        axes[1].plot(x, u, color=color, label=label)
        rows.extend(
            np.column_stack(
                [
                    np.full_like(x, p.discriminant),
                    np.full_like(x, {"elliptic": -1, "parabolic": 0, "hyperbolic": 1}[p.regime]),
                    x,
                    sig,
                    diffusion(x, p),
                    geom,
                    fin,
                    u,
                ]
            )
        )
    axes[0].set_ylabel("$\\sigma(x)$")
    axes[1].set_xlabel("Log state $x$")
    axes[1].set_ylabel("$U_{\\rm eff}(x)$")
    for ax in axes:
        ax.grid(True, linestyle="--", alpha=0.55)
        ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(figures_dir / "regime_potential_landscapes.pdf", bbox_inches="tight")
    plt.close(fig)

    save_csv(
        tables_dir / "regime_potential_landscapes.csv",
        "discriminant,regime_code,x,sigma,diffusion,U_geom,U_fin,U_eff",
        np.asarray(rows, dtype=np.float64),
    )


def run_spectrum_study(figures_dir: Path, tables_dir: Path) -> None:
    x_grid = np.linspace(-3.0, 3.0, 501)
    rows = []
    fig, ax = plt.subplots()
    for p in default_regimes():
        pricer = SpectralPricer(x_grid, p)
        spec = pricer.spectrum(24)
        modes = np.arange(spec.eigenvalues.size)
        ax.plot(modes, spec.eigenvalues, marker="o", markersize=3, label=f"{p.regime}, $\\Delta={p.discriminant:.3g}$")
        rows.extend(np.column_stack([np.full_like(modes, p.discriminant, dtype=float), modes, spec.eigenvalues]))
    ax.set_xlabel("Mode index")
    ax.set_ylabel("Eigenvalue")
    ax.grid(True, linestyle="--", alpha=0.55)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(figures_dir / "spectrum_by_regime.pdf", bbox_inches="tight")
    plt.close(fig)
    save_csv(tables_dir / "spectrum_by_regime.csv", "discriminant,mode,eigenvalue", np.asarray(rows))


def run_smile_by_regime(figures_dir: Path, tables_dir: Path) -> None:
    x_grid = np.linspace(-3.0, 3.0, 501)
    strikes = np.linspace(-1.0, 1.0, 25)
    tau = 1.0
    rows = []
    fig, ax = plt.subplots()
    colors = plt.cm.magma(np.linspace(0.18, 0.88, 3))
    for color, p in zip(colors, default_regimes()):
        ivs = []
        for k in strikes:
            price = fd_price_at(x_grid, p, tau, float(k), steps=240)
            iv = implied_vol(price, float(k), tau, p.r)
            ivs.append(iv)
            rows.append([p.discriminant, {"elliptic": -1, "parabolic": 0, "hyperbolic": 1}[p.regime], k, price, iv])
        ax.plot(strikes, ivs, marker="o", markersize=3, color=color, label=f"{p.regime}, $\\Delta={p.discriminant:.3g}$")
    ax.set_xlabel("$\\log(K/S_0)$")
    ax.set_ylabel("Implied volatility")
    ax.grid(True, linestyle="--", alpha=0.55)
    ax.legend(frameon=True)
    fig.tight_layout()
    fig.savefig(figures_dir / "smile_by_regime.pdf", bbox_inches="tight")
    plt.close(fig)
    save_csv(tables_dir / "smile_by_regime.csv", "discriminant,regime_code,k,fd_price,implied_vol", np.asarray(rows))


def run_quadratic_space(figures_dir: Path, tables_dir: Path, n: int = 3000) -> None:
    rng = np.random.default_rng(123)
    alpha = rng.uniform(0.005, 0.06, n)
    beta = rng.uniform(-0.18, 0.22, n)
    gamma = rng.uniform(0.12, 0.45, n)
    delta = beta * beta - 4.0 * alpha * gamma
    x = np.linspace(-3.0, 3.0, 101)
    min_sigma = np.array([np.min(a * x * x + b * x + g) for a, b, g in zip(alpha, beta, gamma)])
    valid = min_sigma > 0.04
    hyperbolic = delta > 0
    data = np.column_stack([alpha, beta, gamma, delta, hyperbolic.astype(float), valid.astype(float)])
    save_csv(tables_dir / "quadratic_surface_space.csv", "alpha,beta,gamma,discriminant,is_hyperbolic,is_positive_on_grid", data)

    fig = plt.figure(figsize=(6.2, 4.7))
    ax = fig.add_subplot(111, projection="3d")
    rejected = ~(hyperbolic & valid)
    accepted = hyperbolic & valid
    ax.scatter(alpha[rejected], beta[rejected], gamma[rejected], c="lightgray", s=4, alpha=0.22)
    sc = ax.scatter(alpha[accepted], beta[accepted], gamma[accepted], c=delta[accepted], cmap="magma", s=10, alpha=0.9)
    fig.colorbar(sc, ax=ax, shrink=0.72, pad=0.1, label="$\\Delta$")
    ax.set_xlabel("$\\alpha$")
    ax.set_ylabel("$\\beta$")
    ax.set_zlabel("$\\gamma$")
    ax.view_init(elev=24, azim=40)
    fig.tight_layout()
    fig.savefig(figures_dir / "quadratic_surface_space.pdf", bbox_inches="tight")
    plt.close(fig)


def run_speed_benchmark(tables_dir: Path) -> None:
    x_grid = np.linspace(-3.0, 3.0, 501)
    p = default_regimes()[-1]
    pricer = SpectralPricer(x_grid, p)
    contracts = [(tau, k) for tau in (0.25, 1.0, 2.0) for k in np.linspace(-1.0, 1.0, 21)]

    t0 = perf_counter()
    fd_prices = [fd_price_at(x_grid, p, tau, float(k), steps=220) for tau, k in contracts]
    fd_seconds = perf_counter() - t0

    t0 = perf_counter()
    pricer.spectrum(96)
    spectral_precompute = perf_counter() - t0

    t0 = perf_counter()
    spectral_prices = [pricer.price_at(0.0, tau, float(k), 96) for tau, k in contracts]
    spectral_online = perf_counter() - t0
    max_abs_gap = float(np.max(np.abs(np.asarray(fd_prices) - np.asarray(spectral_prices))))

    data = np.array(
        [
            [
                len(contracts),
                fd_seconds,
                spectral_precompute,
                spectral_online,
                spectral_precompute + spectral_online,
                fd_seconds / spectral_online,
                fd_seconds / (spectral_precompute + spectral_online),
                max_abs_gap,
            ]
        ]
    )
    save_csv(
        tables_dir / "speed_benchmark.csv",
        "n_contracts,fd_seconds,spectral_precompute_seconds,spectral_online_seconds,spectral_first_run_seconds,online_speedup,first_run_speedup,max_abs_price_gap",
        data,
    )


def run_all(output_dir: Path) -> None:
    configure_plots()
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    run_regime_summary(tables_dir)
    run_regime_potentials(figures_dir, tables_dir)
    run_spectrum_study(figures_dir, tables_dir)
    run_smile_by_regime(figures_dir, tables_dir)
    run_quadratic_space(figures_dir, tables_dir)
    run_speed_benchmark(tables_dir)


if __name__ == "__main__":
    run_all(Path("results"))
