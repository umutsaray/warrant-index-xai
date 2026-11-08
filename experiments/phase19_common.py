"""Shared deterministic synthetic helpers for Phase 19 diagnostics."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)
BETA0 = np.r_[np.full(8, 0.05), np.zeros(6), np.full(6, 0.008), np.full(5, 0.05), np.full(5, 0.03)]
BETA1 = BETA0.copy()
BETA1[20:25] = -0.05


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def ridge(x: np.ndarray, y: np.ndarray, lam: float = 0.1) -> np.ndarray:
    return np.linalg.solve(x.T @ x + lam * np.eye(x.shape[1]), x.T @ y)


def simulate(
    seed: int,
    *,
    n: int = 700,
    t: int = 50,
    skew_d: bool = False,
) -> dict[str, np.ndarray | float]:
    """Generate one independent deterministic Phase 19 synthetic replicate."""
    rng = np.random.default_rng(seed)
    p = 30
    covariance = np.full((p, p), 0.3)
    np.fill_diagonal(covariance, 1.0)
    x = [rng.multivariate_normal(np.zeros(p), covariance, n) for _ in range(2)]
    if skew_d:
        # Replace only D columns by standardized correlated log-normal margins.
        for site in range(2):
            raw = np.exp(x[site][:, 20:25])
            x[site][:, 20:25] = (raw - raw.mean(axis=0)) / raw.std(axis=0, ddof=0)
    beta = [BETA0, BETA1]
    noise_sd = 0.1
    y = [x[d] @ beta[d] + rng.normal(0.0, noise_sd, n) for d in range(2)]
    phi = np.empty((2, t, n, p))
    coefficients = np.empty((2, t, p))
    for site in range(2):
        centered = x[site] - x[site].mean(axis=0)
        for replica in range(t):
            sampled = rng.integers(0, n, n)
            coefficient = ridge(x[site][sampled], y[site][sampled])
            coefficients[site, replica] = coefficient
            phi[site, replica] = coefficient * centered
    return {
        "x0": x[0],
        "x1": x[1],
        "y0": y[0],
        "y1": y[1],
        "phi": phi,
        "coefficients": coefficients,
        "noise_sd": noise_sd,
    }


def simulate_light(seed: int, *, n: int = 700, t: int = 50) -> dict[str, np.ndarray | float]:
    """Same RNG/DGP as ``simulate`` without materializing the 4-D phi tensor."""
    rng = np.random.default_rng(seed)
    p = 30
    covariance = np.full((p, p), 0.3)
    np.fill_diagonal(covariance, 1.0)
    x = [rng.multivariate_normal(np.zeros(p), covariance, n) for _ in range(2)]
    beta = [BETA0, BETA1]
    noise_sd = 0.1
    y = [x[d] @ beta[d] + rng.normal(0.0, noise_sd, n) for d in range(2)]
    coefficients = np.empty((2, t, p))
    for site in range(2):
        for replica in range(t):
            sampled = rng.integers(0, n, n)
            coefficients[site, replica] = ridge(x[site][sampled], y[site][sampled])
    return {"x0": x[0], "x1": x[1], "y0": y[0], "y1": y[1],
            "coefficients": coefficients, "noise_sd": noise_sd}


def scale_and_n(phi: np.ndarray, kappa: float = 0.25) -> tuple[float, np.ndarray, np.ndarray]:
    magnitude = np.abs(phi[0]).mean(axis=(0, 1))
    scale = float(np.median(magnitude[magnitude > 0]))
    site_n = np.stack([
        np.minimum(1.0, np.abs(phi[site]).mean(axis=1) / (kappa * scale)).mean(axis=0)
        for site in range(2)
    ])
    return scale, magnitude, site_n.min(axis=0)


def quantile_signed(phi: np.ndarray, tolerance: np.ndarray | float) -> np.ndarray:
    means = phi.mean(axis=1)
    ordered = np.sort(means, axis=1)
    return 1.0 - np.minimum(1.0, np.abs(ordered[1] - ordered[0]) / tolerance).mean(axis=0)


def pooled_auc(positive: np.ndarray, negative: np.ndarray) -> float:
    difference = positive[:, None] - negative[None, :]
    return float((np.count_nonzero(difference > 0) + 0.5 * np.count_nonzero(difference == 0)) / difference.size)
