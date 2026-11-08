#!/usr/bin/env python3
"""Phase 18 DGP scale and component-level decomposition of Phase 17 synthetic data."""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
OUT.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from _shared import CALIBRATION

COMPONENTS = ("F", "S", "U", "G_claim", "R")
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def point_components(phi: np.ndarray, a: float, b: float, u: float, eta: float, shared, pilot) -> np.ndarray:
    original = {key: CALIBRATION[key] for key in ("a", "b", "u", "eta")}
    CALIBRATION.update(a=a, b=b, u=u, eta=eta)
    try:
        sites, t, _, p = phi.shape
        f_by_rep = 1 - (np.minimum(1.0, np.abs((phi - phi) / a)) ** 2).mean(axis=2)
        f_point = f_by_rep.mean(axis=1).min(axis=0)
        lower, upper = pilot.importance_ci(phi, CALIBRATION["alpha_U"])
        u_by_rep = 1 - np.minimum(1.0, (upper - lower) / u)
        u_point = u_by_rep.mean(axis=1).min(axis=0)
        importance = np.abs(phi).mean(axis=2)
        targets = np.arange(p)
        s_loss = pilot.pair_losses(phi, targets, a)
        r_loss = pilot.order_pair_losses(importance, targets, eta)
        s_point = (1 - s_loss.sum(axis=(-2, -1)) / (t * (t - 1))).min(axis=0)
        r_point = (1 - r_loss.sum(axis=(-2, -1)) / (t * (t - 1))).min(axis=0)
        difference = importance[1] - importance[0]
        g_magnitude = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / b)
        full_weights = np.ones((1, t), dtype=np.int16)
        g_signed = shared.signed_scores(phi, full_weights)[0]
        g_claim = np.minimum(g_magnitude, g_signed)
        return np.column_stack([f_point, s_point, u_point, g_claim, r_point])
    finally:
        CALIBRATION.update(original)


def main() -> None:
    started = time.perf_counter()
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_breakdown_syn")
    shared = load("12_rsf_certification.py", "phase18_breakdown_shared")
    pilot = load("08_certification_pilot.py", "phase18_breakdown_pilot")
    rng = np.random.default_rng(20261108)
    n, p, t, draws = 700, 30, 50, 999
    beta0 = np.r_[np.full(8, 0.05), np.zeros(6), np.full(6, 0.008), np.full(5, 0.05), np.full(5, 0.03)]
    beta1 = beta0.copy()
    beta1[20:25] = -0.05
    covariance = np.full((p, p), 0.3)
    np.fill_diagonal(covariance, 1)
    component_rows = []
    scale_rows = []
    frozen = pd.read_csv(ROOT / "results/phase17/PHASE17_ADIM5_SYNTH.csv").set_index(["repeat", "feature"])

    for repeat in range(20):
        x = [rng.multivariate_normal(np.zeros(p), covariance, n) for _ in range(2)]
        y = [x[0] @ beta0 + rng.normal(0, 0.1, n), x[1] @ beta1 + rng.normal(0, 0.1, n)]
        phi = np.empty((2, t, n, p))
        for site in range(2):
            for replica in range(t):
                sampled = rng.integers(0, n, n)
                coefficient = syn.ridge(x[site][sampled], y[site][sampled])
                phi[site, replica] = coefficient * (x[site] - x[site].mean(0))
        # Preserve the exact Phase 17 RNG chain even though point calculations
        # do not otherwise need the inner bootstrap indices.
        rng.integers(0, t, size=(draws, t))
        magnitude = np.abs(phi[0]).mean(axis=(0, 1))
        scale = float(np.median(magnitude[magnitude > 0]))
        v1 = point_components(phi, 0.5, 0.5, 3.0, 0.02, shared, pilot)
        v2 = point_components(phi, scale, scale, 6 * scale, scale, shared, pilot)
        n_point = np.stack([
            np.minimum(1.0, np.abs(phi[site]).mean(axis=1) / (0.25 * scale)).mean(axis=0)
            for site in range(2)
        ]).min(axis=0)
        np.testing.assert_allclose(v1.min(axis=1), frozen.loc[repeat, "v1_WI"].to_numpy(), atol=1e-12, rtol=0)
        np.testing.assert_allclose(np.column_stack([v2, n_point]).min(axis=1), frozen.loc[repeat, "v2_WI6"].to_numpy(), atol=1e-12, rtol=0)
        for truth_class in "ABCDE":
            mask = CLASSES == truth_class
            scale_rows.append(
                {
                    "dataset": "synthetic",
                    "repeat": repeat,
                    "truth_class": truth_class,
                    "s_m_sim": scale,
                    "class_median_M_j": float(np.median(magnitude[mask])),
                    "v1_a_over_s_m": 0.5 / scale,
                    "v1_u_over_s_m": 3.0 / scale,
                    "v1_eta_over_s_m": 0.02 / scale,
                    "v2_a_over_s_m": 1.0,
                    "v2_u_over_s_m": 6.0,
                    "v2_eta_over_s_m": 1.0,
                }
            )
        for feature in range(p):
            values_v1 = v1[feature]
            values_v2 = np.r_[v2[feature], n_point[feature]]
            component_rows.append(
                {
                    "repeat": repeat,
                    "feature": feature,
                    "truth_class": CLASSES[feature],
                    "M_j": magnitude[feature],
                    "s_m_sim": scale,
                    **{f"v1_{name}": values_v1[q] for q, name in enumerate(COMPONENTS)},
                    "v1_WI5": values_v1.min(),
                    "v1_bottleneck": COMPONENTS[int(values_v1.argmin())],
                    **{f"v2_{name}": values_v2[q] for q, name in enumerate(COMPONENTS + ("N",))},
                    "v2_WI6": values_v2.min(),
                    "v2_bottleneck": (COMPONENTS + ("N",))[int(values_v2.argmin())],
                }
            )
        print(f"breakdown {repeat + 1}/20", flush=True)

    real_scales = {"CoxPH-v2": 0.023934215222165295, "RSF": 0.01707069437342852, "DeepSurv": 0.012307440008774696}
    for model, scale in real_scales.items():
        scale_rows.append(
            {
                "dataset": "real_fold0",
                "repeat": np.nan,
                "truth_class": model,
                "s_m_sim": scale,
                "class_median_M_j": np.nan,
                "v1_a_over_s_m": 0.5 / scale,
                "v1_u_over_s_m": 3.0 / scale,
                "v1_eta_over_s_m": 0.02 / scale,
                "v2_a_over_s_m": 1.0,
                "v2_u_over_s_m": 6.0,
                "v2_eta_over_s_m": 1.0,
            }
        )
    pd.DataFrame(scale_rows).to_csv(OUT / "PHASE18_ADIM1_DGP_SCALE.csv", index=False)
    pd.DataFrame(component_rows).to_csv(OUT / "PHASE18_ADIM1_COMPONENT_BREAKDOWN.csv", index=False)
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
