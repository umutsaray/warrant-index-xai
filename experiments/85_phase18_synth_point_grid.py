#!/usr/bin/env python3
"""Phase 18 CAL-seed synthetic point grid and matched-scale comparison."""
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

KS = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)
KAPPAS = (0.1, 0.25, 0.5)
CAL_SEEDS = tuple(range(18_000_001, 18_000_061))
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def simulate(seed: int, syn) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n, p, t = 700, 30, 50
    beta0 = np.r_[np.full(8, 0.05), np.zeros(6), np.full(6, 0.008), np.full(5, 0.05), np.full(5, 0.03)]
    beta1 = beta0.copy()
    beta1[20:25] = -0.05
    covariance = np.full((p, p), 0.3)
    np.fill_diagonal(covariance, 1)
    x = [rng.multivariate_normal(np.zeros(p), covariance, n) for _ in range(2)]
    y = [x[0] @ beta0 + rng.normal(0, 0.1, n), x[1] @ beta1 + rng.normal(0, 0.1, n)]
    phi = np.empty((2, t, n, p))
    for site in range(2):
        for replica in range(t):
            sampled = rng.integers(0, n, n)
            coefficient = syn.ridge(x[site][sampled], y[site][sampled])
            phi[site, replica] = coefficient * (x[site] - x[site].mean(0))
    return phi


def main() -> None:
    started = time.perf_counter()
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_grid_syn")
    breakdown = load("84_phase18_synth_breakdown.py", "phase18_grid_breakdown")
    shared = load("12_rsf_certification.py", "phase18_grid_shared")
    pilot = load("08_certification_pilot.py", "phase18_grid_pilot")
    checkpoint = OUT / "PHASE18_ADIM2_POINT_GRID_CHECKPOINT.csv"
    matched_checkpoint = OUT / "PHASE18_ADIM2_MATCHED_SCALE_CHECKPOINT.csv"
    rows = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    matched_rows = pd.read_csv(matched_checkpoint).to_dict("records") if matched_checkpoint.exists() else []
    done = {int(row["seed"]) for row in matched_rows}

    for seed in CAL_SEEDS:
        if seed in done:
            continue
        phi = simulate(seed, syn)
        magnitude = np.abs(phi[0]).mean(axis=(0, 1))
        scale = float(np.median(magnitude[magnitude > 0]))
        n_by_kappa = {
            kappa: np.stack([
                np.minimum(1.0, np.abs(phi[site]).mean(axis=1) / (kappa * scale)).mean(axis=0)
                for site in range(2)
            ]).min(axis=0)
            for kappa in KAPPAS
        }
        v1 = breakdown.point_components(phi, 0.5, 0.5, 3.0, 0.02, shared, pilot)
        for feature in range(30):
            values = v1[feature]
            rows.append(
                {
                    "seed": seed, "feature": feature, "truth_class": CLASSES[feature],
                    "variant": "v1_absolute", "k": np.nan, "rho": np.nan, "kappa": np.nan,
                    "s_m_sim": scale, **{name: values[q] for q, name in enumerate(COMPONENTS)},
                    "N": np.nan, "WI": values.min(), "point_adequate": values.min() >= 0.8,
                    "bottleneck": COMPONENTS[int(values.argmin())],
                }
            )
        for k in KS:
            values5 = breakdown.point_components(phi, k * scale, k * scale, 6 * k * scale, scale, shared, pilot)
            for feature in range(30):
                values = values5[feature]
                rows.append(
                    {
                        "seed": seed, "feature": feature, "truth_class": CLASSES[feature],
                        "variant": "v2_WI5", "k": k, "rho": 1.0, "kappa": np.nan,
                        "s_m_sim": scale, **{name: values[q] for q, name in enumerate(COMPONENTS)},
                        "N": np.nan, "WI": values.min(), "point_adequate": values.min() >= 0.8,
                        "bottleneck": COMPONENTS[int(values.argmin())],
                    }
                )
                for kappa in KAPPAS:
                    values6 = np.r_[values, n_by_kappa[kappa][feature]]
                    names6 = COMPONENTS + ("N",)
                    rows.append(
                        {
                            "seed": seed, "feature": feature, "truth_class": CLASSES[feature],
                            "variant": "v2_WI6", "k": k, "rho": 1.0, "kappa": kappa,
                            "s_m_sim": scale, **{name: values[q] for q, name in enumerate(COMPONENTS)},
                            "N": n_by_kappa[kappa][feature], "WI": values6.min(),
                            "point_adequate": values6.min() >= 0.8,
                            "bottleneck": names6[int(values6.argmin())],
                        }
                    )

        k_v1 = 0.5 / scale
        rho_v1 = 0.02 / scale
        matched = breakdown.point_components(phi, 0.5, 0.5, 3.0, 0.02, shared, pilot)
        np.testing.assert_allclose(matched, v1, atol=1e-12, rtol=0)
        for feature in range(30):
            values5 = matched[feature]
            values6 = np.r_[values5, n_by_kappa[0.25][feature]]
            matched_rows.append(
                {
                    "seed": seed, "feature": feature, "truth_class": CLASSES[feature],
                    "s_m_sim": scale, "k_v1": k_v1, "rho_v1": rho_v1,
                    "v1_WI5": v1[feature].min(), "matched_v2_WI5": values5.min(),
                    "matched_v2_WI6": values6.min(),
                    "v1_point": v1[feature].min() >= 0.8,
                    "matched_v2_WI5_point": values5.min() >= 0.8,
                    "matched_v2_WI6_point": values6.min() >= 0.8,
                }
            )
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
        pd.DataFrame(matched_rows).to_csv(matched_checkpoint, index=False)
        print(f"CAL point grid {len(done) + 1}/60 seed={seed}", flush=True)
        done.add(seed)

    result = pd.DataFrame(rows)
    matched = pd.DataFrame(matched_rows)
    result.to_csv(OUT / "PHASE18_ADIM2_POINT_GRID.csv", index=False)
    matched.to_csv(OUT / "PHASE18_ADIM2_MATCHED_SCALE.csv", index=False)
    print(result.groupby(["variant", "k", "kappa", "truth_class"], dropna=False)["point_adequate"].mean().to_string())
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
