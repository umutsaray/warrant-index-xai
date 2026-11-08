#!/usr/bin/env python3
"""Phase 18 preregistered k-star selection, EVAL validation, and real-data transfer."""
from __future__ import annotations

import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]
EVAL_SEEDS = tuple(range(18_010_001, 18_010_061))
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    started = time.perf_counter()
    cal = pd.read_csv(OUT / "PHASE18_ADIM2_POINT_GRID.csv")
    eligible = cal[(cal.variant == "v2_WI6") & (cal.kappa == 0.25) & (cal.truth_class == "B")]
    rates = eligible.groupby("k").point_adequate.mean().sort_index()
    candidates = rates[rates <= 0.05]
    k_star = float(candidates.index.min()) if len(candidates) else np.nan
    if not np.isfinite(k_star):
        rows = [{"dataset": "selection", "status": "rule_not_satisfied", "k_star": np.nan}]
        pd.DataFrame(rows).to_csv(OUT / "PHASE18_ADIM4_KSTAR.csv", index=False)
        print("k-star rule not satisfied")
        return

    rows: list[dict[str, object]] = [
        {
            "dataset": "CAL_summary", "status": "selected" if k == k_star else "candidate",
            "seed": np.nan, "model": np.nan, "feature": np.nan, "truth_class": "B",
            "k": k, "kappa": 0.25, "k_star": k_star,
            "WI6": np.nan, "point_adequate": np.nan,
            "class_point_rate": rate, "v1_BH": np.nan,
        }
        for k, rate in rates.items()
    ]
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_kstar_syn")
    grid = load("85_phase18_synth_point_grid.py", "phase18_kstar_grid")
    breakdown = load("84_phase18_synth_breakdown.py", "phase18_kstar_breakdown")
    shared = load("12_rsf_certification.py", "phase18_kstar_shared")
    pilot = load("08_certification_pilot.py", "phase18_kstar_pilot")
    checkpoint = OUT / "PHASE18_ADIM4_KSTAR_EVAL_CHECKPOINT.csv"
    eval_rows = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {int(row["seed"]) for row in eval_rows}
    for seed in EVAL_SEEDS:
        if seed in done:
            continue
        phi = grid.simulate(seed, syn)
        magnitude = np.abs(phi[0]).mean(axis=(0, 1))
        scale = float(np.median(magnitude[magnitude > 0]))
        points = breakdown.point_components(
            phi, k_star * scale, k_star * scale, 6 * k_star * scale, scale, shared, pilot
        )
        n_point = np.stack([
            np.minimum(1.0, np.abs(phi[site]).mean(axis=1) / (0.25 * scale)).mean(axis=0)
            for site in range(2)
        ]).min(axis=0)
        for feature in range(30):
            wi6 = min(points[feature].min(), n_point[feature])
            eval_rows.append(
                {
                    "dataset": "EVAL_feature", "status": "held_out",
                    "seed": seed, "model": np.nan, "feature": feature,
                    "truth_class": CLASSES[feature], "k": k_star, "kappa": 0.25,
                    "k_star": k_star, "WI6": wi6, "point_adequate": wi6 >= 0.8,
                    "class_point_rate": np.nan, "v1_BH": np.nan,
                }
            )
        pd.DataFrame(eval_rows).to_csv(checkpoint, index=False)
        done.add(seed)
        print(f"EVAL k-star {len(done)}/60", flush=True)
    rows.extend(eval_rows)

    real = pd.read_csv(OUT / "PHASE18_ADIM5_KPROFILE.csv")
    real = real[(real.k == k_star) & (real.kappa == 0.25)].copy()
    frozen = pd.read_csv(OUT / "PHASE18_ADIM5_KENTER_81.csv")[["model", "feature", "v1_BH"]]
    real = real.merge(frozen, on=["model", "feature"], validate="one_to_one")
    for item in real.itertuples(index=False):
        rows.append(
            {
                "dataset": "real_fold0_feature", "status": "exploratory_transfer",
                "seed": np.nan, "model": item.model, "feature": item.feature,
                "truth_class": np.nan, "k": k_star, "kappa": 0.25,
                "k_star": k_star, "WI6": item.WI6,
                "point_adequate": item.point_adequate6,
                "class_point_rate": np.nan, "v1_BH": item.v1_BH,
            }
        )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "PHASE18_ADIM4_KSTAR.csv", index=False)
    eval_frame = pd.DataFrame(eval_rows)
    print("k_star", k_star, "CAL B rates", rates.to_dict())
    print("EVAL", eval_frame.groupby("truth_class").point_adequate.mean().to_dict())
    print("real point", real.groupby("model").point_adequate6.sum().to_dict())
    print("v1 BH retained", real[real.v1_BH].groupby("model").point_adequate6.sum().to_dict())
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
