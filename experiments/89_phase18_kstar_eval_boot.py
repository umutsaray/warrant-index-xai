#!/usr/bin/env python3
"""Phase 18 held-out EVAL B=999 certification at preregistered k-star."""
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
SEEDS = tuple(range(18_010_001, 18_010_021))
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)
K_STAR = 0.25


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    started = time.perf_counter()
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_eval_boot_syn")
    grid = load("85_phase18_synth_point_grid.py", "phase18_eval_boot_grid")
    shared = load("12_rsf_certification.py", "phase18_eval_boot_shared")
    remedy = load("11_rsf_remedy_calibration.py", "phase18_eval_boot_remedy")
    checkpoint = OUT / "PHASE18_ADIM4_KSTAR_EVAL_BOOT_CHECKPOINT.csv"
    rows = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {int(row["seed"]) for row in rows}
    for seed in SEEDS:
        if seed in done:
            continue
        phi = grid.simulate(seed, syn)
        magnitude = np.abs(phi[0]).mean(axis=(0, 1))
        scale = float(np.median(magnitude[magnitude > 0]))
        indices = np.random.default_rng(seed + 199_000).integers(0, 50, size=(999, 50))
        points, pvalues, bh, by = syn.stats(
            phi, indices, shared, remedy,
            K_STAR * scale, K_STAR * scale, 6 * K_STAR * scale, scale,
            add_n=True, sm=scale,
        )
        for feature in range(30):
            rows.append(
                {
                    "dataset": "EVAL_boot_feature", "status": "held_out_B999",
                    "seed": seed, "model": np.nan, "feature": feature,
                    "truth_class": CLASSES[feature], "k": K_STAR, "kappa": 0.25,
                    "k_star": K_STAR, "WI6": points[feature].min(),
                    "point_adequate": points[feature].min() >= 0.8,
                    "class_point_rate": np.nan, "v1_BH": np.nan,
                    "p_joint": pvalues[feature].max(),
                    "certified_BH": bh[feature], "certified_BY": by[feature],
                    "inner_draws": 999, "inner_seed": seed + 199_000,
                }
            )
        pd.DataFrame(rows).to_csv(checkpoint, index=False)
        done.add(seed)
        print(f"EVAL bootstrap {len(done)}/20", flush=True)
    base = pd.read_csv(OUT / "PHASE18_ADIM4_KSTAR.csv")
    base = base[base.dataset != "EVAL_boot_feature"]
    result = pd.concat([base, pd.DataFrame(rows)], ignore_index=True, sort=False)
    result.to_csv(OUT / "PHASE18_ADIM4_KSTAR.csv", index=False)
    boot = pd.DataFrame(rows)
    print(boot.groupby("truth_class")[["certified_BH", "certified_BY"]].mean())
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
