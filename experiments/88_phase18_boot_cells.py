#!/usr/bin/env python3
"""Phase 18 limited B=999 certification cells for the first 20 CAL seeds."""
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
SEEDS = tuple(range(18_000_001, 18_000_021))
KS = (1.0, 2.0, 4.0)
CLASSES = np.array(["A"] * 8 + ["B"] * 6 + ["C"] * 6 + ["D"] * 5 + ["E"] * 5)
COMPONENTS = ("F", "S", "U", "G_claim", "R", "N")


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def main() -> None:
    started = time.perf_counter()
    syn = load("80_phase17_synthetic_feature_level.py", "phase18_boot_syn")
    grid = load("85_phase18_synth_point_grid.py", "phase18_boot_grid")
    shared = load("12_rsf_certification.py", "phase18_boot_shared")
    remedy = load("11_rsf_remedy_calibration.py", "phase18_boot_remedy")
    checkpoint = OUT / "PHASE18_ADIM2_BOOT_CELLS_CHECKPOINT.csv"
    rows = pd.read_csv(checkpoint).to_dict("records") if checkpoint.exists() else []
    done = {(int(row["seed"]), float(row["k"])) for row in rows}
    for seed in SEEDS:
        phi = grid.simulate(seed, syn)
        magnitude = np.abs(phi[0]).mean(axis=(0, 1))
        scale = float(np.median(magnitude[magnitude > 0]))
        indices = np.random.default_rng(seed + 99_000).integers(0, 50, size=(999, 50))
        for k in KS:
            if (seed, k) in done:
                continue
            cell_started = time.perf_counter()
            points, pvalues, bh, by = syn.stats(
                phi, indices, shared, remedy,
                k * scale, k * scale, 6 * k * scale, scale,
                add_n=True, sm=scale,
            )
            for feature in range(30):
                rows.append(
                    {
                        "seed": seed, "feature": feature, "truth_class": CLASSES[feature],
                        "k": k, "rho": 1.0, "kappa": 0.25, "s_m_sim": scale,
                        **{name: points[feature, q] for q, name in enumerate(COMPONENTS)},
                        **{f"p_{name}": pvalues[feature, q] for q, name in enumerate(COMPONENTS)},
                        "WI6": points[feature].min(), "p_joint": pvalues[feature].max(),
                        "certified_BH": bh[feature], "certified_BY": by[feature],
                        "cell_seconds": time.perf_counter() - cell_started,
                        "inner_draws": 999, "inner_seed": seed + 99_000,
                    }
                )
            done.add((seed, k))
            pd.DataFrame(rows).to_csv(checkpoint, index=False)
            print(f"boot cell seed={seed} k={k:g} done={len(done)}/60", flush=True)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "PHASE18_ADIM2_BOOT_CELLS.csv", index=False)
    print(result.groupby(["k", "truth_class"])[["certified_BH", "certified_BY"]].mean().to_string())
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
