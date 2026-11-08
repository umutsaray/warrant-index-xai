#!/usr/bin/env python3
"""Phase 18 threshold-free class-vs-null discrimination AUCs."""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase18"
BOOTSTRAP_DRAWS = 9_999
BOOTSTRAP_SEED = 18_180_001


def pair_auc_matrix(cell: pd.DataFrame, target: str, seeds: list[int]) -> tuple[np.ndarray, float]:
    positives = {seed: cell[(cell.seed == seed) & (cell.truth_class == target)].WI.to_numpy() for seed in seeds}
    negatives = {seed: cell[(cell.seed == seed) & (cell.truth_class == "B")].WI.to_numpy() for seed in seeds}
    matrix = np.empty((len(seeds), len(seeds)))
    for i, left in enumerate(seeds):
        for j, right in enumerate(seeds):
            difference = positives[left][:, None] - negatives[right][None, :]
            matrix[i, j] = (np.count_nonzero(difference > 0) + 0.5 * np.count_nonzero(difference == 0))
    n_positive = len(positives[seeds[0]])
    n_negative = len(negatives[seeds[0]])
    point = matrix.sum() / (n_positive * n_negative * len(seeds) ** 2)
    return matrix, point


def main() -> None:
    started = time.perf_counter()
    frame = pd.read_csv(OUT / "PHASE18_ADIM2_POINT_GRID.csv")
    seeds = sorted(frame.seed.unique().astype(int).tolist())
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    sampled = rng.integers(0, len(seeds), size=(BOOTSTRAP_DRAWS, len(seeds)))
    counts = np.zeros((BOOTSTRAP_DRAWS, len(seeds)), dtype=np.int16)
    np.add.at(counts, (np.arange(BOOTSTRAP_DRAWS)[:, None], sampled), 1)
    rows = []
    settings = frame[["variant", "k", "rho", "kappa"]].drop_duplicates()
    for setting in settings.itertuples(index=False):
        mask = frame.variant.eq(setting.variant)
        for column in ("k", "rho", "kappa"):
            value = getattr(setting, column)
            mask &= frame[column].isna() if pd.isna(value) else frame[column].eq(value)
        cell = frame[mask]
        for target in "ACDE":
            matrix, point = pair_auc_matrix(cell, target, seeds)
            n_positive = int((cell.truth_class == target).groupby(cell.seed).sum().iloc[0])
            n_negative = int((cell.truth_class == "B").groupby(cell.seed).sum().iloc[0])
            numerator = np.einsum("bi,ij,bj->b", counts, matrix, counts, optimize=True)
            bootstrap = numerator / (n_positive * n_negative * len(seeds) ** 2)
            rows.append(
                {
                    "variant": setting.variant,
                    "k": setting.k,
                    "rho": setting.rho,
                    "kappa": setting.kappa,
                    "positive_class": target,
                    "negative_class": "B",
                    "auc_feature_repeat": point,
                    "bootstrap95_lower": np.quantile(bootstrap, 0.025),
                    "bootstrap95_upper": np.quantile(bootstrap, 0.975),
                    "bootstrap_draws": BOOTSTRAP_DRAWS,
                    "bootstrap_seed": BOOTSTRAP_SEED,
                    "repeat_clusters": len(seeds),
                }
            )
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "PHASE18_ADIM3_AUC.csv", index=False)
    print(result.to_string(index=False))
    print("elapsed", time.perf_counter() - started)


if __name__ == "__main__":
    main()
