"""Validate IUT/BH/BY on known-null synthetic p-values (Section 2.4).

This simulation tests the multiple-testing machinery, NOT validity of bootstrap
p-values estimated from a single WI dataset. Those are deliberately not invented.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd
from scipy.stats import beta

from warrantxai.certification import certify, intersection_union
from _shared import markdown_table, save_run


def run(studies: int = 10000, output: Path = Path("results")) -> pd.DataFrame:
    """Measure Type I and FDR under explicit known nulls (Section 2.4).

    Parameters
    ----------
    studies : int
        Independent family-level simulation runs per dependence condition.
    output : Path
        Artifact directory.

    Returns
    -------
    DataFrame
        Empirical family FDP/power, with true null status known by generation.
    """
    if studies < 1:
        raise ValueError("studies must be positive")
    started = perf_counter()
    rng = np.random.default_rng(20260916)
    n = 200000
    p = np.zeros((n, 5))
    p[:, 0] = rng.uniform(size=n)
    failures = int(np.count_nonzero(intersection_union(p) <= 0.05))
    interval = (float(beta.ppf(0.025, failures, n - failures + 1)), float(beta.ppf(0.975, failures + 1, n - failures)))
    rows = []
    for dependence in ["independent", "shared_rotated_uniform"]:
        for study in range(studies):
            # First 10 claims are null because ONE required criterion is null;
            # other 10 claims are strong alternatives. All null marginals valid.
            pvals = np.zeros((20, 5))
            nulls = rng.uniform(size=10) if dependence == "independent" else (rng.uniform() + np.arange(10) / 10) % 1
            pvals[np.arange(10), np.arange(10) % 5] = nulls
            pvals[10:] = rng.uniform(0, 0.001, size=(10, 5))
            methods = ["bh", "by"] if dependence == "independent" else ["by"]
            for method in methods:
                result = certify(pvals, alpha=0.05, method=method)
                false = int(result.rejected[:10].sum())
                total = int(result.rejected.sum())
                rows.append(
                    {
                        "dependence": dependence,
                        "method": method,
                        "study": study,
                        "FDP": false / max(total, 1),
                        "power": result.rejected[10:].mean(),
                    }
                )
    frame = pd.DataFrame(rows)
    summary = (
        frame.groupby(["dependence", "method"])
        .agg(empirical_FDR=("FDP", "mean"), FDP_std=("FDP", "std"), power=("power", "mean"), runs=("FDP", "size"))
        .reset_index()
    )
    summary["FDR_MC_SE"] = summary.FDP_std / np.sqrt(summary.runs)
    report = (
        "# Deney 05 — certification kontrolü\n\n"
        f"IUT least-favourable union null: {failures}/{n} false positives, rate={failures / n:.6f}; "
        f"exact binomial 95% CI=[{interval[0]:.6f},{interval[1]:.6f}], alpha=0.05.\n\n"
        + markdown_table(summary)
        + "\n\n"
        "assumption: null p-values exact Uniform(0,1); bir feature'ın tek gerekli bileşeni null, diğerleri alternative. "
        "Feature başına max-p kullanılır. Rotated common uniforms arbitrary dependence stresidir; BH için PRDS iddiası kurulmadığından yalnız BY çalıştırılır. "
        "Monte Carlo sonuçları teorik ispat değildir; binomial interval alpha'yı içerebilir. "
        "Bu deney bir gerçek WI dataset'inden üretilmiş p-value'nun geçerliliğini göstermez.\n"
    )
    save_run(
        output, "05_certification", frame, report, started, {"studies": studies, "iut_studies": n, "seed": 20260916}
    )
    print(f"IUT false-positive rate={failures / n:.6f}, binomial CI={interval}")
    print(summary.to_string(index=False))
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--studies", type=int, default=10000)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.studies, args.output)
