"""Benchmark exact vs incomplete S/R on increasing p (Section E)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter
import tracemalloc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from warrantxai.components import importance, rank_reproducibility, stability
from warrantxai.complexity import incomplete_rank, incomplete_stability
from _shared import markdown_table, save_run


def run(repeats: int = 5, output: Path = Path("results")) -> pd.DataFrame:
    """Measure wall time and approximation error, not assumed speedup (Section E).

    Parameters
    ----------
    repeats : int
        Repeated timing and independent subsampling runs.
    output : Path
        Results directory.

    Returns
    -------
    DataFrame
        Timing/error/budget records for each feature dimension.
    """
    if repeats < 1:
        raise ValueError("repeats must be positive")
    started = perf_counter()
    rows = []
    for p in [10, 50, 200, 500]:
        rng = np.random.default_rng(902 + p)
        n, t, sites = 128, 12, 2
        x = rng.normal(size=(n, p))
        beta = rng.normal(size=(sites, t, p))
        phi = x[None, None, :, :] * beta[:, :, None, :]
        imp = importance(phi)
        for repeat in range(repeats):
            for component, exact_fn, approx_fn in [
                ("S", lambda: stability(phi, 2), lambda: incomplete_stability(phi, 2, 12, repeat)),
                (
                    "R",
                    lambda: rank_reproducibility(imp, 0.02),
                    lambda: incomplete_rank(imp, 0.02, 12, min(16, p - 1), repeat),
                ),
            ]:
                start = perf_counter()
                exact = exact_fn()
                exact_s = perf_counter() - start
                start = perf_counter()
                approx = approx_fn()
                approx_s = perf_counter() - start
                # Separate memory run, so tracing overhead does not bias timing.
                tracemalloc.start()
                exact_fn()
                _, peak_exact = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                tracemalloc.start()
                approx_fn()
                _, peak_approx = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                rows.append(
                    {
                        "p": p,
                        "component": component,
                        "repeat": repeat,
                        "exact_seconds": exact_s,
                        "approx_seconds": approx_s,
                        "speedup": exact_s / approx_s,
                        "mean_abs_error": float(np.abs(exact - approx).mean()),
                        "max_abs_error": float(np.abs(exact - approx).max()),
                        "threshold_0_8_disagreement": float(np.mean((exact > 0.8) != (approx > 0.8))),
                        "exact_peak_traced_bytes": peak_exact,
                        "approx_peak_traced_bytes": peak_approx,
                        "cache_bytes": phi.nbytes,
                        "all_training_pairs": t * (t - 1) // 2,
                        "sampled_training_pairs": 12,
                        "comparisons_per_feature": min(16, p - 1) if component == "R" else p - 1,
                    }
                )
    frame = pd.DataFrame(rows)
    summary = (
        frame.groupby(["p", "component"])[
            [
                "exact_seconds",
                "approx_seconds",
                "speedup",
                "mean_abs_error",
                "max_abs_error",
                "threshold_0_8_disagreement",
            ]
        ]
        .median()
        .reset_index()
    )
    report = (
        "# Deney 04 — incomplete U-statistics maliyet/doğruluk\n\n"
        "sites=2,T=12,n=128; exact 66 training pair, approximation 12 pair; R için feature başına en fazla 16 comparison. "
        f"{repeats} timing/subsampling tekrarı; tabloda median. Model eğitimi ve explainer maliyeti dahil değil.\n\n"
        + markdown_table(summary)
        + "\n\n"
        "assumption: cached fixed tensor üzerinde random pair seçimi; site mean unbiased olabilir fakat worst-site min unbiased değildir. "
        "Threshold disagreement 0.8 çizgisinde ölçülür, sıfır çıkması başka eşiklerde garanti vermez. "
        "tracemalloc peak Python/NumPy traced allocations ölçer; process RSS veya GPU VRAM değildir. "
        "Timing tek makineye özgüdür; small-p sampling overhead hız kazancını ortadan kaldırabilir.\n"
    )
    save_run(output, "04_complexity", frame, report, started, {"repeats": repeats, "p_grid": [10, 50, 200, 500]})
    print(summary.to_string(index=False))
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.repeats, args.output)
