"""Run the fixed n/rho oracle grid (Section C.3), from the project root.

Example: python experiments/01_calibration_smallp.py --seeds 5
assumption: tolerances are preregistered illustrative constants, not fitted on
this sweep. 'Calibration' here checks numerical behaviour, not threshold tuning.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from warrantxai.components import component_scores, importance_ci
from warrantxai.dgp import correlated_sem, exact_shapley
from _shared import CALIBRATION, markdown_table, save_run, score_rows


def run(seeds: int = 5, output: Path = Path("results")) -> pd.DataFrame:
    """Run 3 anchor sizes x 3 correlations x independent seeds (Section C.3).

    Parameters
    ----------
    seeds : int
        Number of independent DGP draws per cell.
    output : Path
        Artifact directory.

    Returns
    -------
    DataFrame
        Feature-level scores for oracle and two controlled corruptions.
    """
    if seeds < 1:
        raise ValueError("seeds must be positive")
    started = perf_counter()
    rows = []
    max_oracle_error = 0.0
    beta = np.array([2.0, -1.5, 1.0, -0.5, 0.25])
    for n in [250, 1000, 4000]:
        for rho in [0.0, 0.5, 0.9]:
            dgp = correlated_sem(beta, rho)
            for seed in range(seeds):
                rng = np.random.default_rng(100000 + seed + int(rho * 1000) + n)
                x, _ = dgp.sample(n, rng)
                oracle = dgp.oracle_attributions(x)
                background = rng.normal(size=(16, len(beta)))
                background = np.concatenate([background, -background])
                exact = exact_shapley(lambda a: a @ beta, x[:8], background)
                max_oracle_error = max(max_oracle_error, float(np.max(np.abs(exact - oracle[:8]))))
                reference = np.broadcast_to(oracle, (2, 4, *oracle.shape)).copy()
                # Two replicas/sites share the exact mechanism: replication is
                # stipulated in this oracle unit experiment, not inferred from data.
                for name in ["oracle_invariant", "stable_but_wrong", "replicate_corruption"]:
                    phi = reference.copy()
                    if name == "stable_but_wrong":
                        phi = -phi  # stable sign-reversed explanation, unchanged model
                    if name == "replicate_corruption":
                        phi[:, 1::2] = -phi[:, 1::2]
                    lo, hi = importance_ci(phi, alpha=CALIBRATION["alpha_U"])
                    z = component_scores(phi, reference, lo, hi, **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta"]})
                    rows.extend(score_rows(z, scenario=name, n=n, rho=rho, seed=seed))
    frame = pd.DataFrame(rows)
    summary = frame.groupby(["scenario", "n"])[["F", "S", "U", "G", "R", "WI", "point_adequate"]].mean().reset_index()
    report = (
        "# Deney 01 — small-p oracle taraması\n\n"
        f"Grid: n={{250,1000,4000}}, rho={{0,0.5,0.9}}, {seeds} bağımsız seed; p=5, sites=2, T=4.\n\n"
        f"Exact empirical-coalition vs closed-form maximum error: {max_oracle_error:.3e}.\n\n"
        "assumption: marginal linear SHAP, antithetic zero-mean background; normal-CLT conditional CIs; "
        "tolerances a=b=0.5,u=3,eta=0.02,tau=0.8 bütün evaluation verilerinden önce sabit. "
        "Bu tarama threshold seçimi veya population certification değildir.\n\n" + markdown_table(summary) + "\n\n"
        "Stable sign reversal'da R ve G yüksek kalabilir: ikisi magnitude importance kullanır. "
        "Değişmeyen modelin yanlış açıklamasını F yakalamalıdır. Replicate corruption, trained-model "
        "instability değil, S formülüne yönelik kontrollü attribution stress testidir.\n"
    )
    save_run(
        output,
        "01_calibration",
        frame,
        report,
        started,
        {"seeds": seeds, "grid_n": [250, 1000, 4000], "grid_rho": [0, 0.5, 0.9]},
    )
    print(summary.to_string(index=False))
    print(f"max exact-oracle error={max_oracle_error:.3e}")
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.seeds, args.output)
