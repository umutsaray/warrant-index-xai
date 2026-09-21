"""Compare claim-level WI and baselines on identical synthetic fixtures (C.3)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pandas as pd

from warrantxai.aggregation import geometric_mean, veto_gated_sum, warrant_index
from warrantxai.baselines import eri_linear_stress, eri_mean, eri_min, single_component
from warrantxai.components import component_scores
from warrantxai.scenarios import SCENARIOS, mechanism_shift
from _shared import CALIBRATION, markdown_table, save_run


def run(seeds: int = 20, n: int = 1000, output: Path = Path("results")) -> pd.DataFrame:
    """Compare scores; do not presume WI's superiority (Section C.3).

    Parameters
    ----------
    seeds, n : int
        Independent repeats and anchor count.
    output : Path
        Artifact directory.

    Returns
    -------
    DataFrame
        Whole-claim scores; J* is ALL five features, fixed in advance.
    """
    if seeds < 1:
        raise ValueError("seeds must be positive")
    started = perf_counter()
    rows = []
    generators = {**SCENARIOS, "mechanism_sign_flip": lambda **kw: mechanism_shift(mode="sign_flip", **kw)}
    for name, generator in generators.items():
        for seed in range(seeds):
            s = generator(n=n, seed=20260916 + seed)
            z = component_scores(
                s.phi,
                s.reference,
                s.lower,
                s.upper,
                **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta"]},
                common_support=s.common_support,
            )
            eri = eri_linear_stress(s, scale=CALIBRATION["a"], seed=991 + seed)
            axes = np.array(list(eri.values()))
            scores = {
                "WI": float(np.min(warrant_index(z))),
                "geometric": float(np.min(geometric_mean(z))),
                "veto_sum": float(np.min(veto_gated_sum(z, 0.8))),
                "arithmetic": float(np.min(z.mean(axis=1))),
                "ERI-min": float(eri_min(axes)),
                "ERI-mean": float(eri_mean(axes)),
                **{f"{c}-only": float(np.min(single_component(z, c))) for c in ["S", "F", "G"]},
            }
            for method, score in scores.items():
                rows.append(
                    {
                        "scenario": name,
                        "seed": seed,
                        "method": method,
                        "claim_score": score,
                        "assessed": bool(np.isfinite(score)),
                        "passes_illustrative_0_8": float(score > 0.8) if np.isfinite(score) else np.nan,
                        **eri,
                    }
                )
    frame = pd.DataFrame(rows)
    summary = frame.pivot_table(
        index="scenario", columns="method", values="claim_score", aggfunc="mean", dropna=False
    ).reset_index()
    # Paired contrasts and uncertainty concern score differences, not superiority.
    paired = frame.pivot(index=["scenario", "seed"], columns="method", values="claim_score")
    paired["WI_minus_ERI_min"] = paired.WI - paired["ERI-min"]
    contrasts = paired.groupby(level=0).WI_minus_ERI_min.agg(["mean", "std", "count"])
    contrasts["approx_95_halfwidth"] = 1.96 * contrasts["std"] / np.sqrt(contrasts["count"])
    contrasts = contrasts.reset_index()
    report = (
        "# Deney 03 — baseline karşılaştırması\n\n"
        "J* bütün beş feature olarak önceden sabit. WI/ablation'lar feature minimumu; ERI whole-explanation skoru. "
        "Bu nedenle skor büyüklükleri aynı semantik ölçek değildir. Ortak 0.8 çizgisi yalnız illustrative; "
        "eşit false-positive-rate calibration yapılmadan superiority iddiası kurulamaz.\n\n"
        + markdown_table(summary)
        + "\n\n## ERI uygulama kapsamı\n\n"
        "Kaynak: [ERI v1 Definitions 5–9 / Appendix G](https://arxiv.org/html/2602.05082v1). "
        "Min/mean aggregation birebir; S/R/M/D/T kontrollü linear stress transforms ile hesaplanır. "
        "assumption: distance RMS((a-b)/0.5), Gaussian perturbation sd=0.02, 32 probes, "
        "redundancy alpha uniform [.9,1), T synthetic smooth trajectory, D frozen source model on native cohorts. "
        "Bu ERI-Bench replikasyonu değildir. S/R/M/D/T yerine WI bileşenleri kullanılmamıştır.\n\n"
        "## Yorum\n\n"
        "- Stabil sign-reversed yanlış explainer, ERI distance'larını oracle ile aynı bırakır; WI F bunu ayırır.\n"
        "- Native covariate shift ERI-D'yi düşürebilir; standardized WI G yüksek kalır. Bu farklı hedeflerdir.\n"
        "- Gerçek mechanism sign flip, WI magnitude-only G tarafından kaçırılır; genel üstünlük gösterilmemiştir.\n"
        "- Zero model için her iki yaklaşım yüksek puan verebilir; informative/causal explanation iddiası değildir.\n"
        "- OOD için WI NA korunur; ERI native-distribution stress sayıları certified external transfer değildir.\n\n"
        "## Paired score farkı (üstünlük testi değil)\n\n" + markdown_table(contrasts) + "\n"
    )
    save_run(output, "03_baselines", frame, report, started, {"seeds": seeds, "n": n, "claim_set": [0, 1, 2, 3, 4]})
    summary.to_csv(output / "03_baseline_summary.csv", index=False)
    print(summary.to_string(index=False))
    return frame


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, default=20)
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    run(args.seeds, args.n, args.output)
