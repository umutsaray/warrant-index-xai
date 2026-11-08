"""Phase 19 step 1: diagnose signed-generalizability coupling."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from phase19_common import CLASSES, pooled_auc, quantile_signed, scale_and_n, simulate

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase19"
OUT.mkdir(parents=True, exist_ok=True)
SEEDS = tuple(range(19000001, 19000021))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def one_seed(seed: int, *, skew_d: bool = False) -> pd.DataFrame:
    sim = simulate(seed, skew_d=skew_d)
    phi = np.asarray(sim["phi"])
    coef = np.asarray(sim["coefficients"])
    x0 = np.asarray(sim["x0"])
    scale, magnitude0, _ = scale_and_n(phi)
    importance = np.abs(phi).mean(axis=(1, 2))
    g = 1.0 - np.minimum(1.0, np.abs(importance[1] - importance[0]) / scale)
    sorted_signed = quantile_signed(phi, scale)

    # Same source anchors for both sites; only fitted coefficients differ.
    centered = x0 - x0.mean(axis=0)
    paired_means = np.stack([
        coef[site].mean(axis=0) * centered for site in range(2)
    ])
    paired_signed = 1.0 - np.minimum(
        1.0, np.abs(paired_means[1] - paired_means[0]) / scale
    ).mean(axis=0)

    return pd.DataFrame({
        "seed": seed,
        "feature": np.arange(30),
        "truth_class": CLASSES,
        "scale": scale,
        "magnitude_source": magnitude0,
        "G": g,
        "G_signed_sorted": sorted_signed,
        "G_claim_sorted": np.minimum(g, sorted_signed),
        "G_signed_paired": paired_signed,
        "G_claim_paired": np.minimum(g, paired_signed),
        "skewed_D": skew_d,
    })


def summarize(frame: pd.DataFrame, label: str) -> pd.DataFrame:
    rows = []
    for cls in "AD":
        cell = frame[frame.truth_class == cls]
        for metric in ["G", "G_signed_sorted", "G_claim_sorted", "G_signed_paired", "G_claim_paired"]:
            values = cell[metric].to_numpy(float)
            rows.append({"dataset": label, "truth_class": cls, "metric": metric,
                         "median": np.median(values), "q1": np.quantile(values, .25),
                         "q3": np.quantile(values, .75), "n": len(values)})
    return pd.DataFrame(rows)


def main() -> None:
    started = time.perf_counter()
    cp = OUT / "PHASE19_ADIM1_CHECKPOINT.csv"
    done = pd.read_csv(cp) if cp.exists() else pd.DataFrame()
    completed = set(done.seed.unique()) if not done.empty else set()
    for seed in SEEDS:
        if seed not in completed:
            done = pd.concat([done, one_seed(seed)], ignore_index=True)
            done.to_csv(cp, index=False)
            print(f"P19 step1 seed {seed} complete", flush=True)
    assert len(done) == 20 * 30
    done.to_csv(OUT / "PHASE19_ADIM1_GSIGNED_ALIGNMENT.csv", index=False)

    skew_cp = OUT / "PHASE19_ADIM1_SKEW_CHECKPOINT.csv"
    skew = pd.read_csv(skew_cp) if skew_cp.exists() else pd.DataFrame()
    completed_skew = set(skew.seed.unique()) if not skew.empty else set()
    for seed in SEEDS:
        if seed not in completed_skew:
            skew = pd.concat([skew, one_seed(seed, skew_d=True)], ignore_index=True)
            skew.to_csv(skew_cp, index=False)
            print(f"P19 step1 skew seed {seed} complete", flush=True)
    assert len(skew) == 20 * 30
    skew.to_csv(OUT / "PHASE19_ADIM1_SKEWED_D.csv", index=False)

    auc_rows = []
    for label, frame in [("Gaussian", done), ("Skewed_D", skew)]:
        b = frame.loc[frame.truth_class == "B"]
        for metric in ["G_signed_sorted", "G_claim_sorted", "G_signed_paired", "G_claim_paired"]:
            for cls in "ACDE":
                pos = frame.loc[frame.truth_class == cls, metric].to_numpy(float)
                neg = b[metric].to_numpy(float)
                auc_rows.append({"dataset": label, "metric": metric, "positive_class": cls,
                                 "negative_class": "B", "auc": pooled_auc(pos, neg),
                                 "n_positive": len(pos), "n_negative": len(neg)})
            a = frame.loc[frame.truth_class == "A", metric].to_numpy(float)
            d = frame.loc[frame.truth_class == "D", metric].to_numpy(float)
            auc_rows.append({"dataset": label, "metric": metric, "positive_class": "A",
                             "negative_class": "D", "auc": pooled_auc(a, d),
                             "n_positive": len(a), "n_negative": len(d)})
    auc = pd.DataFrame(auc_rows)
    auc.to_csv(OUT / "PHASE19_ADIM1_PAIRED_AUC.csv", index=False)
    summary = pd.concat([summarize(done, "Gaussian"), summarize(skew, "Skewed_D")], ignore_index=True)
    summary.to_csv(OUT / "PHASE19_ADIM1_SUMMARY.csv", index=False)

    a_med = summary.query("dataset == 'Gaussian' and truth_class == 'A' and metric == 'G_signed_paired'")["median"].iloc[0]
    d_med = summary.query("dataset == 'Gaussian' and truth_class == 'D' and metric == 'G_signed_paired'")["median"].iloc[0]
    ad_auc = auc.query("dataset == 'Gaussian' and metric == 'G_signed_paired' and positive_class == 'A' and negative_class == 'D'").auc.iloc[0]
    meta = {
        "rows": {"alignment": len(done), "skewed": len(skew), "auc": len(auc)},
        "expected_rows": {"alignment": 600, "skewed": 600, "auc": 40},
        "H1": {"paired_A_vs_D_auc": ad_auc, "paired_A_median": a_med,
               "paired_D_median": d_med, "median_gap": a_med - d_med,
               "supported": bool(ad_auc >= .9 and a_med - d_med >= .3)},
        "implementation_alignment": {
            "real_pipeline": "experiments/12_rsf_certification.py:signed_scores sorts query axis after replicate averaging",
            "synthetic": "experiments/84_phase18_synth_breakdown.py:point_components calls the same signed_scores",
            "diagnostic": "paired variant uses identical source anchors without modifying the frozen function"
        },
        "formula_sha256": {name: sha(ROOT / "warrantxai" / name) for name in ["components.py", "aggregation.py", "certification.py"]},
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "PHASE19_ADIM1_METADATA.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(summary.to_string(index=False))
    print(auc.to_string(index=False))
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
