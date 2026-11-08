#!/usr/bin/env python3
"""Phase 19 Step 4: component AUCs and diagnostic relative-G tolerance."""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "phase19"
OUT.mkdir(parents=True, exist_ok=True)
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from phase19_common import CLASSES, load, pooled_auc, quantile_signed, scale_and_n, simulate

P19_SEEDS = tuple(range(19_000_001, 19_000_021))
CAL_SEEDS = tuple(range(18_000_001, 18_000_061))
CS = (0.25, 0.5, 1.0)
COMPONENTS = ("F", "S", "U", "G_claim", "R")
FORMULA_FILES = ("components.py", "aggregation.py", "certification.py")


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest().upper()


def checkpoint_rows(path: Path) -> list[dict[str, object]]:
    return pd.read_csv(path).to_dict("records") if path.exists() else []


def write_checkpoint(rows: list[dict[str, object]], path: Path) -> None:
    pd.DataFrame(rows).to_csv(path, index=False)


def auc_rows(frame: pd.DataFrame, dataset: str) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    metrics = ("F", "S", "U", "G_claim", "R", "N", "M_j", "WI6")
    for metric in metrics:
        for positive_class in "ACDE":
            positive = frame.loc[frame.truth_class == positive_class, metric].to_numpy(float)
            negative = frame.loc[frame.truth_class == "B", metric].to_numpy(float)
            rows.append({
                "dataset": dataset,
                "metric": metric,
                "positive_class": positive_class,
                "negative_class": "B",
                "auc": pooled_auc(positive, negative),
                "n_positive": len(positive),
                "n_negative": len(negative),
            })
        a = frame.loc[frame.truth_class == "A", metric].to_numpy(float)
        d = frame.loc[frame.truth_class == "D", metric].to_numpy(float)
        rows.append({
            "dataset": dataset,
            "metric": metric,
            "positive_class": "A",
            "negative_class": "D",
            "auc": pooled_auc(a, d),
            "n_positive": len(a),
            "n_negative": len(d),
        })
    return rows


def main() -> None:
    started = time.perf_counter()
    breakdown = load("84_phase18_synth_breakdown.py", "phase19_step4_breakdown")
    shared = load("12_rsf_certification.py", "phase19_step4_shared")
    pilot = load("08_certification_pilot.py", "phase19_step4_pilot")

    p19_checkpoint = OUT / "PHASE19_ADIM4_P19_CHECKPOINT.csv"
    p19_rows = checkpoint_rows(p19_checkpoint)
    done_p19 = {int(row["seed"]) for row in p19_rows}
    rel_rows: list[dict[str, object]] = []
    if (OUT / "PHASE19_ADIM4_REL_CHECKPOINT.csv").exists():
        rel_rows = checkpoint_rows(OUT / "PHASE19_ADIM4_REL_CHECKPOINT.csv")
    done_rel = {int(row["seed"]) for row in rel_rows}

    for seed in P19_SEEDS:
        if seed in done_p19 and seed in done_rel:
            continue
        sim = simulate(seed)
        phi = np.asarray(sim["phi"])
        scale, magnitude, n_point = scale_and_n(phi)
        values5 = breakdown.point_components(phi, scale, scale, 6 * scale, scale, shared, pilot)
        if seed not in done_p19:
            for feature in range(30):
                values = values5[feature]
                p19_rows.append({
                    "seed": seed,
                    "feature": feature,
                    "truth_class": CLASSES[feature],
                    "s_m_sim": scale,
                    "M_j": magnitude[feature],
                    **{name: values[q] for q, name in enumerate(COMPONENTS)},
                    "N": n_point[feature],
                    "WI6": min(values.min(), n_point[feature]),
                })
            write_checkpoint(p19_rows, p19_checkpoint)
        if seed not in done_rel:
            mbar = np.abs(phi).mean(axis=(1, 2)).mean(axis=0)
            absolute_g = values5[:, 3]
            for c in CS:
                g_rel = quantile_signed(phi, c * np.maximum(mbar, np.finfo(float).tiny))
                for feature in range(30):
                    rel_rows.append({
                        "seed": seed,
                        "feature": feature,
                        "truth_class": CLASSES[feature],
                        "c": c,
                        "Mbar_j": mbar[feature],
                        "G_absolute_claim": absolute_g[feature],
                        "G_rel": g_rel[feature],
                        "absolute_pass": bool(absolute_g[feature] >= 0.8),
                        "relative_pass": bool(g_rel[feature] >= 0.8),
                    })
            write_checkpoint(rel_rows, OUT / "PHASE19_ADIM4_REL_CHECKPOINT.csv")
        print(f"P19 step4 seed {seed} complete", flush=True)

    p19 = pd.DataFrame(p19_rows).sort_values(["seed", "feature"]).reset_index(drop=True)
    assert len(p19) == 20 * 30 == 600
    assert p19.seed.nunique() == 20

    # Reuse frozen CAL components; simulate only M_j, which was not stored in the grid.
    cal_m_checkpoint = OUT / "PHASE19_ADIM4_CAL_M_CHECKPOINT.csv"
    cal_m_rows = checkpoint_rows(cal_m_checkpoint)
    done_cal = {int(row["seed"]) for row in cal_m_rows}
    for seed in CAL_SEEDS:
        if seed in done_cal:
            continue
        phi = np.asarray(simulate(seed)["phi"])
        _, magnitude, _ = scale_and_n(phi)
        cal_m_rows.extend({
            "seed": seed,
            "feature": feature,
            "truth_class": CLASSES[feature],
            "M_j": magnitude[feature],
        } for feature in range(30))
        write_checkpoint(cal_m_rows, cal_m_checkpoint)
        print(f"CAL M seed {seed} complete", flush=True)

    cal_m = pd.DataFrame(cal_m_rows)
    assert len(cal_m) == 60 * 30 == 1800
    grid = pd.read_csv(ROOT / "results/phase18/PHASE18_ADIM2_POINT_GRID.csv")
    cal = grid[
        (grid.variant == "v2_WI6")
        & (grid.k == 1.0)
        & (grid.rho == 1.0)
        & (grid.kappa == 0.25)
    ].copy()
    assert len(cal) == 60 * 30 == 1800
    cal = cal.merge(cal_m[["seed", "feature", "M_j"]], on=["seed", "feature"], validate="one_to_one")
    cal = cal.rename(columns={"WI": "WI6"})

    legacy = pd.read_csv(ROOT / "results/phase18/PHASE18_ADIM1_COMPONENT_BREAKDOWN.csv").rename(
        columns={
            "repeat": "seed",
            "v2_F": "F",
            "v2_S": "S",
            "v2_U": "U",
            "v2_G_claim": "G_claim",
            "v2_R": "R",
            "v2_N": "N",
        }
    )
    legacy["WI6"] = legacy["v2_WI6"]

    auc = pd.DataFrame(
        auc_rows(legacy, "PHASE18_original20")
        + auc_rows(cal, "CAL60")
        + auc_rows(p19, "P19_first20")
    )
    auc.to_csv(OUT / "PHASE19_ADIM4_AUC_COMPONENTS.csv", index=False)

    rel = pd.DataFrame(rel_rows).sort_values(["seed", "feature", "c"]).reset_index(drop=True)
    assert len(rel) == 20 * 30 * 3 == 1800
    rel.to_csv(OUT / "PHASE19_ADIM4_REL_G.csv", index=False)
    rel_summary = []
    for c in CS:
        cell = rel[rel.c == c]
        for truth_class in "ABCDE":
            values = cell.loc[cell.truth_class == truth_class, "G_rel"].to_numpy(float)
            b = cell.loc[cell.truth_class == "B", "G_rel"].to_numpy(float)
            rel_summary.append({
                "c": c,
                "truth_class": truth_class,
                "median_G_rel": float(np.median(values)),
                "q1_G_rel": float(np.quantile(values, 0.25)),
                "q3_G_rel": float(np.quantile(values, 0.75)),
                "pass_rate_G_rel": float(np.mean(values >= 0.8)),
                "auc_vs_B": np.nan if truth_class == "B" else pooled_auc(values, b),
            })
    rel_summary_frame = pd.DataFrame(rel_summary)
    rel_summary_frame.to_csv(OUT / "PHASE19_ADIM4_REL_G_SUMMARY.csv", index=False)

    legacy_expected = {
        ("G_claim", "A"): 0.235,
        ("G_claim", "C"): 0.335,
        ("G_claim", "D"): 0.230,
        ("G_claim", "E"): 0.269,
        ("N", "A"): 1.000,
        ("N", "D"): 1.000,
        ("N", "E"): 1.000,
        ("N", "C"): 0.806,
        ("M_j", "A"): 1.000,
        ("M_j", "D"): 1.000,
        ("M_j", "E"): 1.000,
        ("M_j", "C"): 0.805,
    }
    # The table also contains the preregistered A-vs-D contrast.  Validate the
    # frozen Phase 18 values on the class-vs-B rows only.
    legacy_auc = auc[
        (auc.dataset == "PHASE18_original20") & (auc.negative_class == "B")
    ].set_index(["metric", "positive_class"]).auc
    for key, expected in legacy_expected.items():
        assert np.isclose(legacy_auc.loc[key], expected, atol=0.002), (key, legacy_auc.loc[key], expected)

    a_abs_rate = float(p19.loc[p19.truth_class == "A", "G_claim"].ge(0.8).mean())
    h5 = {}
    for c in CS:
        a_cell = rel[(rel.c == c) & (rel.truth_class == "A")]
        class_auc = rel_summary_frame[(rel_summary_frame.c == c) & (rel_summary_frame.truth_class != "B")]
        h5[str(c)] = {
            "all_class_auc_ge_0_5": bool(class_auc.auc_vs_B.ge(0.5).all()),
            "A_relative_pass_rate": float(a_cell.relative_pass.mean()),
            "A_absolute_pass_rate": a_abs_rate,
            "A_pass_rate_higher": bool(a_cell.relative_pass.mean() > a_abs_rate),
        }

    metadata = {
        "p19_rows": len(p19),
        "cal_rows": len(cal),
        "relative_rows": len(rel),
        "expected": {"p19_rows": 600, "cal_rows": 1800, "relative_rows": 1800},
        "h5_by_c": h5,
        "completed_p19_seeds": sorted(p19.seed.unique().astype(int).tolist()),
        "completed_cal_seeds": sorted(cal_m.seed.unique().astype(int).tolist()),
        "formula_sha256": {name: sha256(ROOT / "warrantxai" / name) for name in FORMULA_FILES},
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "PHASE19_ADIM4_METADATA.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(auc.to_string(index=False))
    print(rel_summary_frame.to_string(index=False))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
