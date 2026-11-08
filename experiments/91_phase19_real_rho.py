#!/usr/bin/env python3
"""Phase 19 Step 3: frozen-tensor rho sweep and RSF N diagnostics."""
from __future__ import annotations

import hashlib
import importlib.util
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

SCALES = {
    "CoxPH-v2": 0.023934215222165295,
    "RSF": 0.01707069437342852,
    "DeepSurv": 0.012307440008774696,
}
RHOS = (1.0, 2.0, 4.0, 8.0, 16.0)
K = 64.0
KAPPA = 0.25
COMPONENTS = ("F", "S", "U", "G_claim", "R", "N")
FORMULA_FILES = ("components.py", "aggregation.py", "certification.py")


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest().upper()


def site_phi(item: dict) -> list[np.ndarray]:
    if "external" in item:
        return [item["phi"][0][:, :, item["targets"]], item["external"]]
    return [item["phi"][0], item["phi"][1]]


def n_details(item: dict, m0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sites = site_phi(item)
    ihat = np.stack([np.abs(site).mean(axis=(0, 1)) for site in sites])
    leaves = np.stack([
        np.minimum(1.0, np.abs(site).mean(axis=1) / m0).mean(axis=0)
        for site in sites
    ])
    return ihat, leaves, leaves.min(axis=0)


def main() -> None:
    started = time.perf_counter()
    phase12 = load("60_phase12_scale_tie_query_audit.py", "phase19_step3_phase12")
    data = phase12.load_inputs()
    rows: list[dict[str, object]] = []
    n_cache: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}

    for model, item in data.items():
        scale = SCALES[model]
        n_cache[model] = n_details(item, KAPPA * scale)
        n_point = n_cache[model][2]
        for rho in RHOS:
            base = phase12.scores(
                item,
                K * scale,
                K * scale,
                6 * K * scale,
                rho * scale,
            )
            for j, feature in enumerate(item["features"]):
                values = np.r_[
                    base.loc[feature, ["F", "S", "U", "G_claim", "R"]].to_numpy(float),
                    n_point[j],
                ]
                rows.append({
                    "model": model,
                    "feature": feature,
                    "k": K,
                    "rho": rho,
                    "kappa": KAPPA,
                    "a": K * scale,
                    "eta": rho * scale,
                    **{name: values[q] for q, name in enumerate(COMPONENTS)},
                    "WI6": float(values.min()),
                    "point_adequate6": bool(values.min() >= 0.8),
                    "bottleneck6": COMPONENTS[int(values.argmin())],
                })

    output = pd.DataFrame(rows)
    assert len(output) == 3 * len(RHOS) * 27 == 405
    output.to_csv(OUT / "PHASE19_ADIM3_RHO_SWEEP.csv", index=False)

    # Frozen Phase 18 reproduction checks.
    base = output[output.rho == 1.0]
    expected_binding = {
        "CoxPH-v2": {"R": 19, "N": 8},
        "RSF": {"N": 17, "R": 10},
        "DeepSurv": {"R": 17, "N": 10},
    }
    for model, expected in expected_binding.items():
        observed = base.loc[base.model == model, "bottleneck6"].value_counts().to_dict()
        assert observed == expected, (model, observed, expected)
    deep_r = base.loc[base.model == "DeepSurv", "R"]
    assert np.isclose(deep_r.min(), 0.625, atol=5e-4)
    assert np.isclose(deep_r.max(), 0.783, atol=5e-4)
    assert int((deep_r >= 0.8).sum()) == 0

    entries = pd.read_csv(ROOT / "results/phase18/PHASE18_ADIM5_KENTER_81.csv")
    ratios: dict[str, dict[str, float]] = {}
    for model in ("CoxPH-v2", "RSF"):
        group = entries[(entries.model == model) & (~entries.right_censored_128)].copy()
        ratio = group.k_enter * SCALES[model] / group.M_j
        ratios[model] = {
            "count": int(len(ratio)),
            "median": float(ratio.median()),
            "minimum": float(ratio.min()),
            "maximum": float(ratio.max()),
        }
    assert np.isclose(ratios["CoxPH-v2"]["median"], 4.41, atol=0.02)
    assert np.isclose(ratios["RSF"]["minimum"], 8.67, atol=0.02)
    assert np.isclose(ratios["RSF"]["maximum"], 43.0, atol=0.1)

    rsf_base = base[base.model == "RSF"].reset_index(drop=True)
    rsf_item = data["RSF"]
    rsf_ihat, rsf_n_leaves, rsf_n = n_cache["RSF"]
    rsf_rows: list[dict[str, object]] = []
    for j, feature in enumerate(rsf_item["features"]):
        if rsf_base.loc[rsf_base.feature == feature, "bottleneck6"].iloc[0] != "N":
            continue
        rsf_rows.append({
            "feature": feature,
            "I_hat_TCGA": rsf_ihat[0, j],
            "I_hat_METABRIC": rsf_ihat[1, j],
            "threshold_kappa_s_m": KAPPA * SCALES["RSF"],
            "N_TCGA": rsf_n_leaves[0, j],
            "N_METABRIC": rsf_n_leaves[1, j],
            "N": rsf_n[j],
            "binding_at_k64_rho1": True,
            "query36_status": "yapılamadı",
            "query36_reason": "Uyumlu, aynı panel ve açıklayıcı bütçeyle saklanmış 36-sorgu attribution tensörü bulunmuyor.",
        })
    rsf_n_table = pd.DataFrame(rsf_rows)
    assert len(rsf_n_table) == 17
    rsf_n_table.to_csv(OUT / "PHASE19_ADIM3_RSF_N.csv", index=False)

    summary = (
        output.groupby(["model", "rho"])
        .agg(point_adequate=("point_adequate6", "sum"))
        .reset_index()
    )
    binding = (
        output.groupby(["model", "rho"])["bottleneck6"]
        .value_counts()
        .rename("count")
        .reset_index()
    )
    metadata = {
        "row_count": len(output),
        "expected_row_count": 405,
        "rsf_n_rows": len(rsf_n_table),
        "ratios": ratios,
        "deep_r_range_rho1": [float(deep_r.min()), float(deep_r.max())],
        "deep_point_at_rho16": int(summary[(summary.model == "DeepSurv") & (summary.rho == 16.0)].point_adequate.iloc[0]),
        "formula_sha256": {name: sha256(ROOT / "warrantxai" / name) for name in FORMULA_FILES},
        "elapsed_seconds": time.perf_counter() - started,
    }
    (OUT / "PHASE19_ADIM3_METADATA.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(binding.to_string(index=False))
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
