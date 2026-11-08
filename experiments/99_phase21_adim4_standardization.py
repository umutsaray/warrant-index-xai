"""Phase 21 ADIM 4: reference-model METABRIC standardization diagnostic."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "phase21_reviewer1"
DATA = ROOT / "results/phase7_retrain"
PANEL = ROOT / "results/phase72_panel/query_matrices.npz"
SEED = 20260917

sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]
from warrantxai.components import generalizability, generalizability_signed


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def raw_external_queries(active: list[str], official_full: np.ndarray, official_q: np.ndarray):
    phase3 = load("06_brca_real.py", "p21a4_phase3")
    phase3.RAW = ROOT / "data/external/metabric_phase72"
    train = pd.read_parquet(DATA / "X_train_selected.parquet").set_index("sample_key")
    external, _, _ = phase3.external_columns(list(train.columns))
    raw = np.column_stack([external[f] for f in active])
    raw = raw[np.isfinite(raw).all(axis=1)]
    idx = np.array([train.columns.get_loc(f) for f in active])
    mu = train.mean().to_numpy()[idx]
    sd = train.std(ddof=0).to_numpy()[idx]
    aligned = mu + (raw - raw.mean(0)) / raw.std(0) * sd
    np.testing.assert_allclose(aligned, official_full, atol=1e-8, rtol=0)
    selected = []
    for row in official_q:
        distance = np.max(np.abs(aligned - row), axis=1)
        hit = int(np.argmin(distance))
        if distance[hit] > 1e-8:
            raise RuntimeError("Stored METABRIC query could not be mapped to raw cohort row")
        selected.append(hit)
    return raw[np.asarray(selected)]


def g_scores(source: np.ndarray, external: np.ndarray):
    # Single-reference-model tensors are (query, feature).
    magnitude = generalizability(
        np.stack([np.abs(source).mean(axis=0), np.abs(external).mean(axis=0)]), .5)
    signed = generalizability_signed(
        np.stack([np.sort(source, axis=0), np.sort(external, axis=0)]), .5)
    return magnitude, signed, np.minimum(magnitude, signed)


def attach_rows(rows, model, active, method, source, external, official):
    gm, gs, gc = g_scores(source, external)
    for j, feature in enumerate(active):
        base = official.loc[feature]
        values = np.array([base.F, base.S, base.U, gc[j], base.R], float)
        rows.append({
            "model": model, "feature": feature, "standardization": method,
            "scope": "single reference model; no bootstrap replicates; F/S/U/R frozen",
            "G_magnitude": gm[j], "G_signed": gs[j], "G_claim": gc[j],
            "F_frozen": base.F, "S_frozen": base.S, "U_frozen": base.U,
            "R_frozen": base.R, "diagnostic_WI": values.min(),
            "point_adequate": values.min() >= .8,
            "bottleneck": ["F", "S", "U", "G_claim", "R"][int(values.argmin())],
        })


def main():
    started = perf_counter()
    OUT.mkdir(exist_ok=True)
    with np.load(PANEL, allow_pickle=True) as z:
        active = z["active_features"].astype(str).tolist()
        source_q = z["source_query12"].astype(np.float32)
        official_q64 = z["external_query12"].astype(float)
        official_q = official_q64.astype(np.float32)
        official_full = z["external_full"]
    raw_q = raw_external_queries(active, official_full, official_q64)
    train = pd.read_parquet(DATA / "X_train_selected.parquet").set_index("sample_key")
    rows = []

    # Cox values were already exactly reconstructed in Phase 15 from all 50 replicas.
    cox = pd.read_csv(ROOT / "results/phase15/PHASE15_ADIM3_G_STANDARDIZATION.csv")
    names = {"official_cohort_standardized": "cohort_specific_official",
             "tcga_train_scaler_only_no_cohort_standardization": "tcga_train_mean_sd_only"}
    cox = cox[cox.method.isin(names)].copy()
    cox["model"] = "CoxPH-v2"
    cox["standardization"] = cox.method.map(names)
    cox["scope"] = "stored N=50 Cox diagnostic from Phase 15"
    cox["diagnostic_WI"] = cox.WI
    cox = cox.rename(columns={"G_magnitude": "G_magnitude", "G_signed": "G_signed"})

    # RSF reference model, frozen sampled-background and 8-path candidate design.
    rsf_helpers = load("07_rsf_pilot.py", "p21a4_rsf_helpers")
    rsf_bundle = joblib.load(DATA / "rsf_baseline.joblib")
    rsf_features = list(rsf_bundle["features"])
    rsf_idx = np.array([rsf_features.index(f) for f in active])
    unique = train.loc[~train.index.astype(str).str[:12].duplicated()]
    rng = np.random.default_rng(SEED + 4200)
    pick16 = np.sort(rng.choice(len(unique), 16, replace=False))
    backgrounds = unique.iloc[pick16].to_numpy(np.float32)
    active_mean = train.mean().to_numpy(float)[rsf_idx]
    rsf_source = rsf_helpers.background_permutation_attributions(
        rsf_bundle["model"], source_q, backgrounds, active_mean, rsf_idx,
        permutations=8, seed=SEED + 10000)
    rsf_official = rsf_helpers.background_permutation_attributions(
        rsf_bundle["model"], official_q, backgrounds, active_mean, rsf_idx,
        permutations=8, seed=SEED + 10100)
    rsf_raw = rsf_helpers.background_permutation_attributions(
        rsf_bundle["model"], raw_q, backgrounds, active_mean, rsf_idx,
        permutations=8, seed=SEED + 10100)
    rsf_report = pd.read_csv(ROOT / "results/phase72_rsf/feature_report.csv").set_index("feature")
    attach_rows(rows, "RSF", active, "cohort_specific_official", rsf_source, rsf_official, rsf_report)
    attach_rows(rows, "RSF", active, "tcga_train_mean_sd_only", rsf_source, rsf_raw, rsf_report)

    # DeepSurv reference model with the frozen k=16 candidate background.
    phase6 = load("20_deepsurv_pilot.py", "p21a4_phase6")
    phase61 = load("21_deepsurv_explainer_revision.py", "p21a4_phase61")
    phase72 = load("41_phase72_deepsurv.py", "p21a4_phase72")
    checkpoint = torch.load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    model = phase6.model_from_state(checkpoint)
    deep_features = list(checkpoint["ordered_columns"])
    deep_idx = np.array([deep_features.index(f) for f in active])
    deep_train = train[deep_features]
    deep_unique = deep_train.loc[~deep_train.index.astype(str).str[:12].duplicated()]
    rng = np.random.default_rng(SEED + 4200)
    pick16 = np.sort(rng.choice(len(deep_unique), 16, replace=False))
    deep_bg = deep_unique.iloc[pick16].to_numpy(np.float32)
    args = (phase61, phase6, model, checkpoint)
    means = np.asarray(checkpoint["means"], np.float32)
    scales = np.asarray(checkpoint["scales"], np.float32)
    deep_source = phase72.deep_attributions_with_relative_gate(
        *args, source_q, deep_bg, deep_idx, means, scales)[0]
    deep_official = phase72.deep_attributions_with_relative_gate(
        *args, official_q, deep_bg, deep_idx, means, scales)[0]
    deep_raw = phase72.deep_attributions_with_relative_gate(
        *args, raw_q.astype(np.float32), deep_bg, deep_idx, means, scales)[0]
    deep_report = pd.read_csv(ROOT / "results/phase72_deepsurv/feature_report.csv").set_index("feature")
    attach_rows(rows, "DeepSurv", active, "cohort_specific_official", deep_source, deep_official, deep_report)
    attach_rows(rows, "DeepSurv", active, "tcga_train_mean_sd_only", deep_source, deep_raw, deep_report)

    result = pd.concat([cox.reindex(columns=[
        "model", "feature", "standardization", "scope", "G_magnitude", "G_signed",
        "G_claim", "diagnostic_WI", "point_adequate", "bottleneck"]), pd.DataFrame(rows)],
        ignore_index=True, sort=False)
    result.to_csv(OUT / "PHASE21_ADIM4_STANDARDIZATION.csv", index=False)
    summary = (result.groupby(["model", "standardization"])
               .agg(n=("feature", "size"), point=("point_adequate", "sum"),
                    median_G_claim=("G_claim", "median"), median_WI=("diagnostic_WI", "median"))
               .reset_index())
    bn = result.groupby(["model", "standardization", "bottleneck"]).size().rename("n").reset_index()
    summary.to_csv(OUT / "PHASE21_ADIM4_STANDARDIZATION_SUMMARY.csv", index=False)
    text = f"""# Faz 21 ADIM 4 — METABRIC standardizasyon duyarlılığı

RSF ve DeepSurv için bootstrap model nesneleri saklı olmadığından, izin verilen **tek referans model / replika-yok** yedeği kullanıldı. Yeni model fit edilmedi. F/S/U/R resmi Faz 7.2 satırlarından sabit tutuldu; yalnız referans modelin aynı 12+12 sorgudaki G_claim değeri değiştirilerek tanısal WI verildi. CoxPH karşılaştırması saklı N=50 Faz 15 yeniden üretimidir. Dolayısıyla modeller arası inferans kapsamları aynı değildir ve bu tablo resmi sertifikasyon değildir.

```text
{summary.to_string(index=False)}
```

Darboğaz sayımları:

```text
{bn.to_string(index=False)}
```

Wall-clock: {perf_counter()-started:.3f} sn. Dondurulmuş sonuç ve formül dosyaları değiştirilmedi.
"""
    (OUT / "PHASE21_ADIM4_STANDARDIZATION.md").write_text(text, encoding="utf-8")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
