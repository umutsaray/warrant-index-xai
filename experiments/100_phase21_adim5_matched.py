"""Phase 21 ADIM 5: matched Cox audit and DeepSurv background sensitivity."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
import torch
from scipy.stats import pearsonr, spearmanr

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "phase21_reviewer1"
DATA = ROOT / "results/phase7_retrain"
SEED, DRAWS = 20260918, 4999
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]

from warrantxai.certification import certify
from warrantxai.components import generalizability, generalizability_signed


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def cox_matched(active, source_q, external_q):
    shared = load("12_rsf_certification.py", "p21a5_shared")
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as z:
        coefficients = z["coefficients"]
        features = z["features"].astype(str).tolist()
    bundle = joblib.load(DATA / "coxph_baseline_v2.joblib")
    scaler = bundle["scaler"]
    idx = np.array([features.index(f) for f in active])
    source_z = (source_q - scaler.mean_[idx]) / scaler.scale_[idx]
    external_z = (external_q - scaler.mean_[idx]) / scaler.scale_[idx]
    coef = coefficients[:, idx]
    phi = np.stack([
        coef[:, None, :] * source_z[None, :, :],
        coef[:, None, :] * external_z[None, :, :],
    ])
    rng = np.random.default_rng(SEED + 725_000)
    indices = rng.integers(0, phi.shape[1], size=(DRAWS, phi.shape[1]))
    points, pvalues = shared.component_statistics(phi, phi, indices)
    bh, by = certify(pvalues, .05, "bh"), certify(pvalues, .05, "by")
    names = ["F", "S", "U", "G_claim", "R"]
    rows = []
    for j, feature in enumerate(active):
        rows.append({
            "feature": feature, **dict(zip(names, points[j])), "WI": points[j].min(),
            "bottleneck": names[int(points[j].argmin())], "point_adequate": points[j].min() >= .8,
            **{f"p_{name}": pvalues[j, k] for k, name in enumerate(names)},
            "p_joint": bh.joint_p[j], "q_BH": bh.adjusted_p[j], "q_BY": by.adjusted_p[j],
            "certified_BH": bh.rejected[j], "certified_BY": by.rejected[j],
        })
    return pd.DataFrame(rows)


def safe_corr(x, y, kind):
    if np.std(x) == 0 or np.std(y) == 0:
        return 1.0 if np.array_equal(x, y) else np.nan
    return float((pearsonr if kind == "pearson" else spearmanr)(x, y).statistic)


def deep_background(active, queries):
    phase6 = load("20_deepsurv_pilot.py", "p21a5_phase6")
    phase61 = load("21_deepsurv_explainer_revision.py", "p21a5_phase61")
    phase72 = load("41_phase72_deepsurv.py", "p21a5_phase72")
    checkpoint = torch.load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    model = phase6.model_from_state(checkpoint)
    features = list(checkpoint["ordered_columns"])
    idx = np.array([features.index(f) for f in active])
    train = pd.read_parquet(DATA / "X_train_selected.parquet").set_index("sample_key")[features]
    unique = train.loc[~train.index.astype(str).str[:12].duplicated()]
    rng = np.random.default_rng(20260917 + 4200)
    first16 = np.sort(rng.choice(len(unique), 16, replace=False))
    remaining = np.setdiff1d(np.arange(len(unique)), first16)
    # Preserve the frozen k=128 set exactly; only k=256 extends it.
    extra112 = np.sort(rng.choice(remaining, 112, replace=False))
    used128 = np.r_[first16, extra112]
    remaining2 = np.setdiff1d(np.arange(len(unique)), used128)
    extra128 = np.sort(rng.choice(remaining2, 128, replace=False))
    nested = np.r_[used128, extra128]
    means = np.asarray(checkpoint["means"], np.float32)
    scales = np.asarray(checkpoint["scales"], np.float32)
    tensors, diagnostics = {}, []
    for k in (16, 64, 128, 256):
        bg = unique.iloc[nested[:k]].to_numpy(np.float32)
        sites = []
        for site, query in enumerate(queries):
            values, error, output_range, seconds, warnings, fallback = \
                phase72.deep_attributions_with_relative_gate(
                    phase61, phase6, model, checkpoint, query.astype(np.float32), bg,
                    idx, means, scales)
            sites.append(values)
            diagnostics.append({"background_k": k, "site": site, "max_abs_additivity_error": error,
                                "output_range": output_range, "relative_additivity_error": error/max(output_range, np.finfo(float).eps),
                                "seconds": seconds, "fallback_relative_gate": fallback,
                                "warning_count": len(warnings)})
        tensors[k] = np.stack(sites)
    ref = tensors[128]
    rows = []
    for k, tensor in tensors.items():
        gm = generalizability(np.abs(tensor).mean(axis=1), .5)
        gs = generalizability_signed(np.stack([np.sort(tensor[0], axis=0), np.sort(tensor[1], axis=0)]), .5)
        gc = np.minimum(gm, gs)
        for j, feature in enumerate(active):
            x, y = tensor[:, :, j].reshape(-1), ref[:, :, j].reshape(-1)
            denom = max(float(np.mean(np.abs(y))), np.finfo(float).eps)
            rows.append({"feature": feature, "background_k": k,
                         "pearson_vs_k128": safe_corr(x, y, "pearson"),
                         "spearman_vs_k128": safe_corr(x, y, "spearman"),
                         "mean_abs_relative_difference_vs_k128": float(np.mean(np.abs(x-y))/denom),
                         "G_magnitude": gm[j], "G_signed": gs[j], "G_claim": gc[j]})
    return pd.DataFrame(rows), pd.DataFrame(diagnostics)


def main():
    started = perf_counter()
    OUT.mkdir(exist_ok=True)
    with np.load(ROOT / "results/phase72_panel/query_matrices.npz", allow_pickle=True) as z:
        active = z["active_features"].astype(str).tolist()
        source = z["source_query12"].astype(float)
        external = z["external_query12"].astype(float)
    cox = cox_matched(active, source, external)
    cox.to_csv(OUT / "PHASE21_ADIM5_COX_MATCHED.csv", index=False)
    deep, diagnostics = deep_background(active, [source, external])
    deep.to_csv(OUT / "PHASE21_ADIM5_DEEPSURV_BACKGROUND.csv", index=False)
    diagnostics.to_csv(OUT / "PHASE21_ADIM5_DEEPSURV_BACKGROUND_DIAGNOSTICS.csv", index=False)
    cox_summary = pd.DataFrame([{
        "family_m": len(cox), "point": int(cox.point_adequate.sum()),
        "BH": int(cox.certified_BH.sum()), "BY": int(cox.certified_BY.sum()),
        "bottlenecks": str(cox.bottleneck.value_counts().to_dict())}])
    deep_summary = (deep.groupby("background_k")
                    .agg(median_pearson=("pearson_vs_k128", "median"),
                         min_pearson=("pearson_vs_k128", "min"),
                         median_spearman=("spearman_vs_k128", "median"),
                         median_MARD=("mean_abs_relative_difference_vs_k128", "median"),
                         max_MARD=("mean_abs_relative_difference_vs_k128", "max"),
                         median_G_claim=("G_claim", "median")).reset_index())
    text = f"""# Faz 21 ADIM 5 — eşleştirilmiş Cox ve DeepSurv background tanıları

## CoxPH-v2: aynı 12+12 sorgu, aynı 27-feature sıralama evreni

Saklı 50 katsayı replikası kullanıldı; yeni model fit edilmedi. F dahil bütün bileşenler iki kohortlu 12+12 tensör üzerinde, R yalnız 27-feature evreninde hesaplandı. Recursive-studentized B=4,999 ve 27-feature BH/BY ailesi kullanıldı.

```text
{cox_summary.to_string(index=False)}
```

## DeepSurv: background k duyarlılığı

Tek dondurulmuş referans modelde aynı sorgular ve deterministik nested k=16/64/128/256 TCGA background'ları kullanıldı. Korelasyonlar ve MARD, özellik başına iki kohortun 24 attribution değerinin k=128 vektörüyle karşılaştırmasıdır; bu bölüm sertifikasyon değildir.

```text
{deep_summary.to_string(index=False)}
```

Wall-clock: {perf_counter()-started:.3f} sn. Dondurulmuş formüller ve resmi sonuçlar değiştirilmedi.
"""
    (OUT / "PHASE21_ADIM5_MATCHED_AND_BACKGROUND.md").write_text(text, encoding="utf-8")
    print(cox_summary.to_string(index=False))
    print(deep_summary.to_string(index=False))


if __name__ == "__main__":
    main()
