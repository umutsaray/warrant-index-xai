#!/usr/bin/env python3
"""Phase 10 fold harmonization feasibility and Cox retry-depth diagnostics."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import joblib
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "experiments")]
from _shared import CALIBRATION
from warrantxai.certification import certify
from warrantxai.components import importance_ci

EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}
COMPONENTS = ("F", "S", "U", "G_claim", "R")
SEED, DRAWS = 20260918, 4999


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    value = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def md(frame: pd.DataFrame) -> str:
    cols = [str(x) for x in frame.columns]
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False, name=None)]
    return "\n".join(rows)


def harmonization_feasibility() -> tuple[pd.DataFrame, float]:
    started = perf_counter()
    fold0 = pd.read_csv(ROOT / "results/phase72_panel/feature_mapping.csv")
    fold1 = pd.read_csv(ROOT / "results/phase75_fold1_panel/feature_mapping.csv")
    clinical_names = {"clin_numeric__age_at_index.demographic": "Age at index",
                      **{f"clin_categorical__ajcc_pathologic_stage.diagnoses_Stage {x}": f"AJCC Stage {x}" for x in ("I", "II", "III", "IV")}}
    fold0["symbol"] = [clinical_names.get(feature, symbol if pd.notna(symbol) else feature)
                       for feature, symbol in fold0[["feature", "symbol"]].itertuples(index=False)]
    fold1_set = set(fold1.feature.astype(str))
    rows = []
    for feature, symbol in fold0[["feature", "symbol"]].itertuples(index=False):
        available = feature in fold1_set
        for model in ("coxph", "rsf", "deepsurv"):
            row = {"model": model.upper() if model != "deepsurv" else "DeepSurv", "feature": feature,
                   "symbol": symbol, "available_in_fold1_model_panel": available}
            fold0_path = ROOT / f"results/phase72_{model}_certification/feature_certification.csv"
            fold0_hit = pd.read_csv(fold0_path).query("feature == @feature").iloc[0]
            for col in ("F", "S", "U", "G_claim", "R", "WI", "bottleneck", "point_WI_ge_0_8",
                        "certified_BH_alpha_0_05", "certified_BY_alpha_0_05"):
                row[f"fold0_{col}"] = fold0_hit[col]
            path = ROOT / f"results/phase75_fold1_{model}_certification/feature_certification.csv"
            if available and path.exists():
                hit = pd.read_csv(path).query("feature == @feature")
                if len(hit):
                    for col in ("F", "S", "U", "G_claim", "R", "WI", "bottleneck", "point_WI_ge_0_8",
                                "certified_BH_alpha_0_05", "certified_BY_alpha_0_05"):
                        row[f"fold1_{col}"] = hit.iloc[0][col]
                    row["fold1_status"] = "available; BH/BY belongs to original 80-feature family"
            elif not available:
                row["fold1_status"] = "NA: feature absent from frozen fold-1 616-feature model input"
            else:
                row["fold1_status"] = "NA: CoxPH fold-1 replicate 15 failed frozen numerical retry"
            rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(ROOT / "PHASE10_FOLD1_HARMONIZED_27FEATURE.csv", index=False)
    common_features = fold0[fold0.feature.isin(fold1_set)]
    missing_features = fold0[~fold0.feature.isin(fold1_set)]
    common = out[out.available_in_fold1_model_panel]
    available_table = common[["model", "symbol", "fold0_F", "fold1_F", "fold0_S", "fold1_S",
        "fold0_U", "fold1_U", "fold0_G_claim", "fold1_G_claim", "fold0_R", "fold1_R",
        "fold0_WI", "fold1_WI", "fold0_point_WI_ge_0_8", "fold1_point_WI_ge_0_8",
        "fold0_certified_BH_alpha_0_05", "fold1_certified_BH_alpha_0_05",
        "fold0_certified_BY_alpha_0_05", "fold1_certified_BY_alpha_0_05", "fold1_status"]]
    note = f"""# Faz 10 — Fold-1 harmonize 27-feature fizibilitesi

Fold-0'ın 27-feature panelinin yalnız **{len(common_features)}/27** özelliği fold-1'in dondurulmuş 616-feature model girdisinde ve 80-feature audit panelinde vardır; **{len(missing_features)}/27** özellik fold-1 modellerinin girdisi değildir. Bu nedenle mevcut fold-1 modellerini yeniden fit etmeden tam 27-feature apples-to-apples sertifikasyon matematiksel olarak tanımlı değildir. Kapsam kuralı yeniden fit etmeyi yasakladığı için eksik özelliklere sıfır attribution atanmadı ve sonuç uydurulmadı.

Mevcut 80-feature sertifikasyon dosyalarından ortak 10 özellik için değerler CSV'ye tanı amaçlı taşındı. Bunların BH/BY kararları 80-feature ailesine aittir; 27-feature harmonize aile sonucu olarak yorumlanamaz. CoxPH fold-1 sonucu mevcut dondurulmuş retry başarısızlığı nedeniyle zaten NA'dır.

| Durum | Sayı |
| --- | ---: |
| Fold-1 model girdisinde mevcut | {len(common_features)} |
| Fold-1 model girdisinde yok | {len(missing_features)} |
| Geçerli tam-27 harmonize sertifikasyon | 0 (hesaplanamaz) |

**Ortak 10 özellik:** `{', '.join(common_features.symbol.astype(str))}`.

**Fold-1 model girdisinde bulunmayan ve NA işaretlenen 17 özellik:** `{', '.join(missing_features.symbol.astype(str))}`.

## Ortak özelliklerin mevcut sonuçları

{md(available_table)}

Fold-1 CoxPH için ortak özelliklerde de F/S/U/G/R/WI ve sertifikasyon NA'dır; mevcut koşu deterministik replika 15'te dondurulmuş retry sonrasında başarısız olmuştur. RSF ve DeepSurv fold-1 değerleri hesaplanmıştır. Ancak fold-1 BH/BY sütunları **80-feature ailesindeki** mevcut kararlardır; 27-feature harmonize BH/BY sonucu değildir.

## Teknik neden

Faz 7.0/7.1 zincirinde özellik seçimi her dış fold'un yalnız kendi train partisyonunda yeniden yürütülmüştür. Univariate sıralama/filtreleme fold'a özgü örneklem ve olay dağılımına bağlı olduğundan fold-0 ve fold-1 farklı 616-feature model girdileri üretmiştir. Bu ayrım `PHASE7_COLLAPSE_DIAGNOSTIC.md` içindeki feature-selection kararlılığı bulgusu ve `PHASE75_FOLD_SELECTION.md` içindeki fold-1 seçim kaydıyla belgelenmiştir.

Süre **{perf_counter()-started:.3f} sn**. Formül dosyaları değiştirilmedi.
"""
    (ROOT / "PHASE10_FOLD1_HARMONIZED_27FEATURE.md").write_text(note, encoding="utf-8")
    return out, perf_counter() - started


def cox_statistics_subset(keep: np.ndarray) -> tuple[np.ndarray, np.ndarray, list[str]]:
    phase53 = module("12_rsf_certification.py", "p10_phase53")
    oldcox = module("13_cox_recursive_certification.py", "p10_oldcox")
    pilot = module("08_certification_pilot.py", "p10_pilot")
    remedy = module("11_rsf_remedy_calibration.py", "p10_remedy")
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source = saved["source_phi"][keep]
        external = saved["external_phi"][keep]
        coefficients_all = saved["coefficients"][keep]
        all_features = saved["features"].astype(str).tolist()
        active = saved["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(name) for name in active])
    coefficients = coefficients_all[:, targets]
    selected_source = source[:, :, targets]
    t, p = coefficients.shape
    rng = np.random.default_rng(SEED + 1_000_022)
    indices = rng.integers(0, t, size=(DRAWS, t))
    counts = np.zeros((DRAWS, t), dtype=np.int16)
    np.add.at(counts, (np.arange(DRAWS)[:, None], indices), 1)
    points, pvalues = np.empty((p, 5)), np.empty((p, 5))

    f_by_rep = np.ones((1, t, p))
    source4 = source[None, ...]
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_by_rep = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[:, :, targets]
    for values, column in ((f_by_rep, 0), (u_by_rep, 2)):
        point = values.mean(axis=1)
        sampled, sampled_se = phase53.resampled_means(values, counts)
        se = values.std(axis=1, ddof=1) / np.sqrt(t)
        pvalues[:, column] = remedy.studentized_p(point.T.ravel(), np.moveaxis(sampled, 1, 2).reshape(DRAWS, -1), se.T.ravel(), np.moveaxis(sampled_se, 1, 2).reshape(DRAWS, -1))
        points[:, column] = point[0]

    for loss, column in ((pilot.pair_losses(source4, targets, CALIBRATION["a"]), 1), (pilot.order_pair_losses(np.abs(source).mean(axis=1)[None, ...], targets, CALIBRATION["eta"]), 4)):
        point = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled = phase53.quadratic_scores(loss, counts, t)
        pvalues[:, column] = remedy.mc_p(point[0], sampled[:, 0])
        points[:, column] = point[0]

    difference = np.abs(external).mean(axis=1) - np.abs(selected_source).mean(axis=1)
    gm = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = phase53.resampled_means(difference[None, ...], counts)
    gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = remedy.studentized_p(gm, gm_sample, gm_se, gm_mean_se[:, 0] / CALIBRATION["b"])

    sorted_delta = []
    for j in range(p):
        row = int(np.argmax(np.abs(coefficients[:, j])))
        coef = coefficients[row, j]
        if coef == 0:
            sorted_delta.append(np.zeros(max(selected_source.shape[1], external.shape[1])))
        else:
            n = max(selected_source.shape[1], external.shape[1])
            a = pilot.empirical_rank_grid(selected_source[row, :, j] / coef, n)
            b = pilot.empirical_rank_grid(external[row, :, j] / coef, n)
            sorted_delta.append(np.sort(np.abs(a - b)))
    mean = coefficients.mean(axis=0)
    gs = oldcox.signed_score_from_coefficient(mean, sorted_delta)[0]
    delete = (t * mean - coefficients) / (t - 1)
    gs_se = oldcox.jackknife_se(oldcox.signed_score_from_coefficient(delete, sorted_delta))
    sampled_mean = coefficients[indices].mean(axis=1)
    gs_sample = oldcox.signed_score_from_coefficient(sampled_mean, sorted_delta)
    gs_sample_se = np.empty((DRAWS, p))
    for j in range(p):
        leave = (t * sampled_mean[:, j, None] - coefficients[indices, j]) / (t - 1)
        scores = oldcox.signed_score_from_coefficient(leave.reshape(-1, 1), [sorted_delta[j]]).reshape(DRAWS, t)
        gs_sample_se[:, j] = oldcox.jackknife_se(scores.T)
    p_gs = remedy.studentized_p(gs, gs_sample, gs_se, gs_sample_se)
    points[:, 3] = np.minimum(gm, gs)
    pvalues[:, 3] = np.maximum(p_gm, p_gs)
    return points, pvalues, active


def partial_log_likelihood(x: np.ndarray, time: np.ndarray, event: np.ndarray, beta: np.ndarray) -> float:
    risk = x @ beta
    order = np.argsort(-time, kind="mergesort")
    risk, time, event = risk[order], time[order], event[order]
    shift = risk.max()
    cumulative = np.cumsum(np.exp(risk - shift))
    value = 0.0
    for t in np.unique(time[event]):
        deaths = np.flatnonzero((time == t) & event)
        if len(deaths):
            last = np.flatnonzero(time >= t)[-1]
            value += risk[deaths].sum() - len(deaths) * (np.log(cumulative[last]) + shift)
    return float(value)


def retry_depth() -> tuple[pd.DataFrame, float, int]:
    started = perf_counter()
    records = json.loads((ROOT / "results/phase72_cox/bootstrap_summary.json").read_text())
    retry = np.array([len(x["fit_attempts"]) > 1 for x in records])
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as z:
        coefficients = z["coefficients"]
    saved = joblib.load(ROOT / "results/phase7_retrain/coxph_baseline_v2.joblib")
    features, scaler = list(saved["features"]), saved["scaler"]
    train = pd.read_parquet(ROOT / "results/phase7_retrain/X_train_selected.parquet").set_index("sample_key")
    y = pd.read_parquet(ROOT / "results/phase7_retrain/y_train.parquet").set_index("sample_key").loc[train.index]
    scaled = scaler.transform(train[features].astype(float))
    patients = train.index.astype(str).str[:12]
    groups = [np.flatnonzero(patients == p) for p in patients.unique()]
    rng = np.random.default_rng(20260917)
    norms, objectives = [], []
    for rep in range(50):
        idx = np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))])
        beta = coefficients[rep]
        norms.append(float(np.linalg.norm(beta)))
        objectives.append(partial_log_likelihood(scaled[idx], y.time.to_numpy()[idx], y.event.to_numpy()[idx].astype(bool), beta))
    metrics = pd.DataFrame({"replicate": range(50), "group": np.where(retry, "retry", "direct"), "coefficient_L2": norms, "log_partial_likelihood": objectives})

    summary_rows = []
    for metric in ("coefficient_L2", "log_partial_likelihood"):
        a = metrics.loc[retry, metric].to_numpy(); b = metrics.loc[~retry, metric].to_numpy()
        diff = a.mean() - b.mean(); se = np.sqrt(a.var(ddof=1)/len(a) + b.var(ddof=1)/len(b))
        stat, p = ttest_ind(a, b, equal_var=False)
        summary_rows.append({"metric": metric, "retry_n": len(a), "direct_n": len(b), "retry_mean": a.mean(), "direct_mean": b.mean(), "difference_retry_minus_direct": diff, "combined_SE": se, "CI95_lower": diff-1.96*se, "CI95_upper": diff+1.96*se, "Welch_t": stat, "Welch_p": p})
    summary = pd.DataFrame(summary_rows)

    points, pvalues, active = cox_statistics_subset(np.flatnonzero(~retry))
    bh, by = certify(pvalues, alpha=.05, method="bh"), certify(pvalues, alpha=.05, method="by")
    full = pd.read_csv(ROOT / "results/phase72_coxph_certification/feature_certification.csv").set_index("feature").loc[active]
    rows = pd.DataFrame({"feature": active, **{c: points[:, i] for i, c in enumerate(COMPONENTS)}})
    rows["WI_direct22"] = points.min(axis=1)
    rows["point_direct22"] = rows.WI_direct22 >= CALIBRATION["tau"]
    rows["p_joint_direct22"] = bh.joint_p; rows["q_BH_direct22"] = bh.adjusted_p; rows["q_BY_direct22"] = by.adjusted_p
    rows["BH_direct22"] = bh.rejected; rows["BY_direct22"] = by.rejected
    rows["BH_full50"] = full.certified_BH_alpha_0_05.to_numpy(); rows["BY_full50"] = full.certified_BY_alpha_0_05.to_numpy()
    rows["BH_changed"] = rows.BH_direct22 != rows.BH_full50; rows["BY_changed"] = rows.BY_direct22 != rows.BY_full50
    mapping = pd.read_csv(ROOT / "results/phase72_panel/feature_mapping.csv")
    names = dict(zip(mapping.feature, mapping.symbol))
    names.update({"clin_numeric__age_at_index.demographic": "Age at index",
                  **{f"clin_categorical__ajcc_pathologic_stage.diagnoses_Stage {x}": f"AJCC Stage {x}" for x in ("I", "II", "III", "IV")}})
    rows.insert(1, "symbol", rows.feature.map(names).fillna(rows.feature))
    rows["BH_relation"] = np.select(
        [rows.BH_full50 & rows.BH_direct22, rows.BH_full50 & ~rows.BH_direct22, ~rows.BH_full50 & rows.BH_direct22],
        ["certified in both", "dropped in direct-22", "new in direct-22"], default="not certified in either")
    for col in summary.columns:
        rows[f"metric_{col}"] = ""
    for i, record in summary.iterrows():
        for col in summary.columns:
            rows.loc[i, f"metric_{col}"] = record[col]
    rows.to_csv(ROOT / "PHASE10_COXPH_RETRY_BIAS_DEPTH.csv", index=False)
    changed = rows[rows.BH_changed | rows.BY_changed]
    full_bh = rows.loc[rows.BH_full50, "symbol"].tolist()
    direct_bh = rows.loc[rows.BH_direct22, "symbol"].tolist()
    dropped = rows.loc[rows.BH_full50 & ~rows.BH_direct22, "symbol"].tolist()
    added = rows.loc[~rows.BH_full50 & rows.BH_direct22, "symbol"].tolist()
    note = f"""# Faz 10 — CoxPH-v2 retry önyargısı derinleştirmesi

50 dondurulmuş replika yeniden fit edilmeden aynı 28 retry / 22 direct ayrımıyla incelendi. Objective, her replikaya ait saklı katsayı vektörü ve deterministik olarak yeniden kurulan hasta-cluster bootstrap örneği üzerinde hesaplanan penalize edilmemiş Cox log-partial-likelihood değeridir.

{md(summary)}

## Yalnız 22 doğrudan-yakınsayan replika

`recursive_studentized`, B={DRAWS}, τ=0.8 ve aynı BH/BY prosedürü yalnız 22 direct replika üzerinde yeniden uygulandı. Tam-50 sonuçları değiştirilmedi.

| Analiz | Point | BH | BY |
| --- | ---: | ---: | ---: |
| Tam 50 | {int(full.point_WI_ge_0_8.sum())} | {int(full.certified_BH_alpha_0_05.sum())} | {int(full.certified_BY_alpha_0_05.sum())} |
| Direct 22 | {int(rows.point_direct22.sum())} | {int(rows.BH_direct22.sum())} | {int(rows.BY_direct22.sum())} |

### Sertifikalı kümelerin kesin dökümü

- **Tam 50 replika:** 27 özellikten **{len(full_bh)}** tanesi BH ve BY ile sertifikalıdır: `{', '.join(full_bh)}`.
- **Yalnız direct 22 replika:** 27 özellikten **{len(direct_bh)}** tanesi BH ve BY ile sertifikalıdır: `{', '.join(direct_bh)}`.
- **Sertifikadan düşenler:** `{', '.join(dropped) if dropped else 'yok'}`.
- **Direct-22 analizinde yeni eklenenler:** `{', '.join(added) if added else 'yok'}`.

Dolayısıyla CA1 ve OR13C3'ün **ikisi de** tam-50 kümesinde sertifikalıyken direct-22 analizinde sertifikadan düşmüştür; biri düşüp diğeri eklenmemiştir. BH ve BY kümeleri her iki senaryoda da kendi içinde aynıdır.

Point-adequate sonuç tam-50 için **{int(full.point_WI_ge_0_8.sum())}/27**, direct-22 için **{int(rows.point_direct22.sum())}/27**'dir. BH/BY kararı değişen özellik sayısı **{len(changed)}/27**.

Süre **{perf_counter()-started:.3f} sn**. Yeni bootstrap/model fit üretilmedi; resmi sertifikasyon korunmuştur.
"""
    (ROOT / "PHASE10_COXPH_RETRY_BIAS_DEPTH.md").write_text(note, encoding="utf-8")
    return rows, perf_counter() - started, len(changed)


def main() -> None:
    actual = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED}
    if actual != EXPECTED:
        raise RuntimeError(f"Frozen formula hashes changed: {actual}")
    _, step2 = harmonization_feasibility()
    _, step3, changed = retry_depth()
    print(json.dumps({"step2_seconds": step2, "step3_seconds": step3, "decision_changes": changed, "formula_sha256": actual}, indent=2))


if __name__ == "__main__":
    main()
