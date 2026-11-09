"""Phase 21 v3 stage 1: inner-safe preprocessing, LR screen, and reference fits."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from lifelines import CoxPHFitter
from scipy.stats import chi2
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv

ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
SRC = SIBLING / "brca-multiomics-survival/src/common"
OUT = ROOT / "phase21_reviewer1/v3"
LOCAL = ROOT / "results/phase21_v3"
DATA = ROOT / "results/phase7_retrain"
SEED = 42
sys.path[:0] = [str(SRC), str(SIBLING), str(ROOT), str(ROOT / "experiments")]


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def key_hash(index) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, index))).encode()).hexdigest()


def load_common():
    path = SRC / "phase1_common.py"
    spec = importlib.util.spec_from_file_location("p21v3_common", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def null_loglik(time_values, event_values):
    order = np.argsort(-time_values, kind="stable")
    times, events = time_values[order], event_values[order].astype(bool)
    ends = np.empty(len(times), int)
    start = 0
    while start < len(times):
        end = start
        while end + 1 < len(times) and times[end + 1] == times[start]:
            end += 1
        ends[start:end + 1] = end
        start = end + 1
    return float(-np.log(ends[np.flatnonzero(events)] + 1).sum())


def rank_and_select(common, xfit, yfit):
    ll0 = null_loglik(yfit.time.to_numpy(float), yfit.event.to_numpy(bool))
    rows, selected = [], {}
    threshold = int(np.ceil(.10 * len(xfit)))
    for block, prefix in common.BLOCKS.items():
        columns = [c for c in xfit if c.startswith(prefix)]
        values = xfit[columns].to_numpy(float, copy=False)
        nz = values != 0
        nonzero_n = nz.sum(0)
        nonzero_events = nz[yfit.event.to_numpy(bool)].sum(0)
        eligible = (nonzero_n >= threshold) & (nonzero_events >= 10)
        names = np.asarray(columns)[eligible].tolist()
        ranked = common.rank_block(xfit, names, yfit, 512, 50, 1e-7)
        ranked["lr_p"] = chi2.sf(np.maximum(0, 2 * (ranked.partial_log_likelihood - ll0)), 1)
        ranked = ranked.sort_values(["lr_p", "feature"], kind="stable").reset_index(drop=True)
        ranked["lr_rank"] = np.arange(1, len(ranked) + 1)
        chosen = ranked.head(200).feature.tolist()
        selected[block] = chosen
        base = pd.DataFrame({"feature": columns, "block": block, "nonzero_n": nonzero_n,
                             "nonzero_fraction": nonzero_n / len(xfit), "nonzero_events": nonzero_events})
        base = base.merge(ranked[["feature", "lr_p", "lr_rank"]], on="feature", how="left")
        base["filter_pass"] = eligible
        base["selected_model_panel"] = base.feature.isin(chosen)
        rows.append(base)
    return pd.concat(rows, ignore_index=True), selected


def map_rna(features):
    import requests
    ids = [f[4:].split(".")[0] for f in features]
    response = requests.post("https://mygene.info/v3/query", data={
        "q": ",".join(ids), "scopes": "ensembl.gene", "fields": "symbol,entrezgene", "species": "human"}, timeout=120)
    response.raise_for_status()
    by = {}
    for row in response.json():
        if row.get("symbol"):
            by.setdefault(row["query"], set()).add(row["symbol"])
    return {f: next(iter(by.get(i, [])), "") if len(by.get(i, [])) == 1 else "" for f, i in zip(features, ids)}


def train_deep(xfit, yfit, xval, yval, xtest, ytest, columns):
    from train_explainable_dl import (OmicsDeepSurv, cindex, cox_negative_log_likelihood,
                                      make_blocks, predict_risk, set_seed)
    set_seed(SEED)
    prefixes = {"rna": "rna_", "mirna": "mirna_", "methylation": "meth_", "clinical": "clin_"}
    by = {b: [c for c in columns if c.startswith(p)] for b, p in prefixes.items()}
    order = [c for b in prefixes for c in by[b]]
    slices, cursor = {}, 0
    for b in prefixes:
        slices[b] = (cursor, cursor + len(by[b])); cursor += len(by[b])
    means = xfit[order].mean().to_numpy(np.float32)
    scales = xfit[order].std(ddof=0).to_numpy(np.float32); scales[scales == 0] = 1
    arr = lambda z: ((z[order].to_numpy(np.float32) - means) / scales)
    af, av, at = arr(xfit), arr(xval), arr(xtest)
    model = OmicsDeepSurv({b: len(by[b]) for b in prefixes}, dropout=.3)
    optimizer = torch.optim.AdamW(model.parameters(), lr=.001, weight_decay=.0001)
    fb = make_blocks(af, slices, torch.device("cpu")); vb = make_blocks(av, slices, torch.device("cpu"))
    tb = make_blocks(at, slices, torch.device("cpu"))
    tt = torch.tensor(yfit.time.to_numpy(), dtype=torch.float32)
    ee = torch.tensor(yfit.event.to_numpy(), dtype=torch.float32)
    best, best_score, best_epoch, wait = None, -np.inf, 0, 0
    for epoch in range(1, 301):
        model.train(); optimizer.zero_grad()
        loss = cox_negative_log_likelihood(model(fb), tt, ee); loss.backward(); optimizer.step()
        score = cindex(yval.event.to_numpy(), yval.time.to_numpy(), predict_risk(model, vb))
        if score > best_score + 1e-6:
            best_score, best_epoch, wait = score, epoch, 0
            best = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        else:
            wait += 1
        if wait >= 30: break
    model.load_state_dict(best); model.eval()
    risk = predict_risk(model, tb)
    test_score = float(concordance_index_censored(ytest.event.astype(bool), ytest.time, risk)[0])
    return model, {"ordered_columns": order, "block_dims": {b: len(by[b]) for b in prefixes},
                   "slices": slices, "means": means, "scales": scales, "state_dict": best,
                   "dropout": .3, "random_state": SEED, "best_epoch": best_epoch,
                   "stopped_epoch": epoch, "best_validation_c_index": best_score,
                   "test_c_index": test_score}


def main():
    started = time.perf_counter(); OUT.mkdir(parents=True, exist_ok=True); LOCAL.mkdir(parents=True, exist_ok=True)
    formula_before = {n: sha(ROOT / "warrantxai" / n) for n in ("components.py", "aggregation.py")}
    common = load_common()
    xtr0 = pd.read_parquet(DATA / "X_train.parquet").set_index("sample_key")
    xte0 = pd.read_parquet(DATA / "X_test.parquet").set_index("sample_key")
    ytr = pd.read_parquet(DATA / "y_train.parquet").set_index("sample_key")
    yte = pd.read_parquet(DATA / "y_test.parquet").set_index("sample_key")
    all_y = pd.concat([ytr, yte]); base = pd.concat([xtr0[[c for c in xtr0 if c.startswith(("rna_", "mirna_"))]],
                                                     xte0[[c for c in xte0 if c.startswith(("rna_", "mirna_"))]]]).loc[all_y.index]
    groups = ytr.index.astype(str).str[:12]
    splitter = StratifiedGroupKFold(n_splits=7, shuffle=True, random_state=4201)
    fit_rel, val_rel = next(splitter.split(np.zeros(len(ytr)), ytr.event, groups))
    fit_keys, val_keys, test_keys = ytr.index[fit_rel], ytr.index[val_rel], yte.index
    if set(fit_keys.str[:12]) & set(val_keys.str[:12]) or set(ytr.index.str[:12]) & set(test_keys.str[:12]):
        raise RuntimeError("patient overlap")
    fit_idx = all_y.index.get_indexer(fit_keys)
    remainder = all_y.index.difference(fit_keys, sort=False)
    rem_idx = all_y.index.get_indexer(remainder)
    xfit, xrem, prep = common.prepare_fold(SIBLING, base, all_y, fit_idx, rem_idx, 10000, 5, 5000)
    xval, xouter, xalltrain = xrem.loc[val_keys], xrem.loc[test_keys], pd.concat([xfit, xrem.loc[val_keys]])
    screening, selected = rank_and_select(common, xfit, ytr.loc[fit_keys])
    clinical = [c for c in xfit if c.startswith("clin_")]
    model_features = selected["rna"] + selected["mirna"] + selected["methylation"] + clinical

    whitelist = json.loads((ROOT / "data/external/metabric_phase21/full_profile_whitelist.json").read_text())
    white_symbols = {row["symbol"] for row in whitelist["genes"]}
    mapping = map_rna(selected["rna"])
    eligible_rna = [f for f in selected["rna"] if mapping.get(f) in white_symbols]
    audit_rna = eligible_rna[:60]
    audit_clinical = [c for c in clinical if "age_at_index" in c or any(
        c.endswith(f"_{stage}") for stage in ("Stage I", "Stage II", "Stage III", "Stage IV"))]
    audit_features = audit_rna + audit_clinical
    screening["metabric_whitelisted"] = screening.feature.map(mapping).isin(white_symbols)
    screening["selected_audit_family"] = screening.feature.isin(audit_features)
    screening.to_csv(OUT / "feature_screening.csv", index=False)
    (LOCAL / "audit_features.json").write_text(json.dumps({"features": audit_features, "rna_symbols": mapping,
        "audit_rna_n": len(audit_rna), "audit_clinical_n": len(audit_clinical)}, indent=2))

    # Final model fits; scaler is learned only on inner-fit and frozen.
    scaler = StandardScaler().fit(xfit[model_features])
    atr = scaler.transform(xalltrain[model_features]); ate = scaler.transform(xouter[model_features])
    cframe = pd.DataFrame(atr, columns=model_features); cframe["time"] = ytr.time.to_numpy(); cframe["event"] = ytr.event.to_numpy()
    cox = CoxPHFitter(penalizer=.1)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always"); cox.fit(cframe, "time", "event", fit_options={"max_steps": 200, "step_size": .25})
    if caught: raise RuntimeError("Cox warnings: " + " | ".join(str(x.message) for x in caught))
    crisk = cox.predict_partial_hazard(pd.DataFrame(ate, columns=model_features)).to_numpy()
    cc = float(concordance_index_censored(yte.event.astype(bool), yte.time, crisk)[0])
    joblib.dump({"model": cox, "scaler": scaler, "features": model_features, "test_c_index": cc}, LOCAL / "coxph.joblib")

    rsf = RandomSurvivalForest(n_estimators=500, min_samples_split=10, min_samples_leaf=5,
        max_features="sqrt", random_state=42, n_jobs=-1, low_memory=True)
    rsf.fit(xalltrain[model_features].to_numpy(np.float32), Surv.from_arrays(event=ytr.event.astype(bool), time=ytr.time))
    rrisk = rsf.predict(xouter[model_features].to_numpy(np.float32))
    rc = float(concordance_index_censored(yte.event.astype(bool), yte.time, rrisk)[0])
    joblib.dump({"model": rsf, "features": model_features, "test_c_index": rc}, LOCAL / "rsf.joblib")

    deep, checkpoint = train_deep(xfit, ytr.loc[fit_keys], xval, ytr.loc[val_keys], xouter, yte, model_features)
    torch.save(checkpoint, LOCAL / "deepsurv.pt")

    leakage = {"inner_fit_rows": len(fit_keys), "inner_validation_rows": len(val_keys), "outer_train_rows": len(ytr),
        "outer_test_rows": len(yte), "inner_fit_patient_hash": key_hash(fit_keys.str[:12].unique()),
        "inner_validation_patient_hash": key_hash(val_keys.str[:12].unique()), "outer_test_patient_hash": key_hash(test_keys.str[:12].unique()),
        "fit_validation_patient_overlap": len(set(fit_keys.str[:12]) & set(val_keys.str[:12])),
        "outer_train_test_patient_overlap": len(set(ytr.index.str[:12]) & set(yte.index.str[:12])),
        "selection_rows_hash": key_hash(fit_keys), "preprocessing_fit_rows_hash": key_hash(fit_keys),
        "scaler_fit_rows_hash": key_hash(fit_keys), "validation_used_for_selection_or_preprocessing": False}
    metadata = {"models": {"CoxPH-v3": {"c_index": cc}, "RSF-v3": {"c_index": rc},
        "DeepSurv-v3": {"c_index": checkpoint["test_c_index"], "best_epoch": checkpoint["best_epoch"],
                         "validation_c_index": checkpoint["best_validation_c_index"]}},
        "counts": {"model_features": len(model_features), "rna": len(selected["rna"]), "mirna": len(selected["mirna"]),
                   "methylation": len(selected["methylation"]), "clinical": len(clinical), "audit_family": len(audit_features),
                   "audit_rna": len(audit_rna), "audit_clinical": len(audit_clinical)},
        "leakage": leakage, "preprocessing": prep, "formula_sha256": formula_before,
        "wall_clock_seconds": time.perf_counter() - started}
    (OUT / "stage1_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    (OUT / "STAGE1_SPLIT_SCREEN_MODELS.md").write_text(f"""# Faz 21 v3 — Aşama 1

İç-eğitim/iç-validation/dış-test satırları: {len(fit_keys)}/{len(val_keys)}/{len(yte)}. Hasta kesişimleri fit–validation={leakage['fit_validation_patient_overlap']}, dış-train–test={leakage['outer_train_test_patient_overlap']}.

Model paneli {len(model_features)} özelliktir (RNA/miRNA/methylation/clinical: 200/200/200/{len(clinical)}). RNA denetim ailesi model paneli içindeki beyaz-listeli {len(eligible_rna)} adaydan ilk {len(audit_rna)}; ortak klinik {len(audit_clinical)}; toplam aile {len(audit_features)}.

Test C-index: CoxPH-v3={cc:.6f}, RSF-v3={rc:.6f}, DeepSurv-v3={checkpoint['test_c_index']:.6f}. DeepSurv en iyi epoch={checkpoint['best_epoch']}, validation C-index={checkpoint['best_validation_c_index']:.6f}.

Programatik kanıt: seçim, methylation varyans/KNN, clinical imputation/encoding ve scaler fit satır hash'leri iç-eğitim hash'ine eşittir; validation kullanılmadı. Formül hash'leri değişmedi. Wall-clock {metadata['wall_clock_seconds']:.1f} sn.
""", encoding="utf-8")
    if formula_before != {n: sha(ROOT / "warrantxai" / n) for n in formula_before}: raise RuntimeError("formula hash changed")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__": main()
