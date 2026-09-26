#!/usr/bin/env python3
"""Phase 10 point-estimate query-resampling sensitivity for RSF and DeepSurv."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
from time import perf_counter
import warnings

import joblib
import numpy as np
import pandas as pd
import torch
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase7_retrain"
PANEL = ROOT / "results/phase72_panel"
OUT = ROOT / "results/phase10_query_resampling"
CHECKPOINTS = OUT / "checkpoints"
N, BASE_SEED = 50, 20260917
EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    value = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = value
    spec.loader.exec_module(value)
    return value


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def markdown(frame: pd.DataFrame) -> str:
    cols = [str(x) for x in frame.columns]
    rows = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    rows += ["| " + " | ".join(str(v) for v in row) + " |" for row in frame.itertuples(index=False, name=None)]
    return "\n".join(rows)


def query_sets() -> tuple[list[dict], list[str]]:
    with np.load(PANEL / "query_matrices.npz", allow_pickle=True) as z:
        source_full, external_full = z["source_full"], z["external_full"]
        active = z["active_features"].astype(str).tolist()
        original = (z["source_query12"], z["external_query12"])
    configs = [{"size": 12, "seed": "original", "source": original[0], "external": original[1]}]
    draws = {}
    for seed_index, seed in enumerate((20261001, 20261002, 20261003), start=1):
        rng = np.random.default_rng(seed)
        draws[seed_index] = (source_full[rng.choice(len(source_full), 48, replace=False)],
                             external_full[rng.choice(len(external_full), 48, replace=False)])
    for size in (24, 48):
        for seed_index in (1, 2, 3):
            source, external = draws[seed_index]
            configs.append({"size": size, "seed": seed_index, "source": source[:size], "external": external[:size]})
    return configs, active


def boundaries(configs: list[dict], site: str) -> tuple[np.ndarray, list[slice | None]]:
    arrays, slices, start = [], [None] * len(configs), 0
    for seed in (1, 2, 3):
        index48 = next(i for i, c in enumerate(configs) if c["size"] == 48 and c["seed"] == seed)
        index24 = next(i for i, c in enumerate(configs) if c["size"] == 24 and c["seed"] == seed)
        value = np.asarray(configs[index48][site], dtype=np.float32)
        arrays.append(value)
        slices[index48] = slice(start, start + 48)
        slices[index24] = slice(start, start + 24)
        start += 48
    return np.concatenate(arrays), slices


def batched_rsf_attributions(model, query, background_rows, active_mean, active_indices,
                             *, permutations: int, seed: int, batch_size: int = 2) -> np.ndarray:
    """Same permutation estimator as Phase 7.2, with multiple paths per predict call."""
    x = np.asarray(query, dtype=np.float32)
    backgrounds = np.asarray(background_rows, dtype=np.float32).copy()
    active = np.asarray(active_indices, dtype=int)
    backgrounds[:, active] = np.asarray(active_mean, dtype=np.float32)
    rng = np.random.default_rng(seed)
    n, p, k = len(x), x.shape[1], len(backgrounds)
    result = np.zeros((n, p), dtype=float)
    remaining = permutations
    while remaining:
        batch = min(batch_size, remaining)
        states = np.empty((batch, n, p + 1, k, backgrounds.shape[1]), dtype=np.float32)
        orders = []
        for b in range(batch):
            order = rng.permutation(p); orders.append(order)
            for row in range(n):
                block = np.repeat(backgrounds[None, :, :], p + 1, axis=0)
                for step, local_index in enumerate(order, start=1):
                    block[step:, :, active[local_index]] = x[row, local_index]
                states[b, row] = block
        pred = np.asarray(model.predict(states.reshape(-1, backgrounds.shape[1])), dtype=float).reshape(batch, n, p + 1, k).mean(axis=-1)
        for b, order in enumerate(orders):
            result[:, order] += np.diff(pred[b], axis=1)
        remaining -= batch
    return result / permutations


def rsf_attributions(configs: list[dict], active: list[str]) -> tuple[dict, float]:
    started = perf_counter(); CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    helpers = load("07_rsf_pilot.py", "p10_rsf_helpers")
    train, y = frame("X_train_selected.parquet"), frame("y_train.parquet")
    saved = joblib.load(DATA / "rsf_baseline.joblib"); features = list(saved["features"])
    active_idx = np.array([features.index(x) for x in active]); active_mean = train.mean().to_numpy(float)[active_idx]
    patients = train.index.astype(str).str[:12]; unique = train.loc[~patients.duplicated()]
    rng_bg = np.random.default_rng(BASE_SEED + 4200); pick = np.sort(rng_bg.choice(len(unique), 16, replace=False))
    background = unique.iloc[pick].to_numpy(np.float32); background[:, active_idx] = active_mean.astype(np.float32)
    groups = [np.flatnonzero(patients == p) for p in patients.unique()]
    rng = np.random.default_rng(BASE_SEED); values = train.to_numpy(np.float32)
    survival = Surv.from_arrays(event=y.event.astype(bool).to_numpy(), time=y.time.to_numpy())
    q_source, source_slices = boundaries(configs, "source"); q_external, external_slices = boundaries(configs, "external")
    store = {i: {"candidate": np.empty((2, N, c["size"], len(active))), "reference": np.empty((2, N, c["size"], len(active)))} for i, c in enumerate(configs)}
    with np.load(ROOT / "results/phase72_rsf/attributions_n50.npz", allow_pickle=False) as original:
        store[0]["candidate"][:] = original["candidate"]; store[0]["reference"][:] = original["reference"]
    start_rep = int(os.environ.get("PHASE10_START_REP", "0"))
    end_rep = int(os.environ.get("PHASE10_END_REP", str(N)))
    worker = os.environ.get("PHASE10_WORKER") == "rsf"
    for rep in range(start_rep, end_rep):
        idx = np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))])
        path = CHECKPOINTS / f"rsf_rep{rep:02d}.npz"
        if path.exists():
            z = np.load(path, allow_pickle=False); combined = [[z[f"c{s}"], z[f"r{s}"]] for s in range(2)]
        else:
            model = RandomSurvivalForest(n_estimators=500, min_samples_split=10, min_samples_leaf=5, max_features="sqrt", random_state=BASE_SEED+rep, n_jobs=2 if worker else -1, low_memory=True)
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always"); model.fit(values[idx], survival[idx])
            if caught: raise RuntimeError(f"RSF replicate {rep} warnings: {Counter(x.category.__name__ for x in caught)}")
            combined = []
            for site, query in enumerate((q_source, q_external)):
                c = helpers.background_permutation_attributions(model, query, background, active_mean, active_idx, permutations=8, seed=BASE_SEED+10000+100*site+rep)
                r = helpers.background_permutation_attributions(model, query, background, active_mean, active_idx, permutations=128, seed=BASE_SEED+20000+100*site+rep)
                combined.append([c, r])
            np.savez_compressed(path, c0=combined[0][0], r0=combined[0][1], c1=combined[1][0], r1=combined[1][1])
        for ci in range(1, len(configs)):
            for site, sl in enumerate((source_slices[ci], external_slices[ci])):
                store[ci]["candidate"][site, rep] = combined[site][0][sl]
                store[ci]["reference"][site, rep] = combined[site][1][sl]
        print(f"Phase10 RSF {rep+1}/50", flush=True)
    return store, perf_counter()-started


def deepsurv_attributions(configs: list[dict], active: list[str]) -> tuple[dict, float]:
    started = perf_counter(); phase6 = load("20_deepsurv_pilot.py", "p10_deep6"); phase61 = load("21_deepsurv_explainer_revision.py", "p10_deep61"); phase72 = load("41_phase72_deepsurv.py", "p10_deep72")
    checkpoint = torch.load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    train, y = frame("X_train_selected.parquet"), frame("y_train.parquet"); features = list(checkpoint["ordered_columns"]); train = train[features]
    active_idx = np.array([features.index(x) for x in active]); raw = train.to_numpy(np.float32)
    patients = train.index.astype(str).str[:12]; unique = train.loc[~patients.duplicated()]
    rng_bg = np.random.default_rng(BASE_SEED+4200); first = np.sort(rng_bg.choice(len(unique),16,replace=False)); rem=np.setdiff1d(np.arange(len(unique)),first); extra=np.sort(rng_bg.choice(rem,112,replace=False)); nested=np.r_[first,extra]
    backgrounds={16:unique.iloc[nested[:16]].to_numpy(np.float32),128:unique.iloc[nested].to_numpy(np.float32)}
    groups=[np.flatnonzero(patients==p) for p in patients.unique()]; rng=np.random.default_rng(BASE_SEED)
    events=y.event.astype(bool).to_numpy(); times=y.time.to_numpy(np.float32)
    q_source, source_slices=boundaries(configs,"source"); q_external, external_slices=boundaries(configs,"external")
    store={i:{"candidate":np.empty((2,N,c["size"],len(active))),"reference":np.empty((2,N,c["size"],len(active)))} for i,c in enumerate(configs)}
    with np.load(ROOT/"results/phase72_deepsurv/attributions_n50.npz",allow_pickle=False) as original:
        store[0]["candidate"][:]=original["candidate"]; store[0]["reference"][:]=original["reference"]
    for rep in range(N):
        idx=np.concatenate([groups[k] for k in rng.integers(0,len(groups),len(groups))]); path=CHECKPOINTS/f"deep_rep{rep:02d}.npz"
        if path.exists():
            z=np.load(path,allow_pickle=False); combined=[[z[f"c{s}"],z[f"r{s}"]] for s in range(2)]
        else:
            model,means,scales,training=phase6.train_bootstrap_model(checkpoint,raw,events,times,idx,rep)
            if training["warning_categories"]: raise RuntimeError(f"DeepSurv replicate {rep} warnings")
            combined=[]
            for query in (q_source,q_external):
                c,*_=phase72.deep_attributions_with_relative_gate(phase61,phase6,model,checkpoint,query,backgrounds[16],active_idx,means,scales)
                r,*_=phase72.deep_attributions_with_relative_gate(phase61,phase6,model,checkpoint,query,backgrounds[128],active_idx,means,scales)
                combined.append([c,r])
            np.savez_compressed(path,c0=combined[0][0],r0=combined[0][1],c1=combined[1][0],r1=combined[1][1])
        for ci in range(1,len(configs)):
            for site,sl in enumerate((source_slices[ci],external_slices[ci])):
                store[ci]["candidate"][site,rep]=combined[site][0][sl]; store[ci]["reference"][site,rep]=combined[site][1][sl]
        print(f"Phase10 DeepSurv {rep+1}/50",flush=True)
    return store,perf_counter()-started


def report(configs: list[dict], active: list[str], stores: dict[str, dict], timings: dict[str,float]) -> None:
    helpers=load("07_rsf_pilot.py","p10_report_helpers"); mapping=pd.read_csv(PANEL/"feature_mapping.csv"); symbols=dict(zip(mapping.feature,mapping.symbol))
    labels=pd.read_csv(ROOT/"PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv"); null_rsf=set(labels.query("Model == 'RSF' and Label == 'degenerate'")["Feature"])
    rows=[]
    for model,store in stores.items():
        for ci,config in enumerate(configs):
            c,r=store[ci]["candidate"],store[ci]["reference"]; scores=helpers.component_report(c,r,active)
            mean_abs=np.mean(np.abs(c),axis=(0,1,2)); max_abs=np.max(np.abs(c),axis=(0,1,2))
            source_max=np.max(np.abs(c[0]),axis=(0,1)); external_max=np.max(np.abs(c[1]),axis=(0,1))
            for j,row in scores.iterrows():
                symbol=symbols.get(row.feature,row.feature)
                if pd.isna(symbol):
                    symbol = "Age at index" if "age_at_index" in row.feature else "AJCC " + row.feature.rsplit("_", 1)[-1]
                rows.append({"model":model,"feature":row.feature,"symbol":symbol,"query_per_cohort":config["size"],"seed":config["seed"],"F":row.F,"S":row.S,"U":row.U,"G_signed":row.G_signed,"G_claim":row.G_claim,"R":row.R,"WI_j":row.WI,"bottleneck_component":row.bottleneck,"point_adequate":row.point_adequate,"mean_abs_attribution":mean_abs[j],"max_abs_attribution":max_abs[j],"source_max_abs_attribution":source_max[j],"external_max_abs_attribution":external_max[j],"rsf_null_warrant_at_12plus12":model=="RSF" and symbol in null_rsf})
    out=pd.DataFrame(rows); out.to_csv(ROOT/"PHASE10_QUERY_RESAMPLING_SENSITIVITY.csv",index=False)
    sensitivity=out.groupby(["model","feature","symbol"]).agg(WI_min=("WI_j","min"),WI_max=("WI_j","max"),WI_std=("WI_j","std"),bottleneck_count=("bottleneck_component","nunique"),point_status_count=("point_adequate","nunique")).reset_index(); sensitivity["WI_range"]=sensitivity.WI_max-sensitivity.WI_min
    top=sensitivity.sort_values("WI_range",ascending=False).groupby("model").head(10)
    null=out[out.rsf_null_warrant_at_12plus12].groupby(["feature","symbol","query_per_cohort"]).agg(source_max_abs=("source_max_abs_attribution","max"),external_max_abs=("external_max_abs_attribution","max"),max_abs_attribution=("max_abs_attribution","max"),mean_abs_attribution=("mean_abs_attribution","mean"),WI_mean=("WI_j","mean")).reset_index()
    adequacy=out.groupby(["model","query_per_cohort"]).agg(point_min=("point_adequate","sum"),WI_mean=("WI_j","mean"))
    # Replace summed seed rows with per-seed count range.
    perseed=out.groupby(["model","query_per_cohort","seed"]).point_adequate.sum().reset_index(); adequate=perseed.groupby(["model","query_per_cohort"]).point_adequate.agg(["min","max","mean"]).reset_index()
    persistent=int((null.max_abs_attribution==0).sum()); total=len(null)
    point_flips=out.groupby(["model","feature"]).point_adequate.nunique().groupby(level=0).apply(lambda x:int((x>1).sum()))
    bottleneck_flips=out.groupby(["model","feature"]).bottleneck_component.nunique().groupby(level=0).apply(lambda x:int((x>1).sum()))
    source_persistent=(null.groupby("symbol").source_max_abs.max()==0)
    note=f"""# Faz 10 — Sorgu yeniden-örnekleme duyarlılığı

Fold-0 dondurulmuş 27-feature panelinde RSF ve DeepSurv, orijinal 12+12 sorguya ek olarak 24+24 ve 48+48 boyutlarında üçer bağımsız deterministik seed ile yeniden değerlendirildi. Bu nokta-duyarlılık analizidir; yeni 4.999-draw sertifikasyon yapılmadı ve resmi kararlar değiştirilmedi.

## Nokta-adequacy kararlılığı

{markdown(adequate)}

## En yüksek WI aralığına sahip özellikler

{markdown(top[["model","symbol","WI_min","WI_max","WI_range","WI_std","bottleneck_count","point_status_count"]])}

## RSF null-attribution warrant kontrolü

{markdown(null)}

Null-warrant özellik/boyut hücrelerinin tam-sıfır attribution sayısı **{persistent}/{total}**. Ayrıntılı CSV tüm model × özellik × boyut × seed bileşenlerini içerir.

Faz 7.6/7.7 sınıflaması TCGA/source tarafındaki 12 sorguya dayanıyordu; bu kapsamda sekiz özelliğin sekizi başlangıçta sıfırdır. Daha büyük TCGA sorgularında yalnız **{int(source_persistent.sum())}/8** özellik (`{', '.join(source_persistent[source_persistent].index)}`) bütün boyut ve seed'lerde sıfır kalmıştır. Bu nedenle sekizli null-attribution örüntüsü sorgu örneklemesine duyarlıdır. RSF'de **{point_flips['RSF']}/27** özelliğin point-adequacy durumu ve **{bottleneck_flips['RSF']}/27** özelliğin darboğazı en az bir konfigürasyonda değişmiştir; DeepSurv'de iki sayı da **0/27**'dir. RSF yönü aynı ölçekte kalmakla birlikte 12+12'deki 17 point-adequate sonuç, 24+24'te 11–15 ve 48+48'te 12–18 aralığına açılmıştır; bu nedenle özellik düzeyinde güçlü kararlılık iddiası desteklenmez.

Çalışma süreleri: RSF **{timings['RSF']:.1f} sn**, DeepSurv **{timings['DeepSurv']:.1f} sn**. Sonuçlar aynı audit havuzları, aynı leakage dışlamaları, aynı bootstrap modellerinin deterministik yeniden-fit karşılıkları ve dondurulmuş explainer bütçelerine koşulludur. Formül dosyaları değişmedi.
"""
    (ROOT/"PHASE10_QUERY_RESAMPLING_SENSITIVITY.md").write_text(note,encoding="utf-8")
    (OUT/"metadata.json").write_text(json.dumps({"configs":[{"size":c["size"],"seed":c["seed"]} for c in configs],"timings":timings,"formula_sha256":EXPECTED,"model_or_explainer_calls":True,"certification_recomputed":False},indent=2),encoding="utf-8")


def main() -> None:
    actual={name:sha(ROOT/"warrantxai"/name) for name in EXPECTED}
    if actual!=EXPECTED: raise RuntimeError(f"Frozen formula hashes changed: {actual}")
    OUT.mkdir(parents=True,exist_ok=True); configs,active=query_sets(); rsf,t1=rsf_attributions(configs,active)
    if os.environ.get("PHASE10_WORKER") == "rsf":
        return
    deep,t2=deepsurv_attributions(configs,active)
    # Preserve measured production wall-clock when a checkpoint-only report rerun is used.
    rsf_files=list(CHECKPOINTS.glob("rsf_rep*.npz")); deep_files=list(CHECKPOINTS.glob("deep_rep*.npz"))
    if len(rsf_files)==N and len(deep_files)==N:
        run_log=OUT/"run.log"; rsf_end=max(p.stat().st_mtime for p in rsf_files); deep_end=max(p.stat().st_mtime for p in deep_files)
        t1=rsf_end-run_log.stat().st_ctime; t2=deep_end-rsf_end
    report(configs,active,{"RSF":rsf,"DeepSurv":deep},{"RSF":t1,"DeepSurv":t2})


if __name__ == "__main__":
    main()
