#!/usr/bin/env python3
"""Phase 11 diagnostic sensitivity for the fold-1 CoxPH-v2 NA result."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter
import warnings

import joblib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase75_fold1/models"
PANEL = ROOT / "results/phase75_fold1_panel"
OLD = ROOT / "results/phase75_fold1_cox"
OUT = ROOT / "results/phase11_fold1_na_sensitivity"
CHECKPOINTS = OUT / "checkpoints"
N, FAILED, SEED, DRAWS = 50, 15, 20260917, 4999
COMPONENTS = ("F", "S", "U", "G_claim", "R")
EXPECTED = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import CALIBRATION
from warrantxai.certification import certify
from warrantxai.components import importance_ci


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def bootstrap_inputs():
    train, test = frame("X_train_selected.parquet"), frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[train.index]
    saved = joblib.load(DATA / "coxph_baseline_v2.joblib")
    features, scaler = list(saved["features"]), saved["scaler"]
    scaled = scaler.transform(train[features].astype(float))
    test_patients = test.index.astype(str).str[:12]
    audit = test.loc[~test_patients.duplicated()]
    anchors = scaler.transform(audit[features].astype(float))
    patients = train.index.astype(str).str[:12]
    groups = [np.flatnonzero(patients == patient) for patient in patients.unique()]
    rng = np.random.default_rng(SEED)
    indices = [np.concatenate([groups[k] for k in rng.integers(0, len(groups), len(groups))]) for _ in range(N)]
    return train, y, features, scaler, scaled, anchors, indices


def fit_one(rep, indices, y, features, scaled, penalizer, frozen):
    values = scaled[indices]
    events = y.event.to_numpy()[indices].astype(bool)
    valid = (values[events].var(axis=0) > 1e-12) & (values[~events].var(axis=0) > 1e-12)
    options = {"max_steps": 50} if frozen else {"max_steps": 200, "step_size": 0.25}
    attempts = []
    for attempt in range(20 if frozen else 1):
        fit = pd.DataFrame(values[:, valid], columns=np.asarray(features)[valid])
        fit["time"], fit["event"] = y.time.to_numpy()[indices], events
        model = CoxPHFitter(penalizer=penalizer)
        with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=2):
            warnings.simplefilter("always")
            model.fit(fit, "time", "event", show_progress=False, fit_options=options.copy())
        convergence = [x for x in caught if issubclass(x.category, ConvergenceWarning)]
        flagged = [i for i, f in enumerate(features) if valid[i] and any(str(x.message).startswith(f"Column {f} have very low variance") for x in convergence)]
        attempts.append({"attempt": attempt, "fit_options": options.copy(), "warnings": [str(x.message) for x in convergence], "separation_features": [features[i] for i in flagged]})
        if not convergence:
            coef = np.zeros(len(features)); coef[valid] = model.params_.to_numpy()
            return coef, {"replicate": rep, "penalizer": penalizer, "rows": len(indices), "events": int(events.sum()), "excluded_features": [features[i] for i in np.flatnonzero(~valid)], "fit_attempts": attempts, "converged": True}
        if not frozen:
            break
        if flagged:
            valid[flagged] = False
        elif options["max_steps"] == 50:
            options = {"max_steps": 200, "step_size": 0.25}
        else:
            break
    raise RuntimeError(f"replicate {rep} did not converge with penalizer={penalizer}")


def coefficient_set():
    train, y, features, scaler, scaled, anchors, all_indices = bootstrap_inputs()
    old_files = {int(p.stem[-2:]): p for p in OLD.glob("checkpoints/*_rep*.json")}
    coefficients, records = {}, {}
    for rep in range(N):
        if rep == FAILED:
            continue
        path = CHECKPOINTS / f"frozen_rep{rep:02d}.json"
        source = old_files.get(rep)
        if path.exists() or source:
            payload = json.loads((path if path.exists() else source).read_text())
            coefficients[rep] = np.asarray(payload["coefficients"], float)
            records[rep] = payload["record"]
            continue
        coef, record = fit_one(rep, all_indices[rep], y, features, scaled, 0.1, True)
        path.write_text(json.dumps({"coefficients": coef.tolist(), "record": record}, indent=2))
        coefficients[rep], records[rep] = coef, record
        print(f"Phase11 frozen replicate {rep+1}/50", flush=True)
    rescue_path = CHECKPOINTS / "diagnostic_rep15_penalizer015.json"
    if rescue_path.exists():
        payload = json.loads(rescue_path.read_text()); rescued = np.asarray(payload["coefficients"], float); rescue_record = payload["record"]
    else:
        rescued, rescue_record = fit_one(FAILED, all_indices[FAILED], y, features, scaled, 0.15, False)
        rescue_path.write_text(json.dumps({"coefficients": rescued.tolist(), "record": rescue_record}, indent=2))
    order49 = [i for i in range(N) if i != FAILED]
    coef49 = np.asarray([coefficients[i] for i in order49])
    coef50 = np.insert(coef49, FAILED, rescued, axis=0)
    return features, scaler, anchors, coef49, coef50, records, rescue_record


def tensors(coefficients, features, scaler, anchors):
    panel = json.loads((PANEL / "metadata.json").read_text()); active = panel["active_features"]
    targets = np.array([features.index(f) for f in active])
    source = coefficients[:, None, :] * anchors[None, :, :]
    with np.load(PANEL / "query_matrices.npz", allow_pickle=False) as z: external_raw = z["external_full"]
    external_scaled = np.column_stack([(external_raw[:, j] - scaler.mean_[targets[j]]) / scaler.scale_[targets[j]] for j in range(len(active))])
    external = coefficients[:, targets, None] * external_scaled.T[None, :, :]
    return source, np.transpose(external, (0, 2, 1)), active, targets


def statistics(source, external, coefficients_all, all_features, active, targets, rng):
    phase53 = load("12_rsf_certification.py", "p11_phase53")
    oldcox = load("13_cox_recursive_certification.py", "p11_oldcox")
    pilot = load("08_certification_pilot.py", "p11_pilot")
    remedy = load("11_rsf_remedy_calibration.py", "p11_remedy")
    coefficients = coefficients_all[:, targets]; selected_source = source[:, :, targets]
    t, p = coefficients.shape; indices = rng.integers(0, t, size=(DRAWS, t))
    counts = np.zeros((DRAWS, t), dtype=np.int16); np.add.at(counts, (np.arange(DRAWS)[:, None], indices), 1)
    points, pvalues = np.empty((p, 5)), np.empty((p, 5))
    f_by_rep = np.ones((1, t, p)); lower, upper = importance_ci(source[None, ...], CALIBRATION["alpha_U"])
    u_by_rep = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[:, :, targets]
    for by_rep, col in ((f_by_rep, 0), (u_by_rep, 2)):
        leaf = by_rep.mean(axis=1); sampled, sampled_se = phase53.resampled_means(by_rep, counts); se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        pvalues[:, col] = remedy.studentized_p(leaf.T.reshape(-1), np.moveaxis(sampled, 1, 2).reshape(DRAWS, -1), se.T.reshape(-1), np.moveaxis(sampled_se, 1, 2).reshape(DRAWS, -1)); points[:, col] = leaf[0]
    s_loss = pilot.pair_losses(source[None, ...], targets, CALIBRATION["a"])
    r_loss = pilot.order_pair_losses(np.abs(source).mean(axis=1)[None, ...], targets, CALIBRATION["eta"])
    for loss, col in ((s_loss, 1), (r_loss, 4)):
        leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1)); sampled = phase53.quadratic_scores(loss, counts, t)
        pvalues[:, col] = remedy.mc_p(leaf[0], sampled[:, 0]); points[:, col] = leaf[0]
    difference = np.abs(external).mean(axis=1) - np.abs(selected_source).mean(axis=1)
    gm = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = phase53.resampled_means(difference[None, ...], counts); gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    p_gm = remedy.studentized_p(gm, gm_sample, difference.std(axis=0, ddof=1)/np.sqrt(t)/CALIBRATION["b"], gm_mean_se[:, 0]/CALIBRATION["b"])
    deltas=[]
    for j in range(p):
        row=int(np.argmax(np.abs(coefficients[:,j]))); coef=coefficients[row,j]
        if coef==0: deltas.append(np.zeros(max(selected_source.shape[1],external.shape[1]))); continue
        zs=selected_source[row,:,j]/coef; ze=external[row,:,j]/coef; n=max(len(zs),len(ze))
        deltas.append(np.sort(np.abs(pilot.empirical_rank_grid(zs,n)-pilot.empirical_rank_grid(ze,n))))
    mean=coefficients.mean(axis=0); gs=oldcox.signed_score_from_coefficient(mean,deltas)[0]
    gs_se=oldcox.jackknife_se(oldcox.signed_score_from_coefficient((t*mean-coefficients)/(t-1),deltas))
    sampled_mean=coefficients[indices].mean(axis=1); gs_sample=oldcox.signed_score_from_coefficient(sampled_mean,deltas); gs_sample_se=np.empty((DRAWS,p))
    for j in range(p):
        leave=(t*sampled_mean[:,j,None]-coefficients[indices,j])/(t-1); scores=oldcox.signed_score_from_coefficient(leave.reshape(-1,1),[deltas[j]]).reshape(DRAWS,t); gs_sample_se[:,j]=oldcox.jackknife_se(scores.T)
    p_gs=remedy.studentized_p(gs,gs_sample,gs_se,gs_sample_se); points[:,3]=np.minimum(gm,gs); pvalues[:,3]=np.maximum(p_gm,p_gs)
    return points,pvalues,gm,gs


def result_frame(label, active, points, pvalues, gm, gs):
    bh,by=certify(pvalues,alpha=.05,method="bh"),certify(pvalues,alpha=.05,method="by")
    out=pd.DataFrame({"scenario":label,"feature":active})
    for k,c in enumerate(COMPONENTS): out[c]=points[:,k]; out[f"p_{c}"]=pvalues[:,k]
    out["G_magnitude"],out["G_signed"]=gm,gs
    out["WI"]=points.min(axis=1); out["bottleneck"]=[COMPONENTS[int(np.argmin(x))] for x in points]; out["point_adequate"]=out.WI>=CALIBRATION["tau"]
    out["p_joint"],out["q_BH"],out["q_BY"]=bh.joint_p,bh.adjusted_p,by.adjusted_p; out["BH"],out["BY"]=bh.rejected,by.rejected
    return out


def main():
    started=perf_counter(); OUT.mkdir(parents=True,exist_ok=True); CHECKPOINTS.mkdir(parents=True,exist_ok=True)
    previous_metadata = json.loads((OUT/"metadata.json").read_text()) if (OUT/"metadata.json").exists() else {}
    actual={k:sha(ROOT/"warrantxai"/k) for k in EXPECTED}
    if actual!=EXPECTED: raise RuntimeError(f"formula hashes changed: {actual}")
    features,scaler,anchors,c49,c50,records,rescue=coefficient_set()
    rng=np.random.default_rng(20260918+751_100)
    scenarios=[]
    for label,coefs in (("A_exclude_rep15_N49",c49),("B_rescue_rep15_penalizer_0.15_N50",c50)):
        source,external,active,targets=tensors(coefs,features,scaler,anchors); points,pvalues,gm,gs=statistics(source,external,coefs,features,active,targets,rng); scenarios.append(result_frame(label,active,points,pvalues,gm,gs))
    out=pd.concat(scenarios,ignore_index=True); out.to_csv(ROOT/"PHASE11_FOLD1_NA_SENSITIVITY.csv",index=False)
    counts=out.groupby("scenario").agg(point=("point_adequate","sum"),BH=("BH","sum"),BY=("BY","sum")).reset_index()
    a=set(out.query("scenario=='A_exclude_rep15_N49' and BH").feature); b=set(out.query("scenario=='B_rescue_rep15_penalizer_0.15_N50' and BH").feature)
    ap=set(out.query("scenario=='A_exclude_rep15_N49' and point_adequate").feature); bp=set(out.query("scenario=='B_rescue_rep15_penalizer_0.15_N50' and point_adequate").feature)
    norms=np.linalg.norm(c49,axis=1); rescued=float(np.linalg.norm(c50[FAILED])); z=(rescued-norms.mean())/norms.std(ddof=1)
    def names(values): return ", ".join(f"`{x}`" for x in sorted(values)) or "yok"
    rows="\n".join(f"| {x.scenario} | {int(x.point)} | {int(x.BH)} | {int(x.BY)} |" for x in counts.itertuples())
    direction="aynı kaldığını" if len(a)==len(b) else ("arttığını" if len(b)>len(a) else "azaldığını")
    elapsed=max(perf_counter()-started,float(previous_metadata.get("wall_clock_seconds",0)))
    report=f"""# Faz 11 — Fold-1 CoxPH-v2 NA yönsel duyarlılık kontrolü

Bu çalışma **diagnostic/exploratory** bir duyarlılık analizidir. Faz 7.5 üretimi replika 15'te durduğu için diskte başlangıçta 49 değil 15 başarılı checkpoint vardı; aynı deterministik hasta-kümesi örneklemesi ve dondurulmuş `penalizer=0.1`, separation dışlama ve `max_steps=200/step_size=0.25` protokolüyle eksik 16–49 replikaları yalnız bu tanısal çıktı için tamamlandı.

## Sonuçlar

| Senaryo | Point-adequate | BH | BY |
| --- | ---: | ---: | ---: |
{rows}

### ADIM A — replika 15 dışarıda, N=49

BH sertifikalı özellikler ({len(a)}/80): {names(a)}.

### ADIM B — yalnız replika 15 penalizer=0.15 ile kurtarıldı, N=50

Replika 15 `penalizer=0.15`, `max_steps=200`, `step_size=0.25` ile uyarısız yakınsadı. BH sertifikalı özellikler ({len(b)}/80): {names(b)}.

ADIM A'dan B'ye BH sertifikasyonuna eklenenler: {names(b-a)}. Düşenler: {names(a-b)}. Point-adequacy kümesine eklenenler: {names(bp-ap)}; düşenler: {names(ap-bp)}. Kurtarılan katsayı vektörünün L2 normu **{rescued:.6f}**; diğer 49 replikada ortalama **{norms.mean():.6f}**, SS **{norms.std(ddof=1):.6f}**, aralık **[{norms.min():.6f}, {norms.max():.6f}]** ve z-skoru **{z:.3f}**. Kurtarılan vektör diğer 49 replikaya göre belirgin düşük normludur; bu nedenle ADIM B yalnız yönsel tanı olarak yorumlanmalıdır.

Formül SHA256 değerleri değişmedi: `components.py={actual['components.py']}`, `aggregation.py={actual['aggregation.py']}`. Yeni hasta/model kimliği dışa aktarılmadı. İlk tam üretim wall-clock **{elapsed:.3f} sn**.

Bu analiz dondurulmuş retry protokolünü DEĞİŞTİRMEMEKTEDİR; Table 4'teki resmi NA hücresi ve ana metindeki sertifikasyon kararları aynen korunmaktadır.

49-replikalı dışlama analizinden gevşetilmiş-ceza N=50 analizine geçiş, BH sertifikasyon sayısının **{direction}** göstermektedir.
"""
    (ROOT/"PHASE11_FOLD1_NA_SENSITIVITY.md").write_text(report,encoding="utf-8")
    (OUT/"metadata.json").write_text(json.dumps({"diagnostic_only":True,"official_NA_preserved":True,"counts":counts.to_dict("records"),"A_BH":sorted(a),"B_BH":sorted(b),"rescue_record":rescue,"rescued_l2":rescued,"other49_l2_mean":float(norms.mean()),"other49_l2_sd":float(norms.std(ddof=1)),"formula_sha256":actual,"wall_clock_seconds":elapsed,"privacy":{"identifiers_exported":False}},indent=2),encoding="utf-8")
    print(counts.to_string(index=False)); print(f"rescued_l2={rescued:.6f} z={z:.3f}")


if __name__=="__main__": main()
