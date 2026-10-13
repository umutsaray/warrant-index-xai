#!/usr/bin/env python3
"""Phase 13: diagnostic RSF tree traversal and extreme-coalition effects."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import perf_counter
import sys

import joblib
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")
SEED = 20260917
FORMULAS = ("components.py", "aggregation.py", "certification.py")
SPECIAL = {"AJCC Stage IV", "OR6C1", "OR4C46", "OR13C3", "KRTAP19-2",
           "KRTAP19-7", "OR5M10", "RBMY1F", "ACCSL", "AJCC Stage I",
           "AJCC Stage II", "AJCC Stage III"}


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def md(df: pd.DataFrame) -> str:
    shown=df.copy()
    for c in shown.select_dtypes(include=[np.number]):
        shown[c]=shown[c].map(lambda x: "NA" if pd.isna(x) else f"{x:.6g}")
    return "\n".join(["| "+" | ".join(shown.columns)+" |",
        "| "+" | ".join(["---"]*len(shown.columns))+" |",
        *("| "+" | ".join(map(str,r))+" |" for r in shown.itertuples(index=False,name=None))])


def load_fold(fold: int):
    if fold == 0:
        data=ROOT/"results/phase7_retrain"; panel=ROOT/"results/phase72_panel/query_matrices.npz"
        attrs=ROOT/"results/phase72_rsf/attributions_n50.npz"
    else:
        data=ROOT/"results/phase75_fold1/models"; panel=ROOT/"results/phase75_fold1_panel/query_matrices.npz"
        attrs=ROOT/"results/phase75_fold1_rsf/attributions_n50.npz"
    saved=joblib.load(data/"rsf_baseline.joblib"); model=saved["model"]; features=list(saved["features"])
    frozen=(model.n_estimators==500 and model.min_samples_split==10 and model.min_samples_leaf==5
            and model.max_features=="sqrt" and model.random_state==42)
    if not frozen: raise RuntimeError(f"fold {fold} frozen RSF hyperparameters changed")
    train=pd.read_parquet(data/"X_train_selected.parquet").set_index("sample_key")[features]
    with np.load(panel,allow_pickle=False) as z:
        active=z["active_features"].astype(str).tolist(); queries=[z["source_query12"],z["external_query12"]]
    with np.load(attrs,allow_pickle=False) as z: candidate=z["candidate"]
    active_idx=np.array([features.index(f) for f in active]); active_mean=train.mean().to_numpy()[active_idx]
    patients=train.index.astype(str).str[:12]; unique=train.loc[~patients.duplicated()]
    rng=np.random.default_rng(SEED+4200); selection=np.sort(rng.choice(len(unique),16,replace=False))
    background=unique.iloc[selection].to_numpy(np.float32); background[:,active_idx]=active_mean.astype(np.float32)
    return dict(fold=fold,model=model,features=features,active=active,active_idx=active_idx,
                active_mean=active_mean,queries=queries,candidate=candidate,background=background)


def split_inventory(bundle):
    rows=[]
    for local,(name,global_j) in enumerate(zip(bundle["active"],bundle["active_idx"],strict=True)):
        thresholds=[]; trees=0
        for tree in bundle["model"].estimators_:
            hit=tree.tree_.feature==global_j
            if hit.any(): trees+=1; thresholds.extend(tree.tree_.threshold[hit].tolist())
        rows.append((name,local,global_j,np.asarray(thresholds),trees))
    return rows


def paths_for_matrix(model, matrix, global_j):
    counts=np.zeros(len(matrix),dtype=int)
    for tree in model.estimators_:
        nodes=tree.decision_path(matrix)
        split_nodes=np.flatnonzero(tree.tree_.feature==global_j)
        if len(split_nodes): counts += (nodes[:,split_nodes].getnnz(axis=1)>0)
    return counts


def analyze(bundle, display, labels, split_reference):
    fold=bundle["fold"]; step1=[]; step2=[]; step3=[]; summary=[]
    inventory=split_inventory(bundle)
    allowed=set(bundle["active"]) if fold==0 else {
        f for f in bundle["active"] if display.get(f,f) in {"AJCC Stage IV","OR13C3"}}
    for feature,local,global_j,thresholds,tree_count in inventory:
        if feature not in allowed: continue
        if fold==0 and feature in split_reference.index:
            expected=int(split_reference.loc[feature,"Trees with ≥1 split"])
            if tree_count!=expected: raise RuntimeError(f"split-tree mismatch {feature}: {tree_count}!={expected}")
        base=float(bundle["active_mean"][local]); total=len(thresholds)
        aggregate={"fold":fold,"feature_id":feature,"feature":display.get(feature,feature),
            "split_nodes":total,"split_trees":tree_count,"threshold_min":float(thresholds.min()) if total else np.nan,
            "threshold_median":float(np.median(thresholds)) if total else np.nan,
            "threshold_max":float(thresholds.max()) if total else np.nan,"base_value":base,
            "warrant_type":labels.get(feature,"NA")}
        for site,cohort in enumerate(("TCGA","METABRIC")):
            q=bundle["queries"][site][:,local]
            separated=[]
            for qi,value in enumerate(q):
                # sklearn/sksurv Tree follows left when x_j <= threshold.
                count=int(np.sum((value<=thresholds)!=(base<=thresholds))) if total else 0
                ratio=count/total if total else 0.0
                separated.append(ratio)
                step1.append({"fold":fold,"feature_id":feature,"feature":display.get(feature,feature),
                    "cohort":cohort,"query":qi,"query_value":value,"base_value":base,
                    "split_nodes":total,"separating_nodes":count,"separating_ratio":ratio,
                    "thresholds":json.dumps(thresholds.tolist())})
            full=np.repeat(bundle["background"],len(q),axis=0)
            base_matrix=full.copy()
            for qi in range(len(q)):
                start,end=qi*16,(qi+1)*16; ix=np.ix_(np.arange(start,end),bundle["active_idx"])
                full[ix]=bundle["queries"][site][qi]
                base_matrix[ix]=bundle["active_mean"]
            visits_full=paths_for_matrix(bundle["model"],full,global_j).reshape(len(q),16)
            visits_base=paths_for_matrix(bundle["model"],base_matrix,global_j).reshape(len(q),16)
            all_minus=full.copy(); all_minus[:,global_j]=base
            base_plus=base_matrix.copy(); base_plus[:,global_j]=np.repeat(q,16)
            pred_full=bundle["model"].predict(full); pred_minus=bundle["model"].predict(all_minus)
            pred_plus=bundle["model"].predict(base_plus); pred_base=bundle["model"].predict(base_matrix)
            dl=(pred_full-pred_minus).reshape(len(q),16); df=(pred_plus-pred_base).reshape(len(q),16)
            for qi in range(len(q)):
                for bi in range(16):
                    step2 += [
                        {"fold":fold,"feature_id":feature,"feature":display.get(feature,feature),"cohort":cohort,
                         "query":qi,"background_row":bi,"state":"all_active_query","visited_trees":visits_full[qi,bi],"visit_rate":visits_full[qi,bi]/500},
                        {"fold":fold,"feature_id":feature,"feature":display.get(feature,feature),"cohort":cohort,
                         "query":qi,"background_row":bi,"state":"all_active_base","visited_trees":visits_base[qi,bi],"visit_rate":visits_base[qi,bi]/500}]
                    step3.append({"fold":fold,"feature_id":feature,"feature":display.get(feature,feature),"cohort":cohort,
                        "query":qi,"background_row":bi,"query_value":q[qi],"base_value":base,
                        "delta_last":dl[qi,bi],"abs_delta_last":abs(dl[qi,bi]),
                        "delta_first":df[qi,bi],"abs_delta_first":abs(df[qi,bi])})
            aggregate.update({
                f"{cohort}_M_j":float(np.abs(bundle["candidate"][site,:,:,local]).mean()),
                f"{cohort}_query_unique":int(len(np.unique(q))),
                f"{cohort}_separating_min":float(np.min(separated)),f"{cohort}_separating_median":float(np.median(separated)),f"{cohort}_separating_max":float(np.max(separated)),
                f"{cohort}_visit_full_mean":float(visits_full.mean()/500),f"{cohort}_visit_base_mean":float(visits_base.mean()/500),
                f"{cohort}_max_abs_delta_last":float(np.abs(dl).max()),f"{cohort}_nonzero_delta_last":int((np.abs(dl)>1e-12).sum()),
                f"{cohort}_max_abs_delta_first":float(np.abs(df).max()),f"{cohort}_nonzero_delta_first":int((np.abs(df)>1e-12).sum())})
        summary.append(aggregate)
    return map(pd.DataFrame,(step1,step2,step3,summary))


def main():
    started=perf_counter(); before={f:sha(ROOT/"warrantxai"/f) for f in FORMULAS}
    table=pd.read_csv(ROOT/"PHASE77_FULL_FEATURE_TABLE_FOLD0.csv")
    rsf=table[table.Model=="RSF"]; display=dict(zip(rsf["Feature ID"],rsf.Feature)); labels=dict(zip(rsf["Feature ID"],rsf["Informative/Degenerate label"]))
    split=pd.read_csv(ROOT/"results/phase76_rsf_split_usage.csv").set_index("Feature ID")
    outputs=[analyze(load_fold(0),display,labels,split)]
    fold1=load_fold(1); outputs.append(analyze(fold1,display,{},split))
    # map objects are materialized once here.
    f0=list(outputs[0]); f1=list(outputs[1])
    combined=[pd.concat([f0[i],f1[i]],ignore_index=True) for i in range(3)]+[f0[3]]
    names=("PHASE13_ADIM1_SPLIT_SIDES.csv","PHASE13_ADIM2_PATH_VISITS.csv","PHASE13_ADIM3_EXTREME_COALITION.csv","PHASE13_FEATURE_SUMMARY.csv")
    for name,frame in zip(names,combined,strict=True): frame.to_csv(ROOT/name,index=False)
    summary=combined[3]; main=summary; fold1_summary=f1[3]; zero=main[main.TCGA_M_j==0]
    same_side=(zero.TCGA_separating_max==0); no_visit=(zero.TCGA_visit_full_mean==0)&(zero.TCGA_visit_base_mean==0)
    zero_ext=(zero.TCGA_nonzero_delta_last==0)&(zero.TCGA_nonzero_delta_first==0)
    h1=same_side|no_visit
    special=summary[summary.feature.isin(SPECIAL)]
    cols=["feature","warrant_type","TCGA_M_j","METABRIC_M_j","split_trees","split_nodes","base_value",
          "TCGA_query_unique","METABRIC_query_unique","TCGA_separating_max","METABRIC_separating_max",
          "TCGA_visit_full_mean","METABRIC_visit_full_mean","TCGA_max_abs_delta_last","TCGA_max_abs_delta_first",
          "METABRIC_max_abs_delta_last","METABRIC_max_abs_delta_first"]
    after={f:sha(ROOT/"warrantxai"/f) for f in FORMULAS}; assert before==after
    repro=("ATLANDI: Saklı checkpoint'ler yalnız candidate/reference dizilerini içeriyor; bunları üreten bootstrap "
           "RSF model nesnesi saklanmamış. Yasaklanan yeniden fit yapılmadan aynı seed ile tekrar üretim mümkün değildir.")
    text=f"""# Faz 13 — RSF ağaç gezinme ve marjinal etki testi

**Statü:** diagnostic/exploratory. Komut: `python experiments/61_phase13_rsf_traversal.py`. Wall-clock **{perf_counter()-started:.3f} sn**.

## ADIM 0 — Kurulum

Fold-0 ve fold-1 modellerinde 500 ağaç, `min_samples_split=10`, `min_samples_leaf=5`, `max_features=sqrt`, `random_state=42` doğrulandı. Açıklanan skaler çıktı ve bu analizde kullanılan çıktı, `RandomSurvivalForest.predict()` risk skorudur. Ağaç yönü sklearn/sksurv uygulamasındaki `x_j <= threshold` sol dal kuralıdır.

Tek-replika attribution sağlaması: **{repro}** Fold-0 baseline üzerinde yeni attribution üretmek saklı bootstrap tensörüyle aynı modeli karşılaştırmayacağı için yapılmadı.

## ADIM 1–3 — Sonuç

Fold-0'da Tablo 8'de raporlanan 12 aktif özelliğin ağaç split sayıları birebir eşleşti; kalan 15 aktif özellik için sayımlar doğrudan ormandan çıkarıldı. Ham sorgu/split oranları, her background satırındaki yol ziyaretleri ve iki uç koalisyonun risk-skoru farkları ayrı CSV'lerdedir.

Özel özellikler:

{md(special[cols])}

TCGA M_j=0 olan **{len(zero)}** RSF özelliğinin **{int(same_side.sum())}** tanesinde sorgu ve baz bütün split eşiklerinin aynı tarafında, **{int(no_visit.sum())}** tanesinde iki yol durumunda ziyaret oranı sıfır, **{int(zero_ext.sum())}** tanesinde iki uç koalisyonun tüm 12×16 hücrelerinde marjinal etki sıfırdır. H1'in aynı-taraf veya hiç-ziyaret koşulu **{int(h1.sum())}/{len(zero)}** özellikte desteklenmiştir. Dışarıda kalanlar: `{', '.join(zero.loc[~h1,'feature']) or 'yok'}`. Uç-koalisyon testi tam Shapley hesabı değildir.

Fold-1 baseline saklı olduğu için iki analog da aynı analizden geçirildi:

{md(fold1_summary[cols])}

Fold-1 sonuçları ham ADIM 1–3 CSV'lerinde `fold=1` ile işaretlidir; 27 satırlık `PHASE13_FEATURE_SUMMARY.csv` yalnız resmi fold-0 panelini içerir.

## Bütünlük

Başlangıç SHA256: `{json.dumps(before)}`. Bitiş SHA256: `{json.dumps(after)}`. Değerler aynıdır; yeni model fit, bootstrap veya attribution tensörü üretilmedi.

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
"""
    (ROOT/"PHASE13_RSF_TRAVERSAL.md").write_text(text,encoding="utf-8")
    print(text)

if __name__=="__main__": main()
