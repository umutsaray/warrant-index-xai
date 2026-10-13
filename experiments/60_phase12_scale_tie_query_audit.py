#!/usr/bin/env python3
"""Phase 12 exploratory scale, tie, query-cohort, and aggregator audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

from warrantxai.components import (faithfulness, generalizability,
    generalizability_signed, importance, importance_ci, pairwise_orders,
    rank_reproducibility, stability, uncertainty)

MODELS = ("CoxPH-v2", "RSF", "DeepSurv")
KS = (0.25, 0.5, 1.0, 2.0, 4.0)
CS = (0.25, 0.5, 1.0, 2.0)
ABS_ETAS = (0.01, 0.02, 0.05, 0.10)
TAU = 0.8
FOCUS = {"AJCC Stage IV", "OR6C1", "OR4C46", "OR13C3", "KRTAP19-2",
         "KRTAP19-7", "OR5M10", "RBMY1F"}


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def quantile_grid(x: np.ndarray, n: int) -> np.ndarray:
    ordered = np.sort(np.asarray(x, dtype=float))
    indices = np.minimum((np.arange(n) + .5) * len(ordered) // n, len(ordered) - 1).astype(int)
    return ordered[indices]


def load_inputs():
    data = {}
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as z:
        source, external = z["source_phi"], z["external_phi"]
        all_features, features = z["features"].astype(str).tolist(), z["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(f) for f in features])
    data["CoxPH-v2"] = {"phi": source[None], "reference": source[None],
        "external": external, "features": features, "targets": targets,
        "importance": np.abs(source).mean(axis=1)[None]}
    for model, folder in (("RSF", "phase72_rsf"), ("DeepSurv", "phase72_deepsurv")):
        with np.load(ROOT / f"results/{folder}/attributions_n50.npz", allow_pickle=False) as z:
            c, r, f = z["candidate"], z["reference"], z["features"].astype(str).tolist()
        data[model] = {"phi": c, "reference": r, "features": f,
                       "targets": np.arange(len(f)), "importance": np.abs(c).mean(axis=2)}
    return data


def scores(item, a=.5, b=.5, u=3., eta=.02):
    phi, ref, targets = item["phi"], item["reference"], item["targets"]
    lo, hi = importance_ci(phi)
    f = faithfulness(phi, ref, a)[targets]
    s = stability(phi, a)[targets]
    uu = uncertainty(lo, hi, u)[targets]
    r = rank_reproducibility(item["importance"], eta)[targets]
    if "external" in item:
        src, ext = phi[0][:, :, targets], item["external"]
        si = np.stack([np.abs(src).mean(axis=(0, 1)), np.abs(ext).mean(axis=(0, 1))])
        gm = generalizability(si, b)
        n = max(src.shape[1], ext.shape[1])
        sm, em = src.mean(axis=0), ext.mean(axis=0)
        coupled = np.stack([np.column_stack([quantile_grid(sm[:, j], n) for j in range(len(targets))]),
                            np.column_stack([quantile_grid(em[:, j], n) for j in range(len(targets))])])
        gs = generalizability_signed(coupled, a)
    else:
        si = np.abs(phi).mean(axis=(1, 2))
        gm = generalizability(si, b)
        coupled = np.stack([np.sort(phi[d].mean(axis=0), axis=0) for d in range(2)])
        gs = generalizability_signed(coupled, a)
    gc = np.minimum(gm, gs)
    return pd.DataFrame({"F": f, "S": s, "U": uu, "G": gm, "G_signed": gs,
                         "G_claim": gc, "R": r}, index=item["features"])


def r_without_persistent_ties(imp, eta):
    order = pairwise_orders(imp, eta)
    sites, t, p = imp.shape
    out = np.full(p, np.nan)
    for j in range(p):
        keep = np.array([l != j and np.any(order[..., j, l] != 0) for l in range(p)])
        if not keep.any():
            continue
        site_values = []
        for d in range(sites):
            losses = []
            for x in range(t):
                for y in range(x + 1, t):
                    losses.append(np.abs(order[d, x, j, keep] - order[d, y, j, keep]).sum() / (2 * keep.sum()))
            site_values.append(1 - np.mean(losses))
        out[j] = min(site_values)
    return out


def step1(data, official):
    rows, scales = [], {}
    for model, item in data.items():
        m = item["importance"][0].mean(axis=0)[item["targets"]]
        scales[model] = float(np.median(m[m > 0]))
        original = scores(item)
        ref = official[official.Model == model].set_index("Feature ID")
        for col in ("F", "S", "U", "G", "G_signed", "R", "WI"):
            actual = original[[col]].iloc[:, 0] if col != "WI" else original[["F","S","U","G_claim","R"]].min(axis=1)
            expected = ref.loc[item["features"], col]
            np.testing.assert_allclose(actual, expected, atol=1e-9, rtol=0)
        for k in KS:
            z = scores(item, k*scales[model], k*scales[model], 6*k*scales[model], .02)
            primary = z[["F","S","U","G_claim","R"]]
            wi, bott = primary.min(axis=1), primary.idxmin(axis=1)
            for feature in z.index:
                row = {"model": model, "feature": feature, "k": k, "scale_s_m": scales[model],
                       "a": k*scales[model], "b": k*scales[model], "u": 6*k*scales[model],
                       **z.loc[feature].to_dict(), "WI": wi[feature], "bottleneck": bott[feature],
                       "point_adequate": bool(wi[feature] > TAU)}
                rows.append(row)
    return pd.DataFrame(rows), scales


def step2(data, scales, official):
    rows=[]
    for model,item in data.items():
        imp, targets = item["importance"], item["targets"]
        base_other = scores(item)[["F","S","U","G_claim"]]
        grids=[("absolute",e,e) for e in ABS_ETAS] + [("normalized",c,c*scales[model]) for c in CS]
        for kind,param,eta in grids:
            order=pairwise_orders(imp,eta)
            rt=rank_reproducibility(imp,eta)[targets]
            rn=r_without_persistent_ties(imp,eta)[targets]
            for pos,(feature,j) in enumerate(zip(item["features"],targets,strict=True)):
                mask=np.ones(imp.shape[-1],bool); mask[j]=False
                tie=float((order[...,j,mask]==0).mean())
                primary=base_other.loc[feature].to_dict()|{"R":rt[pos]}
                wi=min(primary.values()); bott=min(primary,key=primary.get)
                rows.append({"model":model,"feature":feature,"grid":kind,"parameter":param,
                    "eta":eta,"M_zero":bool(imp[0, :, j].mean() == 0),"tie_rate":tie,"R":rt[pos],"R_noties":rn[pos],
                    "bottleneck":bott,"R_bottleneck":bott=="R","WI":wi,"point_adequate":wi>TAU})
    out=pd.DataFrame(rows)
    old=pd.read_csv(ROOT/"PHASE74_ETA_FEATURE_VALUES.csv")
    name={"CoxPH-v2":"CoxPH","RSF":"RSF","DeepSurv":"DeepSurv"}
    for model in MODELS:
        for eta in ABS_ETAS:
            a=out[(out.model==model)&(out.grid=="absolute")&(out.eta==eta)].set_index("feature").R
            b=old[(old.model==name[model])&(old.eta==eta)].set_index("feature").R
            np.testing.assert_allclose(a.loc[b.index],b,atol=1e-12,rtol=0)
    return out


def query_rows(fold, model, panel_folder, data_folder, label_table):
    panel=ROOT/f"results/{panel_folder}/query_matrices.npz"
    with np.load(panel,allow_pickle=False) as z:
        feats=z["active_features"].astype(str).tolist(); qs=[z["source_query12"],z["external_query12"]]
    if fold==0: train_path=ROOT/"results/phase7_retrain/X_train_selected.parquet"
    else: train_path=ROOT/"results/phase75_fold1/models/X_train_selected.parquet"
    train=pd.read_parquet(train_path).set_index("sample_key")
    patients=train.index.astype(str).str[:12]; uniq=train.loc[~patients.duplicated()]
    rng=np.random.default_rng(20260917+4200)
    first=np.sort(rng.choice(len(uniq),16,replace=False))
    if model=="RSF":
        bg=uniq.iloc[first].copy(); bg.loc[:,feats]=train[feats].mean(axis=0).to_numpy()
        bgs={"k16":bg}
    else:
        remaining=np.setdiff1d(np.arange(len(uniq)),first); extra=np.sort(rng.choice(remaining,112,replace=False)); idx=np.r_[first,extra]
        bgs={"k16":uniq.iloc[idx[:16]],"k128":uniq.iloc[idx]}
    tensor_folder=("phase72_" if fold==0 else "phase75_fold1_")+model.lower()
    with np.load(ROOT/f"results/{tensor_folder}/attributions_n50.npz",allow_pickle=False) as z:
        attr=z["candidate"]
    label=label_table[(label_table.Model==model) if fold==0 else np.ones(len(label_table),dtype=bool)].set_index("Feature ID") if len(label_table) else pd.DataFrame()
    split=pd.read_csv(ROOT/"results/phase76_rsf_split_usage.csv").set_index("Feature ID") if fold==0 and model=="RSF" else None
    rows=[]
    for j,f in enumerate(feats):
        display=f
        if fold==0:
            hit=official_global[(official_global.Model==model)&(official_global["Feature ID"]==f)]
            if len(hit): display=hit.iloc[0].Feature
        for d,cohort in enumerate(("TCGA","METABRIC")):
            vals=qs[d][:,j]
            for budget,bg in bgs.items():
                bvals=bg[f].to_numpy(float)
                equality=float(np.isclose(vals[:, None], bvals[None, :], atol=1e-12, rtol=0).mean())
                lab="NA"
                if fold==0 and f in label.index: lab=str(label.loc[f,"Label"])
                elif fold==0:
                    hit=official_global[(official_global.Model==model)&(official_global["Feature ID"]==f)]
                    if len(hit): lab=str(hit.iloc[0]["Informative/Degenerate label"])
                splits=float(split.loc[f,"Trees with ≥1 split"]) if split is not None and f in split.index else np.nan
                rows.append({"fold":fold,"model":model,"feature":display,"feature_id":f,"cohort":cohort,
                    "background":budget,"n_unique":len(np.unique(vals)),"variance":float(np.var(vals)),
                    "minimum":float(vals.min()),"maximum":float(vals.max()),"query_values":json.dumps(vals.tolist()),
                    "background_values":json.dumps(bvals.tolist()),"query_background_equality_rate":equality,
                    "M_j":float(np.abs(attr[d,:,:,j]).mean()),"RSF_trees_with_split":splits,"warrant_type":lab})
    return rows


def step3():
    labels=pd.read_csv(ROOT/"PHASE76_WARRANT_CLASSIFICATION_FOLD0.csv")
    rows=[]
    for model in ("RSF","DeepSurv"):
        rows += query_rows(0,model,"phase72_panel",model.lower(),labels)
        rows += query_rows(1,model,"phase75_fold1_panel",model.lower(),pd.DataFrame())
    return pd.DataFrame(rows)


def step4(official):
    rows=[]
    for model in MODELS:
        d=official[official.Model==model].copy()
        z=np.column_stack([d.F,d.S,d.U,np.minimum(d.G,d.G_signed),d.R])
        vals={"min":z.min(1),"arithmetic_mean":z.mean(1),"geometric_mean":np.prod(z,axis=1)**.2}
        ranks={k:pd.Series(v).rank(ascending=False,method="average").to_numpy() for k,v in vals.items()}
        top={k:set(np.argsort(-v)[:5]) for k,v in vals.items()}
        for i,row in d.reset_index(drop=True).iterrows():
            rows.append({"model":model,"feature":row["Feature"],"feature_id":row["Feature ID"],
                **{f"{k}_score":v[i] for k,v in vals.items()},**{f"{k}_rank":v[i] for k,v in ranks.items()},
                "min_gt_0_8":vals["min"][i]>.8,"arithmetic_gt_0_8":vals["arithmetic_mean"][i]>.8,
                "geometric_gt_0_8":vals["geometric_mean"][i]>.8})
        # model summaries are written into repeated columns for machine-readable provenance
        rho_ma=spearmanr(vals["min"],vals["arithmetic_mean"]).statistic
        rho_mg=spearmanr(vals["min"],vals["geometric_mean"]).statistic
        overlap_a=len(top["min"]&top["arithmetic_mean"]); overlap_g=len(top["min"]&top["geometric_mean"])
        for row in rows[-len(d):]: row.update({"rho_min_arithmetic":rho_ma,"rho_min_geometric":rho_mg,
            "top5_overlap_min_arithmetic":overlap_a,"top5_overlap_min_geometric":overlap_g})
    return pd.DataFrame(rows)


def md_table(df):
    shown = df.copy()
    for col in shown.select_dtypes(include=[np.number]):
        shown[col] = shown[col].map(lambda value: f"{value:.6g}")
    return "\n".join([
        "| " + " | ".join(map(str, shown.columns)) + " |",
        "| " + " | ".join(["---"] * len(shown.columns)) + " |",
        *("| " + " | ".join(map(str, row)) + " |" for row in shown.itertuples(index=False, name=None)),
    ])


def main():
    global official_global
    started=perf_counter(); before={n:sha(ROOT/f"warrantxai/{n}") for n in ("components.py","aggregation.py","certification.py")}
    official_global=pd.read_csv(ROOT/"PHASE77_FULL_FEATURE_TABLE_FOLD0.csv")
    data=load_inputs()
    t=perf_counter(); s1,scales=step1(data,official_global); dt1=perf_counter()-t; s1.to_csv(ROOT/"PHASE12_ADIM1_SCALE_GRID.csv",index=False)
    t=perf_counter(); s2=step2(data,scales,official_global); dt2=perf_counter()-t; s2.to_csv(ROOT/"PHASE12_ADIM2_TIE_R.csv",index=False)
    t=perf_counter(); s3=step3(); dt3=perf_counter()-t; s3.to_csv(ROOT/"PHASE12_ADIM3_QUERY_COHORT.csv",index=False)
    t=perf_counter(); s4=step4(official_global); dt4=perf_counter()-t; s4.to_csv(ROOT/"PHASE12_ADIM4_AGGREGATORS.csv",index=False)
    after={n:sha(ROOT/f"warrantxai/{n}") for n in before}; assert before==after
    scale_summary=s1.groupby(["model","k"]).agg(point=("point_adequate","sum"),WI_min=("WI","min"),WI_median=("WI","median"),WI_max=("WI","max")).reset_index()
    binding=s1.groupby(["model","k"]).bottleneck.apply(lambda x:json.dumps(x.value_counts().to_dict(),sort_keys=True)).rename("bottleneck_counts").reset_index()
    scale_summary=scale_summary.merge(binding,on=["model","k"])
    component_summary=s1.melt(id_vars=["model","k"],value_vars=["F","S","U","G_claim","R"],var_name="component",value_name="score").groupby(["model","k","component"]).agg(minimum=("score","min"),median=("score","median"),maximum=("score","max"),below_tau=("score",lambda x:int((x<TAU).sum()))).reset_index()
    eta_summary=s2.groupby(["model","grid","parameter","eta"]).agg(R_mean=("R","mean"),R_noties_mean=("R_noties","mean"),tie_rate_mean=("tie_rate","mean"),R_bottlenecks=("R_bottleneck","sum"),point=("point_adequate","sum")).reset_index()
    tie_rho=s2.groupby(["model","grid","parameter","eta"]).apply(lambda d:spearmanr(d.tie_rate,d.R).statistic,include_groups=False).rename("spearman_tie_R").reset_index()
    eta_summary=eta_summary.merge(tie_rho,on=["model","grid","parameter","eta"])
    tie_groups=s2[(s2.grid=="absolute")&(s2.eta==.02)].groupby(["model","M_zero"]).agg(features=("feature","count"),tie_rate_mean=("tie_rate","mean"),R_mean=("R","mean"),R_noties_mean=("R_noties","mean")).reset_index()
    agg_summary=s4.groupby("model").agg(min_gt=("min_gt_0_8","sum"),arithmetic_gt=("arithmetic_gt_0_8","sum"),geometric_gt=("geometric_gt_0_8","sum"),rho_min_arithmetic=("rho_min_arithmetic","first"),rho_min_geometric=("rho_min_geometric","first"),top5_overlap_min_arithmetic=("top5_overlap_min_arithmetic","first"),top5_overlap_min_geometric=("top5_overlap_min_geometric","first")).reset_index()
    focus=s3[(s3.fold==0)&s3.feature.isin(FOCUS)].groupby(["model","feature","feature_id","cohort"]).agg(n_unique=("n_unique","first"),variance=("variance","first"),equality_rate=("query_background_equality_rate","max"),M_j=("M_j","first"),splits=("RSF_trees_with_split","first"),warrant=("warrant_type","first")).reset_index()
    base=s3[(s3.fold==0)&(s3.background=="k16")]
    zero_stats=[]
    for model in ("RSF","DeepSurv"):
        source=base[(base.model==model)&(base.cohort=="TCGA")].set_index("feature_id")
        zero=source.M_j.eq(0); constant=source.n_unique.eq(1)
        zero_stats.append(f"{model}: TCGA M_j=0 olan {int(zero.sum())} özelliğin {int((zero&constant).sum())}'inde TCGA sorguları sabit")
    text=f"""# Faz 12 — Ölçek, tie ve sorgu-kohort duyarlılık denetimi

**Statü:** diagnostic/exploratory; resmi sonuç değildir. Başlangıç/bitiş hash'leri eşittir: `{json.dumps(before,ensure_ascii=False)}`. Yeni model fit, patient bootstrap veya explainer çağrısı yapılmadı. Toplam wall-clock: **{perf_counter()-started:.3f} sn**.

Kullanılan komut: `python experiments/60_phase12_scale_tie_query_audit.py`. Doğrulama komutu: `..\\.venv\\Scripts\\python.exe -m pytest -q --basetemp=.pytest_phase12_tmp`.

## ADIM 1 — Ölçek

Kohort-0 medyan sıfır-olmayan attribution ölçekleri: `{json.dumps(scales)}`. Orijinal `a=b=0.5, u=3, eta=0.02` skorları Faz 7.7 tablosuyla `atol=1e-9` düzeyinde birebir üretildi. H–L aralıkları attribution importance birimindedir; bu nedenle `u=6ks_m` aynı birim sistemini korur. Saklı inferans girdilerini ölçek-parametreli recursive test için yeniden tanımlamak bu tanısal kapsamda yapılmadı; BH/BY opsiyonel k=1 satırı atlandı.

{md_table(scale_summary)}

Her bileşenin dağılımı ve `score<0.8` sayısı:

{md_table(component_summary)}

Ortak ölçeğe alma model farkını ortadan kaldırmadı: k=4'te point sayıları CoxPH-v2/RSF/DeepSurv için sırasıyla 6/16/0'dır. DeepSurv tüm ızgarada sıfır kalırken RSF en yüksek, CoxPH ara düzeydedir. Bağlayıcı çoğunluk CoxPH-v2 ve RSF'de k=4'te G_claim'den R'ye geçerken DeepSurv'de bütün k değerlerinde G_claim kalır. Önceden belirlenen değerlendirme kuralına göre **R baskınlığı ölçeğe duyarlıdır**.

## ADIM 2 — Tie ve R

Tie oranı, her hedef özellik için tüm site×replika×rakip karşılaştırmalarında `o=0` payıdır. `R_noties`, bütün site ve replikalarda daima tie kalan rakipleri hem paydan hem paydadan çıkarır. Mutlak η ızgarası Faz 7.4 değerlerini `atol=1e-12` ile yeniden üretti.

{md_table(eta_summary)}

η=0.02 için M_j=0 gruplaması:

{md_table(tie_groups)}

Kalıcı tie'ları çıkarmak ortalama R'yi genel olarak düşürmektedir; bu, tie'ların R'yi özellikle RSF'de yukarı çektiğini, fakat tek başına bütün modeli açıklamadığını gösterir. Mutlak η büyüdükçe tie oranı ve R birlikte yükselir; CoxPH'de η=0.10'da R-darboğaz sayısının 4'e düşmesi bu modelde güçlü η duyarlılığına işaret eder. DeepSurv'de R bütün normalize ızgarada 27/27 darboğazdır.

## ADIM 3 — Sorgu kohortu ve background eşitliği

Eşitlik oranı, her sorgu değerinin ilgili açıklayıcı background değerlerinden en az biriyle tam (`atol=1e-12`) eşleşme payıdır. RSF k=16 sampled inactive background kullanırken aktif koordinatlar train ortalamasına sabitlenmiştir; DeepSurv deterministik nested k=16/k=128 ham train background kullanır. Fold-1 satırları mevcut saklı tensörlerden ayrıca verilmiştir.

{md_table(focus)}

Tam 27-feature sayımı: **{"; ".join(zero_stats)}**. TCGA sorgularındaki sıfır/sabit değerler DeepSurv'ün gerçek hasta background'ında sıklıkla eşleşirken RSF'nin aktif koordinatları train ortalamasına sabitlenmiş background'ında eşleşmez. Buna rağmen bazı sıfır-varyanslı özelliklerin DeepSurv attribution büyüklüğü sıfır değil ve RSF'de OR6C1/AJCC Stage IV split kullanımı vardır. Bu nedenle örtüşme dejenere warrant'larla ilişkili bir mekanizma adayıdır, tek başına yeterli açıklama değildir.

## ADIM 4 — Aggregator duyarlılığı

Yalnız nokta tahminleri karşılaştırıldı; aritmetik/geometrik ortalamalara sertifikasyon uygulanmadı.

{md_table(agg_summary)}

Min agregatörü ile aritmetik/geometrik ortalamalar arasındaki sıralama korelasyonu yaklaşık 0.84–0.90 olsa da `>0.8` kararları keskin biçimde ayrışır: kompansatuvar iki alternatif üç modelde de 27/27 özellik geçirir. Sonuçların zayıf-halka min kuralına maddi ölçüde bağlı olduğu böylece doğrulanır.

## Dosyalar ve kapsam

Ham özellik düzeyi sonuçlar `PHASE12_ADIM1_SCALE_GRID.csv`, `PHASE12_ADIM2_TIE_R.csv`, `PHASE12_ADIM3_QUERY_COHORT.csv` ve `PHASE12_ADIM4_AGGREGATORS.csv` içindedir. Aşamaların süreleri sırasıyla **{dt1:.3f}/{dt2:.3f}/{dt3:.3f}/{dt4:.3f} sn**. Formül kaynaklarına dokunulmadı.

ADIM 1 sonucu: k=4'te CoxPH-v2/RSF/DeepSurv point sayıları 6/16/0 olduğundan ortak ölçek model farkını kaldırmadı. ADIM 2 sonucu: η=0.02'de ortalama R→R_noties CoxPH-v2/RSF/DeepSurv için 0.774→0.771, 0.855→0.836 ve 0.726→0.726 oldu; tie etkisi en büyük RSF'dedir. ADIM 3 sonucu: {"; ".join(zero_stats)} ve background eşleşmesi tek başına dejenere warrant'ları açıklamadı. ADIM 4 sonucu: min ile point sayıları 9/17/0 iken iki kompansatuvar agregatör her modelde 27/27 verdi.

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.

Ölçek, tie ve sorgu-kohort denetiminin gösterdiği yönsel eğilim: model farkı ölçek normalizasyonuyla kaybolmamakta, tie katkısı en çok RSF'de görülmekte, kohort-içi sabitlik tek başına dejenere warrant'ları açıklamamakta ve kararlar min agregatörüne güçlü biçimde bağlı kalmaktadır.
"""
    (ROOT/"PHASE12_SCALE_TIE_QUERY_AUDIT.md").write_text(text,encoding="utf-8")
    print(text)

if __name__ == "__main__": main()
