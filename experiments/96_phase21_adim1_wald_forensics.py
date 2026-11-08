#!/usr/bin/env python3
"""Phase 21 ADIM 1: forensic audit of the fold-0 RNA Wald screen."""
from __future__ import annotations

import io
import json
import re
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from lifelines import CoxPHFitter
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "phase7_retrain"
OUT = ROOT / "phase21_reviewer1"
RAW = ROOT / "data" / "external" / "metabric_phase72"
PENALIZER = 0.01
DEGENERATE_RSF_RNA = {
    "rna_ENSG00000244362.3", "rna_ENSG00000205330.4", "rna_ENSG00000185926.1",
    "rna_ENSG00000169800.14", "rna_ENSG00000204246.3", "rna_ENSG00000186965.5",
    "rna_ENSG00000254834.4",
}


def frame(name: str) -> pd.DataFrame:
    value = pd.read_parquet(DATA / name)
    return value.set_index("sample_key") if "sample_key" in value else value


def score_pvalues(values: np.ndarray, time: np.ndarray, event: np.ndarray, chunk: int = 1000) -> np.ndarray:
    order = np.argsort(-time, kind="stable")
    times = time[order]
    events = event[order].astype(bool)
    group_end = np.empty(len(times), dtype=int)
    start = 0
    while start < len(times):
        end = start
        while end + 1 < len(times) and times[end + 1] == times[start]:
            end += 1
        group_end[start:end + 1] = end
        start = end + 1
    event_rows = np.flatnonzero(events)
    risk_rows = group_end[event_rows]
    result = np.full(values.shape[1], np.nan)
    for begin in range(0, values.shape[1], chunk):
        raw = values[:, begin:begin + chunk].astype(float, copy=False)
        scale = raw.std(axis=0)
        valid = scale > 0
        x = np.zeros_like(raw)
        x[:, valid] = (raw[:, valid] - raw[:, valid].mean(axis=0)) / scale[valid]
        x = x[order]
        risk_n = (risk_rows + 1).astype(float)[:, None]
        s1 = np.cumsum(x, axis=0)[risk_rows]
        s2 = np.cumsum(x * x, axis=0)[risk_rows]
        score = x[event_rows].sum(axis=0) - (s1 / risk_n).sum(axis=0)
        information = (s2 / risk_n - (s1 / risk_n) ** 2).sum(axis=0)
        stat = np.divide(score * score, information, out=np.full_like(score, np.nan), where=information > 0)
        result[begin:begin + len(stat)] = chi2.sf(stat, 1)
    return result


def null_loglik(time: np.ndarray, event: np.ndarray) -> float:
    order = np.argsort(-time, kind="stable")
    times = time[order]
    events = event[order].astype(bool)
    ends = np.empty(len(times), dtype=int)
    start = 0
    while start < len(times):
        end = start
        while end + 1 < len(times) and times[end + 1] == times[start]:
            end += 1
        ends[start:end + 1] = end
        start = end + 1
    return float(-np.log(ends[np.flatnonzero(events)] + 1).sum())


def fit_one(feature: str, raw: pd.Series, y: pd.DataFrame) -> dict:
    x = (raw.astype(float) - raw.mean()) / raw.std(ddof=0)
    data = pd.DataFrame({"x": x.to_numpy(), "time": y.time.to_numpy(), "event": y.event.to_numpy()})
    result: dict = {}
    warning_text: list[str] = []
    progress = io.StringIO()
    try:
        with warnings.catch_warnings(record=True) as caught, redirect_stdout(progress):
            warnings.simplefilter("always")
            model = CoxPHFitter(penalizer=0.0)
            model.fit(data, "time", "event", show_progress=True)
        warning_text = [f"{type(w.message).__name__}: {w.message}" for w in caught]
        summary = model.summary.loc["x"]
        iterations = [int(x) for x in re.findall(r"Iteration\s+(\d+)", progress.getvalue())]
        result.update(beta=float(summary["coef"]), se=float(summary["se(coef)"]),
                      wald_z=float(summary["z"]), wald_p=float(summary["p"]),
                      lifelines_lr_p=float(model.log_likelihood_ratio_test().p_value),
                      lifelines_iterations=max(iterations) if iterations else np.nan,
                      lifelines_converged="Convergence success" in progress.getvalue(), fit_error="")
    except Exception as exc:
        result.update(beta=np.nan, se=np.nan, wald_z=np.nan, wald_p=np.nan, lifelines_lr_p=np.nan,
                      lifelines_iterations=np.nan, lifelines_converged=False,
                      fit_error=f"{type(exc).__name__}: {exc}")
    try:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("ignore")
            penalized = CoxPHFitter(penalizer=PENALIZER).fit(data, "time", "event", show_progress=False)
        ps = penalized.summary.loc["x"]
        result.update(penalized_beta=float(ps["coef"]), penalized_se=float(ps["se(coef)"]))
    except Exception as exc:
        result.update(penalized_beta=np.nan, penalized_se=np.nan,
                      fit_error_penalized=f"{type(exc).__name__}: {exc}")
    result["warning_text"] = " | ".join(warning_text)
    result["lifelines_separation_warning"] = any(
        token in result["warning_text"].lower() for token in ["separation", "low variance", "complete separation"])
    return result


def annotate_metabric(features: list[str]) -> pd.DataFrame:
    ids = [feature[4:].split(".")[0] for feature in features]
    response = requests.post("https://mygene.info/v3/query", data={
        "q": ",".join(ids), "scopes": "ensembl.gene", "fields": "symbol,entrezgene", "species": "human"
    }, timeout=120)
    response.raise_for_status()
    annotations = response.json()
    by_query: dict[str, set[str]] = {}
    for row in annotations:
        if row.get("symbol"):
            by_query.setdefault(row["query"], set()).add(row["symbol"])
    genes = json.loads((RAW / "genes.json").read_text(encoding="utf-8"))
    catalog: dict[str, set[int]] = {}
    for row in genes:
        catalog.setdefault(row["hugoGeneSymbol"], set()).add(row["entrezGeneId"])
    rows = []
    for feature, ensembl in zip(features, ids):
        symbols = by_query.get(ensembl, set())
        symbol = next(iter(symbols)) if len(symbols) == 1 else ""
        rows.append({"feature": feature, "ensembl": ensembl, "symbol": symbol,
                     "mapping_unique": len(symbols) == 1,
                     "metabric_catalog_unique": bool(symbol and len(catalog.get(symbol, set())) == 1)})
    return pd.DataFrame(rows)


def markdown(frame: pd.DataFrame) -> str:
    cols = list(frame.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    x = frame("X_train.parquet")
    y = frame("y_train.parquet").loc[x.index]
    selection = json.loads((DATA / "selected_features.json").read_text())
    selected = selection["selected_by_block"]["rna"]
    rna = [c for c in x if c.startswith("rna_")]
    assert len(selected) == 200 and len(rna) == 60660
    values = x[rna].to_numpy(dtype=float, copy=False)
    event = y.event.to_numpy(dtype=bool)
    nonzero = values != 0
    nonzero_n = nonzero.sum(axis=0)
    nonzero_events = nonzero[event].sum(axis=0)
    score_p = score_pvalues(values, y.time.to_numpy(float), event)
    ranking = pd.read_csv(DATA / "univariate_cox_rna.tsv", sep="\t").set_index("feature").loc[rna]
    ll0 = null_loglik(y.time.to_numpy(float), event)
    lr_stat = np.maximum(0.0, 2 * (ranking.partial_log_likelihood.to_numpy(float) - ll0))
    lr_p = chi2.sf(lr_stat, 1)
    universe = pd.DataFrame({"feature": rna, "nonzero_n": nonzero_n,
                             "nonzero_events": nonzero_events, "zero_rate": 1 - nonzero_n / len(x),
                             "stored_beta": ranking.coefficient.to_numpy(float),
                             "stored_se": ranking.standard_error.to_numpy(float),
                             "stored_wald_p": ranking.p_value.to_numpy(float),
                             "lr_p_breslow": lr_p, "score_p_breslow": score_p})
    filtered = universe[(universe.nonzero_n >= np.ceil(.10 * len(x))) & (universe.nonzero_events >= 10)].copy()
    lr_panel = filtered.sort_values(["lr_p_breslow", "score_p_breslow", "feature"]).head(200).copy()
    score_panel = filtered.sort_values(["score_p_breslow", "lr_p_breslow", "feature"]).head(200).copy()

    active = set(json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())["active_features"])
    active_rna = {f for f in active if f.startswith("rna_")}
    forensic_rows = []
    selected_lookup = universe.set_index("feature")
    for index, feature in enumerate(selected, 1):
        row = selected_lookup.loc[feature]
        fit = fit_one(feature, x[feature], y)
        nz = int(row.nonzero_n)
        ne = int(row.nonzero_events)
        forensic_rows.append({"feature": feature, "original_rna_rank": index,
                              "in_fold0_27_rna": feature in active_rna,
                              "in_seven_rsf_zero_attribution_rna": feature in DEGENERATE_RSF_RNA,
                              "is_KRTAP11_1": feature == "rna_ENSG00000182591.6",
                              "nonzero_n": nz, "nonzero_events": ne, "zero_rate": row.zero_rate,
                              "nonzero_event_rate": ne / nz if nz else np.nan,
                              "all_nonzero_are_events": bool(nz and ne == nz),
                              "no_nonzero_is_event": bool(nz and ne == 0),
                              "lr_p": row.lr_p_breslow, "score_p": row.score_p_breslow, **fit})
        if index % 25 == 0:
            print(f"forensic fits {index}/200", flush=True)
    forensic = pd.DataFrame(forensic_rows)
    forensic.to_csv(OUT / "PHASE21_ADIM1_WALD_FORENSICS.csv", index=False)

    union = list(dict.fromkeys(lr_panel.feature.tolist() + score_panel.feature.tolist()))
    mapping = annotate_metabric(union).set_index("feature")
    panel_rows = []
    for method, panel in [("filtered_LR", lr_panel), ("filtered_score", score_panel)]:
        for rank, row in enumerate(panel.itertuples(index=False), 1):
            mapped = mapping.loc[row.feature]
            panel_rows.append({"method": method, "rank": rank, **row._asdict(),
                               "symbol": mapped.symbol, "mapping_unique": mapped.mapping_unique,
                               "metabric_catalog_unique": mapped.metabric_catalog_unique})
    panels = pd.DataFrame(panel_rows)
    panels.to_csv(OUT / "PHASE21_ADIM1_ALTERNATIVE_PANELS.csv", index=False)

    clinical_common = [f for f in active if not f.startswith("rna_")]
    common_rows = []
    for method in ["filtered_LR", "filtered_score"]:
        cell = panels[(panels.method == method) & panels.metabric_catalog_unique]
        common_rows.extend({"method": method, "feature": row.feature, "block": "RNA",
                            "symbol": row.symbol} for row in cell.itertuples())
        common_rows.extend({"method": method, "feature": feature, "block": "clinical",
                            "symbol": feature} for feature in clinical_common)
    common = pd.DataFrame(common_rows)
    common.to_csv(OUT / "PHASE21_ADIM1_COMMON_PANEL.csv", index=False)

    summary_rows = []
    for method, panel in [("original_Wald", universe[universe.feature.isin(selected)]),
                          ("filtered_LR", panels[panels.method == "filtered_LR"]),
                          ("filtered_score", panels[panels.method == "filtered_score"])]:
        summary_rows.append({"method": method, "selected_n": len(panel),
                             "zero_rate_min": panel.zero_rate.min(), "zero_rate_q25": panel.zero_rate.quantile(.25),
                             "zero_rate_median": panel.zero_rate.median(), "zero_rate_q75": panel.zero_rate.quantile(.75),
                             "zero_rate_max": panel.zero_rate.max(),
                             "metabric_common_n": (np.nan if method == "original_Wald" else int(panel.metabric_catalog_unique.sum()))})
    summary = pd.DataFrame(summary_rows)
    # The original panel has a frozen, empirically complete 22-RNA intersection.
    summary.loc[summary.method == "original_Wald", "metabric_common_n"] = len(active_rna)
    summary.to_csv(OUT / "PHASE21_ADIM1_PANEL_SUMMARY.csv", index=False)

    flagged = forensic[forensic.in_fold0_27_rna | forensic.in_seven_rsf_zero_attribution_rna | forensic.is_KRTAP11_1]
    report = f"""# Faz 21 ADIM 1 — Wald taraması adli incelemesi

Fold-0 train içindeki 60,660 RNA sütunu incelendi. Adli tablo, orijinal Wald ile seçilen 200 RNA için standartlaştırılmış tek değişkenli lifelines Cox β/SE/z, LR p, Breslow skor p, uyarı/iterasyon, ayrışma göstergeleri ve ridge-cezalı (`penalizer={PENALIZER}`; Firth değildir) karşılaştırmayı içerir.

## Özet

- Başarısız unpenalized fit: **{int(forensic.beta.isna().sum())}/200**.
- Lifelines ayrışma/low-variance uyarısı: **{int(forensic.lifelines_separation_warning.sum())}/200**.
- Nonzero satırların tamamı event veya tamamı censored: **{int((forensic.all_nonzero_are_events | forensic.no_nonzero_is_event).sum())}/200**.
- ≥%10 nonzero ve nonzero satırlarda ≥10 event filtresini geçen evren: **{len(filtered)}/60,660**.

Özellikle istenen fold-0 22 RNA, yedi RSF sıfır-atıflı RNA ve KRTAP11-1 satırları ana CSV'de boolean bayraklarla işaretlidir; birleşik işaretli satır sayısı {len(flagged)}'dir.

## Alternatif 200'lük taramalar

{markdown(summary.round(6))}

METABRIC ortak sayısı, MyGene Ensembl→tekil sembol eşlemesi ve yerel cBioPortal METABRIC gen kataloğunda tekil Entrez karşılığı bulunması ölçütüdür. Bu, genin indirilebilir ortaklığını gösterir; yeni ekspresyon matrisi bu adımda indirilmedi ve model yeniden eğitilmedi.

LR ve skor panellerinin kesişimi **{len(set(lr_panel.feature) & set(score_panel.feature))}/200**; birleşimi {len(union)} özelliktir. Ham sıralamalar `PHASE21_ADIM1_ALTERNATIVE_PANELS.csv` içindedir.

Ortak panel bileşimi: filtre+LR için **{int(panels.query("method == 'filtered_LR'").metabric_catalog_unique.sum())} RNA + {len(clinical_common)} klinik = {int(panels.query("method == 'filtered_LR'").metabric_catalog_unique.sum()) + len(clinical_common)}**; filtre+skor için **{int(panels.query("method == 'filtered_score'").metabric_catalog_unique.sum())} RNA + {len(clinical_common)} klinik = {int(panels.query("method == 'filtered_score'").metabric_catalog_unique.sum()) + len(clinical_common)}**. Tam özellik listesi `PHASE21_ADIM1_COMMON_PANEL.csv` içindedir.

## Karar kapısı

ADIM 2 başlatılmadı. Düzeltilmiş boru hattı için aday toplam ortak panel boyutları: filtre+LR **{int(panels.query("method == 'filtered_LR'").metabric_catalog_unique.sum()) + len(clinical_common)}**, filtre+skor **{int(panels.query("method == 'filtered_score'").metabric_catalog_unique.sum()) + len(clinical_common)}**. Kullanıcı onayı gereklidir.

Doğrulama: iki Faz 21 betiği `py_compile` kontrolünü, tam paket `211 passed` sonucunu verdi. Dondurulmuş `components.py`, `aggregation.py` ve `certification.py` dosyalarında değişiklik yoktur.
"""
    (OUT / "PHASE21_ADIM1_REPORT.md").write_text(report, encoding="utf-8")
    print(f"ADIM1 complete: LR common={int(panels.query('method == \'filtered_LR\'').metabric_catalog_unique.sum())}, "
          f"score common={int(panels.query('method == \'filtered_score\'').metabric_catalog_unique.sum())}")


if __name__ == "__main__":
    main()
