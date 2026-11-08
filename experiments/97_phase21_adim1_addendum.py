#!/usr/bin/env python3
"""Phase 21 ADIM 1 addendum: BH diagnostics and targeted METABRIC validation."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "phase21_reviewer1"
RAW_OLD = ROOT / "data/external/metabric_phase72"
RAW = ROOT / "data/external/metabric_phase21"
BASE = "https://www.cbioportal.org/api"
STUDY = "brca_metabric"
PROFILE = "brca_metabric_mrna"
SAMPLE_LIST = "brca_metabric_mrna"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    OUT.mkdir(exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    panels = pd.read_csv(OUT / "PHASE21_ADIM1_ALTERNATIVE_PANELS.csv")
    forensic = pd.read_csv(OUT / "PHASE21_ADIM1_WALD_FORENSICS.csv")
    original = pd.read_csv(ROOT / "results/phase7_retrain/univariate_cox_rna.tsv", sep="\t").set_index("feature")

    bh_rows = []
    for method, pcol in [("filtered_LR", "lr_p_breslow"), ("filtered_score", "score_p_breslow")]:
        cell = panels[panels.method == method].copy()
        # Reconstruct the full filtered family from the same fold-0 matrix.
        x = pd.read_parquet(ROOT / "results/phase7_retrain/X_train.parquet").set_index("sample_key")
        y = pd.read_parquet(ROOT / "results/phase7_retrain/y_train.parquet").set_index("sample_key").loc[x.index]
        rna = [c for c in x if c.startswith("rna_")]
        values = x[rna].to_numpy(float, copy=False)
        nz = values != 0
        eligible = (nz.sum(0) >= np.ceil(.10 * len(x))) & (nz[y.event.astype(bool).to_numpy()].sum(0) >= 10)
        eligible_names = np.asarray(rna)[eligible]
        if method == "filtered_LR":
            ll0_order = np.argsort(-y.time.to_numpy(float), kind="stable")
            times = y.time.to_numpy(float)[ll0_order]
            events = y.event.to_numpy(bool)[ll0_order]
            ends = np.empty(len(times), int)
            start = 0
            while start < len(times):
                end = start
                while end + 1 < len(times) and times[end + 1] == times[start]: end += 1
                ends[start:end + 1] = end
                start = end + 1
            ll0 = -np.log(ends[np.flatnonzero(events)] + 1).sum()
            p_all = pd.Series(np.nan, index=rna)
            stored = original.reindex(rna)
            p_all[:] = __import__("scipy").stats.chi2.sf(
                np.maximum(0, 2 * (stored.partial_log_likelihood.to_numpy(float) - ll0)), 1)
        else:
            # The exact score p-values for the selected family are preserved in the panel CSV;
            # recompute the full family via the Phase 21 helper.
            import importlib.util, sys
            path = ROOT / "experiments/96_phase21_adim1_wald_forensics.py"
            spec = importlib.util.spec_from_file_location("phase21_wald_helper", path)
            mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod
            assert spec.loader; spec.loader.exec_module(mod)
            p_all = pd.Series(mod.score_pvalues(values, y.time.to_numpy(float), y.event.to_numpy(bool)), index=rna)
        family_p = p_all.loc[eligible_names].to_numpy(float)
        q_all = multipletests(family_p, method="fdr_bh")[1]
        q_lookup = pd.Series(q_all, index=eligible_names)
        top_q = q_lookup.loc[cell.feature].to_numpy(float)
        for row, q in zip(cell.itertuples(), top_q):
            bh_rows.append({"method": method, "feature": row.feature, "rank": row.rank,
                            "p_value": getattr(row, pcol), "bh_q_full_filtered_family": q,
                            "q_lt_0_05": q < .05, "q_lt_0_10": q < .10,
                            "family_m": len(family_p)})
    bh = pd.DataFrame(bh_rows)
    bh.to_csv(OUT / "PHASE21_ADIM1_BH_QVALUES.csv", index=False)

    failed = forensic[forensic.beta.isna()].copy()
    failed["original_wald_beta"] = original.loc[failed.feature, "coefficient"].to_numpy(float)
    failed["ridge_sign_matches_original_wald"] = (
        np.sign(failed.penalized_beta.to_numpy(float)) == np.sign(failed.original_wald_beta.to_numpy(float)))
    failed.to_csv(OUT / "PHASE21_ADIM1_FAILED87.csv", index=False)

    targets = panels[panels.metabric_catalog_unique].drop_duplicates("feature")
    genes = json.loads((RAW_OLD / "genes.json").read_text(encoding="utf-8"))
    catalog = {row["hugoGeneSymbol"]: row["entrezGeneId"] for row in genes}
    target_ids = sorted({int(catalog[s]) for s in targets.symbol if s in catalog})
    if len(target_ids) != targets.symbol.nunique():
        raise RuntimeError("Catalog-matched target count changed before expression download")

    batch_records = []
    molecular = []
    recorded = datetime.now(timezone.utc).isoformat()
    for start in range(0, len(target_ids), 25):
        ids = target_ids[start:start + 25]
        path = RAW / f"mrna_batch_{start // 25:02d}.json"
        if not path.exists():
            response = requests.post(
                BASE + f"/molecular-profiles/{PROFILE}/molecular-data/fetch",
                params={"projection": "DETAILED"},
                json={"sampleListId": SAMPLE_LIST, "entrezGeneIds": ids}, timeout=120)
            if response.status_code != 200:
                raise RuntimeError(f"METABRIC expression download failed: batch={start//25}, HTTP={response.status_code}")
            response.json()
            path.write_bytes(response.content)
        body = json.loads(path.read_text(encoding="utf-8"))
        molecular.extend(body)
        batch_records.append({"file": path.name, "sha256": sha(path), "bytes": path.stat().st_size,
                              "entrez_count": len(ids), "record_count": len(body)})

    mol = pd.DataFrame(molecular)
    if mol.empty or mol.duplicated(["sampleId", "entrezGeneId"]).any():
        raise RuntimeError("Downloaded METABRIC expression is empty or has duplicate keys")
    expression = mol.pivot(index="sampleId", columns="entrezGeneId", values="value")
    samples = pd.DataFrame(json.loads((RAW_OLD / "samples.json").read_text(encoding="utf-8"))).set_index("sampleId")
    patients_raw = pd.DataFrame(json.loads((RAW_OLD / "patient_clinical.json").read_text(encoding="utf-8")))
    patients = patients_raw.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    linked = samples.loc[expression.index, "patientId"]
    time = pd.to_numeric(patients.reindex(linked)["OS_MONTHS"], errors="coerce").set_axis(expression.index)
    status = patients.reindex(linked)["OS_STATUS"].set_axis(expression.index)
    eligible = time.notna() & time.gt(0) & status.isin(["0:LIVING", "1:DECEASED"])
    if int(eligible.sum()) != 1979:
        raise RuntimeError(f"Expected same 1,979-patient cohort, found {int(eligible.sum())}")
    expression = expression.loc[eligible]

    validation_rows = []
    for row in targets.itertuples():
        entrez = int(catalog[row.symbol])
        exists = entrez in expression.columns
        values = expression[entrez] if exists else pd.Series(dtype=float)
        missing = int(values.isna().sum()) if exists else 1979
        complete = int(values.notna().sum()) if exists else 0
        sd = float(values.std(ddof=0)) if complete else np.nan
        validation_rows.append({"feature": row.feature, "ensembl": row.feature[4:].split(".")[0], "symbol": row.symbol,
                                "entrez_gene_id": entrez, "metabric_expression_present": exists,
                                "missing_value_count": missing, "complete_patient_count": complete,
                                "sd": sd, "nonconstant_sd_gt_0": bool(np.isfinite(sd) and sd > 0),
                                "real_expression_eligible": bool(exists and missing == 0 and sd > 0)})
    validation = pd.DataFrame(validation_rows)
    validation.to_csv(OUT / "PHASE21_ADIM1_METABRIC_EXPRESSION_VALIDATION.csv", index=False)
    # Catalog presence is broader than the METABRIC microarray profile. Genes
    # absent from the molecular response are expected platform non-coverage,
    # not a download failure; they remain explicitly ineligible below.

    enriched = panels.merge(validation[["feature", "real_expression_eligible"]], on="feature", how="left")
    real_counts = enriched.groupby("method").real_expression_eligible.sum().astype(int)
    bh_summary = bh.groupby("method").agg(q_lt_0_05=("q_lt_0_05", "sum"), q_lt_0_10=("q_lt_0_10", "sum"),
                                           q_min=("bh_q_full_filtered_family", "min"),
                                           q_median=("bh_q_full_filtered_family", "median"),
                                           q_max=("bh_q_full_filtered_family", "max")).reset_index()
    fail_summary = pd.DataFrame([{
        "failed_n": len(failed), "zero_rate_min": failed.zero_rate.min(),
        "zero_rate_median": failed.zero_rate.median(), "zero_rate_max": failed.zero_rate.max(),
        "nonzero_n_min": failed.nonzero_n.min(), "nonzero_n_median": failed.nonzero_n.median(),
        "nonzero_n_max": failed.nonzero_n.max(), "nonzero_events_min": failed.nonzero_events.min(),
        "nonzero_events_median": failed.nonzero_events.median(), "nonzero_events_max": failed.nonzero_events.max(),
        "complete_separation_n": int((failed.all_nonzero_are_events | failed.no_nonzero_is_event).sum()),
        "ridge_sign_agreement_n": int(failed.ridge_sign_matches_original_wald.sum()),
        "ridge_sign_agreement_rate": float(failed.ridge_sign_matches_original_wald.mean())}])
    fail_summary.to_csv(OUT / "PHASE21_ADIM1_FAILED87_SUMMARY.csv", index=False)
    provenance = {"study_id": STUDY, "molecular_profile_id": PROFILE, "sample_list_id": SAMPLE_LIST,
                  "download_recorded_utc": recorded, "eligible_expression_cohort": 1979,
                  "raw_directory_gitignored": str(RAW.relative_to(ROOT)), "batches": batch_records}
    (RAW / "download_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")

    report = f"""# Faz 21 ADIM 1 eki

## BH q-değerleri

{bh_summary.to_csv(index=False)}
Q-değerleri 43,877 özelliklik filtrelenmiş ailenin tamamı üzerinde BH ile hesaplandı; özet ilk 200 içindir.

## Başarısız 87 Wald özelliği

{fail_summary.to_csv(index=False)}
Başarısızlıkların yalnız {int(fail_summary.complete_separation_n.iloc[0])}/87'si basit tam ayrışma göstergesine sahiptir; baskın desen çok yüksek sıfır oranı ve az nonzero/event desteğidir. Ridge işareti, orijinal tarama katsayısıyla {int(fail_summary.ridge_sign_agreement_n.iloc[0])}/87 özellikte aynıdır.

## Gerçek METABRIC ekspresyon doğrulaması

- Study/profile/sample list: `{STUDY}` / `{PROFILE}` / `{SAMPLE_LIST}`.
- İndirme kaydı (UTC): `{recorded}`.
- Klinik-uygun ekspresyon kohortu: 1,979 hasta.
- Ham hasta-düzeyi yanıtlar yalnız git-ignore'lı `{RAW.relative_to(ROOT)}` altında; commit edilmez.
- Batch SHA-256 ve boyutlar yerel `download_provenance.json` içinde kayıtlıdır.

Gerçek ekspresyon ölçütü: gen sütunu var, 1,979 hastada eksik yok ve SD>0. Filtre+LR gerçek ortak RNA **{real_counts['filtered_LR']}**, filtre+skor gerçek ortak RNA **{real_counts['filtered_score']}**. Gen düzeyi ayrıntılar `PHASE21_ADIM1_METABRIC_EXPRESSION_VALIDATION.csv` içindedir.
"""
    (OUT / "PHASE21_ADIM1_ADDENDUM.md").write_text(report, encoding="utf-8")
    print(f"ADIM1 addendum complete: LR={real_counts['filtered_LR']}, score={real_counts['filtered_score']}")


if __name__ == "__main__":
    main()
