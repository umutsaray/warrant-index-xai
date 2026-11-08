#!/usr/bin/env python3
"""Phase 21 ADIM 0: read-only endpoint, cohort, history and artifact inventory."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter

ROOT = Path(__file__).resolve().parents[1]
SIBLING = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
DATA = ROOT / "results" / "phase7_retrain"
RAW = ROOT / "data" / "external" / "metabric_phase72"
OUT = ROOT / "phase21_reviewer1"


def frame(path: Path) -> pd.DataFrame:
    value = pd.read_parquet(path)
    return value.set_index("sample_key") if "sample_key" in value else value


def reverse_km(time: pd.Series, event: pd.Series) -> float:
    km = KaplanMeierFitter().fit(time.astype(float), event_observed=~event.astype(bool))
    return float(km.median_survival_time_)


def endpoint_row(cohort: str, time: pd.Series, event: pd.Series, unit: str) -> dict:
    return {
        "section": "survival_endpoint", "cohort": cohort, "n": len(time),
        "events": int(event.sum()), "event_rate": float(event.mean()),
        "censor_rate": float(1 - event.mean()), "time_unit": unit,
        "median_followup_reverse_km": reverse_km(time, event),
        "endpoint": "overall survival", "time_origin": "diagnosis",
        "event_coding": "1=death, 0=right-censored",
    }


def markdown(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in frame.itertuples(index=False, name=None):
        values = [f"{value:.6g}" if isinstance(value, float) else str(value) for value in row]
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    train_y, test_y = frame(DATA / "y_train.parquet"), frame(DATA / "y_test.parquet")
    train_x, test_x = frame(DATA / "X_train.parquet"), frame(DATA / "X_test.parquet")
    rows = [
        endpoint_row("TCGA fold-0 train", train_y.time, train_y.event, "days"),
        endpoint_row("TCGA fold-0 test", test_y.time, test_y.event, "days"),
    ]

    patients_raw = pd.DataFrame(json.loads((RAW / "patient_clinical.json").read_text(encoding="utf-8")))
    patients = patients_raw.pivot(index="patientId", columns="clinicalAttributeId", values="value")
    samples = pd.DataFrame(json.loads((RAW / "samples.json").read_text(encoding="utf-8"))).set_index("sampleId")
    molecular = []
    for path in sorted(RAW.glob("mrna_batch_*.json")):
        molecular.extend(json.loads(path.read_text(encoding="utf-8")))
    expression_samples = pd.Index(pd.DataFrame(molecular).sampleId.unique())
    linked = samples.loc[expression_samples, "patientId"]
    os_time = pd.to_numeric(patients.reindex(linked)["OS_MONTHS"], errors="coerce")
    os_status = patients.reindex(linked)["OS_STATUS"]
    eligible = os_time.notna() & os_time.gt(0) & os_status.isin(["0:LIVING", "1:DECEASED"])
    rows.append(endpoint_row("METABRIC expression cohort", os_time[eligible].reset_index(drop=True),
                             os_status[eligible].eq("1:DECEASED").reset_index(drop=True), "months"))

    all_index = train_x.index.append(test_x.index).astype(str)
    patient = all_index.str[:12]
    sample_code = all_index.str[13:15]
    label = pd.Series(np.select([sample_code == "01", sample_code == "06"],
                                ["primary solid tumor", "metastatic"], default="other"))
    for kind, count in label.value_counts().items():
        rows.append({"section": "tcga_sample_type", "cohort": "TCGA retained 842 rows",
                     "category": kind, "n": int(count)})
    counts = pd.Series(patient).value_counts()
    rows.extend([
        {"section": "tcga_repeated_patients", "cohort": "TCGA retained 842 rows",
         "category": "unique patients", "n": int(counts.size)},
        {"section": "tcga_repeated_patients", "cohort": "TCGA retained 842 rows",
         "category": "patients with >1 retained sample", "n": int((counts > 1).sum())},
        {"section": "tcga_repeated_patients", "cohort": "TCGA retained 842 rows",
         "category": "rows belonging to repeated patients", "n": int(counts[counts > 1].sum())},
    ])

    clinical_raw = pd.read_csv(SIBLING / "TCGA-BRCA.clinical.tsv.gz", sep="\t", low_memory=False)
    clinical_raw["sample_key"] = clinical_raw["sample"].astype(str).str[:15]
    clinical_raw = clinical_raw.drop_duplicates("sample_key", keep="first").set_index("sample_key")
    clinical = clinical_raw.reindex(all_index)
    clinical_cols = ["age_at_index.demographic", "gender.demographic", "race.demographic",
                     "ethnicity.demographic", "ajcc_pathologic_stage.diagnoses"]
    for col in clinical_cols:
        rows.append({"section": "clinical_missing_before_imputation", "cohort": "TCGA 842 retained rows",
                     "category": col, "n": int(clinical[col].isna().sum())})
    rows.extend([
        {"section": "exclusion_imputation", "cohort": "TCGA", "category": "rows excluded for missing survival flags", "n": 11},
        {"section": "exclusion_imputation", "cohort": "TCGA fold-0 train", "category": "methylation missing cells before KNN", "n": 729119},
        {"section": "exclusion_imputation", "cohort": "TCGA fold-0 test", "category": "methylation missing cells before KNN", "n": 182284},
        {"section": "exclusion_imputation", "cohort": "TCGA fold-0 train/test", "category": "missing cells after fitted-train imputers", "n": 0},
    ])

    diagnostic = json.loads((ROOT / "results/phase70b_collapse_diagnostic/deepsurv_diagnostic.json").read_text())
    for key in ["row_level_validation", "patient_grouped_validation"]:
        value = diagnostic[key]
        rows.append({"section": "deepsurv_history", "cohort": "fold 0", "category": key,
                     "n": len(value["history"]), "best_epoch": value["best_epoch"],
                     "validation_c_index": value["best_validation_c_index"],
                     "test_c_index": value["best_checkpoint_test_c_index"],
                     "patient_overlap": value["fit_validation_patient_overlap"]})
    rows.append({"section": "deepsurv_history", "cohort": "folds 1-4", "category": "epoch histories preserved", "n": 0})

    for model, directory, pattern in [
        ("RSF", ROOT / "results/phase72_rsf/checkpoints", "*.npz"),
        ("DeepSurv", ROOT / "results/phase72_deepsurv/checkpoints", "*.npz"),
    ]:
        files = list(directory.glob(pattern))
        rows.append({"section": "bootstrap_artifacts", "cohort": model,
                     "category": "serialized bootstrap model objects", "n": 0,
                     "note": f"{len(files)} compact attribution/checkpoint arrays exist; model objects are not serialized"})
        rows.append({"section": "bootstrap_artifacts", "cohort": model,
                     "category": "deterministic reproduction", "n": 50,
                     "note": "seed 20260917 and frozen Phase 7.2 script/input hashes are recorded"})

    output = pd.DataFrame(rows)
    output.to_csv(OUT / "PHASE21_ADIM0_INVENTORY.csv", index=False)
    endpoint = output[output.section == "survival_endpoint"]
    types = output[output.section == "tcga_sample_type"]
    repeated = output[output.section == "tcga_repeated_patients"]
    clinical_missing = output[output.section == "clinical_missing_before_imputation"]
    report = f"""# Faz 21 ADIM 0 — veri ve sağkalım tanımı

Bu adım salt okunur envanterdir. Mevcut `results/phase21/` eski Faz 2.1 artefaktlarını içerdiği için üzerine yazılmadı; yeni dosyalar `phase21_reviewer1/` altındadır.

## Sağkalım tanımı

Sonlanım **overall survival (OS)**; başlangıç tanıdır. TCGA `OS.time` gün, METABRIC `OS_MONTHS` ay birimindedir. Olay 1=ölüm, 0=sağ sansürdür.

{markdown(endpoint[['cohort','n','events','censor_rate','time_unit','median_followup_reverse_km']])}

## TCGA örnekleri ve eksik veri

{markdown(types[['category','n']])}

{markdown(repeated[['category','n']])}

Model matrisi aynı hastanın birden çok satırını dışlamaz; hasta-gruplu dış split bütün satırları aynı fold'da tutar. Fold-0 audit sorgu havuzu ayrıca hasta başına ilk deterministik satırı tutar. Birleşik 853 satırdan survival kontrol bayrağı taşıyan 11 satır çıkarılmıştır. Methylation KNN (k=5), klinik sayısal medyan, kategorik en-sık değer imputasyonu yalnız train'de fit edilip test'e uygulanmıştır; son eksik hücre sayısı sıfırdır.

Klinik ham eksiklik sayıları:

{markdown(clinical_missing[['category','n']])}

## DeepSurv epoch kayıtları

Epoch-bazlı validation/test eğrileri yalnız fold-0 Faz 7.0b tanı koşusunda saklanmıştır: satır-düzeyi geçmiş 68 epoch, hasta-gruplu tanı geçmişi 61 epoch içerir. Fold 1–4 için yalnız seçilen checkpoint/özet değerleri saklanmış, epoch-bazlı validation/test eğrileri saklanmamıştır. Mevcut fold-0 hasta-gruplu tanı kaydı en iyi validation C-index 0.942184, eşzamanlı test C-index 0.488716 ve sıfır fit/validation hasta kesişimi verir.

## Bootstrap model artefaktları

RSF ve DeepSurv için 50 bootstrap model nesnesi diskte saklı değildir. Kompakt checkpoint/atıf dizileri ve `bootstrap_summary.json` kayıtları vardır. Her iki koşu da seed `20260917`, girdi hash'leri ve protokol hash'leriyle deterministik olarak yeniden üretilebilir; yeniden üretim model fit gerektirir.
"""
    (OUT / "PHASE21_ADIM0_REPORT.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
