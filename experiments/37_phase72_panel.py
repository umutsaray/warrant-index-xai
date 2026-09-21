#!/usr/bin/env python3
"""Rebuild the Phase 7.2 TCGA/METABRIC audit panel from clean fold-0 data."""

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


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "phase7_retrain"
OUT = ROOT / "results" / "phase72_panel"
RAW = ROOT / "data" / "external" / "metabric_phase72"
SEED = 20260917


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_phase3():
    path = ROOT / "experiments" / "06_brca_real.py"
    spec = importlib.util.spec_from_file_location("phase72_phase3_helpers", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def main() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    formula = {
        "components.py": sha(ROOT / "warrantxai/components.py"),
        "aggregation.py": sha(ROOT / "warrantxai/aggregation.py"),
    }
    expected = {
        "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    if formula != expected:
        raise RuntimeError(f"Frozen formula hashes changed: {formula}")

    train, test = frame("X_train_selected.parquet"), frame("X_test_selected.parquet")
    cox = joblib.load(DATA / "coxph_baseline_v2.joblib")
    rsf = joblib.load(DATA / "rsf_baseline.joblib")
    deep = __import__("torch").load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    cox_features = list(cox["features"])
    rsf_features = list(rsf["features"])
    deep_features = list(deep["ordered_columns"])
    if rsf_features != list(train.columns) or deep_features != list(train.columns):
        raise RuntimeError("Fold-0 model feature order differs from clean selected matrix")

    phase3 = load_phase3()
    phase3.RAW = RAW
    external, mapping, external_info = phase3.external_columns(list(train.columns))
    active_features = [feature for feature in cox_features if feature in external]
    if not active_features:
        raise RuntimeError("No common CoxPH/RSF/DeepSurv/METABRIC features")
    external_raw = np.column_stack([external[feature] for feature in active_features])
    complete = np.isfinite(external_raw).all(axis=1)
    external_raw = external_raw[complete]
    if len(external_raw) < 12:
        raise RuntimeError("Fewer than 12 complete METABRIC query rows")
    external_std = external_raw.std(axis=0, ddof=0)
    if np.any(external_std <= 0):
        constants = [active_features[i] for i in np.flatnonzero(external_std <= 0)]
        raise RuntimeError(f"Common external panel contains constants: {constants}")

    test_patients = test.index.astype(str).str[:12]
    audit = test.loc[~test_patients.duplicated()].copy()
    audit_patients = audit.index.astype(str).str[:12]
    train_patients = set(train.index.astype(str).str[:12])
    if train_patients & set(audit_patients):
        raise RuntimeError("Patient overlap remains in clean fold-0 audit pool")

    active_indices = np.array([train.columns.get_loc(feature) for feature in active_features])
    train_mean = train.mean(axis=0).to_numpy(float)
    train_std = train.std(axis=0, ddof=0).to_numpy(float)
    active_mean, active_std = train_mean[active_indices], train_std[active_indices]
    if np.any(active_std <= 0):
        raise RuntimeError("Common TCGA panel contains constant train features")
    external_z = (external_raw - external_raw.mean(axis=0)) / external_std
    external_aligned = active_mean + external_z * active_std
    source_pool = audit[active_features].to_numpy(float)

    rng = np.random.default_rng(SEED + 4000)
    source_selection = np.sort(rng.choice(len(source_pool), 12, replace=False))
    external_selection = np.sort(rng.choice(len(external_aligned), 12, replace=False))
    source_query = source_pool[source_selection]
    external_query = external_aligned[external_selection]
    np.savez_compressed(
        OUT / "query_matrices.npz",
        source_full=source_pool,
        external_full=external_aligned,
        source_query12=source_query,
        external_query12=external_query,
        active_features=np.asarray(active_features),
        model_features=np.asarray(train.columns),
        cox_features=np.asarray(cox_features),
    )
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)
    payload = {
        "old": {"common_features": 77, "tcga_audit_patients": 141, "metabric_complete_patients": 1466},
        "new": {
            "fold0_test_rows": len(test),
            "fold0_test_unique_patients": int(test_patients.nunique()),
            "tcga_audit_patients_after_one_sample_per_patient": len(audit),
            "common_features": len(active_features),
            "metabric_complete_patients": len(external_raw),
            "tcga_query_rows": 12,
            "metabric_query_rows": 12,
        },
        "panel_definition": "intersection of CoxPH-v2 fitted features, shared RSF/DeepSurv inputs, and strict METABRIC mapping",
        "active_features": active_features,
        "query_selection": {"seed": SEED + 4000, "deterministic_without_replacement": True},
        "hashes": {
            "query_matrices.npz": sha(OUT / "query_matrices.npz"),
            "coxph_baseline_v2.joblib": sha(DATA / "coxph_baseline_v2.joblib"),
            "rsf_baseline.joblib": sha(DATA / "rsf_baseline.joblib"),
            "explainable_dl_model.pt": sha(DATA / "explainable_dl_model.pt"),
        },
        "formula_sha256": formula,
        "external_inventory": external_info,
        "privacy": {"identifiers_exported": False},
        "wall_clock_seconds": perf_counter() - started,
    }
    (OUT / "metadata.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    report = f"""# Faz 7.2 ADIM 0 — Audit paneli

| Ölçü | Eski | Yeni |
|---|---:|---:|
| Ortak feature | 77 | **{len(active_features)}** |
| TCGA audit hastası | 141 | **{len(audit)}** |
| Tam METABRIC query hastası | 1466 | **{len(external_raw)}** |

Yeni fold-0 test setindeki {len(test)} satır {test_patients.nunique()} benzersiz hastaya aittir. Train ile hasta kesişimi sıfırdır; eski 23-hasta çıkarma adımı uygulanmadı. Hasta başına ilk deterministik satır tutularak audit havuzu {len(audit)} hastaya indirildi.

Panel, CoxPH-v2'nin gerçekten fit edilen feature'ları ile RSF/DeepSurv ortak girdileri ve katı METABRIC eşlemesinin kesişimidir. RSF/DeepSurv için seed={SEED + 4000} ile yeni havuzlardan 12+12 satır kimliksiz ve replacement olmadan seçildi. Query matrisleri `results/phase72_panel/query_matrices.npz` içinde saklıdır.

Formül SHA256 doğrulaması geçti. Gerçek kimlik dışa aktarılmadı. Wall-clock **{payload['wall_clock_seconds']:.3f} saniye**.
"""
    (ROOT / "PHASE72_AUDIT_PANEL.md").write_text(report, encoding="utf-8")
    print(json.dumps(payload["new"], indent=2))


if __name__ == "__main__":
    main()
