#!/usr/bin/env python3
"""Phase 7.5 fold 1 CoxPH-v2 N=50 audit on the rebuilt clean patient panel."""

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
from scipy.stats import norm
from threadpoolctl import threadpool_limits


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "phase75_fold1" / "models"
PANEL = ROOT / "results" / "phase75_fold1_panel"
OUT = ROOT / "results" / "phase75_fold1_cox"
CHECKPOINTS = OUT / "checkpoints"
N, SEED = 50, 20260917

sys.path.insert(0, str(ROOT))
from warrantxai.aggregation import adequate, warrant_index
from warrantxai.components import (
    faithfulness,
    generalizability,
    generalizability_signed,
    importance,
    importance_ci,
    rank_reproducibility,
    stability,
    uncertainty,
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def load_phase3():
    path = ROOT / "experiments" / "06_brca_real.py"
    spec = importlib.util.spec_from_file_location("phase75_fold1_cox_phase3", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    invoked = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    formula = {
        "components.py": sha(ROOT / "warrantxai/components.py"),
        "aggregation.py": sha(ROOT / "warrantxai/aggregation.py"),
    }
    expected = json.loads((PANEL / "metadata.json").read_text())["formula_sha256"]
    if formula != expected:
        raise RuntimeError("Frozen formula hash changed")

    train, test = frame("X_train_selected.parquet"), frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[train.index]
    saved = joblib.load(DATA / "coxph_baseline_v2.joblib")
    features, scaler = list(saved["features"]), saved["scaler"]
    scaled = scaler.transform(train[features].astype(float))
    test_patients = test.index.astype(str).str[:12]
    audit = test.loc[~test_patients.duplicated()]
    anchors = scaler.transform(audit[features].astype(float))
    train_patients = train.index.astype(str).str[:12]
    groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    rng = np.random.default_rng(SEED)
    coefficients, records = [], []
    input_hashes = {
        name: sha(DATA / name)
        for name in ["X_train_selected.parquet", "y_train.parquet", "coxph_baseline_v2.joblib"]
    }
    protocol = hashlib.sha256(
        json.dumps(
            {"inputs": input_hashes, "N": N, "seed": SEED, "protocol": "cluster-cox-separation-retry200-step025-v2"},
            sort_keys=True,
        ).encode()
    ).hexdigest()

    for rep in range(N):
        sampled_groups = rng.integers(0, len(groups), len(groups))
        indices = np.concatenate([groups[index] for index in sampled_groups])
        path = CHECKPOINTS / f"{protocol}_rep{rep:02d}.json"
        if path.exists():
            stored = json.loads(path.read_text())
            coefficients.append(np.asarray(stored["coefficients"], dtype=float))
            records.append(stored["record"])
            print(f"CoxPH 7.2 replicate {rep + 1}/50: verified checkpoint", flush=True)
            continue
        started = perf_counter()
        values = scaled[indices]
        events = y.event.to_numpy()[indices].astype(bool)
        valid = (values[events].var(axis=0) > 1e-12) & (values[~events].var(axis=0) > 1e-12)
        fit_options = {"max_steps": 50}
        attempts = []
        for attempt in range(20):
            fit = pd.DataFrame(values[:, valid], columns=np.asarray(features)[valid])
            fit["time"] = y.time.to_numpy()[indices]
            fit["event"] = events
            model = CoxPHFitter(penalizer=0.1)
            with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=2):
                warnings.simplefilter("always")
                model.fit(fit, "time", "event", show_progress=False, fit_options=fit_options.copy())
            convergence = [item for item in caught if issubclass(item.category, ConvergenceWarning)]
            flagged = [
                index
                for index, feature in enumerate(features)
                if valid[index]
                and any(str(item.message).startswith(f"Column {feature} have very low variance") for item in convergence)
            ]
            attempts.append(
                {
                    "attempt": attempt,
                    "fit_options": fit_options.copy(),
                    "warning_categories": dict(Counter(item.category.__name__ for item in caught)),
                    "separation_features": [features[index] for index in flagged],
                }
            )
            if not convergence:
                break
            if flagged:
                valid[flagged] = False
                continue
            if fit_options["max_steps"] == 50:
                fit_options = {"max_steps": 200, "step_size": 0.25}
                continue
            failure = {
                "replicate": rep,
                "rows": len(indices),
                "events": int(events.sum()),
                "excluded_features": [features[index] for index in np.flatnonzero(~valid)],
                "fit_attempts": attempts,
                "final_warning_messages": [str(item.message) for item in caught],
                "elapsed_seconds": perf_counter() - started,
            }
            (OUT / f"failed_replicate_{rep:02d}.json").write_text(
                json.dumps(failure, indent=2), encoding="utf-8"
            )
            raise RuntimeError(f"Cox replicate {rep} failed numerical retry")
        else:
            raise RuntimeError(f"Cox replicate {rep} exhausted remediation")
        coef = np.zeros(len(features))
        coef[valid] = model.params_.to_numpy()
        record = {
            "replicate": rep,
            "rows": len(indices),
            "events": int(events.sum()),
            "excluded_features": [features[index] for index in np.flatnonzero(~valid)],
            "fit_attempts": attempts,
            "final_warning_categories": dict(Counter(item.category.__name__ for item in caught)),
            "elapsed_seconds": perf_counter() - started,
        }
        if record["final_warning_categories"]:
            raise RuntimeError(f"Cox replicate {rep} retained warnings")
        path.write_text(json.dumps({"coefficients": coef.tolist(), "record": record}, indent=2))
        coefficients.append(coef)
        records.append(record)
        print(
            f"CoxPH 7.2 replicate {rep + 1}/50: events={record['events']} exclusions={len(record['excluded_features'])} time={record['elapsed_seconds']:.1f}s",
            flush=True,
        )

    coefficients = np.asarray(coefficients)
    source_phi = coefficients[:, None, :] * anchors[None, :, :]
    tensor = source_phi[None, ...]
    lower, upper = importance_ci(tensor)
    f = faithfulness(tensor, tensor, 0.5)
    s = stability(tensor, 0.5)
    u = uncertainty(lower, upper, 3.0)
    r = rank_reproducibility(importance(tensor), 0.02)

    panel = json.loads((PANEL / "metadata.json").read_text())
    active = panel["active_features"]
    active_cox = np.array([features.index(feature) for feature in active])
    with np.load(PANEL / "query_matrices.npz", allow_pickle=False) as matrices:
        external_aligned = matrices["external_full"]
    external_scaled = np.column_stack(
        [
            (external_aligned[:, out] - scaler.mean_[active_cox[out]]) / scaler.scale_[active_cox[out]]
            for out in range(len(active))
        ]
    )
    source_active = source_phi[:, :, active_cox]
    external_phi = coefficients[:, active_cox, None] * external_scaled.T[None, :, :]
    external_phi = np.transpose(external_phi, (0, 2, 1))
    gm = generalizability(
        np.stack([np.mean(np.abs(source_active)), np.mean(np.abs(external_phi))])[:, None], 0.5
    )
    # Compute per-feature magnitudes rather than the scalar convenience above.
    site_importance = np.stack(
        [np.mean(np.abs(source_active), axis=(0, 1)), np.mean(np.abs(external_phi), axis=(0, 1))]
    )
    gm = generalizability(site_importance, 0.5)
    phase3 = load_phase3()
    mean_coef = coefficients[:, active_cox].mean(axis=0)
    grid = max(source_active.shape[1], external_phi.shape[1])
    source_signed = np.column_stack(
        [phase3.empirical_rank_grid(mean_coef[j] * anchors[:, active_cox[j]], grid) for j in range(len(active))]
    )
    external_signed = np.column_stack(
        [phase3.empirical_rank_grid(mean_coef[j] * external_scaled[:, j], grid) for j in range(len(active))]
    )
    gs = generalizability_signed(np.stack([source_signed, external_signed]), 0.5)
    gc = np.minimum(gm, gs)

    full = pd.DataFrame({"feature": features, "F": f, "S": s, "U": u, "R": r})
    full.insert(1, "block", [feature.split("_")[0] for feature in features])
    full["G_magnitude"] = np.nan
    full["G_signed"] = np.nan
    full["G_claim"] = np.nan
    lookup = {feature: index for index, feature in enumerate(active)}
    for row, feature in enumerate(features):
        if feature in lookup:
            index = lookup[feature]
            full.loc[row, ["G_magnitude", "G_signed", "G_claim"]] = [gm[index], gs[index], gc[index]]
    z = full[["F", "S", "U", "G_claim", "R"]].to_numpy()
    full["WI"] = warrant_index(z)
    full["point_adequate"] = adequate(z, 0.8)
    full["bottleneck"] = [
        "not assessable" if np.isnan(row).any() else ["F", "S", "U", "G_claim", "R"][int(np.argmin(row))]
        for row in z
    ]
    primary = full[full.feature.isin(active)].set_index("feature").loc[active].reset_index()
    primary.to_csv(OUT / "feature_report.csv", index=False)
    full.to_csv(OUT / "feature_report_all_model_features.csv", index=False)
    np.savez_compressed(
        OUT / "attributions_n50.npz",
        source_phi=source_phi,
        external_phi=external_phi,
        coefficients=coefficients,
        features=np.asarray(features),
        active_features=np.asarray(active),
    )
    (OUT / "bootstrap_summary.json").write_text(json.dumps(records, indent=2))
    elapsed = perf_counter() - invoked
    metadata = {
        "model": "CoxPH-v2",
        "replicates": N,
        "seed": SEED,
        "audit_patients": len(audit),
        "assessed_features": len(primary),
        "point_adequate": int((primary.point_adequate == 1).sum()),
        "bottleneck_counts": primary.bottleneck.value_counts().to_dict(),
        "fit_replicates_with_multiple_attempts": sum(len(row["fit_attempts"]) > 1 for row in records),
        "final_warning_replicates": sum(bool(row["final_warning_categories"]) for row in records),
        "formula_sha256": formula,
        "model_sha256": sha(DATA / "coxph_baseline_v2.joblib"),
        "wall_clock_seconds": elapsed,
        "audit_pipeline_stage": "point estimates only; calibration/certification pending",
        "privacy": {"identifiers_exported": False},
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2))
    report = f"""# Faz 7.5 fold 1 ADIM 1 — CoxPH-v2 N=50 audit

Yeni {len(audit)}-hasta audit havuzu ve {len(primary)}-feature ortak ailede N=50 patient-cluster bootstrap tamamlandı. **{metadata['point_adequate']}/{len(primary)}** feature nokta tahmininde WI≥0.8 geçti. Bottleneck dağılımı `{metadata['bottleneck_counts']}`.

CoxPH candidate ve reference aynı kapalı-form attribution olduğundan F sayısal öz-tutarlılık ölçümüdür. {metadata['fit_replicates_with_multiple_attempts']} replika remediation/retry gerektirdi; nihai warning taşıyan replika sayısı {metadata['final_warning_replicates']}.

Formül SHA256 doğrulandı. Kalibrasyon ve sertifikasyon henüz yapılmadı. Wall-clock **{elapsed:.3f} saniye**.
"""
    (ROOT / "PHASE75_FOLD1_COXPH.md").write_text(report, encoding="utf-8")
    print(json.dumps(metadata, indent=2), flush=True)


if __name__ == "__main__":
    main()
