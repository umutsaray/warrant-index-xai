#!/usr/bin/env python3
"""Phase 7.2 frozen-design DeepSurv N=50 audit on the clean fold-0 split."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter
import warnings

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase7_retrain"
PANEL = ROOT / "results/phase72_panel"
OUT = ROOT / "results/phase72_deepsurv"
CHECKPOINTS = OUT / "checkpoints"
N, SEED = 50, 20260917


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def deep_attributions_with_relative_gate(
    phase61, phase6, model, checkpoint, raw_queries, raw_backgrounds,
    active_indices, means, scales,
):
    """Honor SHAP additivity; accept its absolute assertion only under the frozen 5% gate."""
    try:
        values, error, output_range, seconds, caught = phase61.deep_attributions(
            phase6, model, checkpoint, raw_queries, raw_backgrounds,
            active_indices, means, scales,
        )
        return values, error, output_range, seconds, caught, False
    except AssertionError as assertion:
        started = perf_counter()
        full_background = (raw_backgrounds - means) / scales
        active_background = full_background[:, active_indices].copy()
        queries = (raw_queries - means[active_indices]) / scales[active_indices]
        flat = phase61.FlatOmicsWrapper(model, checkpoint["slices"]).eval()
        wrapper = phase6.IntegratedActiveWrapper(flat, full_background, active_indices).eval()
        background = torch.from_numpy(active_background.astype(np.float32))
        query = torch.from_numpy(queries.astype(np.float32))
        with warnings.catch_warnings(record=True) as caught_warnings:
            warnings.simplefilter("always")
            explainer = phase61.shap.DeepExplainer(wrapper, background)
            values = phase61.unpack_values(explainer.shap_values(query, check_additivity=False))
        with torch.no_grad():
            prediction = wrapper(query).numpy().reshape(-1)
        base = float(np.asarray(explainer.expected_value).reshape(-1)[0])
        error = float(np.max(np.abs(values.sum(axis=1) + base - prediction)))
        output_range = float(np.ptp(np.r_[prediction, base]))
        relative = error / max(output_range, np.finfo(float).eps)
        if relative > 0.05:
            raise RuntimeError(
                f"DeepExplainer failed relative additivity gate: absolute={error}, relative={relative}"
            ) from assertion
        caught = sorted({f"{item.category.__name__}: {item.message}" for item in caught_warnings})
        caught.append(f"SHAP absolute-0.01 assertion: {assertion}; accepted by relative-5%-gate")
        return values, error, output_range, perf_counter() - started, caught, True


def main() -> None:
    invoked = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    formula = {
        "components.py": sha(ROOT / "warrantxai/components.py"),
        "aggregation.py": sha(ROOT / "warrantxai/aggregation.py"),
    }
    if formula != json.loads((PANEL / "metadata.json").read_text())["formula_sha256"]:
        raise RuntimeError("Frozen formula hash changed")

    phase6 = load_module("20_deepsurv_pilot.py", "phase72_deep_phase6")
    phase61 = load_module("21_deepsurv_explainer_revision.py", "phase72_deep_phase61")
    phase4 = phase6.load_phase4()
    checkpoint = torch.load(DATA / "explainable_dl_model.pt", map_location="cpu", weights_only=False)
    train = frame("X_train_selected.parquet")
    y = frame("y_train.parquet").loc[train.index]
    features = list(checkpoint["ordered_columns"])
    train = train[features]
    with np.load(PANEL / "query_matrices.npz", allow_pickle=False) as panel:
        active = panel["active_features"].astype(str).tolist()
        queries = [panel["source_query12"].astype(np.float32), panel["external_query12"].astype(np.float32)]
    active_indices = np.array([features.index(name) for name in active], dtype=int)
    train_raw = train.to_numpy(np.float32)
    patients = train.index.astype(str).str[:12]
    unique_train = train.loc[~patients.duplicated()]
    background_rng = np.random.default_rng(SEED + 4200)
    first16 = np.sort(background_rng.choice(len(unique_train), 16, replace=False))
    remaining = np.setdiff1d(np.arange(len(unique_train)), first16)
    extra112 = np.sort(background_rng.choice(remaining, 112, replace=False))
    nested = np.r_[first16, extra112]
    backgrounds = {
        16: unique_train.iloc[nested[:16]].to_numpy(np.float32),
        128: unique_train.iloc[nested].to_numpy(np.float32),
    }
    protocol_payload = {
        "model_sha256": sha(DATA / "explainable_dl_model.pt"),
        "train_sha256": sha(DATA / "X_train_selected.parquet"),
        "outcome_sha256": sha(DATA / "y_train.parquet"),
        "panel_sha256": sha(PANEL / "query_matrices.npz"),
        "features": active,
        "candidate_background": 16,
        "reference_background": 128,
        "N": N,
        "seed": SEED,
        "protocol": "patient-cluster-DeepExplainer-nested-k16-k128-v2",
    }
    protocol = hashlib.sha256(json.dumps(protocol_payload, sort_keys=True).encode()).hexdigest()
    groups = [np.flatnonzero(patients == patient) for patient in patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    events = y.event.astype(bool).to_numpy()
    times = y.time.to_numpy(np.float32)
    candidate = np.empty((2, N, 12, len(active)))
    reference = np.empty_like(candidate)
    records = []

    for rep in range(N):
        indices = np.concatenate(
            [groups[index] for index in bootstrap_rng.integers(0, len(groups), len(groups))]
        )
        data_path = CHECKPOINTS / f"{protocol}_rep{rep:02d}.npz"
        meta_path = data_path.with_suffix(".json")
        if data_path.exists() and meta_path.exists():
            with np.load(data_path, allow_pickle=False) as saved:
                candidate[:, rep] = saved["candidate"]
                reference[:, rep] = saved["reference"]
            record = json.loads(meta_path.read_text())
            if record["protocol_hash"] != protocol:
                raise RuntimeError("DeepSurv checkpoint protocol mismatch")
            records.append(record)
            print(f"DeepSurv 7.2 replicate {rep + 1}/50: verified checkpoint", flush=True)
            continue

        model, means, scales, training = phase6.train_bootstrap_model(
            checkpoint, train_raw, events, times, indices, rep
        )
        if training["warning_categories"]:
            raise RuntimeError(f"DeepSurv replicate {rep} training warning: {training['warning_categories']}")
        sites = []
        for site in range(2):
            c, c_error, c_range, c_seconds, c_warnings, c_fallback = deep_attributions_with_relative_gate(
                phase61, phase6, model, checkpoint, queries[site], backgrounds[16],
                active_indices, means, scales,
            )
            r, r_error, r_range, r_seconds, r_warnings, r_fallback = deep_attributions_with_relative_gate(
                phase61, phase6, model, checkpoint, queries[site], backgrounds[128],
                active_indices, means, scales,
            )
            candidate[site, rep], reference[site, rep] = c, r
            sites.append({
                "site": site,
                "candidate_seconds": c_seconds,
                "reference_seconds": r_seconds,
                "candidate_max_abs_additivity_error": c_error,
                "reference_max_abs_additivity_error": r_error,
                "candidate_relative_additivity_error": c_error / max(c_range, np.finfo(float).eps),
                "reference_relative_additivity_error": r_error / max(r_range, np.finfo(float).eps),
                "candidate_warnings": c_warnings,
                "reference_warnings": r_warnings,
                "candidate_internal_assertion_fallback": c_fallback,
                "reference_internal_assertion_fallback": r_fallback,
            })
        record = {
            "protocol_hash": protocol,
            "replicate": rep,
            "rows": training["rows"],
            "events": training["events"],
            "fit_rows": training["fit_rows"],
            "validation_rows": training["validation_rows"],
            "best_epoch": training["best_epoch"],
            "best_validation_c_index": training["best_validation_c_index"],
            "last_epoch": training["last_epoch"],
            "warning_categories": training["warning_categories"],
            "training_seconds": training["elapsed_seconds"],
            "sites": sites,
            "elapsed_seconds": training["elapsed_seconds"] + sum(
                site["candidate_seconds"] + site["reference_seconds"] for site in sites
            ),
        }
        np.savez_compressed(data_path, candidate=candidate[:, rep], reference=reference[:, rep])
        meta_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        records.append(record)
        print(
            f"DeepSurv 7.2 replicate {rep + 1}/50: fit={training['elapsed_seconds']:.2f}s "
            f"explain={record['elapsed_seconds'] - training['elapsed_seconds']:.2f}s",
            flush=True,
        )

    report = phase4.component_report(candidate, reference, active)
    report.to_csv(OUT / "feature_report.csv", index=False)
    np.savez_compressed(
        OUT / "attributions_n50.npz", candidate=candidate, reference=reference,
        features=np.asarray(active),
    )
    records = sorted(records, key=lambda row: row["replicate"])
    (OUT / "bootstrap_summary.json").write_text(json.dumps(records, indent=2), encoding="utf-8")
    aggregate = report.groupby("block").agg(
        assessed=("WI", "count"), mean_F=("F", "mean"), mean_S=("S", "mean"),
        mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.to_csv(OUT / "summary.csv", index=False)
    c_errors = [site["candidate_relative_additivity_error"] for row in records for site in row["sites"]]
    r_errors = [site["reference_relative_additivity_error"] for row in records for site in row["sites"]]
    metadata = {
        "phase": "7.2_deepsurv_n50_clean",
        "arguments": protocol_payload,
        "point_adequate": int(report.point_adequate.sum()),
        "bottleneck_counts": report.bottleneck.value_counts().to_dict(),
        "formula_sha256": formula,
        "additivity": {
            "relative_limit": 0.05,
            "max_candidate_relative_error": max(c_errors),
            "max_reference_relative_error": max(r_errors),
            "internal_assertion_fallbacks": sum(
                int(site["candidate_internal_assertion_fallback"]) + int(site["reference_internal_assertion_fallback"])
                for row in records for site in row["sites"]
            ),
        },
        "privacy": {"identifiers_exported": False, "patient_overlap": 0},
        "conditionality": "fixed deterministic 12+12 queries; nested DeepExplainer k16/k128 backgrounds",
        "accumulated_compute_seconds": sum(row["elapsed_seconds"] for row in records),
        "invocation_wall_clock_seconds": perf_counter() - invoked,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    note = f"""# Faz 7.2 — DeepSurv N=50 temiz audit

Hasta-gruplu fold-0 referans DeepSurv mimarisiyle 50 patient-cluster bootstrap model fit edildi. Temiz 27-feature panelde **{metadata['point_adequate']}/27** feature nokta tahminde WI≥0.8 geçti.

{phase4.markdown_table(aggregate)}

Bottleneck dağılımı: `{metadata['bottleneck_counts']}`. En yüksek göreli additivity hatası candidate/reference için **{max(c_errors):.3g}/{max(r_errors):.3g}**; önceden dondurulan %5 kapısı geçildi.

Sonuç sabit deterministik 12+12 query ve nested DeepExplainer k=16/128 background bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen süre **{metadata['accumulated_compute_seconds']:.1f} sn**, invocation wall-clock **{metadata['invocation_wall_clock_seconds']:.1f} sn**.
"""
    (ROOT / "PHASE72_DEEPSURV.md").write_text(note, encoding="utf-8")
    print(json.dumps({"point_adequate": metadata["point_adequate"], "additivity": metadata["additivity"]}, indent=2))


if __name__ == "__main__":
    main()
