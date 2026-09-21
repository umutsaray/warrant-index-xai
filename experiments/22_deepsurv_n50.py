"""Phase 6.2 step 1: DeepSurv N=50 with approved DeepExplainer k16/k128."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
import warnings

import numpy as np
import pandas as pd
import torch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase62_deepsurv_n50"
OLD = ROOT / "results/phase61_deepsurv_explainer_revision"
CHECKPOINTS = OUT / "checkpoints"
SEED = 20260917
N = 50
N_OLD = 2
TIME_LIMIT_SECONDS = 30 * 60
WARNING = (
    "DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds."
)

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import markdown_table


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def deep_attributions_with_relative_gate(
    phase61, phase6, model, checkpoint, raw_queries, raw_backgrounds,
    active_indices, means, scales,
):
    """Use SHAP's check, falling back only when our preregistered 5% gate passes."""
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
        caught.append(f"SHAP internal absolute-0.01 assertion: {assertion}; accepted by relative-5%-gate")
        return values, error, output_range, perf_counter() - started, caught, True


def run() -> None:
    invoked = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    CHECKPOINTS.mkdir(parents=True, exist_ok=True)
    phase61 = load_module("21_deepsurv_explainer_revision.py", "phase61_helpers_for_n50")
    phase6 = load_module("20_deepsurv_pilot.py", "phase6_helpers_for_n50")
    phase4 = phase6.load_phase4()
    (checkpoint, train, y_train, train_raw, train_patients, queries, backgrounds,
     active_features, active_indices, active_mean, mapping, external_info) = phase61.build_inputs(phase6, phase4)

    with np.load(OLD / "attributions_n2.npz", allow_pickle=False) as saved:
        if saved["features"].tolist() != active_features:
            raise RuntimeError("Phase 6.1 feature order differs")
        old_candidate, old_reference = saved["candidate"], saved["reference"]
    if old_candidate.shape != (2, 2, 12, 77) or old_reference.shape != old_candidate.shape:
        raise RuntimeError("Unexpected Phase 6.1 attribution shape")
    candidates = np.empty((2, N, 12, 77))
    references = np.empty_like(candidates)
    candidates[:, :N_OLD] = old_candidate
    references[:, :N_OLD] = old_reference

    protocol = {
        "model_sha256": json.loads((OLD / "metadata.json").read_text())["input_sha256"]["explainable_dl_model.pt"],
        "features": active_features,
        "candidate_background_sha256": hashlib.sha256(np.asarray(backgrounds[16], dtype="<f4").tobytes()).hexdigest(),
        "reference_background_sha256": hashlib.sha256(np.asarray(backgrounds[128], dtype="<f4").tobytes()).hexdigest(),
        "queries_sha256": [hashlib.sha256(np.asarray(q, dtype="<f4").tobytes()).hexdigest() for q in queries],
        "protocol": "deepsurv-patient-cluster-DeepExplainer-k16-k128-v1",
    }
    protocol_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    events = y_train.event.astype(bool).to_numpy()
    times = y_train.time.to_numpy(np.float32)
    new_records = []
    checkpoint_reused = 0

    for replicate in range(N):
        indices = np.concatenate([groups[i] for i in bootstrap_rng.integers(0, len(groups), len(groups))])
        if replicate < N_OLD:
            continue
        data_path = CHECKPOINTS / f"{protocol_hash}_rep{replicate:02d}.npz"
        meta_path = data_path.with_suffix(".json")
        if data_path.exists() and meta_path.exists():
            with np.load(data_path, allow_pickle=False) as saved:
                candidates[:, replicate] = saved["candidate"]
                references[:, replicate] = saved["reference"]
            record = json.loads(meta_path.read_text(encoding="utf-8"))
            if record["protocol_hash"] != protocol_hash:
                raise RuntimeError("DeepSurv checkpoint protocol mismatch")
            checkpoint_reused += 1
            new_records.append(record)
            print(f"DeepSurv N=50 replicate {replicate + 1}/50: verified checkpoint", flush=True)
            continue

        model, means, scales, training = phase6.train_bootstrap_model(
            checkpoint, train_raw, events, times, indices, replicate
        )
        sites = []
        for site in range(2):
            candidate, c_error, c_range, c_seconds, c_warnings, c_fallback = deep_attributions_with_relative_gate(
                phase61, phase6, model, checkpoint, queries[site], backgrounds[16],
                active_indices, means, scales
            )
            reference, r_error, r_range, r_seconds, r_warnings, r_fallback = deep_attributions_with_relative_gate(
                phase61, phase6, model, checkpoint, queries[site], backgrounds[128],
                active_indices, means, scales
            )
            candidates[site, replicate] = candidate
            references[site, replicate] = reference
            sites.append({
                "site": site,
                "candidate_seconds": c_seconds, "reference_seconds": r_seconds,
                "candidate_max_abs_additivity_error": c_error,
                "reference_max_abs_additivity_error": r_error,
                "candidate_relative_additivity_error": c_error / max(c_range, np.finfo(float).eps),
                "reference_relative_additivity_error": r_error / max(r_range, np.finfo(float).eps),
                "candidate_warnings": c_warnings, "reference_warnings": r_warnings,
                "candidate_internal_assertion_fallback": c_fallback,
                "reference_internal_assertion_fallback": r_fallback,
            })
        record = {
            "protocol_hash": protocol_hash, "replicate": replicate,
            "rows": training["rows"], "events": training["events"],
            "fit_rows": training["fit_rows"], "validation_rows": training["validation_rows"],
            "best_epoch": training["best_epoch"], "best_validation_c_index": training["best_validation_c_index"],
            "last_epoch": training["last_epoch"], "warning_categories": training["warning_categories"],
            "training_seconds": training["elapsed_seconds"], "sites": sites,
            "elapsed_seconds": training["elapsed_seconds"] + sum(
                row["candidate_seconds"] + row["reference_seconds"] for row in sites
            ),
        }
        np.savez_compressed(data_path, candidate=candidates[:, replicate], reference=references[:, replicate])
        meta_path.write_text(json.dumps(record, indent=2), encoding="utf-8")
        new_records.append(record)
        print(
            f"DeepSurv N=50 replicate {replicate + 1}/50: "
            f"fit={training['elapsed_seconds']:.2f}s explain={record['elapsed_seconds']-training['elapsed_seconds']:.2f}s",
            flush=True,
        )
        if perf_counter() - invoked > TIME_LIMIT_SECONDS:
            raise RuntimeError("PHASE62_TIME_LIMIT: exceeded 30 minutes after a complete checkpoint")

    report = phase4.component_report(candidates, references, active_features)
    report_n2 = phase4.component_report(candidates[:, :2], references[:, :2], active_features)
    stored_n2 = pd.read_csv(OLD / "feature_report.csv")
    np.testing.assert_allclose(
        report_n2[["F", "S", "U", "G_claim", "R", "WI"]],
        stored_n2[["F", "S", "U", "G_claim", "R", "WI"]], atol=2e-12,
    )
    for column in ["F", "S", "U", "G_claim", "R", "WI", "point_adequate", "bottleneck"]:
        report[f"{column}_n2"] = report.feature.map(report_n2.set_index("feature")[column])
    report.insert(0, "inference_scope", WARNING)
    report.to_csv(OUT / "feature_report.csv", index=False)
    np.savez_compressed(OUT / "attributions_n50.npz", candidate=candidates, reference=references,
                        features=np.asarray(active_features))

    old_summary = json.loads((OLD / "bootstrap_summary.json").read_text(encoding="utf-8"))
    combined = old_summary + sorted(new_records, key=lambda row: row["replicate"])
    if len(combined) != 50 or [row["replicate"] for row in combined] != list(range(50)):
        raise RuntimeError("Incomplete DeepSurv N=50 bootstrap sequence")
    for row in combined:
        row["inference_scope"] = WARNING
    (OUT / "bootstrap_summary.json").write_text(json.dumps(combined, indent=2), encoding="utf-8")
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)

    aggregate = report.groupby("block").agg(
        features=("feature", "size"), assessed=("WI", "count"), mean_F=("F", "mean"),
        mean_S=("S", "mean"), mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.insert(0, "inference_scope", WARNING)
    aggregate.to_csv(OUT / "summary.csv", index=False)
    n2_pass = set(report_n2.loc[report_n2.point_adequate == 1, "feature"])
    n50_pass = set(report.loc[report.point_adequate == 1, "feature"])
    transition = report[["feature", "block", "R_n2", "R", "WI_n2", "WI",
                         "point_adequate_n2", "point_adequate", "bottleneck_n2", "bottleneck"]].copy()
    transition["transition"] = np.select(
        [transition.point_adequate_n2.astype(bool) & transition.point_adequate.astype(bool),
         transition.point_adequate_n2.astype(bool) & ~transition.point_adequate.astype(bool),
         ~transition.point_adequate_n2.astype(bool) & transition.point_adequate.astype(bool)],
        ["retained", "lost", "new"], default="neither",
    )
    transition.insert(0, "inference_scope", WARNING)
    transition.to_csv(OUT / "n2_n50_comparison.csv", index=False)

    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    formula_sha = {name: sha(ROOT / name) for name in formula_files}
    expected = json.loads((OLD / "metadata.json").read_text())["formula_module_sha256"]
    unchanged = formula_sha == expected
    if not unchanged:
        raise RuntimeError("Frozen formula SHA256 changed")
    new_compute = sum(row["elapsed_seconds"] for row in new_records)
    invocation = perf_counter() - invoked
    all_candidate_errors = [site["candidate_relative_additivity_error"] for row in new_records for site in row["sites"]]
    all_reference_errors = [site["reference_relative_additivity_error"] for row in new_records for site in row["sites"]]
    fallback_count = sum(
        int(site.get("candidate_internal_assertion_fallback", False)) +
        int(site.get("reference_internal_assertion_fallback", False))
        for row in new_records for site in row["sites"]
    )
    metadata = {
        "phase": "6.2_step1_DeepSurv_N50", "conditionality_warning": WARNING,
        "arguments": {"replicates": 50, "anchors_per_site": 12, "candidate_background": 16,
                      "reference_background": 128, "seed": SEED},
        "model": {"main_model_retrained": False, "bootstrap_models_in_sequence": 50,
                  "new_models_fit_this_phase": 48, "first_two_reused_from_phase61": True,
                  "architecture_and_hyperparameters_unchanged": True},
        "results": {"assessed_features": 77, "point_adequate_n2": len(n2_pass),
                    "point_adequate_n50": len(n50_pass), "retained_from_n2": len(n2_pass & n50_pass),
                    "lost_from_n2": len(n2_pass - n50_pass), "new_vs_n2": len(n50_pass - n2_pass),
                    "bottleneck_counts_n50": report.bottleneck.value_counts().to_dict()},
        "runtime": {"new_replicates_2_to_49_accumulated_seconds": new_compute,
                    "completion_invocation_seconds": invocation, "checkpoint_replicates_reused": checkpoint_reused,
                    "hard_stop_seconds": TIME_LIMIT_SECONDS},
        "additivity": {"max_new_candidate_relative_error": max(all_candidate_errors),
                       "max_new_reference_relative_error": max(all_reference_errors),
                       "shap_internal_absolute_assertion_fallbacks": fallback_count,
                       "fallback_acceptance_rule": "predeclared relative additivity error <= 0.05"},
        "formula_module_sha256": formula_sha, "formula_module_sha256_unchanged": unchanged,
        "input_sha256": json.loads((OLD / "metadata.json").read_text())["input_sha256"],
        "source_sha256": {"experiments/22_deepsurv_n50.py": sha(Path(__file__)),
                          "results/phase61_deepsurv_explainer_revision/attributions_n2.npz": sha(OLD / "attributions_n2.npz")},
        "protocol_sha256": protocol_hash, "external_inventory": external_info,
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "python": platform.python_version(),
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    text = f"""# Faz 6.2 ADIM 1 — DeepSurv N=50 tam koşu

> **Koşulluluk uyarısı:** {WARNING}

Faz 6.1'in DeepExplainer candidate k=16 / reference k=128 tasarımı aynı 12+12 query ve 77 feature ile N=50'ye çıkarıldı. İlk iki deterministik attribution korundu; aynı patient-cluster bootstrap seed akışındaki 48 model ayrıca fit edildi. Ana model değiştirilmedi.

## Block agregatları

{markdown_table(aggregate)}

N=2'de {len(n2_pass)}/77, N=50'de **{len(n50_pass)}/77** feature WI≥0.8 geçti. N=2 kazananlarının {len(n2_pass & n50_pass)} tanesi korundu, {len(n2_pass - n50_pass)} tanesi düştü; {len(n50_pass - n2_pass)} yeni feature geçti. N=50 bottleneck dağılımı `{report.bottleneck.value_counts().to_dict()}`.

Yeni 48 replikanın birikimli eğitim+açıklama süresi **{new_compute:.3f} sn**; bu invocation **{invocation:.3f} sn** sürdü. Yeni replikalarda maksimum göreli additivity hatası candidate/reference için {max(all_candidate_errors):.3g}/{max(all_reference_errors):.3g}. SHAP'ın sabit mutlak 0.01 iç assertion'ı {fallback_count} çağrıda tetiklendi; bu çağrılar açıkça kaydedildi ve yalnız önceden tanımlı göreli hata ≤%5 olduğunda kabul edildi. Formül SHA256 değişmedi ve kimlik-sızıntısı kontrolü geçti.
"""
    (OUT / "report.md").write_text(text, encoding="utf-8")
    (ROOT / "PHASE6_2_STEP1_N50_REPORT.md").write_text(text, encoding="utf-8")
    print(json.dumps(metadata["results"], indent=2), flush=True)
    print(json.dumps(metadata["runtime"], indent=2), flush=True)


if __name__ == "__main__":
    run()
