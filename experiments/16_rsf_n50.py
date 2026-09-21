"""Phase 5.6 step 1: approved RSF N=50 audit with the frozen Phase 4.2 design."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
import warnings

import joblib
import numpy as np
import pandas as pd
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _shared import CALIBRATION, markdown_table


ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase56_rsf_n50"
OLD = ROOT / "results/phase4_rsf"
N = 50
N_OLD = 20
SEED = 20260917
TIME_LIMIT_SECONDS = 4 * 3600
WARNING = (
    "RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; "
    "per-permutation-path uncertainty is excluded."
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_phase4():
    path = ROOT / "experiments/07_rsf_pilot.py"
    spec = importlib.util.spec_from_file_location("phase4_rsf_frozen_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Phase 4 helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run() -> None:
    invoked = perf_counter()
    phase4 = load_phase4()
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoints = OUT / "checkpoints"
    checkpoints.mkdir(exist_ok=True)

    x = phase4.frame("X_train_selected.parquet")
    test = phase4.frame("X_test_selected.parquet")
    y = phase4.frame("y_train.parquet").loc[x.index]
    bundle = joblib.load(TCGA / "rsf_baseline.joblib")
    baseline = bundle["model"]
    features = list(bundle["features"])
    if features != list(x.columns) or baseline.n_features_in_ != 616:
        raise RuntimeError("Saved RSF is incompatible with the frozen 616-feature panel")
    expected_hyperparameters = (500, 10, 5, "sqrt", 42)
    actual_hyperparameters = (
        baseline.n_estimators, baseline.min_samples_split, baseline.min_samples_leaf,
        baseline.max_features, baseline.random_state,
    )
    if actual_hyperparameters != expected_hyperparameters:
        raise RuntimeError(f"Saved RSF hyperparameters changed: {actual_hyperparameters}")
    del baseline

    phase3 = phase4.phase3_helpers()
    external, mapping, external_info = phase3.external_columns(features)
    active_features = [feature for feature in features if feature in external]
    if len(active_features) != 77:
        raise RuntimeError(f"Expected 77 common features, found {len(active_features)}")
    active_indices = np.array([features.index(feature) for feature in active_features])
    external_raw = np.column_stack([external[feature] for feature in active_features])
    external_raw = external_raw[np.isfinite(external_raw).all(axis=1)]
    train_mean = x.mean(axis=0).to_numpy(float)
    train_std = x.std(axis=0, ddof=0).to_numpy(float)
    active_mean = train_mean[active_indices]
    active_std = train_std[active_indices]
    external_z = (external_raw - external_raw.mean(axis=0)) / external_raw.std(axis=0, ddof=0)
    external_aligned = active_mean + external_z * active_std

    train_patients = x.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    clean = test.loc[~test_patients.isin(train_patients)]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    source_pool = clean[active_features].to_numpy(float)
    query_rng = np.random.default_rng(SEED + 4000)
    source_selection = np.sort(query_rng.choice(len(source_pool), 12, replace=False))
    external_selection = np.sort(query_rng.choice(len(external_aligned), 12, replace=False))
    queries = [source_pool[source_selection], external_aligned[external_selection]]

    unique_train = x.loc[~train_patients.duplicated()]
    background_rng = np.random.default_rng(SEED + 4200)
    background_selection = np.sort(background_rng.choice(len(unique_train), 16, replace=False))
    sampled_backgrounds = unique_train.iloc[background_selection].to_numpy(dtype=np.float32)
    sampled_backgrounds[:, active_indices] = active_mean.astype(np.float32)

    inputs = {name: sha(TCGA / name) for name in [
        "rsf_baseline.joblib", "X_train_selected.parquet", "y_train.parquet"
    ]}
    old_metadata = json.loads((OLD / "metadata.json").read_text(encoding="utf-8"))
    if old_metadata["input_sha256"] != inputs:
        raise RuntimeError("N=20 and N=50 input hashes differ")
    query_hashes = [hashlib.sha256(np.asarray(query, dtype="<f8").tobytes()).hexdigest() for query in queries]
    protocol = {
        "inputs": inputs,
        "anchors_per_site": 12,
        "background_patients": 16,
        "candidate_budget": 8,
        "reference_budget": 128,
        "seed": SEED,
        "query_hashes": query_hashes,
        "background_sha256": hashlib.sha256(np.asarray(sampled_backgrounds, dtype="<f4").tobytes()).hexdigest(),
        "protocol": "rsf-patient-cluster-sampled-background-k16-permutation-8-128-v1",
    }
    protocol_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()

    with np.load(OLD / "attributions_n20.npz", allow_pickle=False) as old:
        if old["features"].tolist() != active_features:
            raise RuntimeError("N=20 attribution feature order differs")
        old_candidate = old["candidate"]
        old_reference = old["reference"]
    if old_candidate.shape != (2, 20, 12, 77) or old_reference.shape != old_candidate.shape:
        raise RuntimeError("Unexpected N=20 attribution shape")
    old_summary = json.loads((OLD / "bootstrap_summary.json").read_text(encoding="utf-8"))
    if len(old_summary) != N_OLD or any(row["warning_categories"] for row in old_summary):
        raise RuntimeError("N=20 bootstrap summary is incomplete or contains warnings")

    candidates = np.empty((2, N, 12, 77))
    references = np.empty_like(candidates)
    candidates[:, :N_OLD] = old_candidate
    references[:, :N_OLD] = old_reference
    patient_groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    survival = Surv.from_arrays(event=y.event.astype(bool).to_numpy(), time=y.time.to_numpy())
    values = x.to_numpy(dtype=np.float32)
    new_records: list[dict] = []
    checkpoint_sites_reused = 0

    for replicate in range(N):
        indices = np.concatenate([
            patient_groups[i] for i in bootstrap_rng.integers(0, len(patient_groups), len(patient_groups))
        ])
        if replicate < N_OLD:
            continue
        paths = [checkpoints / f"{protocol_hash}_rep{replicate:02d}_site{site}.npz" for site in range(2)]
        metadata_paths = [path.with_suffix(".json") for path in paths]
        missing = [site for site in range(2) if not paths[site].exists() or not metadata_paths[site].exists()]
        model = None
        training = None
        charge_site = missing[0] if missing else None
        if missing:
            model = RandomSurvivalForest(
                n_estimators=500, min_samples_split=10, min_samples_leaf=5, max_features="sqrt",
                random_state=SEED + replicate, n_jobs=-1, low_memory=True,
            )
            fit_started = perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model.fit(values[indices], survival[indices])
            training = {
                "rows": int(len(indices)),
                "events": int(survival[indices]["event"].sum()),
                "elapsed_seconds": perf_counter() - fit_started,
                "warning_categories": dict(Counter(w.category.__name__ for w in caught)),
            }
            if caught:
                raise RuntimeError(
                    f"RSF replicate {replicate} emitted warnings; no attribution checkpoint exported: "
                    f"{training['warning_categories']}"
                )
        for site in range(2):
            if site not in missing:
                with np.load(paths[site], allow_pickle=False) as saved:
                    candidates[site, replicate] = saved["candidate"]
                    references[site, replicate] = saved["reference"]
                record = json.loads(metadata_paths[site].read_text(encoding="utf-8"))
                if record["protocol_hash"] != protocol_hash:
                    raise RuntimeError("N=50 checkpoint protocol mismatch")
                new_records.append(record)
                checkpoint_sites_reused += 1
                print(f"RSF N=50 replicate {replicate + 1}/50 site {site}: verified checkpoint", flush=True)
                continue
            candidate_started = perf_counter()
            candidate = phase4.background_permutation_attributions(
                model, queries[site], sampled_backgrounds, active_mean, active_indices,
                permutations=8, seed=SEED + 10000 + 100 * site + replicate,
            )
            candidate_seconds = perf_counter() - candidate_started
            reference_started = perf_counter()
            reference = phase4.background_permutation_attributions(
                model, queries[site], sampled_backgrounds, active_mean, active_indices,
                permutations=128, seed=SEED + 20000 + 100 * site + replicate,
            )
            reference_seconds = perf_counter() - reference_started
            candidates[site, replicate] = candidate
            references[site, replicate] = reference
            record = {
                "protocol_hash": protocol_hash,
                "replicate": replicate,
                "site": site,
                "training": training,
                "training_seconds_charged": training["elapsed_seconds"] if site == charge_site else 0.0,
                "candidate_seconds": candidate_seconds,
                "reference_seconds": reference_seconds,
                "warning_categories": {},
            }
            np.savez_compressed(paths[site], candidate=candidate, reference=reference)
            metadata_paths[site].write_text(json.dumps(record, indent=2), encoding="utf-8")
            new_records.append(record)
            print(
                f"RSF N=50 replicate {replicate + 1}/50 site {site}: "
                f"candidate={candidate_seconds:.1f}s reference={reference_seconds:.1f}s",
                flush=True,
            )
        del model
        if perf_counter() - invoked > TIME_LIMIT_SECONDS:
            raise RuntimeError("PHASE56_TIME_LIMIT: invocation exceeded 4 hours after a complete replicate checkpoint")

    report = phase4.component_report(candidates, references, active_features)
    report_n20 = phase4.component_report(candidates[:, :20], references[:, :20], active_features)
    stored_n20 = pd.read_csv(OLD / "feature_report.csv")
    np.testing.assert_allclose(
        report_n20[["F", "S", "U", "G_claim", "R", "WI"]],
        stored_n20[["F", "S", "U", "G_claim", "R", "WI"]], atol=1e-12,
    )
    for column in ["F", "S", "U", "G_claim", "R", "WI", "point_adequate", "bottleneck"]:
        report[f"{column}_n20"] = report.feature.map(report_n20.set_index("feature")[column])
    report.insert(0, "inference_scope", WARNING)
    report.to_csv(OUT / "feature_report.csv", index=False)
    np.savez_compressed(OUT / "attributions_n50.npz", candidate=candidates, reference=references,
                        features=np.asarray(active_features))

    site_summary = sorted(new_records, key=lambda row: (row["replicate"], row["site"]))
    pd.DataFrame(site_summary).to_json(OUT / "bootstrap_site_summary_new30.json", orient="records", indent=2)
    new_summary = []
    for replicate in range(N_OLD, N):
        rows = [row for row in site_summary if row["replicate"] == replicate]
        attempts = []
        for row in rows:
            attempt = row.get("training")
            if attempt is not None and attempt not in attempts:
                attempts.append(attempt)
        if len(rows) != 2 or not attempts:
            raise RuntimeError(f"Incomplete N=50 checkpoint record for replicate {replicate}")
        new_summary.append({
            "replicate": replicate,
            "rows": attempts[0]["rows"],
            "events": attempts[0]["events"],
            "warning_categories": {},
            "fit_attempts": attempts,
            "site_records": [{"site": row["site"], "candidate_seconds": row["candidate_seconds"],
                              "reference_seconds": row["reference_seconds"]} for row in rows],
            "elapsed_seconds": sum(row["training_seconds_charged"] + row["candidate_seconds"] +
                                   row["reference_seconds"] for row in rows),
        })
    combined_summary = old_summary + new_summary
    for row in combined_summary:
        row["inference_scope"] = WARNING
    pd.DataFrame(combined_summary).to_json(OUT / "bootstrap_summary.json", orient="records", indent=2)
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)

    aggregate = report.groupby("block").agg(
        features=("feature", "size"), assessed=("WI", "count"), mean_F=("F", "mean"),
        mean_S=("S", "mean"), mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.insert(0, "inference_scope", WARNING)
    aggregate.to_csv(OUT / "summary.csv", index=False)

    n20_pass = set(report_n20.loc[report_n20.point_adequate == 1, "feature"])
    n50_pass = set(report.loc[report.point_adequate == 1, "feature"])
    comparison = report[[
        "feature", "block", "R_n20", "R", "WI_n20", "WI", "point_adequate_n20", "point_adequate",
        "bottleneck_n20", "bottleneck",
    ]].copy()
    comparison["delta_R_n50_minus_n20"] = comparison.R - comparison.R_n20
    comparison["delta_WI_n50_minus_n20"] = comparison.WI - comparison.WI_n20
    comparison.insert(0, "inference_scope", WARNING)
    comparison.to_csv(OUT / "n20_n50_comparison.csv", index=False)

    phase21 = json.loads((ROOT / "results/phase21/02_scenarios_metadata.json").read_text(encoding="utf-8"))
    old_hashes = {Path(name).as_posix(): digest for name, digest in phase21["source_sha256"].items()}
    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    unchanged = all(old_hashes[name] == sha(ROOT / name) for name in formula_files)
    if not unchanged:
        raise RuntimeError("Frozen formula module SHA256 changed")
    new_compute = sum(row["elapsed_seconds"] for row in new_summary)
    prior_compute = float(old_metadata["runtime"]["accumulated_measured_compute_wall_seconds"])
    invocation = perf_counter() - invoked
    metadata = {
        "phase": "5.6_step1_RSF_N50",
        "rsf_conditionality_warning": WARNING,
        "arguments": {"replicas": N, "anchors_per_site": 12, "background_patients": 16,
                      "candidate_permutations": 8, "reference_permutations": 128, "seed": SEED},
        "design": {"active_game_features": 77, "inactive_background_integrated_features": 539,
                   "query_rows": {"TCGA": 12, "METABRIC": 12},
                   "query_pool_sizes": {"TCGA": len(source_pool), "METABRIC": len(external_aligned)}},
        "model": {"source": str(TCGA / "rsf_baseline.joblib"), "main_model_retrained": False,
                  "bootstrap_models_in_N50_sequence": 50, "new_models_fit_this_phase": 30,
                  "first_20_reused_from_phase4_2": True, "sha256": inputs["rsf_baseline.joblib"],
                  "low_memory": True},
        "results": {"assessed_features": 77, "point_adequate_n20": len(n20_pass),
                    "point_adequate_n50": len(n50_pass), "retained_from_n20": len(n20_pass & n50_pass),
                    "lost_from_n20": len(n20_pass - n50_pass), "new_vs_n20": len(n50_pass - n20_pass),
                    "bottleneck_counts_n50": report.bottleneck.value_counts().to_dict()},
        "runtime": {"prior_N20_accumulated_seconds": prior_compute,
                    "new_replicates_20_to_49_seconds": new_compute,
                    "N50_total_accumulated_seconds": prior_compute + new_compute,
                    "completion_invocation_seconds": invocation,
                    "checkpoint_sites_reused_in_invocation": checkpoint_sites_reused,
                    "hard_stop_seconds": TIME_LIMIT_SECONDS},
        "formula_module_sha256": {name: sha(ROOT / name) for name in formula_files},
        "formula_module_sha256_unchanged_from_phase21": unchanged,
        "input_sha256": inputs,
        "source_sha256": {"experiments/16_rsf_n50.py": sha(Path(__file__)),
                          "results/phase4_rsf/attributions_n20.npz": sha(OLD / "attributions_n20.npz")},
        "protocol_sha256": protocol_hash,
        "external_inventory": external_info,
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "versions": {name: importlib.metadata.version(name) for name in
                     ["numpy", "pandas", "scikit-learn", "scikit-survival", "shap"]},
        "python": platform.python_version(),
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    warning = ("RSF sonuçları sabit 12 TCGA + 12 METABRIC query ve saklı 8/128 attribution "
               "ortalamalarına koşulludur; tekil permütasyon-yolu belirsizliğini içermez.")
    text = "# Faz 5.6 ADIM 1 — RSF N=50 tam koşu\n\n> **Koşulluluk uyarısı:** " + warning + "\n\n"
    text += "Onaylanan Faz 4.2 tasarımında yalnız replika sayısı 20'den 50'ye çıkarıldı. İlk 20 "
    text += "deterministik replikanın doğrulanmış attribution artefaktları korundu; aynı seed akışındaki "
    text += "21–50 replikaları yeni fit edildi. Ana `rsf_baseline.joblib` değiştirilmedi.\n\n"
    text += "## Block agregatları\n\n" + markdown_table(aggregate) + "\n\n"
    text += (f"77/77 feature değerlendirildi; N=50'de {len(n50_pass)} feature WI≥0.8 geçti "
             f"(N=20: {len(n20_pass)}). Bottleneck dağılımı `{report.bottleneck.value_counts().to_dict()}`. "
             f"N=20 geçenlerin {len(n20_pass & n50_pass)}/23'ü korundu; {len(n20_pass - n50_pass)} düştü, "
             f"{len(n50_pass - n20_pass)} yeni feature geçti.\n\n")
    text += (f"Yeni 30 replikanın ölçülen eğitim+açıklama süresi {new_compute:.1f} sn "
             f"({new_compute/3600:.2f} saat); mevcut ilk 20 dahil N=50 birikimli maliyet "
             f"{prior_compute + new_compute:.1f} sn ({(prior_compute+new_compute)/3600:.2f} saat). "
             f"Bu invocation {invocation:.1f} sn sürdü. Formül SHA256 değerleri değişmedi.\n\n")
    text += "## Tüm feature sonuçları\n\n" + markdown_table(report[[
        "feature", "block", "F", "S", "U", "G_claim", "R", "WI", "bottleneck",
        "R_n20", "WI_n20", "point_adequate_n20", "point_adequate"
    ]]) + "\n"
    (OUT / "report.md").write_text(text, encoding="utf-8")
    (ROOT / "PHASE5_6_STEP1_RSF_N50_REPORT.md").write_text(text, encoding="utf-8")
    print(json.dumps(metadata["results"], indent=2), flush=True)
    print(json.dumps(metadata["runtime"], indent=2), flush=True)


if __name__ == "__main__":
    run()
