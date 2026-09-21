#!/usr/bin/env python3
"""Phase 7.2 frozen-design RSF N=50 audit on the clean fold-0 split."""

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
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results/phase7_retrain"
PANEL = ROOT / "results/phase72_panel"
OUT = ROOT / "results/phase72_rsf"
CHECKPOINTS = OUT / "checkpoints"
N, SEED = 50, 20260917


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def frame(name: str) -> pd.DataFrame:
    return pd.read_parquet(DATA / name).set_index("sample_key")


def load_pilot():
    path = ROOT / "experiments/07_rsf_pilot.py"
    spec = importlib.util.spec_from_file_location("phase72_rsf_helpers", path)
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
    if formula != json.loads((PANEL / "metadata.json").read_text())["formula_sha256"]:
        raise RuntimeError("Frozen formula hash changed")

    train = frame("X_train_selected.parquet")
    y = frame("y_train.parquet").loc[train.index]
    saved = joblib.load(DATA / "rsf_baseline.joblib")
    features = list(saved["features"])
    baseline = saved["model"]
    if features != list(train.columns):
        raise RuntimeError("RSF feature order differs from the clean matrix")
    frozen = (
        baseline.n_estimators == 500
        and baseline.min_samples_split == 10
        and baseline.min_samples_leaf == 5
        and baseline.max_features == "sqrt"
        and baseline.random_state == 42
    )
    if not frozen:
        raise RuntimeError("Frozen RSF hyperparameters changed")
    del baseline

    with np.load(PANEL / "query_matrices.npz", allow_pickle=False) as panel:
        active = panel["active_features"].astype(str).tolist()
        queries = [panel["source_query12"], panel["external_query12"]]
    active_indices = np.array([features.index(name) for name in active], dtype=int)
    active_mean = train.mean(axis=0).to_numpy(float)[active_indices]
    patients = train.index.astype(str).str[:12]
    unique_train = train.loc[~patients.duplicated()]
    background_rng = np.random.default_rng(SEED + 4200)
    selection = np.sort(background_rng.choice(len(unique_train), 16, replace=False))
    backgrounds = unique_train.iloc[selection].to_numpy(dtype=np.float32)
    backgrounds[:, active_indices] = active_mean.astype(np.float32)

    inputs = {
        name: sha(DATA / name)
        for name in ["rsf_baseline.joblib", "X_train_selected.parquet", "y_train.parquet"]
    }
    protocol_payload = {
        "inputs": inputs,
        "panel": sha(PANEL / "query_matrices.npz"),
        "N": N,
        "seed": SEED,
        "background_k": 16,
        "candidate": 8,
        "reference": 128,
        "protocol": "patient-cluster-rsf-sampled-background-permutation-v2",
    }
    protocol = hashlib.sha256(json.dumps(protocol_payload, sort_keys=True).encode()).hexdigest()
    groups = [np.flatnonzero(patients == patient) for patient in patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    survival = Surv.from_arrays(event=y.event.astype(bool).to_numpy(), time=y.time.to_numpy())
    values = train.to_numpy(dtype=np.float32)
    candidate = np.empty((2, N, 12, len(active)))
    reference = np.empty_like(candidate)
    records: list[dict] = []
    helpers = load_pilot()

    for rep in range(N):
        indices = np.concatenate(
            [groups[index] for index in bootstrap_rng.integers(0, len(groups), len(groups))]
        )
        paths = [CHECKPOINTS / f"{protocol}_rep{rep:02d}_site{site}.npz" for site in range(2)]
        metas = [path.with_suffix(".json") for path in paths]
        missing = [site for site in range(2) if not paths[site].exists() or not metas[site].exists()]
        model = None
        training = None
        charge = missing[0] if missing else None
        if missing:
            model = RandomSurvivalForest(
                n_estimators=500,
                min_samples_split=10,
                min_samples_leaf=5,
                max_features="sqrt",
                random_state=SEED + rep,
                n_jobs=-1,
                low_memory=True,
            )
            fit_started = perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model.fit(values[indices], survival[indices])
            training = {
                "rows": len(indices),
                "events": int(survival[indices]["event"].sum()),
                "elapsed_seconds": perf_counter() - fit_started,
                "warning_categories": dict(Counter(item.category.__name__ for item in caught)),
            }
            if caught:
                raise RuntimeError(f"RSF replicate {rep} emitted warnings: {training['warning_categories']}")
        for site in range(2):
            if site not in missing:
                meta = json.loads(metas[site].read_text())
                if meta["protocol_hash"] != protocol:
                    raise RuntimeError("RSF checkpoint protocol mismatch")
                with np.load(paths[site], allow_pickle=False) as saved_site:
                    candidate[site, rep] = saved_site["candidate"]
                    reference[site, rep] = saved_site["reference"]
                records.append(meta)
                print(f"RSF 7.2 replicate {rep + 1}/50 site {site}: verified checkpoint", flush=True)
                continue
            c_started = perf_counter()
            c = helpers.background_permutation_attributions(
                model, queries[site], backgrounds, active_mean, active_indices,
                permutations=8, seed=SEED + 10000 + 100 * site + rep,
            )
            c_seconds = perf_counter() - c_started
            r_started = perf_counter()
            r = helpers.background_permutation_attributions(
                model, queries[site], backgrounds, active_mean, active_indices,
                permutations=128, seed=SEED + 20000 + 100 * site + rep,
            )
            r_seconds = perf_counter() - r_started
            candidate[site, rep], reference[site, rep] = c, r
            meta = {
                "protocol_hash": protocol,
                "replicate": rep,
                "site": site,
                "training": training,
                "training_seconds_charged": training["elapsed_seconds"] if site == charge else 0.0,
                "candidate_seconds": c_seconds,
                "reference_seconds": r_seconds,
                "warning_categories": {},
            }
            np.savez_compressed(paths[site], candidate=c, reference=r)
            metas[site].write_text(json.dumps(meta, indent=2), encoding="utf-8")
            records.append(meta)
            print(
                f"RSF 7.2 replicate {rep + 1}/50 site {site}: "
                f"candidate={c_seconds:.1f}s reference={r_seconds:.1f}s",
                flush=True,
            )
        del model

    report = helpers.component_report(candidate, reference, active)
    report.to_csv(OUT / "feature_report.csv", index=False)
    np.savez_compressed(
        OUT / "attributions_n50.npz", candidate=candidate, reference=reference,
        features=np.asarray(active),
    )
    ordered = sorted(records, key=lambda row: (row["replicate"], row["site"]))
    bootstrap = []
    for rep in range(N):
        rows = [row for row in ordered if row["replicate"] == rep]
        attempts = []
        for row in rows:
            if row.get("training") is not None and row["training"] not in attempts:
                attempts.append(row["training"])
        bootstrap.append({
            "replicate": rep,
            "rows": attempts[0]["rows"],
            "events": attempts[0]["events"],
            "warning_categories": {},
            "fit_attempts": attempts,
            "site_records": [
                {key: row[key] for key in ["site", "candidate_seconds", "reference_seconds"]}
                for row in rows
            ],
            "elapsed_seconds": sum(
                row["training_seconds_charged"] + row["candidate_seconds"] + row["reference_seconds"]
                for row in rows
            ),
        })
    (OUT / "bootstrap_summary.json").write_text(json.dumps(bootstrap, indent=2), encoding="utf-8")
    aggregate = report.groupby("block").agg(
        assessed=("WI", "count"), mean_F=("F", "mean"), mean_S=("S", "mean"),
        mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.to_csv(OUT / "summary.csv", index=False)
    metadata = {
        "phase": "7.2_rsf_n50_clean",
        "arguments": protocol_payload,
        "point_adequate": int(report.point_adequate.sum()),
        "bottleneck_counts": report.bottleneck.value_counts().to_dict(),
        "formula_sha256": formula,
        "privacy": {"identifiers_exported": False, "patient_overlap": 0},
        "conditionality": "fixed deterministic 12+12 queries; unsaved permutation paths; sampled k=16 inactive background",
        "accumulated_compute_seconds": sum(row["elapsed_seconds"] for row in bootstrap),
        "invocation_wall_clock_seconds": perf_counter() - invoked,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    table = helpers.markdown_table(aggregate)
    note = f"""# Faz 7.2 — RSF N=50 temiz audit

Hasta-gruplu fold-0 referans RSF modeli ve temiz 27-feature TCGA–METABRIC paneli üzerinde 50 patient-cluster bootstrap replika tamamlandı. **{metadata['point_adequate']}/27** feature nokta tahminde WI≥0.8 geçti.

{table}

Bottleneck dağılımı: `{metadata['bottleneck_counts']}`.

Bu sonuçlar sabit deterministik 12+12 query kümesine, saklanmamış permütasyon yollarına, sampled inactive background k=16 ve 8/128 permütasyon bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen hesap süresi **{metadata['accumulated_compute_seconds']:.1f} sn**, bu invocation wall-clock **{metadata['invocation_wall_clock_seconds']:.1f} sn**.
"""
    (ROOT / "PHASE72_RSF.md").write_text(note, encoding="utf-8")
    print(json.dumps({"point_adequate": metadata["point_adequate"], "bottlenecks": metadata["bottleneck_counts"]}, indent=2))


if __name__ == "__main__":
    main()
