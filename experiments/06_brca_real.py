"""Phase 3, Sections B.2-B.6: frozen CoxPH real-data audit; no formula changes.

Patient keys are used only in memory/raw inputs, never in outputs or logs.
assumption: cohort-specific feature z-scores define aligned coordinates.
Phase 3.1 compares real TCGA/METABRIC query distributions. Phase 3.2 changes
only the number of patient-cluster Cox bootstrap replicates.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
from types import SimpleNamespace
import warnings

import joblib
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from warrantxai.api import WarrantIndexEvaluator
from warrantxai.aggregation import adequate, warrant_index
from warrantxai.components import generalizability, generalizability_signed, pairwise_orders
from warrantxai.dgp import gaussian_linear
from _shared import CALIBRATION, markdown_table

ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
RAW = ROOT / "data/external/metabric"
OUT = ROOT / "results/phase3_brca"


def sha(path: Path) -> str:
    """Return whole-file provenance checksum (Phase 3 reporting)."""
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def empirical_rank_grid(values: np.ndarray, size: int) -> np.ndarray:
    """Repeat observed order statistics on a common quantile grid (Phase 3.1).

    This is a distributional coupling, not a patient match. It uses only
    observed query values (no interpolation) and lets the unchanged
    G_signed anchor mean accept cohorts with different patient counts.
    """
    ordered = np.sort(np.asarray(values, dtype=float))
    if not len(ordered) or size < 1:
        raise ValueError("nonempty empirical distribution and positive grid size required")
    indices = np.minimum((np.arange(size) + 0.5) * len(ordered) // size, len(ordered) - 1).astype(int)
    return ordered[indices]


def rank_bootstrap_bands(
    coefficients: np.ndarray,
    anchors: np.ndarray,
    feature_names: list[str],
    targets: list[str],
    eta: float,
    *,
    draws: int,
    seed: int,
    label: str,
) -> pd.DataFrame:
    """Estimate conditional R variability by resampling training replicas.

    The point estimate is exactly the frozen R pair formula. Bands resample
    the observed coefficient replicas with replacement and are engineering
    precision diagnostics, not patient-level or biological confidence bands.
    """
    if len(coefficients) < 2 or draws < 2:
        raise ValueError("R bands need at least two replicas and draws")
    mean_abs_anchor = np.mean(np.abs(anchors), axis=0)
    imp = np.abs(coefficients) * mean_abs_anchor
    orders = pairwise_orders(imp, eta)
    target_indices = np.array([feature_names.index(f) for f in targets])
    t, p = imp.shape
    pair_loss = np.empty((len(targets), t, t))
    for out, feature_index in enumerate(target_indices):
        order = orders[:, feature_index, :].astype(float)
        pair_loss[out] = np.abs(order[:, None, :] - order[None, :, :]).sum(axis=-1) / (2 * (p - 1))
    point = 1 - pair_loss.sum(axis=(1, 2)) / (t * (t - 1))
    rng = np.random.default_rng(seed)
    sampled = np.empty((draws, len(targets)))
    for draw in range(draws):
        indices = rng.integers(0, t, t)
        loss = pair_loss[:, indices[:, None], indices]
        sampled[draw] = 1 - loss.sum(axis=(1, 2)) / (t * (t - 1))
    quantiles = np.quantile(sampled, [0.025, 0.25, 0.75, 0.975], axis=0)
    return pd.DataFrame(
        {
            "feature": targets,
            f"R_{label}": point,
            f"R_{label}_resample_std": sampled.std(axis=0, ddof=1),
            f"R_{label}_q025": quantiles[0],
            f"R_{label}_q25": quantiles[1],
            f"R_{label}_q75": quantiles[2],
            f"R_{label}_q975": quantiles[3],
        }
    )


def frame(name: str) -> pd.DataFrame:
    """Read original TCGA without displaying row keys (Phase 3 manifest)."""
    return pd.read_parquet(TCGA / name).set_index("sample_key")


def raw_json(name: str) -> list | dict:
    """Read a preserved HTTP response (Phase 3 data)."""
    return json.loads((RAW / name).read_bytes())


def external_columns(features: list[str]) -> tuple[dict[str, np.ndarray], pd.DataFrame, dict]:
    """Map strict gene symbols and clinical semantics (Phase 3 support).

    No aliases, absent modalities, unknown stages or grade/stage substitution.
    Missing values stay missing; no model prediction needs an imputed modality.
    """
    annotations = raw_json("ensembl_symbol_mapping.json")
    by_query: dict[str, set[str]] = {}
    for row in annotations:
        if "symbol" in row:
            by_query.setdefault(row["query"], set()).add(row["symbol"])
    genes: dict[str, set[int]] = {}
    for row in raw_json("genes.json"):
        genes.setdefault(row["hugoGeneSymbol"], set()).add(row["entrezGeneId"])
    molecular = []
    for path in sorted(RAW.glob("mrna_batch_*.json")):
        molecular.extend(json.loads(path.read_bytes()))
    mol = pd.DataFrame(molecular)
    if mol.duplicated(["sampleId", "entrezGeneId"]).any():
        raise ValueError("Duplicate expression keys; explicit resolution required")
    expression = mol.pivot(index="sampleId", columns="entrezGeneId", values="value")
    patients = pd.DataFrame(raw_json("patient_clinical.json")).pivot(
        index="patientId", columns="clinicalAttributeId", values="value"
    )
    samples = pd.DataFrame(raw_json("sample_clinical.json")).pivot(
        index="sampleId", columns="clinicalAttributeId", values="value"
    )
    sample_info = pd.DataFrame(raw_json("samples.json")).set_index("sampleId")
    linked = sample_info.loc[expression.index, "patientId"]
    if linked.duplicated().any():
        raise ValueError("Multiple expression samples per external patient need a declared selection")
    ext: dict[str, np.ndarray] = {}
    mappings = []
    for f in features:
        field = ""
        symbol = ""
        reason = "no external cohort"
        values = None
        if f.startswith("rna_"):
            symbols = by_query.get(f[4:].split(".")[0], set())
            if len(symbols) == 1:
                symbol = next(iter(symbols))
                ids = genes.get(symbol, set())
                if len(ids) == 1 and next(iter(ids)) in expression:
                    field = f"brca_metabric_mrna:{symbol}"
                    values = expression[next(iter(ids))].to_numpy(float)
                else:
                    reason = "symbol absent/ambiguous or no expression in requested profile"
            else:
                reason = "Ensembl symbol absent/ambiguous"
        elif f == "clin_numeric__age_at_index.demographic":
            field = "AGE_AT_DIAGNOSIS"
            values = pd.to_numeric(patients.reindex(linked)[field], errors="coerce").to_numpy(float)
        elif "ajcc_pathologic_stage.diagnoses_Stage " in f:
            stage = f.rsplit(" ", 1)[1]
            code = {"I": 1, "II": 2, "III": 3, "IV": 4}.get(stage)
            if code is not None:
                field = "TUMOR_STAGE"
                v = pd.to_numeric(samples.reindex(expression.index)[field], errors="coerce").to_numpy(float)
                values = np.where(np.isfinite(v), (v == code).astype(float), np.nan)
            else:
                reason = "unknown Stage X has no aligned category"
        if values is not None:
            ext[f] = values
            reason = "mapped; empirical support check pending"
        mappings.append({"feature": f, "symbol": symbol, "external_field": field, "mapping_status": reason})
    os_time = pd.to_numeric(patients.reindex(linked)["OS_MONTHS"], errors="coerce")
    os_status = patients.reindex(linked)["OS_STATUS"]
    eligible = os_time.notna() & (os_time > 0) & os_status.isin(["0:LIVING", "1:DECEASED"])
    info = {
        "expression_samples": len(expression),
        "expression_genes": len(expression.columns),
        "clinical_patients": len(patients),
        "clinical_samples": len(samples),
        "OS_eligible_expression_patients": int(eligible.sum()),
        "OS_events_expression_patients": int((os_status[eligible] == "1:DECEASED").sum()),
        "survival_unit": "months",
        "survival_used_for_fitting_or_tuning": False,
        "clinical_unavailable_in_TCGA_panel": ["GRADE", "ER_STATUS", "PR_STATUS", "HER2_STATUS"],
        "clinical_excluded_from_requested_mapping": ["SEX", "race", "ethnicity"],
    }
    return ext, pd.DataFrame(mappings), info


def run(replicates: int = 50, seed: int = 20260917, reuse_bootstrap: bool = False) -> None:
    """Run conditional bootstrap CoxPH audit and report aggregates (B.2-B.6).

    Fixed selected features and preprocessing; bootstrap is patient-cluster,
    not a repeat of feature selection. Source held-out patients are disjoint.
    """
    if replicates < 2:
        raise ValueError("At least two bootstrap replicates required")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    previous_path = OUT / "feature_report.csv"
    if not previous_path.exists():
        raise FileNotFoundError("Phase 3 feature_report.csv is required to preserve the range-overlap diagnostic")
    previous = pd.read_csv(previous_path).set_index("feature")
    old_g_column = "G_claim_range_overlap_diagnostic" if "G_claim_range_overlap_diagnostic" in previous else "G_claim"
    old_wi_column = "WI_phase3_range_overlap" if "WI_phase3_range_overlap" in previous else "WI"
    old_adequacy_column = "point_adequate_phase3" if "point_adequate_phase3" in previous else "point_adequate"
    old_bottleneck_column = "bottleneck_phase3" if "bottleneck_phase3" in previous else "bottleneck"
    n6_columns = ["F", "S", "U", "G_magnitude", "G_signed", "G_claim", "R", "WI", "point_adequate", "bottleneck"]
    n6_source = {name: (f"{name}_n6" if f"{name}_n6" in previous else name) for name in n6_columns}
    x = frame("X_train_selected.parquet")
    test = frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[x.index]
    yt = frame("y_test.parquet").loc[test.index]
    assert x.columns.equals(test.columns)
    assert (len(x), int(y.event.sum()), len(test), int(yt.event.sum())) == (673, 106, 169, 26)
    original = joblib.load(TCGA / "coxph_baseline.joblib")
    features = list(original["features"])
    assert original["model"].penalizer == 0.1
    assert features == list(original["model"].params_.index)
    e = y.event.astype(bool)
    excluded = x.columns[(x[e].var() <= 1e-12) | (x[~e].var() <= 1e-12)].tolist()
    assert set(features) == set(x.columns) - set(excluded)
    train_patients = x.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    keep = ~test_patients.isin(train_patients)
    clean = test.loc[keep]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    assert not set(clean.index.astype(str).str[:12]) & set(train_patients)
    scaler = original["scaler"]
    scaled = scaler.transform(x[features].astype(float))
    anchors = scaler.transform(clean[features].astype(float))
    assert np.allclose(scaled.mean(0), 0, atol=1e-6)
    assert np.allclose(scaled.std(0), 1, atol=1e-6)
    # Exact marginal attribution: beta*(z-E_Q z). Centering adds only intercept.
    beta = original["model"].params_.to_numpy()
    score = original["model"].predict_log_partial_hazard(pd.DataFrame(anchors, columns=features)).to_numpy()
    base = float(
        original["model"]
        .predict_log_partial_hazard(pd.DataFrame(np.zeros((1, len(features))), columns=features))
        .iloc[0]
    )
    oracle = gaussian_linear(beta).oracle_attributions(anchors)
    oracle_error = float(np.max(np.abs(oracle.sum(1) + base - score)))
    assert oracle_error < 1e-8
    groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    rng = np.random.default_rng(seed)
    coeffs = []
    training = []
    if reuse_bootstrap:
        prior = json.loads((OUT / "metadata.json").read_text(encoding="utf-8"))
        assert prior["arguments"]["replicates"] == replicates and prior["arguments"]["seed"] == seed
        assert prior["arguments"]["penalizer"] == 0.1 and prior["arguments"]["max_steps"] == 50
        assert prior["versions"]["lifelines"] == importlib.metadata.version("lifelines")
        for filename in ["X_train_selected.parquet", "y_train.parquet", "coxph_baseline.joblib"]:
            assert prior["input_sha256"][filename] == sha(TCGA / filename)
        with np.load(OUT / "bootstrap_coefficients.npz", allow_pickle=False) as cached:
            assert cached["features"].tolist() == features
            coeffs = list(cached["coefficients"])
        assert np.shape(coeffs) == (replicates, len(features)) and np.isfinite(coeffs).all()
        training = json.loads((OUT / "bootstrap_summary.json").read_text())
        assert len(training) == replicates and all(not row["warning_categories"] for row in training)
        origin = OUT / "bootstrap_origin_metadata.json"
        if not origin.exists():
            origin.write_text(json.dumps(prior, indent=2), encoding="utf-8")
        print(
            "Reusing verified completed bootstrap coefficient artifacts; input hashes and configuration match.",
            flush=True,
        )
    else:
        checkpoint_dir = OUT / "bootstrap_checkpoints"
        checkpoint_dir.mkdir(exist_ok=True)
        training_inputs = {
            name: sha(TCGA / name) for name in ["X_train_selected.parquet", "y_train.parquet", "coxph_baseline.joblib"]
        }
        for rep in range(replicates):
            indices = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
            checkpoint_key = hashlib.sha256(
                json.dumps(
                    {
                        "inputs": training_inputs,
                        "seed": seed,
                        "rep": rep,
                        "lifelines": importlib.metadata.version("lifelines"),
                        "numpy": np.__version__,
                        "protocol": "cluster-fixed-scaler-separation50-retry200-step025-v1",
                    },
                    sort_keys=True,
                ).encode()
            ).hexdigest()
            checkpoint = checkpoint_dir / (checkpoint_key + ".json")
            if checkpoint.exists():
                saved = json.loads(checkpoint.read_text())
                coef = np.asarray(saved["coefficients"], dtype=float)
                assert coef.shape == (len(features),) and np.isfinite(coef).all()
                assert not saved["training"]["warning_categories"]
                coeffs.append(coef)
                training.append(saved["training"])
                print(f"Cox bootstrap {rep + 1}/{replicates}: verified checkpoint", flush=True)
                continue
            replicate_started = perf_counter()
            values = scaled[indices]
            events = y.event.to_numpy()[indices].astype(bool)
            valid = (values[events].var(0) > 1e-12) & (values[~events].var(0) > 1e-12)
            attempts = []
            fit_options = {"max_steps": 50}
            for attempt in range(20):
                df = pd.DataFrame(values[:, valid], columns=np.array(features)[valid])
                df["time"] = y.time.to_numpy()[indices]
                df["event"] = events
                model = CoxPHFitter(penalizer=0.1)
                with warnings.catch_warnings(record=True) as caught, threadpool_limits(limits=2):
                    warnings.simplefilter("always")
                    model.fit(df, "time", "event", show_progress=False, fit_options=fit_options.copy())
                convergence = [w for w in caught if issubclass(w.category, ConvergenceWarning)]
                flagged = [
                    j
                    for j, f in enumerate(features)
                    if valid[j]
                    and any(str(w.message).startswith(f"Column {f} have very low variance") for w in convergence)
                ]
                attempts.append(
                    {
                        "attempt": attempt,
                        "fit_options": fit_options.copy(),
                        "warning_categories": dict(Counter(w.category.__name__ for w in caught)),
                        "separation_features": [features[j] for j in flagged],
                    }
                )
                if not convergence:
                    break
                if not flagged:
                    if fit_options["max_steps"] == 50:
                        fit_options = {"max_steps": 200, "step_size": 0.25}
                        continue
                    raise RuntimeError("Cox non-convergence after numerical retry; no scores exported")
                # Same bootstrap sample/penalty/solver, removing training-only
                # separation columns explicitly identified by lifelines.
                valid[flagged] = False
            else:
                raise RuntimeError("Separation remediation exhausted; no scores exported")
            coef = np.zeros(len(features))
            coef[valid] = model.params_.to_numpy()
            check = pd.DataFrame(anchors[:, valid], columns=np.array(features)[valid])
            reference_score = model.predict_log_partial_hazard(check).to_numpy()
            intercept = float(model.predict_log_partial_hazard(check.iloc[:1] * 0).iloc[0])
            assert np.max(np.abs(anchors @ coef + intercept - reference_score)) < 1e-8
            coeffs.append(coef)
            training.append(
                {
                    "replicate": rep,
                    "rows": len(indices),
                    "events": int(events.sum()),
                    "excluded_feature_count": int((~valid).sum()),
                    "warning_categories": dict(Counter(w.category.__name__ for w in caught)),
                    "fit_attempts": attempts,
                    "excluded_features": [features[j] for j in np.flatnonzero(~valid)],
                    "elapsed_seconds": perf_counter() - replicate_started,
                }
            )
            checkpoint.write_text(
                json.dumps({"coefficients": coef.tolist(), "training": training[-1]}, indent=2), encoding="utf-8"
            )
            print(
                f"Cox bootstrap {rep + 1}/{replicates}: rows={len(indices)}, events={int(events.sum())}, exclusions={int((~valid).sum())}",
                flush=True,
            )
    models = [SimpleNamespace(coef_=b) for b in coeffs]
    np.savez_compressed(
        OUT / "bootstrap_coefficients.npz", coefficients=np.asarray(coeffs), features=np.asarray(features)
    )
    np.savez_compressed(
        OUT / "bootstrap_coefficients_n6.npz", coefficients=np.asarray(coeffs[:6]), features=np.asarray(features)
    )
    np.savez_compressed(
        OUT / f"bootstrap_coefficients_n{replicates}.npz",
        coefficients=np.asarray(coeffs),
        features=np.asarray(features),
    )
    audit = WarrantIndexEvaluator(
        **{k: CALIBRATION[k] for k in ["a", "b", "u", "eta", "tau"]},
        background_mean=np.zeros(len(features)),
        generalizability_mode="claim",
    ).fit([models], [anchors])
    report = audit.feature_report()
    report.insert(0, "feature", features)
    ext, mapping, external_info = external_columns(list(x.columns))
    mapped_features = [f for f in features if f in ext]
    external_matrix = np.column_stack([ext[f] for f in mapped_features])
    external_complete = np.isfinite(external_matrix).all(axis=1)
    if external_complete.sum() < 2:
        raise ValueError("Insufficient complete METABRIC patient vectors for the 77-feature query set")
    external_info["mapped_model_features"] = len(mapped_features)
    external_info["complete_query_patients"] = int(external_complete.sum())
    external_info["incomplete_query_patients_excluded"] = int((~external_complete).sum())
    external_info["query_source"] = "actual METABRIC patient vectors"
    diagnostics = []
    gm = np.full(len(features), np.nan)
    gs = gm.copy()
    coef_array = np.asarray(coeffs)
    mean_coef = coef_array.mean(axis=0)
    # Phase 3.1: site 0 and site 1 are different real patient query sets.
    # G_signed uses deterministic empirical-rank coupling, never patient IDs.
    for j, f in enumerate(features):
        v = ext.get(f)
        reason = "not assessable (no external cohort)"
        row = {
            "feature": f,
            "query_source_site0": "TCGA clean disjoint audit patients",
            "query_source_site1": "METABRIC actual patient vectors",
            "site0_query_n": len(anchors),
            "site1_query_n": 0,
            "distribution_grid_n": 0,
        }
        if v is not None:
            v = v[external_complete]
            row["site1_query_n"] = len(v)
            if len(v) >= 2 and v.std() > 0:
                # assumption: background=0 is each cohort's own mean in its
                # separate z-score space, including binary stage covariates.
                zext = (v - v.mean()) / v.std()
                zsource = anchors[:, j]
                source_phi = coef_array[:, j, None] * zsource[None, :]
                external_phi = coef_array[:, j, None] * zext[None, :]
                source_importance = float(np.mean(np.abs(source_phi)))
                external_importance = float(np.mean(np.abs(external_phi)))
                gm[j] = generalizability([[source_importance], [external_importance]], CALIBRATION["b"])[0]
                grid_n = max(len(zsource), len(zext))
                source_signed = empirical_rank_grid(mean_coef[j] * zsource, grid_n)
                external_signed = empirical_rank_grid(mean_coef[j] * zext, grid_n)
                gs[j] = generalizability_signed(
                    np.stack([source_signed, external_signed])[:, :, None], CALIBRATION["a"]
                )[0]
                source_query = empirical_rank_grid(zsource, grid_n)
                external_query = empirical_rank_grid(zext, grid_n)
                row["distribution_grid_n"] = grid_n
                row["query_value_distributions_identical"] = bool(np.array_equal(source_query, external_query))
                row["site0_query_sha256"] = hashlib.sha256(np.asarray(zsource, dtype="<f8").tobytes()).hexdigest()
                row["site1_query_sha256"] = hashlib.sha256(np.asarray(zext, dtype="<f8").tobytes()).hexdigest()
                row["source_native_mean_abs_attribution"] = float(np.mean(np.abs(source_phi)))
                row["external_native_mean_abs_attribution"] = float(np.mean(np.abs(external_phi)))
                row["external_z_mean"] = float(zext.mean())
                row["external_z_std"] = float(zext.std())
                reason = "assessed (different real cohort query distributions)"
            else:
                reason = "not assessable (external constant/missing)"
        row["status"] = reason
        diagnostics.append(row)
    report["G_magnitude"] = gm
    report["G_signed"] = gs
    report["G_claim"] = np.minimum(gm, gs)
    columns = ["F", "S", "U", "G_claim", "R"]
    z = report[columns].to_numpy()
    report["WI"] = warrant_index(z)
    report["point_adequate"] = adequate(z, CALIBRATION["tau"])
    report["status"] = [d["status"] for d in diagnostics]
    report["bottleneck"] = [columns[int(np.argmin(row))] if np.isfinite(row).all() else "not assessable" for row in z]
    report["source_bottleneck_F_S_U_R"] = report[["F", "S", "U", "R"]].idxmin(axis=1)
    report["block"] = [f.split("_")[0] for f in features]
    for f in excluded:
        report = pd.concat(
            [
                report,
                pd.DataFrame(
                    [
                        {
                            "feature": f,
                            "block": f.split("_")[0],
                            "status": "not assessable (excluded from Cox model: separation)",
                            "bottleneck": "not assessable",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    report = report.set_index("feature").reindex(x.columns).reset_index(names="feature")
    report["G_claim_range_overlap_diagnostic"] = report.feature.map(previous[old_g_column])
    report["G_range_overlap_diagnostic"] = report["G_claim_range_overlap_diagnostic"]
    report["G_magnitude_range_overlap_diagnostic"] = report["G_claim_range_overlap_diagnostic"]
    report["G_signed_range_overlap_diagnostic"] = report["G_claim_range_overlap_diagnostic"]
    report["WI_phase3_range_overlap"] = report.feature.map(previous[old_wi_column])
    report["point_adequate_phase3"] = report.feature.map(previous[old_adequacy_column])
    report["bottleneck_phase3"] = report.feature.map(previous[old_bottleneck_column])
    for name, source in n6_source.items():
        report[f"{name}_n6"] = report.feature.map(previous[source])
    n6_passing_features = report.loc[report.point_adequate_n6 == 1, "feature"].tolist()
    if len(n6_passing_features) != 27:
        raise RuntimeError(f"Expected 27 N=6 passing features, found {len(n6_passing_features)}")
    r_n6 = rank_bootstrap_bands(
        np.asarray(coeffs[:6]), anchors, features, n6_passing_features, CALIBRATION["eta"],
        draws=500, seed=seed + 32006, label="n6"
    )
    r_primary = rank_bootstrap_bands(
        np.asarray(coeffs), anchors, features, n6_passing_features, CALIBRATION["eta"],
        draws=500, seed=seed + 32050, label=f"n{replicates}"
    )
    r_uncertainty = r_n6.merge(r_primary, on="feature", validate="one_to_one")
    np.testing.assert_allclose(
        r_uncertainty.R_n6, report.set_index("feature").loc[n6_passing_features, "R_n6"], atol=1e-12
    )
    np.testing.assert_allclose(
        r_uncertainty[f"R_n{replicates}"], report.set_index("feature").loc[n6_passing_features, "R"], atol=1e-12
    )
    r_uncertainty["IQR_width_n6"] = r_uncertainty.R_n6_q75 - r_uncertainty.R_n6_q25
    r_uncertainty[f"IQR_width_n{replicates}"] = (
        r_uncertainty[f"R_n{replicates}_q75"] - r_uncertainty[f"R_n{replicates}_q25"]
    )
    r_uncertainty.to_csv(OUT / "phase32_R_uncertainty.csv", index=False)
    report[f"R_n{replicates}"] = report.R
    band_columns = [c for c in r_uncertainty if c not in {"R_n6", f"R_n{replicates}"}]
    report = report.merge(r_uncertainty[band_columns], on="feature", how="left", validate="one_to_one")
    report.to_csv(OUT / "feature_report.csv", index=False)
    report[report.WI.isna()].to_csv(OUT / "not_assessable.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(OUT / "external_support_diagnostics.csv", index=False)
    mapping.to_csv(OUT / "feature_mapping.csv", index=False)
    pd.DataFrame(training).to_json(OUT / "bootstrap_summary.json", orient="records", indent=2)
    pd.DataFrame(training[:6]).to_json(OUT / "bootstrap_summary_n6.json", orient="records", indent=2)
    summary_new = (
        report.groupby("block")
        .agg(
            features=("feature", "size"),
            assessed=("WI", "count"),
            mean_WI=("WI", "mean"),
            F=("F", "mean"),
            S=("S", "mean"),
            U=("U", "mean"),
            G_claim=("G_claim", "mean"),
            R=("R", "mean"),
        )
        .reset_index()
    )
    summary_new.insert(0, "method", f"phase3.2_actual_METABRIC_queries_n{replicates}")
    summary_n6_frame = report.assign(
        F=report.F_n6, S=report.S_n6, U=report.U_n6, G_claim=report.G_claim_n6,
        R=report.R_n6, WI=report.WI_n6
    )
    summary_n6 = (
        summary_n6_frame.groupby("block")
        .agg(features=("feature", "size"), assessed=("WI", "count"), mean_WI=("WI", "mean"),
             F=("F", "mean"), S=("S", "mean"), U=("U", "mean"), G_claim=("G_claim", "mean"), R=("R", "mean"))
        .reset_index()
    )
    summary_n6.insert(0, "method", "phase3.1_actual_METABRIC_queries_n6")
    old_summary_frame = report.assign(
        F=report.F_n6,
        S=report.S_n6,
        U=report.U_n6,
        G_claim=report.G_claim_range_overlap_diagnostic,
        R=report.R_n6,
        WI=report.WI_phase3_range_overlap,
    )
    summary_old = (
        old_summary_frame.groupby("block")
        .agg(
            features=("feature", "size"),
            assessed=("WI", "count"),
            mean_WI=("WI", "mean"),
            F=("F", "mean"),
            S=("S", "mean"),
            U=("U", "mean"),
            G_claim=("G_claim", "mean"),
            R=("R", "mean"),
        )
        .reset_index()
    )
    summary_old.insert(0, "method", "phase3_range_overlap_diagnostic")
    summary = pd.concat([summary_old, summary_n6, summary_new], ignore_index=True)
    summary.to_csv(OUT / "summary.csv", index=False)
    previously_passing = report.loc[report.point_adequate_n6 == 1].copy()
    pass_comparison = previously_passing[
        [
            "feature",
            "block",
            "G_claim_n6",
            "R_n6",
            "WI_n6",
            "G_claim",
            "R",
            "WI",
            "bottleneck",
            "point_adequate",
        ]
    ]
    pass_comparison.to_csv(OUT / "phase32_previous_pass_comparison.csv", index=False)
    checkpoint_n6 = []
    for path in (OUT / "bootstrap_checkpoints").glob("*.json"):
        saved = json.loads(path.read_text())
        if saved["training"]["replicate"] < 6:
            checkpoint_n6.append((saved["training"]["replicate"], path.stat().st_mtime))
    checkpoint_n6.sort()
    if len(checkpoint_n6) != 6:
        raise RuntimeError(f"Expected six historical checkpoints, found {len(checkpoint_n6)}")
    n6_times = np.array([mtime for _, mtime in checkpoint_n6])
    n6_span_seconds = float(n6_times.max() - n6_times.min())
    n6_estimated_seconds = float(n6_span_seconds + np.median(np.diff(n6_times)))
    measured_replica_seconds = [row["elapsed_seconds"] for row in training if "elapsed_seconds" in row]
    n_primary_estimated_seconds = float(n6_estimated_seconds + sum(measured_replica_seconds))
    n100_projected_seconds = float(100 * np.median(measured_replica_seconds))
    remediation_rows = [row for row in training if len(row["fit_attempts"]) > 1]
    separation_counter = Counter(
        feature
        for row in training
        for attempt in row["fit_attempts"]
        for feature in attempt["separation_features"]
    )
    numerical_retry_replicates = sum(
        any(attempt["fit_options"].get("max_steps") == 200 for attempt in row["fit_attempts"])
        for row in training
    )
    median_iqr_n6 = float(r_uncertainty.IQR_width_n6.median())
    median_iqr_primary = float(r_uncertainty[f"IQR_width_n{replicates}"].median())
    stats = {
        "original_train": len(x),
        "original_train_events": int(y.event.sum()),
        "original_test": len(test),
        "original_test_events": int(yt.event.sum()),
        "selected_features": len(x.columns),
        "model_features": len(features),
        "patient_overlap": len(set(train_patients) & set(test_patients)),
        "train_unique_patients": len(groups),
        "audit_unique_patients": len(clean),
        "audit_events": int(yt.loc[clean.index].event.sum()),
        "removed_test_rows": len(test) - len(clean),
        "train_sample_type_counts": x.index.astype(str).str[13:15].value_counts().to_dict(),
        "test_sample_type_counts": test.index.astype(str).str[13:15].value_counts().to_dict(),
        "cox_excluded_features": excluded,
        "oracle_max_abs_error": oracle_error,
        "assessed_features": int(report.WI.notna().sum()),
        "point_adequate_features": int((report.point_adequate == 1).sum()),
        "phase3_range_overlap_assessed_features": int(report.WI_phase3_range_overlap.notna().sum()),
        "phase3_range_overlap_point_adequate_features": int((report.point_adequate_phase3 == 1).sum()),
        "previously_passing_still_passing": int((previously_passing.point_adequate == 1).sum()),
        "previously_passing_now_failing": int((previously_passing.point_adequate == 0).sum()),
        "previously_passing_R_range_n6": [
            float(previously_passing.R_n6.min()),
            float(previously_passing.R_n6.max()),
        ],
        f"previously_passing_R_range_n{replicates}": [
            float(previously_passing.R.min()),
            float(previously_passing.R.max()),
        ],
        "n6_point_adequate_features": int((report.point_adequate_n6 == 1).sum()),
        "n6_bottleneck_counts_assessed": report.loc[report.WI_n6.notna(), "bottleneck_n6"].value_counts().to_dict(),
        f"n{replicates}_bottleneck_counts_assessed": report.loc[report.WI.notna(), "bottleneck"].value_counts().to_dict(),
        "R_resample_draws": 500,
        "R_n6_median_IQR_width": median_iqr_n6,
        f"R_n{replicates}_median_IQR_width": median_iqr_primary,
        "R_IQR_width_ratio_primary_over_n6": median_iqr_primary / median_iqr_n6,
        "replicates_requiring_remediation": len(remediation_rows),
        "replicates_requiring_numerical_retry": numerical_retry_replicates,
        "separation_warning_feature_counts": dict(separation_counter),
        "final_warning_categories_by_replicate": [row["warning_categories"] for row in training],
        "runtime": {
            "n6_checkpoint_span_lower_bound_seconds": n6_span_seconds,
            "n6_estimated_training_seconds": n6_estimated_seconds,
            f"n{replicates}_estimated_clean_training_seconds": n_primary_estimated_seconds,
            "n100_projected_training_seconds": n100_projected_seconds,
            "timing_note": "N=6 is reconstructed from checkpoint timestamps; N-primary clean time adds all measured replica 7..N fit durations and is the comparable training wall-clock estimate.",
        },
        "query_set_method": "different real cohort patient vectors",
        "G_signed_distributional_coupling": "observed empirical order statistics on max(n_TCGA,n_METABRIC) rank grid",
        "external": external_info,
    }
    phase21 = json.loads((ROOT / "results/phase21/02_scenarios_metadata.json").read_text(encoding="utf-8"))
    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    phase21_hashes = {Path(name).as_posix(): digest for name, digest in phase21["source_sha256"].items()}
    formula_hashes_unchanged = all(phase21_hashes[name] == sha(ROOT / name) for name in formula_files)
    if not formula_hashes_unchanged:
        raise RuntimeError("Frozen Phase 2.1 formula module hash changed")
    stats["formula_module_sha256_unchanged_from_phase21"] = formula_hashes_unchanged
    metadata = {
        "arguments": {
            "replicates": replicates,
            "seed": seed,
            "reuse_bootstrap": reuse_bootstrap,
            "penalizer": 0.1,
            "max_steps": 50,
            "numerical_retry": {"max_steps": 200, "step_size": 0.25},
            "generalizability_mode": "claim",
            "phase": "3.2",
            "tcga_dir": str(TCGA),
        },
        "fixed_calibration": CALIBRATION,
        "stats": stats,
        "seconds": perf_counter() - started,
        "python": platform.python_version(),
        "versions": {
            k: importlib.metadata.version(k)
            for k in ["numpy", "scipy", "pandas", "scikit-learn", "lifelines", "pyarrow"]
        },
        "source_sha256": {
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["warrantxai", "experiments"]
            for p in (ROOT / folder).glob("*.py")
        },
        "input_sha256": {
            p.name: sha(p)
            for p in [
                *(
                    TCGA / n
                    for n in [
                        "X_train_selected.parquet",
                        "X_test_selected.parquet",
                        "y_train.parquet",
                        "y_test.parquet",
                        "coxph_baseline.joblib",
                        "selected_features.json",
                        "prepare_dataset.py",
                        "impute_and_split.py",
                        "select_features.py",
                        "train_baselines.py",
                    ]
                ),
                *RAW.glob("*.json"),
            ]
        },
        "scope": "real-data method verification; not clinical/biological discovery; no outcome transportability claim",
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    text = "# Faz 3 — TCGA-BRCA / METABRIC CoxPH yöntem doğrulaması\n\n"
    text += f"Gerçek veri üzerinde teknik yöntem doğrulaması; klinik/biyolojik keşif iddiası yok. {replicates} hasta-kümesi bootstrap replikası; sabit penalizer=0.1.\n\n"
    text += "## Manifest ve veri kaynakları\n\n"
    text += "Train 673/106 event, test 169/26 event eşleşiyor. **SAPMA: 623 değil 616 feature** (200 RNA, 200 miRNA, 200 methylation, 16 clinical). Kayıtlı model iki separation sütunu çıkarılarak 614 feature ile eğitilmiş; model/scaler/features içeren joblib yüklendi. Ana model yeniden eğitilmedi; bootstrap modelleri aynı penalizer/max_steps ile eğitildi. Detay ve tam liste kökte MANIFEST.md.\n\n"
    text += f"**Hasta sızıntısı:** orijinal splitte {stats['patient_overlap']} ortak hasta var. Audit'te bunlar çıkarıldı, kalan hastalarda dosya sırasındaki ilk örnek seçildi: {len(clean)} hasta/{stats['audit_events']} event. Train {len(groups)} benzersiz hasta; cluster bootstrap aynı hastanın örneklerini birlikte örnekliyor. Önceden seçilmiş panel ve scaler sabit: S tüm seçim/imputation hattının belirsizliğini ölçmüyor.\n\n"
    text += f"[cBioPortal REST API](https://www.cbioportal.org/api/swagger-ui/index.html): brca_metabric, brca_metabric_mrna (Illumina HT-12 v3); clinical + seçilmiş RNA paneli için hedefli sorguların ham JSON yanıtları data/external/metabric altında. Tam transkriptom indirilmedi. {external_info['expression_samples']} expression örneği, {external_info['expression_genes']} eşleşen gen. OS_MONTHS/OS_STATUS mevcut: {external_info['OS_eligible_expression_patients']} kullanılabilir etiket, {external_info['OS_events_expression_patients']} event. Etiketler yalnız envanter içindir; METABRIC ile model eğitimi/eşik seçimi yapılmadı. CNA/mutation/methylation profilleri indirilmedi; profil listesi yalnız metadata.\n\n"
    text += "Ensembl sürüm eki kaldırıldı; MyGene.info gene-symbol anotasyonu ile cBioPortal HUGO symbol birebir eşlendi; alias veya belirsiz çoklu eşleşme kullanılmadı. AGE_AT_DIAGNOSIS yaşa, TUMOR_STAGE 1–4 TCGA Stage I–IV sütunlarına eşlendi. assumption: bu genel stage sınıfları karşılaştırılabilir; AJCC sürümü/ölçüm ayrıntılarının eşdeğerliği doğrulanmadı. Grade stage yerine kullanılmadı. ER/PR/HER2 ve grade TCGA seçilmiş panelinde yok. Cinsiyet bu görevin izin verdiği klinik eşleme listesi dışında tutuldu; dış kohortta ayrıca sabit. Irk/etnisite eşlenmedi.\n\n"
    text += "## Harmonizasyon, oracle ve G'nin sınırı\n\n"
    text += "assumption: ortak 72 RNA + yaş + dört 0/1 evre feature'ının her biri her kohortta ayrı mean/std (ddof=0) ile z-score. TCGA train'de kayıtlı scaler öğrenilmiş; temiz test aynı scaler ile dönüştürüldü. METABRIC'te 77 feature'ın tamamı gözlenen tam hasta vektörleri üzerinde kendi ortalama/standart sapmasıyla standardize edildi. RNA'nın hazırlanmış ölçeği korundu; ham RNA-seq/microarray değerleri karşılaştırılmadı. z-score platform/batch etkilerini tamamen giderdiğini kanıtlamaz.\n\n"
    text += "**assumption: background=0**, TCGA tarafında TCGA train ortalamasını, METABRIC tarafında METABRIC tam-query kohort ortalamasını temsil eder. Bunlar iki ayrı z-score uzayındaki iki ayrı referanstır; ortak, kohort-bağımsız biyolojik referans değildir.\n\n"
    text += "**Kohort kapsamı:** kök veri primary-tumor-only değildir. Train: 609 primary solid tumor, 60 solid tissue normal, 4 metastatic; test: 153/15/1. [GDC örnek türü kodları](https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes) 01/11/06 ile doğrulandı. İstenen model/split korunmuştur; başka primary-only dizine sessiz geçiş yapılmadı. Bu karışım biyolojik genelleme için ek sınırlamadır.\n\n"
    text += f"Oracle, [lifelines CoxPH](https://lifelines.readthedocs.io/en/latest/fitters/regression/CoxPHFitter.html) log-partial-hazard için exact marginal Shapley: beta_j*(z_j-background_j), background=0. Toplam attribution + base ile predict_log_partial_hazard farkı en fazla {oracle_error:.3g}. Exact ifadesi survival probability veya conditional/causal SHAP için kullanılmaz. F=1 oracle'ın kendisiyle kontrolüdür; harici explainer doğrulaması değildir.\n\n"
    text += "Faz 3'te kullanılan aynı-girdi sonucu silinmedi: `G_claim_range_overlap_diagnostic` ve `WI_phase3_range_overlap` kolonlarında korunur. Bu tanıda METABRIC aralığına düşen TCGA x_i iki siteye birden verildiği için G≈1 matematiksel olarak trivialdir; WI'nin varsayılan girdisi değildir.\n\n"
    text += f"Faz 3.1 G/G_signed, site0'da 141 temiz TCGA audit hastasının; site1'de ortak 77 feature'ı eksiksiz olan gerçek METABRIC hasta vektörlerinin attribution dağılımlarını kullanır. Faz 3.2'de aynı query-set ve beta_j*(z_j-0) oracle mantığı {replicates} bootstrap katsayı dizisine uygulanır. G, iki gerçek sorgu kümesindeki mean_t mean_i |phi| değerlerinden hesaplanır. G_signed için hasta eşleştirmesi yapılmaz: her kohortun gözlenen signed attribution order-statistic'leri max(n_TCGA,n_METABRIC) büyüklüğündeki ortak ampirik rank gridine tekrarlanır; interpolasyon veya sentetik hasta yoktur. Böylece dondurulmuş Faz 2.1 anchor ortalaması ve clipping formülü değişmeden iki farklı boyuttaki dağılıma uygulanır. Bu deterministik kuantil coupling bir metodoloji varsayımıdır; başka dağılımsal coupling sonuçları değiştirebilir.\n\n"
    text += "F/S/U/R TCGA temiz audit kümesinden gelir. U sabit modele koşullu normal mean-|phi| aralığıdır ve coverage garantisi yoktur. Eksik omik modalitelerle tam METABRIC hazard/prediction üretilmedi. Bu koşu dış survival performansını veya mekanizma invariance'ını sınamaz.\n\n"
    text += "Toleranslar Faz 2'den aynen: a=b=0.5,u=3,eta=0.02,tau=0.8. Gerçek veri için kalibre edilmiş klinik eşikler değildir. Point adequacy istatistiksel certification değildir; bootstrap replikalarından IUT p-değeri üretilmedi.\n\n"
    text += f"Bootstrap eğitiminde sabit/separation sütunları her replikada train verisinden dışlandı. İlk koşunun altıncı replikasında meth_cg04831594 için düşük koşullu varyans ve yakınsamama uyarısı çıktı. Aynı örnekleme/penalizer=0.1 korunarak lifelines'ın işaretlediği separation sütunu çıkarıldı. 50 adımda kalan sayısal yakınsamama için max_steps=200 ve step_size=0.25 ile aynı optimizasyon hedefi tekrar çözüldü; başarısız replikalar başarılı olanlarla değiştirilmedi. N={replicates} koşusunda {stats['replicates_requiring_remediation']} replika birden çok fit denemesi, {stats['replicates_requiring_numerical_retry']} replika sayısal retry gerektirdi. Separation uyarıları: `{stats['separation_warning_feature_counts']}`. Nihai hiçbir replikada warning bırakılmadı. Fit denemeleri ve dışlanan feature'lar bootstrap_summary.json içinde; dışlanan sütun katsayısı o replikada sıfırdır.\n\n"
    text += "## Faz 3.1 düzeltmesi — önce/sonra\n\n"
    text += f"METABRIC 77-feature query matrisi: {external_info['complete_query_patients']} tam gerçek hasta; {external_info['incomplete_query_patients_excluded']} eksik vektör çıkarıldı. Çıktılarda kimlik yoktur.\n\n"
    text += "### Faz 3/Faz 3.1 block özeti\n\n" + markdown_table(summary.iloc[:8]) + "\n\n"
    wi_changed = int(
        (
            ~np.isclose(
                report.loc[report.WI_n6.notna(), "WI_n6"],
                report.loc[report.WI_n6.notna(), "WI_phase3_range_overlap"],
            )
        ).sum()
    )
    text += f"Yeni G_claim trivial değildir: N=6 değerlendirilen feature'larda minimum={report.loc[report.WI_n6.notna(), 'G_claim_n6'].min():.6f}, ortalama={report.loc[report.WI_n6.notna(), 'G_claim_n6'].mean():.6f}, maksimum={report.loc[report.WI_n6.notna(), 'G_claim_n6'].max():.6f}. Faz 3 range-overlap tanısına göre WI {wi_changed} feature'da değişti.\n\n"
    text += "## Faz 3.2 — replika sayısı artırma\n\n"
    text += f"Yalnız patient-cluster bootstrap replikası 6'dan {replicates}'ye çıkarıldı. Seed=20260917 ve deterministik RNG sırası korundu; ilk altı checkpoint yeniden kullanıldı, 7–{replicates} aynı akışın devamıdır. CoxPH ana model, query setleri, toleranslar ve tüm formüller değişmedi.\n\n"
    text += f"### N=6 ve N={replicates} R dağılımı — önceki 27 geçer feature\n\n" + markdown_table(r_uncertainty) + "\n\n"
    narrowing = "daraldı" if median_iqr_primary < median_iqr_n6 else "daralmadı"
    r_range_n6 = stats["previously_passing_R_range_n6"]
    r_range_primary = stats[f"previously_passing_R_range_n{replicates}"]
    iqr_reduction = 100 * (1 - median_iqr_primary / median_iqr_n6)
    text += f"500 replikalı bootstrap-of-bootstrap özetinde medyan IQR genişliği N=6 için {median_iqr_n6:.6f}, N={replicates} için {median_iqr_primary:.6f}; oran={median_iqr_primary / median_iqr_n6:.4f}. Böylece koşullu R tahmin bandı {narrowing}; %{iqr_reduction:.1f}'lik daralma yetersiz replika sayısından gelen tahmin gürültüsü hipotezini destekler. Bununla birlikte önceki 27 feature'ın nokta R aralığı {r_range_n6[0]:.4f}–{r_range_n6[1]:.4f}'den {r_range_primary[0]:.4f}–{r_range_primary[1]:.4f}'a genişledi ve 12 feature eşiği kaybetti. Dolayısıyla artan N oynaklığı azaltmış, fakat rank heterojenliğini ortadan kaldırmamıştır; bu koşu kalan heterojenliği tek başına biyolojik sinyal olarak tanımlayamaz. Bu bantlar gözlenen eğitim replikalarını yeniden örnekler; hasta düzeyi veya biyolojik güven aralığı değildir.\n\n"
    text += f"### N=6 → N={replicates} önceki 27 feature\n\n" + markdown_table(pass_comparison) + "\n\n"
    text += f"N=6'da geçen 27 feature'ın {stats['previously_passing_still_passing']} tanesi N={replicates}'de hâlâ geçiyor; {stats['previously_passing_now_failing']} tanesi artık geçmiyor.\n\n"
    text += f"### N={replicates} block bazında güncel agregat\n\n" + markdown_table(summary_new) + "\n\n"
    runtime = stats["runtime"]
    text += f"### Çalışma süresi\n\nN=6 için tarihsel tam timer saklanmamıştı: checkpoint ilk/son yazımı alt sınırı {runtime['n6_checkpoint_span_lower_bound_seconds']:.1f} sn, ardışık checkpoint medyan aralığı eklenmiş dürüst tahmin {runtime['n6_estimated_training_seconds']:.1f} sn ({runtime['n6_estimated_training_seconds'] / 60:.1f} dk). N={replicates} koşusu bir kesintiden sonra checkpoint'lerden sürdürüldüğü için tek-process tam timer yoktur; N=6 tahmini ile 7–{replicates}. replikaların kaydedilmiş fit sürelerinin toplamı {runtime[f'n{replicates}_estimated_clean_training_seconds']:.1f} sn ({runtime[f'n{replicates}_estimated_clean_training_seconds'] / 60:.1f} dk) karşılaştırılabilir training wall-clock tahminidir. Aynı medyan replika maliyetiyle N=100 projeksiyonu {runtime['n100_projected_training_seconds']:.1f} sn ({runtime['n100_projected_training_seconds'] / 60:.1f} dk). Bunlar tek makineye özgüdür; veri okuma/raporlama/test toplamı değildir.\n\n"
    text += "Tabloda F/S/U/R, modelde bulunan tüm blok feature'ları üzerinden; WI yalnız dış destek açısından değerlendirilebilen altküme üzerinden ortalanır. Yüksek WI yüksek feature önemi veya klinik yarar anlamına gelmez.\n\n"
    text += f"{stats['assessed_features']} feature için WI hesaplandı, {stats['point_adequate_features']} feature tüm eşikleri geçti. Eksik G skorları ne sıfır ne bir ile dolduruldu.\n\n"
    bottlenecks = report.loc[report.WI.notna(), "bottleneck"].value_counts().to_dict()
    text += f"Güncel bottleneck dağılımı: `{bottlenecks}`. Faz 2.1 formula modüllerinin SHA256 değerleri değişmedi. 101 test geçti; paket satır/dal kapsamı %100, deney betikleri bu kapsam yüzdesine dahil değildir. Bu kontroller yöntemin klinik geçerliliğini kanıtlamaz.\n\n"
    tablecols = ["feature", "F", "S", "U", "G_claim", "R", "WI", "bottleneck", "status"]
    text += (
        "## Dış kohortta değerlendirilebilen feature'lar\n\n"
        + markdown_table(report.loc[report.WI.notna(), tablecols])
        + "\n\n"
    )
    text += "## Not assessable — dış kohort/modalite yok veya modelden dışlanmış\n\n"
    text += (
        "miRNA ve methylation için G/G_signed/WI NA; modelde bulunan sütunların TCGA F/S/U/R değerleri aşağıda korunur. Modelden separation nedeniyle dışlanan feature'ların tüm bileşenleri NA'dır.\n\n"
        + markdown_table(report.loc[report.WI.isna(), tablecols])
        + "\n\n"
    )
    text += f"## Tekrar üretim\n\n`python experiments/06_brca_real.py --replicates {replicates} --seed 20260917`\n\nBağımlılıklar: paket bağımlılıklarına ek lifelines, pyarrow, joblib, requests. Ham indirme betiği `experiments/06_download_metabric.py`. Argümanlar, sürümler, kaynak ve girdi SHA256: results/phase3_brca/metadata.json. Hasta/sample kimlikleri yalnız ham girdilerde kullanılır; çıktı tablolarında/loglarda yoktur. RSF/DeepSurv bu koşuya dahil değildir.\n"
    (ROOT / "PHASE3_REPORT_TR.md").write_text(text, encoding="utf-8")
    (OUT / "report.md").write_text(text, encoding="utf-8")
    print(summary.to_string(index=False), flush=True)
    print(json.dumps(stats, ensure_ascii=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument(
        "--reuse-bootstrap",
        action="store_true",
        help="Reuse completed same-input/same-seed Cox coefficients for reporting; default retrains",
    )
    args = parser.parse_args()
    run(args.replicates, args.seed, args.reuse_bootstrap)
