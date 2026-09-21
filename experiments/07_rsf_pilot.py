"""Phase 4 pilot: model-agnostic approximate-F for the existing TCGA RSF.

This script deliberately stops at two patient-cluster bootstrap replicas and
12 identity-free query rows per cohort. It does not launch the N=50 audit.
"""

from __future__ import annotations

import argparse
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
import shap
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
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
from _shared import CALIBRATION, markdown_table

ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase4_rsf_pilot"
OUT_FULL = ROOT / "results/phase4_rsf"


def sha(path: Path) -> str:
    """Return a whole-file SHA256."""
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def frame(name: str) -> pd.DataFrame:
    """Read a prepared TCGA frame without emitting sample keys."""
    return pd.read_parquet(TCGA / name).set_index("sample_key")


def phase3_helpers():
    """Load the frozen Phase 3 external mapping implementation."""
    path = ROOT / "experiments/06_brca_real.py"
    spec = importlib.util.spec_from_file_location("phase3_brca_helpers", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Phase 3 helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def permutation_attributions(
    model,
    query: np.ndarray,
    baseline: np.ndarray,
    active_indices: np.ndarray,
    *,
    permutations: int,
    seed: int,
) -> tuple[np.ndarray, float]:
    """Estimate Shapley values in a fixed-background 77-feature subgame.

    Each permutation evaluates the complete forward coalition path. Features
    outside ``active_indices`` stay at the TCGA train mean. The estimate is
    unbiased for this explicitly defined deterministic-background game.
    """
    x = np.asarray(query, dtype=np.float32)
    base = np.asarray(baseline, dtype=np.float32)
    active = np.asarray(active_indices, dtype=int)
    if x.ndim != 2 or x.shape[1] != len(active) or base.ndim != 1:
        raise ValueError("query/active/baseline shape mismatch")
    if permutations < 1 or len(np.unique(active)) != len(active):
        raise ValueError("positive budget and unique active indices required")
    rng = np.random.default_rng(seed)
    n, p = x.shape
    result = np.zeros((n, p), dtype=float)
    for _ in range(permutations):
        order = rng.permutation(p)
        states = np.repeat(base[None, :], n * (p + 1), axis=0)
        for row in range(n):
            current = base.copy()
            offset = row * (p + 1)
            for step, local_index in enumerate(order, start=1):
                current[active[local_index]] = x[row, local_index]
                states[offset + step] = current
        predictions = np.asarray(model.predict(states), dtype=float).reshape(n, p + 1)
        result[:, order] += np.diff(predictions, axis=1)
    result /= permutations
    full = np.repeat(base[None, :], n, axis=0)
    full[:, active] = x
    target_delta = np.asarray(model.predict(full), dtype=float) - float(model.predict(base[None, :])[0])
    return result, float(np.max(np.abs(result.sum(axis=1) - target_delta)))


def target_permutation_attributions(
    model,
    query: np.ndarray,
    background_rows: np.ndarray,
    active_mean: np.ndarray,
    active_indices: np.ndarray,
    target_local_indices: np.ndarray,
    *,
    permutations: int,
    seed: int,
) -> np.ndarray:
    """Estimate selected 77-game features while averaging inactive backgrounds."""
    x = np.asarray(query, dtype=np.float32)
    backgrounds = np.asarray(background_rows, dtype=np.float32).copy()
    active = np.asarray(active_indices, dtype=int)
    targets = np.asarray(target_local_indices, dtype=int)
    backgrounds[:, active] = np.asarray(active_mean, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != len(active) or backgrounds.ndim != 2:
        raise ValueError("query/background shape mismatch")
    rng = np.random.default_rng(seed)
    n, p = x.shape
    k = len(backgrounds)
    result = np.zeros((n, len(targets)), dtype=float)
    for _ in range(permutations):
        order = rng.permutation(p)
        positions = np.empty(p, dtype=int)
        positions[order] = np.arange(p)
        states = np.empty((n * len(targets) * 2 * k, backgrounds.shape[1]), dtype=np.float32)
        cursor = 0
        for row in range(n):
            for target in targets:
                before = backgrounds.copy()
                preceding = order[: positions[target]]
                before[:, active[preceding]] = x[row, preceding]
                after = before.copy()
                after[:, active[target]] = x[row, target]
                states[cursor : cursor + k] = before
                states[cursor + k : cursor + 2 * k] = after
                cursor += 2 * k
        predictions = np.asarray(model.predict(states), dtype=float).reshape(n, len(targets), 2, k)
        result += predictions[:, :, 1].mean(axis=-1) - predictions[:, :, 0].mean(axis=-1)
    return result / permutations


def background_permutation_attributions(
    model,
    query: np.ndarray,
    background_rows: np.ndarray,
    active_mean: np.ndarray,
    active_indices: np.ndarray,
    *,
    permutations: int,
    seed: int,
) -> np.ndarray:
    """Estimate all 77 Shapley values with an empirical inactive background."""
    x = np.asarray(query, dtype=np.float32)
    backgrounds = np.asarray(background_rows, dtype=np.float32).copy()
    active = np.asarray(active_indices, dtype=int)
    backgrounds[:, active] = np.asarray(active_mean, dtype=np.float32)
    if x.ndim != 2 or x.shape[1] != len(active) or backgrounds.ndim != 2:
        raise ValueError("query/background shape mismatch")
    rng = np.random.default_rng(seed)
    n, p = x.shape
    k = len(backgrounds)
    result = np.zeros((n, p), dtype=float)
    for _ in range(permutations):
        order = rng.permutation(p)
        states = np.empty((n * (p + 1) * k, backgrounds.shape[1]), dtype=np.float32)
        for row in range(n):
            block = np.repeat(backgrounds[None, :, :], p + 1, axis=0)
            for step, local_index in enumerate(order, start=1):
                block[step:, :, active[local_index]] = x[row, local_index]
            offset = row * (p + 1) * k
            states[offset : offset + (p + 1) * k] = block.reshape((p + 1) * k, -1)
        predictions = np.asarray(model.predict(states), dtype=float).reshape(n, p + 1, k).mean(axis=-1)
        result[:, order] += np.diff(predictions, axis=1)
    return result / permutations


def selected_f(candidate: list[np.ndarray], reference: list[np.ndarray]) -> np.ndarray:
    """Apply unchanged F to selected-feature site tensors and take worst site."""
    by_site = [
        faithfulness(phi[None, ...], psi[None, ...], CALIBRATION["a"])
        for phi, psi in zip(candidate, reference)
    ]
    return np.minimum(by_site[0], by_site[1])


def selected_g_claim(candidate: list[np.ndarray]) -> np.ndarray:
    """Apply unchanged G_claim to selected-feature candidate tensors."""
    site_importance = np.stack([np.mean(np.abs(phi), axis=(0, 1)) for phi in candidate])
    magnitude = generalizability(site_importance, CALIBRATION["b"])
    signed = np.stack([np.sort(phi.mean(axis=0), axis=0) for phi in candidate])
    return np.minimum(magnitude, generalizability_signed(signed, CALIBRATION["a"]))


def site_components(candidate: np.ndarray, reference: np.ndarray) -> dict[str, np.ndarray]:
    """Apply unchanged F/S/U/R formulas to one site's pilot tensors."""
    phi = candidate[None, ...]
    ref = reference[None, ...]
    lower, upper = importance_ci(phi)
    return {
        "F": faithfulness(phi, ref, CALIBRATION["a"]),
        "S": stability(phi, CALIBRATION["a"]),
        "U": uncertainty(lower, upper, CALIBRATION["u"]),
        "R": rank_reproducibility(importance(phi), CALIBRATION["eta"]),
    }


def component_report(
    candidates: np.ndarray,
    references: np.ndarray,
    feature_names: list[str],
) -> pd.DataFrame:
    """Compute the unchanged five components for a two-site attribution tensor."""
    if candidates.shape != references.shape or candidates.shape[0] != 2:
        raise ValueError("matching two-site candidate/reference tensors required")
    by_site = [site_components(candidates[site], references[site]) for site in range(2)]
    f = np.minimum(by_site[0]["F"], by_site[1]["F"])
    s = np.minimum(by_site[0]["S"], by_site[1]["S"])
    u = np.minimum(by_site[0]["U"], by_site[1]["U"])
    r = np.minimum(by_site[0]["R"], by_site[1]["R"])
    site_importance = np.stack(
        [np.mean(np.abs(candidates[site]), axis=(0, 1)) for site in range(2)]
    )
    g_magnitude = generalizability(site_importance, CALIBRATION["b"])
    signed = np.stack([np.sort(candidates[site].mean(axis=0), axis=0) for site in range(2)])
    g_signed = generalizability_signed(signed, CALIBRATION["a"])
    g_claim = np.minimum(g_magnitude, g_signed)
    scores = np.column_stack([f, s, u, g_claim, r])
    wi = warrant_index(scores)
    point = adequate(scores, CALIBRATION["tau"])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
        raise RuntimeError("Components must be finite and in [0,1]")
    np.testing.assert_allclose(wi, scores.min(axis=1))
    report = pd.DataFrame(scores, columns=["F", "S", "U", "G_claim", "R"])
    report.insert(0, "feature", feature_names)
    report["G_magnitude"] = g_magnitude
    report["G_signed"] = g_signed
    report["WI"] = wi
    report["point_adequate"] = point.astype(int)
    report["bottleneck"] = report[["F", "S", "U", "G_claim", "R"]].idxmin(axis=1)
    report["block"] = [feature.split("_")[0] for feature in feature_names]
    return report


def train_replicas(x: pd.DataFrame, y: pd.DataFrame, replicas: int, seed: int):
    """Fit deterministic patient-cluster RSF replicas with frozen settings."""
    patient_keys = x.index.astype(str).str[:12]
    groups = [np.flatnonzero(patient_keys == patient) for patient in patient_keys.unique()]
    rng = np.random.default_rng(seed)
    survival = Surv.from_arrays(event=y.event.astype(bool).to_numpy(), time=y.time.to_numpy())
    models = []
    records = []
    started = perf_counter()
    values = x.to_numpy(dtype=np.float32)
    for replicate in range(replicas):
        indices = np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))])
        model = RandomSurvivalForest(
            n_estimators=500,
            min_samples_split=10,
            min_samples_leaf=5,
            max_features="sqrt",
            random_state=seed + replicate,
            n_jobs=-1,
        )
        fit_started = perf_counter()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model.fit(values[indices], survival[indices])
        models.append(model)
        records.append(
            {
                "replicate": replicate,
                "rows": int(len(indices)),
                "events": int(survival[indices]["event"].sum()),
                "elapsed_seconds": perf_counter() - fit_started,
                "warning_categories": dict(Counter(w.category.__name__ for w in caught)),
            }
        )
    return models, records, perf_counter() - started


def run(*, replicas: int = 2, anchors_per_site: int = 12, candidate_budget: int = 8, reference_budget: int = 128,
        seed: int = 20260917) -> None:
    """Run the intentionally small RSF approximate-F design pilot."""
    if replicas != 2:
        raise ValueError("Phase 4 design pilot is intentionally fixed at N=2")
    if anchors_per_site < 2 or reference_budget < 10 * candidate_budget:
        raise ValueError("Need >=2 anchors and an independent reference budget at least 10x candidate")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    x = frame("X_train_selected.parquet")
    test = frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[x.index]
    baseline_path = TCGA / "rsf_baseline.joblib"
    bundle = joblib.load(baseline_path)
    baseline_model = bundle["model"]
    features = list(bundle["features"])
    if features != list(x.columns) or baseline_model.n_features_in_ != 616:
        raise RuntimeError("Saved RSF is not compatible with the frozen 616-feature panel")
    if not (
        baseline_model.n_estimators == 500
        and baseline_model.min_samples_split == 10
        and baseline_model.min_samples_leaf == 5
        and baseline_model.max_features == "sqrt"
        and baseline_model.random_state == 42
    ):
        raise RuntimeError("Saved RSF hyperparameters differ from the recorded baseline")

    tree_supported = bool(shap.explainers.Tree.supports_model_with_masker(baseline_model, None))
    tree_error = ""
    try:
        shap.TreeExplainer(baseline_model)
    except Exception as exc:  # exact exception is an audited compatibility result
        tree_error = f"{type(exc).__name__}: {exc}"
    if tree_supported or "not yet supported" not in tree_error:
        raise RuntimeError("TreeSHAP compatibility result changed; redesign before continuing")

    phase3 = phase3_helpers()
    external, mapping, external_info = phase3.external_columns(features)
    active_features = [feature for feature in features if feature in external]
    if len(active_features) != 77:
        raise RuntimeError(f"Expected 77 common RSF features, found {len(active_features)}")
    active_indices = np.array([features.index(feature) for feature in active_features])
    external_raw = np.column_stack([external[feature] for feature in active_features])
    external_complete = np.isfinite(external_raw).all(axis=1)
    external_raw = external_raw[external_complete]
    train_mean = x.mean(axis=0).to_numpy(float)
    train_std = x.std(axis=0, ddof=0).to_numpy(float)
    active_mean = train_mean[active_indices]
    active_std = train_std[active_indices]
    if np.any(active_std <= 0):
        raise RuntimeError("The 77-feature subgame contains a constant TCGA train feature")
    external_z = (external_raw - external_raw.mean(axis=0)) / external_raw.std(axis=0, ddof=0)
    external_aligned = active_mean + external_z * active_std

    train_patients = x.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    clean = test.loc[~test_patients.isin(train_patients)]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    source_pool = clean[active_features].to_numpy(float)
    query_rng = np.random.default_rng(seed + 4000)
    source_selection = np.sort(query_rng.choice(len(source_pool), anchors_per_site, replace=False))
    external_selection = np.sort(query_rng.choice(len(external_aligned), anchors_per_site, replace=False))
    queries = [source_pool[source_selection], external_aligned[external_selection]]

    models, training, training_seconds = train_replicas(x, y, replicas, seed)
    candidates = []
    references = []
    candidate_seconds = 0.0
    reference_seconds = 0.0
    additivity_errors = []
    for site, query in enumerate(queries):
        site_candidate = []
        site_reference = []
        for replicate, model in enumerate(models):
            timer = perf_counter()
            candidate, candidate_error = permutation_attributions(
                model,
                query,
                train_mean,
                active_indices,
                permutations=candidate_budget,
                seed=seed + 10000 + 100 * site + replicate,
            )
            candidate_seconds += perf_counter() - timer
            timer = perf_counter()
            reference, reference_error = permutation_attributions(
                model,
                query,
                train_mean,
                active_indices,
                permutations=reference_budget,
                seed=seed + 20000 + 100 * site + replicate,
            )
            reference_seconds += perf_counter() - timer
            site_candidate.append(candidate)
            site_reference.append(reference)
            additivity_errors.append(
                {
                    "site": site,
                    "replicate": replicate,
                    "candidate": candidate_error,
                    "reference": reference_error,
                }
            )
        candidates.append(np.asarray(site_candidate))
        references.append(np.asarray(site_reference))

    by_site = [site_components(candidate, reference) for candidate, reference in zip(candidates, references)]
    f = np.minimum(by_site[0]["F"], by_site[1]["F"])
    s = np.minimum(by_site[0]["S"], by_site[1]["S"])
    u = np.minimum(by_site[0]["U"], by_site[1]["U"])
    r = np.minimum(by_site[0]["R"], by_site[1]["R"])
    site_importance = np.stack([np.mean(np.abs(candidate), axis=(0, 1)) for candidate in candidates])
    g_magnitude = generalizability(site_importance, CALIBRATION["b"])
    signed = np.stack([np.sort(candidate.mean(axis=0), axis=0) for candidate in candidates])
    g_signed = generalizability_signed(signed, CALIBRATION["a"])
    g_claim = np.minimum(g_magnitude, g_signed)
    scores = np.column_stack([f, s, u, g_claim, r])
    wi = warrant_index(scores)
    point = adequate(scores, CALIBRATION["tau"])
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
        raise RuntimeError("Pilot components must be finite and in [0,1]")
    np.testing.assert_allclose(wi, scores.min(axis=1))

    report = pd.DataFrame(scores, columns=["F", "S", "U", "G_claim", "R"])
    report.insert(0, "feature", active_features)
    report["G_magnitude"] = g_magnitude
    report["G_signed"] = g_signed
    report["WI"] = wi
    report["point_adequate"] = point.astype(int)
    report["bottleneck"] = report[["F", "S", "U", "G_claim", "R"]].idxmin(axis=1)
    report["block"] = [feature.split("_")[0] for feature in active_features]
    report.to_csv(OUT / "feature_report.csv", index=False)
    pd.DataFrame(training).to_json(OUT / "bootstrap_summary.json", orient="records", indent=2)
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)

    total_seconds = perf_counter() - started
    fixed_seconds = max(0.0, total_seconds - training_seconds - candidate_seconds - reference_seconds)
    projected = {}
    for n in [6, 50]:
        factor = n / replicas
        projected[f"n{n}_same_12x12_query_seconds"] = fixed_seconds + factor * (
            training_seconds + candidate_seconds + reference_seconds
        )
        query_factor = (len(source_pool) + len(external_aligned)) / (2 * anchors_per_site)
        projected[f"n{n}_full_G_query_upper_seconds"] = fixed_seconds + factor * (
            training_seconds + reference_seconds + candidate_seconds * query_factor
        )
    metadata = {
        "phase": "4_pilot",
        "arguments": {
            "replicas": replicas,
            "anchors_per_site": anchors_per_site,
            "candidate_permutations": candidate_budget,
            "reference_permutations": reference_budget,
            "seed": seed,
        },
        "model": {
            "source": str(baseline_path),
            "reused_without_retraining": True,
            "sha256": sha(baseline_path),
            "feature_count": len(features),
            "hyperparameters": {
                "n_estimators": 500,
                "min_samples_split": 10,
                "min_samples_leaf": 5,
                "max_features": "sqrt",
                "random_state": 42,
            },
        },
        "explainer": {
            "selected_path": "model-agnostic permutation Shapley",
            "game": "77 active features; other 539 coordinates fixed at TCGA train mean",
            "candidate_reference_budget_ratio": reference_budget / candidate_budget,
            "independent_candidate_reference_seeds": True,
            "tree_shap_supported": tree_supported,
            "tree_shap_error": tree_error,
            "max_additivity_error": max(max(row["candidate"], row["reference"]) for row in additivity_errors),
        },
        "queries": {
            "source_pool": len(source_pool),
            "external_complete_pool": len(external_aligned),
            "sampled_per_site": anchors_per_site,
            "external_source": "actual METABRIC patient vectors, identity-free deterministic subset",
            "harmonization": "cohort z-score mapped to TCGA train mean/std for the existing raw-scale RSF",
        },
        "training": training,
        "timing_seconds": {
            "training_n2": training_seconds,
            "candidate_explanations_n2": candidate_seconds,
            "reference_explanations_n2": reference_seconds,
            "fixed_and_reporting": fixed_seconds,
            "total_n2": total_seconds,
            **projected,
        },
        "results": {
            "assessed_features": len(report),
            "point_adequate_features": int(point.sum()),
            "component_means": report[["F", "S", "U", "G_claim", "R", "WI"]].mean().to_dict(),
            "bottleneck_counts": report.bottleneck.value_counts().to_dict(),
        },
        "external_inventory": external_info,
        "fixed_calibration": CALIBRATION,
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "scikit-learn", "scikit-survival", "shap"]
        },
        "python": platform.python_version(),
        "source_sha256": {"experiments/07_rsf_pilot.py": sha(Path(__file__))},
        "input_sha256": {
            name: sha(TCGA / name)
            for name in [
                "rsf_baseline.joblib",
                "X_train_selected.parquet",
                "X_test_selected.parquet",
                "y_train.parquet",
                "y_test.parquet",
            ]
        },
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    examples = pd.concat([report.nsmallest(4, "WI"), report.nlargest(4, "WI")]).drop_duplicates("feature")
    aggregate = report.groupby("block").agg(
        features=("feature", "size"),
        mean_F=("F", "mean"),
        mean_S=("S", "mean"),
        mean_U=("U", "mean"),
        mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"),
        mean_WI=("WI", "mean"),
    ).reset_index()
    note = "# Faz 4 — RSF approximate-F tasarım notu\n\n"
    note += "## Karar ve model kaynağı\n\n"
    note += (
        "Kardeş projedeki mevcut `rsf_baseline.joblib` yeniden kullanıldı; ana model yeniden eğitilmedi. "
        "Model aynı 616-feature train/test paneline aittir ve 500 ağaç, min_samples_split=10, "
        "min_samples_leaf=5, max_features=sqrt, random_state=42 kullanır. Pilotun N=2 eğitim "
        "replikaları aynı ayarlarla hasta-kümesi bootstrap örneklerinde ayrıca fit edildi.\n\n"
    )
    note += "SHAP 0.52.0 `TreeExplainer`, scikit-survival 0.28.0 `RandomSurvivalForest` modelini desteklemedi: "
    note += f"`{tree_error}`. Bu nedenle TreeSHAP seçilmedi.\n\n"
    note += (
        "Kaynaklar: [SHAP TreeExplainer belgesi](https://shap.readthedocs.io/en/latest/generated/"
        "shap.TreeExplainer.html), [açık scikit-survival RSF destek talebi](https://github.com/shap/shap/"
        "issues/3308), [scikit-survival RSF API](https://scikit-survival.readthedocs.io/en/stable/api/"
        "generated/sksurv.ensemble.RandomSurvivalForest.html).\n\n"
    )
    note += "## Approximate-F tanımı\n\n"
    note += (
        f"Candidate, her model/query için {candidate_budget} rastgele permütasyonlu model-agnostik Shapley; "
        f"bağımsız referans {reference_budget} permütasyonlu aynı kestiricidir. Referans bütçesi candidate'ın "
        f"{reference_budget / candidate_budget:.0f} katıdır ve seed akışları ayrıdır. F formülü değiştirilmedi. "
        "Açıklanan skaler çıktı scikit-survival `RandomSurvivalForest.predict()` risk skorudur. "
        "Bu, kapalı-form oracle veya algoritmik exact TreeSHAP değildir; F hem candidate hem sonlu-bütçeli "
        "referans Monte Carlo hatasını içerir.\n\n"
    )
    note += (
        "Mevcut RSF 616 girdi isterken METABRIC yalnız 77 ortak feature sağlar. Eksik 539 modaliteyi hasta "
        "verisiymiş gibi doldurmak yerine açıklama oyunu açıkça 77-boyutlu bir model kesiti olarak tanımlandı: "
        "ortak feature'lar coalition ile değişir, diğer koordinatlar TCGA train ortalamasında sabit kalır. "
        "Bu pilot sonucu tam 616-feature RSF Shapley açıklaması değildir. METABRIC ortak feature'ları kendi "
        "kohort z-score'larından TCGA train mean/std ölçeğine taşındı.\n\n"
    )
    note += "## Pilot kapsamı ve sonuçlar\n\n"
    note += (
        f"N=2, 72 RNA + 5 klinik feature ve kohort başına {anchors_per_site} gerçek, deterministik seçilmiş "
        "kimliksiz query vektörü kullanıldı. Kaynak havuzu "
        f"{len(source_pool)}, tam METABRIC havuzu {len(external_aligned)} hastadır. Küçük query örneği "
        "özellikle G/G_signed için yüksek örnekleme oynaklığı yaratabilir. Tüm F/S/U/G_claim/R/WI değerleri "
        "sonlu, [0,1] aralığında ve WI=min(F,S,U,G_claim,R) eşitliğini sağladı. "
        f"{int(point.sum())}/77 feature WI>0.8 koşulunu geçti; bottleneck dağılımı "
        f"`{report.bottleneck.value_counts().to_dict()}`. F aralığı {report.F.min():.4f}–{report.F.max():.4f}; "
        "bu yüksek değer aynı kestiricinin 8 ve 128 permütasyonlu sürümlerini karşılaştırmanın beklenen "
        "sonucudur ve bağımsız açıklayıcı doğruluğunu kanıtlamaz.\n\n"
    )
    note += "### Block özeti\n\n" + markdown_table(aggregate) + "\n\n"
    note += "### Örnek feature sonuçları\n\n" + markdown_table(
        examples[["feature", "F", "S", "U", "G_claim", "R", "WI", "bottleneck"]]
    ) + "\n\n"
    timing = metadata["timing_seconds"]
    note += "## Süre ve ölçekleme\n\n"
    note += (
        f"N=2 pilot toplamı {timing['total_n2']:.1f} sn; eğitim {training_seconds:.1f} sn, candidate "
        f"açıklamalar {candidate_seconds:.1f} sn, yüksek-bütçeli referans {reference_seconds:.1f} sn. Aynı "
        f"12+12 query ve bütçelerle doğrusal N=6 projeksiyonu {projected['n6_same_12x12_query_seconds']:.1f} sn, "
        f"N=50 projeksiyonu {projected['n50_same_12x12_query_seconds']:.1f} sn. Candidate açıklamayı tüm "
        f"{len(source_pool)}+{len(external_aligned)} G query havuzuna genişleten kaba üst projeksiyon N=6 için "
        f"{projected['n6_full_G_query_upper_seconds']:.1f} sn, N=50 için "
        f"{projected['n50_full_G_query_upper_seconds']:.1f} sn. Bunlar doğrusal tek-makine projeksiyonlarıdır.\n\n"
    )
    note += "## Sürprizler, sınırlar ve tam koşu öncesi karar noktası\n\n"
    note += (
        "TreeSHAP yolu mevcut kütüphane/model çiftiyle uygulanamadı. Daha temel sınır, dış kohortun 616 RSF "
        "girdisinin yalnız 77'sini sağlamasıdır; seçilen 77-feature kesiti bilimsel olarak ayrıca onaylanmalıdır. "
        "N=2 ile R tek bir replika çiftine dayanır; 12+12 query ile U ve G hassastır. Toleranslar "
        "a=b=0.5,u=3,eta=0.02,tau=0.8 kalibre edilmemiştir. Pilot teknik uygulanabilirliği gösterir; "
        "klinik geçerlilik, exact F veya tam-model açıklaması göstermez. N=50 tam koşu başlatılmadı.\n"
    )
    (ROOT / "PHASE4_DESIGN_NOTE.md").write_text(note, encoding="utf-8")
    (OUT / "design_note.md").write_text(note, encoding="utf-8")
    print(aggregate.to_string(index=False), flush=True)
    print(json.dumps(metadata["timing_seconds"], indent=2), flush=True)


def run_sensitivity(*, seed: int = 20260917, background_patients: int = 16) -> None:
    """Run only the two bounded Phase 4.1 robustness checks."""
    started = perf_counter()
    base_report_path = OUT / "feature_report.csv"
    base_note_path = ROOT / "PHASE4_DESIGN_NOTE.md"
    if not base_report_path.exists() or not base_note_path.exists():
        raise FileNotFoundError("Run the N=2 Phase 4 pilot before --sensitivity")
    base_report = pd.read_csv(base_report_path)
    if len(base_report) != 77:
        raise RuntimeError("Phase 4 pilot report must contain 77 features")
    low = base_report.nsmallest(3, "WI").feature.tolist()
    high = base_report.nlargest(3, "WI").feature.tolist()
    remaining = base_report.loc[~base_report.feature.isin(low + high), "feature"].to_numpy()
    selection_rng = np.random.default_rng(seed + 4100)
    random_features = selection_rng.choice(remaining, 2, replace=False).tolist()
    selected_features = low + high + random_features

    x = frame("X_train_selected.parquet")
    test = frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[x.index]
    bundle = joblib.load(TCGA / "rsf_baseline.joblib")
    features = list(bundle["features"])
    if features != list(x.columns):
        raise RuntimeError("Saved RSF feature order changed")
    phase3 = phase3_helpers()
    external, _, _ = phase3.external_columns(features)
    active_features = [feature for feature in features if feature in external]
    active_indices = np.array([features.index(feature) for feature in active_features])
    target_local_indices = np.array([active_features.index(feature) for feature in selected_features])
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
    query_rng = np.random.default_rng(seed + 4000)
    source_selection = np.sort(query_rng.choice(len(source_pool), 12, replace=False))
    external_selection = np.sort(query_rng.choice(len(external_aligned), 12, replace=False))
    queries = [source_pool[source_selection], external_aligned[external_selection]]

    unique_train = x.loc[~train_patients.duplicated()]
    background_rng = np.random.default_rng(seed + 4200)
    background_selection = np.sort(
        background_rng.choice(len(unique_train), background_patients, replace=False)
    )
    sampled_backgrounds = unique_train.iloc[background_selection].to_numpy(dtype=np.float32)
    mean_background = train_mean[None, :]
    models, training, training_seconds = train_replicas(x, y, 2, seed)

    configurations = {
        "mean_candidate_8": (mean_background, 8, 10000),
        "mean_reference_256": (mean_background, 256, 20000),
        "sampled_candidate_8": (sampled_backgrounds, 8, 10000),
        "sampled_reference_128": (sampled_backgrounds, 128, 20000),
    }
    values: dict[str, list[np.ndarray]] = {name: [] for name in configurations}
    timings = {name: 0.0 for name in configurations}
    for site, query in enumerate(queries):
        per_site = {name: [] for name in configurations}
        for replicate, model in enumerate(models):
            for name, (background, budget, seed_offset) in configurations.items():
                timer = perf_counter()
                attribution = target_permutation_attributions(
                    model,
                    query,
                    background,
                    active_mean,
                    active_indices,
                    target_local_indices,
                    permutations=budget,
                    seed=seed + seed_offset + 100 * site + replicate,
                )
                timings[name] += perf_counter() - timer
                per_site[name].append(attribution)
        for name in configurations:
            values[name].append(np.asarray(per_site[name]))

    mean_candidate = values["mean_candidate_8"]
    mean_reference_256 = values["mean_reference_256"]
    sampled_candidate = values["sampled_candidate_8"]
    sampled_reference = values["sampled_reference_128"]
    # The first 128 RNG permutations of the 256 run exactly reproduce the
    # pilot reference seed stream; compute them explicitly for the comparison.
    mean_reference_128 = []
    timer = perf_counter()
    for site, query in enumerate(queries):
        rows = []
        for replicate, model in enumerate(models):
            rows.append(
                target_permutation_attributions(
                    model,
                    query,
                    mean_background,
                    active_mean,
                    active_indices,
                    target_local_indices,
                    permutations=128,
                    seed=seed + 20000 + 100 * site + replicate,
                )
            )
        mean_reference_128.append(np.asarray(rows))
    timings["mean_reference_128"] = perf_counter() - timer

    f_mean_128 = selected_f(mean_candidate, mean_reference_128)
    f_mean_256 = selected_f(mean_candidate, mean_reference_256)
    f_sampled_128 = selected_f(sampled_candidate, sampled_reference)
    g_mean = selected_g_claim(mean_candidate)
    g_sampled = selected_g_claim(sampled_candidate)
    candidate_mean_abs_delta = np.mean(
        np.abs(np.concatenate(sampled_candidate, axis=1) - np.concatenate(mean_candidate, axis=1)), axis=(0, 1)
    )
    candidate_max_abs_delta = np.max(
        np.abs(np.concatenate(sampled_candidate, axis=1) - np.concatenate(mean_candidate, axis=1)), axis=(0, 1)
    )
    reference_mean_abs_delta = np.mean(
        np.abs(np.concatenate(sampled_reference, axis=1) - np.concatenate(mean_reference_128, axis=1)), axis=(0, 1)
    )
    reference_max_abs_delta = np.max(
        np.abs(np.concatenate(sampled_reference, axis=1) - np.concatenate(mean_reference_128, axis=1)), axis=(0, 1)
    )
    result = pd.DataFrame(
        {
            "feature": selected_features,
            "selection_reason": ["low_WI"] * 3 + ["high_WI"] * 3 + ["seeded_random"] * 2,
            "pilot_WI": base_report.set_index("feature").loc[selected_features, "WI"].to_numpy(),
            "candidate_attr_mean_abs_delta_sampled_vs_mean": candidate_mean_abs_delta,
            "candidate_attr_max_abs_delta_sampled_vs_mean": candidate_max_abs_delta,
            "reference_attr_mean_abs_delta_sampled_vs_mean": reference_mean_abs_delta,
            "reference_attr_max_abs_delta_sampled_vs_mean": reference_max_abs_delta,
            "F_mean_bg_ref128": f_mean_128,
            "F_sampled_bg_ref128": f_sampled_128,
            "F_background_delta": f_sampled_128 - f_mean_128,
            "G_claim_mean_bg": g_mean,
            "G_claim_sampled_bg": g_sampled,
            "G_claim_background_delta": g_sampled - g_mean,
            "F_mean_bg_ref256": f_mean_256,
            "F_budget_delta_256_minus_128": f_mean_256 - f_mean_128,
            "G_claim_budget_delta_256_minus_128": np.zeros(len(selected_features)),
        }
    )
    result.to_csv(OUT / "sensitivity_report.csv", index=False)
    max_f_background = float(result.F_background_delta.abs().max())
    max_g_background = float(result.G_claim_background_delta.abs().max())
    max_f_budget = float(result.F_budget_delta_256_minus_128.abs().max())
    background_small = max(max_f_background, max_g_background) < 0.01
    budget_small = max_f_budget < 0.01
    original = json.loads((OUT / "metadata.json").read_text())
    base_timing = original["timing_seconds"]
    sampled_multiplier = background_patients
    reference_multiplier = 2
    recommendation = (
        "77-feature fixed-background subgame and 8/128 budget are acceptable for N=50"
        if background_small and budget_small
        else (
            "Use sampled background (k=16) and 8/128 for N=50"
            if not background_small and budget_small
            else "Use sampled background (k=16) and 8/256 for N=50"
        )
    )
    measured_candidate_multiplier = timings["sampled_candidate_8"] / timings["mean_candidate_8"]
    measured_reference_multiplier = timings["sampled_reference_128"] / timings["mean_reference_128"]
    revised_n50_seconds = base_timing["n50_same_12x12_query_seconds"]
    conservative_n50_upper_seconds = revised_n50_seconds
    if not background_small:
        revised_n50_seconds = 25 * (
            base_timing["training_n2"]
            + measured_candidate_multiplier * base_timing["candidate_explanations_n2"]
            + measured_reference_multiplier * base_timing["reference_explanations_n2"]
        )
        conservative_n50_upper_seconds = 25 * (
            base_timing["training_n2"]
            + sampled_multiplier * base_timing["candidate_explanations_n2"]
            + sampled_multiplier * base_timing["reference_explanations_n2"]
        )
    if not budget_small:
        revised_n50_seconds += 25 * reference_multiplier * base_timing["reference_explanations_n2"]
        conservative_n50_upper_seconds += 25 * reference_multiplier * base_timing["reference_explanations_n2"]
    sensitivity_metadata = {
        "phase": "4.1_sensitivity",
        "arguments": {
            "replicas": 2,
            "query_rows_per_site": 12,
            "selected_features": len(selected_features),
            "background_patients": background_patients,
            "candidate_permutations": 8,
            "reference_budgets": [128, 256],
            "seed": seed,
        },
        "selected_features": selected_features,
        "training": training,
        "timing_seconds": {"training": training_seconds, **timings, "total": perf_counter() - started},
        "maximum_absolute_score_changes": {
            "F_background": max_f_background,
            "G_claim_background": max_g_background,
            "F_reference_128_to_256": max_f_budget,
            "G_claim_reference_128_to_256": 0.0,
        },
        "decision_threshold": "maximum absolute F and G_claim score change < 0.01",
        "background_small": background_small,
        "reference_budget_small": budget_small,
        "recommendation": recommendation,
        "projected_n50_seconds_under_recommendation": revised_n50_seconds,
        "projected_n50_conservative_linear_upper_seconds": conservative_n50_upper_seconds,
        "measured_sampled_over_mean_runtime_multipliers": {
            "candidate": measured_candidate_multiplier,
            "reference": measured_reference_multiplier,
        },
        "method_note": "Frozen G_claim uses candidate attributions only; changing only reference budget cannot change G_claim.",
        "source_sha256": {"experiments/07_rsf_pilot.py": sha(Path(__file__))},
    }
    (OUT / "sensitivity_metadata.json").write_text(
        json.dumps(sensitivity_metadata, indent=2), encoding="utf-8"
    )
    original["source_sha256"]["experiments/07_rsf_pilot.py"] = sha(Path(__file__))
    (OUT / "metadata.json").write_text(json.dumps(original, indent=2), encoding="utf-8")

    base_note = base_note_path.read_text(encoding="utf-8").split("\n## Faz 4.1 — sağlamlık kontrolü")[0]
    note = base_note.rstrip() + "\n\n## Faz 4.1 — sağlamlık kontrolü\n\n"
    note += (
        f"Aynı N=2 modeller, aynı 12+12 query, aynı seed akışı ve sekiz feature kullanıldı: üç düşük WI, "
        f"üç yüksek WI, iki deterministik rastgele feature. Background duyarlılığında pasif 539 girdi "
        f"{background_patients} benzersiz TCGA train hastasından alındı ve her coalition model çıktısı bu "
        "hastalar üzerinden ortalandı.\n\n"
    )
    note += "### Background sabitleme duyarlılığı\n\n" + markdown_table(
        result[
            [
                "feature",
                "selection_reason",
                "candidate_attr_mean_abs_delta_sampled_vs_mean",
                "reference_attr_mean_abs_delta_sampled_vs_mean",
                "F_background_delta",
                "G_claim_background_delta",
            ]
        ]
    ) + "\n\n"
    note += (
        f"Maksimum |ΔF|={max_f_background:.6f}, maksimum |ΔG_claim|={max_g_background:.6f}. "
        "Attribution farkları risk-skoru birimindedir; F/G farkları [0,1] skor birimindedir.\n\n"
    )
    note += "### Referans bütçe yakınsaklığı\n\n" + markdown_table(
        result[
            [
                "feature",
                "F_mean_bg_ref128",
                "F_mean_bg_ref256",
                "F_budget_delta_256_minus_128",
                "G_claim_budget_delta_256_minus_128",
            ]
        ]
    ) + "\n\n"
    note += (
        f"Maksimum |ΔF|={max_f_budget:.6f}. G_claim dondurulmuş tanımda candidate attribution'dan "
        "hesaplanır; candidate=8 sabitken yalnız referansı 128'den 256'ya çıkarmak G_claim'i değiştirmez, "
        "bu nedenle ΔG_claim=0'dır.\n\n"
    )
    note += "### Öneri\n\n"
    note += f"**{recommendation}.** Karar eşiği her kontrolde maksimum mutlak skor farkı <0.01 idi. "
    note += (
        f"Ölçülen sampled/mean süre oranları candidate için {measured_candidate_multiplier:.2f}×, referans "
        f"için {measured_reference_multiplier:.2f}×. Bu oranlarla kaba N=50 maliyet projeksiyonu "
        f"{revised_n50_seconds:.1f} sn ({revised_n50_seconds / 3600:.2f} saat); background maliyetinin "
        f"k=16 ile tam doğrusal arttığı muhafazakâr üst projeksiyon {conservative_n50_upper_seconds:.1f} sn "
        f"({conservative_n50_upper_seconds / 3600:.2f} saat). N=50 bu kontrolde başlatılmadı; ayrıca onay "
        "bekleniyor.\n"
    )
    base_note_path.write_text(note, encoding="utf-8")
    (OUT / "design_note.md").write_text(note, encoding="utf-8")
    print(result.to_string(index=False), flush=True)
    print(json.dumps(sensitivity_metadata, indent=2), flush=True)


def run_full_n20(*, seed: int = 20260917) -> None:
    """Run the approved 77-feature, sampled-background N=20 RSF audit."""
    replicas = 20
    anchors_per_site = 12
    background_patients = 16
    candidate_budget = 8
    reference_budget = 128
    invocation_started = perf_counter()
    OUT_FULL.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = OUT_FULL / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    x = frame("X_train_selected.parquet")
    test = frame("X_test_selected.parquet")
    y = frame("y_train.parquet").loc[x.index]
    bundle = joblib.load(TCGA / "rsf_baseline.joblib")
    baseline_model = bundle["model"]
    features = list(bundle["features"])
    if features != list(x.columns) or baseline_model.n_features_in_ != 616:
        raise RuntimeError("Saved RSF is incompatible with the frozen 616-feature panel")
    if not (
        baseline_model.n_estimators == 500
        and baseline_model.min_samples_split == 10
        and baseline_model.min_samples_leaf == 5
        and baseline_model.max_features == "sqrt"
        and baseline_model.random_state == 42
    ):
        raise RuntimeError("Saved RSF hyperparameters changed")
    del baseline_model

    phase3 = phase3_helpers()
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
    query_rng = np.random.default_rng(seed + 4000)
    source_selection = np.sort(query_rng.choice(len(source_pool), anchors_per_site, replace=False))
    external_selection = np.sort(query_rng.choice(len(external_aligned), anchors_per_site, replace=False))
    queries = [source_pool[source_selection], external_aligned[external_selection]]
    unique_train = x.loc[~train_patients.duplicated()]
    background_rng = np.random.default_rng(seed + 4200)
    background_selection = np.sort(
        background_rng.choice(len(unique_train), background_patients, replace=False)
    )
    sampled_backgrounds = unique_train.iloc[background_selection].to_numpy(dtype=np.float32)
    sampled_backgrounds[:, active_indices] = active_mean.astype(np.float32)

    inputs = {
        name: sha(TCGA / name)
        for name in ["rsf_baseline.joblib", "X_train_selected.parquet", "y_train.parquet"]
    }
    query_hashes = [hashlib.sha256(np.asarray(query, dtype="<f8").tobytes()).hexdigest() for query in queries]
    protocol = {
        "inputs": inputs,
        "replicas": replicas,
        "anchors_per_site": anchors_per_site,
        "background_patients": background_patients,
        "candidate_budget": candidate_budget,
        "reference_budget": reference_budget,
        "seed": seed,
        "query_hashes": query_hashes,
        "background_sha256": hashlib.sha256(np.asarray(sampled_backgrounds, dtype="<f4").tobytes()).hexdigest(),
        "protocol": "rsf-patient-cluster-sampled-background-k16-permutation-8-128-v1",
    }
    protocol_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    patient_groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    bootstrap_rng = np.random.default_rng(seed)
    survival = Surv.from_arrays(event=y.event.astype(bool).to_numpy(), time=y.time.to_numpy())
    values = x.to_numpy(dtype=np.float32)
    candidates = np.empty((2, replicas, anchors_per_site, len(active_features)))
    references = np.empty_like(candidates)
    checkpoint_records = []
    checkpoints_reused = 0
    for replicate in range(replicas):
        indices = np.concatenate(
            [patient_groups[i] for i in bootstrap_rng.integers(0, len(patient_groups), len(patient_groups))]
        )
        paths = [
            checkpoint_dir / f"{protocol_hash}_rep{replicate:02d}_site{site}.npz"
            for site in range(2)
        ]
        metadata_paths = [path.with_suffix(".json") for path in paths]
        missing_sites = [site for site in range(2) if not paths[site].exists() or not metadata_paths[site].exists()]
        model = None
        training_record = None
        training_charge_site = missing_sites[0] if missing_sites else None
        if missing_sites:
            model = RandomSurvivalForest(
                n_estimators=500,
                min_samples_split=10,
                min_samples_leaf=5,
                max_features="sqrt",
                random_state=seed + replicate,
                n_jobs=-1,
                low_memory=True,
            )
            fit_started = perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model.fit(values[indices], survival[indices])
            training_record = {
                "rows": int(len(indices)),
                "events": int(survival[indices]["event"].sum()),
                "elapsed_seconds": perf_counter() - fit_started,
                "warning_categories": dict(Counter(w.category.__name__ for w in caught)),
            }
            if caught:
                raise RuntimeError(
                    f"RSF replicate {replicate} emitted warnings; no checkpoint exported: "
                    f"{training_record['warning_categories']}"
                )
        for site in range(2):
            if site not in missing_sites:
                with np.load(paths[site], allow_pickle=False) as saved:
                    candidates[site, replicate] = saved["candidate"]
                    references[site, replicate] = saved["reference"]
                record = json.loads(metadata_paths[site].read_text())
                if record["protocol_hash"] != protocol_hash:
                    raise RuntimeError("RSF checkpoint protocol mismatch")
                checkpoint_records.append(record)
                checkpoints_reused += 1
                print(f"RSF N=20 replicate {replicate + 1}/20 site {site}: verified checkpoint", flush=True)
                continue
            candidate_started = perf_counter()
            candidate = background_permutation_attributions(
                model,
                queries[site],
                sampled_backgrounds,
                active_mean,
                active_indices,
                permutations=candidate_budget,
                seed=seed + 10000 + 100 * site + replicate,
            )
            candidate_seconds = perf_counter() - candidate_started
            reference_started = perf_counter()
            reference = background_permutation_attributions(
                model,
                queries[site],
                sampled_backgrounds,
                active_mean,
                active_indices,
                permutations=reference_budget,
                seed=seed + 20000 + 100 * site + replicate,
            )
            reference_seconds = perf_counter() - reference_started
            candidates[site, replicate] = candidate
            references[site, replicate] = reference
            record = {
                "protocol_hash": protocol_hash,
                "replicate": replicate,
                "site": site,
                "training": training_record,
                "training_seconds_charged": (
                    training_record["elapsed_seconds"] if site == training_charge_site else 0.0
                ),
                "candidate_seconds": candidate_seconds,
                "reference_seconds": reference_seconds,
                "warning_categories": {},
            }
            np.savez_compressed(paths[site], candidate=candidate, reference=reference)
            metadata_paths[site].write_text(json.dumps(record, indent=2), encoding="utf-8")
            checkpoint_records.append(record)
            print(
                f"RSF N=20 replicate {replicate + 1}/20 site {site}: "
                f"candidate={candidate_seconds:.1f}s reference={reference_seconds:.1f}s",
                flush=True,
            )
        del model

    report = component_report(candidates, references, active_features)
    n2_report = component_report(candidates[:, :2], references[:, :2], active_features)
    report["WI_n2_sampled_background"] = report.feature.map(n2_report.set_index("feature").WI)
    report["point_adequate_n2_sampled_background"] = report.feature.map(
        n2_report.set_index("feature").point_adequate
    )
    pilot = pd.read_csv(OUT / "feature_report.csv").set_index("feature")
    report["WI_n2_fixed_background_pilot"] = report.feature.map(pilot.WI)
    report["point_adequate_n2_fixed_background_pilot"] = report.feature.map(pilot.point_adequate)
    report.to_csv(OUT_FULL / "feature_report.csv", index=False)
    site_summary = sorted(checkpoint_records, key=lambda row: (row["replicate"], row["site"]))
    pd.DataFrame(site_summary).to_json(OUT_FULL / "bootstrap_site_summary.json", orient="records", indent=2)
    bootstrap_summary = []
    for replicate in range(replicas):
        rows = [row for row in site_summary if row["replicate"] == replicate]
        attempts = []
        for row in rows:
            attempt = row.get("training")
            if attempt is not None and attempt not in attempts:
                attempts.append(attempt)
        bootstrap_summary.append(
            {
                "replicate": replicate,
                "rows": attempts[0]["rows"],
                "events": attempts[0]["events"],
                "warning_categories": {},
                "fit_attempts": attempts,
                "site_records": [
                    {
                        "site": row["site"],
                        "candidate_seconds": row["candidate_seconds"],
                        "reference_seconds": row["reference_seconds"],
                    }
                    for row in rows
                ],
                "elapsed_seconds": sum(
                    row["training_seconds_charged"] + row["candidate_seconds"] + row["reference_seconds"]
                    for row in rows
                ),
            }
        )
    pd.DataFrame(bootstrap_summary).to_json(
        OUT_FULL / "bootstrap_summary.json", orient="records", indent=2
    )
    mapping[mapping.feature.isin(active_features)].to_csv(OUT_FULL / "feature_mapping.csv", index=False)
    np.savez_compressed(
        OUT_FULL / "attributions_n20.npz",
        candidate=candidates,
        reference=references,
        features=np.asarray(active_features),
    )

    aggregate = report.groupby("block").agg(
        features=("feature", "size"),
        assessed=("WI", "count"),
        mean_F=("F", "mean"),
        mean_S=("S", "mean"),
        mean_U=("U", "mean"),
        mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"),
        mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.to_csv(OUT_FULL / "summary.csv", index=False)
    n2_passing = set(n2_report.loc[n2_report.point_adequate == 1, "feature"])
    n20_passing = set(report.loc[report.point_adequate == 1, "feature"])
    comparison = report[
        [
            "feature",
            "block",
            "WI_n2_fixed_background_pilot",
            "WI_n2_sampled_background",
            "WI",
            "point_adequate_n2_fixed_background_pilot",
            "point_adequate_n2_sampled_background",
            "point_adequate",
        ]
    ].copy()
    comparison.to_csv(OUT_FULL / "n2_n20_comparison.csv", index=False)

    phase21 = json.loads((ROOT / "results/phase21/02_scenarios_metadata.json").read_text(encoding="utf-8"))
    phase21_hashes = {Path(name).as_posix(): digest for name, digest in phase21["source_sha256"].items()}
    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    formula_unchanged = all(phase21_hashes[name] == sha(ROOT / name) for name in formula_files)
    if not formula_unchanged:
        raise RuntimeError("Frozen formula module SHA256 changed")
    accumulated_compute_seconds = sum(
        row["training_seconds_charged"] + row["candidate_seconds"] + row["reference_seconds"]
        for row in checkpoint_records
    )
    invocation_seconds = perf_counter() - invocation_started
    metadata = {
        "phase": "4.2_rsf_n20",
        "arguments": {
            "replicas": replicas,
            "anchors_per_site": anchors_per_site,
            "background_patients": background_patients,
            "candidate_permutations": candidate_budget,
            "reference_permutations": reference_budget,
            "seed": seed,
        },
        "design": {
            "background": "16 deterministic unique TCGA train patients; inactive predictions averaged",
            "active_game_features": 77,
            "inactive_fixed_feature_count": 539,
            "query_rows": {"TCGA": 12, "METABRIC": 12},
            "query_pool_sizes": {"TCGA": len(source_pool), "METABRIC": len(external_aligned)},
        },
        "model": {
            "source": str(TCGA / "rsf_baseline.joblib"),
            "main_model_retrained": False,
            "bootstrap_models_fit": replicas,
            "sha256": sha(TCGA / "rsf_baseline.joblib"),
            "bootstrap_prediction_storage": "low_memory=True",
            "low_memory_equivalence_check_max_abs_risk_difference": 7.105427357601002e-15,
            "low_memory_equivalence_check_tree_structures_identical": True,
        },
        "results": {
            "assessed_features": len(report),
            "point_adequate_n2_fixed_background_pilot": int(pilot.point_adequate.sum()),
            "point_adequate_n2_sampled_background": int(n2_report.point_adequate.sum()),
            "point_adequate_n20_sampled_background": int(report.point_adequate.sum()),
            "n2_sampled_passing_still_passing_n20": len(n2_passing & n20_passing),
            "n2_sampled_passing_now_failing_n20": len(n2_passing - n20_passing),
            "n20_newly_passing_vs_n2_sampled": len(n20_passing - n2_passing),
            "bottleneck_counts": report.bottleneck.value_counts().to_dict(),
        },
        "runtime": {
            "accumulated_measured_compute_wall_seconds": accumulated_compute_seconds,
            "completion_invocation_wall_seconds": invocation_seconds,
            "checkpoint_sites_reused_in_completion_invocation": checkpoints_reused,
            "timing_note": "Accumulated time sums measured sequential fit/candidate/reference stages across checkpoints.",
        },
        "formula_module_sha256_unchanged_from_phase21": formula_unchanged,
        "formula_module_sha256": {name: sha(ROOT / name) for name in formula_files},
        "fixed_calibration": CALIBRATION,
        "external_inventory": external_info,
        "versions": {
            name: importlib.metadata.version(name)
            for name in ["numpy", "pandas", "scikit-learn", "scikit-survival", "shap"]
        },
        "python": platform.python_version(),
        "source_sha256": {"experiments/07_rsf_pilot.py": sha(Path(__file__))},
        "input_sha256": inputs,
        "protocol_sha256": protocol_hash,
    }
    (OUT_FULL / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    retained = metadata["results"]["n2_sampled_passing_still_passing_n20"]
    lost = metadata["results"]["n2_sampled_passing_now_failing_n20"]
    gained = metadata["results"]["n20_newly_passing_vs_n2_sampled"]
    winner_pattern = "benzer bir eşik gerilemesi görüldü" if lost > gained else "aynı yönde belirgin bir gerileme görülmedi"
    report_text = "# Faz 4.2 — RSF N=20 tam 77-feature koşu raporu\n\n"
    report_text += (
        "Onaylanan tasarım uygulandı: sampled background k=16, candidate/reference bütçesi 8/128, "
        "N=20 patient-cluster bootstrap ve pilotla aynı 12 TCGA + 12 METABRIC query. Ana "
        "`rsf_baseline.joblib` yeniden eğitilmedi; 20 bootstrap modeli ayrıca fit edildi. Tahmin "
        "maliyetini azaltmak için bootstrap modellerinde scikit-survival `low_memory=True` kullanıldı: "
        "aynı seed/veriyle 500 ağaç yapısı birebir aynı, standart ve low-memory risk skorları arasındaki "
        "maksimum fark 7.11e-15 bulundu. Bu bir model/metodoloji değişikliği değil, çıktı saklama "
        "optimizasyonudur.\n\n"
    )
    report_text += "## Block bazında agregat\n\n" + markdown_table(aggregate) + "\n\n"
    report_text += (
        f"{len(report)} feature değerlendirildi; {int(report.point_adequate.sum())} feature WI>0.8 koşulunu "
        f"geçti. Bottleneck dağılımı `{report.bottleneck.value_counts().to_dict()}`.\n\n"
    )
    report_text += "## N=2 pilot ve N=20 karşılaştırması\n\n"
    report_text += (
        f"Eski sabit-background N=2 pilotunda 43/77 feature geçmişti; Faz 4.1 bu background'ı reddettiği "
        f"için doğrudan replika etkisi karşılaştırması N=20 dizisinin ilk iki replikasıyla aynı sampled "
        f"background altında yapıldı. Sampled-background N=2'de {int(n2_report.point_adequate.sum())}/77, "
        f"N=20'de {int(report.point_adequate.sum())}/77 geçti. N=2 geçenlerin {retained} tanesi kaldı, "
        f"{lost} tanesi düştü; {gained} feature yeni geçti. CoxPH Faz 3.2'deki winner's-curse deseniyle "
        f"karşılaştırıldığında {winner_pattern}.\n\n"
    )
    report_text += "Feature bazındaki üçlü karşılaştırma `n2_n20_comparison.csv` içindedir.\n\n"
    report_text += "## Sınırlamalar ve süre\n\n"
    report_text += (
        "Query kümesi yalnız 12+12 gerçek vektördür; CoxPH auditindeki 141 TCGA + 1466 METABRIC "
        "havuzundan çok daha dardır. Bu nedenle özellikle G/G_signed örnekleme oynaklığı daha yüksek "
        "olabilir. Sonuç tam 616-feature Shapley oyunu değildir: 77 ortak feature coalition'a girer, "
        "diğer 539 girdi 16 TCGA background hastası üzerinden entegre edilir. Toleranslar kalibre "
        "edilmemiştir; klinik geçerlilik iddiası yoktur.\n\n"
    )
    report_text += (
        f"Ölçülen birikimli eğitim+açıklama wall-clock süresi {accumulated_compute_seconds:.1f} sn "
        f"({accumulated_compute_seconds / 3600:.2f} saat). Completion invocation süresi "
        f"{invocation_seconds:.1f} sn; {checkpoints_reused} site checkpoint'i yeniden kullanıldı. "
        "Formül modüllerinin Faz 2.1 SHA256 değerleri değişmedi. N=50 çalıştırılmadı.\n\n"
    )
    report_text += "## Tüm feature sonuçları\n\n" + markdown_table(
        report[["feature", "block", "F", "S", "U", "G_claim", "R", "WI", "bottleneck"]]
    ) + "\n"
    (ROOT / "PHASE4_FULL_RUN_REPORT.md").write_text(report_text, encoding="utf-8")
    (OUT_FULL / "report.md").write_text(report_text, encoding="utf-8")
    design_path = ROOT / "PHASE4_DESIGN_NOTE.md"
    design = design_path.read_text(encoding="utf-8").split("\n## Faz 4.2 — N=20 sonucu")[0].rstrip()
    design += "\n\n## Faz 4.2 — N=20 sonucu\n\n"
    design += (
        f"Onaylanan sampled-background k=16 ve 8/128 tasarımla 77 feature/N=20 koşusu tamamlandı. "
        f"{int(report.point_adequate.sum())}/77 feature WI>0.8 geçti; bottleneck dağılımı "
        f"`{report.bottleneck.value_counts().to_dict()}`. Ayrıntılar `PHASE4_FULL_RUN_REPORT.md` içinde. "
        "N=50 başlatılmadı.\n"
    )
    design_path.write_text(design, encoding="utf-8")
    print(aggregate.to_string(index=False), flush=True)
    print(json.dumps(metadata["results"], indent=2), flush=True)
    print(json.dumps(metadata["runtime"], indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicas", type=int, default=2)
    parser.add_argument("--anchors-per-site", type=int, default=12)
    parser.add_argument("--candidate-budget", type=int, default=8)
    parser.add_argument("--reference-budget", type=int, default=128)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--sensitivity", action="store_true")
    parser.add_argument("--full-n20", action="store_true")
    args = parser.parse_args()
    if args.sensitivity and args.full_n20:
        parser.error("choose only one of --sensitivity and --full-n20")
    if args.full_n20:
        run_full_n20(seed=args.seed)
    elif args.sensitivity:
        run_sensitivity(seed=args.seed)
    else:
        run(
            replicas=args.replicas,
            anchors_per_site=args.anchors_per_site,
            candidate_budget=args.candidate_budget,
            reference_budget=args.reference_budget,
            seed=args.seed,
        )
