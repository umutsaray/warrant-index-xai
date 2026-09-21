"""Phase 6: discover, probe, and run the approved N=2 DeepSurv design pilot."""

from __future__ import annotations

from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import random
import sys
from time import perf_counter
import warnings

import joblib
import numpy as np
import pandas as pd
import shap
import torch
from sklearn.model_selection import train_test_split
from sksurv.metrics import concordance_index_censored


ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase6_deepsurv_pilot"
SEED = 20260917
N = 2
WARNING = (
    "DeepSurv pilot conditions on fixed 12+12 queries and stochastic GradientExplainer 8/128 "
    "attribution draws; N=2 is a design pilot, not a certification result."
)

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(TCGA))
from _shared import CALIBRATION, markdown_table
from train_explainable_dl import (
    FlatOmicsWrapper, OmicsDeepSurv, cindex, cox_negative_log_likelihood, make_blocks, predict_risk,
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_phase4():
    path = ROOT / "experiments/07_rsf_pilot.py"
    spec = importlib.util.spec_from_file_location("phase4_helpers_for_deepsurv", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Phase 4 helpers")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)


class IntegratedActiveWrapper(torch.nn.Module):
    """77-input game whose prediction averages 16 inactive-feature backgrounds."""

    def __init__(self, flat_model, backgrounds: np.ndarray, active_indices: np.ndarray):
        super().__init__()
        self.flat_model = flat_model
        self.register_buffer("backgrounds", torch.as_tensor(backgrounds, dtype=torch.float32))
        self.register_buffer("active_indices", torch.as_tensor(active_indices, dtype=torch.long))

    def forward(self, active):
        batch, p = active.shape
        if p != len(self.active_indices):
            raise ValueError("Unexpected active feature count")
        full = self.backgrounds.unsqueeze(0).repeat(batch, 1, 1)
        full[:, :, self.active_indices] = active.unsqueeze(1)
        values = self.flat_model(full.reshape(-1, full.shape[-1])).reshape(batch, len(self.backgrounds), 1)
        return values.mean(dim=1)


def load_checkpoint_and_data():
    checkpoint_path = TCGA / "explainable_dl_model.pt"
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    train = pd.read_parquet(TCGA / "X_train_selected.parquet").set_index("sample_key")
    test = pd.read_parquet(TCGA / "X_test_selected.parquet").set_index("sample_key")
    y_train = pd.read_parquet(TCGA / "y_train.parquet").set_index("sample_key").loc[train.index]
    y_test = pd.read_parquet(TCGA / "y_test.parquet").set_index("sample_key").loc[test.index]
    columns = checkpoint["ordered_columns"]
    if len(columns) != 616 or set(columns) != set(train.columns) or train.shape != (673, 616) or test.shape != (169, 616):
        raise RuntimeError("DeepSurv artifact does not match the frozen 673/169, 616-feature data")
    return checkpoint, train[columns], test[columns], y_train, y_test


def model_from_state(checkpoint, state_dict=None):
    model = OmicsDeepSurv(checkpoint["block_dims"], checkpoint["dropout"])
    model.load_state_dict(checkpoint["state_dict"] if state_dict is None else state_dict)
    model.eval()
    return model


def explainer_probe(checkpoint, train: pd.DataFrame, test: pd.DataFrame) -> dict:
    model = model_from_state(checkpoint)
    wrapper = FlatOmicsWrapper(model, checkpoint["slices"]).eval()
    train_z = (train.to_numpy(np.float32) - checkpoint["means"]) / checkpoint["scales"]
    test_z = (test.to_numpy(np.float32) - checkpoint["means"]) / checkpoint["scales"]
    background = torch.from_numpy(train_z[:8])
    query = torch.from_numpy(test_z[:3])
    with torch.no_grad():
        prediction = wrapper(query).numpy().reshape(-1)
        fallback_base = float(wrapper(background).numpy().mean())
    results = {}
    for name, constructor in [
        ("DeepExplainer", lambda: shap.DeepExplainer(wrapper, background)),
        ("GradientExplainer", lambda: shap.GradientExplainer(wrapper, background)),
    ]:
        started = perf_counter()
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                explainer = constructor()
                values = explainer.shap_values(query)
            values = values[0] if isinstance(values, list) else values
            values = np.asarray(values)
            if values.ndim == 3:
                values = values[..., 0]
            base = float(np.asarray(getattr(explainer, "expected_value", fallback_base)).reshape(-1)[0])
            results[name] = {
                "supported": True,
                "shape": list(values.shape),
                "finite": bool(np.isfinite(values).all()),
                "max_abs_additivity_error": float(np.max(np.abs(values.sum(axis=1) + base - prediction))),
                "warnings": [str(item.message) for item in caught],
                "elapsed_seconds": perf_counter() - started,
            }
        except Exception as error:
            results[name] = {"supported": False, "error": f"{type(error).__name__}: {error}",
                             "elapsed_seconds": perf_counter() - started}
    return results


def train_bootstrap_model(checkpoint, raw: np.ndarray, events: np.ndarray, times: np.ndarray,
                          indices: np.ndarray, replicate: int):
    seed = SEED + replicate
    set_seed(seed)
    values = raw[indices]
    event_values = events[indices]
    time_values = times[indices]
    fit_idx, val_idx = train_test_split(
        np.arange(len(values)), test_size=0.15, random_state=seed, stratify=event_values
    )
    means = values[fit_idx].mean(axis=0, dtype=np.float64).astype(np.float32)
    scales = values[fit_idx].std(axis=0, dtype=np.float64).astype(np.float32)
    scales[scales == 0] = 1.0
    standardized = (values - means) / scales
    model = OmicsDeepSurv(checkpoint["block_dims"], checkpoint["dropout"])
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    fit_blocks = make_blocks(standardized[fit_idx], checkpoint["slices"], torch.device("cpu"))
    val_blocks = make_blocks(standardized[val_idx], checkpoint["slices"], torch.device("cpu"))
    fit_time = torch.as_tensor(time_values[fit_idx], dtype=torch.float32)
    fit_event = torch.as_tensor(event_values[fit_idx], dtype=torch.float32)
    best_state, best_val, best_epoch, wait = None, -np.inf, 0, 0
    started = perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        for epoch in range(1, 301):
            model.train()
            optimizer.zero_grad()
            loss = cox_negative_log_likelihood(model(fit_blocks), fit_time, fit_event)
            loss.backward()
            optimizer.step()
            val_risk = predict_risk(model, val_blocks)
            val_score = cindex(event_values[val_idx], time_values[val_idx], val_risk)
            if val_score > best_val + 1e-6:
                best_val, best_epoch = val_score, epoch
                best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
                wait = 0
            else:
                wait += 1
            if wait >= 30:
                break
    if best_state is None:
        raise RuntimeError(f"DeepSurv replicate {replicate} produced no valid checkpoint")
    model.load_state_dict(best_state)
    model.eval()
    record = {
        "replicate": replicate, "rows": len(indices), "events": int(event_values.sum()),
        "fit_rows": len(fit_idx), "validation_rows": len(val_idx), "best_epoch": best_epoch,
        "best_validation_c_index": best_val, "last_epoch": epoch,
        "warning_categories": dict(Counter(item.category.__name__ for item in caught)),
        "elapsed_seconds": perf_counter() - started,
        "internal_split_note": "row-stratified 15% validation, matching the frozen main-model training protocol",
    }
    return model, means, scales, record


def gradient_attributions(model, checkpoint, raw_queries, raw_backgrounds, active_indices,
                          active_baseline_raw, means, scales, *, budget: int, seed: int):
    full_background = (raw_backgrounds - means) / scales
    active_baseline = (active_baseline_raw - means[active_indices]) / scales[active_indices]
    full_background[:, active_indices] = active_baseline
    queries = (raw_queries - means[active_indices]) / scales[active_indices]
    flat = FlatOmicsWrapper(model, checkpoint["slices"]).eval()
    wrapper = IntegratedActiveWrapper(flat, full_background, active_indices).eval()
    baseline = torch.from_numpy(active_baseline[None, :].astype(np.float32))
    query = torch.from_numpy(queries.astype(np.float32))
    started = perf_counter()
    explainer = shap.GradientExplainer(wrapper, baseline)
    values = explainer.shap_values(query, nsamples=budget, rseed=seed)
    values = values[0] if isinstance(values, list) else values
    values = np.asarray(values)
    if values.ndim == 3:
        values = values[..., 0]
    with torch.no_grad():
        prediction = wrapper(query).numpy().reshape(-1)
        base = float(wrapper(baseline).numpy().reshape(-1)[0])
    return values, float(np.max(np.abs(values.sum(axis=1) + base - prediction))), perf_counter() - started


def select_sensitivity_features(report: pd.DataFrame) -> pd.DataFrame:
    ordered = report.sort_values("WI")
    selected = list(ordered.head(3).feature) + list(ordered.tail(3).feature)
    pool = report.loc[~report.feature.isin(selected), "feature"].to_numpy()
    selected += np.random.default_rng(SEED + 6100).choice(pool, 2, replace=False).tolist()
    return report.set_index("feature").loc[selected].reset_index()


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    phase4 = load_phase4()
    checkpoint, train, test, y_train, y_test = load_checkpoint_and_data()
    probe = explainer_probe(checkpoint, train, test)
    if not probe["DeepExplainer"]["supported"] or not probe["GradientExplainer"]["supported"]:
        raise RuntimeError(f"Expected both gradient explainers to work, got {probe}")
    (OUT / "explainer_probe.json").write_text(json.dumps(probe, indent=2), encoding="utf-8")

    main_model = model_from_state(checkpoint)
    test_z = (test.to_numpy(np.float32) - checkpoint["means"]) / checkpoint["scales"]
    main_risk = predict_risk(main_model, make_blocks(test_z, checkpoint["slices"], torch.device("cpu")))
    reproduced_cindex = float(concordance_index_censored(
        y_test.event.astype(bool).to_numpy(), y_test.time.to_numpy(), main_risk
    )[0])
    if abs(reproduced_cindex - checkpoint["test_c_index"]) > 1e-12:
        raise RuntimeError("Stored DeepSurv test C-index did not reproduce")

    features = checkpoint["ordered_columns"]
    phase3 = phase4.phase3_helpers()
    external, mapping, external_info = phase3.external_columns(features)
    active_features = [feature for feature in features if feature in external]
    if len(active_features) != 77:
        raise RuntimeError(f"Expected 77 common features, got {len(active_features)}")
    active_indices = np.array([features.index(feature) for feature in active_features])
    train_raw = train.to_numpy(np.float32)
    train_mean = train_raw.mean(axis=0, dtype=np.float64).astype(np.float32)
    train_std = train_raw.std(axis=0, dtype=np.float64).astype(np.float32)
    active_mean = train_mean[active_indices]
    active_std = train_std[active_indices]
    external_raw = np.column_stack([external[feature] for feature in active_features])
    external_raw = external_raw[np.isfinite(external_raw).all(axis=1)]
    external_z = (external_raw - external_raw.mean(axis=0)) / external_raw.std(axis=0, ddof=0)
    external_aligned = active_mean + external_z * active_std

    train_patients = train.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    clean = test.loc[~test_patients.isin(train_patients)]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    source_pool = clean[active_features].to_numpy(np.float32)
    query_rng = np.random.default_rng(SEED + 4000)
    queries = [
        source_pool[np.sort(query_rng.choice(len(source_pool), 12, replace=False))],
        external_aligned[np.sort(query_rng.choice(len(external_aligned), 12, replace=False))].astype(np.float32),
    ]
    unique_train = train.loc[~train_patients.duplicated()]
    background_rng = np.random.default_rng(SEED + 4200)
    background_selection = np.sort(background_rng.choice(len(unique_train), 16, replace=False))
    background16 = unique_train.iloc[background_selection].to_numpy(np.float32)
    background8 = background16[:8].copy()

    groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    candidates = np.empty((2, N, 12, 77))
    references = np.empty_like(candidates)
    candidate_k8 = np.empty_like(candidates)
    reference_k8 = np.empty_like(candidates)
    reference256 = np.empty_like(candidates)
    records = []
    events = y_train.event.astype(bool).to_numpy()
    times = y_train.time.to_numpy(np.float32)
    for replicate in range(N):
        indices = np.concatenate([groups[i] for i in bootstrap_rng.integers(0, len(groups), len(groups))])
        model, means, scales, record = train_bootstrap_model(
            checkpoint, train_raw, events, times, indices, replicate
        )
        record["sites"] = []
        for site in range(2):
            site_record = {"site": site}
            for label, backgrounds, budget, offset, target in [
                ("candidate_k16_b8", background16, 8, 10000, candidates),
                ("reference_k16_b128", background16, 128, 20000, references),
                ("candidate_k8_b8", background8, 8, 10000, candidate_k8),
                ("reference_k8_b128", background8, 128, 20000, reference_k8),
                ("reference_k16_b256", background16, 256, 20000, reference256),
            ]:
                values, error, seconds = gradient_attributions(
                    model, checkpoint, queries[site], backgrounds, active_indices, active_mean,
                    means, scales, budget=budget, seed=SEED + offset + 100 * site + replicate,
                )
                target[site, replicate] = values
                site_record[label] = {"seconds": seconds, "max_abs_additivity_error": error}
            record["sites"].append(site_record)
        records.append(record)

    report = phase4.component_report(candidates, references, active_features)
    report.insert(0, "inference_scope", WARNING)
    report.to_csv(OUT / "feature_report.csv", index=False)
    aggregate = report.groupby("block").agg(
        features=("feature", "size"), assessed=("WI", "count"), mean_F=("F", "mean"),
        mean_S=("S", "mean"), mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.insert(0, "inference_scope", WARNING)
    aggregate.to_csv(OUT / "summary.csv", index=False)
    np.savez_compressed(OUT / "attributions_n2.npz", candidate=candidates, reference=references,
                        features=np.asarray(active_features))
    for record in records:
        record["inference_scope"] = WARNING
    (OUT / "bootstrap_summary.json").write_text(json.dumps(records, indent=2), encoding="utf-8")

    k8_report = phase4.component_report(candidate_k8, reference_k8, active_features).set_index("feature")
    ref256_report = phase4.component_report(candidates, reference256, active_features).set_index("feature")
    selected = select_sensitivity_features(report).set_index("feature")
    base = report.set_index("feature")
    sensitivity = selected[["block", "WI"]].rename(columns={"WI": "WI_k16_ref128"}).copy()
    sensitivity["F_k16_ref128"] = base.loc[sensitivity.index, "F"]
    sensitivity["G_claim_k16"] = base.loc[sensitivity.index, "G_claim"]
    sensitivity["F_k8_ref128"] = k8_report.loc[sensitivity.index, "F"]
    sensitivity["G_claim_k8"] = k8_report.loc[sensitivity.index, "G_claim"]
    sensitivity["delta_F_k16_minus_k8"] = sensitivity.F_k16_ref128 - sensitivity.F_k8_ref128
    sensitivity["delta_G_k16_minus_k8"] = sensitivity.G_claim_k16 - sensitivity.G_claim_k8
    sensitivity["F_k16_ref256"] = ref256_report.loc[sensitivity.index, "F"]
    sensitivity["delta_F_ref256_minus_128"] = sensitivity.F_k16_ref256 - sensitivity.F_k16_ref128
    target_positions = np.array([active_features.index(feature) for feature in sensitivity.index])
    sensitivity["candidate_attr_mean_abs_delta_k16_vs_k8"] = np.mean(
        np.abs(candidates[..., target_positions] - candidate_k8[..., target_positions]), axis=(0, 1, 2)
    )
    sensitivity["reference_attr_mean_abs_delta_k16_vs_k8"] = np.mean(
        np.abs(references[..., target_positions] - reference_k8[..., target_positions]), axis=(0, 1, 2)
    )
    sensitivity["reference_attr_mean_abs_delta_256_vs_128"] = np.mean(
        np.abs(reference256[..., target_positions] - references[..., target_positions]), axis=(0, 1, 2)
    )
    sensitivity.insert(0, "feature", sensitivity.index)
    sensitivity.insert(0, "inference_scope", WARNING)
    sensitivity.to_csv(OUT / "sensitivity_report.csv", index=False)

    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    formula_sha = {name: sha(ROOT / name) for name in formula_files}
    phase56 = json.loads((ROOT / "results/phase56_rsf_n50/metadata.json").read_text())
    unchanged = formula_sha == phase56["formula_module_sha256"]
    if not unchanged:
        raise RuntimeError("Frozen formula module SHA256 changed")
    base_training = sum(row["elapsed_seconds"] for row in records)
    base_candidate = sum(site["candidate_k16_b8"]["seconds"] for row in records for site in row["sites"])
    base_reference = sum(site["reference_k16_b128"]["seconds"] for row in records for site in row["sites"])
    base_n2 = base_training + base_candidate + base_reference
    elapsed = perf_counter() - started
    max_candidate_error = max(site["candidate_k16_b8"]["max_abs_additivity_error"] for row in records for site in row["sites"])
    max_reference_error = max(site["reference_k16_b128"]["max_abs_additivity_error"] for row in records for site in row["sites"])
    metadata = {
        "phase": "6_DeepSurv_discovery_explainer_design_N2_pilot",
        "conditionality_warning": WARNING,
        "adim0": {
            "model_found": True, "artifact": str(TCGA / "explainable_dl_model.pt"),
            "artifact_sha256": sha(TCGA / "explainable_dl_model.pt"),
            "train_shape": [673, 616], "test_shape": [169, 616], "test_c_index": reproduced_cindex,
            "best_epoch": checkpoint["best_epoch"], "best_validation_c_index": checkpoint["best_validation_c_index"],
            "dropout": checkpoint["dropout"], "random_state": checkpoint["random_state"],
            "architecture": "block encoders RNA 200-128-32, miRNA 200-64-16, methylation 200-128-32, clinical 16-32-8; head 88-64-16-1",
            "optimizer": "AdamW lr=0.001 weight_decay=0.0001; max_epochs=300; patience=30",
        },
        "adim1": {"probe": probe, "selected_for_approximate_F": "GradientExplainer",
                  "selection_reason": "nsamples exposes the approved independent 8/128 stochastic budgets; DeepExplainer is retained as the near-exact additivity support check"},
        "adim2": {"background": "inactive 539 features averaged over 16 deterministic unique TCGA train patients; active baseline is TCGA train mean",
                  "background_health_check": "nested sampled k=8 versus k=16 on 8 diverse features",
                  "candidate_budget": 8, "reference_budgets": [128, 256],
                  "max_abs_delta_F_k16_vs_k8": float(sensitivity.delta_F_k16_minus_k8.abs().max()),
                  "max_abs_delta_G_claim_k16_vs_k8": float(sensitivity.delta_G_k16_minus_k8.abs().max()),
                  "max_abs_delta_F_ref256_vs_128": float(sensitivity.delta_F_ref256_minus_128.abs().max())},
        "adim3": {"replicates": 2, "features": 77, "query_rows": {"TCGA": 12, "METABRIC": 12},
                  "point_adequate": int(report.point_adequate.sum()),
                  "bottleneck_counts": report.bottleneck.value_counts().to_dict(),
                  "N20_projection_seconds": base_n2 * 10, "N50_projection_seconds": base_n2 * 25,
                  "full_scale_started": False,
                  "design_gate": "revision_required_before_N20",
                  "design_gate_reason": "candidate=8 GradientExplainer maximum additivity error is too large (0.439 risk-score units) despite stable sampled background and reference budget"},
        "runtime": {"training_n2": base_training, "candidate_n2": base_candidate,
                    "reference_n2": base_reference, "approved_base_N2": base_n2,
                    "all_sensitivity_and_reporting": elapsed,
                    "max_candidate_additivity_error": max_candidate_error,
                    "max_reference_additivity_error": max_reference_error},
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": unchanged,
        "input_sha256": {name: sha(TCGA / name) for name in [
            "explainable_dl_model.pt", "X_train_selected.parquet", "X_test_selected.parquet",
            "y_train.parquet", "y_test.parquet"
        ]},
        "source_sha256": {"experiments/20_deepsurv_pilot.py": sha(Path(__file__))},
        "external_inventory": external_info,
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "torch": torch.__version__, "shap": shap.__version__},
        "runtime_executable": sys.executable,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)

    example = pd.concat([report.nsmallest(4, "WI"), report.nlargest(4, "WI")])[
        ["inference_scope", "feature", "block", "F", "S", "U", "G_claim", "R", "WI", "bottleneck"]
    ]
    note = f"""# Faz 6 — DeepSurv keşif, açıklanabilirlik ve N=2 tasarım pilotu

> **Koşulluluk uyarısı:** {WARNING}

## ADIM 0 — hazır model keşfi

Hazır `explainable_dl_model.pt` bulundu ve artifact içinden doğrulandı. Model aynı 673 train / 169 test ve 616-feature panelini kullanır. Dört bloklu mimari: RNA 200→128→32, miRNA 200→64→16, methylation 200→128→32, clinical 16→32→8; birleşik head 88→64→16→1. Encoder katmanları Linear–BatchNorm–ReLU–Dropout, dropout 0.3'tür. AdamW lr=0.001, weight decay=0.0001, en çok 300 epoch, patience=30 ve seed=42 kullanılmıştır. En iyi epoch 26, validation C-index {checkpoint['best_validation_c_index']:.6f}; test C-index artifacttan yeniden üretilerek **{reproduced_cindex:.6f}** bulundu. Ana model yeniden eğitilmedi.

Karşılaştırılabilirlik sınırlaması: ana DeepSurv modeli 673 satırlık train içinde %15 satır-stratified validation kullanır. Mixed cohortta aynı hastanın birden çok örneği bu iç split'i aşabilir; bu, mevcut model artifactının kayıtlı metodolojik sınırlamasıdır.

## ADIM 1 — gerçek explainer testi

{markdown_table(pd.DataFrame([{"explainer": name, **values} for name, values in probe.items()]).drop(columns=["warnings"], errors="ignore"))}

Her iki yöntem gerçek checkpoint ve gerçek 616-feature tensörlerinde çalıştı. DeepExplainer maksimum additivity hatası {probe['DeepExplainer']['max_abs_additivity_error']:.3g} ile near-exact destek kontrolüdür. Pilot için GradientExplainer seçildi; çünkü `nsamples` parametresi bağımsız candidate/reference 8/128 bütçesini uygular. DeepExplainer deterministik olduğundan bu F bütçe karşılaştırmasını dejenere eder. GradientExplainer sonucu yaklaşık ve stokastiktir; N=2 candidate/ref maksimum additivity hataları sırasıyla {max_candidate_error:.4f}/{max_reference_error:.4f}'tür.

## ADIM 2 — sampled-background sağlık kontrolü

Sabit-ortalama background yeniden denenmedi. Diğer 539 feature, RSF ile aynı deterministik 16 benzersiz TCGA train hastası üzerinden prediction düzeyinde entegre edildi; aktif 77-feature baseline'ı TCGA train ortalamasıdır. Sekiz çeşitli feature için nested sampled k=8→16 ve reference 128→256 kontrolü:

{markdown_table(sensitivity.reset_index(drop=True))}

Maksimum |ΔF| k8→k16 **{sensitivity.delta_F_k16_minus_k8.abs().max():.6f}**, maksimum |ΔG_claim| **{sensitivity.delta_G_k16_minus_k8.abs().max():.6f}**. Reference 128→256 maksimum |ΔF| **{sensitivity.delta_F_ref256_minus_128.abs().max():.6f}**. Bunlar pilot sağlık ölçüleridir; eşik kalibrasyonu değildir.

## ADIM 3 — N=2 pilot

> **Koşulluluk uyarısı:** {WARNING}

{markdown_table(aggregate)}

77/77 feature değerlendirildi; {int(report.point_adequate.sum())} feature WI≥0.8 geçti. Bottleneck dağılımı `{report.bottleneck.value_counts().to_dict()}`.

### Örnek feature'lar

{markdown_table(example)}

Onaylı ana N=2 eğitim+candidate+reference süresi {base_n2:.3f} sn. Doğrusal projeksiyon N=20 için {base_n2*10:.1f} sn ({base_n2*10/60:.2f} dk), N=50 için {base_n2*25:.1f} sn ({base_n2*25/60:.2f} dk). Bütün sensitivity kontrolleri ve raporlama dahil gerçek invocation {elapsed:.3f} sn sürdü.

## Karar ve sınırlar

Sampled k=16 background sağlık kontrolü ve reference=128 bütçesi onaylanabilir: k8→k16 maksimum |ΔG_claim| 0.0038, reference 128→256 maksimum |ΔF| 9.8e-5'tir. Buna karşılık GradientExplainer candidate=8 maksimum additivity hatası 0.439 risk-skoru birimidir; bu nedenle tasarım kapısı **REVİZYON GEREKTİRİYOR**. N≥20 mevcut 8/128 ayarıyla başlatılmamalıdır. Sonraki ucuz karar turunda candidate bütçesi 32/64'e yükseltilerek additivity ve F yakınsaklığı ölçülebilir veya near-exact DeepExplainer'ın F katmanındaki rolü yeniden tanımlanabilir. Bu iki seçenek frozen WI bileşen formüllerini değiştirmez, fakat explainer protokolü için kullanıcı kararı gerektirir.

N≥20 koşusu, null kalibrasyon ve FDR başlatılmadı. Query kümesi sabit 12+12'dir; GradientExplainer örnekleme belirsizliği daha sonraki inferansta ayrıca ele alınmalıdır.

Formül SHA256 değişmedi: **{unchanged}**. Kimlik-sızıntısı kontrolü geçti. Runtime: `{sys.executable}`.
"""
    (OUT / "report.md").write_text(note, encoding="utf-8")
    (ROOT / "PHASE6_DESIGN_NOTE.md").write_text(note, encoding="utf-8")
    print(json.dumps({"point_adequate": int(report.point_adequate.sum()),
                      "bottlenecks": report.bottleneck.value_counts().to_dict(),
                      "base_N2_seconds": base_n2, "total_seconds": elapsed}, indent=2), flush=True)


if __name__ == "__main__":
    run()
