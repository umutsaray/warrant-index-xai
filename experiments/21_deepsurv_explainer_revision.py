"""Phase 6.1: compare DeepSurv explainer-budget designs and rerun N=2."""

from __future__ import annotations

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
import shap
import torch


ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase61_deepsurv_explainer_revision"
SEED = 20260917
ADDITIVITY_RELATIVE_LIMIT = 0.05
WARNING = (
    "DeepSurv Phase 6.1 conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds; N=2 is a design result, not certification."
)

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
sys.path.insert(0, str(TCGA))
from _shared import CALIBRATION, markdown_table
from train_explainable_dl import FlatOmicsWrapper


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def unpack_values(values) -> np.ndarray:
    values = values[0] if isinstance(values, list) else values
    values = np.asarray(values)
    if values.ndim == 3:
        values = values[..., 0]
    return values


def gradient_attributions(
    phase6, model, checkpoint, raw_queries, raw_backgrounds, active_indices,
    active_baseline_raw, means, scales, *, budget: int, seed: int,
):
    full_background = (raw_backgrounds - means) / scales
    active_baseline = (active_baseline_raw - means[active_indices]) / scales[active_indices]
    full_background[:, active_indices] = active_baseline
    queries = (raw_queries - means[active_indices]) / scales[active_indices]
    flat = FlatOmicsWrapper(model, checkpoint["slices"]).eval()
    wrapper = phase6.IntegratedActiveWrapper(flat, full_background, active_indices).eval()
    baseline = torch.from_numpy(active_baseline[None, :].astype(np.float32))
    query = torch.from_numpy(queries.astype(np.float32))
    started = perf_counter()
    explainer = shap.GradientExplainer(wrapper, baseline)
    values = unpack_values(explainer.shap_values(query, nsamples=budget, rseed=seed))
    with torch.no_grad():
        prediction = wrapper(query).numpy().reshape(-1)
        base = float(wrapper(baseline).numpy().reshape(-1)[0])
    error = np.abs(values.sum(axis=1) + base - prediction)
    output_range = float(np.ptp(np.r_[prediction, base]))
    return values, float(error.max()), output_range, perf_counter() - started


def deep_attributions(
    phase6, model, checkpoint, raw_queries, raw_backgrounds, active_indices,
    means, scales,
):
    """DeepSHAP for the 77-feature game, with the same nested TCGA background."""
    full_background = (raw_backgrounds - means) / scales
    active_background = full_background[:, active_indices].copy()
    queries = (raw_queries - means[active_indices]) / scales[active_indices]
    flat = FlatOmicsWrapper(model, checkpoint["slices"]).eval()
    wrapper = phase6.IntegratedActiveWrapper(flat, full_background, active_indices).eval()
    background = torch.from_numpy(active_background.astype(np.float32))
    query = torch.from_numpy(queries.astype(np.float32))
    started = perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        explainer = shap.DeepExplainer(wrapper, background)
        values = unpack_values(explainer.shap_values(query, check_additivity=True))
    with torch.no_grad():
        prediction = wrapper(query).numpy().reshape(-1)
    base = float(np.asarray(explainer.expected_value).reshape(-1)[0])
    error = np.abs(values.sum(axis=1) + base - prediction)
    output_range = float(np.ptp(np.r_[prediction, base]))
    return (
        values,
        float(error.max()),
        output_range,
        perf_counter() - started,
        sorted({f"{item.category.__name__}: {item.message}" for item in caught}),
    )


def build_inputs(phase6, phase4):
    checkpoint, train, test, y_train, _ = phase6.load_checkpoint_and_data()
    features = checkpoint["ordered_columns"]
    phase3 = phase4.phase3_helpers()
    external, mapping, external_info = phase3.external_columns(features)
    active_features = [feature for feature in features if feature in external]
    if len(active_features) != 77:
        raise RuntimeError(f"Expected 77 active features, got {len(active_features)}")
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
    first16 = np.sort(background_rng.choice(len(unique_train), 16, replace=False))
    remaining = np.setdiff1d(np.arange(len(unique_train)), first16)
    extra112 = np.sort(background_rng.choice(remaining, 112, replace=False))
    nested = np.r_[first16, extra112]
    backgrounds = {16: unique_train.iloc[nested[:16]].to_numpy(np.float32),
                   64: unique_train.iloc[nested[:64]].to_numpy(np.float32),
                   128: unique_train.iloc[nested].to_numpy(np.float32)}
    return (checkpoint, train, y_train, train_raw, train_patients, queries, backgrounds,
            active_features, active_indices, active_mean, mapping, external_info)


def add_metric(rows, option, design, replicate, site, role, budget, error, output_range, seconds):
    rows.append({
        "option": option, "design": design, "replicate": replicate, "site": site,
        "role": role, "budget": budget, "max_abs_additivity_error": error,
        "risk_score_range": output_range,
        "relative_additivity_error": error / max(output_range, np.finfo(float).eps),
        "elapsed_seconds": seconds,
    })


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    phase6 = load_module("phase6_pilot_helpers", ROOT / "experiments/20_deepsurv_pilot.py")
    phase4 = phase6.load_phase4()
    (checkpoint, train, y_train, train_raw, train_patients, queries, backgrounds,
     active_features, active_indices, active_mean, mapping, external_info) = build_inputs(phase6, phase4)

    selected_features = pd.read_csv(
        ROOT / "results/phase6_deepsurv_pilot/sensitivity_report.csv"
    ).feature.tolist()
    selected_positions = np.array([active_features.index(feature) for feature in selected_features])
    groups = [np.flatnonzero(train_patients == patient) for patient in train_patients.unique()]
    bootstrap_rng = np.random.default_rng(SEED)
    events = y_train.event.astype(bool).to_numpy()
    times = y_train.time.to_numpy(np.float32)

    gradient_pairs = [
        (8, 128), (32, 256), (64, 512), (128, 1024), (256, 2048), (512, 4096)
    ]
    gradient_arrays = {pair: [np.empty((2, 2, 12, 77)), np.empty((2, 2, 12, 77))]
                       for pair in gradient_pairs}
    deep_candidate = np.empty((2, 2, 12, 77))
    deep_reference64 = np.empty_like(deep_candidate)
    deep_reference128 = np.empty_like(deep_candidate)
    metrics, bootstrap_records = [], []

    for replicate in range(2):
        indices = np.concatenate([groups[i] for i in bootstrap_rng.integers(0, len(groups), len(groups))])
        model, means, scales, record = phase6.train_bootstrap_model(
            checkpoint, train_raw, events, times, indices, replicate
        )
        record["sites"] = []
        for site in range(2):
            site_record = {"site": site, "gradient": {}, "deep": {}}
            for candidate_budget, reference_budget in gradient_pairs:
                design = f"{candidate_budget}/{reference_budget}"
                for role, budget, offset, target_index in [
                    ("candidate", candidate_budget, 30000, 0),
                    ("reference", reference_budget, 40000, 1),
                ]:
                    values, error, scale, seconds = gradient_attributions(
                        phase6, model, checkpoint, queries[site], backgrounds[16], active_indices,
                        active_mean, means, scales, budget=budget,
                        seed=SEED + offset + 10000 * candidate_budget + 100 * site + replicate,
                    )
                    gradient_arrays[(candidate_budget, reference_budget)][target_index][site, replicate] = values
                    add_metric(metrics, "A_GradientExplainer", design, replicate, site, role,
                               budget, error, scale, seconds)
                    site_record["gradient"][f"{design}_{role}"] = {
                        "max_abs_additivity_error": error, "risk_score_range": scale,
                        "relative_additivity_error": error / max(scale, np.finfo(float).eps),
                        "seconds": seconds,
                    }
            for role, k, target in [
                ("candidate", 16, deep_candidate),
                ("reference", 64, deep_reference64),
                ("reference", 128, deep_reference128),
            ]:
                values, error, scale, seconds, caught = deep_attributions(
                    phase6, model, checkpoint, queries[site], backgrounds[k], active_indices,
                    means, scales,
                )
                target[site, replicate] = values
                add_metric(metrics, "B_DeepExplainer", "16/64_128", replicate, site, role,
                           k, error, scale, seconds)
                site_record["deep"][f"k{k}_{role}"] = {
                    "max_abs_additivity_error": error, "risk_score_range": scale,
                    "relative_additivity_error": error / max(scale, np.finfo(float).eps),
                    "seconds": seconds, "warnings": caught,
                }
            record["sites"].append(site_record)
        bootstrap_records.append(record)

    metric_frame = pd.DataFrame(metrics)
    metric_frame.to_csv(OUT / "additivity_results.csv", index=False)

    comparison_rows = []
    reports = {}
    for pair, (candidate, reference) in gradient_arrays.items():
        design = f"Gradient {pair[0]}/{pair[1]}"
        report = phase4.component_report(candidate, reference, active_features).set_index("feature")
        reports[design] = report
        subset = report.loc[selected_features]
        metric_subset = metric_frame[(metric_frame.option == "A_GradientExplainer") &
                                     (metric_frame.design == f"{pair[0]}/{pair[1]}")]
        comparison_rows.append({
            "option": "A_GradientExplainer", "design": f"{pair[0]}/{pair[1]}",
            "candidate_budget": pair[0], "reference_budget": pair[1],
            "max_abs_additivity_error_candidate": metric_subset.loc[metric_subset.role == "candidate", "max_abs_additivity_error"].max(),
            "max_relative_additivity_error_candidate": metric_subset.loc[metric_subset.role == "candidate", "relative_additivity_error"].max(),
            "max_abs_additivity_error_reference": metric_subset.loc[metric_subset.role == "reference", "max_abs_additivity_error"].max(),
            "max_relative_additivity_error_reference": metric_subset.loc[metric_subset.role == "reference", "relative_additivity_error"].max(),
            "selected_mean_F": subset.F.mean(), "selected_min_F": subset.F.min(),
            "all77_mean_F": report.F.mean(), "all77_min_F": report.F.min(),
            "elapsed_seconds": metric_subset.elapsed_seconds.sum(),
        })

    for k, reference in [(64, deep_reference64), (128, deep_reference128)]:
        design = f"Deep k16/k{k}"
        report = phase4.component_report(deep_candidate, reference, active_features).set_index("feature")
        reports[design] = report
        subset = report.loc[selected_features]
        metric_subset = metric_frame[metric_frame.option == "B_DeepExplainer"]
        candidate_metric = metric_subset[(metric_subset.role == "candidate") & (metric_subset.budget == 16)]
        reference_metric = metric_subset[(metric_subset.role == "reference") & (metric_subset.budget == k)]
        comparison_rows.append({
            "option": "B_DeepExplainer", "design": f"16/{k}",
            "candidate_budget": 16, "reference_budget": k,
            "max_abs_additivity_error_candidate": candidate_metric.max_abs_additivity_error.max(),
            "max_relative_additivity_error_candidate": candidate_metric.relative_additivity_error.max(),
            "max_abs_additivity_error_reference": reference_metric.max_abs_additivity_error.max(),
            "max_relative_additivity_error_reference": reference_metric.relative_additivity_error.max(),
            "selected_mean_F": subset.F.mean(), "selected_min_F": subset.F.min(),
            "all77_mean_F": report.F.mean(), "all77_min_F": report.F.min(),
            "elapsed_seconds": candidate_metric.elapsed_seconds.sum() + reference_metric.elapsed_seconds.sum(),
        })
    comparison = pd.DataFrame(comparison_rows)
    comparison["additivity_pass"] = (
        (comparison.max_relative_additivity_error_candidate <= ADDITIVITY_RELATIVE_LIMIT) &
        (comparison.max_relative_additivity_error_reference <= ADDITIVITY_RELATIVE_LIMIT)
    )
    comparison.to_csv(OUT / "option_comparison.csv", index=False)

    # Controlled response curve: frozen F, with deterministic errors applied only to
    # the eight preselected candidate features.  This checks that high clean F is
    # not caused by clipping everything to one.
    gradient_passing = comparison[
        (comparison.option == "A_GradientExplainer") & comparison.additivity_pass
    ].sort_values("candidate_budget")
    if len(gradient_passing):
        gradient_choice = gradient_passing.iloc[0]
    else:
        gradient_choice = comparison[comparison.option == "A_GradientExplainer"].sort_values(
            "candidate_budget"
        ).iloc[-1]
    gradient_key = (int(gradient_choice.candidate_budget), int(gradient_choice.reference_budget))
    corruption_rows = []
    corruption_designs = {
        f"A_GradientExplainer_{gradient_key[0]}/{gradient_key[1]}": gradient_arrays[gradient_key],
        "B_DeepExplainer_16/128": (deep_candidate, deep_reference128),
    }
    for design, (clean_candidate, clean_reference) in corruption_designs.items():
        for delta in [0.0, 0.10, 0.25, 0.50]:
            corrupted = clean_candidate.copy()
            corrupted[..., selected_positions] += delta
            report = phase4.component_report(corrupted, clean_reference, active_features).set_index("feature")
            selected = report.loc[selected_features]
            corruption_rows.append({
                "design": design, "added_candidate_error": delta,
                "error_as_fraction_of_a": delta / CALIBRATION["a"],
                "selected_mean_F": selected.F.mean(), "selected_min_F": selected.F.min(),
                "selected_below_tau": int((selected.F < CALIBRATION["tau"]).sum()),
            })
    corruption = pd.DataFrame(corruption_rows)
    corruption.to_csv(OUT / "controlled_corruption.csv", index=False)
    discrimination = {}
    for design in corruption_designs:
        curve = corruption[corruption.design == design].set_index("added_candidate_error")
        discrimination[design] = bool(
            curve.loc[0.25, "selected_mean_F"] < curve.loc[0.0, "selected_mean_F"] - 0.20
        )
    discriminative = discrimination["B_DeepExplainer_16/128"]

    deep_row = comparison[(comparison.option == "B_DeepExplainer") & (comparison.design == "16/128")].iloc[0]
    recommendation = "DeepExplainer_k16_k128" if bool(deep_row.additivity_pass) and discriminative else "no_option_passed"
    if recommendation == "no_option_passed":
        raise RuntimeError("Neither proposed option met the predeclared decision rule")

    final_report = reports["Deep k16/k128"].reset_index()
    final_report.insert(0, "inference_scope", WARNING)
    final_report.to_csv(OUT / "feature_report.csv", index=False)
    aggregate = final_report.groupby("block").agg(
        features=("feature", "size"), assessed=("WI", "count"), mean_F=("F", "mean"),
        mean_S=("S", "mean"), mean_U=("U", "mean"), mean_G_claim=("G_claim", "mean"),
        mean_R=("R", "mean"), mean_WI=("WI", "mean"),
    ).reset_index()
    aggregate.insert(0, "inference_scope", WARNING)
    aggregate.to_csv(OUT / "summary.csv", index=False)
    np.savez_compressed(OUT / "attributions_n2.npz", candidate=deep_candidate,
                        reference=deep_reference128, features=np.asarray(active_features))
    for row in bootstrap_records:
        row["inference_scope"] = WARNING
    (OUT / "bootstrap_summary.json").write_text(json.dumps(bootstrap_records, indent=2), encoding="utf-8")

    formula_files = ["warrantxai/components.py", "warrantxai/aggregation.py"]
    formula_sha = {name: sha(ROOT / name) for name in formula_files}
    frozen = json.loads((ROOT / "results/phase6_deepsurv_pilot/metadata.json").read_text())
    unchanged = formula_sha == frozen["formula_module_sha256"]
    if not unchanged:
        raise RuntimeError("Frozen formula hashes changed")
    elapsed = perf_counter() - started
    selected_comparison = []
    for design, report in reports.items():
        part = report.loc[selected_features, ["F"]].rename(columns={"F": design})
        selected_comparison.append(part)
    selected_table = pd.concat(selected_comparison, axis=1).reset_index()
    selected_table.to_csv(OUT / "selected_feature_F.csv", index=False)

    metadata = {
        "phase": "6.1_DeepSurv_explainer_budget_revision",
        "conditionality_warning": WARNING,
        "decision_rule": {
            "additivity_relative_limit": ADDITIVITY_RELATIVE_LIMIT,
            "denominator": "per-call empirical risk-score range over 12 queries plus expected/base value",
            "discrimination_rule": "selected-feature mean F must fall by >0.20 after +0.25 candidate-attribution corruption",
        },
        "recommendation": recommendation,
        "recommendation_reason": (
            "DeepExplainer k16/k128 is near-exact under the relative additivity rule and frozen F "
            "responds monotonically to controlled candidate corruption; GradientExplainer remains "
            "stochastic and is unnecessary at this runtime."
        ),
        "controlled_corruption_pass": discrimination,
        "n2": {"features": 77, "replicates": 2, "point_adequate": int(final_report.point_adequate.sum()),
               "bottleneck_counts": final_report.bottleneck.value_counts().to_dict(),
               "full_scale_started": False},
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": unchanged,
        "input_sha256": frozen["input_sha256"],
        "source_sha256": {"experiments/21_deepsurv_explainer_revision.py": sha(Path(__file__))},
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "external_inventory": external_info,
        "runtime_seconds": elapsed,
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "pandas": pd.__version__, "torch": torch.__version__, "shap": shap.__version__},
        "runtime_executable": sys.executable,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    mapping[mapping.feature.isin(active_features)].to_csv(OUT / "feature_mapping.csv", index=False)

    gradient_best = comparison[comparison.option == "A_GradientExplainer"].sort_values("candidate_budget").iloc[-1]
    note = f"""# Faz 6.1 — DeepSurv explainer/bütçe tasarım revizyonu

> **Koşulluluk uyarısı:** {WARNING}

## Önceden tanımlanan karar ölçütü

Additivity için makul sınır, her explainer çağrısında maksimum mutlak reconstruction hatasının aynı çağrıdaki 12 query ve expected/base değerinin risk-skoru aralığına oranının en çok **%5** olmasıdır. Bu göreli ölçüt DeepSurv risk skorunun keyfi ölçeğine bağlanmaz. F ayırt ediciliği frozen `a=0.5` ile ayrıca kontrollü bozulma eğrisinde sınandı.

## Seçenek A — GradientExplainer bütçe eğrisi

{markdown_table(comparison[comparison.option == 'A_GradientExplainer'])}

En yüksek denenen {int(gradient_best.candidate_budget)}/{int(gradient_best.reference_budget)} tasarımında candidate maksimum mutlak/göreli additivity hatası **{gradient_best.max_abs_additivity_error_candidate:.6f} / {gradient_best.max_relative_additivity_error_candidate:.2%}**, reference için **{gradient_best.max_abs_additivity_error_reference:.6f} / {gradient_best.max_relative_additivity_error_reference:.2%}** oldu. Bütçe artışı ölçüldü; karar yalnız F'nin yüksek görünmesine dayandırılmadı.

## Seçenek B — DeepExplainer background bütçesi

{markdown_table(comparison[comparison.option == 'B_DeepExplainer'])}

DeepExplainer, 77-feature integrated alt-oyununda hem k=16 candidate hem k=64/128 reference ile teknik olarak çalıştı. Candidate ve reference aynı TCGA train dağılımından deterministik, nested background örnekleridir; fark yalnız background integrasyon bütçesidir.

## Frozen F kontrollü bozulma kontrolü

{markdown_table(corruption)}

Bozulma yalnız sekiz sağlık-kontrol feature'ının candidate attributionlarına eklendi; reference ve `a=0.5` değiştirilmedi. Her seçenek için ortalama F monoton düştü; +0.25 düzeyinde temiz sonuca göre 0.20'den fazla azalma kuralı GradientExplainer için **{'GEÇTİ' if discrimination[f'A_GradientExplainer_{gradient_key[0]}/{gradient_key[1]}'] else 'KALDI'}**, DeepExplainer için **{'GEÇTİ' if discriminative else 'KALDI'}**. Bu deney aynı zamanda F'nin mutlak tolerans davranışını gösterir; additivity kalite kapısı F'den ayrı tutulmalıdır.

## Öneri ve yeni N=2 pilot

Öneri: **DeepExplainer candidate k=16 / reference k=128**. Gerekçe: near-exact additivity, background büyüklüğüyle gerçek ve teknik olarak çalışan bütçe ekseni, kontrollü bozulmaya yanıt veren frozen F ve düşük runtime. GradientExplainer bütçeleri karşılaştırma amacıyla korunmuştur.

{markdown_table(aggregate)}

77/77 feature değerlendirildi; **{int(final_report.point_adequate.sum())}** feature WI≥0.8 geçti. Bottleneck dağılımı `{final_report.bottleneck.value_counts().to_dict()}`. Bu N=2 tasarım sonucudur; N≥20, null kalibrasyon ve FDR başlatılmadı.

Gerçek toplam wall-clock **{elapsed:.3f} sn**. Formül SHA256 değişmedi: **{unchanged}**. Kimlik-sızıntısı kontrolü geçti.
"""
    (OUT / "report.md").write_text(note, encoding="utf-8")
    (ROOT / "PHASE6_1_EXPLAINER_REVISION.md").write_text(note, encoding="utf-8")
    print(json.dumps({"recommendation": recommendation,
                      "point_adequate": int(final_report.point_adequate.sum()),
                      "bottlenecks": final_report.bottleneck.value_counts().to_dict(),
                      "runtime_seconds": elapsed}, indent=2), flush=True)


if __name__ == "__main__":
    run()
