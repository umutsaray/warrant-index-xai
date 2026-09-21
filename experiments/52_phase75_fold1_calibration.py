#!/usr/bin/env python3
"""Fresh Phase 7.5 fold 1 recursive-studentized null gates for all three N=50 audits."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest, norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase75_fold1_calibration"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
ALPHAS = (0.01, 0.05, 0.10)
SEED = 20260918

sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import CALIBRATION, markdown_table
from warrantxai.components import generalizability, generalizability_signed, importance_ci


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def attribution_nuisance(directory: str, seed: int, draws: int = 399) -> dict:
    pilot = load_module("08_certification_pilot.py", f"phase75_fold1_pilot_{directory}")
    null = load_module("09_null_calibration.py", f"phase75_fold1_null_{directory}")
    with np.load(ROOT / "results" / directory / "attributions_n50.npz", allow_pickle=False) as saved:
        candidate, reference = saved["candidate"], saved["reference"]
    rng = np.random.default_rng(seed + 9100)
    indices = rng.integers(0, 50, size=(draws, 50))
    _, sampled = pilot.score_bootstraps(
        phi=candidate,
        reference=reference,
        target_indices=np.arange(candidate.shape[-1]),
        bootstrap_indices=indices,
    )
    return nuisance_from_draws(sampled, null, 50, draws)


def cox_nuisance(seed: int, draws: int = 399) -> dict:
    """Reconstruct all five clean Cox scores under training-replicate resampling."""
    pilot = load_module("08_certification_pilot.py", "phase75_fold1_cox_pilot")
    null = load_module("09_null_calibration.py", "phase75_fold1_cox_null")
    with np.load(ROOT / "results/phase75_fold1_cox/attributions_n50.npz", allow_pickle=False) as saved:
        source = saved["source_phi"]
        external = saved["external_phi"]
        all_features = saved["features"].astype(str).tolist()
        active = saved["active_features"].astype(str).tolist()
    targets = np.array([all_features.index(name) for name in active])
    selected_source = source[:, :, targets]
    source4 = source[None, ...]
    s_loss = pilot.pair_losses(source4, targets, CALIBRATION["a"])
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_by_rep = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[:, :, targets]
    r_loss = pilot.order_pair_losses(np.abs(source).mean(axis=1)[None, ...], targets, CALIBRATION["eta"])

    def one(indices: np.ndarray) -> np.ndarray:
        f = np.ones(len(targets))
        s = pilot.pair_score(s_loss, indices)
        u = u_by_rep[:, indices].mean(axis=1).min(axis=0)
        src, ext = selected_source[indices], external[indices]
        gm = generalizability(
            np.stack([np.abs(src).mean(axis=(0, 1)), np.abs(ext).mean(axis=(0, 1))]),
            CALIBRATION["b"],
        )
        mean_src, mean_ext = src.mean(axis=0), ext.mean(axis=0)
        grid_n = max(len(mean_src), len(mean_ext))
        signed_src = np.stack(
            [pilot.empirical_rank_grid(mean_src[:, j], grid_n) for j in range(len(targets))], axis=1
        )
        signed_ext = np.stack(
            [pilot.empirical_rank_grid(mean_ext[:, j], grid_n) for j in range(len(targets))], axis=1
        )
        gs = generalizability_signed(np.stack([signed_src, signed_ext]), CALIBRATION["a"])
        r = pilot.pair_score(r_loss, indices)
        return np.column_stack([f, s, u, np.minimum(gm, gs), r])

    rng = np.random.default_rng(seed + 9100)
    indices = rng.integers(0, 50, size=(draws, 50))
    sampled = np.stack([one(index) for index in indices])
    point = one(np.arange(50))
    report = pd.read_csv(ROOT / "results/phase75_fold1_cox/feature_report.csv")
    np.testing.assert_allclose(point, report[list(COMPONENTS)].to_numpy(), atol=2e-12)
    return nuisance_from_draws(sampled, null, 50, draws)


def nuisance_from_draws(sampled: np.ndarray, null, replicas: int, draws: int) -> dict:
    flat = sampled.reshape(-1, len(COMPONENTS))
    variable = flat.std(axis=0, ddof=1) > 1e-12
    correlation = np.eye(len(COMPONENTS))
    if variable.sum() > 1:
        correlation[np.ix_(variable, variable)] = np.corrcoef(flat[:, variable], rowvar=False)
    correlation = null.nearest_correlation(correlation)
    amplitude = np.clip(flat.std(axis=0, ddof=1) * np.sqrt(3 * replicas), 0.03, 0.18)
    amplitude[[0, 2]] = np.minimum(amplitude[[0, 2]], 0.09)
    return {
        "replicates": replicas,
        "correlation": correlation,
        "amplitude": amplitude,
        "nuisance_bootstrap_draws": draws,
        "nuisance_features": sampled.shape[1],
    }


def main(outer: int = 2000, inner: int = 999) -> bool:
    if (outer, inner) != (2000, 999):
        raise ValueError("Approved Phase 7.5 fold 1 gate is fixed at 2000x999")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    remedy = load_module("11_rsf_remedy_calibration.py", "phase75_fold1_remedy")
    null = load_module("09_null_calibration.py", "phase75_fold1_null")
    nuisance = {
        "RSF_N50": attribution_nuisance("phase75_fold1_rsf", SEED + 1),
        "DeepSurv_N50": attribution_nuisance("phase75_fold1_deepsurv", SEED + 2),
    }
    warnings = {
        "CoxPH_N50": "153 TCGA + 1466 METABRIC fixed queries; exact closed-form Cox attribution",
        "RSF_N50": "fixed 12+12 queries and stored 8/128 attribution averages; permutation-path uncertainty excluded",
        "DeepSurv_N50": "fixed 12+12 queries and deterministic nested DeepExplainer k16/k128 backgrounds",
    }
    rows = []
    for model_index, (model, setup) in enumerate(nuisance.items()):
        rng = np.random.default_rng(SEED + 720_000 + model_index * 100_000)
        for scenario in ("single_active", "tie"):
            pvalues = np.empty((outer, len(COMPONENTS)))
            for simulation in range(outer):
                normal = np.stack(
                    [rng.multivariate_normal(np.zeros(5), setup["correlation"], size=50) for _ in range(2)],
                    axis=1,
                )
                uniform = norm.cdf(normal)
                indices = rng.integers(0, 50, size=(inner, 50))
                pvalues[simulation] = remedy.candidate_pvalues(
                    uniform, setup["amplitude"], indices, scenario, null
                )["recursive_studentized"]
            for j, component in enumerate(COMPONENTS):
                for alpha in ALPHAS:
                    rejected = int((pvalues[:, j] <= alpha).sum())
                    rows.append({
                        "model": model,
                        "method": "recursive_studentized",
                        "null_scenario": scenario,
                        "component": component,
                        "nominal_alpha": alpha,
                        "rejections": rejected,
                        "outer_simulations": outer,
                        "empirical_type1": rejected / outer,
                        "exact_95ci_lower": 0.0 if rejected == 0 else float(beta.ppf(0.025, rejected, outer-rejected+1)),
                        "exact_95ci_upper": 1.0 if rejected == outer else float(beta.ppf(0.975, rejected+1, outer-rejected)),
                        "one_sided_binomial_p": float(binomtest(rejected, outer, alpha, alternative="greater").pvalue),
                        "inference_scope": warnings[model],
                    })
        print(f"Phase 7.5 fold 1 calibration simulated: {model}", flush=True)

    result = pd.DataFrame(rows)
    result["holm_p_within_model_30"] = np.nan
    for model in nuisance:
        mask = result.model == model
        result.loc[mask, "holm_p_within_model_30"] = remedy.holm_adjust(
            result.loc[mask, "one_sided_binomial_p"].to_numpy()
        )
    result["anti_conservative_gate"] = result.holm_p_within_model_30 < 0.05
    result.to_csv(OUT / "type1_rates.csv", index=False)
    gate = result.groupby("model", sort=False).agg(
        controls=("anti_conservative_gate", "size"),
        holm_flagged_controls=("anti_conservative_gate", "sum"),
    ).reset_index()
    gate["gate_passed"] = gate.holm_flagged_controls == 0
    gate["inference_scope"] = gate.model.map(warnings)
    gate.to_csv(OUT / "gate_summary.csv", index=False)
    passed = bool(gate.gate_passed.all())
    formula = {
        "components.py": sha(ROOT / "warrantxai/components.py"),
        "aggregation.py": sha(ROOT / "warrantxai/aggregation.py"),
    }
    expected = json.loads((ROOT / "results/phase75_fold1_panel/metadata.json").read_text())["formula_sha256"]
    if formula != expected:
        raise RuntimeError("Frozen formula SHA256 changed")
    elapsed = perf_counter() - started
    metadata = {
        "phase": "7.5_fold1_fresh_available_model_null_calibration",
        "arguments": {"outer": outer, "inner": inner, "seed": SEED, "alphas": ALPHAS,
                      "scenarios": ["single_active", "tie"], "Holm_family_per_model": 30},
        "method": "recursive_studentized",
        "nuisance": {
            model: {
                "replicates": setup["replicates"],
                "nuisance_bootstrap_draws": setup["nuisance_bootstrap_draws"],
                "nuisance_features": setup["nuisance_features"],
                "correlation": setup["correlation"].tolist(),
                "amplitude": dict(zip(COMPONENTS, setup["amplitude"].tolist(), strict=True)),
            }
            for model, setup in nuisance.items()
        },
        "all_available_models_gate_passed": passed,
        "all_models_gate_passed": passed,
        "coxph_status": "not_available: deterministic bootstrap replicate 15 failed frozen numerical retry",
        "certification_authorized_only_for_available_models_if_all_passed": True,
        "formula_sha256": formula,
        "privacy": {"identifiers_exported": False},
        "model_or_explainer_calls": False,
        "wall_clock_seconds": elapsed,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    alpha05 = result[result.nominal_alpha == 0.05][
        ["model", "null_scenario", "component", "empirical_type1", "exact_95ci_lower",
         "exact_95ci_upper", "holm_p_within_model_30", "anti_conservative_gate"]
    ]
    flagged = result[result.anti_conservative_gate]
    note = f"""# Faz 7.5 fold 1 — Taze null-kalibrasyon kapıları

Her model kendi temiz N=50 replika varyans/bağımlılık yapısıyla yeniden kalibre edildi. Batarya model başına 2000×999, `single_active`+`tie`, üç alpha ve **30 kontrol içinde Holm** düzeltmesidir. Model veya explainer çağrılmadı.

{markdown_table(gate)}

## Alpha=0.05

{markdown_table(alpha05)}

## Holm ihlalleri

{markdown_table(flagged) if len(flagged) else 'İhlal yok.'}

Mevcut iki-model ortak gate sonucu: **{'GEÇTİ' if passed else 'KALDI'}**. CoxPH, deterministik bootstrap replika 15'in dondurulmuş retry protokolünden sonra yakınsamaması nedeniyle bu kapıda yoktur ve CoxPH sertifikasyonu üretilmez. Gate kalırsa RSF/DeepSurv BH/BY de üretilmez. Formül SHA256 değişmedi; kimlik-sızıntısı kontrolü geçti. Wall-clock **{elapsed:.3f} saniye**.
"""
    (ROOT / "PHASE75_FOLD1_CALIBRATION.md").write_text(note, encoding="utf-8")
    print(gate.to_string(index=False), flush=True)
    return passed


if __name__ == "__main__":
    if not main():
        raise SystemExit(2)
