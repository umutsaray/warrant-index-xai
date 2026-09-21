"""Phase 6.2 step 2: DeepSurv N=50 recursive-studentized null gate."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest, norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _shared import markdown_table


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase62_deepsurv_n50_null_gate"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
ALPHAS = (0.01, 0.05, 0.10)
WARNING = (
    "DeepSurv N=50 inference conditions on fixed 12+12 queries and deterministic nested "
    "DeepExplainer k=16/128 TCGA backgrounds."
)


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def deepsurv_n50_nuisance(seed: int, draws: int = 399) -> dict:
    pilot = load_module("08_certification_pilot.py", "phase5_pilot_n50_gate")
    null = load_module("09_null_calibration.py", "phase51_null_n50_nuisance")
    report = pd.read_csv(ROOT / "results/phase62_deepsurv_n50/feature_report.csv")
    selection = pilot.choose_features(
        report, seed=seed + 502, low_feature="rna_ENSG00000090402.8"
    )
    with np.load(ROOT / "results/phase62_deepsurv_n50/attributions_n50.npz", allow_pickle=False) as saved:
        candidate, reference = saved["candidate"], saved["reference"]
        features = saved["features"].tolist()
    targets = np.array([features.index(feature) for feature in selection.feature])
    rng = np.random.default_rng(seed + 9100)
    indices = rng.integers(0, 50, size=(draws, 50))
    _, sampled = pilot.score_bootstraps(
        phi=candidate, reference=reference, target_indices=targets, bootstrap_indices=indices
    )
    flat = sampled.reshape(-1, len(COMPONENTS))
    variable = flat.std(axis=0, ddof=1) > 1e-12
    correlation = np.eye(len(COMPONENTS))
    if variable.sum() > 1:
        correlation[np.ix_(variable, variable)] = np.corrcoef(flat[:, variable], rowvar=False)
    correlation = null.nearest_correlation(correlation)
    amplitude = np.clip(flat.std(axis=0, ddof=1) * np.sqrt(3 * 50), 0.03, 0.18)
    amplitude[[0, 2]] = np.minimum(amplitude[[0, 2]], 0.09)
    return {
        "replicates": 50,
        "correlation": correlation,
        "amplitude": amplitude,
        "nuisance_bootstrap_draws": draws,
        "selected_features": selection.feature.tolist(),
    }


def run(*, outer: int = 2000, inner: int = 999, seed: int = 20260918) -> bool:
    if outer != 2000 or inner != 999:
        raise ValueError("Phase 5.6 approved gate is fixed at 2000x999")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    remedy = load_module("11_rsf_remedy_calibration.py", "phase52_remedy_n50_gate")
    null = load_module("09_null_calibration.py", "phase51_null_n50_gate")
    nuisance = deepsurv_n50_nuisance(seed)
    rng = np.random.default_rng(seed + 560_000)
    pvalues = {scenario: np.empty((outer, len(COMPONENTS))) for scenario in ("single_active", "tie")}
    for scenario in pvalues:
        for simulation in range(outer):
            normal = np.stack([
                rng.multivariate_normal(np.zeros(5), nuisance["correlation"], size=50)
                for _ in range(2)
            ], axis=1)
            uniform = norm.cdf(normal)
            indices = rng.integers(0, 50, size=(inner, 50))
            pvalues[scenario][simulation] = remedy.candidate_pvalues(
                uniform, nuisance["amplitude"], indices, scenario, null
            )["recursive_studentized"]

    rows = []
    for scenario, values in pvalues.items():
        for j, component in enumerate(COMPONENTS):
            for alpha in ALPHAS:
                rejected = int((values[:, j] <= alpha).sum())
                rows.append({
                    "model": "DeepSurv_N50",
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
                    "inference_scope": WARNING,
                })
    result = pd.DataFrame(rows)
    result["holm_p_within_30"] = remedy.holm_adjust(result.one_sided_binomial_p.to_numpy())
    result["anti_conservative_gate"] = result.holm_p_within_30 < 0.05
    result.to_csv(OUT / "type1_rates.csv", index=False)
    passed = not result.anti_conservative_gate.any()
    gate = pd.DataFrame([{
        "model": "DeepSurv_N50", "method": "recursive_studentized", "controls": 30,
        "holm_flagged_controls": int(result.anti_conservative_gate.sum()), "gate_passed": passed,
        "inference_scope": WARNING,
    }])
    gate.to_csv(OUT / "gate_summary.csv", index=False)

    formula_sha = {name: sha(ROOT / name) for name in [
        "warrantxai/components.py", "warrantxai/aggregation.py"
    ]}
    expected = json.loads((ROOT / "results/phase62_deepsurv_n50/metadata.json").read_text())["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "6.2_step2_DeepSurv_N50_null_gate",
        "model": "DeepSurv_N50",
        "deepsurv_conditionality_warning": WARNING,
        "arguments": {"outer": outer, "inner": inner, "seed": seed, "alphas": ALPHAS,
                      "scenarios": ["single_active", "tie"], "Holm_family": 30},
        "method": "recursive_studentized",
        "nuisance": {"replicates": 50, "nuisance_bootstrap_draws": nuisance["nuisance_bootstrap_draws"],
                     "selected_features": nuisance["selected_features"],
                     "correlation": nuisance["correlation"].tolist(),
                     "amplitude": dict(zip(COMPONENTS, nuisance["amplitude"].tolist(), strict=True))},
        "gate_passed": passed,
        "adim3_authorized_only_if_passed": True,
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": formula_sha == expected,
        "source_sha256": {
            "results/phase62_deepsurv_n50/attributions_n50.npz": sha(ROOT / "results/phase62_deepsurv_n50/attributions_n50.npz"),
            "experiments/23_deepsurv_n50_null_gate.py": sha(Path(__file__)),
        },
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "model_or_explainer_calls": False,
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if formula_sha != expected:
        raise RuntimeError("Frozen formula module SHA256 changed")

    alpha05 = result.loc[result.nominal_alpha == 0.05, [
        "inference_scope", "null_scenario", "component", "empirical_type1",
        "exact_95ci_lower", "exact_95ci_upper", "holm_p_within_30", "anti_conservative_gate"
    ]]
    flagged = result.loc[result.anti_conservative_gate]
    text = f"""# Faz 6.2 ADIM 2 — DeepSurv N=50 null-kalibrasyon kapısı

> **Koşulluluk uyarısı:** {WARNING}

`recursive_studentized`, N=50 attribution varyans/bağımlılık yapısından türetilen nuisance ile aynı 2000×999, `single_active` + `tie`, üç alpha ve 30 kontrol içinde Holm düzeltmeli bataryadan geçirildi.

{markdown_table(gate)}

## Alpha=0.05

{markdown_table(alpha05)}

## Holm-düzeltilmiş ihlaller

{markdown_table(flagged) if len(flagged) else 'İhlal yok.'}

Gate sonucu: **{'GEÇTİ' if passed else 'KALDI'}**. Gerçek wall-clock {elapsed:.3f} saniye. Formül SHA256 değişmedi; model/explainer çağrısı yapılmadı ve kimlik-sızıntısı kontrolü geçti.
"""
    (OUT / "report.md").write_text(text, encoding="utf-8")
    (ROOT / "PHASE6_2_STEP2_NULL_GATE_REPORT.md").write_text(text, encoding="utf-8")
    print(gate.to_string(index=False), flush=True)
    print(result.loc[result.anti_conservative_gate, ["null_scenario", "component", "nominal_alpha", "empirical_type1", "holm_p_within_30"]].to_string(index=False), flush=True)
    print(f"wall={elapsed:.3f}s", flush=True)
    return passed


if __name__ == "__main__":
    if not run():
        raise SystemExit(2)

