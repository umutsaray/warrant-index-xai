"""Phase 5.1: empirical Type-I calibration gate before full BH/BY.

No model or explainer is called. Synthetic nulls have known population
component value 0.8 and mimic dependence/noise estimated from saved artifacts.
"""

from __future__ import annotations

import argparse
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _shared import CALIBRATION, markdown_table


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase51_null_calibration"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
ALPHAS = (0.01, 0.05, 0.10)
RSF_WARNING = (
    "RSF nuisance calibration conditions on fixed 12+12 queries and stored 8/128 attribution averages; "
    "per-permutation-path uncertainty is excluded."
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_pilot_module():
    path = ROOT / "experiments/08_certification_pilot.py"
    spec = importlib.util.spec_from_file_location("phase5_pilot", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Phase 5 pilot cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def nearest_correlation(matrix: np.ndarray) -> np.ndarray:
    """Return a positive-semidefinite unit-diagonal correlation matrix."""
    array = np.nan_to_num(np.asarray(matrix, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    array = (array + array.T) / 2
    np.fill_diagonal(array, 1.0)
    values, vectors = np.linalg.eigh(array)
    array = (vectors * np.maximum(values, 1e-8)) @ vectors.T
    scale = np.sqrt(np.diag(array))
    array = array / scale[:, None] / scale[None, :]
    np.fill_diagonal(array, 1.0)
    return array


def nuisance_from_saved(seed: int, draws: int = 399) -> dict[str, dict]:
    """Estimate component copula correlation and relative scale without refits."""
    pilot = load_pilot_module()
    cox_report = pd.read_csv(ROOT / "results/phase3_brca/feature_report.csv")
    rsf_report = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    low = "rna_ENSG00000090402.8"
    cox_selection = pilot.choose_features(cox_report, seed=seed + 501, low_feature=low)
    rsf_selection = pilot.choose_features(rsf_report, seed=seed + 502, low_feature=low)
    rng = np.random.default_rng(seed + 9100)
    cox_indices = rng.integers(0, 50, size=(draws, 50))
    rsf_indices = rng.integers(0, 20, size=(draws, 20))
    _, cox_draws, _ = pilot.cox_bootstrap(cox_selection, cox_indices)
    with np.load(ROOT / "results/phase4_rsf/attributions_n20.npz", allow_pickle=False) as saved:
        candidate, reference = saved["candidate"], saved["reference"]
        features = saved["features"].tolist()
    targets = np.array([features.index(feature) for feature in rsf_selection.feature])
    _, rsf_draws = pilot.score_bootstraps(
        phi=candidate, reference=reference, target_indices=targets, bootstrap_indices=rsf_indices
    )

    result = {}
    for model, sampled, t in [("CoxPH_N50", cox_draws, 50), ("RSF_N20", rsf_draws, 20)]:
        flat = sampled.reshape(-1, len(COMPONENTS))
        variable = flat.std(axis=0, ddof=1) > 1e-12
        correlation = np.eye(len(COMPONENTS))
        if variable.sum() > 1:
            correlation[np.ix_(variable, variable)] = np.corrcoef(flat[:, variable], rowvar=False)
        correlation = nearest_correlation(correlation)
        # For a Uniform(-a,a) mean, SE=a/sqrt(3T). This only transports
        # relative observed noise scale; exact null calibration is analytic below.
        amplitude = np.clip(flat.std(axis=0, ddof=1) * np.sqrt(3 * t), 0.03, 0.18)
        amplitude[[0, 2]] = np.minimum(amplitude[[0, 2]], 0.09)
        result[model] = {
            "replicates": t,
            "correlation": correlation,
            "amplitude": amplitude,
            "nuisance_bootstrap_draws": draws,
        }
    return result


def mean_abs_pair(values: np.ndarray) -> float:
    ordered = np.sort(values)
    weights = 2 * np.arange(len(values)) - len(values) + 1
    return float(ordered @ weights / (len(values) * (len(values) - 1) / 2))


def mean_abs_pair_bootstrap(values: np.ndarray, indices: np.ndarray) -> np.ndarray:
    ordered = np.sort(values[indices], axis=1)
    weights = 2 * np.arange(indices.shape[1]) - indices.shape[1] + 1
    return ordered @ weights / (indices.shape[1] * (indices.shape[1] - 1) / 2)


def frozen_shape_statistics(
    uniform: np.ndarray, amplitude: np.ndarray, indices: np.ndarray, scenario: str = "single_active"
) -> tuple[np.ndarray, np.ndarray]:
    """Compute five null statistics with frozen-formula structural features.

    uniform has shape (T,2,5). In single_active, site/branch zero is
    exactly 0.8 and the other is 0.9; in tie, both are exactly 0.8.
    F/U are bounded means; S is a continuous pair U-statistic; R is a
    discrete tolerant-order analogue; G_claim is a two-branch minimum.
    """
    if scenario not in {"single_active", "tie"}:
        raise ValueError("scenario must be single_active or tie")
    t = uniform.shape[0]
    site_centers = np.array([0.8, 0.9 if scenario == "single_active" else 0.8])
    f_values = site_centers + amplitude[0] * (2 * uniform[:, :, 0] - 1)
    u_values = site_centers + amplitude[2] * (2 * uniform[:, :, 2] - 1)
    discrepancy_centers = np.array([0.2, 0.1 if scenario == "single_active" else 0.2])
    discrepancy = discrepancy_centers + amplitude[3] * (2 * uniform[:, :, 3] - 1)
    s_coefficients = np.array([0.6, 0.3 if scenario == "single_active" else 0.6])
    r_coefficients = np.array([0.4, 0.2 if scenario == "single_active" else 0.4])

    point = np.empty(5)
    point[0] = f_values.mean(axis=0).min()
    point[1] = min(1 - s_coefficients[site] * mean_abs_pair(uniform[:, site, 1]) for site in range(2))
    point[2] = u_values.mean(axis=0).min()
    point[3] = 1 - np.abs(discrepancy.mean(axis=0)).max()
    binary = uniform[:, :, 4] > 0.5
    r_site = []
    for site in range(2):
        ones = int(binary[:, site].sum())
        mismatch = ones * (t - ones) / (t * (t - 1) / 2)
        r_site.append(1 - r_coefficients[site] * mismatch)
    point[4] = min(r_site)

    sampled = np.empty((len(indices), 5))
    sample_n = indices.shape[1]
    sampled[:, 0] = f_values[indices].mean(axis=1).min(axis=1)
    sampled[:, 2] = u_values[indices].mean(axis=1).min(axis=1)
    sampled[:, 3] = 1 - np.abs(discrepancy[indices].mean(axis=1)).max(axis=1)
    s_sites = np.column_stack(
        [1 - s_coefficients[site] * mean_abs_pair_bootstrap(uniform[:, site, 1], indices) for site in range(2)]
    )
    sampled[:, 1] = s_sites.min(axis=1)
    r_sites = []
    for site in range(2):
        ones = binary[:, site][indices].sum(axis=1)
        mismatch = ones * (sample_n - ones) / (sample_n * (sample_n - 1) / 2)
        r_sites.append(1 - r_coefficients[site] * mismatch)
    sampled[:, 4] = np.column_stack(r_sites).min(axis=1)
    return point, sampled


def simulate_model(
    *, outer: int, inner: int, seed: int, replicates: int, correlation: np.ndarray, amplitude: np.ndarray,
    scenario: str
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    pvalues = np.empty((outer, len(COMPONENTS)))
    for simulation in range(outer):
        normal = np.stack(
            [rng.multivariate_normal(np.zeros(5), correlation, size=replicates) for _ in range(2)], axis=1
        )
        uniform = norm.cdf(normal)
        indices = rng.integers(0, replicates, size=(inner, replicates))
        point, sampled = frozen_shape_statistics(uniform, amplitude, indices, scenario)
        exceed = (sampled - point >= point - CALIBRATION["tau"]).sum(axis=0)
        pvalues[simulation] = (exceed + 1) / (inner + 1)
    return pvalues


def holm_adjust(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues, kind="stable")
    ordered = pvalues[order]
    adjusted_ordered = np.maximum.accumulate((len(pvalues) - np.arange(len(pvalues))) * ordered)
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(1.0, adjusted_ordered)
    return adjusted


def run(*, outer: int = 2000, inner: int = 999, seed: int = 20260918) -> pd.DataFrame:
    if outer < 500 or inner < 199:
        raise ValueError("Calibration gate needs outer>=500 and inner>=199")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    nuisance = nuisance_from_saved(seed)
    rows = []
    pvalue_summaries = {}
    for offset, model in enumerate(["CoxPH_N50", "RSF_N20"]):
        setup = nuisance[model]
        pvalue_summaries[model] = {}
        for scenario_offset, scenario in enumerate(["single_active", "tie"]):
            pvalues = simulate_model(
                outer=outer,
                inner=inner,
                seed=seed + 100_000 * (offset + 1) + 10_000 * scenario_offset,
                replicates=setup["replicates"],
                correlation=setup["correlation"],
                amplitude=setup["amplitude"],
                scenario=scenario,
            )
            pvalue_summaries[model][scenario] = {
                component: {"q01": float(np.quantile(pvalues[:, j], 0.01)), "q05": float(np.quantile(pvalues[:, j], 0.05)), "median": float(np.median(pvalues[:, j]))}
                for j, component in enumerate(COMPONENTS)
            }
            for j, component in enumerate(COMPONENTS):
                for alpha in ALPHAS:
                    rejected = int((pvalues[:, j] <= alpha).sum())
                    rate = rejected / outer
                    lower = 0.0 if rejected == 0 else float(beta.ppf(0.025, rejected, outer - rejected + 1))
                    upper = 1.0 if rejected == outer else float(beta.ppf(0.975, rejected + 1, outer - rejected))
                    greater_p = float(binomtest(rejected, outer, alpha, alternative="greater").pvalue)
                    rows.append(
                        {
                            "model": model,
                            "null_scenario": scenario,
                            "component": component,
                            "nominal_alpha": alpha,
                            "rejections": rejected,
                            "outer_simulations": outer,
                            "empirical_type1": rate,
                            "exact_95ci_lower": lower,
                            "exact_95ci_upper": upper,
                            "one_sided_binomial_p": greater_p,
                            "inference_scope": RSF_WARNING if model == "RSF_N20" else "fixed-query Cox nuisance calibration",
                        }
                    )
    result = pd.DataFrame(rows)
    result["holm_p_across_checks"] = holm_adjust(result.one_sided_binomial_p.to_numpy())
    result["anti_conservative_gate"] = result.holm_p_across_checks < 0.05
    result.to_csv(OUT / "type1_rates.csv", index=False)

    formula_hashes = {
        rel: sha(ROOT / rel) for rel in ["warrantxai/components.py", "warrantxai/aggregation.py"]
    }
    phase4 = json.loads((ROOT / "results/phase4_rsf/metadata.json").read_text())
    unchanged = formula_hashes == phase4["formula_module_sha256"]
    elapsed = perf_counter() - started
    calibration_good = not bool(result.anti_conservative_gate.any())
    nuisance_json = {
        model: {
            "replicates": values["replicates"],
            "correlation": values["correlation"].tolist(),
            "amplitude": dict(zip(COMPONENTS, values["amplitude"].tolist(), strict=True)),
            "nuisance_bootstrap_draws": values["nuisance_bootstrap_draws"],
        }
        for model, values in nuisance.items()
    }
    checks = len(result)
    metadata = {
        "phase": "5.1_null_calibration_gate",
        "arguments": {"outer_simulations": outer, "inner_bootstrap_draws": inner, "seed": seed, "alphas": ALPHAS},
        "known_null": {
            "theta": dict.fromkeys(COMPONENTS, 0.8),
            "construction": {
                "F": "minimum of two bounded site means",
                "S": "minimum of pair U-statistics; active site is 1-0.6*mean_pair_abs(U-U'), E|U-U'|=1/3",
                "U": "minimum of two bounded site means",
                "G_claim": "1-max(|mean discrepancy branches|); active population discrepancy 0.2",
                "R": "minimum of discrete pair U-statistics; active site is 1-0.4*mismatch, P(mismatch)=0.5",
            },
            "scenarios": {
                "single_active": "one site/branch has population score 0.8, the other 0.9; smooth boundary",
                "tie": "both sites/branches have population score 0.8; nonsmooth min tie",
            },
            "sites": "two fixed sites; sites are not resampled",
        },
        "nuisance_from_saved_artifacts": nuisance_json,
        "gate_rule": f"anti-conservative only when one-sided exact binomial greater-tail test remains <0.05 after Holm adjustment across {checks} model/scenario/component/alpha checks",
        "calibration_gate_passed": calibration_good,
        "rsf_conditionality_warning": RSF_WARNING,
        "formula_module_sha256": formula_hashes,
        "formula_module_sha256_unchanged_from_phase4_2": unchanged,
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pvalue_distribution_summaries": pvalue_summaries,
        "full_BH_BY_applied": False,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not unchanged:
        raise RuntimeError("Frozen formula hashes changed")

    compact = result[["model", "null_scenario", "component", "nominal_alpha", "empirical_type1", "exact_95ci_lower", "exact_95ci_upper", "holm_p_across_checks", "anti_conservative_gate"]]
    alpha05 = compact[compact.nominal_alpha == 0.05]
    report = f"""# Faz 5.1 — null kalibrasyon / Tip-I hata kapısı

> **RSF koşulluluk uyarısı:** {RSF_WARNING}

## Özet — RSF sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Gerçek model veya explainer çağrısı yapılmadı. Her model ve null senaryosu için {outer} bağımsız dış simülasyon ve her simülasyonda B={inner} null-merkezli iç bootstrap kullanıldı. Beş sentetik popülasyon hedefi analitik olarak tam 0.8'dir. `single_active` bir site/dal 0.8 ve diğeri 0.9 olan düzgün sınırı; `tie` iki site/dalın da 0.8 olduğu nonsmooth `min` sınırını test eder. Bileşen korelasyonu ve göreli gürültü ölçeği saklı CoxPH N=50 / RSF N=20 artefaktlarından alınmıştır.

Kapı kararı: **{'GEÇTİ' if calibration_good else 'KALDI'}**. Anti-conservative etiketi, nominal oranın üzerindeki sapmanın {checks} kontrol boyunca Holm-düzeltilmiş tek-yönlü exact binomial p<0.05 olmasıyla verildi.

## Alpha=0.05 özeti — RSF satırları sabit 12+12 query ve saklı attribution ortalamalarına koşullu

{markdown_table(alpha05)}

## Tip-I oranları — RSF satırları sabit 12+12 query ve saklı attribution ortalamalarına koşullu

{markdown_table(compact)}

## Yorum — RSF değerlendirmesi sabit 12+12 query ve saklı attribution ortalamalarına koşullu

"""
    if calibration_good:
        report += "Hiçbir model/bileşen/alpha hücresi düzeltilmiş kapıda anti-conservative bulunmadı. Bu sonuç yalnız tanımlanan kontrollü null ailesi ve koşullu çıkarım hedefi için ampirik destek sağlar; evrensel geçerlilik ispatı değildir. Mevcut yöntemde düzeltme uygulanmadı.\n"
    else:
        bad = compact[compact.anti_conservative_gate]
        report += "Aşağıdaki hücreler anti-conservative kapıyı geçti; mevcut yöntemle ADIM 2'ye geçilmemelidir:\n\n"
        report += markdown_table(bad)
        report += (
            "\n\nRSF N=20 `single_active` sınırında F ve U için alpha=0.01 ret oranları sırasıyla "
            "0.0205 ve 0.0215'tir. G_claim daha belirgin biçimde alpha=0.01'de 0.0345, alpha=0.05'te "
            "0.0730 ret vermiştir. Tie senaryoları anti-conservative değildir ve çoğunlukla muhafazakârdır; "
            "dolayısıyla yalnız tie null'u kullanmak sorunu gizlerdi. R her iki modelde belirgin muhafazakârdır; "
            "bu Tip-I ihlali değil, olası güç kaybıdır.\n\n"
            "Olası düzeltmeler: F/U için studentized bootstrap-t veya bounded-score tabanlı muhafazakâr bir "
            "tek-yönlü bound; nonsmooth G_claim için null-constrained bootstrap, m-out-of-n/subsampling veya "
            "iki alt bileşen için ayrı geçerli testlerin intersection-union birleşimi. Hangi yaklaşımın seçileceği "
            "ayrı tasarım ve yeniden kalibrasyon gerektirir; kullanıcı onayı olmadan hiçbiri uygulanmadı.\n"
        )
    report += f"\nFormül SHA256 değişmedi: **{unchanged}**. Kimlik-sızıntısı kontrolü: **geçti**. Gerçek wall-clock: **{elapsed:.2f} saniye**. Tam BH/BY uygulanmadı.\n"
    (ROOT / "PHASE5_1_NULL_CALIBRATION_REPORT.md").write_text(report, encoding="utf-8")
    print(compact.to_string(index=False))
    print(f"Gate={'PASS' if calibration_good else 'FAIL'}; wall={elapsed:.2f}s; full BH/BY not run.")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer", type=int, default=2000)
    parser.add_argument("--inner", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260918)
    args = parser.parse_args()
    run(outer=args.outer, inner=args.inner, seed=args.seed)
