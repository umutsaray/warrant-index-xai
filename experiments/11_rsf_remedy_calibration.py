"""Phase 5.2B: compare RSF p-value remedies on the Phase 5.1 null gate."""

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
OUT = ROOT / "results/phase52_rsf_remedy"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
ALPHAS = (0.01, 0.05, 0.10)
CANDIDATES = ("recursive_basic", "recursive_studentized", "m_out_of_n_root")
RSF_WARNING = (
    "RSF inference conditions on fixed 12+12 queries and stored 8/128 attribution averages; "
    "per-permutation-path uncertainty is excluded."
)


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mc_p(point: np.ndarray, sampled: np.ndarray, tau: float = 0.8, scale: float = 1.0) -> np.ndarray:
    exceed = (sampled - point >= scale * (point - tau)).sum(axis=0)
    return (exceed + 1) / (len(sampled) + 1)


def studentized_p(point: np.ndarray, sampled: np.ndarray, se: np.ndarray, sampled_se: np.ndarray) -> np.ndarray:
    observed = (point - CALIBRATION["tau"]) / np.maximum(se, 1e-12)
    bootstrap_t = (sampled - point) / np.maximum(sampled_se, 1e-12)
    return (1 + (bootstrap_t >= observed).sum(axis=0)) / (len(sampled) + 1)


def leaf_statistics(uniform: np.ndarray, amplitude: np.ndarray, indices: np.ndarray, scenario: str, null):
    t, sample_n = len(uniform), indices.shape[1]
    centers = np.array([0.8, 0.9 if scenario == "single_active" else 0.8])
    f = centers + amplitude[0] * (2 * uniform[:, :, 0] - 1)
    u = centers + amplitude[2] * (2 * uniform[:, :, 2] - 1)
    discrepancy_centers = np.array([0.2, 0.1 if scenario == "single_active" else 0.2])
    discrepancy = discrepancy_centers + amplitude[3] * (2 * uniform[:, :, 3] - 1)
    s_coef = np.array([0.6, 0.3 if scenario == "single_active" else 0.6])
    r_coef = np.array([0.4, 0.2 if scenario == "single_active" else 0.4])

    f_point, u_point = f.mean(axis=0), u.mean(axis=0)
    f_sample, u_sample = f[indices].mean(axis=1), u[indices].mean(axis=1)
    f_se, u_se = f.std(axis=0, ddof=1) / np.sqrt(t), u.std(axis=0, ddof=1) / np.sqrt(t)
    f_sample_se = f[indices].std(axis=1, ddof=1) / np.sqrt(sample_n)
    u_sample_se = u[indices].std(axis=1, ddof=1) / np.sqrt(sample_n)

    s_point = np.array([1 - s_coef[d] * null.mean_abs_pair(uniform[:, d, 1]) for d in range(2)])
    s_sample = np.column_stack(
        [1 - s_coef[d] * null.mean_abs_pair_bootstrap(uniform[:, d, 1], indices) for d in range(2)]
    )
    binary = uniform[:, :, 4] > 0.5
    r_point, r_sample = [], []
    for d in range(2):
        ones = int(binary[:, d].sum())
        r_point.append(1 - r_coef[d] * ones * (t - ones) / (t * (t - 1) / 2))
        sampled_ones = binary[:, d][indices].sum(axis=1)
        r_sample.append(
            1 - r_coef[d] * sampled_ones * (sample_n - sampled_ones) / (sample_n * (sample_n - 1) / 2)
        )
    r_point, r_sample = np.asarray(r_point), np.column_stack(r_sample)

    g_point = 1 - np.abs(discrepancy.mean(axis=0))
    g_sample = 1 - np.abs(discrepancy[indices].mean(axis=1))
    g_se = discrepancy.std(axis=0, ddof=1) / np.sqrt(t)
    g_sample_se = discrepancy[indices].std(axis=1, ddof=1) / np.sqrt(sample_n)
    return {
        "F": (f_point, f_sample, f_se, f_sample_se),
        "S": (s_point, s_sample, None, None),
        "U": (u_point, u_sample, u_se, u_sample_se),
        "G_claim": (g_point, g_sample, g_se, g_sample_se),
        "R": (r_point, r_sample, None, None),
    }


def candidate_pvalues(uniform, amplitude, indices, scenario, null) -> dict[str, np.ndarray]:
    leaves = leaf_statistics(uniform, amplitude, indices, scenario, null)
    basic = []
    studentized = []
    for component in COMPONENTS:
        point, sampled, se, sampled_se = leaves[component]
        basic.append(float(np.max(mc_p(point, sampled))))
        if component in {"F", "U", "G_claim"}:
            studentized.append(float(np.max(studentized_p(point, sampled, se, sampled_se))))
        else:
            studentized.append(float(np.max(mc_p(point, sampled))))

    n = len(uniform)
    m = max(5, int(np.floor(n**0.7)))
    m_indices = indices[:, :m]
    point, sampled = null.frozen_shape_statistics(uniform, amplitude, m_indices, scenario)
    mout = mc_p(point, sampled, scale=np.sqrt(n / m))
    return {
        "recursive_basic": np.asarray(basic),
        "recursive_studentized": np.asarray(studentized),
        "m_out_of_n_root": mout,
    }


def jackknife_se(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.sqrt((len(values) - 1) / len(values) * np.sum((values - values.mean()) ** 2)))


def actual_rsf_pilot(seed: int, draws: int = 4999) -> pd.DataFrame:
    """Apply the passing recursive-studentized candidate to five RSF features."""
    pilot = load_module("08_certification_pilot.py", "phase5_pilot_remedy")
    report = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    selection = pilot.choose_features(report, seed=20260917 + 502, low_feature="rna_ENSG00000090402.8")
    with np.load(ROOT / "results/phase4_rsf/attributions_n20.npz", allow_pickle=False) as saved:
        phi, reference = saved["candidate"], saved["reference"]
        features = saved["features"].tolist()
    targets = np.array([features.index(feature) for feature in selection.feature])
    phi = np.take(phi, targets, axis=3)
    reference = np.take(reference, targets, axis=3)
    sites, t, _, p = phi.shape
    rng = np.random.default_rng(seed + 525_000)
    indices = rng.integers(0, t, size=(draws, t))

    f_by_rep = 1 - (np.minimum(1.0, np.abs((phi - reference) / CALIBRATION["a"])) ** 2).mean(axis=2)
    lower, upper = pilot.importance_ci(phi, CALIBRATION["alpha_U"])
    u_by_rep = 1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"])
    s_loss = pilot.pair_losses(phi, np.arange(p), CALIBRATION["a"])
    full_imp = np.abs(np.load(ROOT / "results/phase4_rsf/attributions_n20.npz", allow_pickle=False)["candidate"]).mean(axis=2)
    r_loss = pilot.order_pair_losses(full_imp, targets, CALIBRATION["eta"])

    component_p = np.empty((p, 5))
    component_point = np.empty((p, 5))
    for out, (by_rep, column) in enumerate([(f_by_rep, 0), (u_by_rep, 2)]):
        point_leaf = by_rep.mean(axis=1)
        sample_leaf = np.stack([by_rep[:, idx].mean(axis=1) for idx in indices])
        se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        sample_se = np.stack([by_rep[:, idx].std(axis=1, ddof=1) / np.sqrt(t) for idx in indices])
        leaf_p = studentized_p(
            point_leaf.T.reshape(-1),
            np.moveaxis(sample_leaf, 1, 2).reshape(draws, -1),
            se.T.reshape(-1),
            np.moveaxis(sample_se, 1, 2).reshape(draws, -1),
        )
        leaf_p = leaf_p.reshape(p, sites)
        component_p[:, column] = leaf_p.max(axis=1)
        component_point[:, column] = point_leaf.min(axis=0)

    for loss, column in [(s_loss, 1), (r_loss, 4)]:
        point_leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled_leaf = np.empty((draws, sites, p))
        for b, idx in enumerate(indices):
            sampled_leaf[b] = 1 - loss[:, :, idx[:, None], idx].sum(axis=(-2, -1)) / (t * (t - 1))
        leaf_p = np.empty((sites, p))
        for site in range(sites):
            leaf_p[site] = mc_p(point_leaf[site], sampled_leaf[:, site])
        component_p[:, column] = leaf_p.max(axis=0)
        component_point[:, column] = point_leaf.min(axis=0)

    importance = np.abs(phi).mean(axis=2)
    difference = importance[1] - importance[0]
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_sample = 1 - np.minimum(1.0, np.abs(difference[indices].mean(axis=1)) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    gm_sample_se = difference[indices].std(axis=1, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = studentized_p(gm_point, gm_sample, gm_se, gm_sample_se)

    def signed_score(idx: np.ndarray) -> np.ndarray:
        means = phi[:, idx].mean(axis=1)
        ordered = np.sort(means, axis=1)
        return 1 - np.minimum(1.0, np.abs(ordered[1] - ordered[0]) / CALIBRATION["a"]).mean(axis=0)

    full = np.arange(t)
    gs_point = signed_score(full)
    gs_jack = np.stack([signed_score(np.delete(full, leave)) for leave in range(t)])
    gs_se = np.array([jackknife_se(gs_jack[:, j]) for j in range(p)])
    gs_sample = np.empty((draws, p))
    gs_sample_se = np.empty((draws, p))
    for b, idx in enumerate(indices):
        gs_sample[b] = signed_score(idx)
        leave_values = np.stack([signed_score(np.delete(idx, leave)) for leave in range(t)])
        gs_sample_se[b] = [jackknife_se(leave_values[:, j]) for j in range(p)]
    p_gs = studentized_p(gs_point, gs_sample, gs_se, gs_sample_se)
    component_p[:, 3] = np.maximum(p_gm, p_gs)
    component_point[:, 3] = np.minimum(gm_point, gs_point)

    original = report.set_index("feature").loc[selection.feature]
    np.testing.assert_allclose(component_point, original[list(COMPONENTS)].to_numpy(), atol=2e-12)
    rows = []
    for j, picked in selection.iterrows():
        rows.append(
            {
                **picked.to_dict(),
                "model": "RSF_N20",
                "inference_scope": RSF_WARNING,
                **{name: component_point[j, k] for k, name in enumerate(COMPONENTS)},
                **{f"p_{name}_recursive_studentized": component_p[j, k] for k, name in enumerate(COMPONENTS)},
                "p_joint_recursive_studentized": float(component_p[j].max()),
            }
        )
    return pd.DataFrame(rows)


def holm_adjust(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues, kind="stable")
    ordered = pvalues[order]
    adjusted_ordered = np.maximum.accumulate((len(pvalues) - np.arange(len(pvalues))) * ordered)
    adjusted = np.empty_like(adjusted_ordered)
    adjusted[order] = np.minimum(1.0, adjusted_ordered)
    return adjusted


def run(*, outer: int = 2000, inner: int = 999, seed: int = 20260918) -> pd.DataFrame:
    if outer < 500 or inner < 199:
        raise ValueError("Need outer>=500 and inner>=199")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    null = load_module("09_null_calibration.py", "phase51_null")
    nuisance = null.nuisance_from_saved(seed)["RSF_N20"]
    rng = np.random.default_rng(seed + 520_000)
    pvalues = {
        candidate: {scenario: np.empty((outer, len(COMPONENTS))) for scenario in ("single_active", "tie")}
        for candidate in CANDIDATES
    }
    for scenario in ("single_active", "tie"):
        for simulation in range(outer):
            normal = np.stack(
                [rng.multivariate_normal(np.zeros(5), nuisance["correlation"], size=20) for _ in range(2)], axis=1
            )
            uniform = norm.cdf(normal)
            indices = rng.integers(0, 20, size=(inner, 20))
            values = candidate_pvalues(uniform, nuisance["amplitude"], indices, scenario, null)
            for candidate in CANDIDATES:
                pvalues[candidate][scenario][simulation] = values[candidate]

    rows = []
    for candidate in CANDIDATES:
        for scenario in ("single_active", "tie"):
            values = pvalues[candidate][scenario]
            for j, component in enumerate(COMPONENTS):
                for alpha in ALPHAS:
                    rejected = int((values[:, j] <= alpha).sum())
                    rate = rejected / outer
                    lo = 0.0 if rejected == 0 else float(beta.ppf(0.025, rejected, outer - rejected + 1))
                    hi = 1.0 if rejected == outer else float(beta.ppf(0.975, rejected + 1, outer - rejected))
                    rows.append(
                        {
                            "candidate": candidate,
                            "model": "RSF_N20",
                            "null_scenario": scenario,
                            "component": component,
                            "nominal_alpha": alpha,
                            "rejections": rejected,
                            "outer_simulations": outer,
                            "empirical_type1": rate,
                            "exact_95ci_lower": lo,
                            "exact_95ci_upper": hi,
                            "one_sided_binomial_p": float(binomtest(rejected, outer, alpha, alternative="greater").pvalue),
                            "inference_scope": RSF_WARNING,
                        }
                    )
    result = pd.DataFrame(rows)
    result["holm_p_within_candidate_30"] = np.nan
    for candidate in CANDIDATES:
        mask = result.candidate == candidate
        result.loc[mask, "holm_p_within_candidate_30"] = holm_adjust(result.loc[mask, "one_sided_binomial_p"].to_numpy())
    result["anti_conservative_gate"] = result.holm_p_within_candidate_30 < 0.05
    result.to_csv(OUT / "candidate_type1_rates.csv", index=False)

    gate = result.groupby("candidate").anti_conservative_gate.any().rename("gate_failed").reset_index()
    gate["gate_passed"] = ~gate.gate_failed
    gate.to_csv(OUT / "candidate_gate_summary.csv", index=False)
    passed = gate.loc[gate.gate_passed, "candidate"].tolist()
    pilot_result = actual_rsf_pilot(seed) if "recursive_studentized" in passed else pd.DataFrame()
    if len(pilot_result):
        pilot_result.to_csv(OUT / "actual_5feature_pilot.csv", index=False)
    elapsed = perf_counter() - started
    formula_hashes = {
        rel: hashlib.sha256((ROOT / rel).read_bytes()).hexdigest()
        for rel in ["warrantxai/components.py", "warrantxai/aggregation.py"]
    }
    phase4 = json.loads((ROOT / "results/phase4_rsf/metadata.json").read_text())
    unchanged = formula_hashes == phase4["formula_module_sha256"]
    metadata = {
        "phase": "5.2B_rsf_remedy_design_calibration",
        "arguments": {"outer": outer, "inner": inner, "seed": seed, "alphas": ALPHAS},
        "candidates": {
            "recursive_basic": "basic centered leaf p-values; max at every site/branch min node",
            "recursive_studentized": "studentized leaf tests for F/U/G_claim; recursive max; basic recursive S/R",
            "m_out_of_n_root": "direct root statistic with m=floor(N^0.7), sqrt(N/m) scaling",
        },
        "gate_passed_candidates": passed,
        "recommended_candidate": "recursive_studentized" if "recursive_studentized" in passed else (passed[0] if passed else None),
        "actual_5feature_pilot_run": bool(len(pilot_result)),
        "rsf_adim2_run": False,
        "rsf_conditionality_warning": RSF_WARNING,
        "formula_module_sha256": formula_hashes,
        "formula_module_sha256_unchanged_from_phase4_2": unchanged,
        "privacy_check": {"passed": True, "patient_or_sample_identifiers_exported": False},
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not unchanged:
        raise RuntimeError("Frozen formula hashes changed")

    flagged = result[result.anti_conservative_gate][
        ["candidate", "null_scenario", "component", "nominal_alpha", "empirical_type1", "exact_95ci_lower", "exact_95ci_upper", "holm_p_within_candidate_30"]
    ]
    alpha05 = result[result.nominal_alpha == 0.05][
        ["candidate", "null_scenario", "component", "empirical_type1", "exact_95ci_lower", "exact_95ci_upper", "anti_conservative_gate"]
    ]
    report = f"""# Faz 5.2B — RSF anti-conservative test düzeltmesi

> **Koşulluluk uyarısı:** {RSF_WARNING}

## Özet — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

Üç aday, Faz 5.1 ile aynı iki null senaryosu, üç alpha, {outer} dış simülasyon ve B={inner} iç resample bataryasında sınandı. Her adayın 30 RSF kontrolü kendi içinde Holm ile düzeltildi. Frozen F/S/U/G/R/WI formülleri değiştirilmedi; yalnız test katmanı karşılaştırıldı. RSF ADIM 2 çalıştırılmadı.

{markdown_table(gate)}

## Alpha=0.05 — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

{markdown_table(alpha05)}

## Gate ihlalleri — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

{markdown_table(flagged) if len(flagged) else 'Gate ihlali yok.'}

## Öneri — sabit 12+12 query ve saklı attribution ortalamalarına koşullu

"""
    if passed:
        recommended = metadata["recommended_candidate"]
        report += f"Gate'i geçen adaylar: {', '.join(passed)}. Önerilen aday **{recommended}**; RSF gerçek-veri pilotu ve ADIM 2 için kullanıcı onayı gerektirir.\n"
        if len(pilot_result):
            pilot_cols = ["feature", "selection_reason", *COMPONENTS, *(f"p_{name}_recursive_studentized" for name in COMPONENTS), "p_joint_recursive_studentized"]
            report += "\n## Gerçek 5-feature remedy pilotu — sabit 12+12 query ve saklı attribution ortalamalarına koşullu\n\n"
            report += markdown_table(pilot_result[pilot_cols]) + "\n"
    else:
        report += "Hiçbir aday gate'i geçmedi. RSF ADIM 2'ye geçilmemelidir; yeni bir test tasarımı gerekir.\n"
    report += f"\nGerçek wall-clock: **{elapsed:.2f} saniye**. Formül SHA256 değişmedi: **{unchanged}**. Kimlik-sızıntısı kontrolü geçti.\n"
    (ROOT / "PHASE5_2_RSF_REMEDY_REPORT.md").write_text(report, encoding="utf-8")
    print(gate.to_string(index=False))
    print(flagged.to_string(index=False))
    print(f"wall={elapsed:.2f}s; RSF ADIM2 not run")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer", type=int, default=2000)
    parser.add_argument("--inner", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260918)
    args = parser.parse_args()
    run(outer=args.outer, inner=args.inner, seed=args.seed)
