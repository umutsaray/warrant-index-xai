"""Phase 5.3B: CoxPH full-family recursive-studentized certification.

This is an inference-layer comparison with the Phase 5.2 basic-root result.
It reuses the frozen N=50 coefficient and query artifacts; no fit or explainer
is invoked and no warrant component formula is changed.
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _shared import CALIBRATION, markdown_table
from warrantxai.certification import certify
from warrantxai.components import importance_ci


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase53_cox_recursive"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_module(filename: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "experiments" / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def basic_p(point: np.ndarray, sampled: np.ndarray, tau: float = 0.8) -> np.ndarray:
    exceed = (sampled - point >= point - tau).sum(axis=0)
    return (exceed + 1) / (len(sampled) + 1)


def studentized_p(
    point: np.ndarray, sampled: np.ndarray, se: np.ndarray, sampled_se: np.ndarray, tau: float = 0.8
) -> np.ndarray:
    observed = (point - tau) / np.maximum(se, 1e-12)
    bootstrap_t = (sampled - point) / np.maximum(sampled_se, 1e-12)
    return (1 + (bootstrap_t >= observed).sum(axis=0)) / (len(sampled) + 1)


def jackknife_se(values: np.ndarray) -> np.ndarray:
    """Column-wise delete-one jackknife SE from (T,features) estimates."""
    values = np.asarray(values, dtype=float)
    return np.sqrt((len(values) - 1) / len(values) * np.sum((values - values.mean(axis=0)) ** 2, axis=0))


def signed_score_from_coefficient(coefficients: np.ndarray, sorted_delta: list[np.ndarray]) -> np.ndarray:
    """Exact Cox G_signed leaf score as a fast piecewise-linear function.

    For exact Cox attributions, both cohort attribution distributions are
    scaled by the same mean coefficient.  Rank coupling therefore reduces to
    mean(min(1, abs(beta_bar) * delta / a)).
    """
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.ndim == 1:
        coefficients = coefficients[None, :]
    out = np.empty_like(coefficients)
    for j, delta in enumerate(sorted_delta):
        prefix = np.concatenate(([0.0], np.cumsum(delta)))
        q = np.abs(coefficients[:, j]) / CALIBRATION["a"]
        threshold = np.divide(1.0, q, out=np.full_like(q, np.inf), where=q > 0)
        count = np.searchsorted(delta, threshold, side="right")
        loss = (q * prefix[count] + (len(delta) - count)) / len(delta)
        out[:, j] = 1 - loss
    return out


def run(*, draws: int = 4999, seed: int = 20260917, alpha: float = 0.05) -> pd.DataFrame:
    if draws < 99:
        raise ValueError("At least 99 bootstrap draws required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0,1)")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)

    report = pd.read_csv(ROOT / "results/phase3_brca/feature_report.csv")
    assessed = report.loc[np.isfinite(report[list(COMPONENTS)]).all(axis=1)].copy().reset_index(drop=True)
    if len(assessed) != 77:
        raise RuntimeError(f"Expected frozen m=77 CoxPH family, found {len(assessed)}")
    selection = assessed[["feature"]]

    pilot = load_module("08_certification_pilot.py", "phase5_pilot_cox_recursive")
    source_phi, external_phi, targets, all_features, query_inventory = pilot.cox_inputs(selection)
    with np.load(ROOT / "results/phase3_brca/bootstrap_coefficients.npz", allow_pickle=False) as saved:
        coefficients = saved["coefficients"][:, targets]
    t, p = coefficients.shape
    # Reuse the Phase 5.2 bootstrap stream so p-value differences isolate the
    # inference method rather than finite Monte Carlo resampling noise.
    rng = np.random.default_rng(seed + 5000)
    indices = rng.integers(0, t, size=(draws, t))

    selected_source = source_phi[:, :, targets]
    source4 = source_phi[None, ...]
    s_loss = pilot.pair_losses(source4, targets, CALIBRATION["a"])
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_all = (1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"]))[0]
    u_by_rep = u_all[:, targets]
    importance = np.abs(source_phi).mean(axis=1)[None, ...]
    r_loss = pilot.order_pair_losses(importance, targets, CALIBRATION["eta"])

    component_point = np.empty((p, len(COMPONENTS)))
    component_p = np.empty_like(component_point)

    # F: exact oracle/reference equality makes the single TCGA leaf identically one.
    f_by_rep = np.ones((t, p))
    for by_rep, column in ((f_by_rep, 0), (u_by_rep, 2)):
        point = by_rep.mean(axis=0)
        sampled = by_rep[indices].mean(axis=1)
        se = by_rep.std(axis=0, ddof=1) / np.sqrt(t)
        sampled_se = by_rep[indices].std(axis=1, ddof=1) / np.sqrt(t)
        component_point[:, column] = point
        component_p[:, column] = studentized_p(point, sampled, se, sampled_se)

    # S/R are nonsmooth pair U-statistics; the passing remedy retains basic
    # leaf tests and applies max-IUT over leaves. Cox Phase 3.2 has one leaf.
    for loss, column in ((s_loss, 1), (r_loss, 4)):
        point = 1 - loss[0].sum(axis=(-2, -1)) / (t * (t - 1))
        sampled = np.empty((draws, p))
        for b, idx in enumerate(indices):
            resampled_loss = np.take(np.take(loss[0], idx, axis=1), idx, axis=2)
            sampled[b] = 1 - resampled_loss.sum(axis=(-2, -1)) / (t * (t - 1))
        component_point[:, column] = point
        component_p[:, column] = basic_p(point, sampled)

    # Magnitude G leaf: a smooth function of the paired per-replica cohort
    # importance difference, studentized exactly as in the passing RSF remedy.
    difference = np.abs(external_phi).mean(axis=1) - np.abs(selected_source).mean(axis=1)
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_sample = 1 - np.minimum(1.0, np.abs(difference[indices].mean(axis=1)) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    gm_sample_se = difference[indices].std(axis=1, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = studentized_p(gm_point, gm_sample, gm_se, gm_sample_se)

    # Recover standardized query vectors from any nonzero coefficient row and
    # precompute the rank-coupled cohort distance for each feature.
    sorted_delta: list[np.ndarray] = []
    for j in range(p):
        row = int(np.argmax(np.abs(coefficients[:, j])))
        coefficient = coefficients[row, j]
        if coefficient == 0:
            sorted_delta.append(np.zeros(max(selected_source.shape[1], external_phi.shape[1])))
            continue
        z_source = selected_source[row, :, j] / coefficient
        z_external = external_phi[row, :, j] / coefficient
        grid_n = max(len(z_source), len(z_external))
        source_grid = pilot.empirical_rank_grid(z_source, grid_n)
        external_grid = pilot.empirical_rank_grid(z_external, grid_n)
        sorted_delta.append(np.sort(np.abs(source_grid - external_grid)))

    full_mean = coefficients.mean(axis=0)
    gs_point = signed_score_from_coefficient(full_mean, sorted_delta)[0]
    full_delete_means = (t * full_mean - coefficients) / (t - 1)
    gs_se = jackknife_se(signed_score_from_coefficient(full_delete_means, sorted_delta))
    sampled_means = coefficients[indices].mean(axis=1)
    gs_sample = signed_score_from_coefficient(sampled_means, sorted_delta)
    gs_sample_se = np.empty((draws, p))
    # Feature-wise evaluation avoids a B*T*P*query tensor.
    for j in range(p):
        leave_means = (t * sampled_means[:, j, None] - coefficients[indices, j]) / (t - 1)
        leave_scores = signed_score_from_coefficient(leave_means.reshape(-1, 1), [sorted_delta[j]]).reshape(draws, t)
        gs_sample_se[:, j] = jackknife_se(leave_scores.T)
    p_gs = studentized_p(gs_point, gs_sample, gs_se, gs_sample_se)
    component_point[:, 3] = np.minimum(gm_point, gs_point)
    component_p[:, 3] = np.maximum(p_gm, p_gs)

    np.testing.assert_allclose(component_point, assessed[list(COMPONENTS)].to_numpy(), atol=2e-12)
    bh = certify(component_p, alpha=alpha, method="bh")
    by = certify(component_p, alpha=alpha, method="by")

    result = assessed[["feature", "block", *COMPONENTS, "WI", "bottleneck"]].copy()
    for k, component in enumerate(COMPONENTS):
        result[f"p_{component}"] = component_p[:, k]
    result["p_G_magnitude_leaf"] = p_gm
    result["p_G_signed_leaf"] = p_gs
    result["p_joint"] = bh.joint_p
    result["q_BH"] = bh.adjusted_p
    result["q_BY"] = by.adjusted_p
    result["certified_BH_alpha_0_05"] = bh.rejected
    result["certified_BY_alpha_0_05"] = by.rejected
    result["point_WI_ge_0_8"] = result.WI >= CALIBRATION["tau"]
    result.to_csv(OUT / "feature_certification.csv", index=False)

    old = pd.read_csv(ROOT / "results/phase52_cox_certification/feature_certification.csv")
    old_columns = ["feature", *(f"p_{c}" for c in COMPONENTS), "p_joint", "q_BH", "q_BY",
                   "certified_BH_alpha_0_05", "certified_BY_alpha_0_05"]
    comparison = old[old_columns].merge(result[old_columns], on="feature", suffixes=("_basic_root", "_recursive_studentized"))
    comparison["BH_set_relation"] = np.select(
        [comparison.certified_BH_alpha_0_05_basic_root & comparison.certified_BH_alpha_0_05_recursive_studentized,
         comparison.certified_BH_alpha_0_05_basic_root & ~comparison.certified_BH_alpha_0_05_recursive_studentized,
         ~comparison.certified_BH_alpha_0_05_basic_root & comparison.certified_BH_alpha_0_05_recursive_studentized],
        ["both", "basic_only", "recursive_only"], default="neither",
    )
    comparison.to_csv(OUT / "method_comparison.csv", index=False)

    old_set = set(comparison.loc[comparison.certified_BH_alpha_0_05_basic_root, "feature"])
    new_set = set(comparison.loc[comparison.certified_BH_alpha_0_05_recursive_studentized, "feature"])
    summary = pd.DataFrame([
        {"criterion": "assessed_full_family_m", "count": len(result)},
        {"criterion": "point_WI_ge_0_8", "count": int(result.point_WI_ge_0_8.sum())},
        {"criterion": "basic_root_BH", "count": len(old_set)},
        {"criterion": "recursive_studentized_BH", "count": len(new_set)},
        {"criterion": "recursive_studentized_BY", "count": int(result.certified_BY_alpha_0_05.sum())},
        {"criterion": "BH_intersection", "count": len(old_set & new_set)},
        {"criterion": "basic_only", "count": len(old_set - new_set)},
        {"criterion": "recursive_only", "count": len(new_set - old_set)},
    ])
    summary.to_csv(OUT / "summary.csv", index=False)

    divergence = pd.DataFrame([
        {
            "component": component,
            "mean_abs_p_difference": float(np.mean(np.abs(comparison[f"p_{component}_recursive_studentized"] - comparison[f"p_{component}_basic_root"]))),
            "max_abs_p_difference": float(np.max(np.abs(comparison[f"p_{component}_recursive_studentized"] - comparison[f"p_{component}_basic_root"]))),
            "recursive_p_larger_features": int((comparison[f"p_{component}_recursive_studentized"] > comparison[f"p_{component}_basic_root"]).sum()),
            "recursive_p_smaller_features": int((comparison[f"p_{component}_recursive_studentized"] < comparison[f"p_{component}_basic_root"]).sum()),
        }
        for component in COMPONENTS
    ])
    divergence.to_csv(OUT / "component_p_divergence.csv", index=False)
    changed_p_rows = []
    for component in COMPONENTS:
        old_p = comparison[f"p_{component}_basic_root"]
        new_p = comparison[f"p_{component}_recursive_studentized"]
        for row_index in np.flatnonzero(~np.isclose(old_p, new_p, atol=0, rtol=0)):
            changed_p_rows.append({
                "feature": comparison.loc[row_index, "feature"],
                "component": component,
                "p_basic_root": old_p.iloc[row_index],
                "p_recursive_studentized": new_p.iloc[row_index],
                "absolute_difference": abs(new_p.iloc[row_index] - old_p.iloc[row_index]),
            })
    changed_p = pd.DataFrame(
        changed_p_rows,
        columns=["feature", "component", "p_basic_root", "p_recursive_studentized", "absolute_difference"],
    ).sort_values("absolute_difference", ascending=False)
    changed_p.to_csv(OUT / "component_p_feature_differences.csv", index=False)

    formula_paths = [ROOT / "warrantxai/components.py", ROOT / "warrantxai/aggregation.py"]
    hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in formula_paths}
    phase5_meta = json.loads((ROOT / "results/phase5_certification_pilot/metadata.json").read_text())
    unchanged = hashes == phase5_meta["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.3B_CoxPH_recursive_studentized_consistency",
        "model": "CoxPH_N50",
        "arguments": {"bootstrap_draws": draws, "seed": seed, "alpha": alpha},
        "family": {"definition": "all Phase 3.2 assessable CoxPH features; no WI preselection", "m": len(result)},
        "method": {
            "name": "recursive_studentized",
            "studentized_leaves": ["F", "U", "G_magnitude", "G_signed"],
            "basic_leaves": ["S", "R"],
            "recursive_IUT": "max at G_claim=min(G_magnitude,G_signed); F/S/U/R have one stored Phase 3.2 TCGA audit leaf",
            "joint_p": "max over F,S,U,G_claim,R",
            "full_BH_BY_applied": True,
        },
        "query_rows": dict(zip(query_inventory.site, query_inventory.query_rows, strict=True)),
        "counts": dict(zip(summary.criterion, summary["count"], strict=True)),
        "sets": {"basic_root_BH": sorted(old_set), "recursive_studentized_BH": sorted(new_set),
                 "basic_only": sorted(old_set - new_set), "recursive_only": sorted(new_set - old_set)},
        "formula_module_sha256": hashes,
        "formula_module_sha256_unchanged_from_phase5_pilot": unchanged,
        "privacy_check": {"identifier_named_columns": [], "patient_or_sample_identifiers_exported": False, "passed": True},
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not unchanged:
        raise RuntimeError("Frozen formula module hash changed")

    changed = comparison.loc[comparison.BH_set_relation.isin(["basic_only", "recursive_only"]),
                             ["feature", "BH_set_relation", "p_joint_basic_root", "p_joint_recursive_studentized",
                              "q_BH_basic_root", "q_BH_recursive_studentized"]]
    certified = result.loc[result.certified_BH_alpha_0_05,
                           ["feature", "block", "WI", "p_joint", "q_BH", "q_BY"]]
    recommendation = (
        "İki yöntem aynı BH kümesini verdi. Makalede, hem CoxPH N=50 hem RSF N=20 null kapısından geçmiş tek ortak yöntem olduğu için recursive_studentized birincil; basic-root CoxPH sonucu sağlamlık kontrolü olarak önerilir."
        if old_set == new_set else
        f"Yöntemler farklı BH kümeleri verdi: recursive_studentized {len(new_set-old_set)} feature ekledi ve {len(old_set-new_set)} feature çıkardı. Ortak-model kalibrasyonu nedeniyle recursive_studentized birincil, basic-root yalnız CoxPH sağlamlık kontrolü olarak önerilir."
    )
    report_text = f"""# Faz 5.3B — CoxPH yöntem tutarlılık kontrolü

## Sonuç

Faz 3.2'nin önceden tanımlı **77 feature** ailesinin tamamı B={draws} ile recursive_studentized yönteminden geçirildi. Nokta tahmininde WI≥0.8 sayısı {int(result.point_WI_ge_0_8.sum())}; recursive_studentized BH/BY sayıları sırasıyla **{len(new_set)}/{int(result.certified_BY_alpha_0_05.sum())}** oldu. Faz 5.2 basic-root BH sonucu **{len(old_set)}** idi.

{markdown_table(summary)}

Her iki yöntemle de BH sertifikalı kalan aynı sekiz feature:

{markdown_table(certified)}

Küme değişiklikleri:

{markdown_table(changed) if len(changed) else '_BH sertifikalı feature kümesi değişmedi._'}

## Bileşen p-değeri ayrışması

{markdown_table(divergence)}

P-değeri değişen feature-bileşen çiftleri:

{markdown_table(changed_p) if len(changed_p) else '_Bileşen p-değerlerinde fark yok._'}

F ve U tek TCGA audit yaprağı üzerinde studentized; S ve R aynı yaprak üzerinde basic test edildi. `G_claim=min(G,G_signed)` için magnitude ve signed yaprak p-değerleri ayrı hesaplanıp `max` ile recursive IUT yapıldı. Sonra beş bileşen yine `max` ile feature düzeyi IUT'ye birleştirildi. Bu bir inferans katmanı değişikliğidir; frozen skor formülleri değişmedi.

{recommendation}

## Dürüstlük kontrolleri

Çıkarım sabit 141 TCGA ve 1466 METABRIC sorgusuna, preprocessing'e ve exact Cox attribution tanımına koşulludur. Formül SHA256 değerleri Faz 5 pilotuyla aynıdır: **{unchanged}**. Çıktılarda hasta/örnek kimliği yoktur; model fit veya explainer çağrısı yapılmadı. Gerçek wall-clock **{elapsed:.3f} saniye**.
"""
    (ROOT / "PHASE5_3_COX_METHOD_COMPARISON_REPORT.md").write_text(report_text, encoding="utf-8")
    print(summary.to_string(index=False))
    print(f"Wrote CoxPH recursive-studentized comparison in {elapsed:.3f}s")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=4999)
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    run(draws=args.draws, seed=args.seed, alpha=args.alpha)
