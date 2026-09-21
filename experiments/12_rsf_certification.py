"""Phase 5.3A: full-family RSF N=20 recursive-studentized certification.

This reads the stored Phase 4.2 attribution averages.  It never fits an RSF
or calls the approximate explainer.  Inference therefore remains conditional
on the fixed queries and the stored permutation-path averages.
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


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase53_rsf_certification"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
WARNING = (
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


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def resampled_means(values: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return bootstrap means and their ordinary within-resample SEs.

    ``values`` is (site,T,feature), and ``counts`` is (B,T).
    """
    t = values.shape[1]
    mean = np.einsum("bt,stp->bsp", counts, values, optimize=True) / t
    second = np.einsum("bt,stp->bsp", counts, values**2, optimize=True)
    variance = np.maximum(0.0, (second - t * mean**2) / (t - 1))
    return mean, np.sqrt(variance / t)


def quadratic_scores(loss: np.ndarray, counts: np.ndarray, t: int) -> np.ndarray:
    """Recompute ordered-pair S/R leaves for every bootstrap multiset."""
    total = np.einsum("bt,sptu,bu->bsp", counts, loss, counts, optimize=True)
    return 1 - total / (t * (t - 1))


def signed_scores(phi: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Quantile-coupled signed score for arbitrary replicate weights.

    ``weights`` is (draw,T).  Output is (draw,feature).
    """
    denominator = weights.sum(axis=1)[:, None, None, None]
    means = np.einsum("bt,stqp->bsqp", weights, phi, optimize=True) / denominator
    ordered = np.sort(means, axis=2)
    return 1 - np.minimum(
        1.0, np.abs(ordered[:, 1] - ordered[:, 0]) / CALIBRATION["a"]
    ).mean(axis=1)


def signed_bootstrap_with_jackknife_se(
    phi: np.ndarray, indices: np.ndarray, *, batch_size: int = 64
) -> tuple[np.ndarray, np.ndarray]:
    """Signed bootstrap values and within-bootstrap jackknife SEs."""
    draws, t = indices.shape
    p = phi.shape[-1]
    scores = np.empty((draws, p))
    ses = np.empty((draws, p))
    eye = np.eye(t, dtype=np.int16)
    for start in range(0, draws, batch_size):
        stop = min(draws, start + batch_size)
        idx = indices[start:stop]
        counts = np.zeros((len(idx), t), dtype=np.int16)
        np.add.at(counts, (np.arange(len(idx))[:, None], idx), 1)
        scores[start:stop] = signed_scores(phi, counts)

        # Deleting a bootstrap *position* is the jackknife used in the approved
        # Phase 5.2 pilot. Duplicate source indices consequently give repeated
        # leave-one-out values, as they should for the empirical bootstrap.
        leave_counts = counts[:, None, :] - eye[idx]
        flat = leave_counts.reshape(-1, t)
        leave = signed_scores(phi, flat).reshape(len(idx), t, p)
        centered = leave - leave.mean(axis=1, keepdims=True)
        ses[start:stop] = np.sqrt((t - 1) / t * np.sum(centered**2, axis=1))
    return scores, ses


def component_statistics(
    phi: np.ndarray, reference: np.ndarray, indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return frozen point estimates and recursive-studentized component p-values."""
    remedy = load_module("11_rsf_remedy_calibration.py", "phase52_rsf_remedy_for_full")
    pilot = load_module("08_certification_pilot.py", "phase5_pilot_for_rsf_full")
    sites, t, _, p = phi.shape
    draws = len(indices)
    counts = np.zeros((draws, t), dtype=np.int16)
    np.add.at(counts, (np.arange(draws)[:, None], indices), 1)
    points = np.empty((p, len(COMPONENTS)))
    pvalues = np.empty_like(points)

    f_by_rep = 1 - (np.minimum(1.0, np.abs((phi - reference) / CALIBRATION["a"])) ** 2).mean(axis=2)
    lower, upper = pilot.importance_ci(phi, CALIBRATION["alpha_U"])
    u_by_rep = 1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"])
    for by_rep, column in ((f_by_rep, 0), (u_by_rep, 2)):
        point_leaf = by_rep.mean(axis=1)
        sampled_leaf, sampled_se = resampled_means(by_rep, counts)
        se = by_rep.std(axis=1, ddof=1) / np.sqrt(t)
        leaf_p = remedy.studentized_p(
            point_leaf.T.reshape(-1),
            np.moveaxis(sampled_leaf, 1, 2).reshape(draws, -1),
            se.T.reshape(-1),
            np.moveaxis(sampled_se, 1, 2).reshape(draws, -1),
        ).reshape(p, sites)
        pvalues[:, column] = leaf_p.max(axis=1)
        points[:, column] = point_leaf.min(axis=0)

    all_targets = np.arange(p)
    importance = np.abs(phi).mean(axis=2)
    s_loss = pilot.pair_losses(phi, all_targets, CALIBRATION["a"])
    r_loss = pilot.order_pair_losses(importance, all_targets, CALIBRATION["eta"])
    for loss, column in ((s_loss, 1), (r_loss, 4)):
        point_leaf = 1 - loss.sum(axis=(-2, -1)) / (t * (t - 1))
        sampled_leaf = quadratic_scores(loss, counts, t)
        leaf_p = np.column_stack(
            [remedy.mc_p(point_leaf[site], sampled_leaf[:, site]) for site in range(sites)]
        )
        pvalues[:, column] = leaf_p.max(axis=1)
        points[:, column] = point_leaf.min(axis=0)

    difference = importance[1] - importance[0]
    gm_point = 1 - np.minimum(1.0, np.abs(difference.mean(axis=0)) / CALIBRATION["b"])
    gm_mean, gm_mean_se = resampled_means(difference[None, ...], counts)
    gm_sample = 1 - np.minimum(1.0, np.abs(gm_mean[:, 0]) / CALIBRATION["b"])
    gm_se = difference.std(axis=0, ddof=1) / np.sqrt(t) / CALIBRATION["b"]
    p_gm = remedy.studentized_p(
        gm_point, gm_sample, gm_se, gm_mean_se[:, 0] / CALIBRATION["b"]
    )

    full_weights = np.ones((1, t), dtype=np.int16)
    gs_point = signed_scores(phi, full_weights)[0]
    full_leave = np.ones((t, t), dtype=np.int16) - np.eye(t, dtype=np.int16)
    gs_jack = signed_scores(phi, full_leave)
    gs_centered = gs_jack - gs_jack.mean(axis=0, keepdims=True)
    gs_se = np.sqrt((t - 1) / t * np.sum(gs_centered**2, axis=0))
    gs_sample, gs_sample_se = signed_bootstrap_with_jackknife_se(phi, indices)
    p_gs = remedy.studentized_p(gs_point, gs_sample, gs_se, gs_sample_se)
    pvalues[:, 3] = np.maximum(p_gm, p_gs)
    points[:, 3] = np.minimum(gm_point, gs_point)
    return points, pvalues


def run(*, draws: int = 4999, seed: int = 20260918, alpha: float = 0.05) -> pd.DataFrame:
    if draws < 99:
        raise ValueError("At least 99 bootstrap draws required")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0,1)")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    stored = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    with np.load(ROOT / "results/phase4_rsf/attributions_n20.npz", allow_pickle=False) as saved:
        phi, reference, features = saved["candidate"], saved["reference"], saved["features"].tolist()
    if len(stored) != 77 or features != stored.feature.tolist():
        raise RuntimeError("Expected the frozen ordered m=77 RSF family")

    rng = np.random.default_rng(seed + 530_000)
    indices = rng.integers(0, phi.shape[1], size=(draws, phi.shape[1]))
    points, component_p = component_statistics(phi, reference, indices)
    np.testing.assert_allclose(points, stored[list(COMPONENTS)].to_numpy(), atol=2e-12)
    bh = certify(component_p, alpha=alpha, method="bh")
    by = certify(component_p, alpha=alpha, method="by")

    result = stored[["feature", "block", *COMPONENTS, "WI", "bottleneck"]].copy()
    result.insert(0, "inference_scope", WARNING)
    for k, component in enumerate(COMPONENTS):
        result[f"p_{component}"] = component_p[:, k]
    result["p_joint"] = bh.joint_p
    result["q_BH"] = bh.adjusted_p
    result["q_BY"] = by.adjusted_p
    result["certified_BH_alpha_0_05"] = bh.rejected
    result["certified_BY_alpha_0_05"] = by.rejected
    result["point_WI_ge_0_8"] = result.WI >= CALIBRATION["tau"]
    result["comparison_BH"] = np.select(
        [result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05,
         result.point_WI_ge_0_8 & ~result.certified_BH_alpha_0_05,
         ~result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05],
        ["both", "point_only", "BH_only"], default="neither",
    )
    result.to_csv(OUT / "feature_certification.csv", index=False)

    counts = [
        ("assessed_full_family_m", len(result)),
        ("point_WI_ge_0_8", int(result.point_WI_ge_0_8.sum())),
        ("BH_certified_alpha_0_05", int(result.certified_BH_alpha_0_05.sum())),
        ("BY_certified_alpha_0_05", int(result.certified_BY_alpha_0_05.sum())),
        ("point_and_BH", int((result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum())),
        ("point_only_not_BH", int((result.point_WI_ge_0_8 & ~result.certified_BH_alpha_0_05).sum())),
        ("BH_only_not_point", int((~result.point_WI_ge_0_8 & result.certified_BH_alpha_0_05).sum())),
    ]
    summary = pd.DataFrame(
        [{"inference_scope": WARNING, "criterion": criterion, "count": count} for criterion, count in counts]
    )
    summary.to_csv(OUT / "summary.csv", index=False)

    formula_paths = [ROOT / "warrantxai/components.py", ROOT / "warrantxai/aggregation.py"]
    formula_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in formula_paths}
    phase52 = json.loads((ROOT / "results/phase52_rsf_remedy/metadata.json").read_text())
    unchanged = formula_hashes == phase52["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.3A_RSF_full_family_certification",
        "model": "RSF_N20",
        "rsf_conditionality_warning": WARNING,
        "arguments": {"bootstrap_draws": draws, "seed": seed, "alpha": alpha, "tau_components": dict.fromkeys(COMPONENTS, 0.8)},
        "family": {"definition": "all Phase 4.2 RSF features; no WI preselection", "m": len(result)},
        "method": {
            "name": "recursive_studentized",
            "studentized_leaves": ["F", "U", "G_magnitude", "G_signed"],
            "basic_bootstrap_leaves": ["S", "R"],
            "recursive_IUT": "max p at every site/branch min node, then max over components",
            "primary_fdr": "BH",
            "robustness_fdr": "BY",
            "full_BH_BY_applied": True,
        },
        "counts": dict(counts),
        "formula_module_sha256": formula_hashes,
        "formula_module_sha256_unchanged_from_phase5_2": unchanged,
        "privacy_check": {"identifier_named_columns": [], "patient_or_sample_identifiers_exported": False, "passed": True},
        "model_or_explainer_calls": False,
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not unchanged:
        raise RuntimeError("Frozen formula hashes changed")
    if any(any(token in column.lower() for token in ("patient", "sample_id", "sample_key")) for column in result.columns):
        raise RuntimeError("Identifier-like output column found")

    passing = result.loc[result.certified_BH_alpha_0_05, ["inference_scope", "feature", "block", "WI", "p_joint", "q_BH", "q_BY", "certified_BY_alpha_0_05"]]
    point_only = result.loc[result.comparison_BH == "point_only", ["inference_scope", "feature", "block", "WI", "bottleneck", "p_joint", "q_BH"]]
    report = f"""# Faz 5.3A — RSF tam-aile recursive_studentized sertifikasyonu

> **Koşulluluk uyarısı:** {WARNING}

## Kapsam ve yöntem — {WARNING}

Faz 4.2'deki **{len(result)} feature'ın tamamı** aileye alındı; WI filtresi uygulanmadı. Gate'i geçen `recursive_studentized` yöntemi B={draws} ile kullanıldı: F/U/G_magnitude/G_signed yapraklarında studentized bootstrap, S/R yapraklarında basic bootstrap ve formül ağacındaki her `min` için recursive `max(p)` IUT. Bileşenler de `p_joint=max_k(p_jk)` ile birleştirildi. BH α={alpha:.2f} birincil, BY bağımlılığa dayanıklı sağlamlık analizidir.

## Nokta tahmini karşılaştırması — {WARNING}

{markdown_table(summary)}

## BH sertifikalı feature'lar — {WARNING}

{markdown_table(passing) if len(passing) else '_Hiçbir feature BH ile sertifikalanmadı._'}

## Nokta-yeterli olup BH ile sertifikalanmayanlar — {WARNING}

{markdown_table(point_only) if len(point_only) else '_Nokta-yeterli olup BH ile sertifikalanmayan feature yok._'}

Bu gruptaki {len(point_only)} feature'ın tamamında nokta-tahmin bottleneck'i `{point_only.bottleneck.mode().iat[0] if len(point_only) else 'yok'}` bileşenidir.

Tam feature-bileşen tablosu `results/phase53_rsf_certification/feature_certification.csv` içindedir. q-değeri feature doğruluk olasılığı değildir.

## Dürüstlük kontrolleri — {WARNING}

Formül SHA256 değerleri Faz 5.2 ile aynıdır: **{unchanged}**. Kimlik-sızıntısı kontrolü geçti. Model yeniden fit edilmedi, explainer çağrılmadı. Gerçek wall-clock **{elapsed:.3f} saniye**.
"""
    (ROOT / "PHASE5_3_RSF_CERTIFICATION_REPORT.md").write_text(report, encoding="utf-8")
    print(summary[["criterion", "count"]].to_string(index=False))
    print(f"Wrote RSF m={len(result)} certification in {elapsed:.3f}s")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=4999)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    run(draws=args.draws, seed=args.seed, alpha=args.alpha)
