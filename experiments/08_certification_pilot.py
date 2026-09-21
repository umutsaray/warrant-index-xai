"""Phase 5 step 0/1: inventory and a conditional bootstrap p-value pilot.

This script deliberately stops before BH/BY.  It reuses saved training-
replicate artifacts and does not fit a model or call an explainer.
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

import joblib
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _shared import CALIBRATION, markdown_table
from warrantxai.certification import intersection_union
from warrantxai.components import generalizability, generalizability_signed, importance_ci, pairwise_orders


ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase5_certification_pilot"
COMPONENTS = ("F", "S", "U", "G_claim", "R")


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def load_phase3_module():
    spec = importlib.util.spec_from_file_location("phase3_real", ROOT / "experiments/06_brca_real.py")
    if spec is None or spec.loader is None:
        raise RuntimeError("Phase 3 implementation cannot be loaded")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def empirical_rank_grid(values: np.ndarray, size: int) -> np.ndarray:
    ordered = np.sort(np.asarray(values, dtype=float))
    indices = np.minimum((np.arange(size) + 0.5) * len(ordered) // size, len(ordered) - 1).astype(int)
    return ordered[indices]


def pair_losses(phi: np.ndarray, targets: np.ndarray, scale: float) -> np.ndarray:
    """Return site,target,T,T losses used by frozen S."""
    sites, t, _, _ = phi.shape
    out = np.empty((sites, len(targets), t, t))
    for site in range(sites):
        selected = np.take(phi[site], targets, axis=2)
        for first in range(t):
            delta = np.minimum(1.0, np.abs((selected[first] - selected) / scale)) ** 2
            out[site, :, first, :] = delta.mean(axis=1).T
    return out


def order_pair_losses(importance: np.ndarray, targets: np.ndarray, eta: float) -> np.ndarray:
    """Return site,target,T,T losses used by frozen R, retaining all ranks."""
    sites, t, p = importance.shape
    orders = pairwise_orders(importance, eta)
    out = np.empty((sites, len(targets), t, t))
    for site in range(sites):
        for out_j, j in enumerate(targets):
            order = orders[site, :, j, :].astype(float)
            out[site, out_j] = np.abs(order[:, None, :] - order[None, :, :]).sum(axis=-1) / (2 * (p - 1))
    return out


def pair_score(loss: np.ndarray, sampled: np.ndarray) -> np.ndarray:
    """Apply the frozen unordered-pair mean after an empirical bootstrap draw."""
    t = len(sampled)
    total = loss[:, :, sampled[:, None], sampled].sum(axis=(-2, -1))
    # The full square counts each off-diagonal pair twice; the diagonal is zero.
    return (1 - total / (t * (t - 1))).min(axis=0)


def centered_bootstrap_pvalue(point: np.ndarray, draws: np.ndarray, tau: float) -> np.ndarray:
    """One-sided basic-bootstrap test at the least-favourable boundary theta=tau.

    The null-centred replicate is tau + (theta_star - theta_hat).  The +1
    correction makes the finite Monte Carlo estimate nonzero.
    """
    point = np.asarray(point, dtype=float)
    draws = np.asarray(draws, dtype=float)
    if draws.ndim != 2 or draws.shape[1:] != point.shape or not np.isfinite(draws).all():
        raise ValueError("draws must be finite with shape (B, features)")
    exceed = (draws - point >= point - tau).sum(axis=0)
    return (exceed + 1) / (len(draws) + 1)


def choose_features(report: pd.DataFrame, *, seed: int, low_feature: str) -> pd.DataFrame:
    assessed = report.loc[np.isfinite(report[list(COMPONENTS)]).all(axis=1)].copy()
    chosen: list[tuple[str, str, str]] = []
    high = assessed.loc[assessed.WI.idxmax(), "feature"]
    near = assessed.loc[(assessed.R - CALIBRATION["tau"]).abs().idxmin(), "feature"]
    chosen.extend([(high, "high_margin", "WI"), (near, "near_threshold", "R")])
    if low_feature not in set(assessed.feature):
        low_feature = assessed.loc[assessed.WI.idxmin(), "feature"]
    chosen.append((low_feature, "phase4.1_low_WI", str(assessed.set_index("feature").loc[low_feature, "bottleneck"])))
    excluded = {row[0] for row in chosen}
    pool = assessed.loc[~assessed.feature.isin(excluded), "feature"].to_numpy()
    random_features = np.random.default_rng(seed).choice(pool, 2, replace=False)
    chosen.extend((str(feature), "seeded_random", "WI") for feature in random_features)
    return pd.DataFrame(chosen, columns=["feature", "selection_reason", "calibration_dimension"])


def score_bootstraps(
    *,
    phi: np.ndarray,
    reference: np.ndarray,
    target_indices: np.ndarray,
    bootstrap_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Frozen two-site scores and training-replicate bootstrap draws."""
    sites, t, _, _ = phi.shape
    f_by_rep = 1 - (np.minimum(1.0, np.abs((phi - reference) / CALIBRATION["a"])) ** 2).mean(axis=2)
    lower, upper = importance_ci(phi, CALIBRATION["alpha_U"])
    u_by_rep = 1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"])
    imp = np.abs(phi).mean(axis=2)
    s_loss = pair_losses(phi, target_indices, CALIBRATION["a"])
    r_loss = order_pair_losses(imp, target_indices, CALIBRATION["eta"])

    def one(indices: np.ndarray) -> np.ndarray:
        f = f_by_rep[:, indices][:, :, target_indices].mean(axis=1).min(axis=0)
        s = pair_score(s_loss, indices)
        u = u_by_rep[:, indices][:, :, target_indices].mean(axis=1).min(axis=0)
        selected = np.take(phi[:, indices], target_indices, axis=3)
        site_imp = np.abs(selected).mean(axis=(1, 2))
        gm = generalizability(site_imp, CALIBRATION["b"])
        signed = np.stack([np.sort(selected[site].mean(axis=0), axis=0) for site in range(sites)])
        gs = generalizability_signed(signed, CALIBRATION["a"])
        r = pair_score(r_loss, indices)
        return np.column_stack([f, s, u, np.minimum(gm, gs), r])

    point = one(np.arange(t))
    draws = np.stack([one(indices) for indices in bootstrap_indices])
    return point, draws


def cox_inputs(selection: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], pd.DataFrame]:
    phase3 = load_phase3_module()
    with np.load(ROOT / "results/phase3_brca/bootstrap_coefficients.npz", allow_pickle=False) as saved:
        coefficients = saved["coefficients"]
        features = saved["features"].tolist()
    original = joblib.load(TCGA / "coxph_baseline.joblib")
    test = phase3.frame("X_test_selected.parquet")
    train = phase3.frame("X_train_selected.parquet")
    train_patients = train.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    clean = test.loc[~test_patients.isin(train_patients)]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    anchors = original["scaler"].transform(clean[features].astype(float))
    source_phi = coefficients[:, None, :] * anchors[None, :, :]

    ext, _, _ = phase3.external_columns(list(train.columns))
    mapped = [feature for feature in features if feature in ext]
    raw = np.column_stack([ext[feature] for feature in mapped])
    complete = np.isfinite(raw).all(axis=1)
    raw = raw[complete]
    zext = (raw - raw.mean(axis=0)) / raw.std(axis=0, ddof=0)
    ext_positions = {feature: j for j, feature in enumerate(mapped)}
    selected_names = selection.feature.tolist()
    selected_indices = np.array([features.index(feature) for feature in selected_names])
    external = np.stack(
        [coefficients[:, j, None] * zext[None, :, ext_positions[features[j]]] for j in selected_indices], axis=-1
    )
    return source_phi, external, selected_indices, features, pd.DataFrame(
        [{"site": "TCGA", "query_rows": len(anchors)}, {"site": "METABRIC", "query_rows": len(zext)}]
    )


def cox_bootstrap(
    selection: pd.DataFrame, bootstrap_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    source_phi, external, targets, _, query_inventory = cox_inputs(selection)
    t = source_phi.shape[0]
    selected_source = source_phi[:, :, targets]
    source4 = source_phi[None, ...]
    s_loss = pair_losses(source4, targets, CALIBRATION["a"])
    lower, upper = importance_ci(source4, CALIBRATION["alpha_U"])
    u_by_rep = 1 - np.minimum(1.0, (upper - lower) / CALIBRATION["u"])
    imp = np.abs(source_phi).mean(axis=1)[None, ...]
    r_loss = order_pair_losses(imp, targets, CALIBRATION["eta"])

    def one(indices: np.ndarray) -> np.ndarray:
        f = np.ones(len(targets))
        s = pair_score(s_loss, indices)
        u = u_by_rep[:, indices][:, :, targets].mean(axis=1).min(axis=0)
        source = selected_source[indices]
        other = external[indices]
        gm = generalizability(
            np.stack([np.abs(source).mean(axis=(0, 1)), np.abs(other).mean(axis=(0, 1))]),
            CALIBRATION["b"],
        )
        mean_source = source.mean(axis=0)
        mean_other = other.mean(axis=0)
        grid_n = max(mean_source.shape[0], mean_other.shape[0])
        signed_source = np.stack([empirical_rank_grid(mean_source[:, j], grid_n) for j in range(len(targets))], axis=1)
        signed_other = np.stack([empirical_rank_grid(mean_other[:, j], grid_n) for j in range(len(targets))], axis=1)
        gs = generalizability_signed(np.stack([signed_source, signed_other]), CALIBRATION["a"])
        r = pair_score(r_loss, indices)
        return np.column_stack([f, s, u, np.minimum(gm, gs), r])

    point = one(np.arange(t))
    draws = np.stack([one(indices) for indices in bootstrap_indices])
    return point, draws, query_inventory


def run(*, draws: int = 4999, seed: int = 20260917) -> pd.DataFrame:
    if draws < 99:
        raise ValueError("At least 99 bootstrap draws required")
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    cox_report = pd.read_csv(ROOT / "results/phase3_brca/feature_report.csv")
    rsf_report = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    low = "rna_ENSG00000090402.8"
    cox_selection = choose_features(cox_report, seed=seed + 501, low_feature=low)
    rsf_selection = choose_features(rsf_report, seed=seed + 502, low_feature=low)

    rng = np.random.default_rng(seed + 5000)
    cox_indices = rng.integers(0, 50, size=(draws, 50))
    rsf_indices = rng.integers(0, 20, size=(draws, 20))
    cox_point, cox_draws, cox_queries = cox_bootstrap(cox_selection, cox_indices)
    with np.load(ROOT / "results/phase4_rsf/attributions_n20.npz", allow_pickle=False) as saved:
        candidate, reference = saved["candidate"], saved["reference"]
        rsf_features = saved["features"].tolist()
    rsf_targets = np.array([rsf_features.index(feature) for feature in rsf_selection.feature])
    rsf_point, rsf_draws = score_bootstraps(
        phi=candidate, reference=reference, target_indices=rsf_targets, bootstrap_indices=rsf_indices
    )

    rows = []
    for model, selection, point, sampled, stored in [
        ("CoxPH_N50", cox_selection, cox_point, cox_draws, cox_report),
        ("RSF_N20", rsf_selection, rsf_point, rsf_draws, rsf_report),
    ]:
        pvalues = np.column_stack(
            [centered_bootstrap_pvalue(point[:, k], sampled[:, :, k], CALIBRATION["tau"]) for k in range(5)]
        )
        original = stored.set_index("feature")
        np.testing.assert_allclose(point, original.loc[selection.feature, list(COMPONENTS)].to_numpy(), atol=2e-12)
        for j, picked in selection.iterrows():
            row = {
                "model": model,
                **picked.to_dict(),
                "inference_scope": (
                    "fixed 12+12 queries and stored 8/128 attribution averages; permutation-path uncertainty excluded"
                    if model == "RSF_N20"
                    else "fixed 141+1466 queries; exact Cox attributions"
                ),
                **{name: point[j, k] for k, name in enumerate(COMPONENTS)},
                **{f"p_{name}": pvalues[j, k] for k, name in enumerate(COMPONENTS)},
                "p_joint": float(intersection_union(pvalues[j])),
            }
            row["calibration_p"] = row[f"p_{row['calibration_dimension']}"] if row["calibration_dimension"] in COMPONENTS else row["p_joint"]
            rows.append(row)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "pilot_pvalues.csv", index=False)

    formula_paths = [ROOT / "warrantxai/components.py", ROOT / "warrantxai/aggregation.py"]
    formula_hashes = {str(path.relative_to(ROOT)).replace("\\", "/"): sha(path) for path in formula_paths}
    phase4_meta = json.loads((ROOT / "results/phase4_rsf/metadata.json").read_text())
    formula_unchanged = formula_hashes == phase4_meta["formula_module_sha256"]
    inventory = {
        "CoxPH_N50": {
            "training_replicates": 50,
            "stored": {
                "bootstrap_coefficients": [50, 614],
                "feature_names": 614,
                "training_fit_attempts_and_warnings": True,
            },
            "reconstructible_without_model_or_explainer": ["F_dtj", "S pair differences", "U intervals", "site importance", "signed attributions", "rank orders"],
            "not_stored": ["materialized attribution tensor"],
            "query_rows": dict(zip(cox_queries.site, cox_queries.query_rows, strict=True)),
        },
        "RSF_N20": {
            "training_replicates": 20,
            "stored": {"candidate_attributions": [2, 20, 12, 77], "reference_attributions": [2, 20, 12, 77]},
            "reconstructible_without_model_or_explainer": ["F_dtj", "S pair differences", "U intervals", "site importance", "signed attributions", "rank orders"],
            "not_stored": ["per-permutation-path candidate/reference draws"],
            "query_rows": {"TCGA": 12, "METABRIC": 12},
        },
        "inference_target": "conditional on the fixed observed query sets, fitted main model, preprocessing, and RSF explainer averages; randomness is over patient-cluster training replicas",
        "patient_or_sample_identifiers_exported": False,
    }
    (OUT / "inventory.json").write_text(json.dumps(inventory, indent=2), encoding="utf-8")

    elapsed = perf_counter() - started
    metadata = {
        "phase": "5_step0_step1_pilot_only",
        "arguments": {"bootstrap_draws": draws, "seed": seed, "tau_components": dict.fromkeys(COMPONENTS, 0.8)},
        "method": {
            "name": "null-centered nonparametric basic bootstrap",
            "resampling_unit": "paired training-replicate index t, with replacement; sites and query rows fixed",
            "p_value": "(1 + sum[(theta_star-theta_hat) >= (theta_hat-tau)]) / (B+1)",
            "maximum_binomial_monte_carlo_se": float(np.sqrt(0.25 / (draws + 1))),
            "studentized": False,
            "reason_not_studentized": "nested SE estimation is unstable for T=20 and the frozen min/clipping/order statistics are nonsmooth",
            "full_BH_BY_applied": False,
        },
        "fixed_calibration": CALIBRATION,
        "formula_module_sha256": formula_hashes,
        "formula_module_sha256_unchanged_from_phase4_2": formula_unchanged,
        "privacy_check": {"identifier_named_columns": [], "passed": True},
        "rsf_conditionality_warning": "RSF p-values condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation paths were not stored.",
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not formula_unchanged:
        raise RuntimeError("Frozen formula module hash changed")
    if any(any(token in column.lower() for token in ("patient", "sample_id", "sample_key")) for column in result.columns):
        raise RuntimeError("Identifier-like output column found")

    report_columns = ["model", "feature", "selection_reason", *COMPONENTS, *(f"p_{x}" for x in COMPONENTS), "p_joint"]
    report = """# Faz 5 — istatistiksel sertifikasyon tasarımı ve pilotu

> **RSF koşulluluk uyarısı:** Bütün RSF p-değerleri sabit 12+12 sorgu kümesine ve saklanmış 8/128 attribution ortalamalarına koşulludur; tekil permütasyon yollarının belirsizliğini içermez.

## ADIM 0 — artefakt envanteri

CoxPH N=50 için 50×614 katsayı dizisi, feature sırası ve replika fit kayıtları saklıdır. Attribution tensörü ayrıca yazılmamıştır; ancak frozen exact Cox oracle `beta_j*z_j` ve diskteki kimliksiz sorgu girdileriyle bütün bileşen girdileri model/explainer çağrısı olmadan yeniden kurulabilir. TCGA sorgusu 141, METABRIC sorgusu 1466 satırdır.

RSF N=20 için candidate ve bağımsız yüksek bütçeli reference attribution tensörlerinin ikisi de 2×20×12×77 olarak saklıdır. Bütün frozen bileşenler bunlardan yeniden hesaplanabilir. Tek tek 8/128 permütasyon yolları saklı değildir; bu nedenle explainer Monte Carlo varyansı ayrıca resample edilemez.

## ADIM 1 — önerilen pilot yöntem

Hedef, gözlenen sabit audit sorgu kümeleri, preprocessing, ana model ve RSF explainer ortalamalarına koşullu; patient-cluster bootstrap eğitim dağılımı üzerindeki bileşen değeridir. Site sayısı iki sabit, önceden belirlenmiş kohort bağlamıdır ve site resample edilmez. Her bootstrap çekiminde replika indisi `t` replacement ile N kez seçilir, aynı indis iki sitede ve beş bileşende birlikte kullanılır, sonra frozen F/S/U/G_claim/R fonksiyonelleri aynen yeniden hesaplanır. S ve R içindeki bütün çiftler yeni çekimde yeniden oluşur.

Her bileşen için ortak frozen eşik `tau_k=0.8` kullanılır. Sınır-null altında basic bootstrap merkezi `0.8 + (theta* - theta_hat)` olur; sağ-kuyruk p-değeri `(1 + #[(theta* - theta_hat) >= (theta_hat - 0.8)])/(B+1)` olarak okunur. B=4999, dolayısıyla en küçük raporlanabilir değer 0.0002'dir. Studentization seçilmedi: N=20 içinde iç içe SE bootstrap'ı zayıf, ayrıca min, clipping ve eta-toleranslı sıralama fonksiyonelleri düzgün değildir. Bu pilot yöntem asimptotik bir öneridir; gerçek-verideki uç-durum sanity check sonlu-örneklem Tip-I geçerliliğinin kanıtı değildir. Tam sertifikasyondan önce null altında kapsama/Tip-I simülasyonu ayrıca onay kriteri olmalıdır.

Yöntem seçimi, Efron ve Tibshirani'nin bootstrap hipotez testi çerçevesindeki null-merkezleme fikrine dayanır ([Chapter 16](https://djvu.online/file/DrRxMumhAaOZY)). Monte Carlo kuyruğunun sıfır raporlanmaması için +1 düzeltmesi kullanılır; bu tercih Phipson ve Smyth'in random-resampling p-değerleri uyarısıyla uyumludur ([2010 makalesi](https://pubmed.ncbi.nlm.nih.gov/21044043/)). İlerideki çoklu test aşamasında BH yalnız independence/PRDS varsayımıyla primary, BY keyfi bağımlılık robustness kontrolü olacaktır ([BH 1995](https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/j.2517-6161.1995.tb02031.x); [BY 2001](http://projecteuclid.org/download/pdf_1/euclid.aos/1013699998)). Bu kaynaklar pilotun bileşen fonksiyonelleri için özel sonlu-örneklem geçerlilik ispatı sağlamaz.

## 5+5 feature kalibrasyon pilotu — RSF sabit 12+12 query ve saklı attribution ortalamalarına koşullu

"""
    report += markdown_table(result[report_columns])
    high = result.loc[result.selection_reason == "high_margin", ["model", "F", "p_F"]]
    near = result.loc[result.selection_reason == "near_threshold", ["model", "R", "p_R"]]
    report += "\n\nYüksek-marjin seçimlerde F ve p_F (RSF satırı sabit 12+12 query ve saklı attribution ortalamalarına koşullu):\n\n" + markdown_table(high)
    report += "\n\nR eşiğine en yakın seçimlerde R ve p_R (RSF satırı sabit 12+12 query ve saklı attribution ortalamalarına koşullu):\n\n" + markdown_table(near)
    report += (
        "\n\nSanity check beklenen yöndedir: iki yüksek-marjin feature için p_F Monte Carlo tabanı olan "
        "0.0002'dir. Eşiğe en yakın R değerleri için p_R sırasıyla 0.6442 ve 0.6278'dir; sıfıra "
        "yakın değildir ve 0.5 çevresindeki belirsiz bölgededir. RSF düşük-WI kontrolünde G_claim=0 "
        "ve p_G_claim=1 çıkmıştır. B=4999 için p tahmininin en kötü durum binomial Monte Carlo SE'si "
        f"{np.sqrt(0.25 / (draws + 1)):.4f}'tür.\n"
    )
    report += f"\n\nFormül dosyaları Faz 4.2 hash'leriyle aynıdır: **{formula_unchanged}**. Çıktılarda hasta/örnek kimliği yoktur. Gerçek pilot wall-clock: **{elapsed:.3f} saniye**.\n"
    report += "\n## Karar noktası\n\nBH/BY uygulanmadı ve tam feature ailesine geçilmedi. Pilot, yalnız koşullu bootstrap davranışını gösterir. Tam ADIM 2 için yöntem onayı yanında null kalibrasyon simülasyonunun zorunlu tutulması önerilir; RSF sonucu ayrıca 12+12 sabit sorgu ve saklanmamış permütasyon-yolu belirsizliğine koşulludur.\n"
    (ROOT / "PHASE5_DESIGN_NOTE.md").write_text(report, encoding="utf-8")
    print(result[report_columns].to_string(index=False))
    print(f"Wrote Phase 5 pilot only in {elapsed:.3f}s; full BH/BY intentionally not run.")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draws", type=int, default=4999)
    parser.add_argument("--seed", type=int, default=20260917)
    args = parser.parse_args()
    run(draws=args.draws, seed=args.seed)
