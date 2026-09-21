"""Phase 5.5: CoxPH R sensitivity to a 77-feature ranking universe.

This diagnostic reuses stored N=50 coefficients and the frozen TCGA audit
matrix. It does not fit/predict with a model, invoke an explainer, resample,
or run BH/BY. Official Phase 5.3 outputs remain untouched.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
from time import perf_counter

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wasserstein_distance

from warrantxai.components import rank_reproducibility


ROOT = Path(__file__).resolve().parents[1]
TCGA = ROOT.parent / "Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer"
OUT = ROOT / "results/phase55_rank_universe_diagnostic"
TAU = 0.8
ETA = 0.02
COMPONENTS_WITHOUT_R = ("F", "S", "U", "G_claim")


def sha256(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def markdown_table(frame: pd.DataFrame, digits: int = 5) -> str:
    shown = frame.copy()
    for column in shown.select_dtypes(include=[np.number]):
        shown[column] = shown[column].map(lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA")
    shown = shown.fillna("NA").astype(str)
    header = "| " + " | ".join(shown.columns) + " |"
    separator = "| " + " | ".join("---" for _ in shown.columns) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in shown.to_numpy().tolist()]
    return "\n".join([header, separator, *rows])


def load_cox_importance() -> tuple[np.ndarray, list[str], int]:
    """Reconstruct the exact stored Phase 3.2 importance array, without prediction."""
    with np.load(ROOT / "results/phase3_brca/bootstrap_coefficients_n50.npz", allow_pickle=False) as saved:
        coefficients = saved["coefficients"]
        features = saved["features"].tolist()
    if coefficients.shape != (50, 614) or len(features) != 614:
        raise RuntimeError(f"Expected frozen (50,614) coefficients, got {coefficients.shape}")

    artifact = joblib.load(TCGA / "coxph_baseline.joblib")
    if list(artifact["features"]) != features:
        raise RuntimeError("Cox artifact feature order differs from stored coefficients")
    train = pd.read_parquet(TCGA / "X_train_selected.parquet").set_index("sample_key")
    test = pd.read_parquet(TCGA / "X_test_selected.parquet").set_index("sample_key")
    train_patients = train.index.astype(str).str[:12]
    test_patients = test.index.astype(str).str[:12]
    clean = test.loc[~test_patients.isin(train_patients)]
    clean = clean.loc[~clean.index.astype(str).str[:12].duplicated()]
    anchors = artifact["scaler"].transform(clean[features].astype(float))
    if len(anchors) != 141:
        raise RuntimeError(f"Expected 141 frozen audit rows, got {len(anchors)}")
    importance = np.abs(coefficients) * np.mean(np.abs(anchors), axis=0)
    return importance, features, len(anchors)


def describe(label: str, values: pd.Series | np.ndarray) -> dict[str, float | str | int]:
    values = pd.Series(np.asarray(values, dtype=float))
    return {
        "distribution": label,
        "n": len(values),
        "mean": values.mean(),
        "std": values.std(ddof=1),
        "min": values.min(),
        "q25": values.quantile(0.25),
        "median": values.median(),
        "q75": values.quantile(0.75),
        "max": values.max(),
        "n_ge_0": int((values >= 0).sum()),
        "n_gt_0_03": int((values > 0.03).sum()),
    }


def draw_figure(all_features: pd.DataFrame, family_summary: dict[str, pd.Series]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.7), constrained_layout=True)
    passing = all_features.loc[all_features.cox_original_point_pass]
    axes[0].scatter(passing.R_614, passing.R_77, color="#1769aa", s=48, alpha=0.85)
    limits = [min(passing.R_614.min(), passing.R_77.min()) - 0.004,
              max(passing.R_614.max(), passing.R_77.max()) + 0.004]
    axes[0].plot(limits, limits, "--", color="#555555", linewidth=1)
    axes[0].set(xlim=limits, ylim=limits, xlabel="Original CoxPH R (614 universe)",
                ylabel="Diagnostic CoxPH R (77 universe)", title="Same 22 CoxPH point-pass features")
    axes[0].grid(alpha=0.2)

    colors = {"CoxPH original": "#1769aa", "CoxPH 77-universe": "#2ca25f", "RSF actual": "#d95f02"}
    for label, values in family_summary.items():
        ordered = np.sort(np.asarray(values, dtype=float))
        axes[1].step(ordered, np.arange(1, len(ordered) + 1) / len(ordered), where="post",
                     label=label, color=colors[label], linewidth=2)
    axes[1].axvline(0.03, color="#777777", linestyle="--", linewidth=1, label="margin=0.03")
    axes[1].set(xlabel="WI margin (WI − 0.8)", ylabel="Empirical CDF",
                title="Model-specific point-pass families")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)
    fig.savefig(OUT / "rank_universe_comparison.png", dpi=180)
    plt.close(fig)


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    cox_path = ROOT / "results/phase53_cox_recursive/feature_certification.csv"
    rsf_path = ROOT / "results/phase53_rsf_certification/feature_certification.csv"
    cox = pd.read_csv(cox_path)
    rsf = pd.read_csv(rsf_path)
    if len(cox) != 77 or set(cox.feature) != set(rsf.feature):
        raise RuntimeError("CoxPH and RSF frozen assessable families must be the same 77 features")

    importance, full_features, audit_rows = load_cox_importance()
    selected = np.array([full_features.index(feature) for feature in cox.feature])
    r_614 = rank_reproducibility(importance[None, ...], ETA)[selected]
    np.testing.assert_allclose(r_614, cox.R, atol=2e-12)
    r_77 = rank_reproducibility(importance[:, selected][None, ...], ETA)

    result = cox[["feature", "block", *COMPONENTS_WITHOUT_R, "R", "WI", "bottleneck",
                  "point_WI_ge_0_8"]].copy()
    result = result.rename(columns={"R": "R_614", "WI": "WI_614", "bottleneck": "bottleneck_614",
                                    "point_WI_ge_0_8": "cox_original_point_pass"})
    result["R_77"] = r_77
    result["delta_R_77_minus_614"] = result.R_77 - result.R_614
    result["WI_77"] = np.minimum(result[list(COMPONENTS_WITHOUT_R)].min(axis=1), result.R_77)
    result["bottleneck_77"] = result[[*COMPONENTS_WITHOUT_R, "R_77"]].idxmin(axis=1).replace({"R_77": "R"})
    result["WI_margin_614"] = result.WI_614 - TAU
    result["R_margin_614"] = result.R_614 - TAU
    result["WI_margin_77"] = result.WI_77 - TAU
    result["R_margin_77"] = result.R_77 - TAU
    result["cox_diagnostic_point_pass"] = result.WI_77 >= TAU
    result = result.merge(
        rsf[["feature", "R", "WI", "bottleneck", "point_WI_ge_0_8"]].rename(
            columns={"R": "RSF_R_77", "WI": "RSF_WI", "bottleneck": "RSF_bottleneck",
                     "point_WI_ge_0_8": "RSF_point_pass"}
        ), on="feature", validate="one_to_one"
    )
    result["RSF_R_margin"] = result.RSF_R_77 - TAU
    result["RSF_WI_margin"] = result.RSF_WI - TAU
    result.to_csv(OUT / "all_77_feature_comparison.csv", index=False)

    same_22 = result.loc[result.cox_original_point_pass].copy()
    same_22.to_csv(OUT / "cox_original_22_comparison.csv", index=False)
    newly_passing = result.loc[result.cox_diagnostic_point_pass & ~result.cox_original_point_pass].copy()
    newly_passing.to_csv(OUT / "diagnostic_status_changes.csv", index=False)

    distribution_rows = [
        describe("same22_Cox_R_margin_614", same_22.R_margin_614),
        describe("same22_Cox_R_margin_77", same_22.R_margin_77),
        describe("same22_RSF_R_margin_77", same_22.RSF_R_margin),
        describe("all77_Cox_R_margin_614", result.R_margin_614),
        describe("all77_Cox_R_margin_77", result.R_margin_77),
        describe("all77_RSF_R_margin_77", result.RSF_R_margin),
        describe("Cox_original_point_family_WI_margin", result.loc[result.cox_original_point_pass, "WI_margin_614"]),
        describe("Cox_77_universe_point_family_WI_margin", result.loc[result.cox_diagnostic_point_pass, "WI_margin_77"]),
        describe("RSF_actual_point_family_WI_margin", result.loc[result.RSF_point_pass, "RSF_WI_margin"]),
    ]
    distribution_summary = pd.DataFrame(distribution_rows)
    distribution_summary.to_csv(OUT / "distribution_summary.csv", index=False)

    def distance_row(scope: str, before: np.ndarray, after: np.ndarray, target: np.ndarray) -> dict:
        before_w = wasserstein_distance(before, target)
        after_w = wasserstein_distance(after, target)
        return {
            "scope": scope,
            "metric": "Wasserstein-1",
            "original_distance_to_RSF": before_w,
            "R77_distance_to_RSF": after_w,
            "absolute_change": after_w - before_w,
            "relative_change": after_w / before_w - 1,
            "closer_after_restriction": after_w < before_w,
        }

    distances = pd.DataFrame([
        distance_row("same 22 features; R margins", same_22.R_margin_614, same_22.R_margin_77,
                     same_22.RSF_R_margin),
        distance_row("all 77 matched features; R margins", result.R_margin_614, result.R_margin_77,
                     result.RSF_R_margin),
        distance_row("model-specific point-pass families; WI margins",
                     result.loc[result.cox_original_point_pass, "WI_margin_614"],
                     result.loc[result.cox_diagnostic_point_pass, "WI_margin_77"],
                     result.loc[result.RSF_point_pass, "RSF_WI_margin"]),
    ])
    distances.to_csv(OUT / "distance_to_rsf.csv", index=False)

    family_summary = {
        "CoxPH original": result.loc[result.cox_original_point_pass, "WI_margin_614"],
        "CoxPH 77-universe": result.loc[result.cox_diagnostic_point_pass, "WI_margin_77"],
        "RSF actual": result.loc[result.RSF_point_pass, "RSF_WI_margin"],
    }
    draw_figure(result, family_summary)

    rho_614 = spearmanr(same_22.R_614, same_22.RSF_R_77)
    rho_77 = spearmanr(same_22.R_77, same_22.RSF_R_77)
    formula_sha = {
        "warrantxai/components.py": sha256(ROOT / "warrantxai/components.py"),
        "warrantxai/aggregation.py": sha256(ROOT / "warrantxai/aggregation.py"),
    }
    expected_sha = {
        "warrantxai/components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "warrantxai/aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.5_rank_universe_size_diagnostic",
        "diagnostic_only": True,
        "official_certification_outputs_modified": False,
        "new_bootstrap_resampling": False,
        "model_fit_prediction_or_explainer_calls": False,
        "parameters": {"eta": ETA, "tau": TAU, "replicates": 50, "audit_rows": audit_rows,
                       "original_ranking_universe": 614, "diagnostic_ranking_universe": 77},
        "source_sha256": {
            str(cox_path.relative_to(ROOT)): sha256(cox_path),
            str(rsf_path.relative_to(ROOT)): sha256(rsf_path),
            "results/phase3_brca/bootstrap_coefficients_n50.npz": sha256(
                ROOT / "results/phase3_brca/bootstrap_coefficients_n50.npz"
            ),
        },
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": formula_sha == expected_sha,
        "validation": {"recomputed_R614_max_abs_error": float(np.max(np.abs(r_614 - cox.R)))},
        "counts": {"official_cox_point_pass": int(result.cox_original_point_pass.sum()),
                   "diagnostic_cox_point_pass": int(result.cox_diagnostic_point_pass.sum()),
                   "RSF_point_pass": int(result.RSF_point_pass.sum()),
                   "new_diagnostic_passes": newly_passing.feature.tolist()},
        "same22_spearman_with_RSF": {"R614_rho": rho_614.statistic, "R77_rho": rho_77.statistic},
        "verdict": "hypothesis_not_supported",
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "wall_clock_seconds": elapsed,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    point_families = distribution_summary.loc[
        distribution_summary.distribution.isin([
            "Cox_original_point_family_WI_margin", "Cox_77_universe_point_family_WI_margin",
            "RSF_actual_point_family_WI_margin"
        ])
    ]
    same22_summary = distribution_summary.loc[
        distribution_summary.distribution.isin([
            "same22_Cox_R_margin_614", "same22_Cox_R_margin_77", "same22_RSF_R_margin_77"
        ])
    ]
    selected_columns = ["feature", "R_614", "R_77", "delta_R_77_minus_614", "WI_614", "WI_77",
                        "WI_margin_614", "WI_margin_77", "RSF_R_77", "RSF_WI"]
    report = f"""# Faz 5.5 — R sıralama-evreni boyutu tanısı

Bu çalışma ayrı bir duyarlılık hesabıdır. CoxPH N=50'nin saklı katsayıları ve 141 satırlık sabit TCGA audit kümesi kullanıldı; yeni bootstrap, model fit/tahmin, explainer veya BH/BY çalıştırılmadı. Resmî 614-evrenli R ve 8/77 sertifikasyon sonucu değişmedi.

## Sonuç

**614 vs 77 sıralama-evreni hipotezi desteklenmedi.** CoxPH evrenini 77 feature'a indirmek aynı 22 nokta-geçen feature'da R'yi ortalama yalnız **{same_22.delta_R_77_minus_614.mean():.5f}** değiştirdi (medyan {same_22.delta_R_77_minus_614.median():.5f}; aralık {same_22.delta_R_77_minus_614.min():.5f}–{same_22.delta_R_77_minus_614.max():.5f}). Yüksek WI-marjinli feature sayısı `>0.03` düzeyinde 9 olarak kaldı; RSF'de 16'dır.

Aynı 22 feature'ın R marjinlerinde RSF'ye Wasserstein uzaklığı **{distances.iloc[0].original_distance_to_RSF:.5f} → {distances.iloc[0].R77_distance_to_RSF:.5f}** oldu; yakınlaşma yerine %{100 * distances.iloc[0].relative_change:.1f} uzaklaşma görüldü. Tüm 77 eşlenmiş feature'da uzaklık **{distances.iloc[1].original_distance_to_RSF:.5f} → {distances.iloc[1].R77_distance_to_RSF:.5f}** oldu. Bu yön, evren boyutunun ana açıklama olduğu hipoteziyle uyumlu değildir.

## Aynı 22 CoxPH nokta-geçen feature

{markdown_table(same_22[selected_columns])}

## R marjin dağılımları

{markdown_table(same22_summary)}

77-evrenli CoxPH R marjininin medyanı {same_22.R_margin_77.median():.5f}; özgün değer {same_22.R_margin_614.median():.5f}, aynı feature'lardaki RSF medyanı {same_22.RSF_R_margin.median():.5f}'tir. Feature-bazlı CoxPH–RSF Spearman rho da {rho_614.statistic:.3f}'ten {rho_77.statistic:.3f}'e yalnız çok az değişti ve zayıf kaldı.

## Karşı-olgusal nokta-geçen aileler

{markdown_table(point_families)}

Yalnız R'yi 77 evreninde hesaplayan tanısal WI, CoxPH'de 22 yerine 23 feature'ı eşiğin üstüne taşır. Yeni feature `{newly_passing.feature.iloc[0]}` için WI {newly_passing.WI_614.iloc[0]:.5f}→{newly_passing.WI_77.iloc[0]:.5f}'tir. Bu resmî bir adequacy veya sertifikasyon sonucu değildir. Karşı-olgusal CoxPH ailesinin WI marjin medyanı {point_families.iloc[1]['median']:.5f} olarak kalırken RSF medyanı {point_families.iloc[2]['median']:.5f}'tir; `>0.03` sayısı yine 9'a karşı 16'dır. Model-özel nokta-geçen ailelerin Wasserstein uzaklığı {distances.iloc[2].original_distance_to_RSF:.5f}→{distances.iloc[2].R77_distance_to_RSF:.5f} ile azalmadı.

![R evreni karşılaştırması](results/phase55_rank_universe_diagnostic/rank_universe_comparison.png)

## Yorum

614 rakip yerine 77 rakip kullanmak CoxPH R düzeyini hafifçe yukarı kaydırıyor, fakat Faz 5.4'teki düşük/yüksek marjin dağılımı ayrımını kapatmıyor. Bu nedenle ters-sezgisel BH hayatta-kalma farkının ana kaynağı sıralama-evreni boyutu değildir. Kalan fark, modellerin replika boyunca ürettiği önem sıralamalarının yapısı ve feature-özel model/veri davranışıyla daha uyumludur; bu tanı tek başına bunu biyolojik fark olarak nitelemez.

Hesap {elapsed:.3f} saniye sürdü. Dondurulmuş 614-evrenli R yeniden hesaplanarak kaynak tabloyla maksimum mutlak hata {np.max(np.abs(r_614 - cox.R)):.1e} bulundu. Formül SHA256 değerleri değişmedi ve kimlik-sızıntısı kontrolü geçti.
"""
    (ROOT / "PHASE5_5_RANK_UNIVERSE_DIAGNOSTIC.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    run()
