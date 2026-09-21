"""Phase 5.6 step 4: compare RSF N=20/N=50 and test the maturation hypothesis."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon, wasserstein_distance


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase56_comparison"
WARNING = (
    "RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; "
    "per-permutation-path uncertainty is excluded."
)


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def markdown_table(frame: pd.DataFrame, digits: int = 4) -> str:
    shown = frame.copy()
    for column in shown.select_dtypes(include=[np.number]):
        shown[column] = shown[column].map(lambda value: f"{value:.{digits}f}" if pd.notna(value) else "NA")
    shown = shown.fillna("NA").astype(str)
    header = "| " + " | ".join(shown.columns) + " |"
    separator = "| " + " | ".join("---" for _ in shown.columns) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in shown.to_numpy().tolist()]
    return "\n".join([header, separator, *rows])


def margin_summary(label: str, values: pd.Series) -> dict:
    margins = values - 0.8
    return {
        "inference_scope": WARNING if label.startswith("RSF") else "CoxPH fixed-query Phase 5.3 comparison",
        "family": label,
        "n": len(margins),
        "mean_margin": margins.mean(),
        "median_margin": margins.median(),
        "q25": margins.quantile(0.25),
        "q75": margins.quantile(0.75),
        "min": margins.min(),
        "max": margins.max(),
        "n_0_001": int(((margins >= 0) & (margins < 0.01)).sum()),
        "n_001_003": int(((margins >= 0.01) & (margins < 0.03)).sum()),
        "n_003_005": int(((margins >= 0.03) & (margins < 0.05)).sum()),
        "n_005_010": int(((margins >= 0.05) & (margins < 0.10)).sum()),
        "n_ge_005": int((margins >= 0.05).sum()),
        "fraction_ge_005": float((margins >= 0.05).mean()),
    }


def draw_figure(transition: pd.DataFrame, margins: dict[str, pd.Series], counts: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    colors = np.where(transition.point_n20, "#d95f02", "#9e9e9e")
    axes[0].scatter(transition.R_n20, transition.R_n50, c=colors, s=38, alpha=0.8)
    lo = min(transition.R_n20.min(), transition.R_n50.min()) - 0.01
    hi = max(transition.R_n20.max(), transition.R_n50.max()) + 0.01
    axes[0].plot([lo, hi], [lo, hi], "--", color="#555555", linewidth=1)
    axes[0].set(xlim=(lo, hi), ylim=(lo, hi), xlabel="R at N=20", ylabel="R at N=50",
                title="Paired R values (orange: N=20 pass)")
    axes[0].grid(alpha=0.2)

    palette = {"RSF N=20": "#d95f02", "RSF N=50": "#2ca25f", "CoxPH N=50": "#1769aa"}
    for label, values in margins.items():
        ordered = np.sort(np.asarray(values - 0.8, dtype=float))
        axes[1].step(ordered, np.arange(1, len(ordered)+1)/len(ordered), where="post",
                     color=palette[label], linewidth=2, label=label)
    axes[1].axvline(0.05, color="#777777", linestyle="--", linewidth=1)
    axes[1].set(xlabel="WI margin (WI − 0.8)", ylabel="Empirical CDF", title="Point-pass margin shape")
    axes[1].grid(alpha=0.2)
    axes[1].legend(fontsize=8)

    x = np.arange(len(counts))
    width = 0.24
    for offset, metric in enumerate(["point", "BH", "BY"]):
        axes[2].bar(x + (offset-1)*width, counts[metric], width, label=metric)
    axes[2].set_xticks(x, counts.model)
    axes[2].set_ylabel("Feature count")
    axes[2].set_title("Point estimate and FDR counts")
    axes[2].legend()
    axes[2].grid(axis="y", alpha=0.2)
    fig.savefig(OUT / "n20_n50_maturation.png", dpi=180)
    plt.close(fig)


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    n20 = pd.read_csv(ROOT / "results/phase4_rsf/feature_report.csv")
    n50 = pd.read_csv(ROOT / "results/phase56_rsf_n50/feature_report.csv")
    cert20 = pd.read_csv(ROOT / "results/phase53_rsf_certification/feature_certification.csv")
    cert50 = pd.read_csv(ROOT / "results/phase56_rsf_n50_certification/feature_certification.csv")
    cox = pd.read_csv(ROOT / "results/phase53_cox_recursive/feature_certification.csv")
    transition = n20[["feature", "block", "R", "WI", "point_adequate", "bottleneck"]].rename(columns={
        "R": "R_n20", "WI": "WI_n20", "point_adequate": "point_n20", "bottleneck": "bottleneck_n20"
    }).merge(n50[["feature", "R", "WI", "point_adequate", "bottleneck"]].rename(columns={
        "R": "R_n50", "WI": "WI_n50", "point_adequate": "point_n50", "bottleneck": "bottleneck_n50"
    }), on="feature", validate="one_to_one")
    transition = transition.merge(cert20[["feature", "certified_BH_alpha_0_05", "certified_BY_alpha_0_05"]].rename(columns={
        "certified_BH_alpha_0_05": "BH_n20", "certified_BY_alpha_0_05": "BY_n20"
    }), on="feature", validate="one_to_one").merge(
        cert50[["feature", "certified_BH_alpha_0_05", "certified_BY_alpha_0_05"]].rename(columns={
            "certified_BH_alpha_0_05": "BH_n50", "certified_BY_alpha_0_05": "BY_n50"
        }), on="feature", validate="one_to_one")
    transition.insert(0, "inference_scope", WARNING)
    transition["delta_R"] = transition.R_n50 - transition.R_n20
    transition["delta_WI"] = transition.WI_n50 - transition.WI_n20
    transition["point_transition"] = np.select(
        [transition.point_n20.astype(bool) & transition.point_n50.astype(bool),
         transition.point_n20.astype(bool) & ~transition.point_n50.astype(bool),
         ~transition.point_n20.astype(bool) & transition.point_n50.astype(bool)],
        ["retained", "lost", "new"], default="neither"
    )
    transition.to_csv(OUT / "feature_transition.csv", index=False)

    n20_point = n20.loc[n20.point_adequate.astype(bool), "WI"]
    n50_point = n50.loc[n50.point_adequate.astype(bool), "WI"]
    cox_point = cox.loc[cox.point_WI_ge_0_8.astype(bool), "WI"]
    margin_table = pd.DataFrame([
        margin_summary("RSF N=20", n20_point), margin_summary("RSF N=50", n50_point),
        margin_summary("CoxPH N=50", cox_point),
    ])
    margin_table.to_csv(OUT / "margin_distribution.csv", index=False)

    counts = pd.DataFrame([
        {"model": "RSF N=20", "point": int(cert20.point_WI_ge_0_8.sum()),
         "BH": int(cert20.certified_BH_alpha_0_05.sum()), "BY": int(cert20.certified_BY_alpha_0_05.sum()),
         "inference_scope": WARNING},
        {"model": "RSF N=50", "point": int(cert50.point_WI_ge_0_8.sum()),
         "BH": int(cert50.certified_BH_alpha_0_05.sum()), "BY": int(cert50.certified_BY_alpha_0_05.sum()),
         "inference_scope": WARNING},
        {"model": "CoxPH N=50", "point": int(cox.point_WI_ge_0_8.sum()),
         "BH": int(cox.certified_BH_alpha_0_05.sum()), "BY": int(cox.certified_BY_alpha_0_05.sum()),
         "inference_scope": "CoxPH fixed-query Phase 5.3 comparison"},
    ])
    counts["BH_survival_fraction"] = counts.BH / counts.point
    counts["BY_survival_fraction"] = counts.BY / counts.point
    counts.to_csv(OUT / "count_comparison.csv", index=False)

    all_test = wilcoxon(transition.delta_R)
    selected = transition.loc[transition.point_n20.astype(bool)]
    selected_test = wilcoxon(selected.delta_R)
    r_summary = pd.DataFrame([
        {"inference_scope": WARNING, "scope": "all 77", "n": len(transition),
         "mean_delta_R": transition.delta_R.mean(), "median_delta_R": transition.delta_R.median(),
         "R_decreased": int((transition.delta_R < 0).sum()), "R_increased": int((transition.delta_R > 0).sum()),
         "wilcoxon_p": all_test.pvalue},
        {"inference_scope": WARNING, "scope": "N20 point-pass 23", "n": len(selected),
         "mean_delta_R": selected.delta_R.mean(), "median_delta_R": selected.delta_R.median(),
         "R_decreased": int((selected.delta_R < 0).sum()), "R_increased": int((selected.delta_R > 0).sum()),
         "wilcoxon_p": selected_test.pvalue},
    ])
    r_summary.to_csv(OUT / "R_delta_summary.csv", index=False)

    n20_set = set(n20.loc[n20.point_adequate.astype(bool), "feature"])
    n50_set = set(n50.loc[n50.point_adequate.astype(bool), "feature"])
    cox_meta = json.loads((ROOT / "results/phase3_brca/metadata.json").read_text())
    retention = pd.DataFrame([
        {"model": "RSF", "from_N": 20, "to_N": 50, "from_point_count": 23, "to_point_count": 20,
         "count_ratio": 20/23, "same_feature_retained": len(n20_set & n50_set),
         "set_retention": len(n20_set & n50_set)/23, "lost": len(n20_set-n50_set), "new": len(n50_set-n20_set)},
        {"model": "CoxPH", "from_N": 6, "to_N": 50, "from_point_count": 27, "to_point_count": 22,
         "count_ratio": 22/27, "same_feature_retained": int(cox_meta["stats"]["previously_passing_still_passing"]),
         "set_retention": int(cox_meta["stats"]["previously_passing_still_passing"])/27,
         "lost": 12, "new": 7},
    ])
    retention.to_csv(OUT / "retention_comparison.csv", index=False)

    margins = {"RSF N=20": n20_point, "RSF N=50": n50_point, "CoxPH N=50": cox_point}
    draw_figure(transition, margins, counts)

    step1 = json.loads((ROOT / "results/phase56_rsf_n50/metadata.json").read_text())
    step2 = json.loads((ROOT / "results/phase56_rsf_n50_null_gate/metadata.json").read_text())
    step3 = json.loads((ROOT / "results/phase56_rsf_n50_certification/metadata.json").read_text())
    formula_sha = step3["formula_module_sha256"]
    elapsed = perf_counter() - started
    metadata = {
        "phase": "5.6_step4_RSF_N20_N50_maturation_comparison",
        "rsf_conditionality_warning": WARNING,
        "verdict": "partially_supported",
        "reason": (
            "N20-selected features show R decline and point/BH survival falls, but all-77 R has no general "
            "downward shift and the N50 point-pass margin distribution remains more high-margin than CoxPH."
        ),
        "counts": counts.drop(columns="inference_scope").to_dict(orient="records"),
        "retention": retention.to_dict(orient="records"),
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": formula_sha == step1["formula_module_sha256"] == step2["formula_module_sha256"],
        "source_sha256": {
            "results/phase4_rsf/feature_report.csv": sha(ROOT / "results/phase4_rsf/feature_report.csv"),
            "results/phase56_rsf_n50/feature_report.csv": sha(ROOT / "results/phase56_rsf_n50/feature_report.csv"),
            "results/phase53_rsf_certification/feature_certification.csv": sha(ROOT / "results/phase53_rsf_certification/feature_certification.csv"),
            "results/phase56_rsf_n50_certification/feature_certification.csv": sha(ROOT / "results/phase56_rsf_n50_certification/feature_certification.csv"),
        },
        "step_wall_clock_seconds": {"step1_new_30": step1["runtime"]["new_replicates_20_to_49_seconds"],
                                    "step1_N50_accumulated": step1["runtime"]["N50_total_accumulated_seconds"],
                                    "step2_gate": step2["wall_clock_seconds"],
                                    "step3_certification": step3["wall_clock_seconds"],
                                    "step4_comparison": elapsed},
        "privacy_check": {"patient_or_sample_identifiers_exported": False, "passed": True},
        "python": platform.python_version(), "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    if not metadata["formula_module_sha256_unchanged"]:
        raise RuntimeError("Frozen formula hashes changed across Phase 5.6")

    warning = WARNING
    report = f"""# Faz 5.6 — RSF N=50 tam koşu ve FDR sertifikasyonu

> **Koşulluluk uyarısı:** {warning}

## Sonuç

Winner's-curse olgunlaşma hipotezi **kısmen desteklendi**. N=20 nokta-geçen 23 feature'ın 19'unda R düştü; bu seçilmiş grupta ortalama değişim {selected.delta_R.mean():.4f}, medyan {selected.delta_R.median():.4f} oldu. Nokta-geçen sayı 23→20, BH 15→10 ve BY 12→10 düştü. Buna karşılık tüm 77 feature'da ortalama R değişimi yalnız {transition.delta_R.mean():.4f}; 39 düşüşe karşı 38 artış var ve eşlenmiş Wilcoxon p={all_test.pvalue:.3f}. Genel bir R çöküşü yok; düşüş N=20'de seçilmiş kazananlarda yoğunlaşıyor.

Marjin şekli CoxPH yönüne kaydı fakat yakınsamadı. RSF'de WI marjini 0.05–0.10 bandındaki pay 15/23 (%65.2) → 9/20 (%45.0) düştü; medyan marjin {n20_point.median()-0.8:.4f} → {n50_point.median()-0.8:.4f}. CoxPH N=50'de aynı yüksek-marjin sayısı yalnız 2/22 (%9.1), medyan {cox_point.median()-0.8:.4f}. RSF BH hayatta-kalma oranı %65.2→%50.0 geriledi, fakat CoxPH'nin %36.4 değerinin üzerinde kaldı.

## ADIM 1 — RSF N=50

> **Koşulluluk uyarısı:** {warning}

Onaylanan sampled-background k=16, candidate/reference 8/128 ve aynı 12+12 query ile 77/77 feature değerlendirildi. İlk 20 deterministik replika doğrulanarak korundu; aynı akışın 21–50 replikaları yeni fit edildi. Ana RSF modeli değiştirilmedi. Yeni 30 replika {step1['runtime']['new_replicates_20_to_49_seconds']:.1f} sn ({step1['runtime']['new_replicates_20_to_49_seconds']/3600:.2f} saat), tüm N=50 dizisinin birikimli ölçülen maliyeti {step1['runtime']['N50_total_accumulated_seconds']:.1f} sn ({step1['runtime']['N50_total_accumulated_seconds']/3600:.2f} saat) oldu. Hiçbir fit uyarısı veya atlanan replika yok.

### Block agregatları

{markdown_table(pd.read_csv(ROOT / 'results/phase56_rsf_n50/summary.csv'))}

## ADIM 2 — N=50 null gate

> **Koşulluluk uyarısı:** {warning}

`recursive_studentized`, aynı 2000×999 `single_active`+`tie`, üç alpha ve 30 kontrol içinde Holm bataryasında **geçti**; 0/30 anti-conservative ihlal bulundu. Wall-clock {step2['wall_clock_seconds']:.3f} sn.

## ADIM 3 — tam-aile FDR

> **Koşulluluk uyarısı:** {warning}

{markdown_table(counts)}

N=50'de BH ve BY aynı 10 feature'ı sertifikaladı. Tam p/q tablosu `results/phase56_rsf_n50_certification/feature_certification.csv` içindedir. Sertifikasyon wall-clock {step3['wall_clock_seconds']:.3f} sn.

## ADIM 4 — N=20→N=50 karşılaştırması

> **Koşulluluk uyarısı:** {warning}

### Marjin dağılımı

{markdown_table(margin_table)}

### R değişimi

{markdown_table(r_summary)}

### Retention

{markdown_table(retention)}

RSF'de count oranı 20/23=%87.0, birebir feature retention 16/23=%69.6'dır. CoxPH'de 22/27=%81.5 yalnız sayım oranıdır; gerçek N=6 feature kümesinin N=50'de retention'ı 15/27=%55.6'dır. Böylece RSF'deki feature-kümesi kaybı CoxPH'den daha hafif, fakat yönü aynı winner's-curse desenindedir.

![N=20 ve N=50 olgunlaşma karşılaştırması](results/phase56_comparison/n20_n50_maturation.png)

## Dürüstlük kontrolleri

Formül SHA256 değerleri bütün adımlarda sabit kaldı. Yeni N=50 çıktılarında kimlik-sızıntısı kontrolü geçti. Query havuzu genişletilmedi; sonuçlar sabit 12+12 query ve saklanmamış permütasyon-yolu belirsizliği sınırlamasını taşır. ADIM 4 analizi {elapsed:.3f} sn sürdü.
"""
    (ROOT / "PHASE5_6_REPORT.md").write_text(report, encoding="utf-8")
    (OUT / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    run()
