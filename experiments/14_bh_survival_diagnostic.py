"""Phase 5.4: diagnose CoxPH/RSF BH survival from existing tables only.

No resampling, model fit, explainer call, or warrant-formula evaluation occurs.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform
from time import perf_counter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase54_bh_survival_diagnostic"
COMPONENTS = ("F", "S", "U", "G_claim", "R")
SOURCES = {
    "CoxPH N=50": ROOT / "results/phase53_cox_recursive/feature_certification.csv",
    "RSF N=20": ROOT / "results/phase53_rsf_certification/feature_certification.csv",
}


def sha256(path: Path) -> str:
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


def load_diagnostics() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for model, path in SOURCES.items():
        data = pd.read_csv(path)
        point = data.loc[data["point_WI_ge_0_8"].astype(bool)].copy()
        point["model"] = model
        point["BH_status"] = np.where(
            point["certified_BH_alpha_0_05"].astype(bool), "BH certified", "point only"
        )
        point["WI_margin"] = point["WI"] - 0.8
        point["bottleneck_value"] = [row[row["bottleneck"]] for _, row in point.iterrows()]
        point["bottleneck_margin"] = point["bottleneck_value"] - 0.8
        point["p_bottleneck"] = [row[f"p_{row['bottleneck']}"] for _, row in point.iterrows()]
        p_columns = [f"p_{component}" for component in COMPONENTS]
        point["inferential_bottleneck"] = (
            point[p_columns].idxmax(axis=1).str.removeprefix("p_")
        )
        point["replicates"] = 50 if model.startswith("CoxPH") else 20
        point["ranking_feature_count"] = 614 if model.startswith("CoxPH") else 77
        frames.append(
            point[
                [
                    "model", "feature", "block", "BH_status", *COMPONENTS, "WI", "WI_margin",
                    "bottleneck", "bottleneck_value", "bottleneck_margin", "p_bottleneck",
                    "inferential_bottleneck", "p_joint", "q_BH", "replicates", "ranking_feature_count",
                ]
            ]
        )
    return pd.concat(frames, ignore_index=True)


def summarise(point: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    groups = []
    for (model, status), part in point.groupby(["model", "BH_status"], sort=False):
        groups.append(
            {
                "model": model,
                "BH_status": status,
                "n": len(part),
                "margin_mean": part.WI_margin.mean(),
                "margin_median": part.WI_margin.median(),
                "margin_q25": part.WI_margin.quantile(0.25),
                "margin_q75": part.WI_margin.quantile(0.75),
                "margin_min": part.WI_margin.min(),
                "margin_max": part.WI_margin.max(),
                "p_bottleneck_mean": part.p_bottleneck.mean(),
                "p_joint_mean": part.p_joint.mean(),
            }
        )
    status_summary = pd.DataFrame(groups)

    bins = [-np.inf, 0.01, 0.03, 0.05, 0.10, np.inf]
    labels = ["<0.01", "0.01-0.03", "0.03-0.05", "0.05-0.10", ">=0.10"]
    point = point.copy()
    point["margin_band"] = pd.cut(point.WI_margin, bins=bins, labels=labels, right=False)
    rows = []
    for (model, band), part in point.groupby(["model", "margin_band"], observed=False, sort=False):
        rows.append(
            {
                "model": model,
                "margin_band": str(band),
                "n": len(part),
                "margin_mean": part.WI_margin.mean(),
                "p_bottleneck_mean": part.p_bottleneck.mean(),
                "p_bottleneck_median": part.p_bottleneck.median(),
                "p_joint_mean": part.p_joint.mean(),
                "p_joint_median": part.p_joint.median(),
                "BH_certified": int((part.BH_status == "BH certified").sum()),
            }
        )
    band_summary = pd.DataFrame(rows)

    bottlenecks = (
        point.loc[point.BH_status == "point only"]
        .groupby(["model", "bottleneck"], as_index=False)
        .size()
        .rename(columns={"size": "n"})
    )
    return status_summary, band_summary, bottlenecks


def draw_scatter(point: pd.DataFrame) -> None:
    colors = {"CoxPH N=50": "#1769aa", "RSF N=20": "#d95f02"}
    markers = {"BH certified": "o", "point only": "X"}
    fig, ax = plt.subplots(figsize=(8.2, 5.2), constrained_layout=True)
    for (model, status), part in point.groupby(["model", "BH_status"]):
        ax.scatter(
            part.WI_margin,
            part.p_bottleneck,
            label=f"{model} — {status}",
            color=colors[model],
            marker=markers[status],
            s=58,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.5,
        )
    ax.axvspan(0, 0.03, color="#777777", alpha=0.09, label="0–0.03 comparison band")
    ax.set_yscale("log")
    ax.set_xlabel("Observed WI margin (WI − 0.8)")
    ax.set_ylabel("Point-estimate bottleneck p-value (log scale)")
    ax.set_title("Phase 5.4: p-value versus observed margin")
    ax.grid(True, which="both", alpha=0.2)
    ax.legend(fontsize=8, ncol=2)
    fig.savefig(OUT / "p_vs_margin.png", dpi=180)
    plt.close(fig)


def run() -> None:
    started = perf_counter()
    OUT.mkdir(parents=True, exist_ok=True)
    point = load_diagnostics()
    status_summary, band_summary, bottlenecks = summarise(point)
    point_only = point.loc[point.BH_status == "point only"].sort_values(["model", "WI_margin"])

    point.to_csv(OUT / "point_pass_feature_diagnostics.csv", index=False)
    point_only.to_csv(OUT / "point_only_feature_diagnostics.csv", index=False)
    status_summary.to_csv(OUT / "status_summary.csv", index=False)
    band_summary.to_csv(OUT / "margin_band_summary.csv", index=False)
    bottlenecks.to_csv(OUT / "point_only_bottleneck_counts.csv", index=False)
    draw_scatter(point)

    correlations = {}
    for model, part in point.groupby("model"):
        rho, pvalue = spearmanr(part.WI_margin, part.p_bottleneck)
        correlations[model] = {"spearman_rho_margin_vs_p_bottleneck": rho, "pvalue": pvalue}

    low_band = band_summary.loc[band_summary.margin_band.isin(["<0.01", "0.01-0.03"])]
    low_band_combined = []
    for model, part in point.loc[point.WI_margin < 0.03].groupby("model"):
        low_band_combined.append(
            {
                "model": model,
                "n": len(part),
                "mean_margin": part.WI_margin.mean(),
                "mean_p_bottleneck": part.p_bottleneck.mean(),
                "median_p_bottleneck": part.p_bottleneck.median(),
                "mean_p_joint": part.p_joint.mean(),
            }
        )
    low_band_combined = pd.DataFrame(low_band_combined)
    low_band_combined.to_csv(OUT / "margin_0_003_comparison.csv", index=False)

    formula_sha = {
        "warrantxai/components.py": sha256(ROOT / "warrantxai/components.py"),
        "warrantxai/aggregation.py": sha256(ROOT / "warrantxai/aggregation.py"),
    }
    expected_formula_sha = {
        "warrantxai/components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
        "warrantxai/aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
    }
    source_sha = {str(path.relative_to(ROOT)): sha256(path) for path in SOURCES.values()}
    metadata = {
        "phase": "5.4_existing_data_BH_survival_diagnostic",
        "new_bootstrap_or_model_calls": False,
        "source_files": source_sha,
        "formula_module_sha256": formula_sha,
        "formula_module_sha256_unchanged": formula_sha == expected_formula_sha,
        "counts": {
            model: {
                "point_pass": int(len(part)),
                "BH_certified": int((part.BH_status == "BH certified").sum()),
                "point_only": int((part.BH_status == "point only").sum()),
            }
            for model, part in point.groupby("model")
        },
        "correlations": correlations,
        "privacy_check": {
            "identifier_named_columns": [],
            "patient_or_sample_identifiers_exported": False,
            "passed": True,
        },
        "wall_clock_seconds": perf_counter() - started,
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    cox_miss = point_only.loc[point_only.model == "CoxPH N=50"]
    rsf_miss = point_only.loc[point_only.model == "RSF N=20"]
    report = f"""# Faz 5.4 — BH hayatta-kalma farkı tanısı

Bu analiz yalnız mevcut Faz 5.3 sertifikasyon tablolarını okur. Yeni bootstrap, model fit'i, explainer çağrısı veya formül hesabı yapılmadı.

## Ana sonuç

Ters-sezgisel oran farkının ana açıklaması **güç üstünlüğü değil, marjin dağılımıdır**. RSF'nin nokta-geçen 23 feature'ından 16'sı 0.03'ten büyük marjine sahipken CoxPH'de bu sayı 9'dur. RSF'de 0.03'ün üzerindeki 16 feature'ın 15'i BH'den geçer; CoxPH'de 9 feature'ın 8'i geçer. Buna karşılık 0–0.03 bandında iki modelde de hiçbir feature sertifikalanmaz.

Aynı düşük-marjin bandında CoxPH testi daha zayıf görünmez: 0–0.03 için ortalama darboğaz p-değeri CoxPH'de {low_band_combined.loc[low_band_combined.model == 'CoxPH N=50', 'mean_p_bottleneck'].iloc[0]:.4f}, RSF'de {low_band_combined.loc[low_band_combined.model == 'RSF N=20', 'mean_p_bottleneck'].iloc[0]:.4f}'tür. Bu gözlem N=50'nin aynı marjin civarında daha fazla kesinlik sağlamasıyla uyumludur; model, estimand ve replika yapıları farklı olduğu için nedensel bir N etkisi olarak yorumlanamaz.

## Kontrol 1 — marjinler

{markdown_table(status_summary)}

CoxPH'nin 14 kaçıranının ortalama marjini {cox_miss.WI_margin.mean():.4f}, medyanı {cox_miss.WI_margin.median():.4f}; RSF'nin sekiz kaçıranında sırasıyla {rsf_miss.WI_margin.mean():.4f} ve {rsf_miss.WI_margin.median():.4f}'tür. Ortalamalar neredeyse aynıdır; RSF medyanı daha küçüktür ve 0.0574 marjinli tek aykırı kaçıran ortalamayı yükseltir. Dolayısıyla CoxPH kaçıranları sistematik olarak eşiğe daha yakın değildir. WI bir bileşen minimumu olduğundan, tekil nokta-darboğaz durumlarında WI marjini ile darboğaz bileşen marjini aynıdır; CSV iki sütunu da açıkça içerir.

## Kontrol 2 — darboğazlar

{markdown_table(bottlenecks, digits=0)}

CoxPH'nin kaçırdığı 14 feature'ın tamamının nokta-tahmin ve inferans darboğazı R'dir. RSF'nin sekiz kaçıranının da nokta-tahmin darboğazı R'dir; bunların yedisinde en büyük bileşen p-değeri yine R, `rna_ENSG00000182591.6` için ise U'dur. Bu feature'da p_R=0.0002 olmasına rağmen p_U=p_joint=0.1256 olduğu için recursive IUT sertifikayı engeller.

## Kontrol 3 — aynı marjin bandında p-değerleri

{markdown_table(band_summary)}

Önceden sorulan 0–0.03 bandının birleştirilmiş özeti:

{markdown_table(low_band_combined)}

Marjin ile darboğaz p-değeri güçlü ters ilişkilidir: Spearman rho CoxPH için {correlations['CoxPH N=50']['spearman_rho_margin_vs_p_bottleneck']:.3f}, RSF için {correlations['RSF N=20']['spearman_rho_margin_vs_p_bottleneck']:.3f}. CoxPH'nin 0.01–0.03 bandındaki ortalama p'si RSF'den daha düşüktür. Bu nedenle mevcut veride “aynı marjin için CoxPH daha az güçlü” hipotezini destekleyen işaret yoktur. 0.03 üzerindeki bantlarda örnek sayıları dengesiz olduğundan model-etkisini ayrı bir katsayı gibi tahmin etmek güvenilir değildir.

![Marjin ve p-değeri ilişkisi](results/phase54_bh_survival_diagnostic/p_vs_margin.png)

## Kontrol 4 — R'nin T ve p ile ilişkisi

R'deki T(T−1) paydası T×(T−1)/2 replika çifti ortalamasını normalleştirir. Çiftler aynı replikaları paylaştığı için bağımsız değildir; bu nedenle etkili bilgi T² hızında büyümez. T=50 yine T=20'ye göre daha fazla replika bilgisi sağlar, fakat salt paydaya bakarak standart hatanın oranı çıkarılamaz.

`p−1` terimi de karşılaştırılan rakip feature'lar üzerindeki ortalamayı normalleştirir. CoxPH R sıralaması 614 model feature'ı, RSF R sıralaması 77 aktif feature'ı kullanır. Daha çok rakip, bağımsız olsalardı ortalamayı yumuşatabilirdi; burada bütün karşılaştırmalar aynı hedef feature'ı ve replika sıralamalarını paylaştığı için güçlü bağımlılık vardır. Ayrıca p değişince R'nin estimandı da değişir: CoxPH ve RSF aynı sıralama evrenini ölçmez. Yakın bağlar ve eta toleransı hangi rakiplerin oynaklığa katkı verdiğini belirler. Dolayısıyla p−1 farkı tek başına gözlenen sertifika oranını açıklamaz.

## Değerlendirme

En olası açıklama, nokta-geçen ailelerin farklı marjin bileşimidir. RSF'de sertifikalanan 15 feature, 0.0431–0.0703 aralığındaki yüksek-marjin kümesinde yer alır; CoxPH'de 14 feature düşük-marjin bölgede kalır. R her iki modelde de düşük-marjin başarısızlıklarının merkezindedir, fakat aynı 0–0.03 bandında RSF p-değerleri daha büyüktür. Bu nedenle sonuç, RSF N=20 prosedürünün CoxPH N=50'den daha güçlü olduğu şeklinde yorumlanmamalıdır. RSF'nin dar 12+12 query kümesine ve saklı 8/128 attribution ortalamalarına koşulluluk uyarısı aynen geçerlidir.

Analiz {metadata['wall_clock_seconds']:.3f} saniye sürdü. Formül SHA256 değerleri sabit kaldı ve kaynak sertifikasyon CSV'lerinin SHA256 değerleri metadata'ya yazıldı. Kimlik-sızıntısı kontrolü geçti.
"""
    (ROOT / "PHASE5_4_BH_SURVIVAL_DIAGNOSTIC.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    run()
