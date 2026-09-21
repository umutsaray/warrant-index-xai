#!/usr/bin/env python3
"""Phase 7.3: format and derive feature-level results from frozen Phase 7.2 artifacts."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import markdown_table

EXPECTED_HASHES = {
    "components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa",
    "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68",
}
EXPECTED_COUNTS = {
    "CoxPH-v2": (27, 9, 8, 8),
    "RSF": (27, 17, 12, 12),
    "DeepSurv": (27, 0, 0, 0),
}
SOURCES = {
    "CoxPH-v2": ROOT / "results/phase72_coxph_certification/feature_certification.csv",
    "RSF": ROOT / "results/phase72_rsf_certification/feature_certification.csv",
    "DeepSurv": ROOT / "results/phase72_deepsurv_certification/feature_certification.csv",
}
OUTPUT_STEMS = {
    "CoxPH-v2": "PHASE73_COXPH_FULLTABLE",
    "RSF": "PHASE73_RSF_FULLTABLE",
    "DeepSurv": "PHASE73_DEEPSURV_FULLTABLE",
}


def sha(path: Path) -> str:
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


def feature_type(feature: str) -> str:
    if feature.startswith("rna_"):
        return "RNA"
    if feature.startswith("clin_numeric__age_at_index"):
        return "clinical (age)"
    if feature.startswith("clin_categorical__ajcc_pathologic_stage"):
        return "clinical (AJCC stage indicator)"
    raise RuntimeError(f"Unclassified Phase 7.2 feature: {feature}")


def yn(value: bool) -> str:
    return "evet" if bool(value) else "hayır"


def label_map() -> dict[str, str]:
    mapping = pd.read_csv(ROOT / "results/phase72_panel/feature_mapping.csv")
    labels = dict(zip(mapping.feature, mapping.symbol, strict=True))
    for feature in json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())["active_features"]:
        if feature.startswith("clin_numeric__age_at_index"):
            labels[feature] = "Age at index"
        elif feature.startswith("clin_categorical__ajcc_pathologic_stage"):
            labels[feature] = feature.rsplit("_", 1)[-1].replace("Stage ", "AJCC Stage ")
    return labels


def validate_sources() -> dict[str, pd.DataFrame]:
    hashes = {name: sha(ROOT / "warrantxai" / name) for name in EXPECTED_HASHES}
    if hashes != EXPECTED_HASHES:
        raise RuntimeError(f"Frozen formula hashes changed: {hashes}")
    tables = {model: pd.read_csv(path) for model, path in SOURCES.items()}
    panel = json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())["active_features"]
    mismatches = []
    for model, table in tables.items():
        observed = (
            len(table),
            int(table.point_WI_ge_0_8.sum()),
            int(table.certified_BH_alpha_0_05.sum()),
            int(table.certified_BY_alpha_0_05.sum()),
        )
        if observed != EXPECTED_COUNTS[model]:
            mismatches.append({"model": model, "expected": EXPECTED_COUNTS[model], "observed": observed})
        if table.feature.tolist() != panel:
            mismatches.append({"model": model, "expected": "ordered Phase 7.2 panel", "observed": "different feature order"})
        component_min = table[["F", "S", "U", "G_claim", "R"]].min(axis=1)
        if not np.allclose(table.WI, component_min, atol=0, rtol=0):
            mismatches.append({"model": model, "expected": "WI=min(components)", "observed": "numerical mismatch"})
        if not np.allclose(table.p_joint, table[["p_F", "p_S", "p_U", "p_G_claim", "p_R"]].max(axis=1)):
            mismatches.append({"model": model, "expected": "p_joint=max(component p)", "observed": "numerical mismatch"})
    if mismatches:
        raise RuntimeError("Phase 7.3 ADIM 1 assertion failed; no outputs written: " + json.dumps(mismatches))
    return tables


def make_full_tables(tables: dict[str, pd.DataFrame], labels: dict[str, str]) -> None:
    for model, source in tables.items():
        out = pd.DataFrame({
            "Feature": source.feature.map(labels),
            "Feature ID": source.feature,
            "Feature type": source.feature.map(feature_type),
            "F": source.F,
            "S": source.S,
            "U": source.U,
            "G^claim": source.G_claim,
            "R": source.R,
            "WI": source.WI,
            "Point adequate": source.point_WI_ge_0_8.map(yn),
            "p_joint": source.p_joint,
            "BH-q": source.q_BH,
            "BY-q": source.q_BY,
            "BH certified": source.certified_BH_alpha_0_05.map(yn),
            "BY certified": source.certified_BY_alpha_0_05.map(yn),
        })
        stem = OUTPUT_STEMS[model]
        out.to_csv(ROOT / f"{stem}.csv", index=False)
        counts = EXPECTED_COUNTS[model]
        note = f"""# Faz 7.3 — {model} tam 27-özellik tablosu

Kaynak: `{SOURCES[model].relative_to(ROOT).as_posix()}`. RNA sembolleri `results/phase72_panel/feature_mapping.csv` dosyasından alınmıştır; `Feature ID` yeniden üretilebilir kanonik panel adıdır.

Satır/nokta-yeterli/BH/BY doğrulaması: **{counts[0]}/{counts[1]}/{counts[2]}/{counts[3]}**. Point-adequate sütunu, Faz 7.2'de saklanan dondurulmuş `WI≥0.8` kararını gösterir.

{markdown_table(out)}
"""
        (ROOT / f"{stem}.md").write_text(note, encoding="utf-8")


def certified_set(table: pd.DataFrame, method: str) -> set[str]:
    column = f"certified_{method}_alpha_0_05"
    return set(table.loc[table[column], "feature"])


def intersection_report(tables: dict[str, pd.DataFrame], labels: dict[str, str]) -> dict:
    cox, rsf = tables["CoxPH-v2"], tables["RSF"]
    sections, payload = [], {}
    for method in ("BH", "BY"):
        cset, rset = certified_set(cox, method), certified_set(rsf, method)
        both, conly, ronly, union = cset & rset, cset - rset, rset - cset, cset | rset
        payload[method] = {
            "cox_only": sorted(conly), "rsf_only": sorted(ronly), "both": sorted(both),
            "union": sorted(union), "counts": {"cox_only": len(conly), "rsf_only": len(ronly),
                                                     "both": len(both), "union": len(union)},
        }
        rows = []
        ci, ri = cox.set_index("feature"), rsf.set_index("feature")
        for feature in sorted(both):
            rows.append({
                "Feature": labels[feature], "Feature ID": feature,
                "CoxPH WI": ci.loc[feature, "WI"], "CoxPH p_joint": ci.loc[feature, "p_joint"],
                "CoxPH q": ci.loc[feature, f"q_{method}"],
                "RSF WI": ri.loc[feature, "WI"], "RSF p_joint": ri.loc[feature, "p_joint"],
                "RSF q": ri.loc[feature, f"q_{method}"],
            })
        table = pd.DataFrame(rows)
        sections.append(
            f"## {method} kesişimi\n\n"
            f"Venn sayıları: **CoxPH-only={len(conly)}, both={len(both)}, RSF-only={len(ronly)}**; "
            f"birleşim={len(union)}.\n\n"
            f"CoxPH-only: {', '.join(labels[x] for x in sorted(conly)) or 'yok'}  \n"
            f"RSF-only: {', '.join(labels[x] for x in sorted(ronly)) or 'yok'}\n\n"
            f"{markdown_table(table)}"
        )
    same = payload["BH"]["both"] == payload["BY"]["both"] and payload["BH"]["union"] == payload["BY"]["union"]
    note = "# Faz 7.3 — CoxPH-v2 ↔ RSF sertifikalı özellik kesişimi\n\n" + "\n\n".join(sections)
    note += f"\n\nBH ve BY sertifikalı kümeleri iki modelde de **{'aynıdır' if same else 'farklıdır'}**.\n"
    (ROOT / "PHASE73_INTERSECTION.md").write_text(note, encoding="utf-8")
    return payload


def type_counts(table: pd.DataFrame, mask: pd.Series | np.ndarray) -> dict[str, int]:
    types = table.loc[mask, "feature"].map(feature_type)
    return {
        "RNA": int((types == "RNA").sum()),
        "clinical (age)": int((types == "clinical (age)").sum()),
        "clinical (AJCC stage indicator)": int((types == "clinical (AJCC stage indicator)").sum()),
    }


def feature_type_report(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    family = next(iter(tables.values()))
    rows = []
    subsets = [("full family", np.ones(len(family), dtype=bool))]
    family_counts = type_counts(family, subsets[0][1])
    for model, table in tables.items():
        for subset, mask in [
            ("point adequate", table.point_WI_ge_0_8),
            ("BH certified", table.certified_BH_alpha_0_05),
        ]:
            counts = type_counts(table, mask)
            total = sum(counts.values())
            clinical = counts["clinical (age)"] + counts["clinical (AJCC stage indicator)"]
            rows.append({
                "Model": model, "Subset": subset, "Total": total, **counts,
                "Clinical total": clinical,
                "RNA %": np.nan if total == 0 else 100 * counts["RNA"] / total,
                "Clinical %": np.nan if total == 0 else 100 * clinical / total,
                "Clinical share minus family pp": np.nan if total == 0 else 100 * clinical / total - 100 * 5 / 27,
            })
    distribution = pd.DataFrame(rows)
    family_table = pd.DataFrame([{"Subset": "27-feature family", "Total": 27, **family_counts,
                                  "Clinical total": 5, "RNA %": 100 * 22/27, "Clinical %": 100 * 5/27}])
    note = f"""# Faz 7.3 — Özellik-tipi dağılımı

## Tam aile

{markdown_table(family_table)}

## Model alt kümeleri

{markdown_table(distribution)}

`Clinical share minus family pp`, alt kümedeki klinik pay ile tam ailenin %18.52 klinik payı arasındaki yüzde-puan farkıdır. Pozitif değer klinik özelliklerin sayısal olarak daha yüksek temsilini, negatif değer daha düşük temsilini gösterir; biyolojik yorum yapılmamıştır. Boş DeepSurv alt kümelerinde oran tanımsızdır.

Sayısal karşılaştırmada klinik pay CoxPH point/BH alt kümelerinde sırasıyla %33.33/%25.00 (+14.81/+6.48 yüzde puan), RSF point/BH alt kümelerinde %23.53/%33.33 (+5.01/+14.81 yüzde puan) oldu. Bu dört dolu alt küme tam aileye göre klinik özellikler yönünde kaymıştır; örneklem sayıları küçük olduğundan bu ifade yalnız betimseldir.
"""
    (ROOT / "PHASE73_FEATURE_TYPE_DISTRIBUTION.md").write_text(note, encoding="utf-8")
    return distribution


def sign_label(value: float) -> str:
    if value > 0:
        return "pozitif"
    if value < 0:
        return "negatif"
    return "sıfır"


def directionality_report(tables: dict[str, pd.DataFrame], labels: dict[str, str]) -> pd.DataFrame:
    with np.load(ROOT / "results/phase72_cox/attributions_n50.npz", allow_pickle=False) as saved:
        cox_phi = saved["source_phi"]
        cox_all = saved["features"].astype(str).tolist()
        cox_active = saved["active_features"].astype(str).tolist()
    with np.load(ROOT / "results/phase72_rsf/attributions_n50.npz", allow_pickle=False) as saved:
        rsf_phi = saved["candidate"]
        rsf_features = saved["features"].astype(str).tolist()
    cox_mean = cox_phi.mean(axis=(0, 1))
    rsf_mean = rsf_phi[0].mean(axis=(0, 1))
    with np.load(ROOT / "results/phase72_panel/query_matrices.npz", allow_pickle=False) as saved:
        source_full, source_query12 = saved["source_full"], saved["source_query12"]
    matched_rows = []
    for query in source_query12:
        matches = np.flatnonzero(np.all(np.isclose(source_full, query, rtol=0, atol=1e-12), axis=1))
        if len(matches) != 1:
            raise RuntimeError(f"Could not uniquely recover a frozen TCGA query row: {matches.tolist()}")
        matched_rows.append(int(matches[0]))
    cox_active_indices = np.array([cox_all.index(feature) for feature in cox_active])
    cox_matched12 = cox_phi[:, matched_rows][:, :, cox_active_indices].mean(axis=(0, 1))
    mean_lookup = {
        "CoxPH-v2": {feature: float(cox_mean[i]) for i, feature in enumerate(cox_all)},
        "RSF": {feature: float(rsf_mean[i]) for i, feature in enumerate(rsf_features)},
    }
    scopes = {
        "CoxPH-v2": "153 TCGA audit patients × 50 replicas",
        "RSF": "stored deterministic 12 TCGA queries × 50 replicas",
    }
    rows = []
    for model in ("CoxPH-v2", "RSF"):
        certified = sorted(certified_set(tables[model], "BH"))
        for feature in certified:
            value = mean_lookup[model][feature]
            rows.append({
                "Model": model, "Feature": labels[feature], "Feature ID": feature,
                "Mean signed attribution": value, "Direction": sign_label(value),
                "Averaging scope": scopes[model],
            })
    direction = pd.DataFrame(rows)
    direction.to_csv(ROOT / "PHASE73_DIRECTIONALITY.csv", index=False)
    both = certified_set(tables["CoxPH-v2"], "BH") & certified_set(tables["RSF"], "BH")
    cross = []
    for feature in sorted(both):
        cv, cm12, rv = (
            mean_lookup["CoxPH-v2"][feature],
            float(cox_matched12[cox_active.index(feature)]),
            mean_lookup["RSF"][feature],
        )
        cross.append({
            "Feature": labels[feature], "Feature ID": feature,
            "CoxPH full-153 mean": cv, "CoxPH full direction": sign_label(cv),
            "CoxPH matched-12 mean": cm12, "CoxPH matched direction": sign_label(cm12),
            "RSF mean": rv, "RSF direction": sign_label(rv),
            "Same sign on matched queries": sign_label(cm12) == sign_label(rv),
        })
    cross_table = pd.DataFrame(cross)
    cross_table.to_csv(ROOT / "PHASE73_DIRECTIONALITY_INTERSECTION.csv", index=False)
    direction_display = direction.copy()
    direction_display["Mean signed attribution"] = direction_display["Mean signed attribution"].map(lambda x: f"{x:.8g}")
    cross_display = cross_table.copy()
    for column in ["CoxPH full-153 mean", "CoxPH matched-12 mean", "RSF mean"]:
        cross_display[column] = cross_display[column].map(lambda x: f"{x:.8g}")
    note = f"""# Faz 7.3 — Sertifikalı özellik yönleri

Yön, saklı candidate attribution'ların TCGA tarafındaki aritmetik ortalamasının işaretidir. CoxPH için tüm 153 fold-0 audit hastası ve 50 replika kullanıldı. RSF artefaktı yalnız dondurulmuş 12 TCGA query için attribution sakladığından RSF yönü bu 12-query alt kümesine ve 50 replikaya koşulludur; 153 hastaya genişletmek yeni explainer çağrısı gerektireceği için yapılmadı. Modeller-arası işaret kontrolünde query farkını gidermek için CoxPH'nin saklı 153-hasta tensöründen RSF ile aynı deterministik 12 satır ayrıca seçildi. Bu seçim yeni attribution hesabı değildir. DeepSurv BH kümesi boş olduğundan yön satırı yoktur. G_signed tek başına yönü saklamaz; bu nedenle ham attribution tensörleri okundu. Tam sıfır ortalama, pozitif veya negatif yön yerine `sıfır` olarak tutuldu.

## BH-sertifikalı özellikler

{markdown_table(direction_display)}

## CoxPH–RSF BH kesişiminde işaret tutarlılığı

{markdown_table(cross_display)}

Aynı 12 query üzerinde aynı işaret sayısı **{int(cross_table['Same sign on matched queries'].sum())}/{len(cross_table)}**. Bu yalnız sayısal yön tutarlılığıdır; biyolojik yorum yapılmamıştır.
"""
    (ROOT / "PHASE73_DIRECTIONALITY.md").write_text(note, encoding="utf-8")
    return cross_table


def main() -> None:
    total_start = perf_counter()
    stage = {}
    start = perf_counter()
    tables = validate_sources()
    labels = label_map()
    make_full_tables(tables, labels)
    stage["ADIM 1 — full tables"] = perf_counter() - start

    start = perf_counter()
    intersections = intersection_report(tables, labels)
    stage["ADIM 2 — intersection"] = perf_counter() - start

    start = perf_counter()
    distribution = feature_type_report(tables)
    stage["ADIM 3 — feature types"] = perf_counter() - start

    start = perf_counter()
    directions = directionality_report(tables, labels)
    stage["ADIM 4 — directionality"] = perf_counter() - start
    total = perf_counter() - total_start

    summary_counts = pd.DataFrame([
        {"Model": model, "Family": values[0], "Point adequate": values[1], "BH": values[2], "BY": values[3]}
        for model, values in EXPECTED_COUNTS.items()
    ])
    timing = pd.DataFrame([{"Step": key, "Seconds": value} for key, value in stage.items()])
    note = f"""# Faz 7.3 — Özet

Faz 7.2'nin saklı per-feature sonuçları yeniden formatlandı; bootstrap, null-kalibrasyon, sertifikasyon, model veya explainer çağrısı yapılmadı.

{markdown_table(summary_counts)}

- Yeni tam tabloların her biri 27 satırdır ve Faz 7.2 sayıları assertion ile birebir doğrulanmıştır.
- CoxPH/RSF BH kesişimi **{intersections['BH']['counts']['both']}**, CoxPH-only **{intersections['BH']['counts']['cox_only']}**, RSF-only **{intersections['BH']['counts']['rsf_only']}** özelliktir. BH ve BY kümeleri aynıdır.
- Aile 22 RNA, 1 yaş ve 4 AJCC stage indicator özelliğinden oluşur.
- CoxPH–RSF BH kesişiminde aynı 12 query üzerindeki attribution yönü **{int(directions['Same sign on matched queries'].sum())}/{len(directions)}** özellikte aynıdır.
- Formül SHA256 değerleri değişmemiştir; çıktı tablolarında hasta veya örnek kimliği yoktur.

## Süre

{markdown_table(timing)}

Toplam wall-clock **{total:.4f} saniye**.
"""
    (ROOT / "PHASE73_SUMMARY.md").write_text(note, encoding="utf-8")
    metadata = {
        "phase": "7.3_feature_level_reporting_only",
        "sources": {model: path.relative_to(ROOT).as_posix() for model, path in SOURCES.items()},
        "asserted_counts": {model: list(values) for model, values in EXPECTED_COUNTS.items()},
        "formula_sha256": EXPECTED_HASHES,
        "model_or_explainer_calls": False,
        "bootstrap_or_calibration_calls": False,
        "privacy": {"patient_or_sample_identifiers_exported": False},
        "stage_seconds": stage,
        "wall_clock_seconds": total,
    }
    (ROOT / "results/phase73_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(summary_counts.to_string(index=False))
    print(f"BH intersection: {intersections['BH']['counts']}")
    print(f"same direction: {int(directions['Same sign on matched queries'].sum())}/{len(directions)}")
    print(f"wall={total:.4f}s")


if __name__ == "__main__":
    main()
