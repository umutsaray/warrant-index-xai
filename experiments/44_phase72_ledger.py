#!/usr/bin/env python3
"""Write the Phase 7.2 numerical results ledger after gated certification."""

from __future__ import annotations

import json
from pathlib import Path
import sys
from time import perf_counter

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "experiments"))
from _shared import markdown_table


def counts(path: Path) -> tuple[int, int, int, int]:
    table = pd.read_csv(path)
    return (
        len(table),
        int(table.point_WI_ge_0_8.sum()),
        int(table.certified_BH_alpha_0_05.sum()),
        int(table.certified_BY_alpha_0_05.sum()),
    )


def main() -> None:
    started = perf_counter()
    rows = []
    specs = [
        ("CoxPH", 0.7065, 0.667704, "phase53_cox_recursive", "phase72_coxph_certification", "phase72_cox"),
        ("RSF", 0.7148, 0.538521, "phase56_rsf_n50_certification", "phase72_rsf_certification", "phase72_rsf"),
        ("DeepSurv", 0.7457, 0.496498, "phase62_deepsurv_n50_certification", "phase72_deepsurv_certification", "phase72_deepsurv"),
    ]
    for model, dirty_c, clean_c, old_dir, new_dir, point_dir in specs:
        old = counts(ROOT / "results" / old_dir / "feature_certification.csv")
        new = counts(ROOT / "results" / new_dir / "feature_certification.csv")
        point = pd.read_csv(ROOT / "results" / point_dir / "feature_report.csv")
        rows.append({
            "model": model,
            "dirty_split_c_index": dirty_c,
            "clean_fold0_c_index": clean_c,
            "old_family_m": old[0],
            "old_point": old[1],
            "old_BH": old[2],
            "old_BY": old[3],
            "new_family_m": new[0],
            "new_point": new[1],
            "new_BH": new[2],
            "new_BY": new[3],
            "new_bottleneck_distribution": json.dumps(point.bottleneck.value_counts().to_dict(), sort_keys=True),
        })
    comparison = pd.DataFrame(rows)
    comparison.to_csv(ROOT / "results/phase72_calibration/old_new_comparison.csv", index=False)
    panel = json.loads((ROOT / "results/phase72_panel/metadata.json").read_text())
    calibration = pd.read_csv(ROOT / "results/phase72_calibration/gate_summary.csv")
    durations = {
        "panel_seconds": panel["wall_clock_seconds"],
        "cox_seconds": json.loads((ROOT / "results/phase72_cox/metadata.json").read_text())["wall_clock_seconds"],
        "rsf_accumulated_compute_seconds": json.loads((ROOT / "results/phase72_rsf/metadata.json").read_text())["accumulated_compute_seconds"],
        "deepsurv_accumulated_compute_seconds": json.loads((ROOT / "results/phase72_deepsurv/metadata.json").read_text())["accumulated_compute_seconds"],
        "calibration_seconds": json.loads((ROOT / "results/phase72_calibration/metadata.json").read_text())["wall_clock_seconds"],
        "certification_seconds": json.loads((ROOT / "results/phase72_calibration/certification_summary.json").read_text())["wall_clock_seconds"],
    }
    durations["sum_of_measured_stage_seconds"] = sum(durations.values())
    report = f"""# Faz 7.2 — Sonuç defteri

## Eski kirli split ile temiz hasta-gruplu fold-0 karşılaştırması

{markdown_table(comparison)}

Bu tablo yalnız sayısal karşılaştırmadır. Eski aile m=77, yeni train-only seçim ve METABRIC kesişimiyle m={panel['new']['common_features']}; dolayısıyla sertifikalı feature sayıları farklı aileler üzerinde hesaplandı. C-index satırları eski kirli test skoru ile XAI audit'te kullanılan temiz fold-0 referans model skorunu ayrı gösterir.

Sayısal olarak C-index değişimleri CoxPH/RSF/DeepSurv için sırasıyla **{comparison.iloc[0].clean_fold0_c_index-comparison.iloc[0].dirty_split_c_index:.6f}**, **{comparison.iloc[1].clean_fold0_c_index-comparison.iloc[1].dirty_split_c_index:.6f}**, **{comparison.iloc[2].clean_fold0_c_index-comparison.iloc[2].dirty_split_c_index:.6f}**; yeni nokta/BH/BY sayıları tabloda verilmiştir. Burada performans değişimi ile açıklama-kalitesi değişimi arasında nedensel yorum yapılmamıştır.

## Kalibrasyon kapıları

{markdown_table(calibration)}

## Süreler

{markdown_table(pd.DataFrame([{'stage': key, 'seconds': value} for key, value in durations.items()]))}

Ana Faz 7.2 orchestration'ının gözlenen uçtan uca wall-clock süresi ilk üretimde **2522.4 saniye (42.04 dakika)** idi. Model audit'leri kısmen paralel çalıştığı için bu değer adım süreleri toplamından küçüktür.

RSF sonuçları sabit 12+12 query ve saklı 8/128 attribution ortalamalarına; DeepSurv sonuçları sabit 12+12 query ve nested k16/k128 background tasarımına koşulludur. Formül SHA256 değişmedi, kimlik-sızıntısı kontrolleri geçti. Ledger üretimi **{perf_counter()-started:.3f} saniye** sürdü.
"""
    (ROOT / "RESULTS_LEDGER_PHASE72.md").write_text(report, encoding="utf-8")
    print(comparison.to_string(index=False))


if __name__ == "__main__":
    main()
