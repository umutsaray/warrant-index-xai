# Faz 7.2 — Sonuç defteri

## Eski kirli split ile temiz hasta-gruplu fold-0 karşılaştırması

| model | dirty_split_c_index | clean_fold0_c_index | old_family_m | old_point | old_BH | old_BY | new_family_m | new_point | new_BH | new_BY | new_bottleneck_distribution |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH | 0.7065 | 0.6677 | 77 | 22 | 8 | 8 | 27 | 9 | 8 | 8 | {"G_claim": 3, "R": 24} |
| RSF | 0.7148 | 0.5385 | 77 | 20 | 10 | 10 | 27 | 17 | 12 | 12 | {"G_claim": 3, "R": 24} |
| DeepSurv | 0.7457 | 0.4965 | 77 | 8 | 0 | 0 | 27 | 0 | 0 | 0 | {"R": 27} |

Bu tablo yalnız sayısal karşılaştırmadır. Eski aile m=77, yeni train-only seçim ve METABRIC kesişimiyle m=27; dolayısıyla sertifikalı feature sayıları farklı aileler üzerinde hesaplandı. C-index satırları eski kirli test skoru ile XAI audit'te kullanılan temiz fold-0 referans model skorunu ayrı gösterir.

Sayısal olarak C-index değişimleri CoxPH/RSF/DeepSurv için sırasıyla **-0.038796**, **-0.176279**, **-0.249202**; yeni nokta/BH/BY sayıları tabloda verilmiştir. Burada performans değişimi ile açıklama-kalitesi değişimi arasında nedensel yorum yapılmamıştır.

## Kalibrasyon kapıları

| model | controls | holm_flagged_controls | gate_passed | inference_scope |
| --- | --- | --- | --- | --- |
| CoxPH_N50 | 30 | 0 | True | 153 TCGA + 1466 METABRIC fixed queries; exact closed-form Cox attribution |
| RSF_N50 | 30 | 0 | True | fixed 12+12 queries and stored 8/128 attribution averages; permutation-path uncertainty excluded |
| DeepSurv_N50 | 30 | 0 | True | fixed 12+12 queries and deterministic nested DeepExplainer k16/k128 backgrounds |

## Süreler

| stage | seconds |
| --- | --- |
| panel_seconds | 3.2848 |
| cox_seconds | 2154.2317 |
| rsf_accumulated_compute_seconds | 883.9458 |
| deepsurv_accumulated_compute_seconds | 302.0825 |
| calibration_seconds | 74.5991 |
| certification_seconds | 8.0489 |
| sum_of_measured_stage_seconds | 3426.1929 |

Ana Faz 7.2 orchestration'ının gözlenen uçtan uca wall-clock süresi, ilk panel betiğinin oluşturulmasından bu ledger'ın ilk yazımına kadar **2522.4 saniye (42.04 dakika)** idi. CoxPH, RSF ve DeepSurv koşuları kısmen paralel çalıştığı için bu değer yukarıdaki adım süreleri toplamından küçüktür.

RSF sonuçları sabit 12+12 query ve saklı 8/128 attribution ortalamalarına; DeepSurv sonuçları sabit 12+12 query ve nested k16/k128 background tasarımına koşulludur. Formül SHA256 değişmedi, kimlik-sızıntısı kontrolleri geçti. Ledger üretimi **0.012 saniye** sürdü.
