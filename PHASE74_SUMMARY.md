# Faz 7.4 — Eşik ve rank-tolerance duyarlılık özeti

## ADIM 0 — Fizibilite

**Yeniden patient-cluster bootstrap gerekmedi.** B=4999 draw matrisi ayrıca saklanmamış olsa da aynı inferential resample akışı, aşağıdaki dondurulmuş T=50 attribution tensörlerinden ve seed'den birebir yeniden üretildi. Model veya explainer çağrısı yapılmadı.

| Model | source_phi | external_phi | coefficients | candidate | reference |
| --- | --- | --- | --- | --- | --- |
| CoxPH | [50, 153, 610] | [50, 1466, 27] | [50, 610] | - | - |
| RSF | - | - | - | [2, 50, 12, 27] | [2, 50, 12, 27] |
| DeepSurv | - | - | - | [2, 50, 12, 27] | [2, 50, 12, 27] |

## Supplementary tablo S1 — τ taraması

Referans satır `reference_primary=True` olan **τ=0.80** satırıdır.

| tau | reference_primary | CoxPH_point | CoxPH_BH | CoxPH_BY | RSF_point | RSF_BH | RSF_BY | DeepSurv_point | DeepSurv_BH | DeepSurv_BY |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.7000 | False | 22 | 17 | 15 | 26 | 25 | 22 | 21 | 6 | 6 |
| 0.7500 | False | 16 | 9 | 9 | 22 | 16 | 12 | 6 | 1 | 1 |
| 0.8000 | True | 9 | 8 | 8 | 17 | 12 | 12 | 0 | 0 | 0 |
| 0.8500 | False | 2 | 1 | 1 | 13 | 12 | 12 | 0 | 0 | 0 |
| 0.9000 | False | 1 | 1 | 1 | 13 | 10 | 10 | 0 | 0 | 0 |

## Supplementary tablo S2 — η taraması

Referans satırlar `reference_primary=True` olan **η=0.02** satırlarıdır.

| eta | reference_primary | model | R_mean | R_std | R_min | R_q25 | R_median | R_q75 | R_max | R_bottleneck_count | point_adequate_tau_0_8 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0100 | False | CoxPH | 0.7343 | 0.1037 | 0.6170 | 0.6483 | 0.6956 | 0.8316 | 0.9833 | 24 | 7 |
| 0.0100 | False | RSF | 0.8362 | 0.0707 | 0.7196 | 0.7724 | 0.8381 | 0.9139 | 0.9159 | 25 | 15 |
| 0.0100 | False | DeepSurv | 0.6661 | 0.0433 | 0.6128 | 0.6401 | 0.6528 | 0.6837 | 0.7858 | 27 | 0 |
| 0.0200 | True | CoxPH | 0.7741 | 0.0756 | 0.6724 | 0.7045 | 0.7674 | 0.8291 | 0.9778 | 24 | 9 |
| 0.0200 | True | RSF | 0.8550 | 0.0784 | 0.7266 | 0.7773 | 0.9038 | 0.9272 | 0.9278 | 24 | 17 |
| 0.0200 | True | DeepSurv | 0.7264 | 0.0359 | 0.6555 | 0.7029 | 0.7301 | 0.7473 | 0.7904 | 27 | 0 |
| 0.0500 | False | CoxPH | 0.8942 | 0.0609 | 0.7602 | 0.8688 | 0.9234 | 0.9322 | 0.9612 | 23 | 23 |
| 0.0500 | False | RSF | 0.8908 | 0.0770 | 0.7301 | 0.8486 | 0.9397 | 0.9441 | 0.9457 | 25 | 21 |
| 0.0500 | False | DeepSurv | 0.8413 | 0.0520 | 0.7398 | 0.8303 | 0.8542 | 0.8802 | 0.9003 | 27 | 21 |
| 0.1000 | False | CoxPH | 0.9679 | 0.0566 | 0.7923 | 0.9828 | 0.9916 | 0.9943 | 0.9964 | 4 | 25 |
| 0.1000 | False | RSF | 0.9158 | 0.0731 | 0.7530 | 0.9313 | 0.9540 | 0.9557 | 0.9621 | 25 | 22 |
| 0.1000 | False | DeepSurv | 0.9191 | 0.0380 | 0.8280 | 0.8855 | 0.9391 | 0.9474 | 0.9572 | 25 | 27 |

η=0.10'da R-darboğaz sayıları CoxPH/RSF/DeepSurv için **4/25/25** oldu. R baskınlığı RSF ve DeepSurv'de geniş η değerinde sürerken CoxPH'de belirgin biçimde azaldı. DeepSurv τ=0.70'te 21/6/6, τ=0.75'te 6/1/1 ve referans τ=0.80'de 0/0/0 olduğundan DeepSurv sertifikasyon sonucu aşağı yönlü τ değişimine duyarlıdır.

## Süre

| Step | Seconds |
| --- | --- |
| ADIM 0 — feasibility | 0.1199 |
| ADIM 1 — tau scan | 40.5169 |
| ADIM 2 — eta scan | 6.7746 |

Toplam wall-clock **47.417 saniye**. Ana Faz 7.2 sonuçları ve artefaktları değiştirilmedi; formül SHA256 değerleri doğrulandı.
