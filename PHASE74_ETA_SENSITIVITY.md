# Faz 7.4 — η rank-bottleneck duyarlılığı

Ana analiz **η=0.02** olarak kalır. R, saklı T=50 candidate attribution tensörlerinden η∈{0.01,0.02,0.05,0.10} için yeniden değerlendirildi; diğer dört bileşen değiştirilmedi.

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

η=0.10'da R-darboğaz sayıları CoxPH/RSF/DeepSurv için **4/25/25**. Tam feature düzeyi R/WI/bottleneck değerleri `PHASE74_ETA_FEATURE_VALUES.csv` içindedir.

R darboğaz baskınlığı η=0.10'da RSF ve DeepSurv için sürerken (25/27 ve 25/27), CoxPH'de 4/27'ye düşer. Bu nedenle η-seçimine dayanıklılık iki nonlinear modelde güçlü, CoxPH'de ise geniş tolerans değerinde zayıftır.
