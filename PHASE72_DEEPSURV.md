# Faz 7.2 — DeepSurv N=50 temiz audit

Hasta-gruplu fold-0 referans DeepSurv mimarisiyle 50 patient-cluster bootstrap model fit edildi. Temiz 27-feature panelde **0/27** feature nokta tahminde WI≥0.8 geçti.

| block | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 0.9993 | 0.9897 | 0.9934 | 0.9831 | 0.7634 | 0.7634 |
| rna | 22 | 0.9959 | 0.9627 | 0.9851 | 0.9389 | 0.7180 | 0.7180 |

Bottleneck dağılımı: `{'R': 27}`. En yüksek göreli additivity hatası candidate/reference için **1.34e-06/8.95e-07**; önceden dondurulan %5 kapısı geçildi.

Sonuç sabit deterministik 12+12 query ve nested DeepExplainer k=16/128 background bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen süre **302.1 sn**, invocation wall-clock **1.4 sn**.


## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m=27 ailenin tamamında B=4999 recursive_studentized uygulandı. Nokta/BH/BY sayıları **0/0/0**.

BH ile sertifikalanan feature yok.
