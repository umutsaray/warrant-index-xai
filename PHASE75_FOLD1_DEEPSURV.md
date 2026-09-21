# Faz 7.5 fold 1 — DeepSurv N=50 temiz audit

Hasta-gruplu fold-1 referans DeepSurv mimarisiyle 50 patient-cluster bootstrap model fit edildi. Temiz 80-feature panelde **7/80** feature nokta tahminde WI≥0.8 geçti.

| block | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 0.9993 | 0.9906 | 0.9944 | 0.9883 | 0.7558 | 0.7558 |
| rna | 75 | 0.9980 | 0.9737 | 0.9863 | 0.9604 | 0.7464 | 0.7464 |

Bottleneck dağılımı: `{'R': 80}`. En yüksek göreli additivity hatası candidate/reference için **4.45e-07/0.00331**; önceden dondurulan %5 kapısı geçildi.

Sonuç sabit deterministik 12+12 query ve nested DeepExplainer k=16/128 background bütçesine koşulludur. Formül SHA256 değişmedi; hasta kimliği dışa aktarılmadı. Birikimli ölçülen süre **400.8 sn**, invocation wall-clock **407.3 sn**.


## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m=80 ailenin tamamında B=4999 recursive_studentized uygulandı. Nokta/BH/BY sayıları **7/0/0**.

BH ile sertifikalanan feature yok.
