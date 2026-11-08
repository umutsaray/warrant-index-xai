# Faz 21 ADIM 4 — METABRIC standardizasyon duyarlılığı

RSF ve DeepSurv için bootstrap model nesneleri saklı olmadığından, izin verilen **tek referans model / replika-yok** yedeği kullanıldı. Yeni model fit edilmedi. F/S/U/R resmi Faz 7.2 satırlarından sabit tutuldu; yalnız referans modelin aynı 12+12 sorgudaki G_claim değeri değiştirilerek tanısal WI verildi. CoxPH karşılaştırması saklı N=50 Faz 15 yeniden üretimidir. Dolayısıyla modeller arası inferans kapsamları aynı değildir ve bu tablo resmi sertifikasyon değildir.

```text
   model          standardization  n  point  median_G_claim  median_WI
CoxPH-v2 cohort_specific_official 27      9        0.958251   0.758629
CoxPH-v2  tcga_train_mean_sd_only 27      3        0.000000   0.000000
DeepSurv cohort_specific_official 27      0        0.983328   0.730094
DeepSurv  tcga_train_mean_sd_only 27      0        0.000000   0.000000
     RSF cohort_specific_official 27     16        0.982562   0.822449
     RSF  tcga_train_mean_sd_only 27     13        0.784267   0.761303
```

Darboğaz sayımları:

```text
   model          standardization bottleneck  n
CoxPH-v2 cohort_specific_official    G_claim  3
CoxPH-v2 cohort_specific_official          R 24
CoxPH-v2  tcga_train_mean_sd_only    G_claim 23
CoxPH-v2  tcga_train_mean_sd_only          R  4
DeepSurv cohort_specific_official          R 27
DeepSurv  tcga_train_mean_sd_only    G_claim 20
DeepSurv  tcga_train_mean_sd_only          R  7
     RSF cohort_specific_official    G_claim  4
     RSF cohort_specific_official          R 23
     RSF  tcga_train_mean_sd_only    G_claim 14
     RSF  tcga_train_mean_sd_only          R 13
```

Wall-clock: 47.956 sn. Dondurulmuş sonuç ve formül dosyaları değiştirilmedi.
