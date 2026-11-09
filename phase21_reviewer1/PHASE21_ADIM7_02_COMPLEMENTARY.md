# Faz 21 ADIM 7.2 — tamamlayıcı sayılar

RSF-v3 BH sertifikalı 14 satırda unexercised (Q_j=0): **1**.

Bağlayıcı dağılım:
```
      model bottleneck  n
   CoxPH-v3    G_claim  2
   CoxPH-v3          R 63
DeepSurv-v3          R 65
     RSF-v3    G_claim  1
     RSF-v3          R 64
```

Cox 7/4/4 tasarımı: 12 TCGA + 12 METABRIC sorgu, 65-feature R evreni ve 65 aile; `matched_three_model.csv` ile aynıdır, tam 153/1,466 veya 610-girdi sonucu değildir.

Standardizasyon:
```
      model                                           method  point  median_G_claim
   CoxPH-v3                     official_cohort_standardized      7        0.973033
   CoxPH-v3 tcga_train_scaler_only_no_cohort_standardization      4        0.742080
DeepSurv-v3                     official_cohort_standardized      0        0.990885
DeepSurv-v3 tcga_train_scaler_only_no_cohort_standardization      0        0.919884
     RSF-v3                     official_cohort_standardized     30        0.959457
     RSF-v3 tcga_train_scaler_only_no_cohort_standardization     16        0.836997
```

DeepSurv validation-test epoch Spearman: 0.3878592791780959.

V1-v3 örtüşme:
```
 v1_features  v1_features_in_v3_family  v1_certified_model_feature_pairs  pairs_feature_in_v3_family  pairs_certified_in_v3
          27                         5                                20                           6                      6
```
