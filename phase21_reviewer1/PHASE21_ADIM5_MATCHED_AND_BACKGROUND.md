# Faz 21 ADIM 5 — eşleştirilmiş Cox ve DeepSurv background tanıları

## CoxPH-v2: aynı 12+12 sorgu, aynı 27-feature sıralama evreni

Saklı 50 katsayı replikası kullanıldı; yeni model fit edilmedi. F dahil bütün bileşenler iki kohortlu 12+12 tensör üzerinde, R yalnız 27-feature evreninde hesaplandı. Recursive-studentized B=4,999 ve 27-feature BH/BY ailesi kullanıldı.

```text
 family_m  point  BH  BY             bottlenecks
       27      3   3   3 {'R': 25, 'G_claim': 2}
```

## DeepSurv: background k duyarlılığı

Tek dondurulmuş referans modelde aynı sorgular ve deterministik nested k=16/64/128/256 TCGA background'ları kullanıldı. Korelasyonlar ve MARD, özellik başına iki kohortun 24 attribution değerinin k=128 vektörüyle karşılaştırmasıdır; bu bölüm sertifikasyon değildir.

```text
 background_k  median_pearson  min_pearson  median_spearman  median_MARD  max_MARD  median_G_claim
           16        0.994700     0.839212         0.946400     0.489467  3.291562        0.983328
           64        0.998931     0.981684         0.986957     0.193681  0.893273        0.984395
          128        1.000000     1.000000         1.000000     0.000000  0.000000        0.985542
          256        0.999467     0.949846         0.993913     0.228719  1.108635        0.985869
```

Wall-clock: 19.832 sn. Dondurulmuş formüller ve resmi sonuçlar değiştirilmedi.
