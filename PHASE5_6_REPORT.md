# Faz 5.6 — RSF N=50 tam koşu ve FDR sertifikasyonu

> **Koşulluluk uyarısı:** RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

## Sonuç

Winner's-curse olgunlaşma hipotezi **kısmen desteklendi**. N=20 nokta-geçen 23 feature'ın 19'unda R düştü; bu seçilmiş grupta ortalama değişim -0.0157, medyan -0.0089 oldu. Nokta-geçen sayı 23→20, BH 15→10 ve BY 12→10 düştü. Buna karşılık tüm 77 feature'da ortalama R değişimi yalnız -0.0006; 39 düşüşe karşı 38 artış var ve eşlenmiş Wilcoxon p=0.743. Genel bir R çöküşü yok; düşüş N=20'de seçilmiş kazananlarda yoğunlaşıyor.

Marjin şekli CoxPH yönüne kaydı fakat yakınsamadı. RSF'de WI marjini 0.05–0.10 bandındaki pay 15/23 (%65.2) → 9/20 (%45.0) düştü; medyan marjin 0.0565 → 0.0454. CoxPH N=50'de aynı yüksek-marjin sayısı yalnız 2/22 (%9.1), medyan 0.0225. RSF BH hayatta-kalma oranı %65.2→%50.0 geriledi, fakat CoxPH'nin %36.4 değerinin üzerinde kaldı.

## ADIM 1 — RSF N=50

> **Koşulluluk uyarısı:** RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

Onaylanan sampled-background k=16, candidate/reference 8/128 ve aynı 12+12 query ile 77/77 feature değerlendirildi. İlk 20 deterministik replika doğrulanarak korundu; aynı akışın 21–50 replikaları yeni fit edildi. Ana RSF modeli değiştirilmedi. Yeni 30 replika 987.9 sn (0.27 saat), tüm N=50 dizisinin birikimli ölçülen maliyeti 3937.2 sn (1.09 saat) oldu. Hiçbir fit uyarısı veya atlanan replika yok.

### Block agregatları

| inference_scope | block | features | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | clin | 5.0000 | 5.0000 | 1.0000 | 0.9790 | 0.9850 | 0.9677 | 0.8262 | 0.8262 |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | rna | 72.0000 | 72.0000 | 0.9998 | 0.9529 | 0.9661 | 0.9113 | 0.7561 | 0.7340 |

## ADIM 2 — N=50 null gate

> **Koşulluluk uyarısı:** RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

`recursive_studentized`, aynı 2000×999 `single_active`+`tie`, üç alpha ve 30 kontrol içinde Holm bataryasında **geçti**; 0/30 anti-conservative ihlal bulundu. Wall-clock 22.842 sn.

## ADIM 3 — tam-aile FDR

> **Koşulluluk uyarısı:** RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

| model | point | BH | BY | inference_scope | BH_survival_fraction | BY_survival_fraction |
| --- | --- | --- | --- | --- | --- | --- |
| RSF N=20 | 23.0000 | 15.0000 | 12.0000 | RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | 0.6522 | 0.5217 |
| RSF N=50 | 20.0000 | 10.0000 | 10.0000 | RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | 0.5000 | 0.5000 |
| CoxPH N=50 | 22.0000 | 8.0000 | 8.0000 | CoxPH fixed-query Phase 5.3 comparison | 0.3636 | 0.3636 |

N=50'de BH ve BY aynı 10 feature'ı sertifikaladı. Tam p/q tablosu `results/phase56_rsf_n50_certification/feature_certification.csv` içindedir. Sertifikasyon wall-clock 10.495 sn.

## ADIM 4 — N=20→N=50 karşılaştırması

> **Koşulluluk uyarısı:** RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded.

### Marjin dağılımı

| inference_scope | family | n | mean_margin | median_margin | q25 | q75 | min | max | n_0_001 | n_001_003 | n_003_005 | n_005_010 | n_ge_005 | fraction_ge_005 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | RSF N=20 | 23.0000 | 0.0440 | 0.0565 | 0.0216 | 0.0622 | 0.0002 | 0.0703 | 5.0000 | 2.0000 | 1.0000 | 15.0000 | 15.0000 | 0.6522 |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | RSF N=50 | 20.0000 | 0.0387 | 0.0454 | 0.0217 | 0.0555 | 0.0034 | 0.0670 | 3.0000 | 6.0000 | 2.0000 | 9.0000 | 9.0000 | 0.4500 |
| CoxPH fixed-query Phase 5.3 comparison | CoxPH N=50 | 22.0000 | 0.0293 | 0.0225 | 0.0121 | 0.0359 | 0.0036 | 0.1503 | 5.0000 | 8.0000 | 7.0000 | 1.0000 | 2.0000 | 0.0909 |

### R değişimi

| inference_scope | scope | n | mean_delta_R | median_delta_R | R_decreased | R_increased | wilcoxon_p |
| --- | --- | --- | --- | --- | --- | --- | --- |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | all 77 | 77.0000 | -0.0006 | -0.0002 | 39.0000 | 38.0000 | 0.7433 |
| RSF results condition on fixed 12+12 queries and stored 8/128 attribution averages; per-permutation-path uncertainty is excluded. | N20 point-pass 23 | 23.0000 | -0.0157 | -0.0089 | 19.0000 | 4.0000 | 0.0035 |

### Retention

| model | from_N | to_N | from_point_count | to_point_count | count_ratio | same_feature_retained | set_retention | lost | new |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RSF | 20.0000 | 50.0000 | 23.0000 | 20.0000 | 0.8696 | 16.0000 | 0.6957 | 7.0000 | 4.0000 |
| CoxPH | 6.0000 | 50.0000 | 27.0000 | 22.0000 | 0.8148 | 15.0000 | 0.5556 | 12.0000 | 7.0000 |

RSF'de count oranı 20/23=%87.0, birebir feature retention 16/23=%69.6'dır. CoxPH'de 22/27=%81.5 yalnız sayım oranıdır; gerçek N=6 feature kümesinin N=50'de retention'ı 15/27=%55.6'dır. Böylece RSF'deki feature-kümesi kaybı CoxPH'den daha hafif, fakat yönü aynı winner's-curse desenindedir.

![N=20 ve N=50 olgunlaşma karşılaştırması](results/phase56_comparison/n20_n50_maturation.png)

## Dürüstlük kontrolleri

Formül SHA256 değerleri bütün adımlarda sabit kaldı. Yeni N=50 çıktılarında kimlik-sızıntısı kontrolü geçti. Query havuzu genişletilmedi; sonuçlar sabit 12+12 query ve saklanmamış permütasyon-yolu belirsizliği sınırlamasını taşır. ADIM 4 analizi 0.353 sn sürdü.
