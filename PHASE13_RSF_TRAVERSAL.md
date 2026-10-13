# Faz 13 — RSF ağaç gezinme ve marjinal etki testi

**Statü:** diagnostic/exploratory. Komut: `python experiments/61_phase13_rsf_traversal.py`. Wall-clock **29.439 sn**.

## ADIM 0 — Kurulum

Fold-0 ve fold-1 modellerinde 500 ağaç, `min_samples_split=10`, `min_samples_leaf=5`, `max_features=sqrt`, `random_state=42` doğrulandı. Açıklanan skaler çıktı ve bu analizde kullanılan çıktı, `RandomSurvivalForest.predict()` risk skorudur. Ağaç yönü sklearn/sksurv uygulamasındaki `x_j <= threshold` sol dal kuralıdır.

Tek-replika attribution sağlaması: **ATLANDI: Saklı checkpoint'ler yalnız candidate/reference dizilerini içeriyor; bunları üreten bootstrap RSF model nesnesi saklanmamış. Yasaklanan yeniden fit yapılmadan aynı seed ile tekrar üretim mümkün değildir.** Fold-0 baseline üzerinde yeni attribution üretmek saklı bootstrap tensörüyle aynı modeli karşılaştırmayacağı için yapılmadı.

## ADIM 1–3 — Sonuç

Fold-0'da Tablo 8'de raporlanan 12 aktif özelliğin ağaç split sayıları birebir eşleşti; kalan 15 aktif özellik için sayımlar doğrudan ormandan çıkarıldı. Ham sorgu/split oranları, her background satırındaki yol ziyaretleri ve iki uç koalisyonun risk-skoru farkları ayrı CSV'lerdedir.

Özel özellikler:

| feature | warrant_type | TCGA_M_j | METABRIC_M_j | split_trees | split_nodes | base_value | TCGA_query_unique | METABRIC_query_unique | TCGA_separating_max | METABRIC_separating_max | TCGA_visit_full_mean | METABRIC_visit_full_mean | TCGA_max_abs_delta_last | TCGA_max_abs_delta_first | METABRIC_max_abs_delta_last | METABRIC_max_abs_delta_first |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| KRTAP19-7 | degenerate | 0 | 0 | 0 | 0 | 0.00154033 | 1 | 12 | 0 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 |
| OR6C1 | degenerate | 0 | 0.00416427 | 10 | 10 | 0.00178377 | 1 | 12 | 0 | 0.3 | 0.0140833 | 0.01425 | 1.06581e-14 | 3.55271e-15 | 0.0683333 | 0.0683333 |
| OR4C46 | degenerate | 0 | 0 | 0 | 0 | 0.000989955 | 1 | 12 | 0 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 |
| RBMY1F | degenerate | 0 | 0 | 0 | 0 | 0.000334774 | 1 | 12 | 0 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 3.55271e-15 | 7.10543e-15 |
| OR13C3 | degenerate | 0 | 0.00250771 | 0 | 0 | 0.00154807 | 1 | 12 | 0 | 0 | 0 | 0 | 1.42109e-14 | 7.10543e-15 | 7.10543e-15 | 2.13163e-14 |
| KRTAP19-2 | degenerate | 0 | 0 | 0 | 0 | 0.00159483 | 1 | 12 | 0 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 |
| ACCSL | informative | 0.00455625 | 0.00689027 | 20 | 21 | 0.0266457 | 2 | 12 | 0.333333 | 0.333333 | 0.0135 | 0.013375 | 0.006 | 0.006 | 0.07 | 0.07 |
| OR5M10 | degenerate | 0 | 0 | 0 | 0 | 0.00177985 | 1 | 12 | 0 | 0 | 0 | 0 | 3.55271e-15 | 1.42109e-14 | 7.10543e-15 | 2.13163e-14 |
| AJCC Stage I | informative | 0.000238564 | 0.000756506 | 2 | 2 | 0.169139 | 2 | 2 | 1 | 1 | 0.000208333 | 0.00025 | 7.10543e-15 | 7.10543e-15 | 1.42109e-14 | 2.13163e-14 |
| AJCC Stage II | informative | 0.00243141 | 0.0032892 | 10 | 10 | 0.55638 | 2 | 2 | 1 | 1 | 0.0015 | 0.0015 | 0.036 | 0.036 | 0.036 | 0.036 |
| AJCC Stage III | informative | 0.00688896 | 0.00337458 | 19 | 19 | 0.255193 | 2 | 2 | 1 | 1 | 0.00504167 | 0.00516667 | 0.092 | 0.092 | 0.00983333 | 0.0378333 |
| AJCC Stage IV | degenerate | 0 | 0 | 25 | 25 | 0.0133531 | 1 | 1 | 0 | 0 | 0.0457083 | 0.0457083 | 1.42109e-14 | 7.10543e-15 | 1.42109e-14 | 3.55271e-15 |

TCGA M_j=0 olan **11** RSF özelliğinin **11** tanesinde sorgu ve baz bütün split eşiklerinin aynı tarafında, **6** tanesinde iki yol durumunda ziyaret oranı sıfır, **11** tanesinde iki uç koalisyonun tüm 12×16 hücrelerinde marjinal etki sıfırdır. H1'in aynı-taraf veya hiç-ziyaret koşulu **11/11** özellikte desteklenmiştir. Dışarıda kalanlar: `yok`. Uç-koalisyon testi tam Shapley hesabı değildir.

Fold-1 baseline saklı olduğu için iki analog da aynı analizden geçirildi:

| feature | warrant_type | TCGA_M_j | METABRIC_M_j | split_trees | split_nodes | base_value | TCGA_query_unique | METABRIC_query_unique | TCGA_separating_max | METABRIC_separating_max | TCGA_visit_full_mean | METABRIC_visit_full_mean | TCGA_max_abs_delta_last | TCGA_max_abs_delta_first | METABRIC_max_abs_delta_last | METABRIC_max_abs_delta_first |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| OR13C3 | NA | 0 | 0 | 0 | 0 | 0.00187875 | 1 | 12 | 0 | 0 | 0 | 0 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 | 7.10543e-15 |
| AJCC Stage IV | NA | 0 | 0 | 60 | 60 | 0.0178042 | 1 | 1 | 0 | 0 | 0.109333 | 0.108625 | 1.42109e-14 | 7.10543e-15 | 1.42109e-14 | 1.42109e-14 |

Fold-1 sonuçları ham ADIM 1–3 CSV'lerinde `fold=1` ile işaretlidir; 27 satırlık `PHASE13_FEATURE_SUMMARY.csv` yalnız resmi fold-0 panelini içerir.

## Bütünlük

Başlangıç SHA256: `{"components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa", "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68", "certification.py": "76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487"}`. Bitiş SHA256: `{"components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa", "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68", "certification.py": "76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487"}`. Değerler aynıdır; yeni model fit, bootstrap veya attribution tensörü üretilmedi.

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
