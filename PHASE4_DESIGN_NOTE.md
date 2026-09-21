# Faz 4 — RSF approximate-F tasarım notu

## Karar ve model kaynağı

Kardeş projedeki mevcut `rsf_baseline.joblib` yeniden kullanıldı; ana model yeniden eğitilmedi. Model aynı 616-feature train/test paneline aittir ve 500 ağaç, min_samples_split=10, min_samples_leaf=5, max_features=sqrt, random_state=42 kullanır. Pilotun N=2 eğitim replikaları aynı ayarlarla hasta-kümesi bootstrap örneklerinde ayrıca fit edildi.

SHAP 0.52.0 `TreeExplainer`, scikit-survival 0.28.0 `RandomSurvivalForest` modelini desteklemedi: `InvalidModelError: Model type not yet supported by TreeExplainer: <class 'sksurv.ensemble.forest.RandomSurvivalForest'>`. Bu nedenle TreeSHAP seçilmedi.

Kaynaklar: [SHAP TreeExplainer belgesi](https://shap.readthedocs.io/en/latest/generated/shap.TreeExplainer.html), [açık scikit-survival RSF destek talebi](https://github.com/shap/shap/issues/3308), [scikit-survival RSF API](https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.ensemble.RandomSurvivalForest.html).

## Approximate-F tanımı

Candidate, her model/query için 8 rastgele permütasyonlu model-agnostik Shapley; bağımsız referans 128 permütasyonlu aynı kestiricidir. Referans bütçesi candidate'ın 16 katıdır ve seed akışları ayrıdır. F formülü değiştirilmedi. Açıklanan skaler çıktı scikit-survival `RandomSurvivalForest.predict()` risk skorudur. Bu, kapalı-form oracle veya algoritmik exact TreeSHAP değildir; F hem candidate hem sonlu-bütçeli referans Monte Carlo hatasını içerir.

Mevcut RSF 616 girdi isterken METABRIC yalnız 77 ortak feature sağlar. Eksik 539 modaliteyi hasta verisiymiş gibi doldurmak yerine açıklama oyunu açıkça 77-boyutlu bir model kesiti olarak tanımlandı: ortak feature'lar coalition ile değişir, diğer koordinatlar TCGA train ortalamasında sabit kalır. Bu pilot sonucu tam 616-feature RSF Shapley açıklaması değildir. METABRIC ortak feature'ları kendi kohort z-score'larından TCGA train mean/std ölçeğine taşındı.

## Pilot kapsamı ve sonuçlar

N=2, 72 RNA + 5 klinik feature ve kohort başına 12 gerçek, deterministik seçilmiş kimliksiz query vektörü kullanıldı. Kaynak havuzu 141, tam METABRIC havuzu 1466 hastadır. Küçük query örneği özellikle G/G_signed için yüksek örnekleme oynaklığı yaratabilir. Tüm F/S/U/G_claim/R/WI değerleri sonlu, [0,1] aralığında ve WI=min(F,S,U,G_claim,R) eşitliğini sağladı. 43/77 feature WI>0.8 koşulunu geçti; bottleneck dağılımı `{'R': 72, 'G_claim': 5}`. F aralığı 0.9947–1.0000; bu yüksek değer aynı kestiricinin 8 ve 128 permütasyonlu sürümlerini karşılaştırmanın beklenen sonucudur ve bağımsız açıklayıcı doğruluğunu kanıtlamaz.

### Block özeti

| block | features | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 0.9996 | 0.9911 | 0.9891 | 0.9739 | 0.8566 | 0.8566 |
| rna | 72 | 0.9998 | 0.9616 | 0.9674 | 0.9169 | 0.7707 | 0.7410 |

### Örnek feature sonuçları

| feature | F | S | U | G_claim | R | WI | bottleneck |
| --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000090402.8 | 0.9993 | 0.6672 | 0.5593 | 0.0000 | 0.9013 | 0.0000 | G_claim |
| rna_ENSG00000133742.14 | 0.9996 | 0.6618 | 0.7237 | 0.0531 | 0.8947 | 0.0531 | G_claim |
| rna_ENSG00000231861.2 | 0.9947 | 0.6745 | 0.8670 | 0.8005 | 0.3947 | 0.3947 | R |
| rna_ENSG00000169432.18 | 0.9991 | 0.6594 | 0.8957 | 0.6691 | 0.4079 | 0.4079 | R |
| rna_ENSG00000266964.6 | 1.0000 | 0.9969 | 0.9917 | 0.9945 | 0.9079 | 0.9079 | R |
| rna_ENSG00000108639.8 | 1.0000 | 0.9999 | 0.9977 | 0.9979 | 0.9079 | 0.9079 | R |
| rna_ENSG00000139874.6 | 1.0000 | 0.9988 | 0.9961 | 0.9917 | 0.9013 | 0.9013 | R |
| rna_ENSG00000187531.14 | 1.0000 | 0.9959 | 0.9945 | 0.9888 | 0.9013 | 0.9013 | R |

## Süre ve ölçekleme

N=2 pilot toplamı 176.7 sn; eğitim 5.6 sn, candidate açıklamalar 10.6 sn, yüksek-bütçeli referans 159.2 sn. Aynı 12+12 query ve bütçelerle doğrusal N=6 projeksiyonu 527.4 sn, N=50 projeksiyonu 4385.3 sn. Candidate açıklamayı tüm 141+1466 G query havuzuna genişleten kaba üst projeksiyon N=6 için 2617.1 sn, N=50 için 21799.3 sn. Bunlar doğrusal tek-makine projeksiyonlarıdır.

## Sürprizler, sınırlar ve tam koşu öncesi karar noktası

TreeSHAP yolu mevcut kütüphane/model çiftiyle uygulanamadı. Daha temel sınır, dış kohortun 616 RSF girdisinin yalnız 77'sini sağlamasıdır; seçilen 77-feature kesiti bilimsel olarak ayrıca onaylanmalıdır. N=2 ile R tek bir replika çiftine dayanır; 12+12 query ile U ve G hassastır. Toleranslar a=b=0.5,u=3,eta=0.02,tau=0.8 kalibre edilmemiştir. Pilot teknik uygulanabilirliği gösterir; klinik geçerlilik, exact F veya tam-model açıklaması göstermez. N=50 tam koşu başlatılmadı.

## Faz 4.1 — sağlamlık kontrolü

Aynı N=2 modeller, aynı 12+12 query, aynı seed akışı ve sekiz feature kullanıldı: üç düşük WI, üç yüksek WI, iki deterministik rastgele feature. Background duyarlılığında pasif 539 girdi 16 benzersiz TCGA train hastasından alındı ve her coalition model çıktısı bu hastalar üzerinden ortalandı.

### Background sabitleme duyarlılığı

| feature | selection_reason | candidate_attr_mean_abs_delta_sampled_vs_mean | reference_attr_mean_abs_delta_sampled_vs_mean | F_background_delta | G_claim_background_delta |
| --- | --- | --- | --- | --- | --- |
| rna_ENSG00000090402.8 | low_WI | 0.0685 | 0.0685 | 0.0003 | 0.0000 |
| rna_ENSG00000133742.14 | low_WI | 0.0676 | 0.0679 | 0.0001 | 0.2504 |
| rna_ENSG00000231861.2 | low_WI | 0.1561 | 0.1515 | 0.0039 | -0.0541 |
| rna_ENSG00000266964.6 | high_WI | 0.0096 | 0.0098 | 0.0000 | 0.0005 |
| rna_ENSG00000108639.8 | high_WI | 0.0109 | 0.0105 | 0.0000 | -0.0072 |
| rna_ENSG00000139874.6 | high_WI | 0.0062 | 0.0062 | -0.0001 | -0.0075 |
| rna_ENSG00000138795.10 | seeded_random | 0.0135 | 0.0137 | 0.0000 | -0.0079 |
| rna_ENSG00000231672.8 | seeded_random | 0.0246 | 0.0264 | -0.0001 | -0.0179 |

Maksimum |ΔF|=0.003897, maksimum |ΔG_claim|=0.250446. Attribution farkları risk-skoru birimindedir; F/G farkları [0,1] skor birimindedir.

### Referans bütçe yakınsaklığı

| feature | F_mean_bg_ref128 | F_mean_bg_ref256 | F_budget_delta_256_minus_128 | G_claim_budget_delta_256_minus_128 |
| --- | --- | --- | --- | --- |
| rna_ENSG00000090402.8 | 0.9993 | 0.9995 | 0.0002 | 0.0000 |
| rna_ENSG00000133742.14 | 0.9996 | 0.9996 | 0.0000 | 0.0000 |
| rna_ENSG00000231861.2 | 0.9947 | 0.9961 | 0.0014 | 0.0000 |
| rna_ENSG00000266964.6 | 1.0000 | 1.0000 | -0.0000 | 0.0000 |
| rna_ENSG00000108639.8 | 1.0000 | 1.0000 | -0.0000 | 0.0000 |
| rna_ENSG00000139874.6 | 1.0000 | 1.0000 | -0.0000 | 0.0000 |
| rna_ENSG00000138795.10 | 1.0000 | 0.9999 | -0.0000 | 0.0000 |
| rna_ENSG00000231672.8 | 0.9998 | 0.9999 | 0.0000 | 0.0000 |

Maksimum |ΔF|=0.001367. G_claim dondurulmuş tanımda candidate attribution'dan hesaplanır; candidate=8 sabitken yalnız referansı 128'den 256'ya çıkarmak G_claim'i değiştirmez, bu nedenle ΔG_claim=0'dır.

### Öneri

**Use sampled background (k=16) and 8/128 for N=50.** Karar eşiği her kontrolde maksimum mutlak skor farkı <0.01 idi. Ölçülen sampled/mean süre oranları candidate için 4.97×, referans için 3.55×. Bu oranlarla kaba N=50 maliyet projeksiyonu 15575.5 sn (4.33 saat); background maliyetinin k=16 ile tam doğrusal arttığı muhafazakâr üst projeksiyon 68039.8 sn (18.90 saat). N=50 bu kontrolde başlatılmadı; ayrıca onay bekleniyor.

## Faz 4.2 — N=20 sonucu

Onaylanan sampled-background k=16 ve 8/128 tasarımla 77 feature/N=20 koşusu tamamlandı. 23/77 feature WI>0.8 geçti; bottleneck dağılımı `{'R': 72, 'G_claim': 5}`. Ayrıntılar `PHASE4_FULL_RUN_REPORT.md` içinde. N=50 başlatılmadı.
