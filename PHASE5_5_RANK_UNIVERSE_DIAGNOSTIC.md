# Faz 5.5 — R sıralama-evreni boyutu tanısı

Bu çalışma ayrı bir duyarlılık hesabıdır. CoxPH N=50'nin saklı katsayıları ve 141 satırlık sabit TCGA audit kümesi kullanıldı; yeni bootstrap, model fit/tahmin, explainer veya BH/BY çalıştırılmadı. Resmî 614-evrenli R ve 8/77 sertifikasyon sonucu değişmedi.

## Sonuç

**614 vs 77 sıralama-evreni hipotezi desteklenmedi.** CoxPH evrenini 77 feature'a indirmek aynı 22 nokta-geçen feature'da R'yi ortalama yalnız **0.00283** değiştirdi (medyan 0.00325; aralık -0.00994–0.00983). Yüksek WI-marjinli feature sayısı `>0.03` düzeyinde 9 olarak kaldı; RSF'de 16'dır.

Aynı 22 feature'ın R marjinlerinde RSF'ye Wasserstein uzaklığı **0.05857 → 0.06076** oldu; yakınlaşma yerine %3.7 uzaklaşma görüldü. Tüm 77 eşlenmiş feature'da uzaklık **0.02174 → 0.02381** oldu. Bu yön, evren boyutunun ana açıklama olduğu hipoteziyle uyumlu değildir.

## Aynı 22 CoxPH nokta-geçen feature

| feature | R_614 | R_77 | delta_R_77_minus_614 | WI_614 | WI_77 | WI_margin_614 | WI_margin_77 | RSF_R_77 | RSF_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000080910.13 | 0.81458 | 0.81904 | 0.00446 | 0.81458 | 0.81904 | 0.01458 | 0.01904 | 0.69969 | 0.69969 |
| rna_ENSG00000133742.14 | 0.84117 | 0.84131 | 0.00014 | 0.84117 | 0.84131 | 0.04117 | 0.04131 | 0.86981 | 0.45451 |
| rna_ENSG00000169800.14 | 0.83309 | 0.82315 | -0.00994 | 0.83309 | 0.82315 | 0.03309 | 0.02315 | 0.85648 | 0.85648 |
| rna_ENSG00000167654.18 | 0.81891 | 0.81669 | -0.00222 | 0.81891 | 0.81669 | 0.01891 | 0.01669 | 0.71984 | 0.71984 |
| rna_ENSG00000270547.6 | 0.82003 | 0.82213 | 0.00210 | 0.82003 | 0.82213 | 0.02003 | 0.02213 | 0.68833 | 0.68833 |
| rna_ENSG00000133800.9 | 0.83671 | 0.83968 | 0.00297 | 0.83671 | 0.83968 | 0.03671 | 0.03968 | 0.82593 | 0.82593 |
| rna_ENSG00000077616.11 | 0.81177 | 0.81294 | 0.00118 | 0.81177 | 0.81294 | 0.01177 | 0.01294 | 0.69346 | 0.69346 |
| rna_ENSG00000109158.11 | 0.83353 | 0.83023 | -0.00330 | 0.83353 | 0.83023 | 0.03353 | 0.03023 | 0.86309 | 0.86309 |
| rna_ENSG00000004846.16 | 0.80561 | 0.81544 | 0.00983 | 0.80561 | 0.81544 | 0.00561 | 0.01544 | 0.86963 | 0.86963 |
| rna_ENSG00000023892.11 | 0.80358 | 0.81038 | 0.00680 | 0.80358 | 0.81038 | 0.00358 | 0.01038 | 0.79557 | 0.79557 |
| rna_ENSG00000183747.12 | 0.80485 | 0.80286 | -0.00199 | 0.80485 | 0.80286 | 0.00485 | 0.00286 | 0.78071 | 0.78071 |
| rna_ENSG00000172901.20 | 0.82527 | 0.83271 | 0.00743 | 0.82527 | 0.83271 | 0.02527 | 0.03271 | 0.74765 | 0.74765 |
| rna_ENSG00000118407.15 | 0.80495 | 0.81454 | 0.00959 | 0.80495 | 0.81454 | 0.00495 | 0.01454 | 0.85059 | 0.85059 |
| rna_ENSG00000182173.13 | 0.82220 | 0.82526 | 0.00306 | 0.82220 | 0.82526 | 0.02220 | 0.02526 | 0.68390 | 0.68390 |
| rna_ENSG00000187758.8 | 0.83713 | 0.83820 | 0.00107 | 0.83713 | 0.83820 | 0.03713 | 0.03820 | 0.84311 | 0.84311 |
| rna_ENSG00000092010.15 | 0.81292 | 0.81798 | 0.00506 | 0.81292 | 0.81798 | 0.01292 | 0.01798 | 0.71357 | 0.71357 |
| rna_ENSG00000204246.3 | 0.83715 | 0.83105 | -0.00610 | 0.83715 | 0.83105 | 0.03715 | 0.03105 | 0.85648 | 0.85648 |
| rna_ENSG00000205186.4 | 0.82277 | 0.82860 | 0.00583 | 0.82277 | 0.82860 | 0.02277 | 0.02860 | 0.75596 | 0.75596 |
| rna_ENSG00000122025.15 | 0.80957 | 0.81514 | 0.00556 | 0.80957 | 0.81514 | 0.00957 | 0.01514 | 0.72808 | 0.72808 |
| clin_numeric__age_at_index.demographic | 0.96901 | 0.97245 | 0.00344 | 0.95035 | 0.95035 | 0.15035 | 0.15035 | 0.69664 | 0.69664 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | 0.86702 | 0.87487 | 0.00785 | 0.86702 | 0.87487 | 0.06702 | 0.07487 | 0.85828 | 0.85828 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.83076 | 0.84020 | 0.00944 | 0.83076 | 0.84020 | 0.03076 | 0.04020 | 0.85935 | 0.85935 |

## R marjin dağılımları

| distribution | n | mean | std | min | q25 | median | q75 | max | n_ge_0 | n_gt_0_03 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| same22_Cox_R_margin_614 | 22.00000 | 0.03012 | 0.03463 | 0.00358 | 0.01206 | 0.02248 | 0.03592 | 0.16901 | 22.00000 | 9.00000 |
| same22_Cox_R_margin_77 | 22.00000 | 0.03295 | 0.03462 | 0.00286 | 0.01575 | 0.02421 | 0.03682 | 0.17245 | 22.00000 | 9.00000 |
| same22_RSF_R_margin_77 | 22.00000 | -0.01563 | 0.07201 | -0.11610 | -0.08486 | -0.01186 | 0.05648 | 0.06981 | 10.00000 | 9.00000 |

77-evrenli CoxPH R marjininin medyanı 0.02421; özgün değer 0.02248, aynı feature'lardaki RSF medyanı -0.01186'tir. Feature-bazlı CoxPH–RSF Spearman rho da 0.283'ten 0.290'e yalnız çok az değişti ve zayıf kaldı.

## Karşı-olgusal nokta-geçen aileler

| distribution | n | mean | std | min | q25 | median | q75 | max | n_ge_0 | n_gt_0_03 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Cox_original_point_family_WI_margin | 22.00000 | 0.02927 | 0.03111 | 0.00358 | 0.01206 | 0.02248 | 0.03592 | 0.15035 | 22.00000 | 9.00000 |
| Cox_77_universe_point_family_WI_margin | 23.00000 | 0.03106 | 0.03005 | 0.00286 | 0.01529 | 0.02315 | 0.03545 | 0.15035 | 23.00000 | 9.00000 |
| RSF_actual_point_family_WI_margin | 23.00000 | 0.04396 | 0.02500 | 0.00017 | 0.02162 | 0.05648 | 0.06217 | 0.07026 | 23.00000 | 16.00000 |

Yalnız R'yi 77 evreninde hesaplayan tanısal WI, CoxPH'de 22 yerine 23 feature'ı eşiğin üstüne taşır. Yeni feature `rna_ENSG00000139874.6` için WI 0.79888→0.81170'tir. Bu resmî bir adequacy veya sertifikasyon sonucu değildir. Karşı-olgusal CoxPH ailesinin WI marjin medyanı 0.02315 olarak kalırken RSF medyanı 0.05648'tir; `>0.03` sayısı yine 9'a karşı 16'dır. Model-özel nokta-geçen ailelerin Wasserstein uzaklığı 0.02291→0.02332 ile azalmadı.

![R evreni karşılaştırması](results/phase55_rank_universe_diagnostic/rank_universe_comparison.png)

## Yorum

614 rakip yerine 77 rakip kullanmak CoxPH R düzeyini hafifçe yukarı kaydırıyor, fakat Faz 5.4'teki düşük/yüksek marjin dağılımı ayrımını kapatmıyor. Bu nedenle ters-sezgisel BH hayatta-kalma farkının ana kaynağı sıralama-evreni boyutu değildir. Kalan fark, modellerin replika boyunca ürettiği önem sıralamalarının yapısı ve feature-özel model/veri davranışıyla daha uyumludur; bu tanı tek başına bunu biyolojik fark olarak nitelemez.

Hesap 2.182 saniye sürdü. Dondurulmuş 614-evrenli R yeniden hesaplanarak kaynak tabloyla maksimum mutlak hata 0.0e+00 bulundu. Formül SHA256 değerleri değişmedi ve kimlik-sızıntısı kontrolü geçti.
