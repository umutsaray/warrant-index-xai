# Faz 4.2 — RSF N=20 tam 77-feature koşu raporu

Onaylanan tasarım uygulandı: sampled background k=16, candidate/reference bütçesi 8/128, N=20 patient-cluster bootstrap ve pilotla aynı 12 TCGA + 12 METABRIC query. Ana `rsf_baseline.joblib` yeniden eğitilmedi; 20 bootstrap modeli ayrıca fit edildi. Tahmin maliyetini azaltmak için bootstrap modellerinde scikit-survival `low_memory=True` kullanıldı: aynı seed/veriyle 500 ağaç yapısı birebir aynı, standart ve low-memory risk skorları arasındaki maksimum fark 7.11e-15 bulundu. Bu bir model/metodoloji değişikliği değil, çıktı saklama optimizasyonudur.

## Block bazında agregat

| block | features | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| clin | 5 | 5 | 0.9999 | 0.9783 | 0.9860 | 0.9694 | 0.8271 | 0.8271 |
| rna | 72 | 72 | 0.9998 | 0.9508 | 0.9655 | 0.9081 | 0.7566 | 0.7320 |

77 feature değerlendirildi; 23 feature WI>0.8 koşulunu geçti. Bottleneck dağılımı `{'R': 72, 'G_claim': 5}`.

## N=2 pilot ve N=20 karşılaştırması

Eski sabit-background N=2 pilotunda 43/77 feature geçmişti; Faz 4.1 bu background'ı reddettiği için doğrudan replika etkisi karşılaştırması N=20 dizisinin ilk iki replikasıyla aynı sampled background altında yapıldı. Sampled-background N=2'de 35/77, N=20'de 23/77 geçti. N=2 geçenlerin 22 tanesi kaldı, 13 tanesi düştü; 1 feature yeni geçti. CoxPH Faz 3.2'deki winner's-curse deseniyle karşılaştırıldığında benzer bir eşik gerilemesi görüldü.

Feature bazındaki üçlü karşılaştırma `n2_n20_comparison.csv` içindedir.

## Sınırlamalar ve süre

Query kümesi yalnız 12+12 gerçek vektördür; CoxPH auditindeki 141 TCGA + 1466 METABRIC havuzundan çok daha dardır. Bu nedenle özellikle G/G_signed örnekleme oynaklığı daha yüksek olabilir. Sonuç tam 616-feature Shapley oyunu değildir: 77 ortak feature coalition'a girer, diğer 539 girdi 16 TCGA background hastası üzerinden entegre edilir. Toleranslar kalibre edilmemiştir; klinik geçerlilik iddiası yoktur.

Ölçülen birikimli eğitim+açıklama wall-clock süresi 2949.3 sn (0.82 saat). Completion invocation süresi 2270.8 sn; 1 site checkpoint'i yeniden kullanıldı. Formül modüllerinin Faz 2.1 SHA256 değerleri değişmedi. N=50 çalıştırılmadı.

## Tüm feature sonuçları

| feature | block | F | S | U | G_claim | R | WI | bottleneck |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rna_ENSG00000121351.8 | rna | 0.9999 | 0.9293 | 0.9794 | 0.8583 | 0.6860 | 0.6860 | R |
| rna_ENSG00000178522.14 | rna | 1.0000 | 0.9918 | 0.9886 | 0.9805 | 0.8002 | 0.8002 | R |
| rna_ENSG00000006116.4 | rna | 1.0000 | 0.9982 | 0.9936 | 0.9858 | 0.8631 | 0.8631 | R |
| rna_ENSG00000080910.13 | rna | 0.9999 | 0.9779 | 0.9871 | 0.9471 | 0.6997 | 0.6997 | R |
| rna_ENSG00000133742.14 | rna | 0.9995 | 0.7751 | 0.8360 | 0.4545 | 0.8698 | 0.4545 | G_claim |
| rna_ENSG00000113889.14 | rna | 0.9997 | 0.9429 | 0.9672 | 0.9434 | 0.7138 | 0.7138 | R |
| rna_ENSG00000182591.6 | rna | 1.0000 | 0.9991 | 0.9962 | 0.9884 | 0.8574 | 0.8574 | R |
| rna_ENSG00000169800.14 | rna | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.8565 | 0.8565 | R |
| rna_ENSG00000198670.12 | rna | 0.9999 | 0.9666 | 0.9712 | 0.9206 | 0.7147 | 0.7147 | R |
| rna_ENSG00000090402.8 | rna | 0.9977 | 0.5920 | 0.6834 | 0.0000 | 0.8260 | 0.0000 | G_claim |
| rna_ENSG00000167654.18 | rna | 0.9999 | 0.9868 | 0.9898 | 0.9717 | 0.7198 | 0.7198 | R |
| rna_ENSG00000233452.7 | rna | 0.9998 | 0.9627 | 0.9758 | 0.9041 | 0.7018 | 0.7018 | R |
| rna_ENSG00000139874.6 | rna | 0.9999 | 0.9633 | 0.9696 | 0.9442 | 0.6947 | 0.6947 | R |
| rna_ENSG00000198099.9 | rna | 0.9999 | 0.9851 | 0.9886 | 0.9824 | 0.8173 | 0.8173 | R |
| rna_ENSG00000169760.17 | rna | 0.9995 | 0.8886 | 0.8889 | 0.7365 | 0.8105 | 0.7365 | G_claim |
| rna_ENSG00000231861.2 | rna | 0.9998 | 0.8598 | 0.9213 | 0.8839 | 0.6905 | 0.6905 | R |
| rna_ENSG00000270069.1 | rna | 0.9994 | 0.9214 | 0.9482 | 0.8525 | 0.6976 | 0.6976 | R |
| rna_ENSG00000174899.12 | rna | 0.9999 | 0.9776 | 0.9780 | 0.9608 | 0.7515 | 0.7515 | R |
| rna_ENSG00000100911.16 | rna | 0.9999 | 0.9830 | 0.9845 | 0.9811 | 0.7242 | 0.7242 | R |
| rna_ENSG00000113389.16 | rna | 0.9998 | 0.9485 | 0.9481 | 0.9138 | 0.6896 | 0.6896 | R |
| rna_ENSG00000169432.18 | rna | 0.9996 | 0.7407 | 0.8436 | 0.5073 | 0.6931 | 0.5073 | G_claim |
| rna_ENSG00000157617.17 | rna | 0.9998 | 0.9564 | 0.9819 | 0.9495 | 0.6742 | 0.6742 | R |
| rna_ENSG00000158122.12 | rna | 0.9999 | 0.9411 | 0.9418 | 0.9014 | 0.7021 | 0.7021 | R |
| rna_ENSG00000168917.9 | rna | 0.9999 | 0.9819 | 0.9855 | 0.9686 | 0.7741 | 0.7741 | R |
| rna_ENSG00000270547.6 | rna | 0.9999 | 0.9609 | 0.9830 | 0.9359 | 0.6883 | 0.6883 | R |
| rna_ENSG00000133800.9 | rna | 0.9999 | 0.9821 | 0.9906 | 0.9794 | 0.8259 | 0.8259 | R |
| rna_ENSG00000077616.11 | rna | 1.0000 | 0.9679 | 0.9737 | 0.9456 | 0.6935 | 0.6935 | R |
| rna_ENSG00000226674.11 | rna | 0.9999 | 0.9106 | 0.9568 | 0.8815 | 0.7052 | 0.7052 | R |
| rna_ENSG00000064042.18 | rna | 0.9977 | 0.7536 | 0.8568 | 0.4337 | 0.7032 | 0.4337 | G_claim |
| rna_ENSG00000138795.10 | rna | 0.9999 | 0.9795 | 0.9845 | 0.9611 | 0.7011 | 0.7011 | R |
| rna_ENSG00000100908.14 | rna | 0.9997 | 0.9325 | 0.9398 | 0.8539 | 0.7920 | 0.7920 | R |
| rna_ENSG00000266964.6 | rna | 1.0000 | 0.9907 | 0.9905 | 0.9945 | 0.7924 | 0.7924 | R |
| rna_ENSG00000109158.11 | rna | 1.0000 | 0.9984 | 0.9957 | 0.9925 | 0.8631 | 0.8631 | R |
| rna_ENSG00000060709.15 | rna | 1.0000 | 0.9737 | 0.9830 | 0.9641 | 0.7077 | 0.7077 | R |
| rna_ENSG00000004846.16 | rna | 1.0000 | 0.9979 | 0.9952 | 0.9920 | 0.8696 | 0.8696 | R |
| rna_ENSG00000023892.11 | rna | 1.0000 | 0.9929 | 0.9913 | 0.9895 | 0.7956 | 0.7956 | R |
| rna_ENSG00000214575.9 | rna | 1.0000 | 0.9787 | 0.9867 | 0.9885 | 0.7629 | 0.7629 | R |
| rna_ENSG00000231925.12 | rna | 0.9997 | 0.9627 | 0.9599 | 0.8839 | 0.7576 | 0.7576 | R |
| rna_ENSG00000183747.12 | rna | 1.0000 | 0.9906 | 0.9880 | 0.9732 | 0.7807 | 0.7807 | R |
| rna_ENSG00000172901.20 | rna | 0.9997 | 0.7981 | 0.9405 | 0.7769 | 0.7476 | 0.7476 | R |
| rna_ENSG00000118407.15 | rna | 1.0000 | 0.9973 | 0.9933 | 0.9910 | 0.8506 | 0.8506 | R |
| rna_ENSG00000187531.14 | rna | 1.0000 | 0.9864 | 0.9871 | 0.9767 | 0.8095 | 0.8095 | R |
| rna_ENSG00000108639.8 | rna | 1.0000 | 0.9975 | 0.9941 | 0.9862 | 0.7980 | 0.7980 | R |
| rna_ENSG00000116663.11 | rna | 1.0000 | 0.9975 | 0.9957 | 0.9949 | 0.8613 | 0.8613 | R |
| rna_ENSG00000124343.14 | rna | 0.9998 | 0.9407 | 0.9278 | 0.8848 | 0.6948 | 0.6948 | R |
| rna_ENSG00000076043.11 | rna | 0.9998 | 0.9467 | 0.9256 | 0.8817 | 0.7389 | 0.7389 | R |
| rna_ENSG00000134531.10 | rna | 1.0000 | 0.9973 | 0.9939 | 0.9892 | 0.8545 | 0.8545 | R |
| rna_ENSG00000182173.13 | rna | 0.9999 | 0.9499 | 0.9565 | 0.9209 | 0.6839 | 0.6839 | R |
| rna_ENSG00000172725.14 | rna | 0.9999 | 0.9513 | 0.9733 | 0.9319 | 0.6957 | 0.6957 | R |
| rna_ENSG00000240108.1 | rna | 0.9999 | 0.9716 | 0.9748 | 0.9108 | 0.6834 | 0.6834 | R |
| rna_ENSG00000187758.8 | rna | 1.0000 | 0.9974 | 0.9958 | 0.9869 | 0.8431 | 0.8431 | R |
| rna_ENSG00000092010.15 | rna | 0.9998 | 0.9886 | 0.9843 | 0.9692 | 0.7136 | 0.7136 | R |
| rna_ENSG00000104327.7 | rna | 0.9999 | 0.9847 | 0.9870 | 0.9288 | 0.7086 | 0.7086 | R |
| rna_ENSG00000204246.3 | rna | 1.0000 | 0.9993 | 0.9986 | 0.9963 | 0.8565 | 0.8565 | R |
| rna_ENSG00000064961.19 | rna | 0.9999 | 0.9454 | 0.9638 | 0.9648 | 0.6952 | 0.6952 | R |
| rna_ENSG00000124440.16 | rna | 0.9999 | 0.9889 | 0.9831 | 0.9408 | 0.6987 | 0.6987 | R |
| rna_ENSG00000205186.4 | rna | 1.0000 | 0.9925 | 0.9905 | 0.9758 | 0.7560 | 0.7560 | R |
| rna_ENSG00000163145.13 | rna | 0.9999 | 0.9490 | 0.9680 | 0.8751 | 0.7632 | 0.7632 | R |
| rna_ENSG00000181788.4 | rna | 1.0000 | 0.9689 | 0.9829 | 0.9612 | 0.7246 | 0.7246 | R |
| rna_ENSG00000246430.7 | rna | 1.0000 | 0.9952 | 0.9930 | 0.9779 | 0.8022 | 0.8022 | R |
| rna_ENSG00000072657.9 | rna | 0.9998 | 0.9307 | 0.9601 | 0.8712 | 0.7006 | 0.7006 | R |
| rna_ENSG00000121871.4 | rna | 1.0000 | 0.9933 | 0.9931 | 0.9875 | 0.8024 | 0.8024 | R |
| rna_ENSG00000169764.16 | rna | 1.0000 | 0.9829 | 0.9879 | 0.9770 | 0.7447 | 0.7447 | R |
| rna_ENSG00000140332.17 | rna | 0.9999 | 0.9470 | 0.9578 | 0.9293 | 0.7123 | 0.7123 | R |
| rna_ENSG00000231672.8 | rna | 0.9997 | 0.8825 | 0.9564 | 0.9167 | 0.7036 | 0.7036 | R |
| rna_ENSG00000163545.11 | rna | 1.0000 | 0.9974 | 0.9949 | 0.9886 | 0.8703 | 0.8703 | R |
| rna_ENSG00000187242.6 | rna | 0.9998 | 0.8484 | 0.9307 | 0.7758 | 0.6899 | 0.6899 | R |
| rna_ENSG00000250334.6 | rna | 1.0000 | 0.9965 | 0.9938 | 0.9907 | 0.8652 | 0.8652 | R |
| rna_ENSG00000187486.5 | rna | 0.9999 | 0.9933 | 0.9927 | 0.9898 | 0.8036 | 0.8036 | R |
| rna_ENSG00000095203.14 | rna | 1.0000 | 0.9408 | 0.9494 | 0.8974 | 0.6852 | 0.6852 | R |
| rna_ENSG00000149557.14 | rna | 1.0000 | 0.9791 | 0.9781 | 0.9588 | 0.7240 | 0.7240 | R |
| rna_ENSG00000122025.15 | rna | 0.9999 | 0.9824 | 0.9885 | 0.9716 | 0.7281 | 0.7281 | R |
| clin_numeric__age_at_index.demographic | clin | 0.9998 | 0.8935 | 0.9373 | 0.8547 | 0.6966 | 0.6966 | R |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 1.0000 | 0.9999 | 0.9993 | 0.9997 | 0.8583 | 0.8583 | R |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II | clin | 1.0000 | 0.9987 | 0.9963 | 0.9966 | 0.8648 | 0.8648 | R |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | clin | 1.0000 | 0.9993 | 0.9970 | 0.9961 | 0.8593 | 0.8593 | R |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | clin | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 0.8565 | 0.8565 | R |
