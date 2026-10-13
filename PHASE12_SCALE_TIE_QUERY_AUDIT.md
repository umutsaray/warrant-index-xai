# Faz 12 — Ölçek, tie ve sorgu-kohort duyarlılık denetimi

**Statü:** diagnostic/exploratory; resmi sonuç değildir. Başlangıç/bitiş hash'leri eşittir: `{"components.py": "7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa", "aggregation.py": "10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68", "certification.py": "76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487"}`. Yeni model fit, patient bootstrap veya explainer çağrısı yapılmadı. Toplam wall-clock: **86.018 sn**.

Kullanılan komut: `python experiments/60_phase12_scale_tie_query_audit.py`. Doğrulama komutu: `..\.venv\Scripts\python.exe -m pytest -q --basetemp=.pytest_phase12_tmp`.

## ADIM 1 — Ölçek

Kohort-0 medyan sıfır-olmayan attribution ölçekleri: `{"CoxPH-v2": 0.023934215222165295, "RSF": 0.01707069437342852, "DeepSurv": 0.012307440008774696}`. Orijinal `a=b=0.5, u=3, eta=0.02` skorları Faz 7.7 tablosuyla `atol=1e-9` düzeyinde birebir üretildi. H–L aralıkları attribution importance birimindedir; bu nedenle `u=6ks_m` aynı birim sistemini korur. Saklı inferans girdilerini ölçek-parametreli recursive test için yeniden tanımlamak bu tanısal kapsamda yapılmadı; BH/BY opsiyonel k=1 satırı atlandı.

| model | k | point | WI_min | WI_median | WI_max | bottleneck_counts |
| --- | --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | 0.25 | 0 | 0 | 0 | 0.415714 | {"G_claim": 26, "S": 1} |
| CoxPH-v2 | 0.5 | 0 | 0 | 0 | 0.660649 | {"G_claim": 26, "S": 1} |
| CoxPH-v2 | 1 | 1 | 0 | 0.254344 | 0.826099 | {"G_claim": 25, "S": 2} |
| CoxPH-v2 | 2 | 3 | 0 | 0.576891 | 0.829058 | {"G_claim": 22, "R": 2, "S": 3} |
| CoxPH-v2 | 4 | 6 | 0.0391485 | 0.709427 | 0.829911 | {"G_claim": 11, "R": 14, "S": 2} |
| DeepSurv | 0.25 | 0 | 0 | 0 | 0.280273 | {"F": 3, "G_claim": 22, "S": 2} |
| DeepSurv | 0.5 | 0 | 0 | 0 | 0.467982 | {"F": 2, "G_claim": 21, "S": 4} |
| DeepSurv | 1 | 0 | 0 | 0 | 0.658687 | {"F": 2, "G_claim": 20, "S": 5} |
| DeepSurv | 2 | 0 | 0 | 0.190974 | 0.790408 | {"F": 1, "G_claim": 19, "R": 1, "S": 6} |
| DeepSurv | 4 | 0 | 0 | 0.568412 | 0.790408 | {"F": 1, "G_claim": 16, "R": 5, "S": 5} |
| RSF | 0.25 | 7 | 0 | 0.121121 | 0.927174 | {"G_claim": 15, "R": 6, "S": 5, "U": 1} |
| RSF | 0.5 | 7 | 0 | 0.403585 | 0.927174 | {"G_claim": 15, "R": 6, "S": 5, "U": 1} |
| RSF | 1 | 10 | 0 | 0.60193 | 0.927174 | {"G_claim": 16, "R": 7, "S": 4} |
| RSF | 2 | 13 | 0 | 0.786542 | 0.927174 | {"G_claim": 16, "R": 8, "S": 2, "U": 1} |
| RSF | 4 | 16 | 0 | 0.812261 | 0.927834 | {"G_claim": 10, "R": 16, "S": 1} |

Her bileşenin dağılımı ve `score<0.8` sayısı:

| model | k | component | minimum | median | maximum | below_tau |
| --- | --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | 0.25 | F | 1 | 1 | 1 | 0 |
| CoxPH-v2 | 0.25 | G_claim | 0 | 0 | 0.793352 | 27 |
| CoxPH-v2 | 0.25 | R | 0.672379 | 0.76741 | 0.977772 | 17 |
| CoxPH-v2 | 0.25 | S | 0.0688296 | 0.342999 | 0.985971 | 23 |
| CoxPH-v2 | 0.25 | U | 0.00495472 | 0.625249 | 1 | 21 |
| CoxPH-v2 | 0.5 | F | 1 | 1 | 1 | 0 |
| CoxPH-v2 | 0.5 | G_claim | 0 | 0 | 0.877235 | 26 |
| CoxPH-v2 | 0.5 | R | 0.672379 | 0.76741 | 0.977772 | 17 |
| CoxPH-v2 | 0.5 | S | 0.135213 | 0.551511 | 0.989799 | 21 |
| CoxPH-v2 | 0.5 | U | 0.0224774 | 0.811779 | 1 | 12 |
| CoxPH-v2 | 1 | F | 1 | 1 | 1 | 0 |
| CoxPH-v2 | 1 | G_claim | 0 | 0.259027 | 0.919177 | 25 |
| CoxPH-v2 | 1 | R | 0.672379 | 0.76741 | 0.977772 | 17 |
| CoxPH-v2 | 1 | S | 0.254344 | 0.747187 | 0.994873 | 19 |
| CoxPH-v2 | 1 | U | 0.111775 | 0.90589 | 1 | 10 |
| CoxPH-v2 | 2 | F | 1 | 1 | 1 | 0 |
| CoxPH-v2 | 2 | G_claim | 0 | 0.581534 | 0.940148 | 23 |
| CoxPH-v2 | 2 | R | 0.672379 | 0.76741 | 0.977772 | 17 |
| CoxPH-v2 | 2 | S | 0.473047 | 0.873969 | 0.998718 | 6 |
| CoxPH-v2 | 2 | U | 0.413685 | 0.952945 | 1 | 2 |
| CoxPH-v2 | 4 | F | 1 | 1 | 1 | 0 |
| CoxPH-v2 | 4 | G_claim | 0.0391485 | 0.784332 | 0.954309 | 15 |
| CoxPH-v2 | 4 | R | 0.672379 | 0.76741 | 0.977772 | 17 |
| CoxPH-v2 | 4 | S | 0.769647 | 0.945473 | 0.99968 | 3 |
| CoxPH-v2 | 4 | U | 0.695433 | 0.976472 | 1 | 1 |
| DeepSurv | 0.25 | F | 0 | 0.285638 | 0.598168 | 27 |
| DeepSurv | 0.25 | G_claim | 0 | 0 | 0.610975 | 27 |
| DeepSurv | 0.25 | R | 0.655463 | 0.730094 | 0.790408 | 27 |
| DeepSurv | 0.25 | S | 0.0212698 | 0.126929 | 0.711516 | 27 |
| DeepSurv | 0.25 | U | 0.00918984 | 0.304394 | 0.976011 | 26 |
| DeepSurv | 0.5 | F | 0.000506022 | 0.462503 | 0.769429 | 27 |
| DeepSurv | 0.5 | G_claim | 0 | 0 | 0.805488 | 26 |
| DeepSurv | 0.5 | R | 0.655463 | 0.730094 | 0.790408 | 27 |
| DeepSurv | 0.5 | S | 0.0435637 | 0.214762 | 0.881249 | 26 |
| DeepSurv | 0.5 | U | 0.0257362 | 0.497571 | 0.988006 | 26 |
| DeepSurv | 1 | F | 0.00886908 | 0.659961 | 0.905684 | 18 |
| DeepSurv | 1 | G_claim | 0 | 0 | 0.902744 | 26 |
| DeepSurv | 1 | R | 0.655463 | 0.730094 | 0.790408 | 27 |
| DeepSurv | 1 | S | 0.0892742 | 0.332041 | 0.967685 | 26 |
| DeepSurv | 1 | U | 0.0576134 | 0.682312 | 0.994003 | 25 |
| DeepSurv | 2 | F | 0.103999 | 0.840703 | 0.964764 | 9 |
| DeepSurv | 2 | G_claim | 0 | 0.227655 | 0.951372 | 24 |
| DeepSurv | 2 | R | 0.655463 | 0.730094 | 0.790408 | 27 |
| DeepSurv | 2 | S | 0.172122 | 0.487653 | 0.991921 | 26 |
| DeepSurv | 2 | U | 0.202002 | 0.819145 | 0.997001 | 12 |
| DeepSurv | 4 | F | 0.351742 | 0.932453 | 0.989478 | 4 |
| DeepSurv | 4 | G_claim | 0 | 0.609777 | 0.975686 | 20 |
| DeepSurv | 4 | R | 0.655463 | 0.730094 | 0.790408 | 27 |
| DeepSurv | 4 | S | 0.338941 | 0.664935 | 0.99798 | 24 |
| DeepSurv | 4 | U | 0.434273 | 0.888171 | 0.998501 | 2 |
| RSF | 0.25 | F | 0.899538 | 0.994313 | 1 | 0 |
| RSF | 0.25 | G_claim | 0 | 0.176511 | 1 | 19 |
| RSF | 0.25 | R | 0.726641 | 0.903752 | 0.927834 | 8 |
| RSF | 0.25 | S | 0.0288509 | 0.740008 | 1 | 15 |
| RSF | 0.25 | U | 0 | 0.559334 | 1 | 20 |
| RSF | 0.5 | F | 0.94622 | 0.996522 | 1 | 0 |
| RSF | 0.5 | G_claim | 0 | 0.512114 | 1 | 18 |
| RSF | 0.5 | R | 0.726641 | 0.903752 | 0.927834 | 8 |
| RSF | 0.5 | S | 0.0530002 | 0.78743 | 1 | 15 |
| RSF | 0.5 | U | 0 | 0.707404 | 1 | 16 |
| RSF | 1 | F | 0.97595 | 0.99913 | 1 | 0 |
| RSF | 1 | G_claim | 0 | 0.732607 | 1 | 16 |
| RSF | 1 | R | 0.726641 | 0.903752 | 0.927834 | 8 |
| RSF | 1 | S | 0.101362 | 0.855771 | 1 | 12 |
| RSF | 1 | U | 0.00684676 | 0.79308 | 1 | 14 |
| RSF | 2 | F | 0.991071 | 0.999783 | 1 | 0 |
| RSF | 2 | G_claim | 0 | 0.853456 | 1 | 13 |
| RSF | 2 | R | 0.726641 | 0.903752 | 0.927834 | 8 |
| RSF | 2 | S | 0.197675 | 0.88446 | 1 | 10 |
| RSF | 2 | U | 0.0209684 | 0.880837 | 1 | 9 |
| RSF | 4 | F | 0.997379 | 0.999946 | 1 | 0 |
| RSF | 4 | G_claim | 0 | 0.926728 | 1 | 10 |
| RSF | 4 | R | 0.726641 | 0.903752 | 0.927834 | 8 |
| RSF | 4 | S | 0.348563 | 0.918616 | 1 | 5 |
| RSF | 4 | U | 0.114539 | 0.937301 | 1 | 6 |

Ortak ölçeğe alma model farkını ortadan kaldırmadı: k=4'te point sayıları CoxPH-v2/RSF/DeepSurv için sırasıyla 6/16/0'dır. DeepSurv tüm ızgarada sıfır kalırken RSF en yüksek, CoxPH ara düzeydedir. Bağlayıcı çoğunluk CoxPH-v2 ve RSF'de k=4'te G_claim'den R'ye geçerken DeepSurv'de bütün k değerlerinde G_claim kalır. Önceden belirlenen değerlendirme kuralına göre **R baskınlığı ölçeğe duyarlıdır**.

## ADIM 2 — Tie ve R

Tie oranı, her hedef özellik için tüm site×replika×rakip karşılaştırmalarında `o=0` payıdır. `R_noties`, bütün site ve replikalarda daima tie kalan rakipleri hem paydan hem paydadan çıkarır. Mutlak η ızgarası Faz 7.4 değerlerini `atol=1e-12` ile yeniden üretti.

| model | grid | parameter | eta | R_mean | R_noties_mean | tie_rate_mean | R_bottlenecks | point | spearman_tie_R |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | absolute | 0.01 | 0.01 | 0.73428 | 0.733025 | 0.285134 | 24 | 7 | 0.0042735 |
| CoxPH-v2 | absolute | 0.02 | 0.02 | 0.774072 | 0.771452 | 0.473234 | 24 | 9 | 0.0409035 |
| CoxPH-v2 | absolute | 0.05 | 0.05 | 0.894185 | 0.882229 | 0.787109 | 23 | 23 | 0.798535 |
| CoxPH-v2 | absolute | 0.1 | 0.1 | 0.967868 | 0.93151 | 0.945954 | 4 | 25 | 0.995116 |
| CoxPH-v2 | normalized | 0.25 | 0.00598355 | 0.723007 | 0.722637 | 0.190007 | 24 | 7 | -0.00671551 |
| CoxPH-v2 | normalized | 0.5 | 0.0119671 | 0.740824 | 0.739239 | 0.326758 | 24 | 8 | -0.0128205 |
| CoxPH-v2 | normalized | 1 | 0.0239342 | 0.792302 | 0.789478 | 0.532146 | 23 | 13 | 0.0903541 |
| CoxPH-v2 | normalized | 2 | 0.0478684 | 0.887895 | 0.877184 | 0.772579 | 24 | 23 | 0.746032 |
| DeepSurv | absolute | 0.01 | 0.01 | 0.666095 | 0.666095 | 0.485698 | 27 | 0 | 0.039072 |
| DeepSurv | absolute | 0.02 | 0.02 | 0.726376 | 0.726376 | 0.644046 | 27 | 0 | 0.611111 |
| DeepSurv | absolute | 0.05 | 0.05 | 0.841315 | 0.839638 | 0.838917 | 27 | 21 | 0.927328 |
| DeepSurv | absolute | 0.1 | 0.1 | 0.919133 | 0.913966 | 0.936154 | 25 | 27 | 0.952629 |
| DeepSurv | normalized | 0.25 | 0.00307686 | 0.625599 | 0.625599 | 0.246353 | 27 | 1 | -0.339795 |
| DeepSurv | normalized | 0.5 | 0.00615372 | 0.639991 | 0.639991 | 0.381481 | 27 | 1 | -0.263777 |
| DeepSurv | normalized | 1 | 0.0123074 | 0.682108 | 0.682108 | 0.533761 | 27 | 0 | 0.206349 |
| DeepSurv | normalized | 2 | 0.0246149 | 0.748951 | 0.748951 | 0.685755 | 27 | 2 | 0.766906 |
| RSF | absolute | 0.01 | 0.01 | 0.836229 | 0.829556 | 0.470114 | 25 | 15 | 0.586381 |
| RSF | absolute | 0.02 | 0.02 | 0.855049 | 0.836346 | 0.60547 | 24 | 17 | 0.60559 |
| RSF | absolute | 0.05 | 0.05 | 0.890792 | 0.86787 | 0.763846 | 25 | 21 | 0.857541 |
| RSF | absolute | 0.1 | 0.1 | 0.91581 | 0.88598 | 0.874786 | 25 | 22 | 0.979626 |
| RSF | normalized | 0.25 | 0.00426767 | 0.825138 | 0.818256 | 0.330399 | 23 | 14 | 0.367128 |
| RSF | normalized | 0.5 | 0.00853535 | 0.832409 | 0.825406 | 0.442536 | 25 | 17 | 0.526856 |
| RSF | normalized | 1 | 0.0170707 | 0.848756 | 0.83861 | 0.574188 | 24 | 15 | 0.6482 |
| RSF | normalized | 2 | 0.0341414 | 0.869933 | 0.848511 | 0.7 | 25 | 18 | 0.89646 |

η=0.02 için M_j=0 gruplaması:

| model | M_zero | features | tie_rate_mean | R_mean | R_noties_mean |
| --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | False | 27 | 0.473234 | 0.774072 | 0.771452 |
| DeepSurv | False | 18 | 0.61812 | 0.729087 | 0.729087 |
| DeepSurv | True | 9 | 0.695897 | 0.720954 | 0.720954 |
| RSF | False | 16 | 0.520649 | 0.833465 | 0.821793 |
| RSF | True | 11 | 0.728846 | 0.886442 | 0.857515 |

Kalıcı tie'ları çıkarmak ortalama R'yi genel olarak düşürmektedir; bu, tie'ların R'yi özellikle RSF'de yukarı çektiğini, fakat tek başına bütün modeli açıklamadığını gösterir. Mutlak η büyüdükçe tie oranı ve R birlikte yükselir; CoxPH'de η=0.10'da R-darboğaz sayısının 4'e düşmesi bu modelde güçlü η duyarlılığına işaret eder. DeepSurv'de R bütün normalize ızgarada 27/27 darboğazdır.

## ADIM 3 — Sorgu kohortu ve background eşitliği

Eşitlik oranı, her sorgu değerinin ilgili açıklayıcı background değerlerinden en az biriyle tam (`atol=1e-12`) eşleşme payıdır. RSF k=16 sampled inactive background kullanırken aktif koordinatlar train ortalamasına sabitlenmiştir; DeepSurv deterministik nested k=16/k=128 ham train background kullanır. Fold-1 satırları mevcut saklı tensörlerden ayrıca verilmiştir.

| model | feature | feature_id | cohort | n_unique | variance | equality_rate | M_j | splits | warrant |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | METABRIC | 1 | 1.69271e-36 | 0 | 0.00119698 | nan | N/A (not certified) |
| DeepSurv | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | KRTAP19-2 | rna_ENSG00000186965.5 | METABRIC | 12 | 0.000459148 | 0 | 0.0588894 | nan | N/A (not certified) |
| DeepSurv | KRTAP19-2 | rna_ENSG00000186965.5 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | KRTAP19-7 | rna_ENSG00000244362.3 | METABRIC | 12 | 0.000347073 | 0 | 0.0430824 | nan | N/A (not certified) |
| DeepSurv | KRTAP19-7 | rna_ENSG00000244362.3 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | OR13C3 | rna_ENSG00000204246.3 | METABRIC | 12 | 0.000395267 | 0 | 0.0269285 | nan | N/A (not certified) |
| DeepSurv | OR13C3 | rna_ENSG00000204246.3 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | OR4C46 | rna_ENSG00000185926.1 | METABRIC | 12 | 5.90828e-05 | 0 | 0.0164054 | nan | N/A (not certified) |
| DeepSurv | OR4C46 | rna_ENSG00000185926.1 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | OR5M10 | rna_ENSG00000254834.4 | METABRIC | 12 | 0.000378812 | 0 | 0.027147 | nan | N/A (not certified) |
| DeepSurv | OR5M10 | rna_ENSG00000254834.4 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| DeepSurv | OR6C1 | rna_ENSG00000205330.4 | METABRIC | 12 | 8.02464e-05 | 0 | 0.0246895 | nan | N/A (not certified) |
| DeepSurv | OR6C1 | rna_ENSG00000205330.4 | TCGA | 1 | 0 | 0.960938 | 0.0203029 | nan | N/A (not certified) |
| DeepSurv | RBMY1F | rna_ENSG00000169800.14 | METABRIC | 12 | 3.85222e-05 | 0 | 0.0199141 | nan | N/A (not certified) |
| DeepSurv | RBMY1F | rna_ENSG00000169800.14 | TCGA | 1 | 0 | 1 | 0 | nan | N/A (not certified) |
| RSF | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | METABRIC | 1 | 1.69271e-36 | 0 | 0 | 25 | degenerate |
| RSF | AJCC Stage IV | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV | TCGA | 1 | 0 | 0 | 0 | 25 | degenerate |
| RSF | KRTAP19-2 | rna_ENSG00000186965.5 | METABRIC | 12 | 0.000459148 | 0 | 0 | 0 | degenerate |
| RSF | KRTAP19-2 | rna_ENSG00000186965.5 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |
| RSF | KRTAP19-7 | rna_ENSG00000244362.3 | METABRIC | 12 | 0.000347073 | 0 | 0 | 0 | degenerate |
| RSF | KRTAP19-7 | rna_ENSG00000244362.3 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |
| RSF | OR13C3 | rna_ENSG00000204246.3 | METABRIC | 12 | 0.000395267 | 0 | 0.00250771 | 0 | degenerate |
| RSF | OR13C3 | rna_ENSG00000204246.3 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |
| RSF | OR4C46 | rna_ENSG00000185926.1 | METABRIC | 12 | 5.90828e-05 | 0 | 0 | 0 | degenerate |
| RSF | OR4C46 | rna_ENSG00000185926.1 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |
| RSF | OR5M10 | rna_ENSG00000254834.4 | METABRIC | 12 | 0.000378812 | 0 | 0 | 0 | degenerate |
| RSF | OR5M10 | rna_ENSG00000254834.4 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |
| RSF | OR6C1 | rna_ENSG00000205330.4 | METABRIC | 12 | 8.02464e-05 | 0 | 0.00416427 | 10 | degenerate |
| RSF | OR6C1 | rna_ENSG00000205330.4 | TCGA | 1 | 0 | 0 | 0 | 10 | degenerate |
| RSF | RBMY1F | rna_ENSG00000169800.14 | METABRIC | 12 | 3.85222e-05 | 0 | 0 | 0 | degenerate |
| RSF | RBMY1F | rna_ENSG00000169800.14 | TCGA | 1 | 0 | 0 | 0 | 0 | degenerate |

Tam 27-feature sayımı: **RSF: TCGA M_j=0 olan 11 özelliğin 11'inde TCGA sorguları sabit; DeepSurv: TCGA M_j=0 olan 9 özelliğin 9'inde TCGA sorguları sabit**. TCGA sorgularındaki sıfır/sabit değerler DeepSurv'ün gerçek hasta background'ında sıklıkla eşleşirken RSF'nin aktif koordinatları train ortalamasına sabitlenmiş background'ında eşleşmez. Buna rağmen bazı sıfır-varyanslı özelliklerin DeepSurv attribution büyüklüğü sıfır değil ve RSF'de OR6C1/AJCC Stage IV split kullanımı vardır. Bu nedenle örtüşme dejenere warrant'larla ilişkili bir mekanizma adayıdır, tek başına yeterli açıklama değildir.

## ADIM 4 — Aggregator duyarlılığı

Yalnız nokta tahminleri karşılaştırıldı; aritmetik/geometrik ortalamalara sertifikasyon uygulanmadı.

| model | min_gt | arithmetic_gt | geometric_gt | rho_min_arithmetic | rho_min_geometric | top5_overlap_min_arithmetic | top5_overlap_min_geometric |
| --- | --- | --- | --- | --- | --- | --- | --- |
| CoxPH-v2 | 9 | 27 | 27 | 0.840049 | 0.857143 | 3 | 3 |
| DeepSurv | 0 | 27 | 27 | 0.879731 | 0.884615 | 5 | 5 |
| RSF | 17 | 27 | 27 | 0.894477 | 0.900031 | 4 | 4 |

Min agregatörü ile aritmetik/geometrik ortalamalar arasındaki sıralama korelasyonu yaklaşık 0.84–0.90 olsa da `>0.8` kararları keskin biçimde ayrışır: kompansatuvar iki alternatif üç modelde de 27/27 özellik geçirir. Sonuçların zayıf-halka min kuralına maddi ölçüde bağlı olduğu böylece doğrulanır.

## Dosyalar ve kapsam

Ham özellik düzeyi sonuçlar `PHASE12_ADIM1_SCALE_GRID.csv`, `PHASE12_ADIM2_TIE_R.csv`, `PHASE12_ADIM3_QUERY_COHORT.csv` ve `PHASE12_ADIM4_AGGREGATORS.csv` içindedir. Aşamaların süreleri sırasıyla **11.390/74.164/0.254/0.007 sn**. Formül kaynaklarına dokunulmadı.

ADIM 1 sonucu: k=4'te CoxPH-v2/RSF/DeepSurv point sayıları 6/16/0 olduğundan ortak ölçek model farkını kaldırmadı. ADIM 2 sonucu: η=0.02'de ortalama R→R_noties CoxPH-v2/RSF/DeepSurv için 0.774→0.771, 0.855→0.836 ve 0.726→0.726 oldu; tie etkisi en büyük RSF'dedir. ADIM 3 sonucu: RSF: TCGA M_j=0 olan 11 özelliğin 11'inde TCGA sorguları sabit; DeepSurv: TCGA M_j=0 olan 9 özelliğin 9'inde TCGA sorguları sabit ve background eşleşmesi tek başına dejenere warrant'ları açıklamadı. ADIM 4 sonucu: min ile point sayıları 9/17/0 iken iki kompansatuvar agregatör her modelde 27/27 verdi.

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.

Ölçek, tie ve sorgu-kohort denetiminin gösterdiği yönsel eğilim: model farkı ölçek normalizasyonuyla kaybolmamakta, tie katkısı en çok RSF'de görülmekte, kohort-içi sabitlik tek başına dejenere warrant'ları açıklamamakta ve kararlar min agregatörüne güçlü biçimde bağlı kalmaktadır.
