# Faz 21 ADIM 0 — veri ve sağkalım tanımı

Bu adım salt okunur envanterdir. Mevcut `results/phase21/` eski Faz 2.1 artefaktlarını içerdiği için üzerine yazılmadı; yeni dosyalar `phase21_reviewer1/` altındadır.

## Sağkalım tanımı

Sonlanım **overall survival (OS)**; başlangıç tanıdır. TCGA `OS.time` gün, METABRIC `OS_MONTHS` ay birimindedir. Olay 1=ölüm, 0=sağ sansürdür.

| cohort | n | events | censor_rate | time_unit | median_followup_reverse_km |
|---|---|---|---|---|---|
| TCGA fold-0 train | 674 | 106 | 0.84273 | days | 984 |
| TCGA fold-0 test | 168 | 26 | 0.845238 | days | 1074 |
| METABRIC expression cohort | 1979 | 1143 | 0.422436 | months | 195.533 |

## TCGA örnekleri ve eksik veri

| category | n |
|---|---|
| primary solid tumor | 762 |
| other | 75 |
| metastatic | 5 |

| category | n |
|---|---|
| unique patients | 769 |
| patients with >1 retained sample | 71 |
| rows belonging to repeated patients | 144 |

Model matrisi aynı hastanın birden çok satırını dışlamaz; hasta-gruplu dış split bütün satırları aynı fold'da tutar. Fold-0 audit sorgu havuzu ayrıca hasta başına ilk deterministik satırı tutar. Birleşik 853 satırdan survival kontrol bayrağı taşıyan 11 satır çıkarılmıştır. Methylation KNN (k=5), klinik sayısal medyan, kategorik en-sık değer imputasyonu yalnız train'de fit edilip test'e uygulanmıştır; son eksik hücre sayısı sıfırdır.

Klinik ham eksiklik sayıları:

| category | n |
|---|---|
| age_at_index.demographic | 0 |
| gender.demographic | 0 |
| race.demographic | 0 |
| ethnicity.demographic | 0 |
| ajcc_pathologic_stage.diagnoses | 6 |

## DeepSurv epoch kayıtları

Epoch-bazlı validation/test eğrileri yalnız fold-0 Faz 7.0b tanı koşusunda saklanmıştır: satır-düzeyi geçmiş 68 epoch, hasta-gruplu tanı geçmişi 61 epoch içerir. Fold 1–4 için yalnız seçilen checkpoint/özet değerleri saklanmış, epoch-bazlı validation/test eğrileri saklanmamıştır. Mevcut fold-0 hasta-gruplu tanı kaydı en iyi validation C-index 0.942184, eşzamanlı test C-index 0.488716 ve sıfır fit/validation hasta kesişimi verir.

## Bootstrap model artefaktları

RSF ve DeepSurv için 50 bootstrap model nesnesi diskte saklı değildir. Kompakt checkpoint/atıf dizileri ve `bootstrap_summary.json` kayıtları vardır. Her iki koşu da seed `20260917`, girdi hash'leri ve protokol hash'leriyle deterministik olarak yeniden üretilebilir; yeniden üretim model fit gerektirir.
