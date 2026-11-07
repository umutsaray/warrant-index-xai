# Faz 14 — Hakem 9 doğrulama paketi

Bu çalışma yalnız tanısal/keşifseldir. Resmî tablolar, saklı tensörler ve karar dosyaları okunmuş; mevcut `results/**` dosyalarının hiçbiri değiştirilmemiştir. Yeni ham çıktılar kökteki `PHASE14_ADIM*.csv` dosyalarıdır.

## ADIM 1 — CoxPH R sıralama evreni

`source_phi` boyutu `(50,153,610)` olarak doğrulandı. `rank_reproducibility(abs(source_phi).mean(axis=1), 0.02)` sonucu 27 aktif özellikte resmî `feature_report.csv` ile `atol=1e-12` düzeyinde birebir eşleşti.

| η | Evren | R min | R medyan | R maks | R darboğazı | Point-adequate |
|---:|---:|---:|---:|---:|---:|---:|
| 0.01 | 27 | 0.6651 | 0.7490 | 0.9541 | 25 | 8 |
| 0.01 | 610 | 0.6170 | 0.6956 | 0.9833 | 24 | 7 |
| **0.02** | **27** | **0.6998** | **0.7983** | **0.9505** | **24** | **11** |
| **0.02** | **610 (resmî)** | **0.6724** | **0.7674** | **0.9778** | **24** | **9** |
| 0.05 | 27 | 0.7648 | 0.9105 | 0.9330 | 24 | 23 |
| 0.05 | 610 | 0.7602 | 0.9234 | 0.9612 | 23 | 23 |
| 0.10 | 27 | 0.8109 | 0.9639 | 0.9744 | 9 | 26 |
| 0.10 | 610 | 0.7923 | 0.9916 | 0.9964 | 4 | 25 |

η=0.02'de `R_610` ile `R_27` Spearman korelasyonu 0.9542'dir; özellik-bazlı `R_27-R_610` farkı min/medyan/maks = -0.0273/+0.0300/+0.0430'dur. 100 evreninde beş tohum point sayısı 9–11, 300 evreninde 9–9'dur. Bütün özellik ve tohum değerleri `PHASE14_ADIM1_R_UNIVERSE.csv` içindedir.

**Sonuç:** R evreni 610→27 olunca Cox'ta R medyanı 0.7674→0.7983, point-adequate sayısı 9→11 olmuştur; model-arası R karşılaştırması evren boyutuna duyarlıdır.

## ADIM 2 — Yinelenen bootstrap replika çiftleri

Kod yolu, `indices = rng.integers(0, T, size=(B,T))` ile yerine koymalı çekim yapar. Pair-loss matrisi çekilen indekslerin sayım ağırlıklarıyla iki taraftan çarpılır; aynı orijinal replikanın çifti toplamda yer alır ve farkı sıfır olduğu için kaybı sıfırdır. Resmî payda yine `T(T-1)`'dir. Farklı-indeks tanısında payda `T²-sum(count²)` olarak ayarlanmıştır.

4,999 çekimde yinelenen sıralı çift oranı 0.019969±0.004071 olup teorik 1/50=0.020000 ile uyumludur.

| Model | Bileşen | Resmî ortalama−nokta | Farklı-indeks ortalama−nokta | Farklı−resmî |
|---|---|---:|---:|---:|
| CoxPH-v2 | S | +0.000201 | -0.000002 | -0.000203 |
| CoxPH-v2 | R | +0.004521 | +0.000010 | -0.004512 |
| RSF | S | +0.000512 | -0.000011 | -0.000523 |
| RSF | R | +0.002449 | -0.000455 | -0.002904 |
| DeepSurv | S | +0.000456 | -0.000193 | -0.000649 |
| DeepSurv | R | +0.005336 | -0.000130 | -0.005466 |

Yinelenen sıfır-kayıplı çiftler S ve özellikle R bootstrap ortalamasını sistematik olarak yukarı iter; farklı-indeks sürümü nokta tahmine yakındır. Feature-bazlı ortalama, SE ve farklar `PHASE14_ADIM2_PAIR_DUPLICATES.csv` içindedir. Ayrı studentized yaprak p-değerleri/BH/BY tanısı yapılmadı: mevcut çıktı çekirdeği farklı-paydalı dağılım için iç studentization saklamıyor; bunu üretmek yeni ve maliyetli iç-iç bootstrap gerektirirdi. Resmî p-değerleri değişmemiştir.

## ADIM 3 — RBMY1F ve düşük ifadeli gen seçimi

Seçim yalnız eğitim bölümünde tek değişkenli Cox ile yapılır ve RNA/miRNA/metilasyon bloklarından p-değerine göre ilk 200 alınır. `select_features.py` varyans/sıfır-oranı/seks-kromozomu filtresi uygulamaz; yalnız fit edilemeyen sütunları eler. Bu nedenle sıfıra yakın genler sıralamaya girebilir.

Fold-0'ın seçilen 200 RNA özelliğinin **166/200**'ünde eğitim sıfır oranı >%90'dır; aktif 22 RNA'nın **12/22**'si bu gruptadır. Tam liste ve fold-1 aktif RNA'ları dahil 270 satır `PHASE14_ADIM3_GENE_SELECTION.csv` içindedir. Yerel eşleme dosyası kromozom anotasyonunu genel olarak içermediği için kromozom alanı tahmin edilmemiştir; RBMY1F'nin yerel adı üzerinden Y bilgisi kaydedilmiştir.

| Gen | Train sıfır | Nonzero (E/K) | Test nonzero | Cox β | p | RNA sıra |
|---|---:|---:|---:|---:|---:|---:|
| RBMY1F | 99.70% | 1/1 | 0 | -0.9476 | 1.36e-55 | 61 |
| OR6C1 | 96.74% | 0/22 | 9 | -1.0163 | 1.84e-78 | 35 |
| OR4C46 | 99.11% | 0/6 | 2 | -1.3672 | 4.81e-63 | 51 |
| OR13C3 | 98.81% | 0/8 | 3 | -1.3333 | 1.15e-43 | 85 |
| OR5M10 | 99.11% | 0/6 | 5 | -0.4596 | 4.89e-11 | 197 |
| KRTAP19-2 | 99.55% | 0/3 | 2 | -1.3219 | 2.82e-17 | 174 |
| KRTAP19-7 | 99.41% | 0/4 | 3 | -1.0567 | 1.05e-118 | 15 |
| KRTAP11-1 | 92.88% | 0/48 | 11 | -0.8626 | 1.07e-37 | 102 |

**Sonuç:** “düşük ifade/az sayıda aykırı hasta” kısmı güçlü biçimde desteklenir; “Y-kromozomu/yalnız erkekler” açıklaması RBMY1F için doğrulanmaz, çünkü iki nonzero eğitim örneğinin biri erkek, biri kadındır.

## ADIM 4 — DeepSurv validasyonu ve C-index belirsizliği

Validasyon, dış fold eğitim satırlarının %15'inin `random_state=42` ile olay durumuna göre stratifiye edilmesiyle kurulur; hasta-gruplu değildir. Early stopping en yüksek validasyon C-index checkpoint'ini seçer, patience=30'dur. Dış train–test hasta örtüşmesi her fold'da 0'dır; fakat iç fit–validasyon hasta örtüşmesi 15,17,19,19,15'tir.

| Fold | Fit satır/hasta/olay | Val satır/hasta/olay | Test satır/hasta/olay | Val C | Test C | Epoch |
|---:|---|---|---|---:|---:|---:|
| 0 | 572/530/90 | 102/101/16 | 168/153/26 | 0.8528 | 0.4965 | 38 |
| 1 | 572/532/90 | 102/101/16 | 168/153/26 | 0.8237 | 0.7192 | 53 |
| 2 | 572/533/90 | 102/102/16 | 168/153/26 | 0.8854 | 0.6500 | 37 |
| 3 | 572/532/89 | 101/101/16 | 169/155/27 | 0.8553 | 0.5319 | 5 |
| 4 | 572/530/89 | 101/99/16 | 169/155/27 | 0.9635 | 0.6628 | 38 |

DeepSurv test C-index hasta-kümeli bootstrap %95 GA'ları sırasıyla [0.3918,0.6204], [0.5778,0.8381], [0.5090,0.7613], [0.3953,0.6575], [0.5427,0.7912]'dir (B=2,000; seed 20261107). CoxPH/RSF için fold 1–4 risk vektörleri veya fit edilmiş model nesneleri saklanmadığından yeni fit olmadan GA hesaplanamadı; CSV'de açık durum alanıyla işaretlendi. 5×3 test C-index ve örneklem sayıları `PHASE14_ADIM4_DEEPSURV_VALIDATION.csv` içindedir.

Epoch-bazlı tam validasyon eğrisi saklanmadığından maksimum-vs-ortalama doğrudan hesaplanamadı. Bununla birlikte checkpoint'in aynı validasyon setindeki maksimumla seçilmesi, yalnız 16 validasyon olayı ve 15–19 hasta örtüşmesi birlikte iyimserlik mekanizmasıyla uyumludur; geniş test GA'ları belirsizliğin de yüksek olduğunu gösterir. Bu sayılar farkın tümünü nedensel olarak kanıtlamaz.

## ADIM 5 — Null kalibrasyonu

Ham 2,000 dış-simülasyon hücrelerinin tamamı, exact %95 GA, tek yönlü binom p ve Holm-gate alanlarıyla `PHASE14_ADIM5_CALIBRATION_TYPE1.csv` dosyasındadır. Aşağıda iki null senaryosundaki maksimum ampirik Tip-I oranları verilmiştir:

| Model/bileşen | α=.01 | α=.05 | α=.10 |
|---|---:|---:|---:|
| Cox F / S / U / G / R | .0065/.0065/.0050/.0090/.0000 | .0425/.0275/.0440/.0420/.0005 | .0955/.0585/.0905/.0950/.0020 |
| RSF F / S / U / G / R | .0060/.0045/.0130/.0080/.0000 | .0480/.0295/.0475/.0530/.0000 | .1035/.0635/.0925/.1020/.0020 |
| Deep F / S / U / G / R | .0095/.0060/.0040/.0110/.0000 | .0530/.0250/.0365/.0525/.0000 | .1030/.0640/.0895/.1035/.0025 |

N=20 RSF remedy arşivinde `recursive_basic` (3 gate), `m_out_of_n_root` (7 gate) ve `recursive_studentized` (0 gate) saklıdır. `direct_root` adlı ayrı bir arşiv hücresi bulunmadığından bu etiketle yeniden tablo üretilemedi; tahmin edilmedi.

## ADIM 6 — Standardizasyon ve Cox G

Panel zincirinde METABRIC ham değerleri önce kendi kohort ortalama/SD'siyle z-skorlanıp TCGA train ortalama/SD ölçeğine taşınır; Cox betiğinde bu `external_full` daha sonra TCGA'da fit edilmiş model scaler'ından geçirilir. Net etki, METABRIC'in kohort-içi standardize edilmiş konumunu Cox ölçeğine taşımaktır; doğrudan ham METABRIC'i yalnız TCGA istatistikleriyle ölçeklemek değildir.

Mevcut yerel ham METABRIC cache aktif 27 özelliğin yalnız 11'ini yeniden eşleyebildi; 16 özellik eksik olduğu için üç-yollu G/G_signed/G_claim/WI karşılaştırması **yapılamadı**. Eksiklerin tam listesi `PHASE14_ADIM6_G_STANDARDIZATION.csv` içindedir. RSF/DeepSurv yeni atıf gerektirdiğinden kapsam dışıdır.

## ADIM 7 — RSF bootstrap ormanlarının yeniden üretimi

Önce replika 0–2, sonra 50 replikanın tamamı aynı patient-cluster örnekleri, RSF ayarları, background ve açıklama tohumlarıyla yeniden üretildi. 2 site × 50 replikanın **100/100**'ü saklı candidate tensörüyle `atol=1e-12` eşleşti; maksimum mutlak fark 1.7764e-15'tir. Ham karşılaştırma `PHASE14_ADIM7_BOOTSTRAP_REPRO.csv`, 2,700 traversal satırı `PHASE14_ADIM7_BOOTSTRAP_TRAVERSAL.csv` içindedir.

| Özellik/site | Ayrıştıran split bulunan replika | Nonzero saklı attribution replika | Toplam split düğümü |
|---|---:|---:|---:|
| OR13C3 / TCGA | 0/50 | 0/50 | 136 |
| OR13C3 / METABRIC | 24/50 | 24/50 | 136 |
| KRTAP11-1 / TCGA | 3/50 | 3/50 | 354 |
| KRTAP11-1 / METABRIC | 38/50 | 36/50 | 354 |
| AJCC Stage IV / TCGA | 0/50 | 0/50 | 875 |
| AJCC Stage IV / METABRIC | 0/50 | 0/50 | 875 |

TCGA'da saklı attribution'ı bütün 50 replikada sıfır olan 11 özelliğin tamamında ayrıştıran split de 0/50'dir. Faz 13 referans ormanında OR13C3 için iki sitede de ayrıştıran split yokken bootstrap ormanları METABRIC'teki nonzero attribution'ın 24 replika boyunca ayrıştıran splitlerle geldiğini gösterir. KRTAP11-1 için TCGA hipotezi de 3/50 replikada birebir doğrulanır. H2 desteklenmiştir.

## Komutlar, bütünlük ve süre

Başlıca komutlar:

```text
python experiments/62_phase14_fast_checks.py
python experiments/63_phase14_pair_duplicates.py
..\.venv\Scripts\python.exe experiments/64_phase14_validation_ci.py
..\.venv\Scripts\python.exe experiments/65_phase14_rsf_reproduce.py --replicates 3
..\.venv\Scripts\python.exe experiments/65_phase14_rsf_reproduce.py --replicates 50
..\.venv\Scripts\python.exe experiments/66_phase14_g_standardization.py
```

RSF 50-replika yeniden üretimi 239 saniye; kalan saklı-veri analizleri toplamı yaklaşık 2 dakika sürdü. Başlangıç ve bitiş SHA256 değerleri aynıdır:

| Dosya | SHA256 |
|---|---|
| `components.py` | `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa` |
| `aggregation.py` | `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` |
| `certification.py` | `76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487` |

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
