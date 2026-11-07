# Faz 15 — Faz 14 sonrası kapanış tanıları

Bu paket yalnız tanısal/keşifseldir. Yeni model eğitilmemiş, özellik paneli yeniden seçilmemiş ve mevcut resmî `results/**` artefaktlarının üzerine yazılmamıştır. Yeni ham çıktılar yalnız `results/phase15/` altındadır.

## ADIM 1 — DeepSurv örtüşmesiz validasyon

Beş fold için eğitim kodundaki `%15`, `random_state=42`, olay-stratifiye satır ayrımı birebir yeniden kuruldu. Validasyon satır/hasta/olay sayıları ve checkpoint validasyon/test C-index'leri Faz 14 değerleriyle `atol=1e-6` eşleşti.

| Fold | Tüm val C | Örtüşmesiz satır/hasta/olay | Örtüşmesiz C (%95 hasta-bootstrap GA) | Örtüşen C | Test C |
|---:|---:|---|---:|---:|---:|
| 0 | 0.8528 | 87/86/9 | 0.8480 [0.6642,0.9544] | 0.7568 | 0.4965 |
| 1 | 0.8237 | 85/84/12 | 0.8576 [0.7370,0.9509] | 0.4545 | 0.7192 |
| 2 | 0.8854 | 83/83/10 | 0.9029 [0.8065,0.9663] | 0.8750 | 0.6500 |
| 3 | 0.8553 | 82/82/7 | 0.8584 [0.7026,0.9706] | 0.8873 | 0.5319 |
| 4 | 0.9635 | 86/84/10 | 0.9672 [0.9033,1.0000] | 0.9474 | 0.6628 |

Tüm-validasyon ortalaması 0.8761, örtüşmesiz ortalama 0.8868, test ortalaması 0.6121'dir. Ortalama val–test farkı 0.2641'den 0.2747'ye yükselmiştir. Epoch-bazlı C-index günlüğü checkpoint'lerde saklanmadığından seçilen epoch/maksimum/ortalama/son-10 karşılaştırması yapılamadı.

**Sonuç:** Örtüşen hastalar çıkarılınca val C-index ortalaması 0.8761→0.8868 yükseliyor; test ile fark 0.2641→0.2747 oluyor. Hasta örtüşmesi farkın ana mekanik açıklaması olarak desteklenmemiştir; yalnız 7–12 örtüşmesiz olay ve checkpoint'in aynı validasyonda seçilmesi önemli kalan sınırlamalardır.

## ADIM 2 — Seyrek RNA kaldıraç tanısı

Seçim kodu her özelliği kendi train ortalama/SD'siyle z-skorlayıp özel Breslow-tie tek-değişkenli Cox Newton döngüsüyle fit eder; sıralama iki yanlı Wald p-değerine, eşitlikte partial log-likelihood'a dayanır. Seçilen 200 RNA'nın resmî β ve p-değerleri kayıtlı dosyayla `rtol=1e-6` eşleşti.

Kaldıraç eşdeğeri olarak en büyük-|z| satırının toplam `sum(z²)` içindeki payı raporlandı. LRT, aynı Breslow partial likelihood'un bağımsız skaler optimizasyonuyla; skor testi β=0'daki skor/bilgiyle; winsorized Wald ise z değerleri önceden sabit ±5'e kırpılarak hesaplandı.

| Gen | Nonzero | maks |z| | maks-z² payı | Wald p / sıra | LRT p / sıra | Skor p / sıra | ±5 Wald p / sıra |
|---|---:|---:|---:|---:|---:|---:|---:|
| RBMY1F | 2 | 24.85 | 91.6% | 1.36e-55 / 61 | 0.0240 / 83 | 7.64e-7 / 28 | 0.0782 / 130 |
| OR6C1 | 22 | 22.57 | 75.6% | 1.84e-78 / 35 | 0.0347 / 114 | 7.95e-4 / 147 | 0.1507 / 144 |
| OR4C46 | 6 | 18.99 | 53.5% | 4.81e-63 / 51 | 0.0311 / 103 | 2.04e-4 / 112 | 0.00776 / 63 |
| OR13C3 | 8 | 23.72 | 83.4% | 1.15e-43 / 85 | 0.0285 / 98 | 3.26e-5 / 69 | 0.0297 / 90 |
| OR5M10 | 6 | 21.85 | 70.8% | 4.89e-11 / 197 | 0.0856 / 189 | 0.00518 / 187 | 0.5344 / 182 |
| KRTAP19-2 | 3 | 24.32 | 87.8% | 2.82e-17 / 174 | 0.0474 / 136 | 9.10e-5 / 85 | 0.3553 / 170 |
| KRTAP19-7 | 4 | 24.48 | 88.9% | 1.05e-118 / 15 | 0.0729 / 173 | 0.00248 / 175 | 0.4006 / 171 |
| KRTAP11-1 | 48 | 20.57 | 62.8% | 1.07e-37 / 102 | 0.0460 / 132 | 0.00346 / 185 | 0.5159 / 178 |

Havuz talimat gereği zaten seçilmiş 200 RNA ile sınırlı olduğundan her özellik her alternatif sıralamada teknik olarak ilk 200 içinde kalır; anlamlı tanı sıra hareketidir, yeni panel seçimi değildir. Nonzero satırların olay ve süre dizileri ham CSV'dedir.

**Sonuç:** Önceden yazılmış hipotez büyük ölçüde desteklenmiştir: resmî Wald p'leri LRT/skor/winsorized p'lerden birçok mertebe küçüktür ve odak genlerin çoğu alternatif sıralamalarda belirgin geriler; ancak RBMY1F skor sırası 61→28 ve KRTAP19-2 skor sırası 174→85 olduğundan etki her kriterde tek yönlü değildir.

## ADIM 3 — METABRIC eşlemesi ve standardizasyon

Faz 14'teki 11/27 eşleme veri eksikliği değildi: tanı betiği `data/external/metabric` eski cache'ine yönelmişti; resmî panel `data/external/metabric_phase72` kullanır. Doğru cache ve özgün `external_columns` fonksiyonuyla 27/27 aktif özellik eşlendi. Yeniden kurulan kohort-standardize matris saklı `external_full` ile `atol=1e-8`, yeniden hesaplanan resmî dış atıf/G/G_signed ise saklı değerlerle `atol=1e-12` eşleşti. Yeni indirme gerekmedi.

| Yöntem | Point-adequate | Darboğaz | Medyan G_claim | Medyan WI |
|---|---:|---|---:|---:|
| Resmî kohort-içi standardizasyon | 9/27 | R:24, G:3 | 0.9583 | 0.7586 |
| Yalnız TCGA-train ölçeği | 3/27 | G:23, R:4 | 0.0000 | 0.0000 |
| Resmî yeniden-hesap kontrolü | 9/27 | R:24, G:3 | 0.9583 | 0.7586 |

“Yalnız TCGA” eksi resmî farkları: G magnitude min/medyan/maks = -0.9933/-0.9651/+0.0295; G_signed = -0.9913/-0.9289/+0.0696; G_claim = -0.9913/-0.9489/+0.0696; WI = -0.8299/-0.7482/+0.0600.

**Sonuç:** CoxPH G ve WI sonuçları kohort-içi konum/ölçek hizalamasına güçlü biçimde duyarlıdır; standardizasyon kaldırıldığında point-adequate sayı 9→3 düşmektedir. Bu tanı resmî yöntemi değiştirmez. RSF/DeepSurv yeni atıf gerektirdiği için kapsam dışıdır.

## ADIM 4 — Yinelenen çift paydası

`recursive_studentized` içinde F/U/G için pivot `(sampled-point)/sampled_se`; S/R için basic Monte Carlo kökü `sampled-point >= point-τ` kullanılır ve site/dal p'leri maksimumla birleştirilir. Yinelenen sıfır-kayıplı çiftlerin `sampled` skorunu yukarı taşıması bootstrap farkını ve dolayısıyla p-değerini büyütür; resmî yöntem daha muhafazakâr yöndedir.

Önceden sabit `outer=1000`, `inner=499`, `seed=20261108` küçültülmüş tanıda:

| Bileşen/senaryo | α | Resmî ret | Farklı-indeks ret |
|---|---:|---:|---:|
| S single_active null | .01/.05/.10 | .003/.034/.058 | .009/.050/.085 |
| S tie null | .01/.05/.10 | .000/.001/.001 | .000/.001/.003 |
| R merkez .85, iki senaryo | .01/.05/.10 | 1.000/1.000/1.000 | 1.000/1.000/1.000 |
| R merkez .90, iki senaryo | .01/.05/.10 | 1.000/1.000/1.000 | 1.000/1.000/1.000 |

Exact %95 binom GA'ları ham CSV'dedir. Bu resmî 2000×999 kapısı değildir. Gerçek-veri düzeltilmiş recursive-studentized p/BH/BY hesabı yapılmadı: farklı payda için her iç çekimde yeniden studentization ve bütün yaprak/IUT ailesi gerekir; mevcut ara değerler saklı değildir ve süre tahmini 60 dakika sınırını aşar.

**Sonuç:** Yinelenenleri çıkarmak single-active S null'unda muhafazakârlığı azaltıp α=.05 ret oranını .034→.050 getirirken nominali aşmadı; tie null hâlâ çok muhafazakârdır. R=.85/.90 alternatiflerinde güç her iki paydada 1.000'dır.

## ADIM 5 — KRTAP11-1 METABRIC 38/36 artığı

Ayrıştıran split bulunup 8-yollu saklı attribution'ı sıfır olan replikalar 3 ve 5'tir. Her iki orman deterministik yeniden üretildi; 128-yollu reference tekrar hesaplandı ve saklı reference ile `atol=1e-12` eşleşti.

| Replika | Özelliği kullanan split | Ayrıştıran ağaç | Ayrıştırılan query | max |candidate-8| | max |reference-128| |
|---:|---:|---:|---:|---:|---:|
| 3 | 4 | 1 | 2/12 | 0 | 0 |
| 5 | 1 | 1 | 4/12 | 0 | 0 |

Replika 3'te ayrışan query değerleri 0.111055 ve 0.095133; replika 5'te 0.111055, 0.047366, 0.095133 ve 0.048892; ortak baz 0.013671'dir. Split varlığı tek başına marjinal tahmin farkını garanti etmez: ağaç katkısı koalisyon durumlarında etkisiz kalabilir veya diğer ağaçlarla tam iptal olabilir.

**Sonuç:** Sıfır atıf 8-yol örnekleme gürültüsüyle uyumsuzdur; 128-yollu referans da tam sıfırdır. 38 split-replika ile 36 nonzero-attribution replika farkı, iki ormanda split varlığının net marjinal etkiye dönüşmemesinden kaynaklanır.

## Komutlar ve bütünlük

```text
..\.venv\Scripts\python.exe experiments/67_phase15_deepsurv_val_nooverlap.py
..\.venv\Scripts\python.exe experiments/68_phase15_sparse_leverage.py
..\.venv\Scripts\python.exe experiments/69_phase15_krtap_residual.py
..\.venv\Scripts\python.exe experiments/70_phase15_metabric_standardization.py
..\.venv\Scripts\python.exe experiments/71_phase15_duplicate_sim.py
```

Başlangıç ve bitiş formül SHA256 değerleri aynıdır:

| Dosya | SHA256 |
|---|---|
| `components.py` | `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa` |
| `aggregation.py` | `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` |
| `certification.py` | `76ab49716cc3ead6dec71656bbd60c6d71c1d76cab8a83b2f0662af6d6a05487` |

Bu analiz dondurulmuş ayarları, Tablo 3–9'u, BH/BY kararlarını, saklı atıf tensörlerini, öznitelik panelini ve CoxPH-v2 retry protokolünü DEĞİŞTİRMEMEKTEDİR.
