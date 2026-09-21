# Hasta-düzeyi veri sızıntısı teşhisi

**Tek cümlelik sonuç: Modeller örnek-düzeyinde ayrık fakat hasta-düzeyinde çakışmalı (kirli) split üzerinde eğitildi ve test edildi.**

Bu teşhis salt okunur yapıldı. Hiçbir model eğitilmedi, hiçbir artefaktın ağırlığı değiştirilmedi ve WI/certification koduna dokunulmadı. Raporda gerçek barcode veya hasta/örnek kimliği yoktur.

## ADIM 0 — Dosya ve kod zinciri

### Kullanılan veri ve model dosyaları

Üç model aynı seçilmiş matris çiftini kullandı:

- Train özellikleri: `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/X_train_selected.parquet` — 673 satır, 616 özellik.
- Test özellikleri: `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/X_test_selected.parquet` — 169 satır, 616 özellik.
- Train sonuçları: `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/y_train.parquet` — 673 satır.
- Test sonuçları: `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/y_test.parquet` — 169 satır.
- Modeller: `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/coxph_baseline.joblib`, `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/rsf_baseline.joblib`, `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/explainable_dl_model.pt`.

`C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:46-49` bu dört parquet dosyasını doğrudan okur. CoxPH train fit'i `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:64-81`, 169 satırlık test prediction/C-index hesabı `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:85-92` içindedir. RSF fit'i `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:97-108`, aynı test dosyasında prediction/C-index hesabı `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:109-113` içindedir. Sonuçlar `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/baseline_results.json:12` ve `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/baseline_results.json:27` satırlarına yazılmıştır.

DeepSurv aynı dört dosyayı `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:144-147` satırlarında okur; test tensörünü `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:175-181` arasında kurar, mevcut test riskini ve C-index'ini `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:227-228` satırlarında hesaplar ve checkpoint'e `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:232-246` arasında kaydeder. Checkpoint içindeki test C-index 0.7456922162804516'dır.

Dosya hashleri daha sonraki audit metadata'sıyla birebir eşleşir: `C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase3_brca/metadata.json` içindeki `input_sha256` CoxPH/veri dosyalarını; `C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase6_deepsurv_pilot/metadata.json` içindeki `input_sha256` DeepSurv/veri dosyalarını doğrular. Tam hash ve zaman damgaları `C:/Users/UMUT/Desktop/ExCred/warrantxai/results/leakage_diagnostic/diagnostic.json` içindedir.

### “23 hasta” kaynağı

Audit script'i aynı dosya çiftini okur: `C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/06_brca_real.py:225-228`. Hasta anahtarını indeksin ilk 12 karakteriyle kurar (`:238-239`), train'de bulunan hastaların test satırlarını atar (`:240-241`) ve kalan test havuzunda hasta başına tek satır tutar (`:242-243`). Yani “23” farklı veya önceden temizlenmiş bir dosya çiftinden gelmez; modellerin kullandığı aynı 673/169 dosyadan gelir. Bu filtre model eğitimine/test C-index hesabına geri uygulanmamış, yalnız sonraki XAI audit query havuzuna uygulanmıştır.

## ADIM 1 — Bağımsız hasta çakışması hesabı

Kaynak üretim kodu `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/prepare_dataset.py:19,29-31` ile `sample_key` değerini TCGA barcode'unun ilk **15** karakteri olarak tanımlar. `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/impute_and_split.py:203-213` yalnız bu 15-karakter örnek anahtarının benzersizliğini kontrol eder; `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/impute_and_split.py:221-230` aynı örnek anahtarlarını satır düzeyinde train/test'e böler ve yalnız tam örnek-anahtarı kesişimini yasaklar. TCGA hasta düzeyi ilk 12 karakterdir. Bu nedenle farklı örnek türlerine ait 15-karakter anahtarlar aynı hastayı iki split'e taşıyabilir.

Orijinal dosyalardan sıfırdan yapılan ilk-12-karakter karşılaştırması:

| Ölçü | Sayı |
|---|---:|
| Train satırı | 673 |
| Test satırı | 169 |
| Train benzersiz hasta | 628 |
| Test benzersiz hasta | 164 |
| Her iki split'te bulunan benzersiz hasta | **23** |
| Ortak hastalara ait train satırı | 25 |
| Ortak hastalara ait test satırı | 23 |
| Ortak hastalar çıkarılınca test satırı | 146 |
| Disjoint test benzersiz hasta | 141 |
| Orijinal / disjoint test olayı | 26 / 19 |

Bağımsız sonuç önceki 23 sayısıyla tam eşleşir. Audit'teki 141 sayısı, 23 ortak-hasta test satırının çıkarılmasından sonra kalan 146 satırda hasta başına tek örnek tutulmasıyla oluşur.

## ADIM 2 — İşlem sırası ve kirli/temiz karar

`C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/select_features.py:173-177,215-216` sample-level split'ten gelen X_train/X_test satırlarını koruyarak 616-feature dosyalarını yazar. CoxPH ve RSF bunları doğrudan `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_baselines.py:46-49` üzerinden; DeepSurv `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:144-147` üzerinden kullanır. Eğitim kodlarında ilk-12-karakter hasta filtresi yoktur.

Dosya zamanları bu sıralamayla uyumludur: seçilmiş train/test dosyaları 2026-08-27 21:12:32 UTC, saklı CoxPH/RSF/DeepSurv modelleri sırasıyla 21:19:26, 21:19:30 ve 21:19:40 UTC'dir. Hasta-düzeyi audit filtresi `C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/06_brca_real.py:238-243` içinde yer alır; script ve Faz 3 metadata'sı 2026-09-17 tarihlidir. Kardeş dizin bir Git deposu değildir; commit geçmişi yoktur, bu nedenle kanıt çalıştırma kodu, SHA256 eşleşmeleri ve dosya zamanlarıyla sınırlıdır.

Ek teşhis: DeepSurv erken-durdurma ayrımı da train satırları üzerinde `C:/Users/UMUT/Desktop/ExCred/Multi-Omics Data Fusion with Explainable Deep Learning for Survival Prediction in Breast Cancer/train_explainable_dl.py:166-172` ile yapılır. Aynı seed ve stratifikasyon yeniden kurulduğunda 572 fit ve 101 validation satırı arasında 12 ortak hasta vardır; bu hastalar her tarafta 12 satıra karşılık gelir. Bu da hasta-gruplu bir validation ayrımı değildir.

## ADIM 3 — Mevcut modellerde disjoint test alt-kümesi

Hiçbir model yeniden fit edilmedi. Saklı scaler/ağırlıklar aynen yüklenerek önce 169 satırlık orijinal test C-index'leri yeniden üretildi; ardından yalnız 23 ortak-hasta test satırı çıkarılarak kalan 146 satırda skorlandı.

| Model | Raporlanan C-index | Yeniden üretilen 169-satır | Disjoint 146-satır | Δ (disjoint−orijinal) |
|---|---:|---:|---:|---:|
| CoxPH | 0.706477 | 0.706477 | 0.700508 | -0.005969 |
| RSF | 0.714795 | 0.714795 | 0.729272 | +0.014477 |
| DeepSurv | 0.745692 | 0.745692 | 0.763959 | +0.018267 |

CoxPH farkı 0.01'in altında; RSF ve DeepSurv mutlak farkları 0.01–0.02 aralığındadır. Disjoint alt-kümede RSF ve DeepSurv skorları düşmek yerine yükselmiştir; dolayısıyla bu tanı, çakışmanın raporlanan C-index'leri yapay olarak yükselttiğine dair ampirik işaret vermez. Bununla birlikte bağımsız test varsayımı 169-satırlık raporlanan metrik için ihlal edilmiştir. Sayısal etki mütevazı olsa da metodolojik sorun makale açısından ihmal edilebilir değildir; tamamen hasta-ayrık performansı belirlemek, bu teşhisin kapsamı dışında kalan hasta-gruplu yeniden split/eğitim gerektirir.

## Bütünlük ve gizlilik kontrolleri

- Dört veri dosyası ve üç model yalnız okundu; model `.fit` çağrısı yapılmadı.
- Orijinal üç C-index 1e-12 toleransında yeniden üretildi.
- Faz 3 ve Faz 6 metadata SHA256 kayıtları mevcut dosyalarla eşleşti.
- `diagnostic.json` ve bu rapor gerçek barcode, patient ID veya sample ID içermez.
