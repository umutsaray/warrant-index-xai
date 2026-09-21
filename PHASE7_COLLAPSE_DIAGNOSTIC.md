# Faz 7.0b — C-index çöküşü kök neden teşhisi

## Sonuç

Faz 7.0 düşüşü tek başına “hasta sızıntısı temizlenince gerçek performans ortaya çıktı” diye yorumlanamaz. Kanıtlar üç ayrı sorunu gösteriyor:

1. **Fold/pipeline oynaklığı:** RSF C-index'i beş geçerli hasta-gruplu fold'da 0.5385–0.7896 arasında değişti ve seçilen fold 0, 5/5 içinde en düşük sonucu verdi. Fold 0–2 satır ve event ölçütlerinde tamamen eşitti; fold 0 yalnız son fold-numarası bağlayıcısıyla seçilmişti.
2. **Feature-selection kararsızlığı:** Eski ve yeni 616 panellerin yalnız 333 özelliği ortaktır (%54.1; Jaccard %37.0). RNA ortaklığı 43/200'dür (%21.5). Beş grouped fold'un ikili panel ortaklığı ortalama 304.9/616'dır ve 277–328 arasında değişir.
3. **Model-spesifik eğitim sorunları:** CoxPH fold 0 fit'i yakınsamamıştı; dondurulmuş `penalizer=0.1` korunup yalnız optimizasyon retry uygulanınca uyarı kayboldu ve C-index 0.6311→0.6677 oldu. DeepSurv ise iç validation yükselirken dış testte kötüleşen güçlü overfitting gösterdi; hasta-gruplu iç validation bu uçurumu gidermedi.

Bu nedenle Faz 7.0 fold 0 skorları tek bir leakage-düzeltilmiş “gerçek performans” tahmini olarak güvenilir değildir. CoxPH çöküşünde yakınsamama ölçülebilir teknik etkendir; RSF'de arbitrar tie-break ile seçilen fold ve yeniden öğrenilen panel zinciri ana etkendir; DeepSurv'de bunlara ek olarak validation tabanlı checkpoint sinyali dış teste ters yönde hareket etmektedir. Feature seçiminin skor kaybındaki bağımsız nedensel payı bu teşhiste ayrıştırılmadı; kararsızlığı doğrudan kanıtlandı, fakat fold bileşimiyle etkileşmektedir.

## 1. Beş grouped fold

Split üretimi `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)` kullanır ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:258)). Faz 7.0b için eklenen `--test-fold` yalnız tanı amacıyla aynı beş adayı ayrı ayrı seçer ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:204)). Her fold'da preprocessing, train-only feature selection ve aynı CoxPH/RSF fit zinciri baştan çalıştırıldı ([30_phase70b_fold_sweep.py](C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/30_phase70b_fold_sweep.py:52)).

| Fold | Train/test | Train/test event | Cox özgün | Cox uyarı | Cox yakınsamış | RSF |
|---:|---:|---:|---:|---:|---:|---:|
| **0 (Faz 7.0)** | 674/168 | .1573/.1548 | 0.6311 | 1 | **0.6677** | **0.5385** |
| 1 | 674/168 | .1573/.1548 | 0.7059 | 0 | 0.7059 | 0.7896 |
| 2 | 674/168 | .1573/.1548 | 0.5971 | 1 | **0.6812** | 0.6362 |
| 3 | 673/169 | .1560/.1598 | 0.6551 | 0 | 0.6551 | 0.7005 |
| 4 | 673/169 | .1560/.1598 | 0.6502 | 0 | 0.6502 | 0.6507 |

Yakınsamış Cox ortalaması 0.6720, aralığı 0.6502–0.7059'dur. RSF ortalaması 0.6631, aralığı 0.5385–0.7896'dır. Fold 0 RSF'de beş fold'un en kötüsüdür ve ortalamanın 0.1246 altındadır. Üstelik fold 0, 1 ve 2 seçim ölçütlerinde aynı satır ve event dengesine sahipti; varsayılan kod eşitliği fold numarasıyla çözer ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:282)). Performans fold seçiminde kullanılmadı; bu metodolojik olarak doğru olsa da tek-fold tahmininin varyansı büyüktür.

Bu deney test alt-kümesinin intrinsik zorluğu ile her fold'da yeniden seçilen panelin etkisini birbirinden ayırmaz. Kanıtlanan sonuç, tüm train→selection→model zincirinin fold'a çok duyarlı olduğu ve deterministik tie-break'in RSF için uçta kötü bir fold'a düştüğüdür. Tam tablo [fold_scores.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/fold_scores.csv) içindedir. Fold taraması 1059.17 saniye sürdü.

## 2. Feature-selection kararlılığı

Seçim kodu her fold'un yalnız `X_train`/`y_train` verisini kullanır ([select_features.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/select_features.py:173), [rank çağrısı](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/select_features.py:186)). Bu leakage açısından doğrudur; sonuçların düşük kararlılığını ortadan kaldırmaz.

| Blok | Eski | Yeni | Ortak | Eskiye göre ortaklık | Jaccard |
|---|---:|---:|---:|---:|---:|
| RNA | 200 | 200 | **43** | **%21.5** | %12.0 |
| miRNA | 200 | 200 | 151 | %75.5 | %60.6 |
| Methylation | 200 | 200 | 123 | %61.5 | %44.4 |
| Clinical | 16 | 16 | 16 | %100 | %100 |
| **Toplam** | **616** | **616** | **333** | **%54.1** | **%37.0** |

Beş yeni grouped fold arasında da ikili toplam ortaklık yalnız 277–328 özellik, ortalama 304.9/616'dır. RNA ikili ortaklığı 22–88/200 aralığındadır. Bu, özellikle univariate RNA sıralamasının split-duyarlı olduğunu doğrular. Sayılar [old_new_feature_stability.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/old_new_feature_stability.csv) ve [summary.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/summary.json) içinde saklıdır.

## 3. CoxPH yakınsama tanısı

Ana eğitim kodu event/censored alt-gruplarında varyansı `≤1e-12` olan sütunları dışlar ([train_baselines.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/train_baselines.py:52)) ve `CoxPHFitter(penalizer=0.1)` kullanır ([train_baselines.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/train_baselines.py:72)). Fold 0'da dışlanan sütunlar:

- `rna_ENSG00000226474.1`
- `rna_ENSG00000200796.1`
- `rna_ENSG00000252433.1`
- `rna_ENSG00000202438.1`
- `mirna_hsa-mir-4681`
- `clin_categorical__race.demographic_american indian or alaska native`

Tam uyarı: “Newton-Raphson failed to converge sufficiently. Please see the following tips in the lifelines documentation: https://lifelines.readthedocs.io/en/latest/Examples.html#problems-with-convergence-in-the-cox-proportional-hazard-model”.

İlk fit'in standardize-kovaryat katsayıları patolojik büyüklükte değildir: β aralığı −0.3205–0.3633, |β| medyanı 0.0565, HR aralığı 0.7258–1.4381 ve |β|>1 katsayı sayısı sıfırdır. Yine de optimizer çözümü tamamlamamıştır. Aynı feature seti ve `penalizer=0.1` ile yalnız `max_steps=200, step_size=0.25` retry ([31_phase70b_cox_diagnostic.py](C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/31_phase70b_cox_diagnostic.py:66)) uyarısız yakınsadı:

| Ölçü | İlk fit | Retry |
|---|---:|---:|
| Test C-index | 0.631128 | **0.667704** |
| β aralığı | −0.3205–0.3633 | −0.1410–0.2425 |
| HR aralığı | 0.7258–1.4381 | 0.8684–1.2744 |
| Yakınsama uyarısı | Var | Yok |

Fold 2 de aynı uyarıyı verdi; retry skoru 0.5971→0.6812 oldu. Bu yüzden Faz 7.0 CoxPH skoru tamamlanmamış optimizasyon nedeniyle 0.0366 aşağıdaydı. Retry sonrası eski 0.7065 skoruna kalan fark 0.0388'dir ve beş-fold dağılımıyla uyumludur. Ayrıntılı katsayı tablosu [cox_largest_coefficients.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/cox_largest_coefficients.csv), retry kayıtları [cox_fold_retries.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/cox_fold_retries.json) içindedir.

## 4. DeepSurv validation–test uçurumu

Ara checkpoint'ler saklanmamıştı. Faz 7.0 deterministik CPU eğitimi 0.74 saniye ölçüldüğünden, aynı mimari/hiperparametre/seed ile yalnız tanı amaçlı yeniden üretildi. Eğitim kodu satır düzeyinde stratified `train_test_split` kullanır ([train_explainable_dl.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/train_explainable_dl.py:167)) ve en yüksek validation C-index checkpoint'ini seçer ([train_explainable_dl.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/train_explainable_dl.py:211)).

| Epoch | Train C-index | Validation C-index | Test C-index |
|---:|---:|---:|---:|
| 1 | 0.7655 | 0.7518 | 0.4918 |
| 10 | 0.9000 | 0.8038 | **0.5533** |
| 26 | 0.9583 | 0.8413 | 0.5105 |
| **38 (seçilen)** | **0.9757** | **0.8528** | **0.4965** |
| 50 | 0.9849 | 0.8268 | 0.5058 |
| 68 (durma) | 0.9918 | 0.8167 | 0.5253 |

Gözlenen en yüksek test C-index 0.5572 ile epoch 8'dedir; test verisi checkpoint seçmek için kullanılamayacağından bu resmi alternatif değildir. Validation ve test C-index'lerinin epoch-bazlı Spearman korelasyonu **−0.711**'dir. Checkpoint seçimi böylece 0.0607 test C-index kaybıyla ilişkilidir, fakat en iyi gözlenen 0.5572 bile eski 0.7457'nin çok altındadır; checkpoint tek başına çöküşü açıklamaz.

İçteki 15 ortak hastayı doğrudan sınamak için aynı eğitim bir kez hasta-gruplu iç validation ile tekrarlandı ([32_phase70b_deepsurv_diagnostic.py](C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/32_phase70b_deepsurv_diagnostic.py:37)). Hasta kesişimi sıfıra indiği halde en iyi validation C-index **0.9422**, bu checkpoint'in test C-index'i **0.4887** oldu; validation–test korelasyonu −0.862'ydi. Dolayısıyla 15 hasta-yayılımı dev uçurumun nedeni olarak desteklenmedi. Her iki iç split dış test dağılımını temsil etmiyor ve model hızlı biçimde train'e aşırı uyuyor.

Tanı koşusu her epoch'ta ek fit/validation/test tahmini yaptığı için 1.56 saniye sürdü; özgün eğitim döngüsü bu ek skorları hesaplamaz ve 0.74 saniyelik ölçüm kodun `perf_counter` sınırlarıyla tutarlıdır ([train_explainable_dl.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/train_explainable_dl.py:229)). Süre bir loglama hatası değildir. Seçili epoch tablosu [deepsurv_selected_epochs.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/deepsurv_selected_epochs.csv), tam eğri [deepsurv_diagnostic.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase70b_collapse_diagnostic/deepsurv_diagnostic.json) içindedir.

## Kapsam ve dürüstlük kontrolleri

- Hiçbir F/S/U/G/G_signed/R, N=50 bootstrap, null kalibrasyon, BH/BY veya METABRIC işlemi çalıştırılmadı.
- `components.py` ve `aggregation.py` değiştirilmedi; SHA256 değerleri sırasıyla `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa` ve `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` olarak kaydedildi.
- Rapor ve özet çıktıları gerçek hasta/örnek kimliği içermez.
- Tanı yeni bir performans-optimizasyon kuralı önermiyor ve fold'u sonucuna bakarak yeniden seçmiyor. Bulgular, Faz 7.1 gibi pahalı bir koşudan önce mevcut tek-fold baseline tasarımının yeterince kararlı olmadığını gösteriyor.
