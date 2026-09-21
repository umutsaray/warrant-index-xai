# Faz 7.0 — Hasta-gruplu split üzerinde model yeniden eğitimi

## Sonuç

Üç model, hasta düzeyinde ayrık yeni 674/168 split üzerinde eski hiperparametreleri, seed'leri ve eğitim prosedürleri korunarak yeniden eğitildi. Yeni test kümesinin 153 hastasının hiçbiri train'deki 616 hastayla kesişmez. Yeni modeller [results/phase7_retrain](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/) altında saklandı; eski artefaktların SHA256 değerleri teşhis kaydıyla birebir eşleşmeye devam ediyor.

| Model | Eski test C-index | Yeni test C-index | Fark (yeni−eski) |
|---|---:|---:|---:|
| CoxPH | 0.706477 | 0.631128 | -0.075348 |
| RSF | 0.714795 | 0.538521 | -0.176274 |
| DeepSurv | 0.745692 | 0.496498 | -0.249194 |

Üçündeki düşüş de büyüktür; özellikle RSF ve DeepSurv yeni testte eski performansı korumadı. Bu sonuç, yeni hasta-gruplu split ve yalnız yeni train'de yeniden öğrenilen feature paneli altında ölçülen baseline performansıdır. Faz 7.0 kapsamında ek hiperparametre araması veya yöntem değişikliği yapılmadı.

## Split ve panel doğrulamaları

- Train/test satır sayısı 674/168, benzersiz hasta sayısı 616/153, örnek ve hasta kesişimi 0/0'dır.
- Event oranı train'de %15.7270, test'te %15.4762'dir.
- Panel yalnız yeni train'de seçildi ve 616 sütundur: 200 RNA, 200 miRNA, 200 methylation, 16 clinical. `selected_features.json` içindeki `selection_scope` değeri `train_only`'dir.
- Klinik kategori sözlüğü önceden sabitlenerek eski 16 girişli klinik mimari korundu; test verisi kategori keşfi veya feature sıralaması için kullanılmadı.

## Model ayrıntıları

**CoxPH.** `lifelines.CoxPHFitter(penalizer=0.1)`, `max_steps=50` ve aynı separation ön kontrolü kullanıldı. Altı sütun dışlandı ve 610 sütunla fit edildi. Newton–Raphson'ın yeterli yakınsamadığı uyarısı üretildi; uyarı sessizce bastırılmadı ve [baseline_results.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/baseline_results.json) içine kaydedildi. Eğitim süresi 24.44 saniyeydi.

**RSF.** 500 ağaç, `min_samples_split=10`, `min_samples_leaf=5`, `max_features="sqrt"`, `random_state=42` kullanıldı. Fit sırasında hata veya uyarı çıkmadı; süre 2.56 saniyeydi.

**DeepSurv.** Blok boyutları RNA 200→128→32, miRNA 200→64→16, methylation 200→128→32 ve clinical 16→32→8; birleşik katmanlar 88→64→16→1 olarak korundu. Dropout 0.3, AdamW `lr=0.001`, `weight_decay=0.0001`, en çok 300 epoch, patience 30 ve seed 42 değişmedi. Eski checkpoint epoch 26 iken yeni en iyi checkpoint epoch 38'de seçildi; validation C-index 0.852814 oldu ve eğitim epoch 68'de erken durdu. Epoch seçiminin belirgin değişmesi ve test C-index'in 0.496498'e düşmesi olağanmış gibi gizlenmemelidir. Eğitim süresi 0.74 saniyeydi.

DeepSurv'ün iç fit/validation ayırma prosedürü talimat gereği değiştirilmedi; hâlâ satır düzeyinde stratified split kullanıyor ve yeni train içinde 15 hasta iki iç parçaya yayılıyor. Bu durum dış test kümesine hasta sızıntısı oluşturmaz; ancak ileride checkpoint seçimi protokolü ayrıca yeniden tasarlanacaksa ele alınması gereken bir sınırlamadır.

## Artefaktlar ve kapsam

Yeni model yolları:

- [coxph_baseline.joblib](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/coxph_baseline.joblib)
- [rsf_baseline.joblib](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/rsf_baseline.joblib)
- [explainable_dl_model.pt](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/explainable_dl_model.pt)
- [validation.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/validation.json)

Eski yedi veri/model artefaktının hash kontrolü 7/7 geçti. Faz 7.0'da F/S/U/G/G_signed/R, N=50 bootstrap, null kalibrasyon, BH/BY, METABRIC veya başka bir audit adımı çalıştırılmadı.
