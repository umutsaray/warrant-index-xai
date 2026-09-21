# Faz 7.1 — XAI audit referans modelleri

## Sabit referans seti

Faz 7.2 XAI audit'i yalnız önceden seçilmiş hasta-gruplu **fold 0** modellerini kullanmalıdır. Fold 0 performansa bakılarak seçilmedi; Faz 7.0'da satır/event dengesi ölçütlerinin deterministik fold-numarası bağlayıcısıyla belirlenmişti. Beş-fold CV sonuçları bu referans seçimini değiştirmek için kullanılmadı.

| Model | Resmi model dosyası | Fold-0 test C-index | Durum |
|---|---|---:|---|
| CoxPH | [coxph_baseline_v2.joblib](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/coxph_baseline_v2.joblib) | **0.667704** | `penalizer=0.1`; `max_steps=200`, `step_size=0.25`; uyarısız retry |
| RSF | [rsf_baseline.joblib](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/rsf_baseline.joblib) | **0.538521** | Faz 7.0 artefaktı değişmeden korundu |
| DeepSurv | [explainable_dl_model.pt](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/explainable_dl_model.pt) | **0.496498** | Faz 7.0 artefaktı değişmeden korundu; checkpoint epoch 38 |

CoxPH v2, aynı 610 feature, aynı StandardScaler ve aynı `penalizer=0.1` ile yalnız doğrulanmış optimizer retry kullanılarak üretildi. Yakınsama uyarısı çıkmadı ve beklenen C-index 0.667704280155642 olarak birebir doğrulandı. Üretim ve koruma kontrolleri [35_phase71_reference_models.py](C:/Users/UMUT/Desktop/ExCred/warrantxai/experiments/35_phase71_reference_models.py) tarafından yapılır; makine-okunur kayıt [reference_models.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/reference_models.json) içindedir.

Eski fold-0 CoxPH v1 silinmedi: [coxph_baseline.joblib](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/coxph_baseline.joblib). V1, RSF ve DeepSurv dosyalarının SHA256 değerleri CoxPH-v2 üretiminden ve CV koşusundan önce/sonra aynıdır. Faz 7.2 model yönlendirmesi CoxPH için açıkça `_v2`, diğer iki model için yukarıdaki mevcut dosyaları kullanmalıdır.

## Kapsam

Bu adım model referanslarını sabitledi. F/S/U/G/G_signed/R, bootstrap, null kalibrasyon, BH/BY ve METABRIC çalıştırılmadı.
