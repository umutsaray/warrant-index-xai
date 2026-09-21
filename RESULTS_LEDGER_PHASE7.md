# Faz 7 sonuç defteri

Bu defter üç farklı performans sayısını birbirinden ayırır. Eski kirli split skorları geçerli genelleme kanıtı değildir. Fold-0 skorları Faz 7.2 audit'inin kullanacağı sabit model örneklerini tanımlar. Beş-fold ortalamaları ise hasta-gruplu genelleme performansının bu üçü arasındaki en güvenilir tahminidir.

| Model | Eski, kirli split C-index | Audit referansı: temiz fold 0 | 5-fold CV genelleme tahmini |
|---|---:|---:|---:|
| CoxPH | 0.706477 | **0.667704** (v2 retry) | **0.672006 [0.650173, 0.705882]** |
| RSF | 0.714795 | **0.538521** | **0.663116 [0.538521, 0.789609]** |
| DeepSurv | 0.745692 | **0.496498** | **0.612080 [0.496498, 0.719193]** |

## Kaynak ayrımı

- **Eski/kirli:** Kardeş projedeki orijinal 673/169 sample-level split 23 hastayı train ve test arasında paylaştırıyordu. Bu skorlar yalnız tarihsel karşılaştırmadır.
- **Audit referansı:** [PHASE7_REFERENCE_MODELS.md](C:/Users/UMUT/Desktop/ExCred/warrantxai/PHASE7_REFERENCE_MODELS.md) içindeki önceden seçilmiş fold-0 modelleridir. Bunlar Faz 7.2'de açıklanacak somut modellerdir; tek-fold skorları genel performans özeti değildir.
- **CV genelleme:** [PHASE7_CV_SUMMARY.md](C:/Users/UMUT/Desktop/ExCred/warrantxai/PHASE7_CV_SUMMARY.md) ve [cv_summary.json](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase71_cv/cv_summary.json) içindeki beş hasta-gruplu fold'un ortalama ve min–maks aralığıdır.

## Dondurulmuş artefakt yönlendirmesi

| Model | Faz 7.2'nin okuyacağı dosya | SHA256 |
|---|---|---|
| CoxPH | `results/phase7_retrain/coxph_baseline_v2.joblib` | `7094e959cc7bd4bf5fa5103db512f8d6071af88ad63f7f78fa3facc2a945915d` |
| RSF | `results/phase7_retrain/rsf_baseline.joblib` | `216ae54ad829d4e3144543a41979f5464ee8851ce4db5282e0500e7f15235e1b` |
| DeepSurv | `results/phase7_retrain/explainable_dl_model.pt` | `435a772042389310e3211b5bbda6ea936acf9c2f0450bb13f1fd9016830da756` |

RSF ve DeepSurv fold-0 artefaktları yeniden eğitilmedi. CoxPH v1 korunmuştur. Bu defter yalnız baseline ve CV sonuçlarını kaydeder; XAI audit sonucu içermez.
