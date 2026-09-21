# Faz 7.0 — Hasta-gruplu split yeniden tasarımı

## Tasarım kararı

Split birimi TCGA `sample_key` yerine ilk 12 karakterden türetilen hasta grubudur. Uygulama `StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=42)` kullanır: [impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:250). Beş test fold'u arasından önce hedef test satır oranı 0.20'ye mutlak uzaklığı, eşitlikte genel event oranına uzaklığı, sonra fold numarasını küçülten fold deterministik olarak seçilir ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:275)).

Bu yaklaşım seçildi çünkü aynı hastanın bütün örneklerini tek fold'da tutarken orijinal yaklaşık 80/20 satır oranını ve event/censoring dengesini birlikte korur. Kullanılan scikit-learn 1.9.0 sürümü `StratifiedGroupKFold` desteğine sahiptir.

Beş adayın agregatları [patient_grouped_fold_candidates.csv](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/patient_grouped_fold_candidates.csv) içindedir. Fold 0–2, 674/168 satır ve %15.7270/%15.4762 train/test event oranı; fold 3–4, 673/169 satır ve %15.6018/%15.9763 oranı verdi. 0.20 hedef satır oranına uzaklık ilk grupta 0.000475, ikinci grupta 0.000713 olduğundan fold 0 seçildi. Fold 0–2 seçim ölçütlerinde eşit olduğundan son bağlayıcı fold numarasıdır.

## Eski ve yeni split karşılaştırması

| Ölçü | Eski sample-level split | Yeni patient-grouped split |
|---|---:|---:|
| Train satırı | 673 | 674 |
| Test satırı | 169 | 168 |
| Train benzersiz hasta | 628 | 616 |
| Test benzersiz hasta | 164 | 153 |
| Train/test ortak hasta | **23** | **0** |
| Train event oranı | %15.7504 | %15.7270 |
| Test event oranı | %15.3846 | %15.4762 |

Yeni sayılar [split_summary.md](C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/split_summary.md:1) tarafından kaydedildi. Kod hem tam örnek anahtarı kesişimini hem ilk-12-karakter hasta kesişimini çalışma zamanında sıfır olmaya zorlar ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:288)).

## Preprocessing ve feature-selection sınırı

Methylation varyans filtresi yeni train fold'unda öğrenildi; KNN imputer train'de fit edilip train/test'e uygulandı; clinical median/mode imputation ve one-hot encoder da yalnız train'de fit edildi. Yeni ham işlenmiş matrisler 674×72564 ve 168×72564'tür. Yeni train'de görülmeyip test'te bulunan tek bir ırk kategorisinin klinik giriş boyutunu 16'dan 15'e indirmemesi için eski model şemasındaki kategori sözlüğü kodda önceden sabitlendi ([impute_and_split.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/impute_and_split.py:41)); encoder yine yalnız train'de fit edildi. Böylece test kategorileri şema keşfi için okunmadan 200+200+200+16=616 giriş korundu.

Feature-selection kodu yeni `X_train.parquet` ile `y_train.parquet` dosyalarını okur ([select_features.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/select_features.py:173)); her omics bloğundaki univariate Cox sıralaması yalnız `x_train` ve `y_train` ile hesaplanır ([select_features.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/select_features.py:183)). Test matrisi yalnız train'de seçilmiş `final_columns` üzerine projekte edilir ve selection istatistiğine girmez ([select_features.py](C:/Users/UMUT/Desktop/ExCred/Multi-Omics%20Data%20Fusion%20with%20Explainable%20Deep%20Learning%20for%20Survival%20Prediction%20in%20Breast%20Cancer/select_features.py:207)). `selected_features.json` ayrıca `selection_scope=train_only` kaydını üretir.

## İzolasyon

Yeni çıktılar `C:/Users/UMUT/Desktop/ExCred/warrantxai/results/phase7_retrain/` altında yazılır. Eski 673/169 parquet dosyaları ve eski üç model artefaktı değiştirilmemiştir. Faz 7.1 audit işlemleri bu adımda çalıştırılmayacaktır.
