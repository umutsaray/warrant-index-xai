# Faz 21 v3 — Aşama 1

İç-eğitim/iç-validation/dış-test satırları: 577/97/168. Hasta kesişimleri fit–validation=0, dış-train–test=0.

Model paneli 615 özelliktir (RNA/miRNA/methylation/clinical: 200/200/200/15). RNA denetim ailesi model paneli içindeki beyaz-listeli 139 adaydan ilk 60; ortak klinik 5; toplam aile 65.

Test C-index: CoxPH-v3=0.515953, RSF-v3=0.384436, DeepSurv-v3=0.697276. DeepSurv en iyi epoch=30, validation C-index=0.742857.

Programatik kanıt: seçim, methylation varyans/KNN, clinical imputation/encoding ve scaler fit satır hash'leri iç-eğitim hash'ine eşittir; validation kullanılmadı. Formül hash'leri değişmedi. Wall-clock 62.4 sn.
