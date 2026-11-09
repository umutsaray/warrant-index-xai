# Faz 21 ADIM 7 — sayısal özet

## 7.1 RSF-v3 C-index

- İç-eğitim: 0.558292; ters işaret: 0.441708.
- İç-validation: 0.670588; ters işaret: 0.329412.
- Dış-eğitim: 0.571492; ters işaret: 0.428508.
- Dış-test: 0.384436; ters işaret: 0.615564.
- Alternatif RSF seed 43/44/45 dış-test C-index: 0.551751 / 0.539300 / 0.536965.
- Yön/kod uyuşmazlığı: HAYIR.

Hasta-bootstrap %95 GA:

| Model | C-index | Alt | Üst |
|---|---:|---:|---:|
| CoxPH-v3 | 0.515953 | 0.374395 | 0.681016 |
| RSF-v3 | 0.384436 | 0.256022 | 0.525717 |
| DeepSurv-v3 | 0.697276 | 0.553722 | 0.829938 |
| CoxPH-v1 | 0.667704 | 0.541912 | 0.781333 |
| RSF-v1 | 0.538521 | 0.400283 | 0.677048 |
| DeepSurv-v1 | 0.496498 | 0.385362 | 0.620431 |

## 7.2 Tamamlayıcı sayılar

- RSF-v3 BH sertifikalı: 14; unexercised (`Q_j=0`): 1.
- Bağlayıcılar: CoxPH `G_claim=2, R=63`; RSF `G_claim=1, R=64`; DeepSurv `R=65`.
- Cox 7/4/4: 12+12 sorgu, 65 özellikli R evreni, 65 özellikli aile; eşleştirilmiş üç-model tasarımıyla aynı.
- 2×2 sayıları `PHASE21_ADIM7_07_FOURCELL.csv` dosyasındadır.
- Resmî/TCGA-ölçekli nokta yeterlilik: CoxPH 7/4; RSF 30/16; DeepSurv 0/0.
- Resmî/TCGA-ölçekli medyan `G_claim`: CoxPH 0.973033/0.742080; RSF 0.959457/0.836997; DeepSurv 0.990885/0.919884.
- DeepSurv en iyi validation epoch: 30; validation/test C-index: 0.742857/0.697276; fark 0.045581; epoch Spearman 0.387859.
- v1 27 özellikten v3 ailesinde bulunan: 5.
- v1 sertifikalı 20 model-özellik çiftinden özelliği v3 ailesinde bulunan: 6; v3'te de sertifikalı: 6.

Formül dosyaları değiştirilmedi.
