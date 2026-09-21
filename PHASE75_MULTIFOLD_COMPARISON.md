# Faz 7.5 — Çoklu-fold karşılaştırması

Fold 0 resmi, önceden belirlenmiş audit referansıdır. Fold 1 ek sağlamlık analizidir; feature aileleri fold'a özgü train-only seçim ve TCGA–METABRIC kesişimi nedeniyle farklı büyüklüktedir.

| Fold | Model | C-index | Panel m | Point WI≥0.8 | BH | BY | Bottleneck | Durum |
|---:|---|---:|---:|---:|---:|---:|---|---|
| 0 | CoxPH-v2 | 0.667704 | 27 | 9 | 8 | 8 | R:24; G_claim:3 | resmi referans |
| 0 | RSF | 0.538521 | 27 | 17 | 12 | 12 | R:24; G_claim:3 | resmi referans |
| 0 | DeepSurv | 0.496498 | 27 | 0 | 0 | 0 | R:27 | resmi referans |
| 1 | CoxPH-v2 | 0.705882 | 80 | NA | NA | NA | NA | bootstrap replika 15 frozen retry sonrası başarısız |
| 1 | RSF | 0.789609 | 80 | 10 | 7 | 7 | R:76; G_claim:4 | gate geçti: 0/30 Holm flag |
| 1 | DeepSurv | 0.719193 | 80 | 7 | 0 | 0 | R:80 | gate geçti: 0/30 Holm flag |

## Sayısal kontroller

- **RSF–CoxPH sertifikasyon karşılaştırması:** Fold 1 CoxPH N=50 tamamlanamadığı için bu fold'da sayısal karşılaştırma **NA**; “RSF daha çok sertifikalı feature üretir” önermesi fold 1 CoxPH karşılığı olmadan değerlendirilemez.
- **DeepSurv:** Dış-test C-index'i 0.496498'den 0.719193'e çıkmasına ve point-adequate sayısı 0/27'den 7/80'e yükselmesine rağmen BH/BY sayısı **0/0** olarak kaldı.
- **R darboğazı:** Fold 1'de RSF için 76/80 (%95), DeepSurv için 80/80 (%100) R darboğazıdır. Tamamlanan nonlinear modellerde fold 0'daki R baskınlığı sürdü; CoxPH fold-bağımsızlığı bu çalıştırmada değerlendirilemedi.

RSF ve DeepSurv sonuçları sabit deterministik 12+12 query setlerine koşulludur. RSF ayrıca saklanmamış permütasyon yolları ile 8/128 bütçesine; DeepSurv nested DeepExplainer k=16/k=128 background tasarımına koşulludur.
