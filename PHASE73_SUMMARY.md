# Faz 7.3 — Özet

Faz 7.2'nin saklı per-feature sonuçları yeniden formatlandı; bootstrap, null-kalibrasyon, sertifikasyon, model veya explainer çağrısı yapılmadı.

| Model | Family | Point adequate | BH | BY |
| --- | --- | --- | --- | --- |
| CoxPH-v2 | 27 | 9 | 8 | 8 |
| RSF | 27 | 17 | 12 | 12 |
| DeepSurv | 27 | 0 | 0 | 0 |

- Yeni tam tabloların her biri 27 satırdır ve Faz 7.2 sayıları assertion ile birebir doğrulanmıştır.
- CoxPH/RSF BH kesişimi **6**, CoxPH-only **2**, RSF-only **6** özelliktir. BH ve BY kümeleri aynıdır.
- Aile 22 RNA, 1 yaş ve 4 AJCC stage indicator özelliğinden oluşur.
- CoxPH–RSF BH kesişiminde aynı 12 query üzerindeki attribution yönü **0/6** özellikte aynıdır.
- Formül SHA256 değerleri değişmemiştir; çıktı tablolarında hasta veya örnek kimliği yoktur.

## Süre

| Step | Seconds |
| --- | --- |
| ADIM 1 — full tables | 0.0139 |
| ADIM 2 — intersection | 0.0021 |
| ADIM 3 — feature types | 0.0021 |
| ADIM 4 — directionality | 0.0887 |

Toplam wall-clock **0.1069 saniye**.
