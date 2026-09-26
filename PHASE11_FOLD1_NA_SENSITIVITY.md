# Faz 11 — Fold-1 CoxPH-v2 NA yönsel duyarlılık kontrolü

Bu çalışma **diagnostic/exploratory** bir duyarlılık analizidir. Faz 7.5 üretimi replika 15'te durduğu için diskte başlangıçta 49 değil 15 başarılı checkpoint vardı; aynı deterministik hasta-kümesi örneklemesi ve dondurulmuş `penalizer=0.1`, separation dışlama ve `max_steps=200/step_size=0.25` protokolüyle eksik 16–49 replikaları yalnız bu tanısal çıktı için tamamlandı.

## Sonuçlar

| Senaryo | Point-adequate | BH | BY |
| --- | ---: | ---: | ---: |
| A_exclude_rep15_N49 | 19 | 13 | 11 |
| B_rescue_rep15_penalizer_0.15_N50 | 20 | 13 | 11 |

### ADIM A — replika 15 dışarıda, N=49

BH sertifikalı özellikler (13/80): `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I`, `clin_numeric__age_at_index.demographic`, `rna_ENSG00000107187.17`, `rna_ENSG00000109158.11`, `rna_ENSG00000157343.9`, `rna_ENSG00000165762.3`, `rna_ENSG00000177182.11`, `rna_ENSG00000186965.5`, `rna_ENSG00000196406.4`, `rna_ENSG00000197038.6`, `rna_ENSG00000204246.3`, `rna_ENSG00000243040.6`, `rna_ENSG00000270547.6`.

### ADIM B — yalnız replika 15 penalizer=0.15 ile kurtarıldı, N=50

Replika 15 `penalizer=0.15`, `max_steps=200`, `step_size=0.25` ile uyarısız yakınsadı. BH sertifikalı özellikler (13/80): `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I`, `clin_numeric__age_at_index.demographic`, `rna_ENSG00000107187.17`, `rna_ENSG00000109158.11`, `rna_ENSG00000157343.9`, `rna_ENSG00000165762.3`, `rna_ENSG00000177182.11`, `rna_ENSG00000186965.5`, `rna_ENSG00000196406.4`, `rna_ENSG00000197038.6`, `rna_ENSG00000204246.3`, `rna_ENSG00000243040.6`, `rna_ENSG00000270547.6`.

ADIM A'dan B'ye BH sertifikasyonuna eklenenler: yok. Düşenler: yok. Point-adequacy kümesine eklenenler: `clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II`; düşenler: yok. Kurtarılan katsayı vektörünün L2 normu **1.086970**; diğer 49 replikada ortalama **1.431849**, SS **0.107485**, aralık **[1.235166, 1.730571]** ve z-skoru **-3.209**. Kurtarılan vektör diğer 49 replikaya göre belirgin düşük normludur; bu nedenle ADIM B yalnız yönsel tanı olarak yorumlanmalıdır.

Formül SHA256 değerleri değişmedi: `components.py=7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`, `aggregation.py=10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68`. Yeni hasta/model kimliği dışa aktarılmadı. İlk tam üretim wall-clock **1606.582 sn**.

Bu analiz dondurulmuş retry protokolünü DEĞİŞTİRMEMEKTEDİR; Table 4'teki resmi NA hücresi ve ana metindeki sertifikasyon kararları aynen korunmaktadır.

49-replikalı dışlama analizinden gevşetilmiş-ceza N=50 analizine geçiş, BH sertifikasyon sayısının **aynı kaldığını** göstermektedir.
