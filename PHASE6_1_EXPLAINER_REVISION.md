# Faz 6.1 — DeepSurv explainer/bütçe tasarım revizyonu

> **Koşulluluk uyarısı:** DeepSurv Phase 6.1 conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds; N=2 is a design result, not certification.

## Önceden tanımlanan karar ölçütü

Additivity için makul sınır, her explainer çağrısında maksimum mutlak reconstruction hatasının aynı çağrıdaki 12 query ve expected/base değerinin risk-skoru aralığına oranının en çok **%5** olmasıdır. Bu göreli ölçüt DeepSurv risk skorunun keyfi ölçeğine bağlanmaz. F ayırt ediciliği frozen `a=0.5` ile ayrıca kontrollü bozulma eğrisinde sınandı.

## Seçenek A — GradientExplainer bütçe eğrisi

| option | design | candidate_budget | reference_budget | max_abs_additivity_error_candidate | max_relative_additivity_error_candidate | max_abs_additivity_error_reference | max_relative_additivity_error_reference | selected_mean_F | selected_min_F | all77_mean_F | all77_min_F | elapsed_seconds | additivity_pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| A_GradientExplainer | 8/128 | 8 | 128 | 0.2635 | 0.3287 | 0.0700 | 0.0411 | 0.9999 | 0.9998 | 0.9998 | 0.9921 | 0.7531 | False |
| A_GradientExplainer | 32/256 | 32 | 256 | 0.1738 | 0.0952 | 0.0941 | 0.0264 | 1.0000 | 1.0000 | 0.9999 | 0.9974 | 1.4587 | False |
| A_GradientExplainer | 64/512 | 64 | 512 | 0.1303 | 0.0838 | 0.0264 | 0.0302 | 1.0000 | 1.0000 | 1.0000 | 0.9987 | 2.7211 | False |
| A_GradientExplainer | 128/1024 | 128 | 1024 | 0.1552 | 0.0683 | 0.0409 | 0.0224 | 1.0000 | 0.9998 | 0.9999 | 0.9985 | 5.2725 | False |
| A_GradientExplainer | 256/2048 | 256 | 2048 | 0.0484 | 0.0265 | 0.0185 | 0.0161 | 1.0000 | 1.0000 | 1.0000 | 0.9999 | 10.2748 | True |
| A_GradientExplainer | 512/4096 | 512 | 4096 | 0.0385 | 0.0295 | 0.0211 | 0.0115 | 1.0000 | 1.0000 | 1.0000 | 0.9998 | 20.3791 | True |

En yüksek denenen 512/4096 tasarımında candidate maksimum mutlak/göreli additivity hatası **0.038499 / 2.95%**, reference için **0.021063 / 1.15%** oldu. Bütçe artışı ölçüldü; karar yalnız F'nin yüksek görünmesine dayandırılmadı.

## Seçenek B — DeepExplainer background bütçesi

| option | design | candidate_budget | reference_budget | max_abs_additivity_error_candidate | max_relative_additivity_error_candidate | max_abs_additivity_error_reference | max_relative_additivity_error_reference | selected_mean_F | selected_min_F | all77_mean_F | all77_min_F | elapsed_seconds | additivity_pass |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| B_DeepExplainer | 16/64 | 16 | 64 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.9958 | 0.9787 | 0.9981 | 0.9787 | 1.6593 | True |
| B_DeepExplainer | 16/128 | 16 | 128 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.9953 | 0.9748 | 0.9983 | 0.9748 | 7.4282 | True |

DeepExplainer, 77-feature integrated alt-oyununda hem k=16 candidate hem k=64/128 reference ile teknik olarak çalıştı. Candidate ve reference aynı TCGA train dağılımından deterministik, nested background örnekleridir; fark yalnız background integrasyon bütçesidir.

## Frozen F kontrollü bozulma kontrolü

| design | added_candidate_error | error_as_fraction_of_a | selected_mean_F | selected_min_F | selected_below_tau |
| --- | --- | --- | --- | --- | --- |
| A_GradientExplainer_256/2048 | 0.0000 | 0.0000 | 1.0000 | 1.0000 | 0 |
| A_GradientExplainer_256/2048 | 0.1000 | 0.2000 | 0.9600 | 0.9599 | 0 |
| A_GradientExplainer_256/2048 | 0.2500 | 0.5000 | 0.7499 | 0.7496 | 8 |
| A_GradientExplainer_256/2048 | 0.5000 | 1.0000 | 0.0003 | 0.0001 | 8 |
| B_DeepExplainer_16/128 | 0.0000 | 0.0000 | 0.9953 | 0.9748 | 0 |
| B_DeepExplainer_16/128 | 0.1000 | 0.2000 | 0.9585 | 0.9311 | 0 |
| B_DeepExplainer_16/128 | 0.2500 | 0.5000 | 0.7524 | 0.6885 | 7 |
| B_DeepExplainer_16/128 | 0.5000 | 1.0000 | 0.0374 | 0.0053 | 8 |

Bozulma yalnız sekiz sağlık-kontrol feature'ının candidate attributionlarına eklendi; reference ve `a=0.5` değiştirilmedi. Her seçenek için ortalama F monoton düştü; +0.25 düzeyinde temiz sonuca göre 0.20'den fazla azalma kuralı GradientExplainer için **GEÇTİ**, DeepExplainer için **GEÇTİ**. Bu deney aynı zamanda F'nin mutlak tolerans davranışını gösterir; additivity kalite kapısı F'den ayrı tutulmalıdır.

## Öneri ve yeni N=2 pilot

Öneri: **DeepExplainer candidate k=16 / reference k=128**. Gerekçe: near-exact additivity, background büyüklüğüyle gerçek ve teknik olarak çalışan bütçe ekseni, kontrollü bozulmaya yanıt veren frozen F ve düşük runtime. GradientExplainer bütçeleri karşılaştırma amacıyla korunmuştur.

| inference_scope | block | features | assessed | mean_F | mean_S | mean_U | mean_G_claim | mean_R | mean_WI |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DeepSurv Phase 6.1 conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds; N=2 is a design result, not certification. | clin | 5 | 5 | 0.9992 | 0.9940 | 0.9917 | 0.9769 | 0.5961 | 0.5961 |
| DeepSurv Phase 6.1 conditions on fixed 12+12 queries and deterministic nested DeepExplainer k=16/128 TCGA backgrounds; N=2 is a design result, not certification. | rna | 72 | 72 | 0.9982 | 0.9724 | 0.9826 | 0.9401 | 0.6761 | 0.6724 |

77/77 feature değerlendirildi; **11** feature WI≥0.8 geçti. Bottleneck dağılımı `{'R': 75, 'G_claim': 2}`. Bu N=2 tasarım sonucudur; N≥20, null kalibrasyon ve FDR başlatılmadı.

Gerçek toplam wall-clock **53.389 sn**. Formül SHA256 değişmedi: **True**. Kimlik-sızıntısı kontrolü geçti.
