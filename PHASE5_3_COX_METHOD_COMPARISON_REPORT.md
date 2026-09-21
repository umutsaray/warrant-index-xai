# Faz 5.3B — CoxPH yöntem tutarlılık kontrolü

## Sonuç

Faz 3.2'nin önceden tanımlı **77 feature** ailesinin tamamı B=4999 ile recursive_studentized yönteminden geçirildi. Nokta tahmininde WI≥0.8 sayısı 22; recursive_studentized BH/BY sayıları sırasıyla **8/8** oldu. Faz 5.2 basic-root BH sonucu **8** idi.

| criterion | count |
| --- | --- |
| assessed_full_family_m | 77 |
| point_WI_ge_0_8 | 22 |
| basic_root_BH | 8 |
| recursive_studentized_BH | 8 |
| recursive_studentized_BY | 8 |
| BH_intersection | 8 |
| basic_only | 0 |
| recursive_only | 0 |

Her iki yöntemle de BH sertifikalı kalan aynı sekiz feature:

| feature | block | WI | p_joint | q_BH | q_BY |
| --- | --- | --- | --- | --- | --- |
| rna_ENSG00000133742.14 | rna | 0.8412 | 0.0002 | 0.0019 | 0.0095 |
| rna_ENSG00000169800.14 | rna | 0.8331 | 0.0002 | 0.0019 | 0.0095 |
| rna_ENSG00000133800.9 | rna | 0.8367 | 0.0002 | 0.0019 | 0.0095 |
| rna_ENSG00000109158.11 | rna | 0.8335 | 0.0002 | 0.0019 | 0.0095 |
| rna_ENSG00000187758.8 | rna | 0.8371 | 0.0002 | 0.0019 | 0.0095 |
| rna_ENSG00000204246.3 | rna | 0.8371 | 0.0002 | 0.0019 | 0.0095 |
| clin_numeric__age_at_index.demographic | clin | 0.9503 | 0.0002 | 0.0019 | 0.0095 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I | clin | 0.8670 | 0.0002 | 0.0019 | 0.0095 |

Küme değişiklikleri:

_BH sertifikalı feature kümesi değişmedi._

## Bileşen p-değeri ayrışması

| component | mean_abs_p_difference | max_abs_p_difference | recursive_p_larger_features | recursive_p_smaller_features |
| --- | --- | --- | --- | --- |
| F | 0.0000 | 0.0000 | 0 | 0 |
| S | 0.0000 | 0.0000 | 0 | 0 |
| U | 0.0000 | 0.0000 | 0 | 0 |
| G_claim | 0.0014 | 0.0676 | 3 | 0 |
| R | 0.0000 | 0.0000 | 0 | 0 |

P-değeri değişen feature-bileşen çiftleri:

| feature | component | p_basic_root | p_recursive_studentized | absolute_difference |
| --- | --- | --- | --- | --- |
| rna_ENSG00000121351.8 | G_claim | 0.0436 | 0.1112 | 0.0676 |
| rna_ENSG00000090402.8 | G_claim | 0.0002 | 0.0384 | 0.0382 |
| rna_ENSG00000178522.14 | G_claim | 0.0002 | 0.0022 | 0.0020 |

F ve U tek TCGA audit yaprağı üzerinde studentized; S ve R aynı yaprak üzerinde basic test edildi. `G_claim=min(G,G_signed)` için magnitude ve signed yaprak p-değerleri ayrı hesaplanıp `max` ile recursive IUT yapıldı. Sonra beş bileşen yine `max` ile feature düzeyi IUT'ye birleştirildi. Bu bir inferans katmanı değişikliğidir; frozen skor formülleri değişmedi.

İki yöntem aynı BH kümesini verdi. Makalede, hem CoxPH N=50 hem RSF N=20 null kapısından geçmiş tek ortak yöntem olduğu için recursive_studentized birincil; basic-root CoxPH sonucu sağlamlık kontrolü olarak önerilir.

## Dürüstlük kontrolleri

Çıkarım sabit 141 TCGA ve 1466 METABRIC sorgusuna, preprocessing'e ve exact Cox attribution tanımına koşulludur. Formül SHA256 değerleri Faz 5 pilotuyla aynıdır: **True**. Çıktılarda hasta/örnek kimliği yoktur; model fit veya explainer çağrısı yapılmadı. Gerçek wall-clock **11.673 saniye**.
