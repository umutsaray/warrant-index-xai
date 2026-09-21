# Faz 7.6 — Fold-0 RSF split-kullanım diagnostik

Ana fold-0 RSF artefaktındaki 500 fitted ağacın `tree_.feature` dizileri doğrudan okundu; model yeniden fit edilmedi. `feature_importances_` bu scikit-survival sınıfında mevcut görünse de çağrı doğrudan `NotImplementedError` üretir; bu yüzden değer uydurulmadı.

| Feature        | Feature ID                                                  | Label       | feature_importances_                                                         | Trees with ≥1 split | Total split nodes |
| --- | --- | --- | --- | --- | --- |
| KRTAP19-7      | rna_ENSG00000244362.3                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| OR6C1          | rna_ENSG00000205330.4                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 10                  | 10                |
| OR4C46         | rna_ENSG00000185926.1                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| RBMY1F         | rna_ENSG00000169800.14                                      | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| OR13C3         | rna_ENSG00000204246.3                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| KRTAP19-2      | rna_ENSG00000186965.5                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| ACCSL          | rna_ENSG00000205126.2                                       | informative | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 20                  | 21                |
| OR5M10         | rna_ENSG00000254834.4                                       | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 0                   | 0                 |
| AJCC Stage I   | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage I   | informative | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 2                   | 2                 |
| AJCC Stage II  | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage II  | informative | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 10                  | 10                |
| AJCC Stage III | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | informative | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 19                  | 19                |
| AJCC Stage IV  | clin_categorical__ajcc_pathologic_stage.diagnoses_Stage IV  | degenerate  | unavailable: scikit-survival RandomSurvivalForest raises NotImplementedError | 25                  | 25                |

Beklenen tam uyum görülmedi. Sekiz degenerate feature'ın altısı hiçbir ağaçta kullanılmadı; **OR6C1 10/500 ağaçta 10 split**, **AJCC Stage IV 25/500 ağaçta 25 split** içinde kullanıldığı halde saklı fixed-query candidate attribution'ları exact sıfırdır. Bu tablo yalnız gözlenen tutarsızlığı kaydeder. Süre **1.097 sn**.
