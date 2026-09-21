# Faz 7.2 ADIM 1 — CoxPH-v2 N=50 audit

Yeni 153 hastalık audit havuzu ve 27 özellikli ortak ailede N=50 patient-cluster bootstrap tamamlandı. **9/27** özellik nokta tahmininde WI≥0.8 geçti. Darboğaz dağılımı `{'R': 24, 'G_claim': 3}`.

CoxPH candidate ve reference aynı kapalı-form attribution ile hesaplandığından F sayısal öz-tutarlılık ölçümüdür. 28 replika önceden belirlenen remediation/retry yolunu gerektirdi; nihai warning taşıyan replika sayısı sıfırdı. N=50 koşusu **2154.232 saniye** sürdü.

## ADIM 5 — Tam-aile sertifikasyon

Taze null kapısı geçtikten sonra m=27 ailenin tamamında B=4999 `recursive_studentized` uygulandı. Nokta/BH/BY sayıları **9/8/8**.

| feature | WI | p_joint | q_BH | q_BY |
| --- | ---: | ---: | ---: | ---: |
| rna_ENSG00000244362.3 | 0.8299 | 0.0002 | 0.0009 | 0.0035 |
| rna_ENSG00000169800.14 | 0.8292 | 0.0002 | 0.0009 | 0.0035 |
| rna_ENSG00000133742.14 | 0.8248 | 0.0038 | 0.0128 | 0.0499 |
| rna_ENSG00000204246.3 | 0.8205 | 0.0006 | 0.0023 | 0.0090 |
| rna_ENSG00000186965.5 | 0.8291 | 0.0002 | 0.0009 | 0.0035 |
| rna_ENSG00000254834.4 | 0.8297 | 0.0002 | 0.0009 | 0.0035 |
| clin_numeric__age_at_index.demographic | 0.9724 | 0.0002 | 0.0009 | 0.0035 |
| clin_categorical__ajcc_pathologic_stage.diagnoses_Stage III | 0.8509 | 0.0002 | 0.0009 | 0.0035 |

Formül SHA256 değerleri doğrulandı. Çıktılarda hasta/örnek kimliği yoktur.
