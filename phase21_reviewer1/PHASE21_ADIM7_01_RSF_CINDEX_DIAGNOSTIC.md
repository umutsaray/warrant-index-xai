# Faz 21 ADIM 7.1 — RSF-v3 C-index doğrulaması

## Yön
```
 model               scope  rows  events  cindex_risk  cindex_negated_risk             risk_definition
RSF-v3     inner_train_577   577      91     0.558292             0.441708 model.predict; higher=worse
RSF-v3 inner_validation_97    97      15     0.670588             0.329412 model.predict; higher=worse
RSF-v3     outer_train_674   674     106     0.571492             0.428508 model.predict; higher=worse
RSF-v3      outer_test_168   168      26     0.384436             0.615564 model.predict; higher=worse
```

## Hasta-bootstrap %95 GA
```
      model  point_cindex  bootstrap_valid  ci95_lower  ci95_upper  bootstrap_seed
   CoxPH-v3      0.515953             1000    0.374395    0.681016        20262107
     RSF-v3      0.384436             1000    0.256022    0.525717        20262107
DeepSurv-v3      0.697276             1000    0.553722    0.829938        20262107
   CoxPH-v1      0.667704             1000    0.541912    0.781333        20262107
     RSF-v1      0.538521             1000    0.400283    0.677048        20262107
DeepSurv-v1      0.496498             1000    0.385362    0.620431        20262107
```

## Alternatif RSF seed'leri
```
 random_state  test_cindex_risk  test_cindex_negated_risk
           43          0.551751                  0.448249
           44          0.539300                  0.460700
           45          0.536965                  0.463035
```

Kod yolu: ortak skor fonksiyonu bu dosyada satır 12'de `concordance_index_censored(event,time,risk)` çağrısıdır. v3 ve v1 riskleri sırasıyla satır 33 ve 35'te üretilip satır 36'da aynı fonksiyona verilir. İlk v3 Cox/RSF yolu `experiments/101_phase21_v3_stage1.py:194-202`; DeepSurv yolu aynı dosyada `train_deep` içindeki satır 139'dur. RSF için her iki sürümde de `model.predict` riski kullanılır.

Yön/kod uyuşmazlığı: HAYIR.
