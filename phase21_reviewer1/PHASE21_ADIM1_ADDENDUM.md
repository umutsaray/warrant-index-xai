# Faz 21 ADIM 1 eki

## BH q-değerleri

method,q_lt_0_05,q_lt_0_10,q_min,q_median,q_max
filtered_LR,113,200,0.003799918712851941,0.0471715101199035,0.061081693862346714
filtered_score,200,200,1.1975339926901468e-17,0.005659959292406913,0.01561499267075438

Q-değerleri 43,877 özelliklik filtrelenmiş ailenin tamamı üzerinde BH ile hesaplandı; özet ilk 200 içindir.

## Başarısız 87 Wald özelliği

failed_n,zero_rate_min,zero_rate_median,zero_rate_max,nonzero_n_min,nonzero_n_median,nonzero_n_max,nonzero_events_min,nonzero_events_median,nonzero_events_max,complete_separation_n,ridge_sign_agreement_n,ridge_sign_agreement_rate
87,0.8427299703264095,0.9925816023738872,0.9970326409495548,2,5.0,106,1,2.0,20,3,31,0.3563218390804598

Başarısızlıkların yalnız 3/87'si basit tam ayrışma göstergesine sahiptir; baskın desen çok yüksek sıfır oranı ve az nonzero/event desteğidir. Ridge işareti, orijinal tarama katsayısıyla 31/87 özellikte aynıdır.

## Gerçek METABRIC ekspresyon doğrulaması

- Study/profile/sample list: `brca_metabric` / `brca_metabric_mrna` / `brca_metabric_mrna`.
- İndirme kaydı (UTC): `2026-11-08T22:54:04.898724+00:00`.
- Klinik-uygun ekspresyon kohortu: 1,979 hasta.
- Ham hasta-düzeyi yanıtlar yalnız git-ignore'lı `data\external\metabric_phase21` altında; commit edilmez.
- Batch SHA-256 ve boyutlar yerel `download_provenance.json` içinde kayıtlıdır.

Gerçek ekspresyon ölçütü: gen sütunu var, 1,979 hastada eksik yok ve SD>0. Filtre+LR gerçek ortak RNA **130**, filtre+skor gerçek ortak RNA **100**. Gen düzeyi ayrıntılar `PHASE21_ADIM1_METABRIC_EXPRESSION_VALIDATION.csv` içindedir.
