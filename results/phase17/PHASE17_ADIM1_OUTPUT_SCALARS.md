# Faz 17 ADIM 1 — Atıf oyunlarının skaler çıktıları

| Model | Atıf oyunundaki skaler | Kod yolu | Birim/ölçek notu |
|---|---|---|---|
| CoxPH-v2 | Log kısmi hazard, kapalı form `β_j(x_j-background_j)` | `experiments/39_phase72_cox.py`: `source_phi = coefficients * anchors`; model doğrulaması `predict_log_partial_hazard` | Boyutsuz log-risk/log-hazard oranı. Candidate ve reference aynı kapalı-form hesap olduğundan F sayısal öz-tutarlılıktır. |
| RSF | `RandomSurvivalForest.predict(states)` risk skoru | `experiments/07_rsf_pilot.py::background_permutation_attributions`: bütün koalisyon durumlarında `model.predict(states)`; permütasyon marjinal farkları | scikit-survival RSF risk skoru, yani tahmin edilen kümülatif hazardın gözlenen olay zamanları üzerindeki toplamı; tek bir zaman noktasındaki sağkalım olasılığı değildir. |
| DeepSurv | Ağın tek skaler Cox log-risk çıktısı | `experiments/41_phase72_deepsurv.py`; `FlatOmicsWrapper`/`IntegratedActiveWrapper` çıktısı DeepExplainer'a verilir. Eğitim kaybı `cox_negative_log_likelihood(risk,...)` | Boyutsuz log-risk. DeepExplainer candidate/reference farkı background örnek büyüklüğüdür (k=16/k=128). |

Bu üç çıktı aynı fiziksel/istatistiksel ölçekte değildir. CoxPH-v2 ve DeepSurv log-risk üretirken RSF'nin skoru kümülatif hazard toplamıdır; bu nedenle mutlak `a,b,u,η` toleranslarının üç modele ortak uygulanması model-çıktı ölçeğine duyarlıdır. v2'nin model-özgü tek `s_m` normalizasyonunun gerekçesi budur; bu belge v1 veya v2 formül sonucunu değiştirmez.
