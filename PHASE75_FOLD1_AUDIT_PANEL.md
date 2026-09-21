# Faz 7.5 fold 1 ADIM 0 — Audit paneli

| Ölçü | Eski | Yeni |
|---|---:|---:|
| Ortak feature | 77 | **80** |
| TCGA audit hastası | 141 | **153** |
| Tam METABRIC query hastası | 1466 | **1466** |

Yeni fold-1 test setindeki 168 satır 153 benzersiz hastaya aittir. Train ile hasta kesişimi sıfırdır; eski 23-hasta çıkarma adımı uygulanmadı. Hasta başına ilk deterministik satır tutularak audit havuzu 153 hastaya indirildi.

Panel, CoxPH-v2'nin gerçekten fit edilen feature'ları ile RSF/DeepSurv ortak girdileri ve katı METABRIC eşlemesinin kesişimidir. RSF/DeepSurv için seed=20264917 ile yeni havuzlardan 12+12 satır kimliksiz ve replacement olmadan seçildi. Query matrisleri `results/phase75_fold1_panel/query_matrices.npz` içinde saklıdır.

Formül SHA256 doğrulaması geçti. Gerçek kimlik dışa aktarılmadı. Wall-clock **4.013 saniye**.
