# Faz 7.5 fold 1 — CoxPH-v2 N=50 audit (TAMAMLANAMADI)

Fold 1 referans CoxPH modeli aynı dondurulmuş `penalizer=0.1` ile yeniden fit edildi; uyarısız test C-index **0.705882** olarak Faz 7.1 sonucunu birebir doğruladı. Fold-spesifik ortak panel 80 feature, TCGA audit havuzu 153 hasta ve METABRIC havuzu 1466 hastadır.

N=50 patient-cluster bootstrap'ın ilk 15 replikası tamamlandı. Sıfır tabanlı replika 15, koşullu sıfır varyans ön taramasında 12 sütun dışlandıktan sonra hem `max_steps=50` ilk denemesinde hem dondurulmuş `max_steps=200, step_size=0.25` retry denemesinde `Newton-Raphson failed to converge sufficiently` uyarısı verdi. Lifelines isimli bir separation sütunu işaretlemedi. Aynı deterministik replika ikinci kez çalıştırılarak sonuç doğrulandı; son deneme 406.7 saniye sürdü.

Replika atlanmadı, başka replikayla değiştirilmedi ve frozen solver protokolü genişletilmedi. Bu nedenle CoxPH için N=50 component/WI sonucu, taze null gate ve BH/BY sertifikasyonu **NA** bırakıldı. Başarısızlık kaydı `results/phase75_fold1_cox/failed_replicate_15.json`, başarıyla tamamlanan önceki replika checkpoint'leri `results/phase75_fold1_cox/checkpoints/` içindedir.

Formül SHA256 değerleri değişmedi; kimlik dışa aktarılmadı.
