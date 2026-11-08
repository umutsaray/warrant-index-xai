# Faz 21 ADIM 1 — Wald taraması adli incelemesi

Fold-0 train içindeki 60,660 RNA sütunu incelendi. Adli tablo, orijinal Wald ile seçilen 200 RNA için standartlaştırılmış tek değişkenli lifelines Cox β/SE/z, LR p, Breslow skor p, uyarı/iterasyon, ayrışma göstergeleri ve ridge-cezalı (`penalizer=0.01`; Firth değildir) karşılaştırmayı içerir.

## Özet

- Başarısız unpenalized fit: **87/200**.
- Lifelines ayrışma/low-variance uyarısı: **10/200**.
- Nonzero satırların tamamı event veya tamamı censored: **4/200**.
- ≥%10 nonzero ve nonzero satırlarda ≥10 event filtresini geçen evren: **43877/60,660**.

Özellikle istenen fold-0 22 RNA, yedi RSF sıfır-atıflı RNA ve KRTAP11-1 satırları ana CSV'de boolean bayraklarla işaretlidir; birleşik işaretli satır sayısı 22'dir.

## Alternatif 200'lük taramalar

| method | selected_n | zero_rate_min | zero_rate_q25 | zero_rate_median | zero_rate_q75 | zero_rate_max | metabric_common_n |
|---|---|---|---|---|---|---|---|
| original_Wald | 200 | 0.0 | 0.9477 | 0.98368 | 0.992582 | 0.997033 | 22.0 |
| filtered_LR | 200 | 0.0 | 0.0 | 0.0 | 0.106825 | 0.87092 | 173.0 |
| filtered_score | 200 | 0.0 | 0.0 | 0.236647 | 0.587908 | 0.893175 | 141.0 |

METABRIC ortak sayısı, MyGene Ensembl→tekil sembol eşlemesi ve yerel cBioPortal METABRIC gen kataloğunda tekil Entrez karşılığı bulunması ölçütüdür. Bu, genin indirilebilir ortaklığını gösterir; yeni ekspresyon matrisi bu adımda indirilmedi ve model yeniden eğitilmedi.

LR ve skor panellerinin kesişimi **112/200**; birleşimi 288 özelliktir. Ham sıralamalar `PHASE21_ADIM1_ALTERNATIVE_PANELS.csv` içindedir.

Ortak panel bileşimi: filtre+LR için **173 RNA + 5 klinik = 178**; filtre+skor için **141 RNA + 5 klinik = 146**. Tam özellik listesi `PHASE21_ADIM1_COMMON_PANEL.csv` içindedir.

## Karar kapısı

ADIM 2 başlatılmadı. Düzeltilmiş boru hattı için aday toplam ortak panel boyutları: filtre+LR **178**, filtre+skor **146**. Kullanıcı onayı gereklidir.

Doğrulama: iki Faz 21 betiği `py_compile` kontrolünü, tam paket `211 passed` sonucunu verdi. Dondurulmuş `components.py`, `aggregation.py` ve `certification.py` dosyalarında değişiklik yoktur.
