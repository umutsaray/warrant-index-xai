# Faz 10 ADIM 1 — sayısal doğrulama

Bu doğrulama yalnız `PHASE10_QUERY_RESAMPLING_SENSITIVITY.csv` dosyasındaki 378 satırı okur. Yeni model fit'i, bootstrap, query örneklemesi veya sertifikasyon çalıştırılmamıştır.

## Point-adequacy değişimi

Her RSF özelliğinin 24+24 ve 48+48 sonuçları, aynı özelliğin özgün 12+12 `point_adequate` değeriyle doğrudan karşılaştırıldı. Adequacy, saklı pipeline çıktısındaki `WI >= 0.8` ve bileşen eşiklerinin birlikte sağlanması tanımını kullanır.

**Nihai sonuç: 27 özelliğin 10'unda en az bir yeni konfigürasyonda point-adequacy durumu değişmiştir.**

| Özellik | 12+12 durumu | Referanstan farklı konfigürasyonlar |
| --- | --- | --- |
| OR6C1 | adequate | 24+24 seed 1, 24+24 seed 2: inadequate |
| MRGPRX1 | inadequate | 24+24 seed 3, 48+48 seed 3: adequate |
| LHFPL3 | inadequate | 24+24 seed 3: adequate |
| CACNG3 | adequate | 24+24 seed 1–3, 48+48 seed 1–2: inadequate |
| STXBP5-AS1 | adequate | 24+24 seed 1: inadequate |
| NPAP1 | adequate | 24+24 seed 1–3, 48+48 seed 1–3: inadequate |
| ACCSL | adequate | 24+24 seed 1–3, 48+48 seed 1–2: inadequate |
| LCE3C | adequate | 24+24 seed 1–3, 48+48 seed 1–2: inadequate |
| Age at index | inadequate | 48+48 seed 3: adequate |
| AJCC Stage IV | adequate | 48+48 seed 1: inadequate |

## Darboğaz değişimi

Her yeni konfigürasyonun `bottleneck_component` değeri aynı özelliğin özgün 12+12 darboğazıyla karşılaştırıldı. Saklı çıktı, genel-genellenebilirlik girdisini `G_claim = min(G, G_signed)` olarak adlandırdığı için tabloda `G_claim` kullanılır.

**Nihai sonuç: 27 özelliğin 5'inde en az bir yeni konfigürasyonda darboğaz değişmiştir.**

| Özellik | 12+12 darboğazı | Referanstan farklı konfigürasyonlar |
| --- | --- | --- |
| MRGPRX1 | G_claim | 24+24 seed 1–2, 48+48 seed 1–2: R |
| SMR3A | R | 24+24 seed 3: G_claim |
| CA1 | G_claim | 24+24 seed 1–3, 48+48 seed 1–3: R |
| STXBP5-AS1 | G_claim | 48+48 seed 3: R |
| Age at index | R | 24+24 seed 1 ve 3, 48+48 seed 1: G_claim |

## Tutarsızlığın kök nedeni

`8/4` sayımı özellik kimliği yerine `symbol` alanına göre gruplanmıştı. Klinik özelliklerin `symbol` alanı ilk ara tabloda `NaN` olduğundan pandas `groupby` işlemi `Age at index` ile `AJCC Stage IV` özelliklerini varsayılan olarak dışladı. `Age at index` hem point hem darboğaz sayımına, `AJCC Stage IV` yalnız point sayımına eklenince sonuç `8/4` yerine **10/5** olur. Seed'leri ayrı özellikler gibi sayma veya başka bir adequacy tanımı söz konusu değildir; `8/4` eksik-gruplama kaynaklı bir özet hatasıdır.

## Manuskript önerisi

Limitations bölümünde kullanılacak doğrulanmış ifade: **“RSF'de 10/27 özelliğin point-adequacy durumu ve 5/27 özelliğin darboğaz bileşeni en az bir sorgu-boyutu/seed konfigürasyonunda özgün 12+12 sonucundan farklıydı.”**

Bu düzeltme WI formülünü, eşikleri veya Table 3–8'deki BH/BY sertifikasyon kararlarını etkilemez. `components.py` SHA256 değeri `7887cf447958726df0fc2fd17ffc5bc4d70b00de92b0888f36d3c4f49057e7aa`, `aggregation.py` SHA256 değeri `10179864a63104eca5aeec8145983f0d25b2539a191a5b2e2095eb4b3bf5ff68` olarak değişmeden doğrulanmıştır.
