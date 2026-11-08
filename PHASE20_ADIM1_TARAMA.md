# Faz 20 — ADIM 1 hasta-düzeyi / yeniden dağıtılamaz veri taraması

## Kapı kararı

**BULGU VAR — ADIM 2–5 başlatılmadı.** Kamuya açık Git ağacında veya Git geçmişinde hasta barkodu, ham kohort matrisi, model/atıf ikilisi ya da açık bir kimlik bilgisi bulunmadı. Bununla birlikte yerel çalışma ağacında Git tarafından izlenmeyen ve `.gitignore` ile dışlanan hasta-düzeyi/yeniden dağıtımı kısıtlı olabilecek dosyalar vardır. Bu dosyalar silinmedi, taşınmadı, içerikleri rapora kopyalanmadı ve Git geçmişi yeniden yazılmadı.

Tarama başlangıcı: commit `41ca39d1ee8afba72b3b8801c8c7a44ae21fff27`, dal `main`, 320 izlenen dosya. Başlangıç çalışma ağacı temizdi.

## Tarama yöntemi

- `git ls-files`, çalışma ağacı ve `git log --all --diff-filter=A --name-only` tarandı.
- İstenen ikili/kap formatları, 5 MB üzerindeki dosyalar, hassas klasör adları ve hasta kimliği örüntüleri kontrol edildi.
- CSV dosyalarında yalnız başlık ve satır sayısı incelendi; satır içeriği rapora alınmadı.
- Parquet dosyalarında yalnız şema/sütun adları; NPZ dosyalarında yalnız anahtar ve dizi biçimleri incelendi.
- Yerel mutlak yol, e-posta, olası anahtar/token ve hakem/öz-hakemlik atıfları hem mevcut Git ağacında hem tüm erişilebilir commit'lerde tarandı. Eşleşme değerleri rapora alınmadı.

## Git tarafından izlenen içerik

| Denetim | Sonuç | Değerlendirme |
|---|---:|---|
| İstenen ikili uzantılar | 0 | İzlenen `.npy/.npz/.pkl/.parquet/.joblib/.pt/...` yok. |
| `TCGA-XX-XXXX` veya `MB-####` örüntüsü | 0 | İzlenen içerikte hasta barkodu eşleşmesi yok. |
| Kimlik benzeri CSV başlığı | 0 / 65 CSV | `patient_id`, `case_id`, `sample_id`, `submitter_id`, `barcode` yok. |
| 5 MB üzeri izlenen dosya | 2 | Aşağıdaki iki dosya türetilmiş, hasta kimliği içermeyen tablolardır. |

| Yol | Boyut | Satır | Başlık özeti | Hasta-düzeyi bilgi |
|---|---:|---:|---|---|
| `results/phase16/PHASE16_ADIM1_FULL_POOL_RANKING.csv` | 14,921,145 bayt | 60,660 | Özellik ve sıralama/istatistik alanları | HAYIR — özellik-düzeyi türetilmiş sıralama tablosu |
| `results/phase18/PHASE18_ADIM2_POINT_GRID.csv` | 8,827,104 bayt | 52,200 | Seed, sentetik truth sınıfı, varyant ve bileşen skorları | HAYIR — sentetik/türetilmiş grid |

Kod içinde `patient_id`/`sample_id` gibi sözcükler 60 konumda kullanılıyor; bunlar deney/test mantığı ve şema kontrolleridir, veri kaydı değildir. İzlenen içerikte gerçek barkod örüntüsü bulunmadığından dosya/satır listesi hassas bulgu olarak raporlanmadı.

## Yerel, Git tarafından izlenmeyen bulgular

### Ham/önbellek kohort verileri

| Yol/grup | Dosya | Toplam boyut | Hasta-düzeyi bilgi | Neden |
|---|---:|---:|---|---|
| `data/external/metabric/` | 15 | 84,854,895 bayt | EVET | Hasta ve örnek klinik JSON'ları, örnek listeleri ve mRNA parti yanıtları içeriyor. |
| `data/external/metabric_phase72/` | 11 | 52,468,277 bayt | EVET | `patient_clinical.json`, `sample_clinical.json`, `samples.json` ve mRNA parti önbellekleri içeriyor. |
| `data/external/metabric_phase75_fold1/` | 12 | 86,643,736 bayt | EVET | Hasta/örnek klinik JSON'ları ve fold-1 mRNA parti önbellekleri içeriyor. |

Bu üç grup `.gitignore` içindeki `/data/` kuralıyla dışlanmaktadır. Dosya içerikleri veya hasta kimlikleri rapora kopyalanmamıştır.

### Hasta matrisleri, modeller ve saklı atıf tensörleri

| Yol/grup | Dosya sayısı | Toplam boyut | Hasta-düzeyi bilgi | Neden |
|---|---:|---:|---|---|
| `results/phase7_retrain/` (`*.parquet`, `*.joblib`, `*.pt`) | 10 | 572,845,260 bayt | EVET/BELİRSİZ | Parquet şemalarında `sample_key` ve hasta satırları var; model dosyaları eğitim verisinden türetilmiş fakat doğrudan kimlik içerip içermediği güvenle kanıtlanmadı. |
| `results/phase70b_collapse_diagnostic/` (fold Parquet/model dosyaları) | 20 | 10,422,398 bayt | EVET/BELİRSİZ | Fold matrisleri `sample_key` içeriyor; model dosyaları için yeniden dağıtım/kimlik riski belirsiz. |
| `results/phase75_fold1/` (Parquet/modeller) | 8 | 343,264,722 bayt | EVET/BELİRSİZ | Hasta-düzeyi fold matrisi ve eğitilmiş model artefaktları. |
| `results/phase72_cox/attributions_n50.npz` | 1 | 32,908,069 bayt | EVET | Sorgu/hasta düzeyinde saklı atıf tensörü; kimlik etiketi taşımadığı görülse de hasta verisinden türetilmiş. |
| `results/phase72_panel/`, `results/phase75_fold1_panel/` | 2 | 1,191,928 bayt | EVET | Sorgu matrisleri; doğrudan kimlik etiketi olmadan hasta-düzeyi değerler içeriyor. |
| Diğer `results/**` NPZ checkpoint/atıf grupları | 566 | yaklaşık 13 MB | EVET/BELİRSİZ | Candidate/reference atıfları, katsayılar veya bootstrap dizileri; çoğu hasta/sorgu düzeyinde türetilmiş, kimlik alanı görünmüyor. |
| `.pytest*` ve Faz 17–18 test önbellekleri | 58 | 12,992 bayt | HAYIR | 224 baytlık sentetik test önbellekleri. |
| `dist/warrantxai-0.1.0.tar.gz` | 1 | 48,322 bayt | HAYIR | Kaynak dağıtım paketi; `/dist/` kuralıyla dışlanmış. |

Tüm `results/**` bulguları `.gitignore` içindeki `/results/` kuralıyla; model ve dizi türleri ayrıca `*.pt`, `*.joblib`, `*.npz`, `*.npy` kurallarıyla dışlanmaktadır. Öneri: bu artefaktları Git'e eklemeyin; yayımlanacak türetilmiş tablolar için ayrı, kimliksiz ve lisans açısından doğrulanmış bir allowlist kullanın. Model artefaktlarının yayınlanması düşünülürse önce model inversion/memorization ve kaynak veri kullanım koşulları bakımından ayrıca değerlendirin.

## Git geçmişi

`git log --all --diff-filter=A --name-only` taramasında istenen ikili uzantılarla daha önce eklenip silinmiş dosya bulunmadı. `data/`, `raw/`, `.cache/` veya `metabric_phase72/` altında geçmişte izlenmiş veri dosyası bulunmadı. Yalnız veri indirme/manifest betik adları eşleşti:

- `experiments/06_download_metabric.py`
- `experiments/06_manifest_tcga.py`
- `experiments/38_phase72_download_metabric.py`
- `experiments/47_phase75_fold1_download_metabric.py`
- `experiments/70_phase15_metabric_standardization.py`

Bunlar kod dosyalarıdır; hasta verisi değildir. Bu nedenle mevcut bulgular için geçmişi yeniden yazma gereği saptanmadı.

## Gizlilik, kimlik bilgisi ve süreç dili taraması

| Örüntü sınıfı | Mevcut izlenen ağaç | Git geçmişi |
|---|---:|---:|
| Yerel mutlak yol (`C:\\Users\\`) | 0 | 0 |
| Kişisel e-posta biçimi | 0 | 0 |
| `ghp_`, `AKIA`, `password`, `secret` | 0 | 0 |
| Hakem/öz-hakemlik/Claude atıf örüntüleri | 0 | 0 |

Yerel mutlak yol (`C:\\Users\\...`) eşleşmesi bulunmadı; sonuç hem mevcut izlenen ağaç hem de erişilebilir Git geçmişi için sıfırdır.

Hakem metni alıntısı veya bir hipotezi/itirazı belirli bir hakeme, öz-hakemliğe ya da Claude'a atfeden ifade bulunmadı; sonuç hem mevcut izlenen ağaç hem de erişilebilir Git geçmişi için sıfırdır.

## Formül dosyası bütünlüğü — başlangıç

| Dosya | SHA256 |
|---|---|
| `warrantxai/components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` |
| `warrantxai/aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` |
| `warrantxai/certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` |

## Önerilen karar

Kamuya açık Git deposunun mevcut izlenen içeriğinde acil kaldırma/geçmiş yeniden yazma gerektiren bir bulgu yoktur. Buna karşın yerel hasta-düzeyi dosyaların varlığı nedeniyle Faz 20 talimatındaki kapı uygulanmıştır. ADIM 2–5'e devam edilmeden önce kullanıcı onayı alınmalıdır. Devam onayı verilirse yalnız izlenen kaynak/dokümantasyon dosyaları değiştirilmeli; `/data/`, `/results/`, saklı atıf tensörleri ve modeller kesinlikle stage edilmemelidir.
