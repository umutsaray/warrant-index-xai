# Faz 20 — public depo hazırlığı kapanış raporu

## Özet

Bu faz bilimsel hesap yürütmeden public depo hijyeni, ikili lisanslama, atıf/Zenodo metadatası, README güncellemesi, test doğrulaması ve sürüm etiketi hazırlığını tamamladı.

## ADIM 1 — hasta-düzeyi / yeniden dağıtılamaz veri taraması

- Ayrıntılı kayıt: [`PHASE20_ADIM1_TARAMA.md`](PHASE20_ADIM1_TARAMA.md).
- İzlenen Git ağacında ve erişilebilir geçmişte hasta barkodu, ham kohort matrisi, hasta-düzeyi ikili artefakt, kişisel e-posta, yerel `C:\\Users\\...` yolu veya anahtar/token bulunmadı.
- Hakem metni alıntısı ya da belirli bir hakem/öz-hakemlik/Claude atfı bulunmadı.
- Yerel çalışma ağacında hasta-düzeyi METABRIC cache'leri, Parquet matrisleri, modeller ve atıf tensörleri bulundu; tamamı `.gitignore` ile dışlanmış ve hiçbiri stage/commit edilmemiştir.
- Geçmişi yeniden yazma gereği saptanmadı; force-push/filter-repo yapılmadı.

## ADIM 2 — lisans ve atıf metadata

- Mevcut `LICENSE` MIT olarak korundu.
- `LICENSE-DATA`, Creative Commons'ın resmî `https://creativecommons.org/licenses/by/4.0/legalcode.txt` metniyle normalize satır sonları dışında birebir eşleşecek şekilde eklendi.
- `CITATION.cff`, Umut Saray ve ORCID `0000-0003-3339-6876` ile eklendi; DOI alanı uydurulmadı ve e-posta eklenmedi.
- `.zenodo.json`, açık erişimli yazılım kaydı taslağı olarak eklendi. Yazılım lisansı MIT, türetilmiş depo çıktılarının lisansı açıklamada CC BY 4.0 olarak belirtildi.

## ADIM 3 — README ve MANIFEST

- README'ye v1 resmî, v2 keşifsel ve Faz 18 tolerans/sentetik analiz katmanları; raporlardan alınan sayısal bulgular; yeniden üretim sırası/süreleri; veri erişimi ve ikili lisans kapsamı eklendi.
- ADIM 1 doğruladığı için depoda izlenen hasta-düzeyi veri bulunmadığı açıkça belirtildi.
- `MANIFEST.md`, Faz 17–20 dosyaları ve SHA256 tablosuyla güncellendi.
- Depoda bağımsız `PHASE19_PREREG.md`, `PHASE19_KAPANIS.md` veya Faz 19 betiği bulunmadığından Phase 19 içeriği uydurulmadı; bu eksiklik README ve MANIFEST'te açıkça kaydedildi.

## ADIM 4 — sürüm etiketi ve Zenodo hazırlığı

- Anotasyonlu etiket: `v1.0-bmi-submission`.
- Etiket mesajı: `Code and derived results accompanying the BioMedInformatics submission`.
- Etiketin işaret ettiği commit SHA: `3904910e9d2a6a6bc8d891652373ba9fc5ca628c`.
- GitHub release oluşturulmadı. Zenodo DOI'si henüz mevcut olmadığından DOI ve README rozeti eklenmedi.

## Testler ve biçim doğrulamaları

- Sistem Python 3.14 ortamında `pytest` kurulu değildi; bu deneme test çalıştırmadı.
- Proje sanal ortamındaki ilk tam test denemesinde 210 test geçti, bir test Windows genel temp klasörü izin hatası nedeniyle setup aşamasında kaldı.
- Aynı paket çalışma alanı içindeki izleme-dışı `--basetemp=.pytest-phase20-final` ile yeniden çalıştırıldı: **211/211 geçti (3.73 s)**.
- `LICENSE-DATA` resmî kaynakla eşleşti; `.zenodo.json` geçerli JSON, `CITATION.cff` geçerli YAML olarak doğrulandı; `git diff --check` temizdi.

## Değişen/eklenen dosyalar

- Değişen: `README.md`, `MANIFEST.md`.
- Eklenen: `LICENSE-DATA`, `CITATION.cff`, `.zenodo.json`, `PHASE20_ADIM1_TARAMA.md`, `PHASE20_KAPANIS.md`.
- Bilimsel raporlar, `results/**`, Tablo 3–9 kaynakları, model/checkpoint dosyaları ve saklı atıf tensörleri değiştirilmedi.

## Formül dosyası SHA256 başlangıç/bitiş

| Dosya | Başlangıç | Bitiş | Durum |
|---|---|---|---|
| `warrantxai/components.py` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` | `7887CF447958726DF0FC2FD17FFC5BC4D70B00DE92B0888F36D3C4F49057E7AA` | değişmedi |
| `warrantxai/aggregation.py` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` | `10179864A63104ECA5AEEC8145983F0D25B2539A191A5B2E2095EB4B3BF5FF68` | değişmedi |
| `warrantxai/certification.py` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` | `76AB49716CC3EAD6DEC71656BBD60C6D71C1D76CAB8A83B2F0662AF6D6A05487` | değişmedi |

## Commit güvenliği

Dosyalar yalnız açık dosya-adı listeleriyle stage edildi. `git add .`, `git add -A`, `git add -f` veya wildcard stage kullanılmadı. Her commit öncesinde staged dosya listesi kontrol edildi; `.parquet`, `.npy`, `.npz`, `.pkl`, model veya cache yolu bulunmadı.

Bu faz bilimsel hesap içermemektedir; dondurulmuş ayarlar, Tablo 3–9, BH/BY kararları, saklı atıf tensörleri ve CoxPH-v2 retry protokolü değiştirilmemiştir.
