# Faz 8 — GitHub deposu özeti

## Depo

- Public URL: https://github.com/umutsaray/warrant-index-survival-xai
- Varsayılan dal: `main`
- İlk yayın commit'i: `f18ab12`
- Lisans: MIT

## Dahil edilenler

- `warrantxai/` WI bileşenleri, agregasyon, sertifikasyon ve API kodu.
- Faz 3–8 dahil tüm `experiments/` analiz scriptleri.
- Birim, regresyon, gizlilik ve artefakt testleri.
- Kimliksiz, feature-düzeyi toplu CSV ve Markdown raporları.
- Yayın şekilleri ve toplu LaTeX tablo parçaları.
- `README.md`, sabit paket sürümlerini içeren `requirements.txt`, `pyproject.toml` ve MIT `LICENSE`.

## Hariç tutulanlar

`.gitignore` ve push-öncesi tracked-file taraması aşağıdaki kategorileri dışlar:

- `data/` ve `results/` altındaki bütün yerel veri/hesap artefaktları.
- Ham TCGA/METABRIC ve cBioPortal yanıtları; hasta/örnek düzeyi parquet dosyaları.
- Fitted model ve checkpoint dosyaları (`.joblib`, `.pkl`, `.pt`, `.pth`).
- Attribution/bootstrap tensörleri (`.npz`, `.npy`).
- Sanal ortamlar, test cache'leri, scratch/temp/log çıktıları.
- `.env`, anahtar, token, credential ve secret isim kalıpları.

İlk commit 210 dosya ve yaklaşık 3.1 MB izlenen içerik taşıdı. Yasak uzantı/yol taraması sıfır eşleşme verdi. GitHub sayfasında public görünürlük, `main` dalı, README ve MIT lisansı push sonrasında doğrulandı.
