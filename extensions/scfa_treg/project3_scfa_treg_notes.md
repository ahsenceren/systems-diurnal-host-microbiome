# Proje 3: SCFA-Treg Metabolik Entegrasyonu — Çalışma Notları

Bu dosya, abstract/ek özet yazımında kullanılmak üzere şimdiye kadar doğrulanmış gerçek bulguları tutuyor. Her madde gerçekten çalıştırılmış bir kontrolün sonucu, tahmin değil.

## Kullanılan model

- **HTimmR** (human T-immuno reconstructor): jenerik insan CD4+ T hücre GEM'i. BioModels erişim kodu **MODEL2101270002**, Cell Reports (Kasım 2021) makalesinin ek materyali. 7558 reaksiyon, 8 kompartman (extracellular, peroxisome, mitochondria, cytosol, lysosome, ER, Golgi, nucleus).
- Bağlam-spesifik türetilmiş, doğrulanmış alt-tip modelleri indirildi ve `cobrapy` ile başarıyla yüklendi:
  - **Thp.xml** (T-naive): 5778 reaksiyon, 6751 metabolit
  - **iTreg.xml**: 5884 reaksiyon, 6855 metabolit
  - **Th17.xml**: 5633 reaksiyon, 7204 metabolit
- Reaksiyon sayıları makalenin bildirdiği HTimmR sayısıyla (7558) birebir eşleşti — dosyaların doğru/bozulmamış olduğunu doğruladık.
- Not: bu dosyalarda reaksiyon `name`/`id` alanları anlamsız ("HTimmR_XXXX" placeholder), gen ilişkileri (GPR) cobrapy tarafından okunamadı (0 gen). Bu yüzden her arama **metabolit isimleri üzerinden** yapılmalı, reaksiyon ismi/ID'si üzerinden değil.

## Bütirat giriş noktası (doğrulandı)

- Metabolit ID'si: **m01410** = butyrate (kompartmanlar: c=sitozol, s=hücre-dışı/sınır).
- iTreg modelinde m01410s ile ilişkili 13 reaksiyon bulundu.
- Gerçek exchange (sınır) reaksiyonu: **HTimmR_8729**, ` <=> m01410s`, `exchange=True`, mevcut sınırlar (-1000, 1000) — kısıtlanmamış. Bunu gerçek topluluk-modeli bütirat çıktımızla (FA/FT/NA, mmol/gDW/h) kısıtlayacağız.
- Asetat (m01252) ve propanoat (m02772) da modelde mevcut — ileride genişletme için not.

## Kritik metodolojik bulgu: biyokütle/büyüme hedefi YOK

- Thp/iTreg/Th17 modellerinin hiçbirinde tanımlı bir amaç fonksiyonu (`objective_coefficient`) ya da "biomass" içeren bir metabolit/reaksiyon bulunamadı — hem `objective_coefficient` kontrolüyle hem metabolit-isim taramasıyla doğrulandı, ikisi de boş çıktı.
- Orijinal makalenin (bioRxiv 10.1101/2021.01.29.428853, tam metin) doğrudan alıntısı: modelin amaç fonksiyonu olarak biyokütle değil, **8 farklı sfingolipid yolağının akısı teker teker maksimize edilmiş** ("flux through 8 different sphingolipid pathways... were maximized one-by-one (as the objective function)"), sonuçlar yüzdeye çevrilip karşılaştırılmış.
- **Sonuç:** Biz de aynı, makalede zaten doğrulanmış yöntemi izleyeceğiz — kendi hipotezimize uygun spesifik reaksiyonları (OXPHOS için bir ATP-sentaz/oksidatif fosforilasyon proxy'si, glikoliz için laktat ihracatı) bulup, bütirat kısıtı altında teker teker maksimize edeceğiz. Rastgele/uydurma bir "maximize biomass" hedefi kullanmayacağız — bu hem yanlış olurdu hem de orijinal makaleyle tutarsız olurdu.

## Sıradaki adım (henüz yapılmadı)

OXPHOS ve glikoliz proxy reaksiyonlarını (laktat, ATP/O2 ile ilişkili metabolitler üzerinden) bulmak için tarama çalıştırıldı, sonucu henüz almadık.
