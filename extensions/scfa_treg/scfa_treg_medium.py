"""
scfa_treg_medium.py
====================
Proje 3 (SCFA-Treg) icin kaldigimiz yerden devam: hedefli besiyeri kalibrasyonu.

DURUM: Model dogrulandi, butirat/laktat/O2 exchange reaksiyonlari bulundu
(bkz. project3_scfa_treg_notes.md). Ancak modelin 460 exchange reaksiyonunun
HEPSI sinirsiz (-1000, 1000) -- gercekci olmayan "acik besiyeri". Bunu
kismen kapatip, hipotezle dogrudan ilgili birkac reaksiyonu (glikoz, O2,
laktat, butirat/asetat/propanoat) gercek literatur degerleriyle
sinirlayacagiz, geri kalanini simdilik acik birakacagiz.

EKSIK: Glikoz ve O2 icin GERCEK, mutlak bir alim hizi (mmol/gDW/h ya da
donusturulebilir bir birimde) bulamadim. Asagidaki GLUCOSE_UPTAKE ve
O2_UPTAKE degerleri BILEREK None birakildi -- uydurma bir sayi
yerlestirmedim. Bu ikisini doldurmadan script "gercek" bir sonuc
uretmeyecek, sadece yapiyi kontrol edecek.

NEREDE ARANACAK (en umut verici, henuz okumadigim kaynak):
  PMC5684100 - "Metabolic Adaptation of Human CD4+ and CD8+ T-Cells to
  T-Cell Receptor-Mediated Stimulation" -- bu makalenin Seahorse
  figurlerinde/tablolarinda mutlak OCR (pmol O2/dk) ve ECAR ya da glikoz
  tuketim degerleri olma ihtimali yuksek. Bulursan, hucre basina degeri
  bir T hucresinin gercek kuru agirligiyla (araniyor, ~150-200 pg
  araliginda tahminler var ama dogrulanmis degil) mmol/gDW/h'ye cevir.

  Rathmell et al. 2001 (PMID 11739504, BioNumbers BNID 102644): dinlenmis
  T hucrelerinde ~2-3 nmol O2 / nmol glikoz orani -- bu, GLUCOSE_UPTAKE ve
  O2_UPTAKE'i doldurduktan sonra bir SAGLAMLIK KONTROLU olarak kullanilabilir
  (hesapladigimiz oranin bu araliga yakin cikip cikmadigina bak).
"""
import warnings
warnings.filterwarnings("ignore")
import json
import cobra
import logging
logging.getLogger("cobra").setLevel(logging.ERROR)

# --- reaksiyon ID'leri (dogrulandi, notes dosyasinda kayitli) ---
BUTYRATE_RXN = "HTimmR_8729"   #  <=> m01410s
LACTATE_RXN  = "HTimmR_8835"   #  <=> m02403s
O2_RXN       = "HTimmR_8752"   #  <=> m02630s  (yon: negatif = alim)

# --- HENUZ DOLDURULMADI: gercek, kaynakli deger bulunmadan calistirma ---
GLUCOSE_UPTAKE = None  # mmol/gDW/h, gercek deger lazim (yukaridaki kaynaga bak)
O2_UPTAKE      = None  # mmol/gDW/h, gercek deger lazim

GEMS_DIR = "/home/aceren/diurnal_host_microbiome/data/gems_immune"
BUTYRATE_JSON = "/home/aceren/diurnal_host_microbiome/butyrate_all_replicates_results.json"


def find_exchange_by_metabolite(model, met_id):
    """Reaksiyon adlari anlamsiz oldugu icin metabolit ID'sinden geriye gidiyoruz."""
    met = model.metabolites.get_by_id(met_id)
    for r in met.reactions:
        if r in model.exchanges:
            return r
    return None


def apply_targeted_medium(model, glucose_uptake=None, o2_uptake=None):
    """
    Sadece hipotezle ilgili reaksiyonlari kisitlar, geri kalan 450+ exchange
    reaksiyonunu ACIK birakir. Bu bilerek yapilan bir basitlestirme --
    tam bir besiyeri kalibrasyonu degil, README/abstract'ta acikca
    "hedefli besiyeri, tam kalibre edilmedi" diye belirtilmeli.
    """
    if glucose_uptake is not None:
        glc = find_exchange_by_metabolite(model, "m01965s")  # glukoz ID'sini dogrula
        if glc:
            glc.lower_bound = -abs(glucose_uptake)
    if o2_uptake is not None:
        o2 = model.reactions.get_by_id(O2_RXN)
        o2.lower_bound = -abs(o2_uptake)


def run_pilot(cell_type, butyrate_value, glucose_uptake, o2_uptake):
    model = cobra.io.read_sbml_model(f"{GEMS_DIR}/{cell_type}.xml")
    apply_targeted_medium(model, glucose_uptake, o2_uptake)

    with model:
        model.reactions.get_by_id(BUTYRATE_RXN).lower_bound = -butyrate_value
        model.objective = LACTATE_RXN
        model.objective_direction = "max"
        lac = model.optimize().objective_value

        model.objective = O2_RXN
        model.objective_direction = "min"
        o2 = model.optimize().objective_value

    return {"cell_type": cell_type, "butyrate": butyrate_value, "lactate_max": lac, "o2_min": o2}


def main():
    if GLUCOSE_UPTAKE is None or O2_UPTAKE is None:
        print("DUR: GLUCOSE_UPTAKE / O2_UPTAKE henuz gercek bir kaynaktan doldurulmadi.")
        print("Once literatur degerini bul, dosyanin basindaki sabitleri guncelle, sonra calistir.")
        return

    with open(BUTYRATE_JSON) as f:
        replicates = json.load(f)

    results = []
    for rep in replicates:
        for cell_type in ["Thp", "iTreg", "Th17"]:
            res = run_pilot(cell_type, rep["max_butyrate_mmol_gDW_h"], GLUCOSE_UPTAKE, O2_UPTAKE)
            res.update({"condition": rep["condition"], "ZT": rep["ZT"], "sample": rep["sample"]})
            results.append(res)

    with open("scfa_treg_results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"{len(results)} sonuc yazildi -> scfa_treg_results.json")


if __name__ == "__main__":
    main()
