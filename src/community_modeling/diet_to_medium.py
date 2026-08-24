"""
diet_to_medium.py
==================

Converts the study's diet/feeding conditions (NA, FA, FT) into a time-resolved
"medium" definition: a dict of {exchange_reaction_id: upper_bound (mmol/gDW/h)}
suitable for constraining a cobrapy/pyTFA community model's model.medium at a
given Zeitgeber time t.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field

# TestDiet 58Y1 -- DIO Purified Diet, 60% kcal fat (fully quantified).
# VERIFIED against user-provided source files (2026-08-14):
#   - Diet_Reference_Table.xlsx, sheet "Diets", rows under
#     "1. TestDiet 58Y1 (HFD, ad libitum / TRF - FA & FT groups)"
#   - Diet for Project.xlsx, sheet "Sayfa1", rows for "TestDiet 58Y1"
# All 14 ingredient percentages below match both files exactly (cross-checked
# programmatically, not from memory).
DIET_58Y1 = {
    "Lard": 31.66,
    "Casein_VitaminTested": 25.845,
    "Maltodextrin": 16.153,
    "Sucrose": 8.847,
    "Powdered_Cellulose": 6.461,
    "Soybean_Oil": 3.231,
    "Potassium_Citrate": 2.132,
    "Calcium_Phosphate": 1.68,
    "DIO_Mineral_Mix": 1.292,
    "AIN76A_Vitamin_Mix": 1.292,
    "Calcium_Carbonate": 0.711,
    "L_Cystine": 0.388,
    "Choline_Bitartrate": 0.258,
    "FDC_Blue1": 0.05,
}

# LabDiet 5001 -- Natural-ingredient chow (NCD). Exact ingredient % is not
# published by the manufacturer (confirmed: "Diet for Project.xlsx" lists
# all ingredients with Amount="N/A"). Using the "Nutrient Basis" quantitative
# proxy instead (protein/fat/fiber/starch+sugar/ash from the Guaranteed
# Analysis table).
# VERIFIED against Diet_Reference_Table.xlsx, sheet "Diets", Section 3
# ("LabDiet 5001 (Nutrient Basis)") -- all 6 values match exactly.
# DECISION: the same reference file's README also offers TestDiet 58Y2 as a
# fully-quantified alternative, but explicitly recommends AGAINST it when
# replicating the original study matters more than having exact quantities
# for every diet. Since our goal is reproducing Ramos et al. 2025's actual
# feeding conditions (which used real LabDiet 5001, not 58Y2), we keep the
# Nutrient Basis proxy -- this is the file's own recommended choice for our
# use case, not an arbitrary pick.
DIET_5001_APPROX = {
    "Protein_mixed": 24.1,
    "Fat_mixed": 5.1,
    "Crude_Fiber": 5.3,
    "Starch_Sugar": 21.9 + 3.25,
    "Ash_Minerals": 7.2,
    "Unresolved_NFE": 100 - (24.1 + 5.1 + 5.3 + 21.9 + 3.25 + 7.2),
}

CLASS_MAP = {
    "Lard": "fat", "Soybean_Oil": "fat", "Fat_mixed": "fat",
    "Casein_VitaminTested": "protein", "L_Cystine": "protein", "Protein_mixed": "protein",
    "Maltodextrin": "digestible_cho", "Sucrose": "digestible_cho",
    "Starch_Sugar": "digestible_cho", "Unresolved_NFE": "digestible_cho",
    "Powdered_Cellulose": "fiber", "Crude_Fiber": "fiber",
    "Potassium_Citrate": "mineral", "Calcium_Phosphate": "mineral",
    "DIO_Mineral_Mix": "mineral", "Calcium_Carbonate": "mineral",
    "Ash_Minerals": "mineral",
    "AIN76A_Vitamin_Mix": "vitamin",
    "Choline_Bitartrate": "other", "FDC_Blue1": "other",
}

ABSORPTION_FRACTION = {
    "fat": 0.92,
    "protein": 0.90,
    "digestible_cho": 0.97,
    "fiber": 0.02,
    "mineral": 0.55,
    "vitamin": 0.70,
    "other": 0.50,
}

CLASS_MW = {
    "fat": 284.0,
    "protein": 110.0,
    "digestible_cho": 180.0,
    "fiber": 162.0,
    "mineral": 60.0,
    "vitamin": 300.0,
    "other": 150.0,
}

# NOTE (fixed): BiGG/CarveMe uses a DOUBLE underscore before the
# stereo-descriptor (e.g. "ala__L", not "ala_L"). The single-underscore
# forms used previously never matched real model exchange IDs -- confirmed
# via empirical exchange-ID inspection of the WoL2 community model.
CLASS_EXCHANGES = {
    "fat": ["EX_hdca_e", "EX_ocdca_e"],
    "protein": [
        "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e",
        "EX_cys__L_e", "EX_glu__L_e", "EX_gln__L_e", "EX_gly_e",
        "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
        "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e",
        "EX_thr__L_e", "EX_trp__L_e", "EX_tyr__L_e", "EX_val__L_e",
    ],
    "digestible_cho": [],
    "fiber": ["EX_cellb_e", "EX_glc__D_e"],
    "mineral": ["EX_ca2_e", "EX_pi_e", "EX_k_e", "EX_fe2_e"],
    "vitamin": ["EX_thm_e", "EX_ribflv_e", "EX_fol_e"],
    "other": [],
}


@dataclass
class FeedingCondition:
    name: str
    diet: dict
    daily_intake_g: float
    trf_window: tuple | None = None
    cecal_biomass_gDW: float = 0.05


CONDITIONS = {
    "NA": FeedingCondition("NA (LabDiet 5001, ad libitum)", DIET_5001_APPROX,
                            daily_intake_g=4.5, trf_window=None),
    "FA": FeedingCondition("FA (TestDiet 58Y1, ad libitum)", DIET_58Y1,
                            daily_intake_g=3.0, trf_window=None),
    "FT": FeedingCondition("FT (TestDiet 58Y1, 8h TRF)", DIET_58Y1,
                            daily_intake_g=3.0, trf_window=(13, 21)),
}


def cecal_class_flux(cond: FeedingCondition) -> dict:
    class_mass_g = {}
    for ingredient, pct in cond.diet.items():
        cls = CLASS_MAP.get(ingredient, "other")
        mass_g = cond.daily_intake_g * (pct / 100.0)
        class_mass_g[cls] = class_mass_g.get(cls, 0.0) + mass_g

    flux = {}
    for cls, mass_g in class_mass_g.items():
        surviving_g = mass_g * (1 - ABSORPTION_FRACTION.get(cls, 0.5))
        mmol = (surviving_g / CLASS_MW.get(cls, 150.0)) * 1000.0
        flux[cls] = mmol / 24.0 / cond.cecal_biomass_gDW
    return flux


def trf_gate(t_zt: float, window: tuple, steepness: float = 4.0) -> float:
    start, end = window
    t = t_zt % 24.0

    def sigmoid(x):
        return 1.0 / (1.0 + math.exp(-steepness * x))

    if start < end:
        return sigmoid(t - start) * sigmoid(end - t)
    else:
        return max(sigmoid(t - start), sigmoid(end - t + 24 if t < start else end - t))


def build_medium(condition_key: str, t_zt: float) -> dict:
    cond = CONDITIONS[condition_key]
    class_flux = cecal_class_flux(cond)
    gate = 1.0 if cond.trf_window is None else trf_gate(t_zt, cond.trf_window)
    medium = {}
    for cls, mmol_gDW_h in class_flux.items():
        for rxn in CLASS_EXCHANGES.get(cls, []):
            medium[rxn] = round(mmol_gDW_h * gate, 5)
    return medium


if __name__ == "__main__":
    for cond_key in CONDITIONS:
        print(f"\n{CONDITIONS[cond_key].name}")
        for zt in [1, 13, 17, 21]:
            m = build_medium(cond_key, zt)
            print(f"  ZT{zt:>2}: n_active_exchanges={sum(1 for v in m.values() if v > 1e-6)}")
