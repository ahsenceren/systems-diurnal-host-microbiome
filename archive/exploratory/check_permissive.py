"""
check_permissive_and_combo.py
================================

No single exchange unlocked growth -- so it's either (a) a genuine GEM
limitation (can't grow under ANY medium, meaning a curation problem), or
(b) a multi-nutrient combination requirement. Test both quickly.
"""

import cobra
import sys

sys.path.insert(0, "/home/aceren/diurnal_host_microbiome")
from diet_to_medium import build_medium

MODELS = {
    "Muribaculaceae": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
    "Faecalibacterium": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
}

CORE_UNIVERSAL_BIGG = {
    "EX_h2o_e": 100.0, "EX_h_e": 100.0, "EX_pi_e": 10.0, "EX_co2_e": 50.0,
    "EX_cl_e": 10.0, "EX_na1_e": 10.0, "EX_mg2_e": 5.0, "EX_mn2_e": 1.0,
    "EX_zn2_e": 1.0, "EX_cu2_e": 1.0, "EX_fe3_e": 1.0, "EX_so4_e": 5.0,
    "EX_nh4_e": 10.0, "EX_cobalt2_e": 0.5, "EX_ni2_e": 0.5,
}
diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}

for name, path in MODELS.items():
    print("=" * 70)
    print(name)
    print("=" * 70)
    model = cobra.io.read_sbml_model(path)

    with model:
        for r in model.exchanges:
            r.lower_bound = -1000
            r.upper_bound = 1000
        g_permissive = model.slim_optimize()
        print(f"  Growth under FULLY PERMISSIVE medium: {g_permissive}")

    with model:
        for r in model.exchanges:
            r.lower_bound = 0
            r.upper_bound = 1000
        ex_ids = [r.id for r in model.exchanges]
        for bigg_id, val in {**CORE_UNIVERSAL_BIGG, **diet_medium_bigg}.items():
            if bigg_id in ex_ids:
                model.reactions.get_by_id(bigg_id).lower_bound = -val * 10
        g_scaled = model.slim_optimize()
        print(f"  Growth with core+diet bounds x10: {g_scaled}")

    with model:
        for r in model.exchanges:
            r.lower_bound = -5.0
            r.upper_bound = 1000
        ex_ids = [r.id for r in model.exchanges]
        for bigg_id, val in {**CORE_UNIVERSAL_BIGG, **diet_medium_bigg}.items():
            if bigg_id in ex_ids:
                model.reactions.get_by_id(bigg_id).lower_bound = -val
        g_combo = model.slim_optimize()
        print(f"  Growth with core+diet + ALL others open at -5: {g_combo}")

    print()
