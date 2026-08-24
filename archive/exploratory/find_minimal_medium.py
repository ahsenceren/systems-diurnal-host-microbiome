"""
find_minimal_medium.py
========================

Uses cobra's own minimal_medium() to find the SMALLEST set of exchange
reactions (beyond our core+diet set) needed for Muribaculaceae and
Faecalibacterium to grow -- algorithmic, not guessed.
"""

import cobra
from cobra.medium import minimal_medium
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
    ex_ids = [r.id for r in model.exchanges]

    for r in model.exchanges:
        r.lower_bound = -5.0
        r.upper_bound = 1000
    for bigg_id, val in {**CORE_UNIVERSAL_BIGG, **diet_medium_bigg}.items():
        if bigg_id in ex_ids:
            model.reactions.get_by_id(bigg_id).lower_bound = -val

    try:
        mm = minimal_medium(model, minimize_components=True, open_exchanges=False)
        print("Minimal medium (all compounds needed for growth):")
        print(mm)
        print()
        already_have = set(CORE_UNIVERSAL_BIGG.keys()) | set(diet_medium_bigg.keys())
        new_needed = {k: v for k, v in mm.items() if k not in already_have}
        print(f"NEW compounds needed (not already in core+diet): {len(new_needed)}")
        for k, v in new_needed.items():
            met = model.reactions.get_by_id(k).reactants[0] if model.reactions.get_by_id(k).reactants else None
            metname = met.name if met else "?"
            print(f"  {k} ({metname}): {v}")
    except Exception as e:
        print(f"minimal_medium failed: {e}")
    print()
