"""
find_blocking_nutrient.py
============================

Muribaculaceae and Faecalibacterium show growth_rate=0 under the current
CORE+diet medium (only Roseburia grows). Instead of guessing what's
missing, this tests EACH currently-closed exchange reaction one at a time
(open it briefly, check if growth becomes possible, close it again) on the
INDIVIDUAL (non-community) GEMs, to find exactly which real nutrient(s)
are blocking growth. Fast: single-species FBA solves are cheap.
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
    ex_ids = [r.id for r in model.exchanges]

    for r in model.exchanges:
        r.lower_bound = 0
        r.upper_bound = 1000

    applied = 0
    for bigg_id, val in {**CORE_UNIVERSAL_BIGG, **diet_medium_bigg}.items():
        if bigg_id in ex_ids:
            model.reactions.get_by_id(bigg_id).lower_bound = -val
            applied += 1

    baseline = model.slim_optimize()
    print(f"  Medium applied: {applied} exchanges open")
    print(f"  Baseline growth on this medium: {baseline}")

    if baseline is not None and baseline > 1e-6:
        print("  Already grows -- no blocker.")
        continue

    print("  Growth is 0 -- testing which closed exchange unlocks it...")
    unlockers = []
    open_set = set(CORE_UNIVERSAL_BIGG.keys()) | set(diet_medium_bigg.keys())
    for r in model.exchanges:
        if r.id in open_set:
            continue
        old_lb = r.lower_bound
        r.lower_bound = -10.0
        g = model.slim_optimize()
        r.lower_bound = old_lb
        if g is not None and g > 1e-6:
            met = list(r.metabolites.keys())[0]
            unlockers.append((r.id, met.name, g))

    print(f"  Found {len(unlockers)} single exchanges that unlock growth:")
    for rid, metname, g in sorted(unlockers, key=lambda x: -x[2])[:15]:
        print(f"    {rid} ({metname}): growth={g:.4f}")
    print()
