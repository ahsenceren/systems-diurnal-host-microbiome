"""
max_butyrate.py
=================
Instead of pFBA (which requires a secondary norm-minimization LP that keeps
hitting numerical infeasibility under CoCo's bound scaling), directly
constrain each taxon's growth near its cooperative-tradeoff-optimal value
(90% fraction) and MAXIMIZE net community butyrate export (EX_but_m) as the
objective. This is a simpler, single-stage LP and gives a well-defined,
interpretable number.
"""

import sys
import pandas as pd

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
sys.path.insert(0, "/home/aceren/diurnal_host_microbiome")

from coco.coco import CoCo
import micom as mc
from diet_to_medium import build_medium

gene_expr = pd.read_pickle("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
cocoGEM_builder = CoCo(gene_expr, default_ub=1000.0)

taxonomy = pd.DataFrame({
    "id": ["Muribaculaceae", "Roseburia", "Faecalibacterium"],
    "genus": ["Muribaculaceae", "Roseburia", "Faecalibacterium"],
    "species": ["Muribaculaceae sp. CAG-873", "Roseburia intestinalis", "Faecalibacterium duncaniae"],
    "file": [
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2_bsh2.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed_bsh2.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
    ],
    "abundance": [1913/2779, 443/2779, 423/2779],
})
com = mc.Community(taxonomy)
exchange_ids = [r.id for r in com.exchanges]

growth_rxns = [r.id for r in com.reactions if "growth" in r.id.lower() or "biomass" in r.id.lower()]
print(f"Found {len(growth_rxns)} growth/biomass reactions: {growth_rxns}")
if len(growth_rxns) < 3:
    print("ABORTING -- could not find one growth reaction per taxon.")
    sys.exit(1)

if "EX_but_m" not in exchange_ids:
    print("ABORTING -- EX_but_m not found in community exchange list.")
    sys.exit(1)

CORE_UNIVERSAL_BIGG = {
    "EX_h2o_e": 100.0, "EX_h_e": 100.0, "EX_pi_e": 10.0, "EX_co2_e": 50.0,
    "EX_cl_e": 10.0, "EX_na1_e": 10.0, "EX_mg2_e": 5.0, "EX_mn2_e": 1.0,
    "EX_zn2_e": 1.0, "EX_cu2_e": 1.0, "EX_fe3_e": 1.0, "EX_so4_e": 5.0,
    "EX_nh4_e": 10.0, "EX_cobalt2_e": 0.5, "EX_ni2_e": 0.5,
}
diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}
BILE_ACID_HOST_BIGG = {"EX_tcholate_e": 0.2, "EX_tmurichol_e": 0.2}

TRACE_BOUND = 2.0
ANALYTE_EXCLUDE_SUBSTR = ["but", "cholate", "murichol", "taurine"]

full_medium = {}
for bigg_id, val in CORE_UNIVERSAL_BIGG.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
for bigg_id, val in diet_medium_bigg.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
for bigg_id, val in BILE_ACID_HOST_BIGG.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
for ex_id in exchange_ids:
    if ex_id in full_medium:
        continue
    if any(s in ex_id.lower() for s in ANALYTE_EXCLUDE_SUBSTR):
        continue
    full_medium[ex_id] = TRACE_BOUND
com.medium = full_medium

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]

print(f"Sample: {sample_name}")
print("=" * 70)

with com as coco_com:
    cocoGEM_builder.build(coco_com, sample_name, delta=5.0, gamma=1.0, meta=True)

    sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=True)
    print(f"Normal solve status: {sol.status}, objective: {sol.objective_value}")
    print(sol.members[["abundance", "growth_rate"]])
    print()

    growth_rates = sol.members["growth_rate"].dropna().to_dict()
    fixed_count = 0
    for member_id, g in growth_rates.items():
        candidates = [r for r in coco_com.reactions
                      if r.id.endswith(f"__{member_id}") and
                      ("growth" in r.id.lower() or "biomass" in r.id.lower())]
        if candidates:
            rxn = candidates[0]
            rxn.lower_bound = max(0, g * 0.9)
            print(f"Fixed {rxn.id}: lower_bound = {rxn.lower_bound:.6f} (90% of {g:.6f})")
            fixed_count += 1
        else:
            print(f"WARNING: could not find growth reaction for {member_id}")

    if fixed_count < 3:
        print("ABORTING -- not all 3 members' growth was fixed.")
        sys.exit(1)

    but_rxn = coco_com.reactions.get_by_id("EX_but_m")
    coco_com.objective = but_rxn
    coco_com.objective_direction = "max"

    sol2 = coco_com.optimize(fluxes=True)
    print()
    print("=" * 70)
    print("MAX BUTYRATE SOLVE (growth constrained to >=90% of optimal per taxon)")
    print("=" * 70)
    print(f"status: {sol2.status}")
    print(f"Max feasible net community butyrate export (EX_but_m): {sol2.fluxes.get('EX_but_m', 'NOT FOUND')}")

    but_cols = [c for c in sol2.fluxes.index if "but" in c.lower()]
    print()
    print("All butyrate-related fluxes at this solution:")
    for c in but_cols:
        print(f"  {c}: {sol2.fluxes[c]}")
