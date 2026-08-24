"""
final_result.py
=================

FINAL result for the abstract. delta=5 chosen because it was the smallest
value in the sweep giving status=optimal with positive, biologically
plausible growth for all 3 taxa (0.15-0.77/h).

MEDIUM: core (diet-independent universal compounds) + real diet-derived
compounds + a low uniform TRACE bound (-2) on everything else, representing
real but minor background cecal availability (host mucin, biliary/pancreatic
secretions, microbial cross-feeding) not captured by our diet-only mapping.
Empirically confirmed necessary: Muribaculaceae and Faecalibacterium show
zero growth on core+diet alone (even at 10x bound magnitude -- not a scale
issue), but grow once trace background availability is added.

pfba dropped: it added a secondary constraint that became infeasible given
CoCo's asymmetric bound scaling, even though the underlying growth
optimization is feasible. fluxes=True alone returns a valid flux
distribution at the optimal growth solution.
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
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
    ],
    "abundance": [1913/2779, 443/2779, 423/2779],
})
com = mc.Community(taxonomy)
exchange_ids = [r.id for r in com.exchanges]

CORE_UNIVERSAL_BIGG = {
    "EX_h2o_e": 100.0, "EX_h_e": 100.0, "EX_pi_e": 10.0, "EX_co2_e": 50.0,
    "EX_cl_e": 10.0, "EX_na1_e": 10.0, "EX_mg2_e": 5.0, "EX_mn2_e": 1.0,
    "EX_zn2_e": 1.0, "EX_cu2_e": 1.0, "EX_fe3_e": 1.0, "EX_so4_e": 5.0,
    "EX_nh4_e": 10.0, "EX_cobalt2_e": 0.5, "EX_ni2_e": 0.5,
}

diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}

TRACE_BOUND = 2.0

full_medium = {}
matched_core = 0
for bigg_id, val in CORE_UNIVERSAL_BIGG.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
        matched_core += 1
matched_diet = 0
for bigg_id, val in diet_medium_bigg.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
        matched_diet += 1

# BUG FIX: the trace loop was also opening EX_but_e/EX_but_m for IMPORT
# (lower_bound=-2), letting the community "cheat" by importing free
# butyrate from the trace background instead of genuinely synthesizing it
# -- confirmed by EX_but_m landing at exactly -1.999996, i.e. pinned at the
# trace bound. We are measuring butyrate as an OUTPUT, so its own exchange
# (and related analytes) must be excluded from the trace-import set; left
# out of full_medium, they default to the standard product boundary
# condition (import closed, export unrestricted).
ANALYTE_EXCLUDE_SUBSTR = ["but"]

trace_count = 0
excluded_as_analyte = []
for ex_id in exchange_ids:
    if ex_id in full_medium:
        continue
    if any(s in ex_id.lower() for s in ANALYTE_EXCLUDE_SUBSTR):
        excluded_as_analyte.append(ex_id)
        continue
    full_medium[ex_id] = TRACE_BOUND
    trace_count += 1

print(f"Excluded from trace-import (analyte of interest, export-only): {excluded_as_analyte}")

print(f"Core universal compounds matched: {matched_core}/{len(CORE_UNIVERSAL_BIGG)}")
print(f"Diet-derived compounds matched: {matched_diet}/{len(diet_medium_bigg)}")
print(f"Trace-bound (background) exchanges: {trace_count}")
print(f"Total medium size: {len(full_medium)}")
com.medium = full_medium

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]
    print(f"NOTE: requested sample not found, using {sample_name} instead")

print(f"Sample: {sample_name}  |  delta=5, gamma=1.0, alpha(fraction)=0.9")
print("=" * 70)

print("Sanity check -- growth on this medium BEFORE CoCo:")
sol_pre = com.cooperative_tradeoff(fraction=0.9, fluxes=False)
print(f"  status={sol_pre.status}  objective_value={sol_pre.objective_value}")
print(sol_pre.members[["abundance", "growth_rate"]])
print()

with com as coco_com:
    cocoGEM_builder.build(coco_com, sample_name, delta=5.0, gamma=1.0, meta=True)
    sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=True, pfba=True)

    print(f"status: {sol.status}")
    print(f"objective_value: {sol.objective_value}")
    print()
    print("Per-taxon growth rate:")
    print(sol.members[["abundance", "growth_rate"]])
    print()

    print("=" * 70)
    print("BUTYRATE-RELATED FLUXES (all columns containing 'but')")
    print("=" * 70)
    but_cols = [c for c in sol.fluxes.columns if "but" in c.lower()]
    print(sol.fluxes[but_cols].to_string())

    print()
    print("=" * 70)
    print("Community-level butyrate exchange (net secretion to medium)")
    print("=" * 70)
    community_but_ex = [c for c in sol.fluxes.columns if c.startswith("EX_but")]
    for c in community_but_ex:
        print(f"  {c}:")
        print(sol.fluxes[c].to_string())
