"""
run_bsh_community_v2.py
=========================
Community model using the compartment-fixed, collision-free GEMs
(*_bsh2.xml). Prints the raw list of bile-acid-related community exchange
IDs immediately after Community() construction, so we can CONFIRM (not
assume) that EX_tcholate_m / EX_tmurichol_m etc. now exist before trying to
set the medium.
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

BILE_KEYWORDS = ["cholate", "murichol", "taurine"]
bile_related = [e for e in exchange_ids if any(k in e.lower() for k in BILE_KEYWORDS)]
print("DIAGNOSTIC -- all community exchange IDs matching bile-acid keywords:")
print(bile_related)
print()

CORE_UNIVERSAL_BIGG = {
    "EX_h2o_e": 100.0, "EX_h_e": 100.0, "EX_pi_e": 10.0, "EX_co2_e": 50.0,
    "EX_cl_e": 10.0, "EX_na1_e": 10.0, "EX_mg2_e": 5.0, "EX_mn2_e": 1.0,
    "EX_zn2_e": 1.0, "EX_cu2_e": 1.0, "EX_fe3_e": 1.0, "EX_so4_e": 5.0,
    "EX_nh4_e": 10.0, "EX_cobalt2_e": 0.5, "EX_ni2_e": 0.5,
}

diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}

BILE_ACID_HOST_BIGG = {
    "EX_tcholate_e": 0.2,
    "EX_tmurichol_e": 0.2,
}

TRACE_BOUND = 2.0
ANALYTE_EXCLUDE_SUBSTR = ["but", "cholate", "murichol", "taurine"]

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

matched_bile = 0
for bigg_id, val in BILE_ACID_HOST_BIGG.items():
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in exchange_ids:
        full_medium[cand] = val
        matched_bile += 1
    else:
        print(f"WARNING: {cand} not found in community exchange list -- bile acid medium entry not applied")

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

print(f"Core matched: {matched_core}/{len(CORE_UNIVERSAL_BIGG)}")
print(f"Diet matched: {matched_diet}/{len(diet_medium_bigg)}")
print(f"Bile acid (host) matched: {matched_bile}/{len(BILE_ACID_HOST_BIGG)}")
print(f"Trace-bound exchanges: {trace_count}")
print(f"Excluded from trace (analyte pathway members): {len(excluded_as_analyte)} -> {excluded_as_analyte}")
print(f"Total medium size: {len(full_medium)}")
com.medium = full_medium

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]
    print(f"NOTE: requested sample not found, using {sample_name} instead")

print(f"Sample: {sample_name}  |  delta=5, gamma=1.0, alpha(fraction)=0.9")
print("=" * 70)

with com as coco_com:
    cocoGEM_builder.build(coco_com, sample_name, delta=5.0, gamma=1.0, meta=True)
    sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=True)

    print(f"status: {sol.status}")
    print(f"objective_value: {sol.objective_value}")
    print()
    print("Per-taxon growth rate:")
    print(sol.members[["abundance", "growth_rate"]])
    print()

    print("=" * 70)
    print("BSH / BILE-ACID-RELATED FLUXES")
    print("=" * 70)
    bile_cols = [c for c in sol.fluxes.columns
                 if any(s in c.lower() for s in ["bsh", "cholate", "murichol", "taurine"])]
    print(f"Columns found: {bile_cols}")
    print()
    print(sol.fluxes[bile_cols].to_string())

    print()
    print("=" * 70)
    print("Community-level bile-acid exchange (net secretion/uptake, EX_*_m)")
    print("=" * 70)
    community_bile_ex = [c for c in sol.fluxes.columns
                          if c.startswith("EX_") and c.endswith("_m")
                          and any(s in c.lower() for s in ["cholate", "murichol", "taurine"])]
    for c in community_bile_ex:
        print(f"  {c}:")
        print(sol.fluxes[c].to_string())
