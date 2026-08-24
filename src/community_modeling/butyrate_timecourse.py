"""
butyrate_timecourse.py
========================
Reruns the validated production-envelope butyrate-maximization approach
across ALL 6 real Zeitgeber timepoints available for the FT condition
(ZT01, ZT05, ZT09, ZT13, ZT17, ZT21; 'a' replicate), confirmed present in
coco_gene_expr.pkl. This replaces our single static butyrate value with a
real, multi-point diurnal profile -- no interpolation/fabrication between
real optimizations, each point is its own real solved LP.
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

CORE_UNIVERSAL_BIGG = {
    "EX_h2o_e": 100.0, "EX_h_e": 100.0, "EX_pi_e": 10.0, "EX_co2_e": 50.0,
    "EX_cl_e": 10.0, "EX_na1_e": 10.0, "EX_mg2_e": 5.0, "EX_mn2_e": 1.0,
    "EX_zn2_e": 1.0, "EX_cu2_e": 1.0, "EX_fe3_e": 1.0, "EX_so4_e": 5.0,
    "EX_nh4_e": 10.0, "EX_cobalt2_e": 0.5, "EX_ni2_e": 0.5,
}
BILE_ACID_HOST_BIGG = {"EX_tcholate_e": 0.2, "EX_tmurichol_e": 0.2}
TRACE_BOUND = 2.0
ANALYTE_EXCLUDE_SUBSTR = ["but", "cholate", "murichol", "taurine"]

ZT_POINTS = [1, 5, 9, 13, 17, 21]
SAMPLE_MAP = {1: "cFT01a_S51", 5: "cFT05a_S52", 9: "cFT09a_S53",
              13: "cFT13a_S48", 17: "cFT17a_S49", 21: "cFT21a_S50"}

results = []

for zt in ZT_POINTS:
    sample_name = SAMPLE_MAP[zt]
    if sample_name not in cocoGEM_builder.samples:
        print(f"ZT{zt}: sample {sample_name} NOT FOUND -- skipping")
        continue

    com = mc.Community(taxonomy)
    exchange_ids = [r.id for r in com.exchanges]

    diet_medium_bigg = build_medium("FT", zt)
    diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}

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

    try:
        with com as coco_com:
            cocoGEM_builder.build(coco_com, sample_name, delta=5.0, gamma=1.0, meta=True)
            sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=True)
            growth_rates = sol.members["growth_rate"].dropna().to_dict()

            fixed_count = 0
            for member_id, g in growth_rates.items():
                candidates = [r for r in coco_com.reactions
                              if r.id.endswith(f"__{member_id}") and "growth" in r.id.lower()]
                if candidates:
                    candidates[0].lower_bound = max(0, g * 0.9)
                    fixed_count += 1

            if fixed_count < 3:
                print(f"ZT{zt}: could not fix all 3 growth reactions -- skipping")
                continue

            but_rxn = coco_com.reactions.get_by_id("EX_but_m")
            coco_com.objective = but_rxn
            coco_com.objective_direction = "max"
            sol2 = coco_com.optimize(fluxes=True)

            max_but = sol2.fluxes.loc["medium", "EX_but_m"] if "medium" in sol2.fluxes.index else None
            print(f"ZT{zt} ({sample_name}): status={sol2.status}, max_butyrate_export={max_but}, growth_status={sol.status}")
            results.append({"ZT": zt, "sample": sample_name, "status": sol2.status,
                             "max_butyrate_mmol_gDW_h": max_but,
                             "growth_rates": growth_rates})
    except Exception as e:
        print(f"ZT{zt}: FAILED -- {type(e).__name__}: {e}")

print()
print("=" * 60)
print("REAL DIURNAL BUTYRATE EXPORT-CAPACITY PROFILE (FT condition)")
print("=" * 60)
for r in results:
    print(f"  ZT{r['ZT']:>2}: {r['max_butyrate_mmol_gDW_h']:.4f} mmol/gDW/h  (status={r['status']})")

import json
with open("/home/aceren/diurnal_host_microbiome/butyrate_timecourse_results.json", "w") as f:
    json.dump([{"ZT": r["ZT"], "sample": r["sample"], "status": r["status"],
                "max_butyrate_mmol_gDW_h": r["max_butyrate_mmol_gDW_h"]} for r in results], f, indent=2)
print("\nSaved: /home/aceren/diurnal_host_microbiome/butyrate_timecourse_results.json")
