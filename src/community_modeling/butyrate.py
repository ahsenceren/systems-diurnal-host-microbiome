"""
butyrate_all_conditions_replicates.py
=======================================
Extends butyrate_all_conditions.py from single-replicate ('a' only, 18
optimizations) to ALL real replicates: NA=3 reps, FA=2 reps, FT=2 reps,
each x 6 ZT = 42 real, independent LP optimizations total. Reports
mean +/- std per (condition, ZT) instead of a single descriptive value,
turning the earlier "descriptive, not statistically tested" pattern into
one with real inter-replicate variance attached.
"""
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
sys.path.insert(0, "/home/aceren/diurnal_host_microbiome")

from coco.coco import CoCo
import micom as mc
from diet_to_medium import build_medium

gene_expr = pd.read_pickle("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
cocoGEM_builder = CoCo(gene_expr, default_ub=1000.0)

taxonomy = pd.DataFrame({
    "id": ["Muribaculaceae", "Roseburia", "Faecalibacterium", "Dubosiella"],
    "genus": ["Muribaculaceae", "Roseburia", "Faecalibacterium", "Dubosiella"],
    "species": ["Muribaculaceae sp. CAG-873", "Roseburia intestinalis", "Faecalibacterium duncaniae", "Dubosiella newyorkensis"],
    "file": [
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2_bsh2.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed_bsh2.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
        "/home/aceren/diurnal_host_microbiome/data/gems_wol2/dubosiella_wol2.xml",
    ],
    "abundance": [0.4053, 0.0975, 0.0955, 0.4017],
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


# Full real replicate map (from list_samples_grouped.py output)
SAMPLE_MAP = {
    "NA": {
        1: ["cNA01a_S4", "cNA01b_S10", "cNA01c_S16"],
        5: ["cNA05a_S5", "cNA05b_S11", "cNA05c_S17"],
        9: ["cNA09a_S6", "cNA09b_S12", "cNA09c_S18"],
        13: ["cNA13a_S1", "cNA13b_S7", "cNA13c_S13"],
        17: ["cNA17a_S2", "cNA17b_S8", "cNA17c_S14"],
        21: ["cNA21a_S3", "cNA21b_S9", "cNA21c_S15"],
    },
    "FA": {
        1: ["cFA01a_S39", "cFA01b_S45"],
        5: ["cFA05a_S40", "cFA05b_S46"],
        9: ["cFA09a_S41", "cFA09b_S47"],
        13: ["cFA13a_S36", "cFA13b_S42"],
        17: ["cFA17a_S37", "cFA17b_S43"],
        21: ["cFA21a_S38", "cFA21b_S44"],
    },
    "FT": {
        1: ["cFT01a_S51", "cFT01b_S57"],
        5: ["cFT05a_S52", "cFT05b_S58"],
        9: ["cFT09a_S53", "cFT09b_S59"],
        13: ["cFT13a_S48", "cFT13b_S54"],
        17: ["cFT17a_S49", "cFT17b_S55"],
        21: ["cFT21a_S50", "cFT21b_S56"],
    },
}

def run_one(condition, zt, sample_name):
    if sample_name not in cocoGEM_builder.samples:
        print(f"{condition} ZT{zt} {sample_name}: NOT FOUND -- skipping")
        return None

    com = mc.Community(taxonomy)
    exchange_ids = [r.id for r in com.exchanges]

    try:
        diet_medium_bigg = build_medium(condition, zt)
    except Exception as e:
        print(f"{condition} ZT{zt} {sample_name}: build_medium FAILED -- {type(e).__name__}: {e}")
        return None
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
                print(f"{condition} ZT{zt} {sample_name}: could not fix all 3 growth reactions -- skipping")
                return None

            but_rxn = coco_com.reactions.get_by_id("EX_but_m")
            coco_com.objective = but_rxn
            coco_com.objective_direction = "max"
            sol2 = coco_com.optimize(fluxes=True)

            max_but = sol2.fluxes.loc["medium", "EX_but_m"] if "medium" in sol2.fluxes.index else None
            print(f"{condition} ZT{zt:>2} {sample_name}: status={sol2.status}, max_butyrate={max_but}")
            return max_but
    except Exception as e:
        print(f"{condition} ZT{zt} {sample_name}: FAILED -- {type(e).__name__}: {e}")
        return None

all_results = []
for condition, zt_map in SAMPLE_MAP.items():
    for zt, sample_list in zt_map.items():
        for sample_name in sample_list:
            val = run_one(condition, zt, sample_name)
            all_results.append({"condition": condition, "ZT": zt, "sample": sample_name,
                                 "max_butyrate_mmol_gDW_h": val})

df = pd.DataFrame(all_results)
df.to_json("/home/aceren/diurnal_host_microbiome/butyrate_all_replicates_results.json",
           orient="records", indent=2)

print()
print("=" * 70)
print("SUMMARY: mean +/- std butyrate export capacity, all real replicates")
print("=" * 70)
summary = df.dropna(subset=["max_butyrate_mmol_gDW_h"]).groupby(["condition", "ZT"])["max_butyrate_mmol_gDW_h"].agg(["mean", "std", "count"])
print(summary.to_string())

print()
print("Overall mean across ZT by condition:")
overall = df.dropna(subset=["max_butyrate_mmol_gDW_h"]).groupby("condition")["max_butyrate_mmol_gDW_h"].agg(["mean", "std", "count"])
print(overall.to_string())

print("\nSaved: /home/aceren/diurnal_host_microbiome/butyrate_all_replicates_results.json")
