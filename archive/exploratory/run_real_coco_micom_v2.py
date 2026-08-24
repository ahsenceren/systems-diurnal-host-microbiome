"""
run_real_coco_micom_v2.py
===========================

Corrected version following the VALIDATED pattern from coco_paper's
bin/simulations.ipynb (Zampieri et al. 2023):

  com = mc.Community(species_df, solver=...)
  com = <apply medium>
  with com as coco_com:
      cocoGEM_builder.build(coco_com, sample, delta, gamma, meta=True)
      sol = coco_com.cooperative_tradeoff(fraction=alpha, fluxes=True)

Differences from their anaerobic-digestion setup:
  - No CPLEX (not licensed) -- solver=None resolved to optlang hybrid interface
  - No set_biochemical_constraints / set_MAG_constraints (AD-specific)
  - Medium: our own diet-derived medium, matched to the REAL community
    exchange naming (trailing "_e" -> "_m", confirmed empirically)
"""

import sys
import time
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")

print("=" * 70)
print("STEP 1: Imports")
print("=" * 70)
from coco.coco import CoCo
import micom as mc
print(f"  micom version: {mc.__version__}")

print()
print("=" * 70)
print("STEP 2: Load real gene expression data")
print("=" * 70)
gene_expr = pd.read_pickle("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
print(f"  Shape: {gene_expr.shape}")

print()
print("=" * 70)
print("STEP 3: Initialize CoCo")
print("=" * 70)
t0 = time.time()
cocoGEM_builder = CoCo(gene_expr, default_ub=1000.0)
print(f"  Done in {time.time()-t0:.2f}s. Taxa: {cocoGEM_builder.taxa}")

print()
print("=" * 70)
print("STEP 4: Build MICOM Community")
print("=" * 70)

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

t0 = time.time()
com = mc.Community(taxonomy)
print(f"  Community built in {time.time()-t0:.2f}s")
print(f"  Solver in use: {com.solver.interface}")
print(f"  Total reactions: {len(com.reactions)}")
print(f"  Total metabolites: {len(com.metabolites)}")

exchange_ids = [r.id for r in com.exchanges]
print(f"  Total exchange reactions: {len(exchange_ids)}")

print()
print("=" * 70)
print("STEP 5: Apply diet-informed medium")
print("=" * 70)

diet_medium_bigg = {
    "EX_glc__D_e": 10.0,
    "EX_cellb_e": 8.0,
    "EX_ala__L_e": 2.0, "EX_arg__L_e": 1.0, "EX_asn__L_e": 0.6, "EX_asp__L_e": 1.5,
    "EX_cys__L_e": 0.5, "EX_gln__L_e": 1.0, "EX_glu__L_e": 2.0, "EX_gly_e": 1.5,
    "EX_his__L_e": 0.3, "EX_ile__L_e": 0.7, "EX_leu__L_e": 1.2, "EX_lys__L_e": 0.9,
    "EX_met__L_e": 0.4, "EX_phe__L_e": 0.8, "EX_pro__L_e": 0.5, "EX_ser__L_e": 0.8,
    "EX_thr__L_e": 0.7, "EX_trp__L_e": 0.2, "EX_tyr__L_e": 0.6, "EX_val__L_e": 0.9,
    "EX_hdca_e": 3.0, "EX_ca2_e": 1.0, "EX_pi_e": 0.5, "EX_k_e": 2.0, "EX_fe2_e": 0.1,
}

medium_dict = {}
matched = 0
unmatched = []
for bigg_id, val in diet_medium_bigg.items():
    candidates = []
    if bigg_id.endswith("_e"):
        candidates.append(bigg_id[:-2] + "_m")
    candidates.append(bigg_id)
    found = False
    for cand in candidates:
        if cand in exchange_ids:
            medium_dict[cand] = val
            matched += 1
            found = True
            break
    if not found:
        unmatched.append(bigg_id)

print(f"  Matched {matched}/{len(diet_medium_bigg)} diet components to real exchange IDs")
if unmatched:
    print(f"  Unmatched: {unmatched}")
if matched == 0:
    print("  *** WARNING: no matches -- CANNOT PROCEED with this medium ***")
    sys.exit(1)
com.medium = medium_dict
print(f"  Medium applied ({len(medium_dict)} components)")

print()
print("=" * 70)
print("STEP 6: Apply CoCo bound scaling + run cooperative_tradeoff (REAL)")
print("=" * 70)

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]
    print(f"  Requested sample not found, using: {sample_name}")
else:
    print(f"  Using sample: {sample_name}")

ALPHA = 0.9

with com as coco_com:
    t0 = time.time()
    cocoGEM_builder.build(coco_com, sample_name, delta=1.0, gamma=1.0, meta=True)
    t1 = time.time()
    print(f"  CoCo.build() applied in {t1-t0:.2f}s")

    print(f"  Running cooperative_tradeoff(fraction={ALPHA})...")
    t0 = time.time()
    sol = coco_com.cooperative_tradeoff(fraction=ALPHA, fluxes=True, pfba=True)
    elapsed = time.time() - t0

    print()
    print(f"  *** ELAPSED TIME FOR 1 REAL SIMULATION: {elapsed:.2f} seconds ***")
    print()
    print(f"  Objective value: {sol.objective_value}")
    print(f"  Growth rates:")
    print(sol.members[["abundance", "growth_rate"]])

    print()
    print("  Butyrate-related fluxes:")
    flux_but_cols = [c for c in sol.fluxes.columns if "but" in c.lower()]
    if flux_but_cols:
        print(sol.fluxes[flux_but_cols])
    else:
        print("  (none found in flux table columns)")

print()
print("=" * 70)
print(f"TIMING SUMMARY: 1 simulation = {elapsed:.2f}s")
print(f"  Full grid (72 sims) serial estimate: {elapsed*72/60:.1f} min")
print(f"  Full grid (72 sims) parallel(20-core) estimate: {elapsed*72/20/60:.1f} min")
print("=" * 70)
