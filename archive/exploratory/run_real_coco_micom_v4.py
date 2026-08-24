"""
run_real_coco_micom_v4.py
===========================

Fixes the growth=0.0 bug from v3: com.medium was being REPLACED with only
27 diet-derived exchanges, closing everything else (water, protons, CO2,
core ions) that the community needs to grow but that isn't diet-limited.

v4 approach:
  1. Start from the community's OWN default medium (264 exchanges, all
     open at CarveMe/MICOM's default bounds) as the "always available"
     gut-lumen baseline.
  2. OVERRIDE only the diet-relevant subset with REAL numbers from the
     (now-verified) diet_to_medium.py build_medium() function.
  3. Everything else keeps its default bound.
"""

import sys
import time
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
sys.path.insert(0, "/home/aceren/diurnal_host_microbiome")

print("=" * 70)
print("STEP 1: Imports")
print("=" * 70)
from coco.coco import CoCo
import micom as mc
from diet_to_medium import build_medium
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

exchange_ids = [r.id for r in com.exchanges]
default_medium = dict(com.medium)
print(f"  Default (baseline) medium size: {len(default_medium)}")

print()
print("=" * 70)
print("STEP 4b: SANITY CHECK -- growth on community's OWN default medium?")
print("=" * 70)
sol_default = com.cooperative_tradeoff(fraction=0.9, fluxes=False)
print(f"  Growth on default medium (objective_value): {sol_default.objective_value}")
print(sol_default.members[["abundance", "growth_rate"]])

print()
print("=" * 70)
print("STEP 5: Build REAL diet medium (FT condition, ZT13)")
print("=" * 70)
diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}
print(f"  Real diet-derived medium ({len(diet_medium_bigg)} nonzero components):")
for k, v in diet_medium_bigg.items():
    print(f"    {k}: {v}")

print()
print("=" * 70)
print("STEP 6: Layer diet medium ON TOP of the community's default medium")
print("=" * 70)
full_medium = dict(default_medium)
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
            full_medium[cand] = val
            matched += 1
            found = True
            break
    if not found:
        unmatched.append(bigg_id)

print(f"  Matched {matched}/{len(diet_medium_bigg)} diet components")
if unmatched:
    print(f"  Unmatched (dropped): {unmatched}")
com.medium = full_medium
print(f"  Final medium size: {len(full_medium)}")

print()
print("=" * 70)
print("STEP 6b: SANITY CHECK -- growth on combined medium (no CoCo yet)?")
print("=" * 70)
sol_diet_only = com.cooperative_tradeoff(fraction=0.9, fluxes=False)
print(f"  Growth on combined diet medium (objective_value): {sol_diet_only.objective_value}")
print(sol_diet_only.members[["abundance", "growth_rate"]])

print()
print("=" * 70)
print("STEP 7: Apply CoCo bound scaling + run cooperative_tradeoff (REAL)")
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
