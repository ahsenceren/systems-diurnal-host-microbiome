"""
diag_coco_bounds.py
=====================

Diagnostic: v4 confirmed the community grows fine on default medium
(23.76) and on the real diet medium alone (2.51), but growth drops to
EXACTLY 0.0 after cocoGEM_builder.build() is applied. This script does
NOT guess at delta/gamma -- it directly compares reaction bounds before
and after build(), so we can see exactly which reactions CoCo is
crushing and whether that's biologically sane or a bug.
"""

import sys
import time
import pandas as pd
import numpy as np

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
sys.path.insert(0, "/home/aceren/diurnal_host_microbiome")

from coco.coco import CoCo
import micom as mc
from diet_to_medium import build_medium

print("Loading gene expression + building community (same as v4)...")
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
default_medium = dict(com.medium)

diet_medium_bigg = build_medium("FT", 13)
diet_medium_bigg = {k: v for k, v in diet_medium_bigg.items() if v > 1e-9}
full_medium = dict(default_medium)
for bigg_id, val in diet_medium_bigg.items():
    candidates = []
    if bigg_id.endswith("_e"):
        candidates.append(bigg_id[:-2] + "_m")
    candidates.append(bigg_id)
    for cand in candidates:
        if cand in exchange_ids:
            full_medium[cand] = val
            break
com.medium = full_medium

print("Growth check before CoCo (should match v4's 2.51):")
sol0 = com.cooperative_tradeoff(fraction=0.9, fluxes=False)
print(f"  objective_value = {sol0.objective_value}")

before = {r.id: (r.lower_bound, r.upper_bound) for r in com.reactions}

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]

print(f"\nApplying CoCo.build() for sample {sample_name}...")
cocoGEM_builder.build(com, sample_name, delta=1.0, gamma=1.0, meta=True)

print("\nGrowth check AFTER CoCo (this is the one that comes out 0.0):")
sol1 = com.cooperative_tradeoff(fraction=0.9, fluxes=False)
print(f"  objective_value = {sol1.objective_value}")
print(f"  status: {getattr(sol1, 'status', 'n/a')}")

after = {r.id: (r.lower_bound, r.upper_bound) for r in com.reactions}

print("\n" + "=" * 70)
print("DIAGNOSTIC 1: Any bounds that became inverted (lower > upper)?")
print("=" * 70)
inverted = [(rid, lb, ub) for rid, (lb, ub) in after.items() if lb > ub]
print(f"  Found {len(inverted)} inverted bound pairs")
for rid, lb, ub in inverted[:20]:
    print(f"    {rid}: lb={lb:.6g} > ub={ub:.6g}  (before: {before[rid]})")

print("\n" + "=" * 70)
print("DIAGNOSTIC 2: Biomass reactions -- bounds before/after")
print("=" * 70)
biomass_rxns = [r.id for r in com.reactions if "biomass" in r.id.lower() or "BIOMASS" in r.id]
for rid in biomass_rxns:
    print(f"  {rid}: before={before.get(rid)} after={after.get(rid)}")

print("\n" + "=" * 70)
print("DIAGNOSTIC 3: Diet-medium exchange reactions -- bounds before/after")
print("=" * 70)
for bigg_id in list(diet_medium_bigg.keys())[:10]:
    cand = bigg_id[:-2] + "_m" if bigg_id.endswith("_e") else bigg_id
    if cand in before:
        print(f"  {cand}: before={before[cand]} after={after[cand]}")

print("\n" + "=" * 70)
print("DIAGNOSTIC 4: Top 20 reactions with the LARGEST upper_bound shrinkage")
print("(as a ratio after/before, excluding already-zero reactions)")
print("=" * 70)
shrinkage = []
for rid in before:
    lb0, ub0 = before[rid]
    lb1, ub1 = after[rid]
    if ub0 > 1e-9:
        ratio = ub1 / ub0
        shrinkage.append((ratio, rid, ub0, ub1, lb0, lb1))
shrinkage.sort(key=lambda x: x[0])
for ratio, rid, ub0, ub1, lb0, lb1 in shrinkage[:20]:
    print(f"  {rid}: ub {ub0:.6g} -> {ub1:.6g} (ratio {ratio:.4g})  lb {lb0:.6g} -> {lb1:.6g}")

print("\n" + "=" * 70)
print("DIAGNOSTIC 5: How many reactions have ub<=0 AND lb>=0 after build()")
print("(i.e. completely closed reactions)")
print("=" * 70)
closed = [rid for rid, (lb, ub) in after.items() if ub <= 1e-9 and lb >= -1e-9]
closed_before = [rid for rid in closed if not (before[rid][1] <= 1e-9 and before[rid][0] >= -1e-9)]
print(f"  Total closed after build(): {len(closed)}")
print(f"  NEWLY closed by build() (were open before): {len(closed_before)}")
print(f"  Examples of newly-closed: {closed_before[:20]}")
