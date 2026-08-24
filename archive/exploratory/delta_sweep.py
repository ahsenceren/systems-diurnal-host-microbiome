"""
delta_sweep.py
================

The diagnostic showed CoCo.build() with delta=1.0 crushes Roseburia's
central glycolysis reactions (GAPD, TPI, etc.) from unbounded (1000) down
to ~1.0-1.5 mmol/gDW/h, killing community growth. This matches exactly why
the ORIGINAL CoCo paper's own notebook explores delta over a grid
(np.arange(1., 11.) in their exploration phase) instead of fixing it at 1.
We skipped that exploration -- doing it now, properly.

Uses the `with com as coco_com:` pattern (same as the real simulations.ipynb)
so each delta trial starts from the SAME diet-constrained base bounds,
not accumulating mutations from previous trials.
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

sample_name = "cFT13a_S48"
if sample_name not in cocoGEM_builder.samples:
    sample_name = cocoGEM_builder.samples[0]

print(f"Sweeping delta at gamma=1.0, alpha(fraction)=0.9, sample={sample_name}")
print("=" * 70)

deltas = [1, 2, 3, 5, 7, 10, 15, 20, 30, 50, 75, 100]
results = []
for d in deltas:
    with com as coco_com:
        cocoGEM_builder.build(coco_com, sample_name, delta=float(d), gamma=1.0, meta=True)
        sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=False)
        obj = sol.objective_value
        results.append((d, obj))
        print(f"  delta={d:>4}: objective_value = {obj}")

print()
print("=" * 70)
print("SUMMARY")
print("=" * 70)
growing = [(d, o) for d, o in results if o is not None and o > 1e-6]
if growing:
    print(f"  Growth-permitting deltas: {[d for d, o in growing]}")
    print(f"  Smallest growth-permitting delta: {growing[0][0]} (objective={growing[0][1]:.4f})")
else:
    print("  *** NO delta value in this range permits growth -- need wider range or gamma/alpha adjustment ***")
