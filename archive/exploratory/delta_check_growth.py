"""
delta_check_growth.py
=======================

The delta sweep's objective_value oscillated (negative at some deltas,
positive at others), which is suspicious for something that should
represent community growth. This checks the ACTUAL per-taxon growth_rate
and solver status at each delta, not just the (possibly secondary-QP-stage)
objective_value.
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

for d in [1, 5, 10, 30, 50]:
    print("=" * 70)
    print(f"delta = {d}")
    print("=" * 70)
    with com as coco_com:
        cocoGEM_builder.build(coco_com, sample_name, delta=float(d), gamma=1.0, meta=True)
        sol = coco_com.cooperative_tradeoff(fraction=0.9, fluxes=False)
        print(f"  status: {sol.status}")
        print(f"  objective_value: {sol.objective_value}")
        print(sol.members[["abundance", "growth_rate"]])
    print()
