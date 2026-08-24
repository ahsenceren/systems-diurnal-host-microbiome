"""
list_samples_grouped.py
========================
Groups all real samples in coco_gene_expr.pkl by condition (NA/FA/FT/NT)
and Zeitgeber timepoint, to identify every real replicate available --
not just the 'a' replicate used in the first pass. This is step 1 toward
averaging the NA/FA/FT butyrate comparison across all real replicates
instead of one.
"""
import sys
import re
import pandas as pd
from collections import defaultdict

sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
from coco.coco import CoCo

gene_expr = pd.read_pickle("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
cocoGEM_builder = CoCo(gene_expr, default_ub=1000.0)

samples = sorted(cocoGEM_builder.samples)
print(f"Total samples: {len(samples)}\n")

# parse pattern like cFT13a_S48 -> condition=FT, zt=13, rep=a
pattern = re.compile(r"^c([A-Z]{2})(\d{2})([a-z])_S\d+$")

grouped = defaultdict(list)
unparsed = []
for s in samples:
    m = pattern.match(s)
    if m:
        cond, zt, rep = m.group(1), int(m.group(2)), m.group(3)
        grouped[(cond, zt)].append((rep, s))
    else:
        unparsed.append(s)

for cond in ["NA", "FA", "FT", "NT"]:
    print(f"--- {cond} ---")
    zts = sorted(set(zt for (c, zt) in grouped if c == cond))
    for zt in zts:
        reps = sorted(grouped[(cond, zt)])
        rep_str = ", ".join(f"{r}={name}" for r, name in reps)
        print(f"  ZT{zt:>2}: {rep_str}")
    print()

if unparsed:
    print("--- Unparsed sample names (check pattern) ---")
    for s in unparsed:
        print(" ", s)
