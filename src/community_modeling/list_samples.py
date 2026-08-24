import sys
import pandas as pd
sys.path.insert(0, "/home/aceren/cycling_dir/coco_paper/bin")
from coco.coco import CoCo

gene_expr = pd.read_pickle("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
cocoGEM_builder = CoCo(gene_expr, default_ub=1000.0)

print(f"Total samples: {len(cocoGEM_builder.samples)}")
print()
for s in sorted(cocoGEM_builder.samples):
    print(" ", s)

print()
ft_samples = [s for s in cocoGEM_builder.samples if s.lower().startswith("cft") or "ft" in s.lower()]
print(f"FT-condition-looking samples: {len(ft_samples)}")
for s in sorted(ft_samples):
    print(" ", s)
