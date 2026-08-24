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

biomass_rxns = [r.id for r in com.reactions if "biomass" in r.id.lower()]
print("All reaction IDs containing 'biomass':")
for r in biomass_rxns:
    print(" ", r)

print()
print("Community objective expression (post-construction, pre-solve):")
print(com.objective.expression)

print()
print("com.variables sample (first 5, to see naming pattern):")
for v in list(com.variables)[:5]:
    print(" ", v.name)
