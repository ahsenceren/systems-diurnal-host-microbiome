import pandas as pd
import numpy as np
import cobra
from pathlib import Path

GENOME_TPM_PATH = "/home/aceren/diurnal_host_microbiome/data/genome_metaT/genome-TPM.tsv"

GEM_PATHS = {
    "Muribaculaceae": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
    "Roseburia": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml",
    "Faecalibacterium": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
}

GENOME_PREFIX = {
    "Muribaculaceae": "G001689415_",
    "Roseburia": "G000156535_",
    "Faecalibacterium": "G000162015_",
}

FLOOR_TPM = 0.01

print("STEP 1: Load genome-TPM.tsv")
tpm_full = pd.read_csv(GENOME_TPM_PATH, sep="\t")
tpm_full = tpm_full.rename(columns={"#FeatureID": "FeatureID"})
print(f"  Full table shape: {tpm_full.shape}")

print("STEP 2: Filter to our 3 genomes + build MultiIndex")
records = []
for taxon, prefix in GENOME_PREFIX.items():
    sub = tpm_full[tpm_full["FeatureID"].str.startswith(prefix)].copy()
    sub["taxon"] = taxon
    print(f"  {taxon}: {len(sub)} genes found (detected)")
    records.append(sub)

detected = pd.concat(records, ignore_index=True)
detected = detected.set_index(["taxon", "FeatureID"])
detected.index.names = ["taxon", "gene_id"]

sample_cols = [c for c in tpm_full.columns if c != "FeatureID"]

print("STEP 3: Add missing genes from each GEM at floor TPM")
all_rows = [detected]
for taxon, gem_path in GEM_PATHS.items():
    model = cobra.io.read_sbml_model(gem_path)
    gem_genes = set(g.id for g in model.genes)
    detected_genes = set(detected.loc[taxon].index) if taxon in detected.index.get_level_values(0) else set()
    missing_genes = gem_genes - detected_genes
    print(f"  {taxon}: {len(gem_genes)} in GEM, {len(detected_genes)} detected, {len(missing_genes)} missing")
    if missing_genes:
        floor_df = pd.DataFrame(
            FLOOR_TPM,
            index=pd.MultiIndex.from_product([[taxon], list(missing_genes)], names=["taxon", "gene_id"]),
            columns=sample_cols,
        )
        all_rows.append(floor_df)

gene_expr = pd.concat(all_rows)
gene_expr = gene_expr[sample_cols]

print(f"Final gene_expr shape: {gene_expr.shape}")
print(f"Taxa: {gene_expr.index.get_level_values('taxon').unique().tolist()}")

out_path = Path("/home/aceren/diurnal_host_microbiome/data/coco_gene_expr.pkl")
gene_expr.to_pickle(out_path)
print(f"Saved: {out_path}")

print("Sample of data:")
print(gene_expr.iloc[:5, :5])
