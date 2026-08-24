import pandas as pd

lineage_map = {}
with open("../data/wol2_taxonomy/lineages.txt") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        acc = parts[0].strip()
        lineage = parts[1]
        species = None
        for token in lineage.split(";"):
            token = token.strip()
            if token.startswith("s__"):
                species = token[3:]
                break
        if species:
            lineage_map[acc] = species

df = pd.read_csv("data/genome_metaG/genome.tsv", sep="\t")
feature_col = df.columns[0]
sample_cols = df.columns[1:]
df["genome_acc"] = df[feature_col].str.replace(r"_\d+$", "", regex=True)
df["species"] = df["genome_acc"].map(lineage_map)
df["row_total"] = df[sample_cols].sum(axis=1)

for sp in ["Roseburia intestinalis", "Faecalibacterium prausnitzii_C", "CAG-873 sp001689415", "Dubosiella newyorkensis"]:
    sub = df[df["species"] == sp]
    n_genomes = sub["genome_acc"].nunique()
    n_rows = len(sub)
    print(f"{sp:35s}: {n_genomes} distinct genome accessions, {n_rows} feature rows, total={sub['row_total'].sum():,.0f}")
    print(f"    per-genome breakdown: {sub.groupby('genome_acc')['row_total'].sum().to_dict()}")
