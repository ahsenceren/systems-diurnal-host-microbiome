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

targets = ["Roseburia intestinalis", "Faecalibacterium prausnitzii_C", "CAG-873 sp001689415", "Dubosiella newyorkensis"]
results = {}
for sp in targets:
    sub = df[df["species"] == sp]
    n_features = len(sub)
    total = sub["row_total"].sum()
    per_feature_mean = total / n_features
    results[sp] = per_feature_mean
    print(f"{sp:35s}: total={total:>8,.0f}  n_features={n_features:>5d}  per-feature mean={per_feature_mean:>8.3f}")

print()
grand_total = sum(results.values())
print("=== Normalized abundance fractions (per-feature mean based) ===")
for sp, val in results.items():
    print(f"{sp:35s}: {val/grand_total:.4f}")
