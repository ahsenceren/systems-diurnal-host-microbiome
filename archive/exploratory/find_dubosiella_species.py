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

print(f"Loaded {len(lineage_map)} genome->species mappings")

df = pd.read_csv("data/genome_metaG/genome.tsv", sep="\t")
feature_col = df.columns[0]
sample_cols = df.columns[1:]

df["genome_acc"] = df[feature_col].str.replace(r"_\d+$", "", regex=True)
df["species"] = df["genome_acc"].map(lineage_map)

df["row_total"] = df[sample_cols].sum(axis=1)
species_totals = df.groupby("species")["row_total"].sum().sort_values(ascending=False)

print()
print("=== Species containing target keywords ===")
for kw in ["Roseburia intestinalis", "Faecalibacterium", "CAG-873", "Dubosiella", "newyorkensis"]:
    matches = species_totals[species_totals.index.str.contains(kw, case=False, na=False)]
    print(f"--- matching '{kw}' ---")
    print(matches.to_string())
    print()
