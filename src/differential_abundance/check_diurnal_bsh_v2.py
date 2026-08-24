import re
import pandas as pd

f = "/home/aceren/diurnal_host_microbiome/birdman_results_combined.txt"

rows = []
malformed = 0
malformed_examples = []
with open(f) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) == 7:
            rows.append(parts)
        elif len(parts) == 6:
            m = re.match(r"^(dark|light)(NAvFA|FAvFT)$", parts[4])
            if m:
                phase, comparison = m.groups()
                fixed = parts[:4] + [phase, comparison] + parts[5:]
                rows.append(fixed)
            else:
                malformed += 1
                if len(malformed_examples) < 5:
                    malformed_examples.append(line)
        else:
            malformed += 1
            if len(malformed_examples) < 5:
                malformed_examples.append(line)

print(f"Parsed rows: {len(rows)}")
print(f"Malformed/unparseable rows skipped: {malformed}")
for ex in malformed_examples:
    print(f"  SKIPPED: {ex}")
print()

df = pd.DataFrame(rows, columns=["taxon_ortholog", "ratio", "min", "max", "phase", "comparison", "method"])
for col in ["ratio", "min", "max"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

n_bad_numeric = df["ratio"].isna().sum()
print(f"Rows where ratio failed to parse as numeric even after fix: {n_bad_numeric}")
df = df.dropna(subset=["ratio"])

targets = ["Roseburia", "Duncaniella", "Paramuribaculum", "Muribaculum"]
mask = df["taxon_ortholog"].str.contains("|".join(targets), case=False, na=False)
sub = df[mask].copy()

sub_ft = sub[(sub["comparison"] == "FAvFT") & (sub["method"] == "targetted")]

print(f"Total target-taxa rows (all comparisons/methods): {len(sub)}")
print(f"FAvFT + targetted rows: {len(sub_ft)}")
print()

pivot = sub_ft.pivot_table(index="taxon_ortholog", columns="phase", values="ratio", aggfunc="first")
if "dark" in pivot.columns and "light" in pivot.columns:
    pivot = pivot.dropna(subset=["dark", "light"])
    pivot["diff_dark_minus_light"] = pivot["dark"] - pivot["light"]
    pivot["sign_flip"] = (pivot["dark"] * pivot["light"] < 0)

    pd.set_option("display.max_rows", 200)
    pd.set_option("display.width", 140)
    print(pivot.sort_values("diff_dark_minus_light", ascending=False).to_string())

    print()
    n_flip = pivot["sign_flip"].sum()
    print(f"Rows with SAME taxon/ortholog present in BOTH phases: {len(pivot)}")
    print(f"Rows where dark/light ratio have OPPOSITE sign: {n_flip} / {len(pivot)}")
    print(f"Mean |dark - light| difference: {pivot['diff_dark_minus_light'].abs().mean():.3f}")
    print(f"Median |dark - light| difference: {pivot['diff_dark_minus_light'].abs().median():.3f}")
else:
    print("Could not find both 'dark' and 'light' columns after pivot -- check phase value spelling.")
    print(pivot)
