"""
butyrate_anova.py
==================
Formal statistical test of the condition (NA/FA/FT) and Zeitgeber-timepoint
effects on community butyrate export capacity, using the real 42-replicate
dataset produced by butyrate.py (butyrate_all_replicates_results.json).

Replaces the earlier eyeballed "means look similar / different" comparison
with a real two-way ANOVA: does condition matter (magnitude), does ZT
matter (diurnal effect pooled across conditions), and critically, does the
condition-by-ZT interaction matter (i.e., is the diurnal SHAPE different
between conditions -- the original, unreplicated claim)?
"""
import json
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols

with open("/home/aceren/diurnal_host_microbiome/butyrate_all_replicates_results.json") as f:
    data = json.load(f)

df = pd.DataFrame(data)
df = df.dropna(subset=["max_butyrate_mmol_gDW_h"])
df["condition"] = df["condition"].astype("category")
df["ZT"] = df["ZT"].astype("category")

print(f"N = {len(df)}")
print(df.groupby(["condition", "ZT"], observed=True).size())
print()

model = ols("max_butyrate_mmol_gDW_h ~ C(condition) + C(ZT) + C(condition):C(ZT)", data=df).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
print(anova_table)
print()
print("Model R-squared:", model.rsquared)

print()
print("=" * 70)
print("Interpretation:")
print(f"  condition main effect p = {anova_table.loc['C(condition)', 'PR(>F)']:.4f}"
      f" ({'significant' if anova_table.loc['C(condition)', 'PR(>F)'] < 0.05 else 'not significant'})")
print(f"  ZT main effect p        = {anova_table.loc['C(ZT)', 'PR(>F)']:.4f}"
      f" ({'significant' if anova_table.loc['C(ZT)', 'PR(>F)'] < 0.05 else 'not significant'})")
print(f"  condition x ZT p        = {anova_table.loc['C(condition):C(ZT)', 'PR(>F)']:.4f}"
      f" ({'significant' if anova_table.loc['C(condition):C(ZT)', 'PR(>F)'] < 0.05 else 'not significant'})")
