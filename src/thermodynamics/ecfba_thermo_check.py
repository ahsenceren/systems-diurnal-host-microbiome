"""
ecfba_thermo_check.py
======================
Real thermodynamic feasibility check for the butyryl-CoA:acetate
CoA-transferase (BCoAT) reaction -- the enzyme confirmed real via
UniProt G2SYC0 (Roseburia hominis) -- using equilibrator-api
(component-contribution method, Noor lab, real published thermodynamic
database, https://equilibrator.weizmann.ac.il). No organism-specific
data needed: standard Gibbs free energies are metabolite-identity-based,
not organism-specific.

Reaction (Rhea:30071, confirmed from UniProt G2SYC0 entry):
  butanoyl-CoA + acetate <=> butanoate + acetyl-CoA
  (KEGG: C00332 + C00033 <=> C00246 + C00024)

This is the reaction our community model uses to produce net butyrate
export. We test: is this reaction thermodynamically favorable in the
direction our FBA optimum requires (net butyrate production), at
gut-like pH and 37C.

Run in WSL (needs internet access to download the equilibrator compound
cache from Zenodo the first time -- ~200MB, one-time download):
  pip install equilibrator-api
  python3 ecfba_thermo_check.py
"""
from equilibrator_api import ComponentContribution, Q_

print("Loading component-contribution model (real thermodynamic database)...")
print("(first run downloads ~200MB compound cache from Zenodo -- be patient)")
cc = ComponentContribution()

# gut/colonic physiological condition: pH ~6.5-7.0, 37C
cc.p_h = Q_(6.8)
cc.ionic_strength = Q_("0.25M")
cc.p_mg = Q_(3.0)
cc.temperature = Q_("310.15K")  # 37C

rxn_str = "kegg:C00332 + kegg:C00033 = kegg:C00246 + kegg:C00024"
rxn = cc.parse_reaction_formula(rxn_str)

print("\nReaction: butanoyl-CoA + acetate <=> butanoate + acetyl-CoA")
print(f"Balanced: {rxn.is_balanced()}")

dg0 = cc.standard_dg_prime(rxn)
print("\nStandard Gibbs free energy (delta-G'-standard, pH=6.8, 37C, I=0.25M):")
print(f"  {dg0}")

dgm = cc.physiological_dg_prime(rxn)
print("\nPhysiological delta-G'm (1 mM reference concentration for all reactants,")
print("standard tFBA 'reversibility index' convention, Jankowski et al. 2008):")
print(f"  {dgm}")

print("\n--- Interpretation ---")
print("If delta-G' is negative (favorable) in the forward direction as written")
print("(butanoyl-CoA + acetate -> butanoate + acetyl-CoA), this is thermodynamic")
print("evidence that our FBA-predicted net butyrate export direction is consistent")
print("with the real Gibbs free energy of the reaction, not just LP-feasible.")
print("If it is near zero or positive, the reaction is close to equilibrium or")
print("reverse-favored at standard conditions, and real gut metabolite")
print("concentrations (acetate/butanoate typically mM-range and abundant,")
print("CoA-thioesters typically sub-mM) are what would need to pull it forward --")
print("worth flagging honestly either way.")
