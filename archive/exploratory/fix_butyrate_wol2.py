import cobra
from cobra import Reaction, Metabolite

def add_butyrate_pathway(model_path, output_path):
    model = cobra.io.read_sbml_model(model_path)

    # Detect compartment suffix used in this model (e.g. _c, _e)
    ref_c = None
    for met in model.metabolites:
        if met.id.endswith("_c") and "atp" in met.id.lower():
            ref_c = met.compartment
            break
    ref_e = None
    for met in model.metabolites:
        if met.id == "glc__D_e":
            ref_e = met.compartment
            break

    # Ensure but_c exists
    if "but_c" not in [m.id for m in model.metabolites]:
        but_c = Metabolite("but_c", formula="C4H7O2", name="Butyrate",
                            compartment=ref_c or "C_c")
        model.add_metabolites([but_c])
    else:
        but_c = model.metabolites.get_by_id("but_c")

    # Ensure but_e exists
    if "but_e" not in [m.id for m in model.metabolites]:
        but_e = Metabolite("but_e", formula="C4H7O2", name="Butyrate",
                            compartment=ref_e or "C_e")
        model.add_metabolites([but_e])
    else:
        but_e = model.metabolites.get_by_id("but_e")

    # BUTCT2: butyryl-CoA + acetate <-> butyrate + acetyl-CoA
    # Requires btcoa_c and ac_c and accoa_c to exist in the model
    needed = ["btcoa_c", "ac_c", "accoa_c"]
    missing = [m for m in needed if m not in [x.id for x in model.metabolites]]
    if missing:
        print(f"  WARNING: missing precursor metabolites for BUTCT2: {missing}")
        print(f"  Skipping synthesis reaction -- model lacks upstream butyryl-CoA pathway")
    else:
        if "BUTCT2" not in [r.id for r in model.reactions]:
            butct2 = Reaction("BUTCT2")
            butct2.name = "Butyryl-CoA:acetate CoA-transferase"
            butct2.lower_bound = -1000
            butct2.upper_bound = 1000
            btcoa = model.metabolites.get_by_id("btcoa_c")
            ac = model.metabolites.get_by_id("ac_c")
            accoa = model.metabolites.get_by_id("accoa_c")
            butct2.add_metabolites({
                btcoa: -1, ac: -1, but_c: 1, accoa: 1
            })
            model.add_reactions([butct2])
            print("  Added BUTCT2 (synthesis)")

    # BUTt: but_c <-> but_e transport
    if "BUTt" not in [r.id for r in model.reactions]:
        butt = Reaction("BUTt")
        butt.name = "Butyrate transport"
        butt.lower_bound = -1000
        butt.upper_bound = 1000
        butt.add_metabolites({but_c: -1, but_e: 1})
        model.add_reactions([butt])
        print("  Added BUTt (transport)")

    # EX_but_e: exchange
    if "EX_but_e" not in [r.id for r in model.reactions]:
        ex_but = Reaction("EX_but_e")
        ex_but.name = "Butyrate exchange"
        ex_but.lower_bound = 0
        ex_but.upper_bound = 1000
        ex_but.add_metabolites({but_e: -1})
        model.add_reactions([ex_but])
        print("  Added EX_but_e (exchange)")

    cobra.io.write_sbml_model(model, output_path)
    print(f"  Saved: {output_path}")

    # Verify with permissive FBA
    with model:
        for rxn in model.exchanges:
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
        model.objective = "EX_but_e"
        model.objective_direction = "max"
        sol = model.optimize()
        print(f"  Verification FBA (permissive): {sol.objective_value:.4f} mmol/gDW/h")


print("=" * 70)
print("Roseburia")
print("=" * 70)
add_butyrate_pathway(
    "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2.xml",
    "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml"
)

print()
print("=" * 70)
print("Faecalibacterium")
print("=" * 70)
add_butyrate_pathway(
    "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2.xml",
    "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml"
)
