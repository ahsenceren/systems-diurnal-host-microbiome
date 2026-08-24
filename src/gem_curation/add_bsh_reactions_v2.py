import cobra
from cobra import Metabolite, Reaction

TARGET_MODELS = {
    "Roseburia": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml",
    "Muribaculaceae": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
}
OUT_SUFFIX = "_bsh2"

C_C = "C_c"
C_E = "C_e"

NEW_MET_IDS = ["tcholate_c", "tcholate_e", "tmurichol_c", "tmurichol_e",
               "cholate_c", "cholate_e", "murichol_c", "murichol_e",
               "taurine_c", "taurine_e"]


def add_bsh(model):
    collisions = [mid for mid in NEW_MET_IDS if mid in model.metabolites]
    if collisions:
        raise RuntimeError(f"ID COLLISION detected, aborting: {collisions} already exist in model")

    added_mets = []
    added_rxns = []

    def create_met(mid, name, compartment, formula=None):
        m = Metabolite(mid, name=name, compartment=compartment, formula=formula)
        model.add_metabolites([m])
        added_mets.append(mid)
        return m

    h2o_c = None
    if "h2o_c" in model.metabolites:
        h2o_c = model.metabolites.get_by_id("h2o_c")
    else:
        h2o_c = create_met("h2o_c", "Water", C_C, "H2O")

    tchol_c = create_met("tcholate_c", "Taurocholate", C_C, "C26H44NO7S")
    tchol_e = create_met("tcholate_e", "Taurocholate", C_E, "C26H44NO7S")
    tmca_c = create_met("tmurichol_c", "Tauro-beta-muricholate", C_C, "C26H44NO8S")
    tmca_e = create_met("tmurichol_e", "Tauro-beta-muricholate", C_E, "C26H44NO8S")
    chol_c = create_met("cholate_c", "Cholate", C_C, "C24H39O5")
    chol_e = create_met("cholate_e", "Cholate", C_E, "C24H39O5")
    mca_c = create_met("murichol_c", "Muricholate", C_C, "C24H39O6")
    mca_e = create_met("murichol_e", "Muricholate", C_E, "C24H39O6")
    taur_c = create_met("taurine_c", "Taurine", C_C, "C2H6NO3S")
    taur_e = create_met("taurine_e", "Taurine", C_E, "C2H6NO3S")

    def add_rxn(rid, name, mets, lb, ub):
        if rid in model.reactions:
            raise RuntimeError(f"REACTION ID COLLISION: {rid} already exists")
        r = Reaction(rid, name=name, lower_bound=lb, upper_bound=ub)
        r.add_metabolites(mets)
        model.add_reactions([r])
        added_rxns.append(rid)

    add_rxn("BSH_TCHOLATE", "Bile salt hydrolase (taurocholate)",
             {tchol_c: -1, h2o_c: -1, chol_c: 1, taur_c: 1}, 0, 1000)
    add_rxn("BSH_TMURICHOL", "Bile salt hydrolase (tauro-beta-muricholate)",
             {tmca_c: -1, h2o_c: -1, mca_c: 1, taur_c: 1}, 0, 1000)

    add_rxn("TCHOLATEt", "Taurocholate transport", {tchol_e: -1, tchol_c: 1}, -1000, 1000)
    add_rxn("TMURICHOLt", "Tauro-beta-muricholate transport", {tmca_e: -1, tmca_c: 1}, -1000, 1000)
    add_rxn("CHOLATEt", "Cholate transport", {chol_c: -1, chol_e: 1}, -1000, 1000)
    add_rxn("MURICHOLt", "Muricholate transport", {mca_c: -1, mca_e: 1}, -1000, 1000)
    add_rxn("TAURINEt", "Taurine transport", {taur_c: -1, taur_e: 1}, -1000, 1000)

    add_rxn("EX_tcholate_e", "Taurocholate exchange", {tchol_e: -1}, 0, 1000)
    add_rxn("EX_tmurichol_e", "Tauro-beta-muricholate exchange", {tmca_e: -1}, 0, 1000)
    add_rxn("EX_cholate_e", "Cholate exchange", {chol_e: -1}, 0, 1000)
    add_rxn("EX_murichol_e", "Muricholate exchange", {mca_e: -1}, 0, 1000)
    add_rxn("EX_taurine_e", "Taurine exchange", {taur_e: -1}, 0, 1000)

    return added_mets, added_rxns


for name, path in TARGET_MODELS.items():
    print("=" * 70)
    print(name)
    print("=" * 70)
    model = cobra.io.read_sbml_model(path)
    print(f"Before: {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")

    added_mets, added_rxns = add_bsh(model)

    print(f"After:  {len(model.metabolites)} metabolites, {len(model.reactions)} reactions")
    print(f"New metabolites added: {len(added_mets)} -> {added_mets}")
    print(f"New reactions added: {len(added_rxns)} -> {added_rxns}")

    for mid in ["tcholate_e", "cholate_e"]:
        m = model.metabolites.get_by_id(mid)
        print(f"  sanity: {mid}.compartment = '{m.compartment}' (expect 'C_e')")

    growth = model.slim_optimize()
    print(f"Sanity check -- growth (default/permissive bounds unchanged): {growth}")

    out_path = path.replace(".xml", f"{OUT_SUFFIX}.xml")
    cobra.io.write_sbml_model(model, out_path)
    print(f"Saved: {out_path}")
    print()
