import cobra

MODELS = {
    "Muribaculaceae": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
    "Roseburia": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2.xml",
    "Faecalibacterium": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2.xml",
}

for name, path in MODELS.items():
    print("=" * 70)
    print(f"{name}")
    print("=" * 70)

    model = cobra.io.read_sbml_model(path)

    has_but_c = "but_c" in [m.id for m in model.metabolites]
    has_but_e = "but_e" in [m.id for m in model.metabolites]
    print(f"  but_c (intracellular) exists: {has_but_c}")
    print(f"  but_e (extracellular) exists: {has_but_e}")

    transport_rxns = []
    if has_but_c and has_but_e:
        but_c_met = model.metabolites.get_by_id("but_c")
        but_e_met = model.metabolites.get_by_id("but_e")
        for rxn in model.reactions:
            mets_in_rxn = set(rxn.metabolites.keys())
            if but_c_met in mets_in_rxn and but_e_met in mets_in_rxn:
                transport_rxns.append(rxn.id)
    print(f"  Transport reactions (but_c<->but_e): {transport_rxns if transport_rxns else 'NONE FOUND'}")

    has_ex_but = "EX_but_e" in [r.id for r in model.reactions]
    print(f"  EX_but_e exchange reaction exists: {has_ex_but}")

    synthesis_candidates = []
    if has_but_c:
        but_c_met = model.metabolites.get_by_id("but_c")
        for rxn in but_c_met.reactions:
            coef = rxn.metabolites.get(but_c_met, 0)
            if coef > 0:
                synthesis_candidates.append((rxn.id, rxn.reaction))
    print(f"  Reactions producing but_c ({len(synthesis_candidates)} found):")
    for rxn_id, rxn_str in synthesis_candidates:
        print(f"    {rxn_id}: {rxn_str}")

    if has_ex_but:
        with model:
            for rxn in model.exchanges:
                rxn.lower_bound = -1000
                rxn.upper_bound = 1000
            model.objective = "EX_but_e"
            model.objective_direction = "max"
            solution = model.optimize()
            print(f"  FBA max butyrate secretion (permissive medium): "
                  f"{solution.objective_value:.4f} mmol/gDW/h "
                  f"(status: {solution.status})")
            if solution.objective_value < 1e-6:
                print(f"  *** WARNING: pathway present but NOT flux-connected even under permissive medium ***")
    else:
        print(f"  *** Cannot test FBA -- no EX_but_e reaction ***")
    print()
