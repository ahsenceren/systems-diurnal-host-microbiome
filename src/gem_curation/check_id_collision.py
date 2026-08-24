import cobra

MODELS = {
    "Roseburia (original, pre-BSH)": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml",
    "Muribaculaceae (original, pre-BSH)": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
    "Faecalibacterium (never touched)": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
}

CHECK_IDS = ["chol_c", "chol_e", "tchol_c", "tchol_e", "taur_c", "taur_e",
             "mca_c", "mca_e", "tmca_c", "tmca_e"]

for name, path in MODELS.items():
    print("=" * 70)
    print(name)
    print("=" * 70)
    model = cobra.io.read_sbml_model(path)
    for cid in CHECK_IDS:
        if cid in model.metabolites:
            m = model.metabolites.get_by_id(cid)
            print(f"  FOUND (pre-existing): {cid}  ->  name='{m.name}'  formula={m.formula}")
        else:
            print(f"  not present: {cid}")
    chol_ex = [r.id for r in model.reactions if "chol" in r.id.lower() and r.id.startswith("EX_")]
    print(f"  All pre-existing EX_* reactions containing 'chol': {chol_ex}")
    for rid in chol_ex:
        r = model.reactions.get_by_id(rid)
        mets = list(r.metabolites.keys())
        print(f"    {rid}: {r.reaction}  |  metabolite name(s): {[m.name for m in mets]}")
    print()
