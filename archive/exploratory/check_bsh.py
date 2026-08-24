"""
check_bsh.py
=============

Checks whether the 3 existing GEMs already have BSH (bile salt hydrolase /
choloylglycine hydrolase, EC 3.5.1.24) reactions.
"""

import cobra

MODELS = {
    "Muribaculaceae": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/muribaculaceae_wol2.xml",
    "Roseburia": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml",
    "Faecalibacterium": "/home/aceren/diurnal_host_microbiome/data/gems_wol2/faecalibacterium_wol2_fixed.xml",
}

for name, path in MODELS.items():
    print("=" * 70)
    print(name)
    print("=" * 70)
    model = cobra.io.read_sbml_model(path)

    bile_mets = []
    for m in model.metabolites:
        mid_l = m.id.lower()
        name_l = (m.name or "").lower()
        if "chol" in mid_l or "bile" in name_l or "chol" in name_l:
            bile_mets.append(m)

    print(f"  Metabolites with 'chol'/'bile' in ID or name: {len(bile_mets)}")
    for m in bile_mets[:30]:
        print(f"    {m.id}  |  {m.name}  |  compartment={m.compartment}")

    bsh_rxns = [r for r in model.reactions if "bsh" in r.id.lower() or "cholylglycine" in (r.name or "").lower()
                or "choloylglycine" in (r.name or "").lower()]
    print(f"  Reactions with BSH-like name: {len(bsh_rxns)}")
    for r in bsh_rxns:
        print(f"    {r.id}: {r.name} | {r.reaction}")

    print()
