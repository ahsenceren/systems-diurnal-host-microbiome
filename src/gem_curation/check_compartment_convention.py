import cobra

model = cobra.io.read_sbml_model("/home/aceren/diurnal_host_microbiome/data/gems_wol2/roseburia_wol2_fixed.xml")

for mid in ["so4_e", "h2o_e", "glc__D_e", "but_e"]:
    if mid in model.metabolites:
        m = model.metabolites.get_by_id(mid)
        print(f"{mid}: compartment attribute = '{m.compartment}'")

print()
print("model.compartments dict:", model.compartments)

print()
for mid in ["h2o_c", "atp_c", "but_c"]:
    if mid in model.metabolites:
        m = model.metabolites.get_by_id(mid)
        print(f"{mid}: compartment attribute = '{m.compartment}'")
