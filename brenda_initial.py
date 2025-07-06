import numpy as np
from brendapyrser import BRENDA

dataFile = "brenda_db.txt"
brenda = BRENDA(dataFile)

ec_number = "4.2.1.1"
carbonic_anhydrase = brenda.reactions.get_by_id(ec_number)

if carbonic_anhydrase:
    print(f"Found enzyme: {carbonic_anhydrase.name} ({ec_number})")

    organism_to_id_map = {}
    for protein_info in carbonic_anhydrase.proteins.values():
        name = protein_info['name']
        protein_id = protein_info['proteinID']
        if name not in organism_to_id_map:
            organism_to_id_map[name] = []
        if protein_id:
            organism_to_id_map[name].append(protein_id)

    print("--- Substrate, Organism, UniProt ID, KM Value (mM) ---")

    km_values_dict = carbonic_anhydrase.KMvalues

    for substrate, measurement_list in km_values_dict.items():
        for measurement_data in measurement_list:
            km_value = measurement_data['value']

            organism_names = measurement_data['species']

            for org_name in organism_names:
                uniprot_ids = organism_to_id_map.get(org_name, ['N/A'])

                print(
                    f"Substrate: {substrate}, "
                    f"Organism: {org_name}, "
                    f"UniProt ID: {', '.join(uniprot_ids)}, "
                    f"KM Value: {km_value}"
                )
else:
    print(f"Enzyme with EC number {ec_number} not found in the database.")
