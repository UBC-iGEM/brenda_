from brendapyrser import BRENDA
from brenda_map import PHdata, TempData, BrendaVariant, BrendaData

brenda = BRENDA("brenda_db.txt")

rxn = brenda.reactions.get_by_id("4.2.1.1")

variants = []
default_organism = "unknown"
default_uniprot = None

organism_to_uniprot = {
    p.get("name").lower(): p.get("proteinID")
    for p in rxn.proteins.values()
    if p.get("name") and p.get("proteinID")
}

for substrate, measurements in rxn.KMvalues.items():
    for entry in measurements:
        km_value = entry.get("value")
        if not isinstance(km_value, (int, float)):
            continue

        species = entry.get("species", [default_organism])
        organism = species[0] if species else default_organism
        uniprot_id = organism_to_uniprot.get(organism.lower(), default_uniprot)

        ph_opt = next(
            (
                e["value"]
                for e in rxn.PH.get("optimum", [])
                if isinstance(e["value"], (int, float))
            ),
            None,
        )
        ph_range = next(
            (
                e["value"]
                for e in rxn.PH.get("range", [])
                if isinstance(e["value"], list)
            ),
            None,
        )

        temp_opt = next(
            (
                e["value"]
                for e in rxn.temperature.get("optimum", [])
                if isinstance(e["value"], (int, float))
            ),
            None,
        )
        temp_range = next(
            (
                e["value"]
                for e in rxn.temperature.get("range", [])
                if isinstance(e["value"], list)
            ),
            None,
        )

        variant = BrendaVariant(
            organism=organism,
            UniprotID=uniprot_id,
            KM=float(km_value),
            ph=PHdata(
                optimum=float(ph_opt) if ph_opt is not None else -1.0,
                range=tuple(float(v) for v in ph_range) if ph_range else (-1.0, -1.0),
            ),
            temp=TempData(
                optimum=float(temp_opt) if temp_opt is not None else -1.0,
                range=(
                    tuple(float(v) for v in temp_range) if temp_range else (-1.0, -1.0)
                ),
            ),
        )

        variants.append(variant)


brenda_data = BrendaData(entries=variants)
brenda_data.to_json(indent=4)

print(f"Saved {len(variants)} entries to 'brenda_serialized.json'")
