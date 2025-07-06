import json
from numpy import float32
from brendapyrser import BRENDA
from brenda_map import PHdata, TempData, BrendaVariant, BrendaData
brenda = BRENDA("brenda_db.txt")

rxn = brenda.reactions.get_by_id("4.2.1.1")

variants = []
default_organism = "unknown"
default_uniprot = None

organism_to_uniprot = {
    p.get("organism", "").lower(): p.get("proteinID")
    for p in rxn.proteins.values()
    if p.get("organism") and p.get("proteinID")
}

for substrate, measurements in rxn.KMvalues.items():
    for entry in measurements:
        km_value = entry.get("value")
        if not isinstance(km_value, (int, float)):
            continue

        species = entry.get("species", [default_organism])
        organism = species[0] if species else default_organism
        uniprot_id = organism_to_uniprot.get(organism, default_uniprot)

        
        ph_opt = next((e["value"] for e in rxn.PH.get("optimum", []) if isinstance(e["value"], (int, float))), None)
        ph_range = next((e["value"] for e in rxn.PH.get("range", []) if isinstance(e["value"], list)), None)

        temp_opt = next((e["value"] for e in rxn.temperature.get("optimum", []) if isinstance(e["value"], (int, float))), None)
        temp_range = next((e["value"] for e in rxn.temperature.get("range", []) if isinstance(e["value"], list)), None)

    
        variant = BrendaVariant(
            organism=organism,
            UniprotID=uniprot_id,
            KM=float32(km_value),
            ph=PHdata(
                optimum=float32(ph_opt) if ph_opt is not None else float32(-1),
                range=tuple(float32(v) for v in ph_range) if ph_range else (float32(-1), float32(-1))
            ),
            temp=TempData(
                optimum=float32(temp_opt) if temp_opt is not None else float32(-1),
                range=tuple(float32(v) for v in temp_range) if temp_range else (float32(-1), float32(-1))
            )
        )

        variants.append(variant)


brenda_data = BrendaData(entries=variants)


with open("brenda_serialized.json", "w") as f:
    def convert(obj):
        if hasattr(obj, "__dict__"):
            return obj.__dict__
        elif isinstance(obj, float32):
            return float(obj)
        return str(obj)

    json.dump(brenda_data, f, indent=2, default=convert)

print(f"Saved {len(variants)} entries to 'brenda_serialized.json'")