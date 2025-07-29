from dataclasses import dataclass
from dataclasses_json import dataclass_json
from typing import Optional, Tuple
from numpy import float32, float64
from brendapyrser import BRENDA

@dataclass_json
@dataclass
class PHdata:
    optimum: Optional[float] = None
    range: Optional[Tuple[float, float]] = None

@dataclass_json
@dataclass
class TempData:
    optimum: Optional[float] = None
    range: Optional[Tuple[float, float]] = None

@dataclass_json
@dataclass
class BrendaVariant:
    UniprotID: Optional[str]
    KM: float
    ph: PHdata
    temp: TempData

@dataclass_json
@dataclass
class BrendaData:
    entries: dict[str, BrendaVariant]

brenda = BRENDA("brenda_db.txt")
rxn = brenda.reactions.get_by_id("4.2.1.1")

organism_to_uniprot = {
    p.get("name").lower(): p.get("proteinID")
    for p in rxn.proteins.values()
    if p.get("name") and p.get("proteinID")
}

entries: dict[str, BrendaVariant] = {}

for substrate, measurements in rxn.KMvalues.items():
    for entry in measurements:
        km_value = entry.get("value")
        if not isinstance(km_value, (int, float, float32, float64)):
            continue

        species = entry.get("species", [])
        if not species:
            continue
        organism = species[0]
        uniprot_id = organism_to_uniprot.get(organism.lower())

        if organism not in entries:
            entries[organism] = BrendaVariant(
                UniprotID=uniprot_id,
                KM=float(km_value),
                ph=PHdata(),
                temp=TempData()
            )

for ph_type in ["optimum", "range"]:
    for entry in rxn.PH.get(ph_type, []):
        value = entry.get("value")
        if value is None:
            continue
        species = entry.get("species", [])
        if not species:
            continue
        organism = species[0]
        if organism not in entries:
            continue
        if ph_type == "optimum" and isinstance(value, (int, float, float32, float64)):
            entries[organism].ph.optimum = float(value)
        elif ph_type == "range" and isinstance(value, list) and len(value) == 2:
            entries[organism].ph.range = tuple(float(v) for v in value)

for temp_type in ["optimum", "range"]:
    for entry in rxn.temperature.get(temp_type, []):
        value = entry.get("value")
        if value is None:
            continue
        species = entry.get("species", [])
        if not species:
            continue
        organism = species[0]
        if organism not in entries:
            continue
        if temp_type == "optimum" and isinstance(value, (int, float, float32, float64)):
            entries[organism].temp.optimum = float(value)
        elif temp_type == "range" and isinstance(value, list) and len(value) == 2:
            entries[organism].temp.range = tuple(float(v) for v in value)

brenda_data = BrendaData(entries=entries)

with open("brenda_serialized.json", "w") as f:
    f.write(brenda_data.to_json(indent=4))

print(f"Saved {len(entries)} organism-specific entries to 'brenda_serialized.json'")
