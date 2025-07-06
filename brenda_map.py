from dataclasses import dataclass
from dataclasses_json import dataclass_json

@dataclass_json
@dataclass
class PHdata:
    optimum: float
    range: tuple[float, float]

@dataclass_json
@dataclass
class TempData:
    optimum: float
    range: tuple[float, float]

@dataclass_json
@dataclass
class BrendaVariant:
    organism: str
    UniprotID: str | None
    KM: float
    ph: PHdata
    temp: TempData

@dataclass_json
@dataclass
class BrendaData:
    entries: list[BrendaVariant]
