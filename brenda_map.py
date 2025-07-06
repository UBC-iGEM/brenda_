from numpy import float32
from dataclasses import dataclass

@dataclass
class PHdata:
    optimum: float32
    range: tuple[float32, float32]

@dataclass
class TempData:
    optimum: float32
    range: tuple[float32, float32]

@dataclass
class BrendaVariant:
    organism: str
    UniprotID: str | None
    KM: float32
    ph: PHdata
    temp: TempData

@dataclass
class BrendaData:
    entries: list[BrendaVariant]
