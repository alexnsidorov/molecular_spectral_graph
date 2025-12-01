from typing import List, Union
from dataclasses import dataclass, field

@dataclass
class MoleculaPeaks:
    Intensity: int
    Peak: int

@dataclass
class Molecula:
    name: str
    smiles: str
    width_peak: int
    peaks: List[MoleculaPeaks]
