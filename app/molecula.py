import numpy as np
from typing import List, Union
from dataclasses import dataclass

from rdkit import Chem
from rdkit.Chem import Draw

import app.default_variable as dv

@dataclass
class Peaks:
    intensity: Union[int, float]
    peak: Union[int, float]

@dataclass
class Params:
    peaks: List[Peaks]
    width: int = 0

@dataclass
class Molecula:
    name: str
    smiles: str
    ir_params: Params
    uv_params: Params

    def get_ir_spectrum(self):
        """Генерация модельного ИК спектра"""
        wavenumber = np.linspace(dv.IR_START, dv.IR_STOP, dv.IR_NUM)
        transmittance = np.ones_like(wavenumber)

        for ir_peak in self.ir_params.peaks:
            transmittance -= ir_peak.intensity * np.exp(-((wavenumber - ir_peak.peak) / 20) ** 2)

        transmittance = np.clip(transmittance, 0, 1)
        return wavenumber, transmittance

    def get_uv_spectrum(self):
        """Генерация модельного УФ спектра"""
        wavelength = np.linspace(dv.UV_START, dv.UV_STOP, dv.UV_NUM)
        absorbance = np.zeros_like(wavelength)

        for peak in self.uv_params.peaks:
            absorbance += peak.intensity * np.exp(-((wavelength - peak.peak) / self.uv_params.width) ** 2)

        return wavelength, absorbance

    def draw_molecule(self, height=300, width=300):
        """Отрисовка молекулы"""

        mol = Chem.MolFromSmiles(self.smiles)
        if mol:
            img = Draw.MolToImage(mol, size=(height, width))
            return img
        return None