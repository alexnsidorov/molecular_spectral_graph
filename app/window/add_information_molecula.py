from PyQt6.QtWidgets import QMessageBox, QDialog

from app.molecula import Molecula, Peaks, Params
from app.window.ui.add_infornation_molecula import UIAddInformationMolecula
from rdkit import Chem
from typing import Tuple, List, Optional


class AddInformationMolecula(UIAddInformationMolecula):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.btn_accept.clicked.connect(self._validate_data)
        self.btn_cancel.clicked.connect(self.reject)

    def _validate_data(self):
        if not self.name_molecula.text():
            QMessageBox.critical(self, "Ошибка валидации", "Поле названия молекулы пустое")
            return

        _is_validate_smale, massage = self.__validate_smiles(self.smiles.text())
        if not _is_validate_smale:
            QMessageBox.critical(self, "Ошибка валидации SMILES", massage)
            return

        if not self.__validate_peaks(self.uv_peaks.toPlainText()):
            QMessageBox.critical(self, "Ошибка валидации", "В поле УФ не правильные данные")
            return

        if not self.__validate_peaks(self.ir_peaks.toPlainText()):
            QMessageBox.critical(self, "Ошибка валидации", "В поле ИК не правильные данные")
            return

        self.accept()

    def show_modal(self) -> Optional[Molecula]:
        if self.exec() == QDialog.DialogCode.Accepted:
            ir_pears = self._parse_peaks(self.ir_peaks.toPlainText())
            uv_pears = self._parse_peaks(self.uv_peaks.toPlainText())

            _molecula = Molecula(
                name=self.name_molecula.text(),
                smiles=self.smiles.text(),
                ir_params=Params(peaks=ir_pears),
                uv_params=Params(peaks=uv_pears, width=self.width_peak.value())
            )
            return _molecula

        return None

    def __validate_peaks(self, text: str) -> bool:
        if not text:
            return False

        for line in text.split("\n"):
            if not line:
                continue

            try:
                if len(tuple(map(float, line.split(" ")))) != 2:
                    return False
            except ValueError:
                return False

        return True

    def __validate_smiles(self, smiles) -> Tuple[bool, str]:
        # Шаг 1: Синтаксическая проверка
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Неверный синтаксис"

        # Chem.SanitizeMol(mol, catchErrors=True)
        # # Шаг 3: Проверка валентностей
        # for atom in mol.GetAtoms():
        #     if not atom.HasProp('valences'):
        #         return False, "Некорректная валентность"

        return True, "SMILES корректен"

    def _parse_peaks(self, text: str) -> List[Peaks]:
        list_peaks = []
        for line in text.split("\n"):
            data = line.strip().split(" ")
            if line and len(data) == 2:
                peak, intensity = data
                list_peaks.append(Peaks(peak=int(peak), intensity=float(intensity)))

        return list_peaks
