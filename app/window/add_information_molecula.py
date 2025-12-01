from PyQt6.QtWidgets import QMessageBox
from app.window.ui.add_infornation_molecula import UIAddInformationMolecula
import traceback


class AddInformationMolecula(UIAddInformationMolecula):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.btn_accept.clicked.connect(self._validate_data)
        self.btn_cancel.clicked.connect(self.reject)

    def _validate_data(self):
        if not self.name_molecula.text():
            QMessageBox.critical(self, "Ошибка валидации", "Поле названия молекулы пустое")
            return

        if not self.smiles.text():
            QMessageBox.critical(self, "Ошибка валидации", "Поле SMILES пустое")
            return

        if not self.__validate_peaks(self.uv_peaks.toPlainText()):
            QMessageBox.critical(self, "Ошибка валидации", "В поле УФ не правильные данные")
            return

        if not self.__validate_peaks(self.ir_peaks.toPlainText()):
            QMessageBox.critical(self, "Ошибка валидации", "В поле ИК не правильные данные")
            return

        self.accept()

    def show_modal(self):
        if self.exec():
            print("Что-то добавили")
        else:
            pass

    def __validate_peaks(self, text: str) -> bool:
        if not text:
            return False

        for line in text.split("\n"):
            if not line:
                continue

            try:
                _numbers = tuple(map(int, line.split(" ")))
                if len(_numbers) != 2:
                    return False
            except ValueError:
                return False

        return True
