import os.path
import sys
import json

from PyQt6.QtCore import QAbstractListModel, Qt, QModelIndex

from typing import List
from app.molecula import Molecula, Params, Peaks
from app.window.ui.main_window import UIMainWindow
from app.window.add_information_molecula import AddInformationMolecula
from pathlib import Path



class ListMoleculaModel(QAbstractListModel):
    def __init__(self):
        super().__init__()

        self._list_molecula: List[Molecula] = []

    def rowCount(self, parent = None):
        return len(self._list_molecula)

    def data(self, index, role = Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DisplayRole:
            _molecula: Molecula = self._list_molecula[index.row()]
            return _molecula.name

        return None

    def append_row(self, _molecula: Molecula):
        self.beginInsertRows(QModelIndex(), self.rowCount(), self.rowCount())
        self._list_molecula.append(_molecula)
        self.endInsertRows()

    def removeRow(self, position, parent=QModelIndex()):
        if 0 <= position < self.rowCount():
            self.beginRemoveRows(parent, position, position)
            del self._list_molecula[position]
            self.endRemoveRows()
            return True
        return False

class MainWindow(UIMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._list_molecula_model = ListMoleculaModel()
        self.list_molecula_widget.list_molecula_view.setModel(self._list_molecula_model)

        self.exit_program.triggered.connect(lambda: sys.exit(0))

        self.list_molecula_widget.add_molecula.clicked.connect(self.show_dialog_add_molecula)
        self.load_default_value()

    def show_dialog_add_molecula(self):
        _add_information_molecula = AddInformationMolecula(self)

        molecula = _add_information_molecula.show_modal()
        if molecula is not None:
            self._list_molecula_model.append_row(molecula)


    def load_default_value(self):
        root_path = Path(os.path.abspath(os.path.dirname(sys.argv[0])))
        path_file = root_path / "example_data.json"


        if path_file.is_file():
            with open(path_file, 'r', encoding="UTF-8") as file:
                _json = json.load(file)

                for key, val in _json.items():
                    molecula = Molecula(
                        name=key,
                        smiles=val['SMILES'],
                        uv_params=Params(width=val["uv_params"]["width"],
                                         peaks = [Peaks(peak=intensity,
                                                        intensity=peak) for peak, intensity in val["uv_params"]["peaks"]]),
                        ir_params=Params(peaks = [Peaks(peak=intensity,
                                                        intensity=peak) for peak, intensity in val["ir_params"]["peaks"]])
                    )

                    self._list_molecula_model.append_row(molecula)
