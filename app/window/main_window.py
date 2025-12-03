import os.path
import sys
import json

from PyQt6.QtCore import QAbstractListModel, Qt, QModelIndex

from typing import List, Optional, Union

from PyQt6.QtWidgets import QMenu, QFileDialog

from app.molecula import Molecula, Params, Peaks
from app.window.ui.main_window import UIMainWindow
from app.window.add_information_molecula import AddInformationMolecula
from pathlib import Path

from app.window.ui.plot_canvas import FullGraphics, IRGraphic, UVGraphic, IRUVCombinate, StructuraImg


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

    def get_molecula(self, index: QModelIndex) -> Optional[Molecula]:
        if 0 <= index.row() < self.rowCount():
            return self._list_molecula[index.row()]
        return None

    def append_row(self, _molecula: Molecula):
        self.beginInsertRows(QModelIndex(), self.rowCount(), self.rowCount())
        self._list_molecula.append(_molecula)
        self.endInsertRows()

    def remove_row(self, index):
        if 0 <= index.row() < self.rowCount():
            self.beginRemoveRows(QModelIndex(), index.row(), index.row())
            molecula = self._list_molecula[index.row()]
            #Сигналим что удалились
            molecula.i_am_removed.emit()
            del self._list_molecula[index.row()]
            self.endRemoveRows()
            return True
        return False



class MainWindow(UIMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._list_molecula_model = ListMoleculaModel()
        self.list_molecula_widget.setModel(self._list_molecula_model)
        self.list_molecula_widget.doubleClicked.connect(self._open_widget_graph)
        self.exit_program.triggered.connect(lambda: sys.exit(0))
        self.load_data_action.triggered.connect(self._show_select_json_file)
        # self.export_data_action.triggered.connect(self._export_molecula)

        self.load_default_value()

        self.list_molecula_widget.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.list_molecula_widget.customContextMenuRequested.connect(self._show_context_menu)

    def _show_select_json_file(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Открыть файл',
                                               './', "JSON файл (*.json);;Все файлы (*)")

        self._load_data_from_file(Path(fname))

    def _show_context_menu(self, pos):
        index = self.list_molecula_widget.currentIndex()
        molecula = self._list_molecula_model.get_molecula(index)
        if molecula:
            menu = QMenu(self)

            menu_graphics = menu.addMenu("Графики")
            open_only_ir = menu_graphics.addAction("ИК")
            open_only_ir.triggered.connect(lambda: self._open_only_ir(molecula))
            open_only_uv = menu_graphics.addAction("УФ")
            open_only_uv.triggered.connect(lambda: self._open_only_uv(molecula))
            combinate_graphics = menu_graphics.addAction("Комбинированый график УФ и ИК")
            combinate_graphics.triggered.connect(lambda: self._open_combinate_graphics(molecula))

            open_struct = menu.addAction("Открыть структуру")
            open_struct.triggered.connect(lambda: self._open_struct(molecula))
            delete_molecula = menu.addAction("Удалить")
            delete_molecula.triggered.connect(lambda: self._list_molecula_model.remove_row(index))
            menu.exec(self.list_molecula_widget.viewport().mapToGlobal(pos))

    def _open_struct(self, molecula):
        ir_graphic = StructuraImg(self, molecula, dpi=500, width=400, height=400)
        self.tab_widget.addTab(ir_graphic, F"Структура {molecula.name}")
        self.tab_widget.setCurrentIndex(self.tab_widget.count() - 1)
        ir_graphic.set_function_close_table(lambda: self.tab_widget.removeTab(self.tab_widget.count() - 1))

    def _open_combinate_graphics(self, molecula):
        ir_graphic = IRUVCombinate(self, molecula)
        self.tab_widget.addTab(ir_graphic, F"Комбинированый график УФ и ИК {molecula.name}")
        self.tab_widget.setCurrentIndex(self.tab_widget.count() - 1)
        ir_graphic.set_function_close_table(lambda: self.tab_widget.removeTab(self.tab_widget.count() - 1))

    def _open_only_uv(self, molecula):
        ir_graphic = UVGraphic(self, molecula)
        self.tab_widget.addTab(ir_graphic, F"УФ спектр {molecula.name}")
        self.tab_widget.setCurrentIndex(self.tab_widget.count() - 1)
        ir_graphic.set_function_close_table(lambda: self.tab_widget.removeTab(self.tab_widget.count() - 1))

    def _open_only_ir(self, molecula):
        ir_graphic = IRGraphic(self, molecula)
        self.tab_widget.addTab(ir_graphic, F"ИК спектр {molecula.name}")
        self.tab_widget.setCurrentIndex(self.tab_widget.count() - 1)
        ir_graphic.set_function_close_table(lambda: self.tab_widget.removeTab(self.tab_widget.count() - 1))

    def show_dialog_add_molecula(self):
        _add_information_molecula = AddInformationMolecula(self)

        molecula = _add_information_molecula.show_modal()
        if molecula is not None:
            self._list_molecula_model.append_row(molecula)


    def _load_data_from_file(self, file: Path):

        if file.is_file():
            with open(file, 'r', encoding="UTF-8") as file:
                _json = json.load(file)

                for key, val in _json.items():
                    molecula = Molecula(
                        name=key,
                        smiles=val['SMILES'],
                        uv_params=Params(width=val["uv_params"]["width"],
                                         peaks = self.__get_peaks(val["uv_params"]["peaks"])),
                        ir_params=Params(peaks = self.__get_peaks(val["ir_params"]["peaks"]))
                    )

                    self._list_molecula_model.append_row(molecula)

    def __get_peaks(self, _list: List[List[Union[int, float]]]) -> List[Peaks]:
        _result: List[Peaks] = []
        for peak, intensity in _list:
            _result.append(Peaks(peak=peak, intensity=intensity))
        return _result

    def  _open_widget_graph(self, index: QModelIndex):
        molecula = self._list_molecula_model.get_molecula(index)
        _plot_canvas = FullGraphics(molecula=molecula)
        self.tab_widget.addTab(_plot_canvas, molecula.name)

        self.tab_widget.setCurrentIndex(self.tab_widget.count() - 1)
        _plot_canvas.set_function_close_table(lambda: self.tab_widget.removeTab(self.tab_widget.count() - 1))

    def load_default_value(self):
        root_path = Path(os.path.abspath(os.path.dirname(sys.argv[0])))
        path_file = root_path / "example_data.json"

        self._load_data_from_file(path_file)