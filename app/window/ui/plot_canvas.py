from PyQt6.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from app.molecula import Molecula
from app.default_variable import *


class PlotCanvas(QWidget):
    def __init__(self, parent=None, molecula: Molecula = None, width=5, height=4, dpi=100):
        super().__init__(parent)
        base_layout = QVBoxLayout()
        self._figure = Figure(figsize=(width, height), dpi=dpi)
        self._canvas = FigureCanvas(self._figure)
        self._molecula = molecula
        # УФ
        self._wavelength, self._absorbance = self._molecula.get_uv_spectrum()
        # ИК
        self._wavenumber, self._transmittance = self._molecula.get_ir_spectrum()

        self._create_ir_spectre()
        self._create_uv_spectre()
        self._create_structure_molecula()
        self._create_combinate()

        base_layout.addWidget(self._canvas)
        self.setLayout(base_layout)

    def _create_ir_spectre(self):
        self._ir_plot = self._figure.add_subplot(2, 2, 2)
        self._ir_plot.plot(self._wavenumber, self._transmittance * 100, 'r-', linewidth=2)
        self._ir_plot.set_title(F"ИК спектр {self._molecula.name}")
        self._ir_plot.set_ylabel("Пропускание (%)")
        self._ir_plot.set_xlabel("Волновое число (см⁻¹)")
        self._ir_plot.grid(True, alpha=0.3)


    def _create_uv_spectre(self):
        self._uv_plot = self._figure.add_subplot(2, 2, 1)
        self._uv_plot.plot(self._wavelength, self._absorbance, 'b-', linewidth=2)
        self._uv_plot.set_title(F"УФ спектр {self._molecula.name}")
        self._uv_plot.set_ylabel("Молярный коэффициент поглощения")
        self._uv_plot.set_xlabel("Длина волны (нм)")
        self._uv_plot.grid(True, alpha=0.3)
        self._ir_plot.set_xlim(UV_START, UV_STOP)

    def _create_structure_molecula(self):
        self._structure_molecula = self._figure.add_subplot(2, 2, 3)
        img = self._molecula.draw_molecule()
        if img is not None:
            self._structure_molecula.imshow(img)
            self._structure_molecula.axis('off')
            self._structure_molecula.set_title(f'Структура {self._molecula.name}')

    def _create_combinate(self):
        self._combinate_plot = self._figure.add_subplot(2, 2, 4)
        self._combinate_plot.plot(self._wavelength, self._absorbance, 'b-', linewidth=2, label='УФ спектр')
        _combinate_twin = self._combinate_plot.twinx()
        _combinate_twin.plot(self._wavenumber, self._transmittance * 100, 'r-', linewidth=2, label='ИК спектр')
        self._combinate_plot.set_xlabel('Длина волны (нм) / Волновое число (см⁻¹)')
        self._combinate_plot.set_ylabel('Молярный коэффициент поглощения', color='b')
        _combinate_twin.set_ylabel('Пропускание (%)', color='r')
        self._combinate_plot.set_title(f'Комбинированные спектры {self._molecula.name}')
        self._combinate_plot.grid(True, alpha=0.3)
        self._combinate_plot.set_xlim(UV_START, UV_STOP)


