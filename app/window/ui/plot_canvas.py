from PyQt6.QtWidgets import QWidget, QVBoxLayout, QToolBar, QLabel, QStatusBar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT
from matplotlib.backends.backend_qt import LocationEvent

from matplotlib.figure import Figure
import numpy as np
from app.molecula import Molecula
from app.default_variable import *

class PlotCreator:
    @staticmethod
    def create_uv_spectre(plot, molecula_name: str, wavelength: np.ndarray, absorbance: np.ndarray):
        plot.plot(wavelength, absorbance, 'b-', linewidth=2)
        plot.set_title(F"УФ спектр {molecula_name}")
        plot.set_ylabel("Молярный коэффициент поглощения")
        plot.set_xlabel("Длина волны (нм)")
        plot.grid(True, alpha=0.3)
        plot.set_xlim(UV_START, UV_STOP)


    @staticmethod
    def create_ir_spectre(plot, molecula_name: str, wavenumber: np.ndarray, transmittance: np.ndarray):
        plot.plot(wavenumber, transmittance * 100, 'r-', linewidth=2)
        plot.set_xlabel('Волновое число (см⁻¹)')
        plot.set_ylabel('Пропускание (%)')
        plot.set_title(f'ИК спектр {molecula_name}')
        plot.grid(True, alpha=0.3)
        plot.set_xlim(IR_START, IR_STOP)
        plot.invert_xaxis()

    @staticmethod
    def create_structure_molecula(plot, molecula: Molecula):
        img = molecula.draw_molecule(500, 500)
        if img is not None:
            plot.imshow(img)
            plot.axis('off')
            plot.set_title(f'Структура {molecula.name}')

    @staticmethod
    def create_combinate_graph(plot, molecula_name: str, wavelength, absorbance, wavenumber, transmittance):
        plot.plot(wavelength, absorbance, 'b-', linewidth=2, label='УФ спектр')
        plot.set_xlabel('Длина волны (нм) / Волновое число (см⁻¹)')
        plot.set_ylabel('Молярный коэффициент поглощения', color='b')

        plot_twin = plot.twinx()
        plot_twin.plot(wavenumber, transmittance * 100, 'r-', linewidth=2, label='ИК спектр')
        plot_twin.set_xlim(IR_START, IR_STOP)
        plot_twin.set_ylabel('Пропускание (%)', color='r')

        plot.set_title(f'Комбинированные спектры {molecula_name}')
        plot.grid(True, alpha=0.3)
        plot.set_xlim(min((UV_START, IR_START)), max((UV_STOP, IR_STOP)))

        lines1, labels1 = plot.get_legend_handles_labels()
        lines2, labels2 = plot_twin.get_legend_handles_labels()
        plot.legend(lines1 + lines2, labels1 + labels2, loc='upper right')


class BaseCanvas(QWidget):
    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent)

        self._remove_tab = None
        base_layout = QVBoxLayout()

        self._figure = Figure(figsize=(width, height), dpi=dpi)
        self._figure.set_frameon(False)
        self._canvas = FigureCanvas(self._figure)
        self._canvas.mpl_connect('motion_notify_event', self._mouse_move)

        navigator = NavigationToolbar2QT(self._canvas, self, coordinates=False)

        self._molecula = molecula
        self._molecula.i_am_removed.connect(self._rem_tab)

        self._status_label = QLabel(self)
        self._tool_bar = QToolBar()
        self._tool_bar.addWidget(self._status_label)

        base_layout.addWidget(navigator)
        base_layout.addWidget(self._canvas)
        base_layout.addWidget(self._tool_bar)

        self.setLayout(base_layout)

    def _mouse_move(self, event: LocationEvent):
        if event.inaxes:
            self._status_label.setText(F"{event.inaxes.get_title()} x: {round(event.xdata, 2)} y: {round(event.ydata, 2)}")

    def set_function_close_table(self, fun):
        self._remove_tab = fun

    def _rem_tab(self):
        if self._remove_tab is not None:
            self._remove_tab()

class FullGraphics(BaseCanvas):

    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent=parent, molecula=molecula, width=width, height=height, dpi=dpi)

        # УФ
        self._wavelength, self._absorbance = self._molecula.get_uv_spectrum()
        # ИК
        self._wavenumber, self._transmittance = self._molecula.get_ir_spectrum()

        ax1 = self._figure.add_subplot(2, 2, 1)
        ax2 = self._figure.add_subplot(2, 2, 2)
        ax3 = self._figure.add_subplot(2, 2, 3)
        ax4 = self._figure.add_subplot(2, 2, 4)

        PlotCreator.create_uv_spectre(ax1, self._molecula.name, self._wavelength, self._absorbance)
        PlotCreator.create_ir_spectre(ax2, self._molecula.name, self._wavenumber, self._transmittance)
        PlotCreator.create_structure_molecula(ax3, self._molecula)
        PlotCreator.create_combinate_graph(ax4, self._molecula.name, self._wavelength, self._absorbance,
                                        self._wavenumber, self._transmittance)


class IRGraphic(BaseCanvas):
    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent=parent, molecula=molecula, width=width, height=height, dpi=dpi)
        PlotCreator.create_ir_spectre(self._figure.add_subplot(1, 1, 1), self._molecula.name, *molecula.get_ir_spectrum())

class UVGraphic(BaseCanvas):
    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent=parent, molecula=molecula, width=width, height=height, dpi=dpi)
        PlotCreator.create_uv_spectre(self._figure.add_subplot(1, 1, 1), self._molecula.name, *molecula.get_uv_spectrum())

class IRUVCombinate(BaseCanvas):
    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent=parent, molecula=molecula, width=width, height=height, dpi=dpi)
        PlotCreator.create_combinate_graph(self._figure.add_subplot(1, 1, 1), self._molecula.name,
                                           *molecula.get_ir_spectrum(), *molecula.get_uv_spectrum())

class StructuraImg(BaseCanvas):
    def __init__(self, parent=None, molecula: Molecula = None, width=16, height=6, dpi=72):
        super().__init__(parent=parent, molecula=molecula, width=width, height=height, dpi=dpi)
        PlotCreator.create_structure_molecula(self._figure.add_subplot(1, 1, 1), molecula)