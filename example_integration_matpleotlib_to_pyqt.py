import sys
from PyQt6.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Четыре графика")
        self.setGeometry(100, 100, 800, 600)

        # Создание центрального виджета
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)

        # Создание макета
        self.layout = QVBoxLayout()
        self.central_widget.setLayout(self.layout)

        # Создание фигуры с 4 графиками
        self.fig = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.fig)
        self.layout.addWidget(self.canvas)

        # Создание подграфиков
        self.ax1 = self.fig.add_subplot(2, 2, 1)
        self.ax2 = self.fig.add_subplot(2, 2, 2)
        self.ax3 = self.fig.add_subplot(2, 2, 3)
        self.ax4 = self.fig.add_subplot(2, 2, 4)

        # Генерация данных
        self.x = np.linspace(0, 10, 100)
        self.y1 = np.sin(self.x)
        self.y2 = np.cos(self.x)
        self.y3 = np.tan(self.x)
        self.y4 = np.exp(-self.x) * np.sin(self.x)

        # Построение графиков
        self.plot_graphs()

    def plot_graphs(self):
        # График 1
        self.ax1.plot(self.x, self.y1, 'r-')
        self.ax1.set_title('Синус')
        self.ax1.set_xlabel('X')
        self.ax1.set_ylabel('sin(x)')

        # График 2
        self.ax2.plot(self.x, self.y2, 'g--')
        self.ax2.set_title('Косинус')
        self.ax2.set_xlabel('X')
        self.ax2.set_ylabel('cos(x)')

        # График 3
        self.ax3.plot(self.x, self.y3, 'b-.')
        self.ax3.set_title('Тангенс')
        self.ax3.set_xlabel('X')
        self.ax3.set_ylabel('tan(x)')

        # График 4
        self.ax4.plot(self.x, self.y4, 'm:')
        self.ax4.set_title('Экспонента')
        self.ax4.set_xlabel('X')
        self.ax4.set_ylabel('exp(-x)*sin(x)')

        # Обновление графика
        self.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())