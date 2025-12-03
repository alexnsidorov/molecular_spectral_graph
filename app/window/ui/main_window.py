from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QMainWindow, QStyle, QSplitter, QHBoxLayout, QTabWidget, QWidget, QListView
from app.default_variable import NAME_PROGRAM
from app.window.ui.list_molecula_view import UIListMoleculaView

class UIMainWindow(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(NAME_PROGRAM)
        self.setWindowState(Qt.WindowState.WindowMaximized)

        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("&File")


        _base_layout = QHBoxLayout()

        _splitter = QSplitter(Qt.Orientation.Horizontal)
        _base_layout.addWidget(_splitter)

        _central_widget = QWidget()
        _central_widget.setLayout(_base_layout)
        self.setCentralWidget(_central_widget)

        close_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_TitleBarCloseButton)
        self.exit_program = file_menu.addAction(close_icon, "&Exit")
        self.load_data_action = menu_bar.addAction("Загрузить молекулы")
        # self.export_data_action = menu_bar.addAction("Выгрузить молекулы")

        self.list_molecula_widget = QListView(self)
        _splitter.addWidget(self.list_molecula_widget)

        self.tab_widget = QTabWidget(self)
        self.tab_widget.setTabsClosable(True)
        self.tab_widget.tabCloseRequested.connect(self._close_tab)
        _splitter.addWidget(self.tab_widget)

        _splitter.setStretchFactor(1, 1)
        _splitter.setStretchFactor(2, 2)

    def _close_tab(self, index):
        self.tab_widget.removeTab(index)
