import sys
from app.window.ui.ui_main_window import UIMainWindow
from app.window.add_information_molecula import AddInformationMolecula


class MainWindow(UIMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.exit_program.triggered.connect(lambda: sys.exit(0))

        self.list_molecula_widget.add_molecula.clicked.connect(self.show_window_add_molecula)

    def show_window_add_molecula(self):
        _add_information_molecula = AddInformationMolecula(self)

        molecula = _add_information_molecula.show_modal()
