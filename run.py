import sys
from PyQt6.QtWidgets import QApplication
from app.window import MainWindow
from rdkit import RDLogger
import traceback

if __name__ == "__main__":
    # Отключим вывод от rdkit
    RDLogger.DisableLog('rdApp.*')

    app = QApplication(sys.argv)

    try:
        main_window = MainWindow()
        main_window.show()
    except Exception as e:
        print(traceback.print_exc())
        print(e)

    sys.exit(app.exec())
