from PyQt6.QtWidgets import QVBoxLayout, QPushButton, QListView, QWidget


class UIListMoleculaView(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)

        base_layout = QVBoxLayout()
        self.setLayout(base_layout)

        self.add_molecula = QPushButton("Добавить молекулу")
        base_layout.addWidget(self.add_molecula)

        self.list_molecula_view = QListView()
        base_layout.addWidget(self.list_molecula_view)