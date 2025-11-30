from PyQt6.QtWidgets import QDialog, QHBoxLayout, QVBoxLayout, QPushButton, QSpinBox, QTextEdit, QLineEdit


class UIAddInformationMolecula(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        base_layout = QVBoxLayout()
        self.setWindowTitle("Добвление молекулы")

        self.name_molecula = QLineEdit()
        self.name_molecula.setPlaceholderText("Названия молекулы")
        base_layout.addWidget(self.name_molecula)

        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText("SMILES")
        self.smiles.setToolTip("SMILES (Simplified Molecular Input Line Entry System) — система однозначной записи структуры молекул.")
        base_layout.addWidget(self.smiles)

        self.width_peak = QSpinBox()
        self.width_peak.setSuffix(" нм")
        self.width_peak.setValue(15)
        self.width_peak.setMinimum(1)
        base_layout.addWidget(self.width_peak)

        button_layout = QHBoxLayout()

        self.btn_accept = QPushButton("Добавить")
        button_layout.addWidget(self.btn_accept)
        self.btn_cancel = QPushButton("Отменить")
        button_layout.addWidget(self.btn_cancel)

        base_layout.addLayout(button_layout)

        self.setLayout(base_layout)

