from PyQt6.QtWidgets import QDialog, QHBoxLayout, QVBoxLayout, QPushButton, QSpinBox, QTextEdit, QLineEdit


class UIAddInformationMolecula(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setWindowTitle("Добвление молекулы")

        base_layout = QVBoxLayout()

        self.name_molecula = QLineEdit()
        self.name_molecula.setPlaceholderText("Названия молекулы")
        base_layout.addWidget(self.name_molecula)

        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText("SMILES")
        self.smiles.setToolTip("SMILES (Simplified Molecular Input Line Entry System) — система однозначной записи структуры молекул.")
        base_layout.addWidget(self.smiles)


        spectrum_layout = QHBoxLayout()

        self.ir_peaks = QTextEdit()
        self.ir_peaks.setPlaceholderText("Пики в инфракрасном спектре")
        self.ir_peaks.setToolTip("Пример ввода:\n700 0.8\n750 0.9\n1000 0.6")
        spectrum_layout.addWidget(self.ir_peaks)

        uv_layout = QVBoxLayout()
        self.width_peak = QSpinBox()
        self.width_peak.setSuffix(" нм")
        self.width_peak.setValue(15)
        self.width_peak.setMinimum(1)
        uv_layout.addWidget(self.width_peak)

        self.uv_peaks = QTextEdit()
        self.uv_peaks.setPlaceholderText("Пики поглощения в ультрафиолетовой области спектра")
        self.uv_peaks.setToolTip("Пример ввода:\n250 18000\n280 12000\n320 8000")
        uv_layout.addWidget(self.uv_peaks)

        spectrum_layout.addLayout(uv_layout)
        base_layout.addLayout(spectrum_layout)


        button_layout = QHBoxLayout()

        self.btn_accept = QPushButton("Добавить")
        button_layout.addWidget(self.btn_accept)
        self.btn_cancel = QPushButton("Отменить")
        button_layout.addWidget(self.btn_cancel)

        base_layout.addLayout(button_layout)

        self.setLayout(base_layout)

