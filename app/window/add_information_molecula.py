
from app.window.ui.add_infornation_molecula import UIAddInformationMolecula

class AddInformationMolecula(UIAddInformationMolecula):
    def __init__(self, parent = None):
        super().__init__(parent)
        self.btn_accept.clicked.connect(self.accept)
        self.btn_cancel.clicked.connect(self.reject)

    def show_modal(self):
        if self.exec():
            print("Что-то добавили")
        else:
            pass