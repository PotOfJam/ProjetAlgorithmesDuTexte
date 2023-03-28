import os
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton
from app.application import Application


def action_btn():
    print("Ok")


os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Application()

    qBtn = QPushButton(window)
    qBtn.setText("Démarrer")
    qBtn.setGeometry(100, 100, 200, 30)
    qBtn.clicked.connect(action_btn)
    window.show()

    app.exec_()
