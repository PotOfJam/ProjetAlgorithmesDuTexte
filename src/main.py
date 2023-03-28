import os
import sys
from PyQt5.QtWidgets import QApplication
from app.application import Application

os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = Application()
    app.exec_()