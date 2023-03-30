import os
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QFileSystemModel
from app.application import Application
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *



os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # Set-up logging
    #logging.basicConfig(filename="application.log", encoding="utf-8", level=logging.DEBUG)

    # Create application
    app = QApplication(sys.argv)
    window = Application()

    # Run application
    app.exec_()
