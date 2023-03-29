import os
import sys
import logging
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton
from app.application import Application


os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # Set-up logging
    logging.basicConfig(filename="application.log", encoding="utf-8", level=logging.DEBUG)

    # Create application
    app = QApplication(sys.argv)
    window = Application()

    # Run application
    app.exec_()
