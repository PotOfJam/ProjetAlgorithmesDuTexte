import os, sys, logging
from PyQt5.QtWidgets import QApplication
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
