import os, sys
from PyQt5.QtWidgets import QApplication
from src import *

# os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # Create application
    _app = QApplication(sys.argv)
    window = Application()

    # Run application
    _app.exec_()