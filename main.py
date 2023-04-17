import sys
from PyQt5.QtWidgets import QApplication
from src import *

if __name__ == "__main__":
    # Create application
    _app = QApplication(sys.argv)
    window = Application()

    # Run application
    sys.exit(_app.exec_())