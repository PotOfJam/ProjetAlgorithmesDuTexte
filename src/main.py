import os
import sys

# GUI
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QApplication, QLabel, QTabWidget, QLineEdit, QComboBox, QPushButton

# Tree
from genbank.tree import updateTree

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class AA(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(AA, self).__init__()

        # Load the ui file
        uic.loadUi("main.ui", self)

        # Define widgets
        self.defineWidgets()

        # Assign fonction
        self.asignWidgetsToFunction()

        # Show the app
        self.show()

    def defineWidgets(self):
        """
        Defines all widget used in the application.
        """
        pass
    
    def asignWidgetsToFunction(self):
        """
        Assign a function to each widget.
        """
        pass

if __name__ == "__main__":
    # Initialise the app
    updateTree()
    app = QApplication(sys.argv)
    window = AA()
    app.exec_()