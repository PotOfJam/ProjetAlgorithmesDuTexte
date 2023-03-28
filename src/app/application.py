import sys

# Tree
sys.path.append("../")
from genbank.tree import updateTree

# GUI
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow

class Application(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(Application, self).__init__()

        # Update Results file tree
        updateTree()

        # Load the ui file
        uic.loadUi("app/application.ui", self)

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