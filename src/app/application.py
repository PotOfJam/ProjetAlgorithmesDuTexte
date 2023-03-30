from PyQt5.QtWidgets import *
from PyQt5.Qt import QStandardItemModel, QStandardItem
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import uic
from PyQt5 import QtWidgets
# from genbank.tree import updateTree
import sys

# Tree
sys.path.append("../")

# GUI
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class Application(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(Application, self).__init__()

        # Update Results file tree
        # updateTree()

        # Load the ui file
        uic.loadUi("app/application.ui", self)

        grid = QGridLayout()
        self.setLayout(grid)

        # treeView.setHeaderHidden(True)
        #grid.addWidget(treeView, 3, 3)
        # self.setCentralWidget(treeView)

        #treeView.resize(100, 100)

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
        model = QFileSystemModel()
        model.setRootPath(QDir.currentPath())
        self.treeView.setModel(model)
        self.treeView.setRootIndex(model.index(QDir.currentPath()+"/../Results"))

        #self.tabs.insertTopLevelItems(None, tree)

    def asignWidgetsToFunction(self):
        """
        Assign a function to each widget.
        """
        pass
