# GUI
from PyQt5 import uic 
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# GenBank functions
import sys
sys.path.append("../")
from genbank.tree import updateTree
import genbank.search, genbank.fetch, genbank.feature_parser

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

        # Define layout
        grid = QGridLayout()
        self.setLayout(grid)
        self.splitter = self.findChild(QSplitter, "splitter")
        print(self.splitter)
        self.splitter.setStretchFactor(1, 10)

        # treeView.setHeaderHidden(True)
        # grid.addWidget(treeView, 3, 3)
        # self.setCentralWidget(treeView)

        # treeView.resize(100, 100)

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
        model.setFilter(QDir.NoDotAndDotDot | QDir.Dirs)
        self.treeView.setModel(model)
        self.treeView.setRootIndex(model.index(QDir.currentPath()+"/../Results"))
        for column in range(1, model.columnCount()):
            self.treeView.hideColumn(column)

        #self.tabs.insertTopLevelItems(None, tree)

    def asignWidgetsToFunction(self):
        """
        Assign a function to each widget.
        """
        # Push button
        self.startbutton = self.findChild(QPushButton, "pushButton")
        self.startbutton.clicked.connect(self.test)

    def test(self):
        # Test
        id = "NC_018416" # For testing purpose, very small organism
        record = genbank.fetch.fetchFromID(id)
        genbank.feature_parser.parseFeatures("", id, "ORGANISME_TEST", record)
        print("FIN DU TEST")
        return