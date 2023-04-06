# GUI
import logging
from PyQt5 import uic 
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# GenBank functions
import os, sys
sys.path.append("../")
import genbank.tree, genbank.search, genbank.fetch, genbank.feature_parser

class Application(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(Application, self).__init__()

        # Update Results file tree
        genbank.tree.updateTree()

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
        self.model = QFileSystemModel()
        self.model.setRootPath(QDir.currentPath())
        self.model.setFilter(QDir.NoDotAndDotDot | QDir.Dirs)
        self.treeView.setModel(self.model)
        self.treeView.setRootIndex(self.model.index(QDir.currentPath() + "/../Results"))
        for column in range(1, self.model.columnCount()):
            self.treeView.hideColumn(column)

        self.path = ""
        self.treeView.clicked.connect(self.on_treeView_clicked)

        self.region_type = []

        # Assign fonction
        self.asignWidgetsToFunction()

        # Show the app
        self.show()


    def asignWidgetsToFunction(self):
        """
        Assign a function to each widget.
        """
        # Push button
        self.startbutton = self.findChild(QPushButton, "pushButton")
        self.startbutton.clicked.connect(self.startParsing)


    def on_treeView_clicked(self, index):
        temp_mod = index.model()
        self.path = temp_mod.filePath(index)
        print(self.path)


    def test(self):

        if os.path.exists("CDS_ORGANISME_TEST_NC_000021.txt"):
            os.remove("CDS_ORGANISME_TEST_NC_000021.txt")
        if os.path.exists("intron_ORGANISME_TEST_NC_000021.txt"):
            os.remove("intron_ORGANISME_TEST_NC_000021.txt")
        print("DEBUT DU TEST")
        region_type = ["CDS", "intron"]
        # id = "NC_018416" # For testing purpose, very small organism
        id = "NC_000021" # For testing purpose, very small organism
        record = genbank.fetch.fetchFromID(id)
        genbank.feature_parser.parseFeatures(region_type, "", id, "ORGANISME_TEST", record)
        print("FIN DU TEST")
        return# id = "NC_018416" # For testing purpose, very small organism
    

    def startParsing(self):

        logging.info("Start parsing")

        if self.path == "" or not os.path.isdir(self.path):
            logging.error("Invalid path")
            return
        else:
            organisms = genbank.tree.findOrganisms(self.path)
            print("ORGANISMES:", organisms)

            for organism in organisms:
                ids = genbank.search.searchID(organism)
                for id in ids:
                    record = genbank.fetch.fetchFromID(id)
                    genbank.feature_parser.parseFeatures(self.region_type, , id, organism, record)
