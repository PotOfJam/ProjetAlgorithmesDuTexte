# GUI
import logging
from PyQt5 import uic 
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

class CustomFormatter(logging.Formatter):
    FORMATS = {
        logging.ERROR:   ("[%(levelname)-8s] %(message)s", QColor("red")),
        logging.DEBUG:   ("[%(levelname)-8s] [%(filename)s:%(lineno)d] %(message)s", "green"),
        logging.INFO:    ("[%(levelname)-8s] %(message)s", "#0000FF"),
        logging.WARNING: ('%(asctime)s - %(name)s - %(levelname)s - %(message)s', QColor(100, 100, 0))
    }

    def format( self, record ):
        last_fmt = self._style._fmt
        opt = CustomFormatter.FORMATS.get(record.levelno)
        if opt:
            fmt, color = opt
            self._style._fmt = "<font color=\"{}\">{}</font>".format(QColor(color).name(),fmt)
        res = logging.Formatter.format( self, record )
        self._style._fmt = last_fmt
        return res

class QPlainTextEditLogger(logging.Handler):
    def __init__(self, parent=None):
        super().__init__()
        self.widget = QPlainTextEdit(parent)
        self.widget.setReadOnly(True)    

    def emit(self, record):
        msg = self.format(record)
        self.widget.appendHtml(msg) 
        # move scrollbar
        scrollbar = self.widget.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())



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
        self.splitter.setStretchFactor(1, 10)
        

        self.logger_box = self.findChild(QFormLayout, "formLayout_6")
        logTextBox = QPlainTextEditLogger()
        self.logger_box.addWidget(logTextBox.widget)
        logging.getLogger().addHandler(logTextBox)
        logTextBox.setFormatter(CustomFormatter())
        logging.getLogger().setLevel(logging.DEBUG)

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

        # check boxes
        self.checkBoxes=[self.CDS,self.CENTRO,self.INTRON,self.MOBILE,self.NC_RNA,self.R_RNA,self.TELOMETRE,self.T_RNA,self.UTR_3,self.UTR_5,self.OTHER]
        self.checkBoxes[0].setChecked(True)

        for k in range(len(self.checkBoxes)):
            self.checkBoxes[k].toggled.connect(self.onChecked)

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

    def onChecked(self):

        self.region_type = []

        if(self.CDS.isChecked()):
            self.region_type.append("CDS")
        if(self.CENTRO.isChecked()):
            self.region_type.append("centromere")
        if(self.INTRON.isChecked()):
            self.region_type.append("intron")
        if(self.MOBILE.isChecked()):
            self.region_type.append("mobile_element")
        if(self.NC_RNA.isChecked()):
            self.region_type.append("ncRNA")
        if(self.R_RNA.isChecked()):
            self.region_type.append("rRNA")
        if(self.TELOMETRE.isChecked()):
            self.region_type.append("telomere")
        if(self.T_RNA.isChecked()):
            self.region_type.append("tRNA")
        if(self.UTR_3.isChecked()):
            self.region_type.append("3'UTR")
        if(self.UTR_5.isChecked()):
            self.region_type.append("5'UTR")

        print("check: region type= ",self.region_type)


        

    def onTreeViewClicked(self, index):
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
        
        logging.info("Initialising parsing")

        if self.path == "" or not os.path.isdir(self.path):
            logging.error("Invalid path")
            return
        else:
            logging.info("Start parsing")
            organisms = genbank.tree.findOrganisms(self.path)

            for organism in organisms:
                ids = genbank.search.searchID(organism)
                for id in ids:
                    record = genbank.fetch.fetchFromID(id)

                    if(self.CDS.isChecked()):
                        t=1
                    
                    genbank.feature_parser.parseFeatures(self.region_type, self.path, id, organism, record)