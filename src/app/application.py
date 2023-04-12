# System
import os, sys, threading

# GUI
import logging
from PyQt5 import uic 
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from app.logger import CustomFormatter, QPlainTextEditLogger

# GenBank functions
sys.path.append("../")
import genbank.tree, genbank.search, genbank.fetch, genbank.feature_parser


class Application(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(Application, self).__init__()

        # Multithreading attributes
        self.nb_max_threads = 8

        # File parsing attributes
        self.path = ""
        self.region_type = []
        self.nb_organisms_to_parse = 0
        self.nb_files_to_parse = 0
        self.nb_parsed_organisms = 0
        self.nb_parsed_files = 0

        # Load the ui file
        uic.loadUi("app/application.ui", self)

        # Define layout
        self.defineLayout()

        # Create widgets
        self.createWidgets()

        # Set-up logger
        self.setUpLogger()

        # Show the app
        self.show()
        
        # Update Results file tree
        genbank.tree.updateTree()

    def setUpLogger(self):
        """
        Set-up application logger.
        """
        self.log_file = "application.log"
        if os.path.exists(self.log_file):
            os.remove(self.log_file)
        logging.basicConfig(filename=self.log_file, encoding="utf-8", level=logging.DEBUG)
        self.logger_box = self.findChild(QFormLayout, "formLayout_6")
        logTextBox = QPlainTextEditLogger()
        self.logger_box.addWidget(logTextBox.widget)
        logging.getLogger().addHandler(logTextBox)
        logTextBox.setFormatter(CustomFormatter())
        logging.getLogger().setLevel(logging.DEBUG)

    def defineLayout(self):
        """
        Define Application layout.
        """
        grid = QGridLayout()
        self.setLayout(grid)
        self.splitter = self.findChild(QSplitter, "splitter")
        self.splitter.setStretchFactor(1, 10)

    def createWidgets(self):
        """
        Create widgets and assign a function to each widget.
        """
        # Model
        self.model = QFileSystemModel()
        self.model.setRootPath(QDir.currentPath())
        self.model.setFilter(QDir.NoDotAndDotDot | QDir.Dirs)

        # Push button
        self.button = self.findChild(QPushButton, "pushButton")
        self.button.clicked.connect(self.onButtonClicked)
        self.button_state = 0

        # Tree view
        self.treeView.setModel(self.model)
        self.treeView.setRootIndex(self.model.index(QDir.currentPath() + "/../Results"))
        for column in range(1, self.model.columnCount()):
            self.treeView.hideColumn(column)
        self.treeView.clicked.connect(self.onTreeViewClicked)

        # Check boxes
        self.checkBoxes=[self.CDS,self.CENTRO,self.INTRON,self.MOBILE,self.NC_RNA,self.R_RNA,self.TELOMETRE,self.T_RNA,self.UTR_3,self.UTR_5,self.OTHER]
        self.allChecked=False
        for k in range(len(self.checkBoxes)):
            self.checkBoxes[k].toggled.connect(self.onChecked)
        self.checkBoxes[0].setChecked(True)
        self.allChecked=False


    
    def progressbarAdvance(self):
        """
            change the progress bar 
        """
        self.progressBar_2.setValue(self.nb_parsed_organisms)
        self.progressBar_2.setMaximum(self.nb_organisms_to_parse+1)
        self.progressBar_2.setFormat("%v / %m")

        self.progressBar.setValue(self.nb_parsed_files)
        self.progressBar.setMaximum(self.nb_files_to_parse+1)
        self.progressBar.setFormat("%v / %m")

    def onButtonClicked(self):
        """
        Function to execute when the button is cliscked.
        """
        if self.button_state == 0:
            self.button_state = 1
            self.button.setText("Arrêter l'analyse")
            self.startParsing()
        else:
            self.button_state = 0
            self.button.setText("Lancer l'analyse")
            self.stopParsing()

    def onTreeViewClicked(self, index):
        """
        Function to execute when the tree view is clicked.

        Args:
            index (...): ???
        """
        temp_mod = index.model()
        self.path = temp_mod.filePath(index)
        logging.info("Select path: %s" % self.path)

    def onChecked(self):
        """
        Function to execute when a check box is clicked.
        """
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
        if(self.OTHER.isChecked()):
            for k in range(len(self.checkBoxes)):
                self.checkBoxes[k].setChecked(True)
            self.allChecked=True
        if(self.OTHER.isChecked()==False and self.allChecked):
            for k in range(len(self.checkBoxes)):
                self.checkBoxes[k].setChecked(False)
            self.allChecked=False

        logging.info("Selected DNA regions: " + str(self.region_type))

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
        return
    
    def singleThreadParsing(self, organisms):
        """
        Parse organisms sequentially (not using multithreading).

        Args:
            organisms (list): Tuples containing the name of the organism and the path to its folder.
        """
        for organism, organism_path in organisms:
            ids = genbank.search.searchID(organism)
            organism_files_to_parse = genbank.tree.needParsing(organism_path, ids)
            if organism_files_to_parse > 0:
                self.nb_organisms_to_parse += 1
                self.nb_files_to_parse += organism_files_to_parse
                for id in ids:
                    record = genbank.fetch.fetchFromID(id)
                    genbank.feature_parser.parseFeatures(self.region_type, organism_path, id, organism, record)
        


    def multiThreadParsing(self, organisms):
        """
        Parse organims sequentially using multithreading.

        Args:
            organisms (list): Tuples containing the name of the organism and the path to its folder.
        """
        parsing_attributes = []
        threads = []

        def threadFunction(organism_path, id, organism):
            logging.info("Start parsing file: %s" % id)
            record = genbank.fetch.fetchFromID(id)
            genbank.feature_parser.parseFeatures(self.region_type, organism_path, id, organism, record)

        for organism, organism_path in organisms:
            logging.info("Start parsing organism: %s" % organism)
            ids = genbank.search.searchID(organism)
            organism_files_to_parse = genbank.tree.needParsing(organism_path, ids)
            logging.info("Organism %s has %d file(s) that need(s) to be parsed" % (organism, organism_files_to_parse))
            if organism_files_to_parse > 0:
                self.nb_organisms_to_parse += 1
                self.nb_files_to_parse += organism_files_to_parse
                for id in ids:
                    parsing_attributes.append((organism_path, id, organism))

            # Create threads
            t = 0
            for attributes in parsing_attributes:
                t += 1
                logging.debug("Creating thread %d" % t)
                threads.append(threading.Thread(target=threadFunction, args=attributes))

            # Start threads
            for thread in threads:
                thread.start()

            # Wait for threads to finish
            for thread in threads:
                thread.join()
                #self.nb_files_to_parse -= 1
                self.nb_parsed_files+=1

            #self.nb_organisms_to_parse -= 1
            self.nb_parsed_organisms += 1
            self.progressbarAdvance()

        logging.info("Fin de l'analyse des fichiers sélectionnés")
        self.onButtonClicked()


    def startParsing(self):
        """
        Start file parsing.
        """
        logging.info("Initialising parsing")

        if self.path == "" or not os.path.isdir(self.path):
            logging.error("Invalid path: " + self.path)
            self.onButtonClicked()
            return
        else:
            logging.info("Start parsing")
            organisms = genbank.tree.findOrganisms(self.path)
            self.multiThreadParsing(organisms)

    def stopParsing(self):
        """
        Stop file parsing.
        """
        pass