# System
import os
from queue import Queue

# GUI
import logging
from PyQt5 import uic
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *

# Log and Thread
from .logger import CustomFormatter, QPlainTextEditLogger
from .parser_thread import *

# GenBank functions
from ..genbank import tree, fetch, feature_parser
from .logger import emitLog, Log


class Application(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(Application, self).__init__()

        # Multithreading attributes
        self.threadpool = QThreadPool.globalInstance()
        self.worker_queue = Queue()
        self.max_nb_threads = self.threadpool.maxThreadCount()
        self.nb_running_threads = 0

        # GUI attributes
        self.all_checked = False
        self.all_unchecked = True

        # File parsing attributes
        self.selected_path = ""
        self.region_type = []
        self.organisms_to_parse = []
        self.nb_organisms_to_parse = 0
        self.nb_files_to_parse = 0
        self.nb_parsed_organisms = 0
        self.nb_parsed_files = 0
        self.start_parsing_label = False

        # Load the ui file
        uic.loadUi("src/app/application.ui", self)

        # Define layout
        self.defineLayout()

        # Create widgets
        self.createWidgets()

        # Set-up logger
        self.setUpLogger()

        # Show the app
        self.show()

        # Update Results file tree
        tree.updateTree()

        emitLog(Log.INFO, "Multithreading with maximum %d threads" % self.threadpool.maxThreadCount())

    def signalHandler(self, frame):
            emitLog(Log.WARNING, "User interupt")
            sys.exit(42)

    def setUpLogger(self):
        """
        Set-up application logger.
        """
        self.log_file = "application.log"
        if os.path.exists(self.log_file):
            os.remove(self.log_file)
        logging.basicConfig(filename=self.log_file,
                            encoding="utf-8", level=logging.DEBUG)
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

        # Proxy Model to sort
        self.sort_proxy_model = QSortFilterProxyModel()
        self.sort_proxy_model.setSourceModel(self.model)
        self.sort_proxy_model.setDynamicSortFilter(True)
        self.sort_proxy_model.sort(0, Qt.AscendingOrder)

        # Push button
        self.button = self.findChild(QPushButton, "pushButton")
        self.button.clicked.connect(self.onButtonClicked)
        self.button_state = 0

        # Tree view
        self.treeView.setModel(self.sort_proxy_model)
        self.treeView.setRootIndex(self.sort_proxy_model.mapFromSource(
            self.model.index(os.path.join(QDir.currentPath(), "Results"))))
        for column in range(1, self.model.columnCount()):
            self.treeView.hideColumn(column)
        self.treeView.clicked.connect(self.onTreeViewClicked)

        # Check boxes
        self.checkboxes = [self.CDS, self.CENTRO, self.INTRON, self.MOBILE, self.NC_RNA,
                           self.R_RNA, self.TELOMETRE, self.T_RNA, self.UTR_3, self.UTR_5, self.ALL, self.NONE]
        self.all_checked = False
        self.none_checked = False
        for checkbox in self.checkboxes:
            checkbox.toggled.connect(self.onChecked)

        self.NONE.toggled.connect(self.onChecked_NONE)
        self.ALL.toggled.connect(self.onChecked_ALL)

        # Readme
        text_read = self.textEdit
        text_read.setReadOnly(True)
        with open('README.md', encoding='utf8') as f:
            markdown = f.read()
            text_read.setMarkdown(markdown)

    def updateProgressBar(self):
        """
        Change the progress bar.
        """
        self.fileProgressBar.setValue(self.nb_parsed_files)
        self.fileProgressBar.setMaximum(self.nb_files_to_parse)
        self.fileProgressBar.setFormat("%v / %m")

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
        mapped_index = self.sort_proxy_model.mapToSource(index)
        self.selected_path = self.model.filePath(mapped_index)
        emitLog(Log.INFO, "Selected path: %s" % self.selected_path)

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

        if(self.all_checked == False and self.none_checked == False):
            emitLog(Log.INFO, "Selected DNA regions: " + str(self.region_type))

    def onChecked_ALL(self):
        if(self.ALL.isChecked() and self.none_checked == False):
            self.all_checked = True
            for k in range(len(self.checkboxes)-2):
                self.checkboxes[k].setChecked(True)
            self.all_checked = False
            emitLog(Log.INFO, "Selected DNA regions: " + str(self.region_type))

    def onChecked_NONE(self):
        if(self.NONE.isChecked() and self.all_checked == False):
            self.none_checked = True
            for checkbox in self.checkboxes:
                checkbox.setChecked(False)
            self.none_checked = False
            emitLog(Log.INFO, "Selected DNA regions: " + str(self.region_type))

    def addWorker(self, worker):
        self.worker_queue.put(worker)
        self.processQueue()

    def processQueue(self):
        while not self.worker_queue.empty() and self.nb_running_threads < self.max_nb_threads - 1:
            worker = self.worker_queue.get()
            self.threadpool.start(worker)
            self.nb_running_threads += 1

    def threadWork(self, progress_callback, parsing_attribute, worker=None):
        organism_path, id, organism, worker = parsing_attribute
        emitLog(Log.INFO, "Start parsing file: %s" % id, worker)
        record = fetch.fetchFromID(id, worker=worker)
        if record is not None:
            feature_parser.parseFeatures(self.region_type, organism_path, id, organism, record, worker=worker)

    def threadLog(self, level, message):
        emitLog(level, message)

    def threadComplete(self):
        emitLog(Log.INFO, "Thread complete")
        self.nb_running_threads -= 1
        self.nb_parsed_files += 1
        self.processQueue()
        self.updateProgressBar()
        if self.worker_queue.empty() and self.nb_running_threads==0:
            logging.info("Fin de l'analyse des fichiers sélectionnés")
            self.button_state = 0
            self.button.setText("Lancer l'analyse")
            self.button.setEnabled(True)

    def threadPreparsing(self, parsing_attributes):
        emitLog(Log.INFO, "Starting workers to parse files...")
        t = 0
        self.nb_files_to_parse = len(parsing_attributes)
        for parsing_attribute in parsing_attributes:
            # Pass the function to execute
            # Any other args, kwargs are passed to the run function
            worker = Worker(self.threadWork, parsing_attribute=parsing_attribute)
            worker.signals.finished.connect(self.threadComplete)
            worker.signals.log.connect(self.threadLog)

            # Start the thread
            self.addWorker(worker)
            logging.info("Starting thread %d" % t)
            t += 1        

    def multiThreadParsing(self, organisms):
        """
        Parse organims sequentially using multithreading.

        Args:
            organisms (list): Tuples containing the name of the organism and the path to its folder.
        """
        emitLog(Log.INFO, "Start looking for files to parse...")
        pre_worker = Preworker(None, organisms=organisms)
        pre_worker.signals.result.connect(self.threadPreparsing)
        pre_worker.signals.log.connect(self.threadLog)
        self.threadpool.start(pre_worker)

    def startParsing(self):
        """
        Start file parsing.
        """
        emitLog(Log.INFO, "Initialising parsing")

        if self.selected_path == "" or not os.path.isdir(self.selected_path):
            emitLog(Log.ERROR, "Invalid path: " + self.selected_path)
            self.button_state = 0
            self.button.setText("Lancer l'analyse")
            return
        elif self.region_type == []:
            emitLog(Log.ERROR, "No DNA region selected")
            self.button_state = 0
            self.button.setText("Lancer l'analyse")
            return
        else:
            emitLog(Log.INFO, "Start parsing")
            logging.info("Start parsing")
            self.fileProgressBar.setValue(0)
            self.organisms_to_parse = tree.findOrganisms(self.selected_path)
            self.multiThreadParsing(self.organisms_to_parse)

    def stopParsing(self):
        """
        Stop file parsing.
        """
        logging.info("Stopping the parsing")
        if(self.nb_running_threads != 0):
            self.button.setEnabled(False)
            self.button.setText("Halting parsing")
        self.worker_queue = Queue()
        self.start_parsing_label = False

        logging.info("Stop parsing")
        return

    def signalHandler(self):
        """
        Delete all organisms files in the selected folder (self.selected_path).
        Note: Also remove files if they were created during a previous parsing. (For instance, files not needing update)
        """
        emitLog(Log.DEBUG, "Received SIGINT")
        if self.organisms_to_parse != [] and self.nb_files_to_parse - self.nb_parsed_files > 0:
            for organism, organism_path in self.organisms_to_parse:
                files = [file for file in os.listdir(organism_path) if os.path.isfile(file)]
                for file in files:
                    if os.path.exists(file):
                        os.remove(file)
                        emitLog(Log.INFO, "Successfully deleted: %s" % file)