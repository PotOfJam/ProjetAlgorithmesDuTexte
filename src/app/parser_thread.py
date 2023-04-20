import sys, traceback
from PyQt5.QtCore import QObject, QRunnable, pyqtSignal, pyqtSlot
from ..genbank.search import searchID
from ..genbank.tree import needParsing
from .logger import Log

class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.

    Supported signals are:

    finished
        No data

    error
        tuple (exctype, value, traceback.format_exc() )

    result
        object data returned from processing, anything

    progress
        int indicating % progress

    """
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(list)
    progress = pyqtSignal(int)
    log = pyqtSignal(Log, str)


class Worker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    """
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs["progress_callback"] = self.signals.progress
        self.kwargs["parsing_attribute"] = kwargs["parsing_attribute"] + (self,)

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """
        # Retrieve args/kwargs here and fire processing using them
        try:

            result = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        finally:
            self.signals.finished.emit()  # Done


class Preworker(QRunnable):
    """
    Worker thread

    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.

    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    """
    def __init__(self, fn, *args, **kwargs):
        super(Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs["progress_callback"] = self.signals.progress
        self.organisms = self.kwargs["Orga"][0]

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """
        # Retrieve args/kwargs here and fire processing using them
        try:
            parsing_attributes = []

            for organism, organism_path in self.organisms:
                #logging.info("Start parsing organism: %s" % organism)
                ids = searchID(organism)
                if ids == []:
                    #logging.warning("Did not find any NC corresponding to organism: %s" % organism)
                    #logging.info("Fin de l'analyse des fichiers sélectionnés")
                    break
                organism_files_to_parse = needParsing(organism_path, ids) # AMELIORABLE ?
                #logging.info("Organism %s has %d file(s) that need(s) to be parsed" % (organism, organism_files_to_parse))
                if organism_files_to_parse > 0:
                    self.nb_organisms_to_parse += 1
                    self.nb_files_to_parse += organism_files_to_parse
                    for id in ids:
                        parsing_attributes.append((organism_path, id, organism))
            
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        finally:
            self.signals.result.emit(parsing_attributes)  # Done