import sys, traceback
from PyQt5.QtCore import QObject, QRunnable, pyqtSignal, pyqtSlot
from ..genbank.search import searchID
from ..genbank.tree import needParsing
from .logger import emitLog, Log

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
        super(Preworker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs["progress_callback"] = self.signals.progress
        self.organisms = self.kwargs["organisms"]

    @pyqtSlot()
    def run(self):
        """
        Initialise the runner function with passed args, kwargs.
        """
        try:
            parsing_attributes = []

            for organism, organism_path in self.organisms:
                emitLog(Log.INFO, "Start parsing organism: %s" % organism, self)
                ids = searchID(organism, worker=self)
                if ids == []:
                    emitLog(Log.WARNING, "Did not find any NC corresponding to organism: %s" % organism, self)
                    emitLog(Log.INFO, "Fin de l'analyse des fichiers sélectionnés", self)
                    break
                organism_files_to_parse = needParsing(organism_path, ids, worker=self)
                emitLog(Log.INFO, "Organism %s has %d file(s) that need(s) to be parsed" % (organism, organism_files_to_parse), self)
                if organism_files_to_parse > 0:
                    for id in ids:
                        parsing_attributes.append((organism_path, id, organism))
            self.signals.result.emit(parsing_attributes)
        except:
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))            