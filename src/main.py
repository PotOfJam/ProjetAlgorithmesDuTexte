import os, sys, logging
from PyQt5.QtWidgets import QApplication
from app.application import Application

# Test
import genbank.search, genbank.fetch, genbank.feature_parser

os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # Set-up logging
    log_file = "application.log"
    if os.path.exists(log_file):
        os.remove(log_file)
    logging.basicConfig(filename=log_file, encoding="utf-8", level=logging.DEBUG)

    # Create application
    app = QApplication(sys.argv)
    window = Application()

    # Test
    id = ["NC_018416"] # For testing purpose, very small organism
    record = genbank.fetch.fetchFromID(id)
    genbank.feature_parser.parseFeatures("", id, "ORGANISME_TEST", record)
    
    # Run application
    app.exec_()