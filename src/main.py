import os, sys, logging
from PyQt5.QtWidgets import QApplication
from app.application import Application

# Test
import genbank.search, genbank.fetch, genbank.feature_parser

os.chdir(os.path.dirname(os.path.abspath(__file__)))

if __name__ == "__main__":
    # Set-up logging
    logging.basicConfig(filename="application.log", encoding="utf-8", level=logging.DEBUG)

    # Create application
    app = QApplication(sys.argv)
    window = Application()

    # Run application
    app.exec_()

    # Test
    id = ["NC_018416"] # For testing purpose, very small organism
    record = genbank.fetch.fetchFromID(id)
    genbank.feature_parser("", id, "ORGANISME_TEST", record)