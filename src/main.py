import os
import sys
from PyQt5 import uic
from PyQt5.QtWidgets import QMainWindow, QApplication, QLabel, QTabWidget, QLineEdit, QComboBox, QPushButton
from math import pi

os.chdir(os.path.dirname(os.path.abspath(__file__)))

class AA(QMainWindow):

    def __init__(self):
        """
        Initialise the application.
        """
        super(AA, self).__init__()

        # Load the ui file
        uic.loadUi("main.ui", self)

        # Define widgets
        self.defineWidgets()

        # Assign fonction
        self.asignWidgetsToFunction()

        #
        self.computed_once = False

        # Show the app
        self.show()

    def defineWidgets(self):
        """
        Defines all widget used in the application.
        """
        # # Tabs
        # self.tabs = self.findChild(QTabWidget, "tabs")

        # # Density
        # self.density_input = self.findChild(QLineEdit, "density_input")
        # self.density_unit = self.findChild(QComboBox, "density_unit")

        # # Mass flow
        # self.mass_flow_input = self.findChild(QLineEdit, "mass_flow_input")
        # self.mass_flow_unit = self.findChild(QComboBox, "mass_flow_unit")

        # # Exterior diameter
        # self.ext_diameter_input = self.findChild(QLineEdit, "ext_diameter_input")
        # self.ext_diameter_unit = self.findChild(QComboBox, "ext_diameter_unit")

        # # Thickness
        # self.thickness_input = self.findChild(QLineEdit, "thickness_input")
        # self.thickness_unit = self.findChild(QComboBox, "thickness_unit")

        # # Result
        # self.result_output = self.findChild(QLineEdit, "result_output")
        # self.result_unit = self.findChild(QComboBox, "result_unit")

        # # Compute
        # self.compute_button = self.findChild(QPushButton, "compute")

        # # Reset
        # self.reset_button = self.findChild(QPushButton, "reset")

    
    def asignWidgetsToFunction(self):
        """
        Assign a function to each widget.
        """

        # # Compute
        # self.compute_button.clicked.connect(self.startComputing)

        # # Reset
        # self.reset_button.clicked.connect(self.resetComputeUI)

    
    def startComputing(self):
        """
        Start the computation.
        """
        # if self.allInputsFilled():
        #     if not self.computed_once:
        #         self.computeSpeed()
        #     else:
        #         if self.density_input.isModified() or self.mass_flow_input.isModified() or self.ext_diameter_input.isModified() or self.thickness_input.isModified():
        #             self.computeSpeed()


    def allInputsFilled(self):
        """
        Checks if all inputs fields are filled.
        Returns:
            bool: True if all inputs are filled.
        """
        # if self.density_input.text() and self.mass_flow_input.text() and self.ext_diameter_input.text() and self.thickness_input.text():
        #     return True
        # return False

    
    def allInputsCorrect(self, density, mass_flow, ext_diameter, thickness):
        """
        Check if all inputs are correct.
        Args:
            density (float): density of the fluid (kg/m3).
            mass_flow (float): mass flow (kg/h).
            ext_diameter (float): outer diameter of the pipe (m).
            thickness (float): thickness of the pipe (m).
        Returns:
            _type_: _description_
        """
        # if density <= 0 or mass_flow <= 0 or ext_diameter <= 0 or thickness <= 0:
        #     return False
        # return True

    def resetComputeUI(self):
        """
        Clear all inputs and outputs in the compute tab.
        """

        # Reset all inputs
        self.density_input.clear()
        self.mass_flow_input.clear()
        self.ext_diameter_input.clear()
        self.thickness_input.clear()

        # Reset output
        self.result_output.clear()

if __name__ == "__main__":
    # Initialise the app
    app = QApplication(sys.argv)
    window = AA()
    app.exec_()