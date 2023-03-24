from PyQt5.QtWidgets import *
from PyQt5.QtCore import *


class Widget(QWidget):
    def __init__(self, *args, **kwargs):

        QWidget.__init__(self, *args, **kwargs)

        hbox = QHBoxLayout(self)

        splitter1 = QSplitter(self)
        splitter1.setOrientation(Qt.Horizontal)

        center = QFrame(splitter1)
        center.setFrameShape(QFrame.StyledPanel)

        splitter2 = QSplitter(splitter1)
        sizePolicy = splitter2.sizePolicy()
        sizePolicy.setHorizontalStretch(1)

        splitter2.setSizePolicy(sizePolicy)
        splitter2.setOrientation(Qt.Vertical)

        top_right = QFrame(splitter2)
        top_right.setFrameShape(QFrame.StyledPanel)
        bottom_right = QFrame(splitter2)
        bottom_right.setFrameShape(QFrame.StyledPanel)

        hbox.addWidget(splitter1)
        self.setGeometry(10, 10, 2000, 1500)
        self.setWindowTitle("GENOME")


if __name__ == "__main__":
    import sys
    app = QApplication(sys.argv)
    app.setStyle("fusion")
    w = Widget()
    w.show()
    sys.exit(app.exec_())
