import sys
from PyQt4.QtCore import QRect
from PyQt4.QtGui import QApplication, QWidget, QTabWidget, QHBoxLayout

class Widget(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        # Create the layout.
        self.h_layout = QHBoxLayout()

        # Create the QTabWidget.
        self.tabWidget = QTabWidget()
        self.tabWidget.setEnabled(True)
        self.tabWidget.setGeometry(QRect(20, 40, 601, 501))
        self.tabWidget.setTabPosition(QTabWidget.North)
        self.tabWidget.setObjectName('tabWidget')

        # Add the QTabWidget to the created layout and set the 
        # layout on the QWidget.
        self.h_layout.addWidget(self.tabWidget)
        self.setLayout(self.h_layout)


if __name__ == '__main__':
  app = QApplication(sys.argv)

  widget = Widget()
  widget.show()

  sys.exit(app.exec_())