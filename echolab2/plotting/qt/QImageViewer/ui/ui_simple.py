# Form implementation generated from reading ui file 'C:\Users\rick.towler\Work\noaa-afsc-mace\MaceFunctions\QImageViewer\ui\simple.ui'
#
# Created by: PyQt6 UI code generator 6.6.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets


class Ui_simple(object):
    def setupUi(self, simple):
        simple.setObjectName("simple")
        simple.resize(919, 627)
        self.centralwidget = QtWidgets.QWidget(parent=simple)
        self.centralwidget.setObjectName("centralwidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.ImageViewLeft = QImageViewer(parent=self.centralwidget)
        self.ImageViewLeft.setObjectName("ImageViewLeft")
        self.horizontalLayout.addWidget(self.ImageViewLeft)
        self.ImageViewRight = QImageViewer(parent=self.centralwidget)
        self.ImageViewRight.setObjectName("ImageViewRight")
        self.horizontalLayout.addWidget(self.ImageViewRight)
        simple.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(parent=simple)
        self.statusbar.setObjectName("statusbar")
        simple.setStatusBar(self.statusbar)

        self.retranslateUi(simple)
        QtCore.QMetaObject.connectSlotsByName(simple)

    def retranslateUi(self, simple):
        _translate = QtCore.QCoreApplication.translate
        simple.setWindowTitle(_translate("simple", "Simple Zoom Example"))
from QImageViewer.QImageViewer import QImageViewer
