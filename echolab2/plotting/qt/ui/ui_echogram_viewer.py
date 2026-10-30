# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'echogram_viewer.ui'
##
## Created by: Qt User Interface Compiler version 6.11.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QApplication, QMainWindow, QSizePolicy, QStatusBar,
    QVBoxLayout, QWidget)

from ..QImageViewer.QEchogramViewer import QEchogramViewer

class Ui_echogram_viewer(object):
    def setupUi(self, echogram_viewer):
        if not echogram_viewer.objectName():
            echogram_viewer.setObjectName(u"echogram_viewer")
        echogram_viewer.resize(1000, 600)
        echogram_viewer.setMinimumSize(QSize(0, 0))
        echogram_viewer.setMaximumSize(QSize(9999, 9999))
        self.actionLoad = QAction(echogram_viewer)
        self.actionLoad.setObjectName(u"actionLoad")
        self.actionExit = QAction(echogram_viewer)
        self.actionExit.setObjectName(u"actionExit")
        self.actionSet_EK60_Raw_data_path = QAction(echogram_viewer)
        self.actionSet_EK60_Raw_data_path.setObjectName(u"actionSet_EK60_Raw_data_path")
        self.centralwidget = QWidget(echogram_viewer)
        self.centralwidget.setObjectName(u"centralwidget")
        self.verticalLayout_4 = QVBoxLayout(self.centralwidget)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.QEchogramViewer = QEchogramViewer(self.centralwidget)
        self.QEchogramViewer.setObjectName(u"QEchogramViewer")
        self.QEchogramViewer.setMinimumSize(QSize(400, 400))

        self.verticalLayout_4.addWidget(self.QEchogramViewer)

        echogram_viewer.setCentralWidget(self.centralwidget)
        self.statusBar = QStatusBar(echogram_viewer)
        self.statusBar.setObjectName(u"statusBar")
        echogram_viewer.setStatusBar(self.statusBar)

        self.retranslateUi(echogram_viewer)

        QMetaObject.connectSlotsByName(echogram_viewer)
    # setupUi

    def retranslateUi(self, echogram_viewer):
        echogram_viewer.setWindowTitle(QCoreApplication.translate("echogram_viewer", u"Qt Echogram Viewer", None))
        self.actionLoad.setText(QCoreApplication.translate("echogram_viewer", u"Load Deployment...", None))
        self.actionExit.setText(QCoreApplication.translate("echogram_viewer", u"Exit", None))
        self.actionSet_EK60_Raw_data_path.setText(QCoreApplication.translate("echogram_viewer", u"Set Echogram properties", None))
    # retranslateUi

