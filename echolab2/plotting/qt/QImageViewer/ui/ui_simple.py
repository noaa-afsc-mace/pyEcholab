# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'simple.ui'
##
## Created by: Qt User Interface Compiler version 6.11.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QHBoxLayout, QMainWindow, QSizePolicy,
    QStatusBar, QWidget)

from ..QImageViewer import QImageViewer

class Ui_simple(object):
    def setupUi(self, simple):
        if not simple.objectName():
            simple.setObjectName(u"simple")
        simple.resize(919, 627)
        self.centralwidget = QWidget(simple)
        self.centralwidget.setObjectName(u"centralwidget")
        self.horizontalLayout = QHBoxLayout(self.centralwidget)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.ImageViewLeft = QImageViewer(self.centralwidget)
        self.ImageViewLeft.setObjectName(u"ImageViewLeft")

        self.horizontalLayout.addWidget(self.ImageViewLeft)

        self.ImageViewRight = QImageViewer(self.centralwidget)
        self.ImageViewRight.setObjectName(u"ImageViewRight")

        self.horizontalLayout.addWidget(self.ImageViewRight)

        simple.setCentralWidget(self.centralwidget)
        self.statusbar = QStatusBar(simple)
        self.statusbar.setObjectName(u"statusbar")
        simple.setStatusBar(self.statusbar)

        self.retranslateUi(simple)

        QMetaObject.connectSlotsByName(simple)
    # setupUi

    def retranslateUi(self, simple):
        simple.setWindowTitle(QCoreApplication.translate("simple", u"Simple Zoom Example", None))
    # retranslateUi

