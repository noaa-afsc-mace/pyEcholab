# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'imageAdjustmentsDlg.ui'
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
from PySide6.QtWidgets import (QApplication, QCheckBox, QDialog, QFormLayout,
    QGroupBox, QHBoxLayout, QLabel, QPushButton,
    QRadioButton, QSizePolicy, QSlider, QSpacerItem,
    QVBoxLayout, QWidget)

class Ui_imageAdjustmentsDlg(object):
    def setupUi(self, imageAdjustmentsDlg):
        if not imageAdjustmentsDlg.objectName():
            imageAdjustmentsDlg.setObjectName(u"imageAdjustmentsDlg")
        imageAdjustmentsDlg.resize(414, 500)
        imageAdjustmentsDlg.setMinimumSize(QSize(316, 500))
        self.verticalLayout_2 = QVBoxLayout(imageAdjustmentsDlg)
        self.verticalLayout_2.setSpacing(4)
        self.verticalLayout_2.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.gbBrightnessContrast = QGroupBox(imageAdjustmentsDlg)
        self.gbBrightnessContrast.setObjectName(u"gbBrightnessContrast")
        font = QFont()
        font.setPointSize(10)
        font.setBold(True)
        self.gbBrightnessContrast.setFont(font)
        self.gbBrightnessContrast.setCheckable(True)
        self.gbBrightnessContrast.setChecked(False)
        self.verticalLayout_4 = QVBoxLayout(self.gbBrightnessContrast)
        self.verticalLayout_4.setSpacing(3)
        self.verticalLayout_4.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_4.setObjectName(u"verticalLayout_4")
        self.label_8 = QLabel(self.gbBrightnessContrast)
        self.label_8.setObjectName(u"label_8")
        font1 = QFont()
        font1.setPointSize(10)
        font1.setBold(False)
        self.label_8.setFont(font1)
        self.label_8.setAlignment(Qt.AlignCenter)

        self.verticalLayout_4.addWidget(self.label_8)

        self.brightnessSlider = QSlider(self.gbBrightnessContrast)
        self.brightnessSlider.setObjectName(u"brightnessSlider")
        self.brightnessSlider.setMinimum(-100)
        self.brightnessSlider.setMaximum(100)
        self.brightnessSlider.setOrientation(Qt.Horizontal)
        self.brightnessSlider.setTickPosition(QSlider.TicksAbove)
        self.brightnessSlider.setTickInterval(50)

        self.verticalLayout_4.addWidget(self.brightnessSlider)

        self.bcAutomatic = QRadioButton(self.gbBrightnessContrast)
        self.bcAutomatic.setObjectName(u"bcAutomatic")
        font2 = QFont()
        font2.setBold(False)
        self.bcAutomatic.setFont(font2)
        self.bcAutomatic.setChecked(False)

        self.verticalLayout_4.addWidget(self.bcAutomatic)

        self.gbAutoBC = QGroupBox(self.gbBrightnessContrast)
        self.gbAutoBC.setObjectName(u"gbAutoBC")
        self.formLayout_4 = QFormLayout(self.gbAutoBC)
        self.formLayout_4.setObjectName(u"formLayout_4")
        self.formLayout_4.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.label_10 = QLabel(self.gbAutoBC)
        self.label_10.setObjectName(u"label_10")
        self.label_10.setFont(font2)

        self.formLayout_4.setWidget(0, QFormLayout.ItemRole.LabelRole, self.label_10)

        self.bcClipLimit = QSlider(self.gbAutoBC)
        self.bcClipLimit.setObjectName(u"bcClipLimit")
        self.bcClipLimit.setMinimum(5)
        self.bcClipLimit.setMaximum(50)
        self.bcClipLimit.setValue(17)
        self.bcClipLimit.setOrientation(Qt.Horizontal)
        self.bcClipLimit.setTickPosition(QSlider.TicksAbove)
        self.bcClipLimit.setTickInterval(5)

        self.formLayout_4.setWidget(0, QFormLayout.ItemRole.FieldRole, self.bcClipLimit)


        self.verticalLayout_4.addWidget(self.gbAutoBC)

        self.bcManual = QRadioButton(self.gbBrightnessContrast)
        self.bcManual.setObjectName(u"bcManual")
        self.bcManual.setFont(font2)
        self.bcManual.setChecked(True)

        self.verticalLayout_4.addWidget(self.bcManual)

        self.gbManualBC = QGroupBox(self.gbBrightnessContrast)
        self.gbManualBC.setObjectName(u"gbManualBC")
        self.verticalLayout_5 = QVBoxLayout(self.gbManualBC)
        self.verticalLayout_5.setObjectName(u"verticalLayout_5")
        self.verticalLayout_6 = QVBoxLayout()
        self.verticalLayout_6.setObjectName(u"verticalLayout_6")
        self.label_7 = QLabel(self.gbManualBC)
        self.label_7.setObjectName(u"label_7")
        self.label_7.setFont(font1)
        self.label_7.setAlignment(Qt.AlignCenter)

        self.verticalLayout_6.addWidget(self.label_7)

        self.contrastSlider = QSlider(self.gbManualBC)
        self.contrastSlider.setObjectName(u"contrastSlider")
        self.contrastSlider.setMinimum(-100)
        self.contrastSlider.setMaximum(100)
        self.contrastSlider.setValue(0)
        self.contrastSlider.setSliderPosition(0)
        self.contrastSlider.setOrientation(Qt.Horizontal)
        self.contrastSlider.setTickPosition(QSlider.TicksAbove)
        self.contrastSlider.setTickInterval(50)

        self.verticalLayout_6.addWidget(self.contrastSlider)


        self.verticalLayout_5.addLayout(self.verticalLayout_6)


        self.verticalLayout_4.addWidget(self.gbManualBC)

        self.horizontalLayout = QHBoxLayout()
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.horizontalSpacer_3 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout.addItem(self.horizontalSpacer_3)

        self.pbBCReset = QPushButton(self.gbBrightnessContrast)
        self.pbBCReset.setObjectName(u"pbBCReset")
        self.pbBCReset.setFont(font1)

        self.horizontalLayout.addWidget(self.pbBCReset)


        self.verticalLayout_4.addLayout(self.horizontalLayout)


        self.verticalLayout_2.addWidget(self.gbBrightnessContrast)

        self.gbColorCorrection = QGroupBox(imageAdjustmentsDlg)
        self.gbColorCorrection.setObjectName(u"gbColorCorrection")
        self.gbColorCorrection.setFont(font)
        self.gbColorCorrection.setCheckable(True)
        self.gbColorCorrection.setChecked(False)
        self.verticalLayout_3 = QVBoxLayout(self.gbColorCorrection)
        self.verticalLayout_3.setSpacing(3)
        self.verticalLayout_3.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_3.setObjectName(u"verticalLayout_3")
        self.cbAWB = QCheckBox(self.gbColorCorrection)
        self.cbAWB.setObjectName(u"cbAWB")
        self.cbAWB.setFont(font2)

        self.verticalLayout_3.addWidget(self.cbAWB)

        self.ccManual = QRadioButton(self.gbColorCorrection)
        self.ccManual.setObjectName(u"ccManual")
        self.ccManual.setFont(font2)
        self.ccManual.setChecked(True)
        self.ccManual.setAutoExclusive(True)

        self.verticalLayout_3.addWidget(self.ccManual)

        self.gbManualCC = QGroupBox(self.gbColorCorrection)
        self.gbManualCC.setObjectName(u"gbManualCC")
        self.verticalLayout_7 = QVBoxLayout(self.gbManualCC)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.formLayout = QFormLayout()
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
        self.label = QLabel(self.gbManualCC)
        self.label.setObjectName(u"label")
        self.label.setFont(font1)

        self.formLayout.setWidget(0, QFormLayout.ItemRole.LabelRole, self.label)

        self.redSlider = QSlider(self.gbManualCC)
        self.redSlider.setObjectName(u"redSlider")
        self.redSlider.setMinimum(-50)
        self.redSlider.setMaximum(50)
        self.redSlider.setOrientation(Qt.Horizontal)
        self.redSlider.setTickPosition(QSlider.TicksAbove)
        self.redSlider.setTickInterval(10)

        self.formLayout.setWidget(0, QFormLayout.ItemRole.FieldRole, self.redSlider)

        self.label_2 = QLabel(self.gbManualCC)
        self.label_2.setObjectName(u"label_2")
        self.label_2.setFont(font1)

        self.formLayout.setWidget(1, QFormLayout.ItemRole.LabelRole, self.label_2)

        self.greenSlider = QSlider(self.gbManualCC)
        self.greenSlider.setObjectName(u"greenSlider")
        self.greenSlider.setMinimum(-50)
        self.greenSlider.setMaximum(50)
        self.greenSlider.setOrientation(Qt.Horizontal)
        self.greenSlider.setTickPosition(QSlider.TicksAbove)
        self.greenSlider.setTickInterval(10)

        self.formLayout.setWidget(1, QFormLayout.ItemRole.FieldRole, self.greenSlider)

        self.label_3 = QLabel(self.gbManualCC)
        self.label_3.setObjectName(u"label_3")
        self.label_3.setFont(font1)

        self.formLayout.setWidget(2, QFormLayout.ItemRole.LabelRole, self.label_3)

        self.blueSlider = QSlider(self.gbManualCC)
        self.blueSlider.setObjectName(u"blueSlider")
        self.blueSlider.setMinimum(-50)
        self.blueSlider.setMaximum(50)
        self.blueSlider.setOrientation(Qt.Horizontal)
        self.blueSlider.setTickPosition(QSlider.TicksAbove)
        self.blueSlider.setTickInterval(10)

        self.formLayout.setWidget(2, QFormLayout.ItemRole.FieldRole, self.blueSlider)


        self.verticalLayout_7.addLayout(self.formLayout)


        self.verticalLayout_3.addWidget(self.gbManualCC)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalSpacer_2 = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_3.addItem(self.horizontalSpacer_2)

        self.pbColorReset = QPushButton(self.gbColorCorrection)
        self.pbColorReset.setObjectName(u"pbColorReset")
        self.pbColorReset.setFont(font1)

        self.horizontalLayout_3.addWidget(self.pbColorReset)


        self.verticalLayout_3.addLayout(self.horizontalLayout_3)


        self.verticalLayout_2.addWidget(self.gbColorCorrection)

        self.horizontalLayout_2 = QHBoxLayout()
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.horizontalSpacer = QSpacerItem(40, 20, QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)

        self.horizontalLayout_2.addItem(self.horizontalSpacer)

        self.pbApply = QPushButton(imageAdjustmentsDlg)
        self.pbApply.setObjectName(u"pbApply")
        font3 = QFont()
        font3.setPointSize(10)
        self.pbApply.setFont(font3)

        self.horizontalLayout_2.addWidget(self.pbApply)

        self.pbCancel = QPushButton(imageAdjustmentsDlg)
        self.pbCancel.setObjectName(u"pbCancel")
        self.pbCancel.setFont(font3)

        self.horizontalLayout_2.addWidget(self.pbCancel)


        self.verticalLayout_2.addLayout(self.horizontalLayout_2)


        self.retranslateUi(imageAdjustmentsDlg)

        QMetaObject.connectSlotsByName(imageAdjustmentsDlg)
    # setupUi

    def retranslateUi(self, imageAdjustmentsDlg):
        imageAdjustmentsDlg.setWindowTitle(QCoreApplication.translate("imageAdjustmentsDlg", u"Image Adjustments:", None))
        self.gbBrightnessContrast.setTitle(QCoreApplication.translate("imageAdjustmentsDlg", u"Brightness and Contrast", None))
        self.label_8.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Brightness", None))
        self.bcAutomatic.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Automatic Contrast (CLAHE)", None))
        self.gbAutoBC.setTitle("")
        self.label_10.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Clip Limit", None))
        self.bcManual.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Manual Contrast", None))
        self.gbManualBC.setTitle("")
        self.label_7.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Contrast", None))
        self.pbBCReset.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Reset", None))
        self.gbColorCorrection.setTitle(QCoreApplication.translate("imageAdjustmentsDlg", u"Color Correction", None))
        self.cbAWB.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Automatic White Balance", None))
        self.ccManual.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Manual", None))
        self.gbManualCC.setTitle("")
        self.label.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"R", None))
        self.label_2.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"G", None))
        self.label_3.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"B", None))
        self.pbColorReset.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Reset", None))
        self.pbApply.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Apply", None))
        self.pbCancel.setText(QCoreApplication.translate("imageAdjustmentsDlg", u"Cancel", None))
    # retranslateUi

