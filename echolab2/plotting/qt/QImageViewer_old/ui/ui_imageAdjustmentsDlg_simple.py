# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'C:\Users\rick.towler\Work\AFSCGit\SurveyApps\MaceFunctions3\QImageViewer\ui\imageAdjustmentsDlg_simple.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_imageAdjustmentsDlg(object):
    def setupUi(self, imageAdjustmentsDlg):
        imageAdjustmentsDlg.setObjectName("imageAdjustmentsDlg")
        imageAdjustmentsDlg.resize(414, 500)
        imageAdjustmentsDlg.setMinimumSize(QtCore.QSize(316, 500))
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(imageAdjustmentsDlg)
        self.verticalLayout_2.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_2.setSpacing(4)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.gbBrightnessContrast = QtWidgets.QGroupBox(imageAdjustmentsDlg)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.gbBrightnessContrast.setFont(font)
        self.gbBrightnessContrast.setCheckable(True)
        self.gbBrightnessContrast.setChecked(False)
        self.gbBrightnessContrast.setObjectName("gbBrightnessContrast")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.gbBrightnessContrast)
        self.verticalLayout_4.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_4.setSpacing(3)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.label_8 = QtWidgets.QLabel(self.gbBrightnessContrast)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_8.setFont(font)
        self.label_8.setAlignment(QtCore.Qt.AlignCenter)
        self.label_8.setObjectName("label_8")
        self.verticalLayout_4.addWidget(self.label_8)
        self.brightnessSlider = QtWidgets.QSlider(self.gbBrightnessContrast)
        self.brightnessSlider.setMinimum(-100)
        self.brightnessSlider.setMaximum(100)
        self.brightnessSlider.setOrientation(QtCore.Qt.Horizontal)
        self.brightnessSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.brightnessSlider.setTickInterval(50)
        self.brightnessSlider.setObjectName("brightnessSlider")
        self.verticalLayout_4.addWidget(self.brightnessSlider)
        self.bcAutomatic = QtWidgets.QRadioButton(self.gbBrightnessContrast)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.bcAutomatic.setFont(font)
        self.bcAutomatic.setChecked(False)
        self.bcAutomatic.setObjectName("bcAutomatic")
        self.verticalLayout_4.addWidget(self.bcAutomatic)
        self.gbAutoBC = QtWidgets.QGroupBox(self.gbBrightnessContrast)
        self.gbAutoBC.setTitle("")
        self.gbAutoBC.setObjectName("gbAutoBC")
        self.formLayout_4 = QtWidgets.QFormLayout(self.gbAutoBC)
        self.formLayout_4.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout_4.setObjectName("formLayout_4")
        self.label_10 = QtWidgets.QLabel(self.gbAutoBC)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label_10)
        self.bcClipLimit = QtWidgets.QSlider(self.gbAutoBC)
        self.bcClipLimit.setMinimum(5)
        self.bcClipLimit.setMaximum(50)
        self.bcClipLimit.setProperty("value", 17)
        self.bcClipLimit.setOrientation(QtCore.Qt.Horizontal)
        self.bcClipLimit.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.bcClipLimit.setTickInterval(5)
        self.bcClipLimit.setObjectName("bcClipLimit")
        self.formLayout_4.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.bcClipLimit)
        self.verticalLayout_4.addWidget(self.gbAutoBC)
        self.bcManual = QtWidgets.QRadioButton(self.gbBrightnessContrast)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.bcManual.setFont(font)
        self.bcManual.setChecked(True)
        self.bcManual.setObjectName("bcManual")
        self.verticalLayout_4.addWidget(self.bcManual)
        self.gbManualBC = QtWidgets.QGroupBox(self.gbBrightnessContrast)
        self.gbManualBC.setTitle("")
        self.gbManualBC.setObjectName("gbManualBC")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.gbManualBC)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.label_7 = QtWidgets.QLabel(self.gbManualBC)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_7.setFont(font)
        self.label_7.setAlignment(QtCore.Qt.AlignCenter)
        self.label_7.setObjectName("label_7")
        self.verticalLayout_6.addWidget(self.label_7)
        self.contrastSlider = QtWidgets.QSlider(self.gbManualBC)
        self.contrastSlider.setMinimum(-100)
        self.contrastSlider.setMaximum(100)
        self.contrastSlider.setProperty("value", 0)
        self.contrastSlider.setSliderPosition(0)
        self.contrastSlider.setOrientation(QtCore.Qt.Horizontal)
        self.contrastSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.contrastSlider.setTickInterval(50)
        self.contrastSlider.setObjectName("contrastSlider")
        self.verticalLayout_6.addWidget(self.contrastSlider)
        self.verticalLayout_5.addLayout(self.verticalLayout_6)
        self.verticalLayout_4.addWidget(self.gbManualBC)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.pbBCReset = QtWidgets.QPushButton(self.gbBrightnessContrast)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.pbBCReset.setFont(font)
        self.pbBCReset.setObjectName("pbBCReset")
        self.horizontalLayout.addWidget(self.pbBCReset)
        self.verticalLayout_4.addLayout(self.horizontalLayout)
        self.verticalLayout_2.addWidget(self.gbBrightnessContrast)
        self.gbColorCorrection = QtWidgets.QGroupBox(imageAdjustmentsDlg)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(True)
        font.setWeight(75)
        self.gbColorCorrection.setFont(font)
        self.gbColorCorrection.setCheckable(True)
        self.gbColorCorrection.setChecked(False)
        self.gbColorCorrection.setObjectName("gbColorCorrection")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.gbColorCorrection)
        self.verticalLayout_3.setContentsMargins(5, 5, 5, 5)
        self.verticalLayout_3.setSpacing(3)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.cbAWB = QtWidgets.QCheckBox(self.gbColorCorrection)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cbAWB.setFont(font)
        self.cbAWB.setObjectName("cbAWB")
        self.verticalLayout_3.addWidget(self.cbAWB)
        self.ccManual = QtWidgets.QRadioButton(self.gbColorCorrection)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.ccManual.setFont(font)
        self.ccManual.setChecked(True)
        self.ccManual.setAutoExclusive(True)
        self.ccManual.setObjectName("ccManual")
        self.verticalLayout_3.addWidget(self.ccManual)
        self.gbManualCC = QtWidgets.QGroupBox(self.gbColorCorrection)
        self.gbManualCC.setTitle("")
        self.gbManualCC.setObjectName("gbManualCC")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.gbManualCC)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.formLayout.setObjectName("formLayout")
        self.label = QtWidgets.QLabel(self.gbManualCC)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.label)
        self.redSlider = QtWidgets.QSlider(self.gbManualCC)
        self.redSlider.setMinimum(-50)
        self.redSlider.setMaximum(50)
        self.redSlider.setOrientation(QtCore.Qt.Horizontal)
        self.redSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.redSlider.setTickInterval(10)
        self.redSlider.setObjectName("redSlider")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.redSlider)
        self.label_2 = QtWidgets.QLabel(self.gbManualCC)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.label_2)
        self.greenSlider = QtWidgets.QSlider(self.gbManualCC)
        self.greenSlider.setMinimum(-50)
        self.greenSlider.setMaximum(50)
        self.greenSlider.setOrientation(QtCore.Qt.Horizontal)
        self.greenSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.greenSlider.setTickInterval(10)
        self.greenSlider.setObjectName("greenSlider")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.greenSlider)
        self.label_3 = QtWidgets.QLabel(self.gbManualCC)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.label_3)
        self.blueSlider = QtWidgets.QSlider(self.gbManualCC)
        self.blueSlider.setMinimum(-50)
        self.blueSlider.setMaximum(50)
        self.blueSlider.setOrientation(QtCore.Qt.Horizontal)
        self.blueSlider.setTickPosition(QtWidgets.QSlider.TicksAbove)
        self.blueSlider.setTickInterval(10)
        self.blueSlider.setObjectName("blueSlider")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.blueSlider)
        self.verticalLayout_7.addLayout(self.formLayout)
        self.verticalLayout_3.addWidget(self.gbManualCC)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.pbColorReset = QtWidgets.QPushButton(self.gbColorCorrection)
        font = QtGui.QFont()
        font.setPointSize(10)
        font.setBold(False)
        font.setWeight(50)
        self.pbColorReset.setFont(font)
        self.pbColorReset.setObjectName("pbColorReset")
        self.horizontalLayout_3.addWidget(self.pbColorReset)
        self.verticalLayout_3.addLayout(self.horizontalLayout_3)
        self.verticalLayout_2.addWidget(self.gbColorCorrection)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem2)
        self.pbApply = QtWidgets.QPushButton(imageAdjustmentsDlg)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pbApply.setFont(font)
        self.pbApply.setObjectName("pbApply")
        self.horizontalLayout_2.addWidget(self.pbApply)
        self.pbCancel = QtWidgets.QPushButton(imageAdjustmentsDlg)
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pbCancel.setFont(font)
        self.pbCancel.setObjectName("pbCancel")
        self.horizontalLayout_2.addWidget(self.pbCancel)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)

        self.retranslateUi(imageAdjustmentsDlg)
        QtCore.QMetaObject.connectSlotsByName(imageAdjustmentsDlg)

    def retranslateUi(self, imageAdjustmentsDlg):
        _translate = QtCore.QCoreApplication.translate
        imageAdjustmentsDlg.setWindowTitle(_translate("imageAdjustmentsDlg", "Image Adjustments:"))
        self.gbBrightnessContrast.setTitle(_translate("imageAdjustmentsDlg", "Brightness and Contrast"))
        self.label_8.setText(_translate("imageAdjustmentsDlg", "Brightness"))
        self.bcAutomatic.setText(_translate("imageAdjustmentsDlg", "Automatic Contrast (CLAHE)"))
        self.label_10.setText(_translate("imageAdjustmentsDlg", "Clip Limit"))
        self.bcManual.setText(_translate("imageAdjustmentsDlg", "Manual Contrast"))
        self.label_7.setText(_translate("imageAdjustmentsDlg", "Contrast"))
        self.pbBCReset.setText(_translate("imageAdjustmentsDlg", "Reset"))
        self.gbColorCorrection.setTitle(_translate("imageAdjustmentsDlg", "Color Correction"))
        self.cbAWB.setText(_translate("imageAdjustmentsDlg", "Automatic White Balance"))
        self.ccManual.setText(_translate("imageAdjustmentsDlg", "Manual"))
        self.label.setText(_translate("imageAdjustmentsDlg", "R"))
        self.label_2.setText(_translate("imageAdjustmentsDlg", "G"))
        self.label_3.setText(_translate("imageAdjustmentsDlg", "B"))
        self.pbColorReset.setText(_translate("imageAdjustmentsDlg", "Reset"))
        self.pbApply.setText(_translate("imageAdjustmentsDlg", "Apply"))
        self.pbCancel.setText(_translate("imageAdjustmentsDlg", "Cancel"))