# -*- coding: utf-8 -*-

#import PySide2
#from PySide2.QtWidgets import QApplication, QMessageBox

import csv
import sys
import PySide2
import adme_predict
import os
from PySide2 import QtWidgets
from PySide2.QtUiTools import QUiLoader
from PySide2.QtWidgets import QApplication,QMainWindow,QFileDialog,QTextEdit,QAction

from ui_adme import Ui_Form
from adme_predict.run import predict
from PySide2 import QtXml
#from PySide2.QtGui import QIcon, QFont
#import os
#QtCore, QtGui,


class Wdl_guiForm(QtWidgets.QWidget,Ui_Form):
    def __init__(self,parent=None):
        super(Wdl_guiForm,self).__init__(parent)
        self.setupUi(self)
        self.initUI()
        self.ui = QUiLoader().load('admee.ui')

    def initUI(self):

#       self.radioButton.setIcon(QIcon('D:\PY\\admet\首页.ICO'))
#        self.radioButton.clicked.connect(self.showDialog1)
#        self.lineEdit_2.clicked.connect(self.predicting)

#        self.toolButton_2.setIcon(QIcon('D:\PY\\admet\诊断.ICO'))
        self.toolButton.clicked.connect(self.showDialog2)

        self.toolButton_2.clicked.connect(self.showDialog3)

        self.pushButton.clicked.connect(self.predicting)

        self.toolButton_3.clicked.connect(self.generate)

        self.toolButton_4.clicked.connect(self.showDialog4)

#    def showDialog1(self,):
#        filename,filetype= (QFileDialog.getOpenFileName(self,'Select Path ','/'))
#        self.lineEdit_2.setText(filename)

    def showDialog2(self):
        filename,filetype = (QFileDialog.getOpenFileName(self, 'Select Path ', ' ', '*.csv'))
        self.lineEdit.setText(filename)

    def showDialog3(self):
        filename = (QFileDialog.getExistingDirectory(self,'Select Path ','/'))
        self.lineEdit_3.setText(filename)

    def showDialog4(self):
        filename = (QFileDialog.getExistingDirectory(self,'Select Path ','/'))
        self.lineEdit_4.setText(filename)


    def predicting(self):

        filename = self.lineEdit.text()
        #print(filename)
        savename = self.lineEdit_3.text() +'\\result.csv'
        print(savename)
        filepath = os.path.abspath(__file__)
        splitdir = filepath.split('\\')
        condaDir = splitdir[:-1]
        if condaDir[0] == 'D:':
		    condaDir[0] = 'D:\\'
        if condaDir[0] == 'C:':
		    condaDir[0] = 'C:\\'
        #print(condaDir)
        a= os.path.join(*condaDir)
        current_path=a+'\\PKL\\'
        #print(current_path)
        predict(filename, savename, current_path)

    def generate(self):
        smis = self.lineEdit_2.text()
        #print(smis)
        savename=self.lineEdit_4.text() +'\smiles.csv'
        print(savename)
        f = open(savename,'w')
        csv_writer = csv.writer(f)
        csv_writer.writerow(["Smiles"])
        csv_writer.writerow([smis])
        f.close()

if __name__ =='__main__':
    app = QApplication(sys.argv)
    mainw = Wdl_guiForm()
    mainw.show()
    sys.exit(app.exec_())


