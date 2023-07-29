"""
Script to save point clicked on image to CSV file with pixel values for later analysis
author: Thomas Van Der Weide
"""

#####--------------------------------------------------------------------------------------------------------------------####
#########------------------------------------------pyQT Zoom -------------------------------------------------------########
#####--------------------------------------------------------------------------------------------------------------------####
#https://stackoverflow.com/questions/50379530/how-to-get-the-coordinate-of-the-loaded-image-and-not-the-one-from-the-display
from PyQt5.QtWidgets import QWidget, QApplication, QSlider, QGraphicsView, QGraphicsScene, QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QGraphicsPixmapItem, QShortcut
from PyQt5.QtGui import QPainter, QColor, QImage
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtOpenGL import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import pandas as pd
import os.path
from skimage import io
from datetime import datetime


class View(QGraphicsView):
    photo_clicked = QtCore.pyqtSignal(QtCore.QPoint)
    
    def __init__(self, parent):
        super(View, self).__init__()       
        self.new_Img()
                
    def new_Img(self):
        global i
        img = imgTarget_df['Img_Name'][i]
        self.image = io.imread(img)
        h,w,z = self.image.shape
        self.photo = QImage(self.image,self.image.shape[1],self.image.shape[0],self.image.strides[0],QImage.Format_RGB888)
        self.pixmap = QtGui.QPixmap()
        x_val = int(imgTarget_df['Img_X'][i])
        y_val = int(imgTarget_df['Img_Y'][i])
        if (x_val == 0) or (y_val == 0) or (x_val == 99999) or (y_val == 99999):
            pass
        else:
            self.photo.setPixelColor(x_val, y_val, QtGui.QColor(0, 0, 0))
        
        self.pixmap = self.pixmap.fromImage(self.photo)
        
        self.imgLabel = QtWidgets.QGraphicsPixmapItem()
        self.imgLabel.setPixmap(self.pixmap)
        
        self.begin = QtCore.QPoint()
        self.end = QtCore.QPoint()
        self.points = []
        
        self.scene = QtWidgets.QGraphicsScene(self)
        self.scene.addItem(self.imgLabel)
        self.setScene(self.scene)
        self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        
        self.resetTransform()
        self.scale(0.3, 0.3)
    
    def blank_Img(self):
        global i
        img = imgTarget_df['Img_Name'][i]
        self.image = io.imread(img)
        self.photo = QImage(self.image.shape[1],self.image.shape[0], QImage.Format_RGB888)
        self.pixmap = QtGui.QPixmap()

        self.pixmap = self.pixmap.fromImage(self.photo)
        
        self.imgLabel = QtWidgets.QGraphicsPixmapItem()
        self.imgLabel.setPixmap(self.pixmap)
        
        self.begin = QtCore.QPoint()
        self.end = QtCore.QPoint()
        self.points = []
        
        self.scene = QtWidgets.QGraphicsScene(self)
        self.scene.addItem(self.imgLabel)
        self.setScene(self.scene)
        self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
        
        self.resetTransform()
        self.scale(0.3, 0.3)
        
        
    def Hand_drag(self):
            self.setDragMode(QtWidgets.QGraphicsView.ScrollHandDrag)
            
    def pixel_pointer(self):
            self.setDragMode(QtWidgets.QGraphicsView.NoDrag)
            
    def mousePressEvent(self, event):
        if self.imgLabel.isUnderMouse():
            p = self.imgLabel.mapToItem(self.imgLabel, self.mapToScene(event.pos()))
            # Remove rounding errors
            p.setX(int(p.x()))
            p.setY(int(p.y()))
            self.pt = p.toPoint()
            # self.center = self.photo_clicked.emit(self.pt)
            if self.parent().btn_pix_info1.isChecked():
                self.center = self.photo_clicked.emit(self.pt)
                self.fillImage()
            else:
                pass
        super(View, self).mousePressEvent(event)
    
    def fillImage(self):
        i = self.pt.x()
        j = self.pt.y()
        print("Fill Img at %d,%d" %(i,j))
        # reload image to get rid of previous point clicked
        self.photo = QImage(self.image,self.image.shape[1],self.image.shape[0],self.image.strides[0],QImage.Format_RGB888)
        self.photo.setPixelColor(i,j, QtGui.QColor(0, 0, 0))
        self.pixmap = self.pixmap.fromImage(self.photo)
        self.imgLabel.setPixmap(self.pixmap)
        self.update()
        
    def show_prev(self):
        global i
        if i > 0:
            i = i -1
            self.new_Img()
        else:
            i = 0
            self.blank_Img()
        
    
    def show_next(self):
        global i
        global img_max
        if i < img_max:
            i = i + 1
            self.new_Img()
        else:
            i = i
            self.blank_Img()

        
    def wheelEvent(self, event):
        # Zoom Factor
        zoomInFactor = 1.25
        zoomOutFactor = 1 / zoomInFactor
    
        # Set Anchors
        self.setTransformationAnchor(View.NoAnchor)
        self.setResizeAnchor(View.NoAnchor)
    
        # Save the scene pos
        oldPos = self.mapToScene(event.pos())
    
        # Zoom
        if event.angleDelta().y() > 0:
            zoomFactor = zoomInFactor
        else:
            zoomFactor = zoomOutFactor
        self.scale(zoomFactor, zoomFactor)
    
        # Get the new position
        newPos = self.mapToScene(event.pos())
    
        # Move scene to old position
        delta = newPos - oldPos
        self.translate(delta.x(), delta.y())
        


class Window(QWidget):    
    def __init__(self):
        super(Window, self).__init__()
        self.view = View(self)
        self.center = None 
        
        self.btn_hand_drag = QtWidgets.QCheckBox("Hand drag", self)
        self.btn_hand_drag.clicked.connect(self.view.Hand_drag)
        self.btn_hand_drag.clicked.connect(self.btn_hand_drag_uncheck_others)

        self.btn_pix_info1 = QtWidgets.QCheckBox("Select Point", self)
        self.btn_pix_info1.clicked.connect(self.view.pixel_pointer)
        self.btn_pix_info1.clicked.connect(self.btn_pix_info1_drag_uncheck_other)
        self.editPixInfo1 = QtWidgets.QLineEdit(self)
        self.editPixInfo1.setReadOnly(True)
        
        self.prev = QtWidgets.QPushButton("Previous", self)
        self.prev.clicked.connect(self.view.show_prev)
        self.prev.clicked.connect(self.set_interval)
        self.prev.clicked.connect(self.set_gcpLabel)
        self.prev.clicked.connect(self.view.Hand_drag)
        self.prev.clicked.connect(self.btn_hand_drag_uncheck_others)
        self.next = QtWidgets.QPushButton("Next", self)
        self.next.clicked.connect(self.view.show_next)
        self.next.clicked.connect(self.set_interval)
        self.next.clicked.connect(self.set_gcpLabel)
        self.next.clicked.connect(self.btn_hand_drag_uncheck_others)
        self.next.clicked.connect(self.view.Hand_drag)
        
        self.intervals = QtWidgets.QLabel("Image Number 0", self)
        self.gcpLabel = QtWidgets.QLabel("gcp", self)
        self.btn_enter = QtWidgets.QPushButton("Enter", self)
        self.btn_enter.clicked.connect(self.save_TargetLoc)
        self.btn_skip = QtWidgets.QPushButton("Skip", self)
        self.btn_skip.clicked.connect(self.set_TargetSkip)

        self.view.photo_clicked.connect(self.photo_clicked)
        
        #Hotkeys
        # Space => Skip
        # Enter and T => Enter
        # S => Pointer for point selection
        # D => Hand for movement
        # A => Prev
        # F => Next
        
        # self.setFocusPolicy(Qt.StrongFocus)
        self._skip_key = QShortcut(QKeySequence(Qt.Key_Space), self)
        self._skip_key.activated.connect(self.btn_skip.click)
        
        self._enter_key = QShortcut(QKeySequence(Qt.Key_Return), self)
        self._enter_key.activated.connect(self.btn_enter.click)
        self._enter_key2 = QShortcut(QKeySequence(Qt.Key_T), self)
        self._enter_key2.activated.connect(self.btn_enter.click)
        
        self._select_key = QShortcut(QKeySequence(Qt.Key_S), self)
        self._select_key.activated.connect(self.btn_pix_info1.click)
        
        self._hand_key = QShortcut(QKeySequence(Qt.Key_D), self)
        self._hand_key.activated.connect(self.btn_hand_drag.click)
        
        self._prev_key = QShortcut(QKeySequence(Qt.Key_A), self)
        self._prev_key.activated.connect(self.prev.click)
        
        self._next_key = QShortcut(QKeySequence(Qt.Key_F), self)
        self._next_key.activated.connect(self.next.click)
        
        layout1 = QHBoxLayout()
        vbox = QVBoxLayout()
        layout3 = QVBoxLayout()
        
        vbox.addWidget(self.view)
        
        layout1.addLayout(vbox)
        layout3.addStretch()
        layout3.addWidget(self.btn_hand_drag)
        layout3.addWidget(self.btn_pix_info1)
        layout3.addWidget(self.editPixInfo1)
        layout3.addWidget(self.intervals)
        layout3.addWidget(self.prev)
        layout3.addWidget(self.next)
        layout3.addWidget(self.gcpLabel)
        layout3.addWidget(self.btn_enter)
        layout3.addWidget(self.btn_skip)
        layout3.addStretch()
        layout1.addLayout(layout3)
        
        self.setLayout(layout1)
        self.setWindowTitle("Image viewer")
        self.showMaximized()
        
        
    def btn_hand_drag_uncheck_others(self):
        self.btn_pix_info1.setChecked(False)

    def btn_pix_info1_drag_uncheck_other(self):
        self.btn_hand_drag.setChecked(False)
        
    def photo_clicked(self, pos):
        print(pos)
        global i
        self.editPixInfo1.setText('%d, %d' % (pos.x(), pos.y()))
        if self.btn_pix_info1.isChecked():
            imgTarget_df.loc[(imgTarget_df['Img_Name']==imgTarget_df['Img_Name'][i]) & (imgTarget_df['Marker_Name']==imgTarget_df['Marker_Name'][i]), ['Img_X']] = pos.x()
            imgTarget_df.loc[(imgTarget_df['Img_Name']==imgTarget_df['Img_Name'][i]) & (imgTarget_df['Marker_Name']==imgTarget_df['Marker_Name'][i]), ['Img_Y']] = pos.y()
            
        
    def set_interval(self):
        global i
        global img_max
        x_valtemp = int(imgTarget_df['Img_X'][i])
        y_valtemp = int(imgTarget_df['Img_Y'][i])
        self.intervals.setText('Image Number %d of %d' % (i, img_max))
        self.editPixInfo1.setText('%d, %d' % (x_valtemp, y_valtemp))
        
    def set_gcpLabel(self):
        global i
        self.gcpLabel.setText(imgTarget_df['Marker_Name'][i])
        
    def save_TargetLoc(self):
        global i
        x_val = int(imgTarget_df['Img_X'][i])
        y_val = int(imgTarget_df['Img_Y'][i])
        
        if (x_val == 0) or (y_val == 0) or (x_val == 99999) or (y_val == 99999):
            # Create textbox
            msg = QMessageBox()
            msg.setWindowTitle("Warning")
            msg.setIcon(QMessageBox.Critical)
            msg.setText("Are you sure (x,y) should be (%d, %d)?" % (x_val,y_val))
            msg.exec_()
        else:
            self.next.click()
            self.btn_hand_drag.click()
            
    def set_TargetSkip(self):
        global i
        imgTarget_df.loc[(imgTarget_df['Img_Name']==imgTarget_df['Img_Name'][i]) & (imgTarget_df['Marker_Name']==imgTarget_df['Marker_Name'][i]), ['Img_X']] = 99999
        imgTarget_df.loc[(imgTarget_df['Img_Name']==imgTarget_df['Img_Name'][i]) & (imgTarget_df['Marker_Name']==imgTarget_df['Marker_Name'][i]), ['Img_Y']] = 99999
        self.next.click()
        self.btn_hand_drag.click()
        


######-------------------------------MAIN-----------------------------------------######
if __name__== "__main__":
    # Enter Field Site, Camera type, and Date
    mainFold = "F:/SnowDrones/"
    fieldSiteList = ["LDP","BogusRidge","BullTrout","Headwall","PoleCat","TableRock","Treeline"]
    fieldSite = fieldSiteList[1]
    imgTypeList = ["RGB","Multispec","Thermal"]
    imgType = imgTypeList[0]
    locFold = mainFold + fieldSite + "/"
    
    AllDates = [dI for dI in sorted(os.listdir(locFold)) if os.path.isdir(os.path.join(locFold,dI))]
    ProcessDate = AllDates[3]
    # ProcessDate = ["02-04-2020"]
    print(ProcessDate)
    
    fn = locFold + ProcessDate + "/" + imgType + "/"
    # Load the Dataframe
    df_fn = fn + "imgTargets_df.feather"
    if os.path.isfile(df_fn):
        imgTarget_df = pd.read_feather(df_fn)
    else:
        print(df_fn)
        print('does not exist.')
    
    # Image File Path
    i = 0
    img_max = int(len(imgTarget_df.index) -1)
    
    if len(imgTarget_df) > 1:
        Cur_Date =  datetime.now().strftime("%m-%d")
        
        # Run the Visual
        app = QApplication.instance()
        if app is None:
                app = QApplication([])
        w = Window()
        w.show()
        w.raise_()
        app.exec_()
        
        save_df = imgTarget_df.copy()
        
        print(save_df.to_string())
        # Drop any columns that don't have entries
        save_df = save_df[(save_df != 0).all(1)]
        save_df = save_df[(save_df != 99999).all(1)]
        
        save_df = save_df.reset_index()
        save_df = save_df.drop(['index'], axis=1)
        
        # Rename Images so they match the names in Agisoft            
        save_dfTempName = save_df['Img_Name'].str.rpartition("/")[2]
        save_df['Img_Name'] = save_dfTempName.str.rpartition(".")[0]
        
        
        print(save_df.to_string())
        
        # Save as comma delim CSV file without index
        save_Fold = fn + "GCPwithImageLocations_11-24.csv"
        save_df.to_csv(save_Fold, encoding='utf-8', index=False)
    
    print("Finished Running Script.")


