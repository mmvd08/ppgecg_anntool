from pyqtgraph.Qt import QtGui, QtCore
import io
import pathlib
import logging
import platform
import pyqtgraph as pg
import numpy as np
from os import getcwd
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QMenuBar, QMenu, QDockWidget
from PyQt5.QtWidgets import QAction
from PyQt5.QtCore import Qt
from pyqtgraph import DataTreeWidget
from pyqtgraph.dockarea import *
from mainwindow import Ui_MainWindow
from scipy.io.matlab import savemat, loadmat
from senales import Senal
from cursor import CursorInfo
from sklearn.cluster import DBSCAN
import scipy.signal as signal
from PyQt5.QtWidgets import QMenuBar, QMenu, QAction
from func import detectRealIdx
logger = logging.getLogger()
import path as camino
#/home/aasl/PycharmProjects/PPG
########################Funciones de aasl ###############################################################################
def ClickedLabel(ev):
    cInfo.setPoint(QtGui.QCursor.pos())
    signal_obj.setLastEvent(ev.opts['name'])
    cmToolMenu.popup(QtGui.QCursor.pos())


def Clicked(ev):
    signal_obj.setLastEvent(ev.opts['name'])
    if cInfo.getCross():
        cInfo.setCross(False)
    else:
        cInfo.setPoint(QtGui.QCursor.pos())
        cmToolMenu.popup(QtGui.QCursor.pos())
        menuC.popup(QtGui.QCursor.pos())


def ClickedArrhythmia(ev):
    cmToolMenu2.popup(QtGui.QCursor.pos())


def mouseMoved1(evt):
    flg = 0
    pos = evt[0]
    if(signal_obj.getLastEvent()=="Peak" or signal_obj.getLastEvent() == "Foot" or signal_obj.getLastEvent() == "Dichr"):
        s = signal_obj.getPPG()
        d2 = area.docks['PPG']
    else:
        s = signal_obj.getECG()
        d2 = area.docks['ECG']
    fs = signal_obj.getfs()
    p21 = d2.widgets[0]
    vb = p21.getViewBox()

    #l=p21.plotItem.items
    vLine = p21.plotItem.items[2]
    hLine = p21.plotItem.items[-1]

    rect = p21.sceneBoundingRect()

    if rect.contains(pos):
        mousePoint = vb.mapSceneToView(pos)
        idx = int(fs * mousePoint.x())
        vLine.setPos(mousePoint.x())
        hLine.setPos(s[idx])


def mouseClicked1(evt):

    if evt[0]._double:
        a = signal_obj.getLastEvent()

        if a == "Peak_R":
            d2 = area.docks['ECG']
        if a == "Peak" :
            d2 = area.docks['PPG']
        if a == "Foot":
            d2 = area.docks['PPG']
        if a == "Dichr":
            d2 = area.docks['PPG']
        p21 = d2.widgets[0]
        vb = p21.getViewBox()
        cInfo.setCross(False)
        cInfo.setPoint(vb.mapSceneToView(evt[0]._scenePos))

def detqrs( y, fs, fac=1.5):      #y=senal
    ########Esto no esparte del codgo original#############
     max1=np.max(y)
     min1=np.min(y)
     d=max1+ np.abs(min1)
     fac=0.3+d
    ######################################################
     yf = y
     yfd = np.diff(yf)
     nterm = np.max(np.abs(yfd))
     nterm2=np.median(np.abs(yfd))
     yfds = np.round(yfd * fac / nterm)
     yfdds = np.diff(yfds)
     k = np.nonzero(yfdds)[0]
     l = len(k)
     r = []
     for i in range(0, (l - 1)):
         if (k[i + 1] - k[i]) > round(0.1 * fs):
              r.append(k[i])
     N_lat = len(r)
     ra = np.array(r)
     rb = ra.copy()
     d = round(0.1 * fs)
     for i in range(1, N_lat - 1):
          a = y[(ra[i] - d):(ra[i] + d)]
          idx = np.argmax(np.diff(a))
          rb[i] = rb[i] + idx - d
     rc = np.unique(rb)
     rr = np.diff(rc)
     q = np.min(rr)
     s = np.argmin(rr)
     lr = list(rc)
     thr = np.mean(rr) - 2 * np.std(rr)
     if q < thr:
        lr.remove(lr[s])
     rd = np.array(lr)
     rdl = rd.size
     wns = np.intp(0.05 * fs)
     for i in range(1, rdl):
         wdw = y[(rd[i] - wns):(rd[i] + wns)]
         if wdw.shape[0] > 0:
            pp = np.argmax(wdw)
            rd[i] = rd[i] + pp - wns
         else:
            continue
     re = np.unique(rd)
     return re

def another0(p):
    if signal_obj.getLastEvent()=="Type":
        cmToolMenu2.setEnabled(True)
        v = p
        #print(p)

    else:
        txt = p.data()
        txtnum = float(txt)

        wdw_st = txtnum - 2.5
        wdw_ed = txtnum + 2.5

        d2 = area.docks['ECG']
        p3 = d2.widgets[0]

        lr3 = p3.plotItem.items[1]
        lr3.setRegion([wdw_st, wdw_ed])

#####################################################################################################################################################
def filtrado_ppg(senal,fs):           #Filtrado para la senal de PPG
    FNyquist1 = fs / 2
    Wn = np.array([0.5, 10])  # Para pasabandas o supresores de banda
    Wn11 = (2 / fs) * Wn
    b1, a1 = signal.butter(4, Wn11, 'bandpass')  # Filtro(orden,frecuencia de corte, tipo)
    #b1,a1=signal.cheby2(4,40,Wn11,'bandpass')
    senalF = signal.filtfilt(b1, a1, senal)
    return senalF


def filtrado_ecg(senal,fs):              #Filtrado pa la senal de ECG:no se usa porque la base de datos del computing ya esta filtrada
    FNyquist1 = fs / 2
    Wn = np.array([0.05, 40])  # Para pasabandas o supresores de banda
    Wn11 = (2 / fs) * Wn
    b1, a1 = signal.butter(4, Wn11, 'bandpass')  # Filtro(orden,frecuencia de corte, tipo)
    senalF = signal.filtfilt(b1, a1, senal)
    return senalF


def plotear(dock, senal,t,puntos, name):      #funcion para plotear puntos
    d1=area.docks['PPG']
    d2 = area.docks['ECG']
    p1 = d1.widgets[0]
    p1.plotItem.clear()
    p2 = d2.widgets[0]
    p2.plotItem.clear()
    ecg=signal_obj.getECG()
    ppg=signal_obj.getPPG()
    p1.plot(t, ppg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name='PPG')
    p2.plot(t, ecg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name='ECG')
    dx=area.docks[dock]
    px = dx.widgets[0]
    if (name=='Peak'):
        color='r'
    if (name=='Foot'):
        color= 'm'
    if (name=='Dichr'):
        color='b'
    if (dock=="PPG"):
        #p2.plot(t, senal, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name=dock)
        p=px.plot(t[puntos], senal[puntos], pen=None, symbolBrush=color, symbolPen='w', name=name)
        p.sigClicked.connect(ClickedLabel)
    else:
        #p2.plot(t, senal, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name=dock)
        e = px.plot(t[puntos], senal[puntos], pen=None, symbolBrush='r', symbolPen='w', name=name)
        e.sigClicked.connect(ClickedLabel)

def menu_arrhythmia(q):   #Funcion del menu para borrar arritmias almacenadas
    #cmToolMenu2.setEnabled(True)

    #str=None
    # if signal_obj.getLastEvent()== "Type":
    #     str = 'ArrhythmiaD'
    # if q.text() == "Remove arrhythmia":
    #     PosI=signal_obj.getPI()
    #     PosF=signal_obj.getPF()
    #     Tipo=signal_obj.getTipo()

    #quit(0)

    if q.text() == "Remove arrhythmia":
        if signal_obj.getArrhythmiaDialog() == True:
            ql01 = win.findChild(QtWidgets.QListWidget, name='ListaA')
            if ql01.currentItem() is not None:
                signal_obj.removeItemArr(ql01.row(ql01.currentItem()))
                ql01.clear()

                piTotal = signal_obj.getPI()
                pfTotal = signal_obj.getPF()
                tipTotal = signal_obj.getTipo()

                l = len(piTotal)
                if l > 0:
                    for i in range(l):
                        ql01.addItem(str(piTotal[i]) + " - " + str(pfTotal[i]) + " : " + tipTotal[i])
                else:
                    frame = win.findChild(QtWidgets.QFrame, name='ArrhFrame')
                    frame.deleteLater()
                    signal_obj.setArrhythmiaDialog(False)


    elif q.text() == "Sort arrhythmia":
        posI=signal_obj.getPI()
        posF=signal_obj.getPF()
        tipo=signal_obj.getTipo()

        signal_obj.clear_PI()
        signal_obj.clear_PF()
        signal_obj.clear_Tipo()

        l1=len(posI)
        for i in range(l1):
            minim=np.argmin(posI)
            signal_obj.setPI(posI[minim])
            signal_obj.setPF(posF[minim])
            signal_obj.setTipo(tipo[minim])

            posI[minim]=1000
            posF[minim]=1000
            tipo[minim]=1000
        ql01 = win.findChild(QtWidgets.QListWidget, name='ListaA')
        piTotal = signal_obj.getPI()
        pfTotal = signal_obj.getPF()
        tipTotal = signal_obj.getTipo()
        ql01.clear()
        for i in range(l1):
            ql01.addItem(str(piTotal[i]) + " - " + str(pfTotal[i]) + " : " + tipTotal[i])



        pass

        #if signal_obj.getArrhythmiaDialog() == True:
        #    ql01 = win.findChild(QtWidgets.QListWidget, name='ListaA')
        #    ql01.sortItems()


def adjust_points(q):       #Funcion del menu para corregir puntos
    fs=signal_obj.getfs()
    t=signal_obj.getTime()
    str = None
    if signal_obj.getLastEvent() == "Peak" or signal_obj.getLastEvent() == "Foot" or signal_obj.getLastEvent() == "Dichr":
        d2 = area.docks['PPG']
        s = signal_obj.getPPG()
        p=signal_obj.getP()
        f=signal_obj.getpie()
        d=signal_obj.getdicr()
        str = 'PPG'
    else:
        s = signal_obj.getECG()
        d2 = area.docks['ECG']
        pR=signal_obj.getR()
        str='ECG'

    if q.text() == "Adjust point":     # si se selecciona e el menu 'Adjust point'
        p21 = d2.widgets[0]
        vb = p21.getViewBox()
        ipoint = cInfo.getPoint() - win.pos()
        point = vb.mapSceneToView(ipoint)
        if signal_obj.getLastEvent() == "Peak":
            puntos=p
        if signal_obj.getLastEvent() == "Foot":
            puntos=f
        if signal_obj.getLastEvent() == "Dichr":
            puntos = d
        if signal_obj.getLastEvent() == "Peak_R":
            puntos = pR
        idx = int(point.x() * fs)
        isFirst, isLast, ind = detectRealIdx(puntos, idx)    #buscando la posicion del punto que se desea ajustar
        rl = list(puntos)
        rl.remove(rl[ind])      #eliminando el punto que se desea ajustar
        rm = np.array(np.intp(rl))
        #p21.plotItem.items[1].sigClicked.connect(Clicked)
        if(signal_obj.getplot_i()==True):
            color='m'
        if (signal_obj.getplot_p()==True):
            color='r'
        #p21hdl1 = p21.plotItem.items[1]
        #p21hdl1.sigClicked.disconnect(Clicked)
        plotear(str, s, t, rm, signal_obj.getLastEvent())
        #p21.removeItem(p21.plotItem.items[-1])
        #p21.plot(t[rm], s[rm], pen=None, symbolBrush=color, symbolPen='w', name=signal_obj.getLastEvent())

        vLine = pg.InfiniteLine(angle=90, movable=False)
        hLine = pg.InfiniteLine(angle=0, movable=False)
        p21.addItem(vLine, ignoreBounds=True)
        p21.addItem(hLine, ignoreBounds=True)
        proxy1 = pg.SignalProxy(p21.scene().sigMouseMoved, rateLimit=60, slot=mouseMoved1)
        proxy2 = pg.SignalProxy(p21.scene().sigMouseClicked, delay=0.5, rateLimit=1, slot=mouseClicked1)
        app.setOverrideCursor(QtGui.QCursor(QtCore.Qt.CrossCursor))
        cInfo.setCross(True)

        while cInfo.getCross():
            app.processEvents()

        app.restoreOverrideCursor()
        proxy1.disconnect()
        proxy2.disconnect()
        mypoint = cInfo.getPoint()
        idx2 = int(mypoint.x() * fs)
        rl.insert(0, idx2)
        rm1 = np.array(np.intp(rl))
        rm2 = np.sort(rm1)
        if signal_obj.getLastEvent() == "Peak":
            signal_obj.setP(rm2)
        if signal_obj.getLastEvent() == "Foot":
            signal_obj.setpie(rm2)
        if signal_obj.getLastEvent() == "Dichr":
            signal_obj.setdicr(rm2)
        if signal_obj.getLastEvent() == "Peak_R":
            signal_obj.setR(rm2)
        if signal_obj.getplot_i()==True:
            color='m'
        if signal_obj.getplot_p() == True:
            color = 'r'
        #p21.plotItem.removeItem(p21.plotItem.items[-1])
        #p21hdl1 = p21.plot(t[rm2], s[rm2], pen=None, symbolBrush=color, symbolPen='w', name=signal_obj.getLastEvent())
        #p21hdl1.sigClicked.connect(Clicked)
        plotear(str, s, t, rm2, signal_obj.getLastEvent())

    elif q.text() == "Remove point":
        p21 = d2.widgets[0]
        vb = p21.getViewBox()
        ipoint = cInfo.getPoint() - win.pos()
        # opoint = win.pos()
        point = vb.mapSceneToView(ipoint)
        idx = int(point.x() * fs)
        if signal_obj.getLastEvent() == "Peak":
            puntos=p
        if signal_obj.getLastEvent() == "Foot":
            puntos=f
        if signal_obj.getLastEvent() == "Dichr":
            puntos = d
        if signal_obj.getLastEvent() == "Peak_R":
            puntos = pR
        isFirst, isLast, ind = detectRealIdx(puntos, idx)
        rl = list(puntos)
        rl.remove(rl[ind])
        rm = np.array(np.intp(rl))
        if signal_obj.getLastEvent() == "Peak":
            signal_obj.setP(rm)
        if signal_obj.getLastEvent() == "Foot":
            signal_obj.setpie(rm)
        if signal_obj.getLastEvent() == "Dichr":
            signal_obj.setdicr(rm)
        if signal_obj.getLastEvent() == "Peak_R":
            signal_obj.setR(rm)
        if (signal_obj.getplot_i() == True):
            color = 'm'
        if (signal_obj.getplot_p() == True):
            color = 'r'

        p21hdl1 = p21.plotItem.items[1]
        #p21hdl1.sigClicked.disconnect(Clicked)
        # p21hdl1.sigPointsClicked.disconnect(ClickedPoints)
        #p21.removeItem(p21.plotItem.items[-1])
        #p21hdl1 = p21.plot(t[rm], s[rm], pen=None, symbolBrush=color, symbolPen='w', name=signal_obj.getLastEvent())
        #p21hdl1.sigClicked.connect(Clicked)
        # p21hdl1.sigPointsClicked.connect(ClickedPoints)
        plotear(str, s, t, rm, signal_obj.getLastEvent())
    elif q.text() == "Add point":
        p21 = d2.widgets[0]
        vb = p21.getViewBox()
        ipoint = cInfo.getPoint() - win.pos()
        point = vb.mapSceneToView(ipoint)
        idx = int(point.x() * fs)
        if signal_obj.getLastEvent() == "Peak":
            puntos = p
        if signal_obj.getLastEvent() == "Foot":
            puntos = f
        if signal_obj.getLastEvent() == "Dichr":
            puntos = d
        if signal_obj.getLastEvent() == "Peak_R":
            puntos = pR
        isFirst, isLast, ind = detectRealIdx(puntos, idx)

        if isFirst:
            newR = np.intp(puntos[0] / 2)
        else:
            newR = np.intp((puntos[ind - 1] + puntos[ind]) / 2)

        rl = list(puntos)
        rl.insert(ind, newR)
        rm = np.array(np.intp(rl))
        if signal_obj.getLastEvent() == "Peak":
            signal_obj.setP(rm)
        if signal_obj.getLastEvent() == "Foot":
            signal_obj.setpie(rm)
        if signal_obj.getLastEvent() == "Dichr":
            signal_obj.setdicr(rm)
        if signal_obj.getLastEvent() == "Peak_R":
            signal_obj.setR(rm)
            #p21hdl1 = p21.plotItem.items[1]
        #p21hdl1.sigClicked.disconnect(Clicked)
        # p21hdl1.sigPointsClicked.disconnect(ClickedPoints)
        if (signal_obj.getplot_i() == True):
            color = 'm'
        if (signal_obj.getplot_p() == True):
            color = 'r'
        plotear(str, s, t, rm, signal_obj.getLastEvent())
        #p21.removeItem(p21.plotItem.items[-1])
         #p21hdl1 = p21.plot(t[rm], s[rm], pen=None, symbolBrush=color, symbolPen='w', name=signal_obj.getLastEvent())
        #p21hdl1.sigClicked.connect(Clicked)


            #if signal_obj.getViewInfoPanell():
            #buildInfoDock(True)
def FileOpen():
    file_name = QtGui.QFileDialog.getOpenFileName(None,
                                                  directory='/home/aasl/SharedXP/Arianna/Maestria/codigo /Bases de datos',
                                                  caption='Open  file...',
                                                  filter="Matlab files (*.mat)")
    v = file_name[0]
    flname = pathlib.PurePath(file_name[0])
    ext = flname.suffix
    pos = file_name[0].find(ext)
    full_header = str(file_name[0][0:pos] + ".hea")
    onlyValECG = 1
    onlyValPPG = 1

    if any(v):
        file_header = io.open(full_header, 'r')
        if file_header:
            m1 = file_header.readline()
            m2 = file_header.readline()
            m3 = file_header.readline()
            m4 = file_header.readline()

            mfs=m1.split(' ')
            fs=float(mfs[2])

            mList2 = m2.split(' ')

            pos_num2 = mList2[2].find('/')
            onlyValECG = mList2[2][0:pos_num2]
            unitsECG = mList2[2][pos_num2 + 1:-1]

            mList2 = m4.split(' ')

            pos_num2 = mList2[2].find('/')
            onlyValPPG = mList2[2][0:pos_num2]
            unitsPPG = mList2[2][pos_num2 + 1:-1]

    gainECG = float(onlyValECG)
    gainPPG = float(onlyValPPG)

    if any(v):
        arch = loadmat(file_name[0])
        datos = arch['val']
        ecg = np.array(datos[0, :] / gainECG)
        ppg = np.array(datos[2, :] / gainPPG)
        # fs = 250.0
        x = 1
        return x, ecg, ppg, fs
    else :
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(Exception).__name__, Exception.args)
        logger.debug(message)
        logger.error("Operation 'file opening' cancelled by user")
        x = 0
        return x, x, x, x


def Archivo(accion):
    if(accion.text()=="Open"):
        x,ecg,ppg,fs=FileOpen()
        if x==1:
            signal_obj.signal_clear()
            area.clear()
            d2 = Dock("ECG", size=(1000, 600))
            d3 = Dock("PPG", size=(1000, 600))

            area.addDock(d2, 'right')  ## place d2 at right edge of dock area
            area.addDock(d3, 'bottom', d2)  ## place d3 at right edge of dock area

            n = len(ecg)
            x = np.linspace(0, n / fs, n, endpoint=False)
            ecg = filtrado_ecg(ecg, fs)    #no es necesarip porque estas senales ya estan filtradas
            p2 = pg.PlotWidget()
            p2hdl1 = p2.plot(x, ecg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="ECG")
            p2.showGrid(x=True, y=True)
            # p2.plotItem.setMenuEnabled(False)

            ppg = filtrado_ppg(ppg, fs)
            p3 = pg.PlotWidget()
            p2hdl1 = p3.plot(x, ppg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="PPG")
            p3.showGrid(x=True, y=True)
            # p3.plotItem.setMenuEnabled(False)

            # vlayout = win.findChild(QtWidgets.QWidget, "centralWidget")

            d2.addWidget(p2, row=0, col=0)
            d3.addWidget(p3, row=1, col=0)

            p2.setXLink(p3)

            signal_obj.setECG(ecg)
            signal_obj.setPPG(ppg)
            signal_obj.setfs(fs)
            signal_obj.setTime(x)

            uimain.actionPeak_2.setEnabled(True)
            uimain.actionFoot_2.setEnabled(True)
            uimain.actionDichrotic_2.setEnabled(True)
            uimain.actionPeak_3.setEnabled(True)
            uimain.actionArrhythmia.setEnabled(True)
            uimain.actionDistance.setEnabled(True)
            uimain.actionHeart_rate.setEnabled(True)
    elif(accion.text()=="Save"):
        name, fltyp = QtGui.QFileDialog.getSaveFileName(caption='Guardar resultados', directory='/home/aasl/PycharmProjects/PPG', filter='MATLAB files(*.mat);; JSON files(*.json)')
        if name:
            signal_obj.save(name)
    else:
        QtGui.QApplication.instance().quit()

    # if (accion.text() == "Exit"):
    #     #quit(0)
    #     QtGui.QApplication.instance().quit(0)


def detect_peak():
    ppg = signal_obj.getPPG()
    fs = signal_obj.getfs()
    m = np.diff(ppg)  # buscando la primera derivada
    picos = list()

    for i in range(0, len(m) - 2):  # buscando los picos en la primera derivada(cambios de signo)
        pos1 = m[i]
        pos2 = m[i + 1]
        if pos2 < 0 and pos1 > 0:
            picos.append(i + 1)

    picos = np.array(picos, dtype=int)
    picobinario = np.zeros((len(m), 1), dtype=int)
    for i in range(0, len(picos) - 1):
        window = ppg[picos[i:i + 10]]
        amp_max = np.max(window)
        umbral = amp_max * 0.40  # estableciendo como umbral el 40% del maximo de la senal
        if ppg[picos[i]] > umbral:
            picobinario[picos[i]] = 1
    refract = np.ceil(fs * 0.25)
    picobinario = np.array(picobinario)

    ###### Algoritmo DBSCAN
    hw_size = refract
    detect = np.array(picobinario, dtype='int').squeeze()
    detect2 = np.array(np.nonzero(detect)).squeeze()
    detect3 = detect2.reshape(-1, 1)
    db = DBSCAN(eps=hw_size, min_samples=1,
                metric='cityblock')  # Agrupar los 1 que estan a una distancia de una muestra hasta 1/4 del periodo refractario
    clusters = db.fit_predict(detect3)
    top = np.max(clusters) + 1
    picosdbscan = list()
    for i in range(top - 1):
        idx1 = clusters == i
        idxint = np.array(idx1, dtype='int')
        ppg1 = idxint * ppg[detect2]
        val = np.argmax(ppg1)
        picosdbscan.append(detect2[val])

    picosdbscanx = np.array(picosdbscan)
    picos = picosdbscanx
    t = signal_obj.getTime()
    a='PPG'
    b=ppg
    c=t
    d=picos
    e='Peak'
    return a,b,c,d,e


def Deteccion(accion):
    cmToolMenu.setEnabled(True)
    if (accion.text() == "Peak R"):    #Para el ECG
        read_R=signal_obj.getR()
        try:
            if read_R==None:
                s = signal_obj.getECG()
                fs = signal_obj.getfs()
                t = signal_obj.getTime()
                p = detqrs(s, fs, fac=1.5)
                plotear('ECG', s, t, p, 'Peak_R')
                signal_obj.setR(p)
                signal_obj.setRD(p)
                uimain.actionSave.setEnabled(True)
        except ValueError:
            s = signal_obj.getECG()
            fs = signal_obj.getfs()
            t = signal_obj.getTime()
            p = read_R
            plotear('ECG', s, t, p, 'Peak_R')
            signal_obj.setR(p)
            uimain.actionSave.setEnabled(True)

    if(accion.text()=="Peak"):         #Para el PPG
        read_p=signal_obj.getP()
        uimain.actionSave.setEnabled(True)
        try:
            if read_p==None:      # Verificando si ya se han detectado los picos anteriormente
                a, b, c, d, e = detect_peak()  # a='PPG', b=ppg, c=t, d=picos, e='Peak'
                signal_obj.setplot_p(True)
                signal_obj.setplot_i(False)
                signal_obj.setP(d)
                signal_obj.setPD(d)
                plotear(a, b, c, d, e)
        except ValueError:
            a='PPG'
            b=signal_obj.getPPG()
            c=signal_obj.getTime()
            d=signal_obj.getP()
            e='Peak'
            signal_obj.setplot_p(True)
            signal_obj.setplot_i(False)
            plotear(a, b, c, d, e)

    if(accion.text()=="Foot"):
        cmToolMenu.setEnabled(True)
        read_f=signal_obj.getpie()
        uimain.actionSave.setEnabled(True)
        try:
            if read_f==None:       # Verificando si ya se han detectado los pies anteriormente
                a, b, c, d, e = detect_peak()  # a='PPG', b=ppg, c=t, d=picos, e='Peak',
                signal_obj.setP(d)
                ppg = b
                picos = d
                t = c
                areax = list()
                area_position = list()
                pie = list()
                areas_min = list()
                suma = 0
                minimo = np.abs(np.min(ppg))
                for i in range(0, len(picos) - 2):  # Moverse por toda la senal
                    for w in range(picos[i],
                                   picos[i + 1] - 20):  # desplazando la ventana del area entre pico y pico de la senal
                        for m in range(0, 20):  # calculando el area en una ventana de 20 muestras
                            suma = suma + minimo + ppg[(w + m)]
                        areax.append(suma)
                        area_position.append((w + 15))
                        suma = 0

                    Area_array = areax.copy()
                    areax.clear()
                    areas_min.append(np.min(Area_array))
                    areamin = np.argmin(Area_array)
                    pie.append(area_position[areamin])
                    area_position.clear()

                foot = np.array(pie, dtype='int')
                signal_obj.setpie(foot)
                signal_obj.setpieD(foot)
                signal_obj.setplot_i(True)
                signal_obj.setplot_p(False)
                w=plotear('PPG', ppg, t, foot, 'Foot')
        except ValueError:
            a = 'PPG'
            b = signal_obj.getPPG()
            c = signal_obj.getTime()
            d = signal_obj.getpie()
            e = 'Foot'
            signal_obj.setplot_i(True)
            signal_obj.setplot_p(False)
            plotear(a, b, c, d, e)

    if(accion.text()=="Dichrotic"):
        cmToolMenu.setEnabled(True)
        uimain.actionSave.setEnabled(True)
        ppg=signal_obj.getPPG()
        fs=signal_obj.getfs()
        read_d = signal_obj.getdicr()
        read_p=signal_obj.getP()
        read_f=signal_obj.getpie()
        try:
            if read_d == None:  # Verificando si ya se han detectado los puntos dicroticos anteriormente
                try:
                    if read_p== None:   # Verificando si ya se han detectado picos anteriormente
                        a, b, c, d, e = detect_peak()  # a='PPG', b=ppg, c=t, d=picos, e='Peak'
                        picos = d
                except ValueError:
                    picos = read_p
                try:
                    if read_f == None:     # Verificando si ya se han detectado los pies anteriormente
                        a, b, c, d, e = detect_peak()  # a='PPG', b=ppg, c=t, d=picos, e='Peak',
                        signal_obj.setP(d)
                        ppg = b
                        picos = d
                        t = c
                        areax = list()
                        area_position = list()
                        pie = list()
                        areas_min = list()
                        suma = 0
                        minimo = np.abs(np.min(ppg))
                        for i in range(0, len(picos) - 2):  # Moverse por toda la senal
                            for w in range(picos[i],
                                           picos[
                                               i + 1] - 20):  # desplazando la ventana del area entre pico y pico de la senal
                                for m in range(0, 20):  # calculando el area en una ventana de 20 muestras
                                    suma = suma + minimo + ppg[(w + m)]
                                areax.append(suma)
                                area_position.append((w + 15))
                                suma = 0

                            Area_array = areax.copy()
                            areax.clear()
                            areas_min.append(np.min(Area_array))
                            areamin = np.argmin(Area_array)
                            pie.append(area_position[areamin])
                            area_position.clear()

                        foot = np.array(pie, dtype='int')
                        signal_obj.setpie(foot)
                        signal_obj.setpieD(foot)
                        signal_obj.setplot_i(True)
                        signal_obj.setplot_p(False)
                        pies = foot
                except ValueError:
                    pies = read_f
                #Detectando los puntos dicroticos
                blank = 10
                l1 = len(picos)
                dichrotic = list()
                for i in range(0, l1 - 2):
                    posI = picos[i] + blank
                    posF = pies[i] - blank
                    posm = int((posF - posI) / 3)
                    PD = np.diff(ppg[posI:posF])
                    x = False
                    cont = 1
                    l = list()
                    while x == False:
                        PD[PD > 0] = 100
                        minimo1 = np.argmin(np.abs(PD))
                        vmin = np.min(np.abs(PD))
                        l.append(minimo1)
                        if cont == 1:
                            save = minimo1
                            cont = 2
                        if minimo1 > posm - 30 and minimo1 < posm + 30 and PD[minimo1 - 1] != 100:
                            x = True
                            dic3 = int(posI + minimo1 + 1)
                        else:
                            PD[minimo1] = 100
                            if vmin == 100:
                                l = np.array(l)
                                r = l - posm
                                r2 = np.min(np.abs(r))
                                dic3 = int(posI + r2 + 1)
                                x = True
                    dichrotic.append(dic3)

                dichrotic = np.array(dichrotic)
                dic=dichrotic
            signal_obj.setdicr(dic)
            signal_obj.setdicrD(dic)

            a = 'PPG'
            c = signal_obj.getTime()
            d = signal_obj.getdicr()
            e = 'Dichr'
            signal_obj.setplot_i(False)
            signal_obj.setplot_p(False)
            signal_obj.setplot_d(True)
            plotear(a, ppg, c, d, e)
        except ValueError:
            a = 'PPG'
            c = signal_obj.getTime()
            d = signal_obj.getdicr()
            e = 'Dichr'
            signal_obj.setplot_i(False)
            signal_obj.setplot_p(False)
            signal_obj.setplot_d(True)
            plotear(a, ppg, c, d, e)

    if (accion.text() == "Arrhythmia"):
        t=signal_obj.getTime()
        ecg=signal_obj.getECG()
        ppg=signal_obj.getPPG()
        fs=signal_obj.getfs()

        d2 = area.docks['ECG']
        p2 = d2.widgets[0]
        p2.plotItem.clear()
        p21=p2.plot(t, ecg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="ECG")
        p21.sigClicked.connect(ClickedLabel)

        lr2 = pg.LinearRegionItem([0, 10])
        lr2.setZValue(-10)
        lr2.sigRegionChanged.connect(ventana1)
        p2.addItem(lr2)

        d3 = area.docks['PPG']
        p3 = d3.widgets[0]
        p3.plotItem.clear()
        ppg = filtrado_ppg(ppg, fs)
        p31=p3.plot(t, ppg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="PPG")
        p31.sigClicked.connect(ClickedLabel)

        lr3 = pg.LinearRegionItem([0, 10])
        lr3.setZValue(-10)
        p3.addItem(lr3)
        lr3.sigRegionChanged.connect(ventana2)

        uimain.actionhhh.setEnabled(True)
        uimain.actionSave.setEnabled(True)

def clasificar(q):
    menuC.popup(QtGui.QCursor.pos())
    cmToolMenu2.setEnabled(True)
    signal_obj.setLastEvent("Type")

def ventana1(object):
    menuC.setEnabled(True)
    x,y=object.getRegion()
    signal_obj.setWindI(np.int(x))
    signal_obj.setWindF(np.int(y))
    d3 = area.docks['PPG']
    p3 = d3.widgets[0]

    list = p3.plotItem.items

    for objx in list:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue


def ventana2(object):
    menuC.setEnabled(True)
    x, y = object.getRegion()
    signal_obj.setWindI(np.int(x))
    signal_obj.setWindF(np.int(y))
    d2 = area.docks['ECG']
    p2 = d2.widgets[0]

    list = p2.plotItem.items

    for objx in list:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue

def ventana3(object):
    menuC.setEnabled(True)
    x,y=object.getRegion()

    d2 = area.docks['ECG']
    p2 = d2.widgets[0]

    d3 = area.docks['PPG']
    p3 = d3.widgets[0]

    listPPG = p3.plotItem.items
    listECG = p2.plotItem.items

    for objx in listPPG:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue

    p2.removeItem(listECG[-1])
    p2.removeItem(listECG[-1])

    lr5 = pg.LineSegmentROI([(x, 1), (y, 1)], movable=False)
    lr5.removable = False
    lr5.rotateAllowed = False
    lr5.aspectLocked = True
    p2.addItem(lr5)

    lblfmt = "Distance:{p:.2f} s"

    lbl = pg.TextItem(lblfmt.format(p = (y - x)))
    lbl.size = '6pt'
    p2.addItem(lbl)
    w = lbl.pixelHeight()

    vb = p2.getViewBox()
    point = QtCore.QPointF(y + 2, 1.02)

    #if ((y - x - w) < 0):
        #point = vb.mapSceneToView(QtCore.QPointF(y+2, 1.10))
    #    point = QtCore.QPointF(y + 2, 1.20)
    #else:
        #point = vb.mapSceneToView(QtCore.QPointF(x + (y - x  - w) / 2, 1.10))
    #    point = QtCore.QPointF(x + (y - x - w) / 2, 1.20)

    lbl.setPos(point.x(), point.y())




     # pg.LabelItem.set

    #lbl = pg.LabelItem(text=str(x + (x+y)/2))
    #pg.LabelItem.set
    #p2.addItem(lbl)
    # for objy in listECG:
    #     try:
    #         objy.scaleSnap = False
    #         objy.setPos([(x,1) (y, 1)])
    #         break
    #     except:
    #         continue


def ventana4(object):
    menuC.setEnabled(True)
    x, y = object.getRegion()

    d2 = area.docks['ECG']
    p2 = d2.widgets[0]

    #line = pg.LineSegmentROI((x,y))
    #p2.addItem(line)

    list = p2.plotItem.items

    for objx in list:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue

def procTest(q):
    print(q.text())
def clasificacion(accion):
    menuC.setEnabled(False)
    if (accion.text()=="Asystole"):
        signal_obj.setClasses("Asystole")
    if (accion.text() == "Bradycardia"):
        signal_obj.setClasses("Bradycardia")
    if (accion.text() == "Tachycardia"):
        signal_obj.setClasses("Tachycardia")
    if (accion.text() == "Ventricular Tachycardia"):
        signal_obj.setClasses("Ventricular Tachycardia")
    if (accion.text() == "Ventricular Fibrillation"):
        signal_obj.setClasses("Ventricular Fibrillation")
    if (accion.text() == "Low signal quality ECG"):
        signal_obj.setClasses("Low signal quality ECG")
    if (accion.text() == "Low signal quality PPG"):
        signal_obj.setClasses("Low signal quality PPG")
    if (accion.text() == "Flutter"):
            signal_obj.setClasses("Flutter")
    if (accion.text() == "Atrial Fibrillation"):
            signal_obj.setClasses("Atrial Fibrillationr")
    if (accion.text() == "Isolated PVC"):
            signal_obj.setClasses("Isolated PVC")
    if (accion.text() == "Isolated PAC"):
            signal_obj.setClasses("Isolated PAC")
    if (accion.text() == "Bigeminy"):
        signal_obj.setClasses("Bigeminy")
    if (accion.text() == "Trigeminy"):
            signal_obj.setClasses("Trigeminy")

    # pi=signal_obj.getWindI()
    piTotal=signal_obj.getPI()
    #
    # pf=signal_obj.getWindF()
    pfTotal=signal_obj.getPF()
    #
    # tip=signal_obj.getClasses()
    tipTotal=signal_obj.getTipo()
    #
    # n=len(piTotal)
    # lista=list()
    #
    # for i in range(0,n):
    #     lista.append(str(piTotal[i]) + " - " + str(pfTotal[i]) + " : " + tipTotal[i])

    if not signal_obj.getArrhythmiaDialog():
        qf01 = QtWidgets.QFrame()
        qf01.setObjectName('ArrhFrame')
        qf01lm = QtWidgets.QGridLayout()

        qf01.setLayout(qf01lm)

        ql01 = QtWidgets.QListWidget()
        ql01.setFixedWidth(150)
        ql01.setFixedHeight(60)
        ql01.setObjectName('ListaA')

        lbl01 = QtWidgets.QLabel("Arrhythmias stored")
        lbl01.setFixedHeight(20)
        lbl01.setFixedWidth(150)

        btndel = QtWidgets.QPushButton()
        btndel.setText("Actions")
        btndel.setMenu(cmToolMenu2)

        qf01lm.addWidget(lbl01)
        qf01lm.addWidget(ql01, 1, 0, Qt.Alignment(36))
        qf01lm.addWidget(btndel, 2, 0, Qt.Alignment(36))
        d2 = area.docks['ECG']
        d2.addWidget(qf01, row=0, col=1)
        signal_obj.setArrhythmiaDialog()
    else:
        try:
            ql01 = win.findChild(QtWidgets.QListWidget, name='ListaA')
        except Exception as e:
            e.printStackTrace()

        ql01.clear()

    #ql01.add

    l=len(piTotal)
    for i in range(l):
        ql01.addItem(str(piTotal[i]) + " - " + str(pfTotal[i]) + " : " + tipTotal[i])

        #ql01.clicked.connect(another0)
        #cmToolMenu2.setEnabled(True)
        #ql01.clicked.connect(menu_arrhythmia)

        #logger.error(c)


def ventana5(object):
    menuC.setEnabled(True)
    x,y=object.getRegion()
    signal_obj.setWindI(np.int(x))
    signal_obj.setWindF(np.int(y))
    d3 = area.docks['ECG']
    p3 = d3.widgets[0]

    list = p3.plotItem.items
    picosR=signal_obj.getR()
    picosRseg=picosR/250
    pp=signal_obj.getP()
    if(y-x)>10:
        p1 = picosRseg < y     #devuelve un arreglo de True or False
        p0 = picosRseg * p1    #para recuperar los picos verdaderos
        p22 = p0 > x            #devuelve un arreglo de True or False
        p4 = p0 * p22           #para recuperar los picos verdaderos
        p5=p4[p4>0]
        if len(p5) > 2:
            picof = p5[-1]
            picoi = p5[0]
            fc = (len(p5)*60)/(picof-picoi)
            p3.removeItem(list[-1])
            lblfmt = "Heart rate:{p:.0f} bpm"
            lbl1 = pg.TextItem(lblfmt.format(p=fc), color=pg.mkColor(255,255,255), fill=pg.mkBrush(pg.mkColor(0,0,0)))
            lbl2 = pg.TextItem(lblfmt.format(p=fc), color=pg.mkColor(255, 255, 255),fill=pg.mkBrush(pg.mkColor(0, 0, 0)))
            w = lbl1.pixelWidth()
            lbl1.size = '6pt'
            p3.addItem(lbl1)
            point = QtCore.QPointF(x + (y -x -w)/2, 1)
            lbl1.setPos(point.x(), point.y())


    for objx in list:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue


def ventana6(object):
    menuC.setEnabled(True)
    x, y = object.getRegion()
    signal_obj.setWindI(np.int(x))
    signal_obj.setWindF(np.int(y))
    d2 = area.docks['ECG']
    p2 = d2.widgets[0]

    list = p2.plotItem.items
    picosP = signal_obj.getP()
    picosPseg = picosP / 250
    if (y-x)>10:
        p1 = picosPseg < y  # devuelve un arreglo de True or False
        p31 = picosPseg * p1  # para recuperar los picos verdaderos
        p0 = p31 > x  # devuelve un arreglo de True or False
        p4 = p31 * p0  # para recuperar los picos verdaderos
        p5 = p4[p4 > 0]
        if len(p5)>2:
            picof = p5[-1]
            picoi = p5[0]
            fc = (len(p5) * 60) / (picof - picoi)
            p2.removeItem(list[-1])
            lblfmt = "Heart rate:{p:.0f} bpm"
            lbl1 = pg.TextItem(lblfmt.format(p=fc), color=pg.mkColor(255, 255, 255),
                              fill=pg.mkBrush(pg.mkColor(0, 0, 0)))
            w = lbl1.pixelWidth()
            lbl1.size = '6pt'
            p2.addItem(lbl1)
            point = QtCore.QPointF(x + (y - x - w) / 2, 1)
            lbl1.setPos(point.x(), point.y())

    for objx in list:
        try:
            objx.setRegion(object.getRegion())
            break
        except:
            continue
def Medicion(accion):
    ecg=signal_obj.getECG()
    ppg=signal_obj.getPPG()
    t=signal_obj.getTime()
    fs=signal_obj.getfs()
    if(accion.text()=="Distance"):
        d2 = area.docks['ECG']
        p2 = d2.widgets[0]
        p2.plotItem.clear()
        p21 = p2.plot(t, ecg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="ECG")
        p21.sigClicked.connect(ClickedLabel)

        lr2 = pg.LinearRegionItem([0, 10])    #Creando la ventana
        lr2.setZValue(-10)
        lr2.sigRegionChanged.connect(ventana3)
        p2.addItem(lr2)

        #lr4 = pg.LineROI(0, 10, width=0)
        #p2.addItem(lr4)

        lr5 = pg.LineSegmentROI([(0, 1), (10, 1)],movable=False)
        lr5.removable = False
        lr5.rotateAllowed = False
        lr5.aspectLocked = True

        p2.addItem(lr5)

        lbl = pg.TextItem("10")
        lbl.size = '6pt'
        p2.addItem(lbl)

        d3 = area.docks['PPG']
        p3 = d3.widgets[0]
        p3.plotItem.clear()
        ppg = filtrado_ppg(ppg, fs)
        p31 = p3.plot(t, ppg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="PPG")
        p31.sigClicked.connect(ClickedLabel)

        lr3 = pg.LinearRegionItem([0, 10])
        lr3.setZValue(-10)
        p3.addItem(lr3)
        lr3.sigRegionChanged.connect(ventana4)

    if (accion.text() == "Heart rate"):
        d2 = area.docks['ECG']
        p2 = d2.widgets[0]
        p2.plotItem.clear()
        p21 = p2.plot(t, ecg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="ECG")
        p21.sigClicked.connect(ClickedLabel)

        lr2 = pg.LinearRegionItem([0, len(ecg)/fs])
        lr2.setZValue(-10)
        lr2.sigRegionChanged.connect(ventana5)
        p2.addItem(lr2)

        picosR=detqrs(ecg, fs, fac=1.5)
        signal_obj.setR(picosR)
        a, b, c, d, e=detect_peak()
        picosP=d
        signal_obj.setP(picosP)
        fc=len(picosR) * 60 / (len(ecg)/ fs)
        #signal_obj.setFC(fc)

        d3 = area.docks['PPG']
        p3 = d3.widgets[0]
        p3.plotItem.clear()
        ppg = filtrado_ppg(ppg, fs)
        p31 = p3.plot(t, ppg, pen=pg.mkPen(pg.mkColor(182, 240, 9)), name="PPG")
        p31.sigClicked.connect(ClickedLabel)

        # lr3 = pg.LinearRegionItem([0, len(ppg)/fs])
        # lr3.setZValue(-10)
        # p3.addItem(lr3)
        # lr3.sigRegionChanged.connect(ventana6)

        lblfmt = "Heart rate:{p:.0f} bpm"
        lbl = pg.TextItem(lblfmt.format(p=fc), color=pg.mkColor(255,255,255), fill=pg.mkBrush(pg.mkColor(0,0,0)))
        lbl.size = '6pt'
        p2.addItem(lbl)
        point = QtCore.QPointF(0.5 * len(ppg)/fs, 1)
        lbl.setPos(point.x(), point.y())

        # lbl1 = pg.TextItem(lblfmt.format(p=fc), color=pg.mkColor(255, 255, 255), fill=pg.mkBrush(pg.mkColor(0, 0, 0)))
        # lbl1.size = '6pt'
        # p3.addItem(lbl1)
        # point = QtCore.QPointF(0.5 * len(ppg) / fs, 1)
        # lbl1.setPos(point.x(), point.y())



# Main Window construction
app = QtGui.QApplication([])

win = QtWidgets.QMainWindow()
uimain = Ui_MainWindow()
uimain.setupUi(win)

menuFile = win.findChild(QMenu,name="menuFile")
menuFile.triggered[QAction].connect(Archivo)
menuDetect = win.findChild(QMenu,name="menuDetect")
menuDetect.triggered[QAction].connect(Deteccion)
menuMeassure = win.findChild(QMenu,name="menuMeassure")
menuMeassure.triggered[QAction].connect(Medicion)

hiddenMenu = QMenu()
hiddenMenu.addAction(uimain.actionhhh)
hiddenMenu.triggered[QAction].connect(clasificar)
uimain.actionhhh.setEnabled(False)
uimain.actionPeak_2.setEnabled(False)
uimain.actionFoot_2.setEnabled(False)
uimain.actionDichrotic_2.setEnabled(False)
uimain.actionPeak_3.setEnabled(False)
uimain.actionArrhythmia.setEnabled(False)
uimain.actionSave.setEnabled(False)
uimain.actionExit.setEnabled(True)
uimain.actionDistance.setEnabled(False)
uimain.actionHeart_rate.setEnabled(False)

#p = win.findChild(QAction, name="actionhhh")

#uimain.actionhhh.
#triggered[uimain.actionhhh].connect(clasificar)    #connect(clasificar)


area = DockArea()
win.setCentralWidget(area)

pg.setConfigOption("background", pg.mkColor(0,0,0))
pg.setConfigOption("foreground", pg.mkColor(72,240,9))

signal_obj=Senal()

ContextMenu = QMenu()
cmToolMenu = ContextMenu.addMenu("Tools")

adjusta = QAction("Adjust point", win)
remova = QAction("Remove point", win)
adda = QAction("Add point", win)

cmToolMenu.addAction(adjusta)
cmToolMenu.addAction(remova)
cmToolMenu.addAction(adda)

cmToolMenu.triggered[QAction].connect(adjust_points)
cmToolMenu.setEnabled(False)

ContextMenu2=QMenu()
cmToolMenu2=ContextMenu2.addMenu("Tools2")
removarrhythmia = QAction("Remove arrhythmia", win)
clickarrhythmia = QAction("Sort arrhythmia", win)
cmToolMenu2.addAction(removarrhythmia)
cmToolMenu2.addAction(clickarrhythmia)
cmToolMenu2.triggered[QAction].connect(menu_arrhythmia)

#cmToolMenu2.setEnabled(False)


MenuClassif = QMenu()
menuC = MenuClassif.addMenu("Arrythmias")
clase1 = QAction("Asystole", win)
clase2 = QAction("Bradycardia", win)
clase3 = QAction("Tachycardia", win)
clase4 = QAction("Ventricular Tachycardia", win)
clase5 = QAction("Ventricular Fibrillation", win)
clase6 = QAction("Atrial Fibrillation", win)
clase7 = QAction("Flutter", win)
clase8 = QAction("Isolated PVC", win)
clase9 = QAction("Isolated PAC", win)
clase10 = QAction("Bigeminy", win)
clase11 = QAction("Trigeminy", win)
clase12 = QAction("Low signal quality ECG", win)
clase13 = QAction("Low signal quality PPG", win)

menuC.addAction(clase1)
menuC.addAction(clase2)
menuC.addAction(clase3)
menuC.addAction(clase4)
menuC.addAction(clase5)
menuC.addAction(clase6)
menuC.addAction(clase7)
menuC.addAction(clase8)
menuC.addAction(clase9)
menuC.addAction(clase10)
menuC.addAction(clase11)
menuC.addAction(clase12)
menuC.addAction(clase13)
menuC.triggered[QAction].connect(clasificacion)
menuC.setEnabled(False)

cInfo = CursorInfo(False, QtCore.QPointF())
win.show()



# open = QAction("Open", win)
# save = QAction("Save", win)
# peak = QAction("Peak", win)
# foot = QAction("Foot", win)
# dichrotic = QAction("Dichrotic", win)
#
#
# win.
# MenuFile.addAction(open)
# MenuFile.addAction(save)
# MenuDetect.addAction(peak)
# MenuDetect.addAction(foot)
# MenuDetect.addAction(dichrotic)
## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()