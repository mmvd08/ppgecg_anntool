
import pathlib
from io import *
from scipy.io.matlab import savemat
from json import JSONEncoder

class Senal:
    def __init__(self, viewSpectrum=False, viewInfoPanel=False):
        self.data = {'fs': {}, 'ECG': {}, 'PPG': [], 'Time': [], 'R': None, 'P': None, 'pie': None,'dic': None,
                     'RR': None, 'PP': None, 'FC': None, 'PI': [], 'PF': [], 'Tipo': [],'PD': None, 'RD': None, 'pieD': None, 'dicrD': None}
        self.Outliers = {'ECG': [], 'PPG': []}
        self.arrythmia=[]
        self.Classes = []
        self.plot= {'picos': False, 'inicio': False, 'dichr': False}
        self.modif = {'picos': False, 'inicio': False, 'dichr': False}
        self.posWindI=[]
        self.posWindF=[]
        self.haveArrhythmiaDialog = False

        self.viewSpectrum = viewSpectrum
        #self.viewInfoPanel = viewInfoPanel
        self.SelFilter = 0
        #self.SelRDet = 0
        #self.IHB = 0
        #self.LeadCode = 0
        self.__isProcessing = False
        self.__lastEvent = None
        self.__Filtered = False

    ########### Funcion de guardar #############
    def save(self,name):
        fs = self.data['fs']
        ecg = self.data['ECG']
        ppg = self.data['PPG']
        time = self.data['Time']
        pi=self.data['PI']
        pf=self.data['PF']
        t=self.data['Tipo']

        try:
            if self.data['R'] == None:
                r= []
        except ValueError:
            r = self.data['R']
        try:
            if self.data['P'] == None:
                p = []
        except ValueError:
            p = self.data['P']
        try:
            if self.data['pie'] == None:
                pie = []
        except ValueError:
            pie = self.data['pie']
        try:
            if self.data['dic'] == None:
                dic = []
        except ValueError:
            dic = self.data['dic']
        try:
            if self.data['RR'] == None:
                rr = []
        except ValueError:
            rr = self.data['RR']
        try:
            if self.data['PP'] == None:
                pp = []
        except ValueError:
            pp = self.data['PP']
        try:
            if self.data['FC'] == None:
                fc = []
        except ValueError:
            fc = self.data['FC']
        try:
            if self.data['PD'] == None:
                pd = []
        except ValueError:
            pd = self.data['PD']
        try:
            if self.data['RD'] == None:
                rd = []
        except ValueError:
            rd = self.data['RD']
        try:
            if self.data['pieD'] == None:
                pieD = []
        except ValueError:
            pieD = self.data['pieD']
        try:
            if self.data['dicrD'] == None:
                dicrD = []
        except ValueError:
            dicrD = self.data['dicrD']

        struct = {'fs': fs, 'ECG': ecg, 'PPG': ppg, 'Time': time, 'R': r, 'P': p, 'pie': pie, 'dic': dic, 'RR': rr,
                  'PP': pp, 'FC': fc, 'PI': pi, 'PF': pf, 'Tipo': t, 'PD': pd, 'RD': rd, 'pieD': pieD, 'dicrD': dicrD}
        savemat(name, struct)

    ########### Funciones de self.plot #############
    def setplot_p(self,valor):
        self.plot['picos']=valor

    def getplot_p(self):
        return self.plot['picos']

    def setplot_i(self,valor):
        self.plot['inicio']=valor

    def getplot_i(self):
        return self.plot['inicio']

    def setplot_d(self,valor):
        self.plot['dichr']=valor

    def getplot_d(self):
        return self.plot['dichr']

    ########## Funciones de self.modif############

    def setmodif_p(self, valor):
        self.modif['picos'] = valor

    def getmodif_p(self):

        return self.modif['picos']

    def setmodif_i(self, valor):
        self.modif['inicio'] = valor

    def getmodif_i(self):
        return self.modif['inicio']

    def setmodif_p(self, valor):
        self.modif['dichr'] = valor

    def getmodif_p(self):
        return self.modif['dichr']

    ########## Funciones de Processing############
    def setProcessing(self, proc=False):
        self.__isProcessing = proc

    def getProcessing(self):
        return self.__isProcessing

    ########## Funciones de self.data############
    def setfs(self, fs):
        self.data['fs'] = fs

    def getfs(self):
        return self.data['fs']

    def setECG(self, asignal):
        self.data['ECG'] = asignal

    def getECG(self):
        return self.data['ECG']

    def setPPG(self, asignal):
        self.data['PPG'] = asignal

    def getPPG(self):
        return self.data['PPG']

    def setTime(self, atime):
        self.data['Time'] = atime

    def getTime(self):
        return self.data['Time']

    def setR(self, R):
        self.data['R'] = R

    def getR(self):
        return self.data['R']

    def setP(self, P):
        self.data['P'] = P

    def getP(self):
        return self.data['P']

    def setpie(self, pie):
        self.data['pie'] = pie

    def getpie(self):
        return self.data['pie']

    def setdicr(self, dic):
        self.data['dic'] = dic

    def getdicr(self):
        return self.data['dic']

    def setRR(self, rr):
        self.data['RR'] = rr

    def getRR(self):
        return self.data['RR']

    def setPP(self, pp):
        self.data['PP'] = pp

    def getPP(self):
        return self.data['PP']

    def setFC(self, fc):
        self.data['FC'] = fc

    def getFC(self):
        return self.data['FC']

    def setPI(self, valor):
        self.data['PI'].append(valor)

    def getPI(self):
        return self.data['PI'].copy()

    def setPF(self, valor):
        self.data['PF'].append(valor)

    def getPF(self):
        return self.data['PF'].copy()

    def setTipo(self, valor):
        self.data['Tipo'].append(valor)

    def getTipo(self):
        return self.data['Tipo'].copy()

    def setPD(self, PD):
        self.data['PD'] = PD

    def getPD(self):
        return self.data['PD']

    def setRD(self, RD):
        self.data['RD'] = RD

    def getRD(self):
        return self.data['RD']

    def setpieD(self, pieD):
        self.data['pieD'] = pieD

    def getpieD(self):
        return self.data['pieD']

    def setdicrD(self, dicrD):
        self.data['dicrD'] = dicrD

    def getdicrD(self):
        return self.data['dicrD']

    ########## Funciones de posicion de ventana############
    def setWindI(self,valor):
        self.posWindI=valor

    def getWindF(self):
        return self.posWindF

    def setWindF(self,valor):
        self.posWindF=valor

    def getWindI(self):
        return self.posWindI

    ########## ############
    def setLastEvent(self, ev):
        self.__lastEvent = ev

    def getLastEvent(self):
        return self.__lastEvent

    def setViewInfoPanel(self, viewInfoPanel=False):
        self.viewInfoPanel = viewInfoPanel

    def getViewInfoPanell(self):
        return self.viewInfoPanel

    def setViewSpectrum(self, viewSpectrum=False):
        self.viewSpectrum = viewSpectrum

    def getViewSpectrum(self):
        return self.viewSpectrum


    def setFiltered(self, filt=False):
        self.__Filtered = filt

    def getFiltered(self):
        return self.__Filtered

    def setOutLier(self, lead, olist):
        self.Outliers[lead] = olist

    def getOutLiers(self):
        return self.Outliers

    def setClasses(self, clase):
        pi=self.getWindI()
        pf=self.getWindF()
        self.setPI(pi)
        self.setPF(pf)
        self.setTipo(clase)
        self.Classes = clase

    def getClasses(self):
        return self.Classes

    def getArrhythmiaDialog(self):
        return self.haveArrhythmiaDialog

    def setArrhythmiaDialog(self, val=True):
        self.haveArrhythmiaDialog = val

  ########## Funcion clear ############
    def clear_PI (self):
        self.data['PI'].clear()
    def clear_PF (self):
        self.data['PF'].clear()
    def clear_Tipo (self):
        self.data['Tipo'].clear()

    def signal_clear (self):
          self.data['fs'] = []
          self.data['ECG'] = []
          self.data['PPG'] = []
          self.data['Time'] = []
          self.data['R'] = None
          self.data['P'] = None
          self.data['pie'] = None
          self.data['dic'] = None
          self.data['RR'] = None
          self.data['PP'] = None
          self.data['FC'] = None
          self.data['PI'].clear()
          self.data['PF'].clear()
          self.data['Tipo'].clear()
          self.data['PD'] = None
          self.data['RD'] = None
          self.data['pieD'] = None
          self.data['dicrD'] = None

    def removeItemArr(self, index):
        self.data['PI'].remove(self.data['PI'][index])
        self.data['PF'].remove(self.data['PF'][index])
        self.data['Tipo'].remove(self.data['Tipo'][index])

