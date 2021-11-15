class CursorInfo:

    def __init__(self, enCross, lastPoint):

        self.enCross = enCross
        self.lastPoint = lastPoint
        self.currentLead = 'A'


    def setCross(self, adata):

        self.enCross = adata


    def getCross(self):

        return self.enCross


    def setPoint(self, adata):

        self.lastPoint = adata


    def getPoint(self):

        return self.lastPoint

    def setCurrentLead(self, lead):

        self.currentLead = lead

    def getCurrentLead(self):

        return self.currentLead