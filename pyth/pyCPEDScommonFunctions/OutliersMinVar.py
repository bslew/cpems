'''
Created on Apr 13, 2017

@author: blew
'''

import numpy as np
from statsmodels import robust
# from pyCPEDScommonFunctions import cpedsPythCommon

class OutliersMinVar():
    '''
    classdocs
    '''


    def __init__(self, data, col=0, Nmin=1):
        '''
        Constructor
        '''
        self.data=data
        self.col=col
        self.dataCol=np.asarray(data[:,col],dtype="float")        
        self.cleaned=self.dataCol
        self.mask=None
        self.thres=None
        self.shist=[]
        self.gaussianProcessRuns=10
        self.dataMinimalLength=Nmin
        
    def getGaussianSample(self, m,s,N):
        data = np.random.normal(0,s,N).reshape((-1,1))
        data+=m-np.mean(data)

    def findOutliers_gaussianProcess(self):

        
        # copy input data
        dataOrig=self.data
        dataColOrig=self.dataCol

        shistory=self.findOutliers(0)
        
        m=np.mean(self.dataCol)
#         m=np.median(self.dataCol)
        s=np.std(self.dataCol)
        for i in range(self.gaussianProcessRuns):
            self.data=dataOrig
            self.dataCol=dataColOrig
            data=self.getGaussianSample(m, s, len(self.data))
            if i==0:
                shistory=self.findOutliers(0)
            else:
                shistory+=self.findOutliers(0)
                
        shistory/=self.gaussianProcessRuns

        

    def findOutliers(self,thres):
        m=np.mean(self.dataCol)
#         m=np.median(self.dataCol)
        s=np.std(self.dataCol)
#         s=robust.mad(self.dataCol)
        
        if s==0: # we only have one data sample
            return
        
        dataCol=self.dataCol

        shist=list()

        cont=True
        outlierCount=0
        if len(dataCol)==self.dataMinimalLength:
            cont=False
        while cont:
            sprev=s
            d=np.abs((dataCol-m))
            
#             print 'data:', dataCol
#             print 'dist: ',d
            
            iMax=d.argmax(axis=0)
#             print 'max d idx:',iMax
            if iMax==0:
                dataCol=dataCol[1:len(dataCol)]
                self.data=self.data[1:len(self.dataCol)]
                outlierCount+=1
            elif iMax==len(dataCol-1):
                dataCol=dataCol[0:iMax]
                self.data=self.data[0:iMax]
                outlierCount+=1
            else:
                dataCol=np.hstack([dataCol[0:iMax],dataCol[iMax+1:]])
                self.data=np.vstack([self.data[0:iMax],self.data[iMax+1:]])
                outlierCount+=1

            
            m=np.mean(dataCol)
#             m=np.median(dataCol)
            s=np.std(dataCol)
#             s=robust.mad(self.dataCol)
            shist.append([outlierCount,m,s,(sprev-s)/sprev])
#             print
#             print 's: ',s
#             print 'sprev: ',sprev
            
            if len(dataCol)==self.dataMinimalLength:
                cont=False
            if np.abs(sprev-s)/sprev<thres:
                cont=False
#             if s>sprev:
#                 cont=False
#                 print 'Removed %i outliers' % outlierCount
            
        self.cleaned=dataCol
        self.shist=shist
#         self.thres=np.median(shist[:,1])
#     def getMaxDevIdx(self):
        
        
        
        
    def getCleanDataCol(self):
        return self.cleaned
    
    
    def getCleanData(self):
        return self.data
    
    
class OutliersMinMAD(OutliersMinVar):
    def __init__(self, data, col=0, Nmin=1):
        '''
        Constructor
        '''
        self.data=data
        self.col=col
        self.dataCol=np.asarray(data[:,col],dtype="float")        
        self.cleaned=self.dataCol
        self.mask=None
        self.thres=None
        self.shist=[]
        self.gaussianProcessRuns=10
        self.dataMinimalLength=Nmin


    def findOutliers(self,thres):
#         m=np.mean(self.dataCol)
        m=np.median(self.dataCol)
#         s=np.std(self.dataCol)
        s=robust.mad(self.dataCol)
        
        if s==0: # we only have one data sample
            return
        
        dataCol=self.dataCol

        shist=list()

        cont=True
        outlierCount=0
        print('Size of data: ',len(dataCol))
        if len(dataCol)==self.dataMinimalLength:
            cont=False
            
        while cont:
            sprev=s
            d=np.abs((dataCol-m))
            
#             print 'data:', dataCol
#             print 'dist: ',d
            
            iMax=d.argmax(axis=0)
#             print 'max d idx:',iMax
            if iMax==0:
                dataCol=dataCol[1:len(dataCol)]
                self.data=self.data[1:len(self.dataCol)]
                outlierCount+=1
            elif iMax==len(dataCol-1):
                dataCol=dataCol[0:iMax]
                self.data=self.data[0:iMax]
                outlierCount+=1
            else:
                dataCol=np.hstack([dataCol[0:iMax],dataCol[iMax+1:]])
                self.data=np.vstack([self.data[0:iMax],self.data[iMax+1:]])
                outlierCount+=1

            
#             m=np.mean(dataCol)
            m=np.median(dataCol)
#             s=np.std(dataCol)
            s=robust.mad(self.dataCol)
            shist.append([outlierCount,m,s,(sprev-s)/sprev])
#             print
#             print 's: ',s
#             print 'sprev: ',sprev
            
            if len(dataCol)==self.dataMinimalLength:
                cont=False
            if np.abs(sprev-s)/sprev<thres:
                cont=False
#             if s>sprev:
#                 cont=False
#                 print 'Removed %i outliers' % outlierCount
            
        self.cleaned=dataCol
        self.shist=shist
#         self.thres=np.median(shist[:,1])
#     def getMaxDevIdx(self):
