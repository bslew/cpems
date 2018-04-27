#!/usr/bin/env python2.7

'''
Module description: 

Created on Sep 18, 2012
@author: blew
'''

import sys
import os
sys.path.append('./')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
#sys.path.append(os.environ['OCRA_TOOLKIT_DIR']+'/scripts/fluxCalibration/')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
from pyCPEDScommonFunctions.cpedsPythCommon import saveHowWeWereCalled
from pyCPEDScommonFunctions import cpedsPythCommon
from pyCPEDScommonFunctions import OutliersMinVar

#from pylab import *

import numpy as np
from scipy.interpolate import interp1d
#import itertools
#import matplotlib.pyplot as plt
#import matplotlib.mlab
#from matplotlib.ticker import FuncFormatter

#from matplotlib.collections import PatchCollection
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
import csv
from optparse import OptionParser



programDescription="""

example usage:



"""

parser = OptionParser(description=programDescription)

#options
parser.add_option("", "--data", dest="data", default="", type="string", help='name of the data file: if gauss3000 then 3000 gaussian samples is generated internally', metavar="STRING")
parser.add_option("", "--col", dest="col", default=0, type="int", help='column in data file to look for outliers', metavar="VAL")
parser.add_option("", "--thres", dest="thres", default=0.1, type="float", help='threshold to remove outliers. Defines fraction by which if the cleaned sample standard deviation drops below, the outlier removal is stopped.', metavar="VAL")
parser.add_option("-m", "--mean", dest="mean", default=0.0, type="float", help='generated sample mean', metavar="VAL")
parser.add_option("-s", "--dev", dest="dev", default=1.0, type="float", help='generated sample deviation', metavar="VAL")
parser.add_option("-N", "--size", dest="N", default=3000, type="int", help='generated sample size', metavar="VAL")


# switches
# parser.add_option("", "--saveClusterPart", action="store_true", dest="saveClusterPart", default=False, help="triggers saving particles positions that form found clusters. This option may produce large amounts of data, but may be useful for control visualizations ")




(option, args) = parser.parse_args()

# check for required options
#if not option.rawDataFile:
#    parser.error("No --xxx option supplied. Please specify ....")

# cpedsPythCommon.saveHowWeWereCalled()



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FUNCTIONS SPACE


class Test():
    
    def __init__(self, data):
        self.data=list(data)
        
        del self.data[-1]


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# GLOBAL VARIABLES SPACE
data=list()
data.append(list([1,2,3]))
data.append(list([1,2,3]))
data.append(list([1,2,3]))

print data
test=Test(data)
print data

sys.exit(0)



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM

if option.data=='gauss':
    data = np.random.normal(option.mean,option.dev,option.N).reshape((-1,1))
#     option.col=0
    data+=option.mean-np.mean(data)
#     s=np.std(data)
#     data*=option.dev/s
    print 'mean: ',np.mean(data)
    print 'std: ',np.std(data)
    print 'sample size: ',len(data)

elif option.data=='uniform':
    data = np.random.uniform(-1,1,option.N).reshape((-1,1))
    s=np.std(data)
    data*=option.dev/s
    data+=option.mean-np.mean(data)
    option.col=0
    print 'mean: ',np.mean(data)
    print 'std: ',np.std(data)
    print 'sample size: ',len(data)
    
else:
    data=np.loadtxt(option.data,dtype="string")
    datacol=np.asarray(data[:,option.col], dtype="float")
    print 'mean: ',np.mean(datacol)
    print 'std: ',np.std(datacol)
    print 'sample size: ',len(datacol)
    
out,shist=cpedsPythCommon.removeOutliersByMinimizingSampleVariance(data,option.col,option.thres)
np.savetxt(option.data+'.clean', out, fmt="%s", header='clean by minimal sample variance method')
np.savetxt(option.data+'.shist', shist, header='clean by minimal sample variance method')

