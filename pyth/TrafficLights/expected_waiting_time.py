#!/usr/bin/env python

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
parser.add_option("", "--ns", dest="ns", default=2, type="int", help='Number of traffic signals', metavar="VALUE")
parser.add_option("", "--nm", dest="nm", default=1, type="int", help='Number of magic wands that can change red to green instantaneously', metavar="VALUE")
parser.add_option("", "--Pred", dest="Pred", default=0.5, type="float", help='Probability of red lights', metavar="VALUE")
parser.add_option("", "--tMax", dest="tMax", default=80.0, type="float", help='Maximal waiting time at red lights', metavar="VALUE")


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
def getRandomCounterValues(cMin,cMax,num):
#     return np.random.uniform(cMin,cMax,num)
    return np.random.random_integers(cMin,cMax,num)

def EX(X):
    return np.mean(X)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# GLOBAL VARIABLES SPACE


sampleSize=10000
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM

#
# choose traffic lights counter values
#

c=getRandomCounterValues(1, option.tMax, option.ns*sampleSize).reshape([option.ns,sampleSize])

ex1=option.Pred*EX(c[0])
print 'EX [s]: ',ex1

Nwands=option.nm
for i in range(option.ns):
    len(c[0][c[0]<ex1])
