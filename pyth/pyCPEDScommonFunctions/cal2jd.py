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
import scipy
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
Converts txt files containing dates to the same files but containing JD in place 
of the time columns

example usage:



"""

parser = OptionParser(description=programDescription)

#options
# parser.add_option("", "--data", dest="data", default="", type="string", help='name of the data file: if gauss3000 then 3000 gaussian samples is generated internally', metavar="STRING")
parser.add_option("-c", "--col", dest="col", default=0, type="int", help='column in data file to look for date. It should have format yyyy-mm-dd H:M:S', metavar="VAL")
parser.add_option("-o", "", dest="outfile", default="out", type="string", help='output file name', metavar="STRING")


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


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# GLOBAL VARIABLES SPACE



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM
# print args
data=np.loadtxt(args[0], dtype="string")

dt=data[:,[option.col,option.col+1]]
dt=map(' '.join, zip(dt[:,0],dt[:,1]))
# print dt
jd=cpedsPythCommon.cal2jd(dt)

data[:,option.col]=jd
data=scipy.delete(data,option.col+1,1)
data=np.array(data,dtype="float")
# print data
np.savetxt(option.outfile, data, fmt="%.6f")
# print data

# data1=None
# data2=None
# if option.col>0:
#     data1=data[:,np.arange(0,option.col)]
# if option.col<len(data[0])-2:
#     data2=data[:,np.arange(option.col+2,len(data[0]))]
# 
# if data1!=None:
#     out=np.hstack([data1,jd])
# if data2!=None:
#     out=np.hstack([data1,jd])
    
