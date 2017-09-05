#!/usr/bin/env python2.7

'''
Module description: matrix operations

Created on Sep 18, 2012
@author: blew
'''

import sys
import os
sys.path.append('./')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
sys.path.append(os.environ['OCRA_TOOLKIT_DIR']+'/scripts/fluxCalibration/')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
from pyCPEDScommonFunctions.cpedsPythCommon import saveHowWeWereCalled
from pyCPEDScommonFunctions import cpedsPythCommon

#from pylab import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab
import os
from matplotlib.ticker import FuncFormatter

from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import csv
from optparse import OptionParser


programDescription="""
program description here please
"""

parser = OptionParser(description=programDescription)

#options
parser.add_option("-o", "--outFile", dest="outFile", default="out.mat", type="string", help='name of the output file', metavar="STRING")
#parser.add_option("-i", "--inFile", dest="dataFile", default="controlSystem.timingInfo", type="string", help='name of the input data file to be processed', metavar="STRING")

# switches
parser.add_option("", "--add", action="store_true", dest="add", default=False, help="triggers calculating sum of two matrices of the same size")
parser.add_option("", "--sub", action="store_true", dest="sub", default=False, help="triggers calculating difference of two matrices of the same size")
#parser.add_option("", "--save", action="store_true", dest="savePlot", default=False, help="triggers saving plot to file")

#parser.add_option("", "--plotSpecial1", action="store_true", dest="plotSpecial1", default=False, help="plot special block1")

(option, args) = parser.parse_args()

# check for required options
#if not option.rawDataFile:
#    parser.error("No --xxx option supplied. Please specify ....")


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


if option.add:
    res=np.array([])
    m=np.array([])
    for i in np.arange(len(args)):
        m=np.loadtxt(args[i])
        if i==0:
            res=m
        else:
            res=res+m
    np.savetxt(option.outFile,res)

if option.sub:
    res=np.array([])
    m=np.array([])
    for i in np.arange(len(args)):
        m=np.loadtxt(args[i])
        if i==0:
            res=m
        else:
            res=res-m
    np.savetxt(option.outFile,res)


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################



