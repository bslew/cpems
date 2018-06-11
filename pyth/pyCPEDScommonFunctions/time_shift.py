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
parser.add_option("", "--outDir", dest="outDir", default="", type="string", help='name of the output directory. If not given then will be given automatically from run parameters', metavar="STRING")


# switches
parser.add_option("", "--saveClusterPart", action="store_true", dest="saveClusterPart", default=False, help="triggers saving particles positions that form found clusters. This option may produce large amounts of data, but may be useful for control visualizations ")




(option, args) = parser.parse_args()

# check for required options
#if not option.rawDataFile:
#    parser.error("No --xxx option supplied. Please specify ....")

cpedsPythCommon.saveHowWeWereCalled()



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



