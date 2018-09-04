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
from pyCPEDScommonFunctions import cpedsRotation

#from pylab import *
import math
import numpy as np
# from scipy.interpolate import interp1d
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
parser.add_option("-p", "--point", dest="point", default="1,0,0", type="string", help='coma-separated x,y,z coordinates of a point in 3D space', metavar="STRING")
parser.add_option("", "--Rx", dest="Rx", default="90", type="float", help='Rx rotation angle [deg]', metavar="VALUE")


# switches
# parser.add_option("", "--saveClusterPart", action="store_true", dest="saveClusterPart", default=False, help="triggers saving particles positions that form found clusters. This option may produce large amounts of data, but may be useful for control visualizations ")




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

# p=cpedsPythCommon.getFloatList(option.point)
# print p
# print cpedsRotation.Rx(p[0],p[1],p[2],option.Rx*math.pi/180.0)

data=np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]],dtype='float')
# data=np.array([1,2,3])
print data
rotated=cpedsRotation.rotXarray(option.Rx,data)
print rotated

