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
Amend the input data file with data from the file being joined using the interpolated
values from the joined file at the values from the first column in the input file.

For example:
The input data is two column data file with
JD in 1st column and some data values on 2nd column

We want to join this file with content of the --join file_to_join that are interpolated
at JD values from the input file.
The file_to_join has columns JD,col2,col3,...

The program loads the input data file and the file_to_join and interpolates all remaining
columns from file_to_join at JD values from input file and adds those columns 
to the output file.

The output file contains columns
JD, col2_in, col2,col3,...


example usage:



"""

parser = OptionParser(description=programDescription)

#options
parser.add_option("-o", "", dest="outfile", default="out", type="string", help='output file name', metavar="STRING")
parser.add_option("-j", "--join", dest="file_to_join", default="", type="string", help='file to join to the input file', metavar="STRING")


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

indata=np.loadtxt(args[0])
interdata=np.loadtxt(option.file_to_join)

x=interdata[:,0]
y=interdata[:,range(1,len(interdata[0]),1)]

f = interp1d(x, y, kind='linear', axis=0)

intercols=f(indata[:,0])

outdata=np.hstack([indata,intercols])

np.savetxt(option.outfile, outdata)


