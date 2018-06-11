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
Converts txt files containing JD dates to the same files but containing date strings in place 
of the JD time columns

example usage:



"""

parser = OptionParser(description=programDescription)

#options
# parser.add_option("", "--data", dest="data", default="", type="string", help='name of the data file: if gauss3000 then 3000 gaussian samples is generated internally', metavar="STRING")
parser.add_option("-c", "--col", dest="col", default=0, type="int", help='column in data file to look for JD date. ', metavar="VAL")
parser.add_option("-o", "", dest="outfile", default="out", type="string", help='output file name', metavar="STRING")
parser.add_option("", "--offset", dest="offset", default=0, type="float", help='time offset to apply to the input JD times [JD]', metavar="VALUE")
parser.add_option("", "--fmt", dest="fmt", default="%Y-%m-%d %H:%M:%S", type="string", help='output file time format (default: %Y-%m-%d %H:%M:%S)', metavar="STRING")


# switches
parser.add_option("", "--test", action="store_true", dest="test", default=False, help="triggers test mode")
parser.add_option("", "--testFmt", action="store_true", dest="testFmt", default=False, help="triggers test mode")




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

if option.test:
    dt='2017-10-27 09:10:11.23'
    val=cpedsPythCommon.cal2jd(date_time_str=dt,offset=option.offset)
    print dt,' UTC is %.15f' % val
    shouldBe=2458053.882074421271682
    print 'should be: %.15f' % shouldBe
    print 'diff: ',val-shouldBe
    sys.exit()

if option.testFmt:
    dt='2017 10 27 09 10 11.23'
    val=cpedsPythCommon.cal2jd(date_time_str=dt,offset=option.offset,DT_FMT=option.fmt)
    print dt,' UTC is %.15f' % val
    shouldBe=2458053.882074421271682
    print 'should be: %.15f' % shouldBe
    print 'diff: ',val-shouldBe
    sys.exit()

from matplotlib.dates import strpdate2num
data=np.loadtxt(args[0], dtype="float")
# print data

dt=cpedsPythCommon.jd2cal(jd=data[:,option.col], DT_FMT=option.fmt).reshape([len(data),1])
print np.shape(data[:,range(option.col+1,len(data[0]),1)])
print np.shape(dt)
data=np.hstack([dt,data[:,range(option.col+1,len(data[0]),1)]])
# data=np.array(data,dtype="float")

np.savetxt(option.outfile, data, fmt="%s")
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
    
