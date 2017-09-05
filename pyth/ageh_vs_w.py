#!/usr/bin/env python2.7

'''
Module description: 

This program calculates age of the Universe as a function of time-independent DE eq. of state parameter

Created on Feb 22, 2012
@author: blew
'''

import sys
sys.path.append('./')

from pylab import *

import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.mlab
import os
# from matplotlib.ticker import FuncFormatter

# from matplotlib.collections import PatchCollection
# import matplotlib.path as mpath
# import matplotlib.patches as mpatches
import csv
from optparse import OptionParser
from pyCPEDScommonFunctions import cpedsPythCommon


programDescription="""
This program calculates age of the Universe as a function of time-independent DE eq. of state parameter

"""

parser = OptionParser(description=programDescription)

#calculations setup
# parser.add_option("-o", "--outFile", dest="outFile", default="WlWm_vs_X.txt", type="string", help='name of the output file to be read by make_mysql_source_database.sh script', metavar="STRING")
# parser.add_option("", "--ZDstart", dest="ZDstart", default=80.0, type="float", help='Initial true ZD [deg]', metavar="VALUE")
# parser.add_option("", "--ZDstop", dest="ZDstop", default=92.0, type="float", help='Final true ZD [deg]', metavar="VALUE")
# parser.add_option("", "--dZD", dest="dZD", default=0.02, type="float", help='ZD step [deg]', metavar="VALUE")
# parser.add_option("", "--moonSize", dest="moonSize", default=0.5, type="float", help='Moon angular diameter [deg]', metavar="VALUE")

# parser.add_option("", "--Wlmax", dest="Wlmax", default=1.0, type="float", help='Maximal value of the Wl0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--dWl", dest="dWl", default=0.1, type="float", help='Step value for the Wl0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--Wmmin", dest="Wmmin", default=0.0, type="float", help='Minimal value of the Wm0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--Wmmax", dest="Wmmax", default=1.0, type="float", help='Maximal value of the Wm0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--dWm", dest="dWm", default=0.1, type="float", help='Step value for the Wm0 parameter for the calculations', metavar="VALUE")


# switches
#parser.add_option("", "--mkCalSrcDwingeloo", action="store_true", dest="mkCalSrcDwingeloo", default=False, help="triggers calculating the dwingeloo flux calibration sources for RT32 receiver bands")

(option, args) = parser.parse_args()


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FNUCTION SPACE


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM

w=-2.0
wEnd=0.0
h=0.71

agew=list()
ageLCDM=0
cmd='cosmocalc --agez --w0  %f --dontSaveRun' % (-1)
data=cpedsPythCommon.getStdOutValues(cmd,['age at redshift'],':')
age=float(data['age at redshift'])
ageLCDM=age

while w<wEnd:
    print 'w: ',w
    cmd='cosmocalc --agez --w0  %f --dontSaveRun' % (w)
    data=cpedsPythCommon.getStdOutValues(cmd,['age at redshift'],':')
#     print data['age at redshift']
    age=float(data['age at redshift'])
    agew.append(np.array([w,age/ageLCDM]))

    w+=0.1

np.savetxt('ageh_vs_w.txt',agew)




