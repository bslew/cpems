#!/usr/bin/env python2.7

'''
Module description: 

This program calculates sigma8 as a function of h and w0 using camb.

Created on Feb 22, 2012
@author: blew
'''

import sys
sys.path.append('./')

from pylab import *

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab
import os
# from matplotlib.ticker import FuncFormatter

# from matplotlib.collections import PatchCollection
# import matplotlib.path as mpath
# import matplotlib.patches as mpatches
import csv
from optparse import OptionParser
from pyCPEDScommonFunctions import cpedsPythCommon


programDescription="""
This program calculates sigma8 as a function of h and w0 using camb.
The calculations are done using the external CPEDS cosmocalc program and CAMB.

"""

parser = OptionParser(description=programDescription)

#calculations setup
parser.add_option("-o", "--outFile", dest="outFile", default="WlWm_vs_X.txt", type="string", help='name of the output file to be read by make_mysql_source_database.sh script', metavar="STRING")
parser.add_option("", "--hmin", dest="hmin", default=0.3, type="float", help='Minimal value of the Wl0 parameter for the calculations', metavar="VALUE")
parser.add_option("", "--hmax", dest="hmax", default=1.1, type="float", help='Maximal value of the Wl0 parameter for the calculations', metavar="VALUE")
parser.add_option("", "--dh", dest="dh", default=0.1, type="float", help='Step value for the Wl0 parameter for the calculations', metavar="VALUE")
parser.add_option("", "--wmin", dest="wmin", default=-2.0, type="float", help='Minimal value of the Wm0 parameter for the calculations', metavar="VALUE")
parser.add_option("", "--wmax", dest="wmax", default=0.0, type="float", help='Maximal value of the Wm0 parameter for the calculations', metavar="VALUE")
parser.add_option("", "--dw", dest="dw", default=0.1, type="float", help='Step value for the Wm0 parameter for the calculations', metavar="VALUE")

#cosmology
parser.add_option("", "--Wb", dest="Wb", default=0.044, type="float", help='Omega_b0', metavar="VAL")
parser.add_option("", "--H0", dest="H0", default=70.4, type="float", help='Hubble parameter [km/s/Mpc]', metavar="VAL")
parser.add_option("", "--w0", dest="w0", default=-1, type="float", help='DE EoS parameter today', metavar="VAL")
parser.add_option("", "--z0", dest="z0", default=0, type="float", help='redshift at which evaluate the quantity', metavar="VAL")

parser.add_option("", "--NP", dest="NP", default=10000, type="int", help='accuracy parameter for cosmocalc program (10000)', metavar="VAL")


# switches
#parser.add_option("", "--mkCalSrcDwingeloo", action="store_true", dest="mkCalSrcDwingeloo", default=False, help="triggers calculating the dwingeloo flux calibration sources for RT32 receiver bands")
# parser.add_option("", "--vsAge", action="store_true", dest="vsAge", default=False, help="triggers calculating the age on a regular grid defined by the WXmin/max parameters")
# parser.add_option("", "--vsdL", action="store_true", dest="vsdL", default=False, help="triggers calculating the luminosity distance on a regular grid defined by the WXmin/max parameters")
# parser.add_option("-q", "--quiet", action="store_true", dest="quiet", default=False, help="triggers small verbocity")
# parser.add_option("", "--plot", action="store_true", dest="plot", default=False, help="triggers plotting the results")

#parser.add_option("", "--plotSpecial1", action="store_true", dest="plotSpecial1", default=False, help="plot special block1")

(option, args) = parser.parse_args()


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FNUCTION SPACE


def saveConfig(fname, mode, plotOption="",extraComment=""):
    conftxt=""
    conftxt=conftxt+"RUN MODE: %s\n" % mode
    conftxt=conftxt+"Wlmin: %f\n" % option.Wlmin
    conftxt=conftxt+"Wlmax: %f\n" % option.Wlmax
    conftxt=conftxt+"dWl: %f\n" % option.dWl
    conftxt=conftxt+"Wmmin: %f\n" % option.Wmmin
    conftxt=conftxt+"Wmmax: %f\n" % option.Wmmax
    conftxt=conftxt+"dWm: %f\n" % option.dWm
    conftxt=conftxt+"H0 [km/s/Mpc]: %f\n" % option.H0
    conftxt=conftxt+"z0: %f\n" % option.z0
    conftxt=conftxt+"w0: %f\n" % option.w0
    conftxt=conftxt+"Wb0: %f\n" % option.Wb
    conftxt=conftxt+"Np: %i\n" % option.NP
    conftxt=conftxt+"\n"
    conftxt=conftxt+"Plot option:\n%s\n" % plotOption
    conftxt=conftxt+"\n"
    conftxt=conftxt+"Extra comments:\n%s\n" % extraComment

    f = open(fname, 'w')
    f.write(conftxt)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM




ageL=list()
Nw=len(arange(option.wmin,option.wmax,option.dw))
Nh=len(arange(option.hmin,option.hmax,option.dh))
N=Nw*Nh
i=0

data=list()
for w in arange(option.wmin,option.wmax,option.dw):
    for h in arange(option.hmin,option.hmax,option.dh):
# calculate sigma8 using camb

        cmd="cat %s/pyth/cosmoPlots/w0h_vs_sigma8-camb-LCDM.template" % (os.environ["CPEDS_DIR"])
        cmd+="|sed -e 's/w              = -1.0/w              = %f/' -e 's/hubble         = 70.0/hubble         = %f/' > camb.param" % (w,h*100.0)
        cpedsPythCommon.sayAndExecute("preparing camb script", cmd, True)
        cmd="camb camb.param"
        res=cpedsPythCommon.getStdOutValues(cmd, ['Age of universe','sigma8'], '=')
        sigma8=res['sigma8'].split(' ')[-1]
#         print res
        data.append(np.array([w,h,float(res['Age of universe']),sigma8]))
        print data[-1]
        a=np.asarray(data, dtype='float')
        np.savetxt('w0h_vs_sigma8_age', a)
# #       calculate the age  
#         cosmoStr='  --Wb0 %f  --Wcdm0 %f --Wl0 %f --h %f  --w0 %f  ' % (option.Wb, Wm-option.Wb, Wl, h, w) 
#         cmd='cosmocalc --agez -z %f --NP %i --script -v 0' % (option.z0, option.NP )  +cosmoStr
#         age=float(getStdOutValue(cmd, "age at redshift",":"))
#         ageL.append(age)
#         i=i+1
#         print "done: %lf %%  (w=%f, h=%f)" % (float(i)/N*100, w, h)
#         
# ageM=array(ageL)
# ageM=ageM.reshape(Nw,Nh)    
# np.savetxt(option.outFile,ageM)
# pltOption=''
# if option.plot:
#     pltOption='plot_matrix.py %s --colorbar --xmin 0 --xmax 1 --ymin 0 --ymax 1 --grid --xlabel "$\Omega_m$" --ylabel "$\Omega_l$" --zlabel "age [Gyr]"  --flipY --contours 10,11,12,13,14,15,16,17,18,19,20 --title "Contours: 10,11,12,13,14,15,16,17,18,19,20" --ticks 20 --save -o %s.png' % ( option.outFile, option.outFile)
#     os.system(pltOption);
# saveConfig(option.outFile+".conf","vsAge",pltOption)

################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################

