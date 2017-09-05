#!/usr/bin/env python2.7

'''
Module description: 

This program calculates shape of the Moon when it sets accounting for refraction

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
This program calculates shape of the Moon when it sets accounting for refraction

"""

parser = OptionParser(description=programDescription)

#calculations setup
# parser.add_option("-o", "--outFile", dest="outFile", default="WlWm_vs_X.txt", type="string", help='name of the output file to be read by make_mysql_source_database.sh script', metavar="STRING")
parser.add_option("", "--ZDstart", dest="ZDstart", default=80.0, type="float", help='Initial true ZD [deg]', metavar="VALUE")
parser.add_option("", "--ZDstop", dest="ZDstop", default=92.0, type="float", help='Final true ZD [deg]', metavar="VALUE")
parser.add_option("", "--dZD", dest="dZD", default=0.02, type="float", help='ZD step [deg]', metavar="VALUE")
parser.add_option("", "--moonSize", dest="moonSize", default=0.5, type="float", help='Moon angular diameter [deg]', metavar="VALUE")

# parser.add_option("", "--Wlmax", dest="Wlmax", default=1.0, type="float", help='Maximal value of the Wl0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--dWl", dest="dWl", default=0.1, type="float", help='Step value for the Wl0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--Wmmin", dest="Wmmin", default=0.0, type="float", help='Minimal value of the Wm0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--Wmmax", dest="Wmmax", default=1.0, type="float", help='Maximal value of the Wm0 parameter for the calculations', metavar="VALUE")
# parser.add_option("", "--dWm", dest="dWm", default=0.1, type="float", help='Step value for the Wm0 parameter for the calculations', metavar="VALUE")


# switches
parser.add_option("", "--testRefr", action="store_true", dest="testRefr", default=False, help="test sla refraction for range of wavelengths at the horizon")
parser.add_option("", "--testRefrPW1", action="store_true", dest="testRefrPW1", default=False, help="test sla refraction vs Patrick Wallace refraction for T=-10, P=1000, H=30%")
parser.add_option("", "--testRefrPW2", action="store_true", dest="testRefrPW2", default=False, help="test sla refraction vs Patrick Wallace refraction for T= 30, P=1000, H=90%")
parser.add_option("", "--testRefrPW3", action="store_true", dest="testRefrPW3", default=False, help="test sla refraction vs Patrick Wallace refraction for T= 30, P=950, H=90%")
parser.add_option("", "--testRefrOptical1", action="store_true", dest="testRefrOptical1", default=False, help="test sla radio refraction vs SLA optical refraction for T=-10, P=1000, H=30%")
parser.add_option("", "--testRefrOptical2", action="store_true", dest="testRefrOptical2", default=False, help="test sla radio refraction vs SLA optical refraction for T= 30, P=1000, H=90%")
parser.add_option("", "--testRefrOptical3", action="store_true", dest="testRefrOptical3", default=False, help="test sla radio refraction vs SLA optical refraction for T= 30, P=950, H=90%")
parser.add_option("", "--encode", action="store_true", dest="encode", default=False, help="encode the pictures")

(option, args) = parser.parse_args()


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FNUCTION SPACE


################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM

if option.testRefrOptical1:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,-10,1000,30,5.5e-5,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestOptical1',r)
    sys.exit()

if option.testRefrOptical2:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,30,1000,90,5.5e-5,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestOptical2',r)
    sys.exit()
    
if option.testRefrOptical3:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,30,950,90,5.5e-5,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestOptical3',r)
    sys.exit()


if option.testRefrPW1:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,-10,1000,30,1,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestRadio1',r)
    sys.exit()

if option.testRefrPW2:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,30,1000,90,1,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestRadio2',r)
    sys.exit()
    
if option.testRefrPW3:
    H=list([90,85,80,75,70,65,60,55,50,45,40,35,30,25,20,15,10,9,8,7,6,5,4,3,2,1])
    ZD=0
    r=list()
    for h in H:
        ZD=90.0-h
        cmd='cosmocalc --refr %f,30,950,90,1,133,53,0.0065,1e-9 --script --dontSaveRun' % (ZD)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([ZD,float(ref)]))
        print 'ZD [deg]:',ZD
        ZD=ZD+1
    np.savetxt('refrtestRadio3',r)
    sys.exit()


if option.testRefr:
    lam=0.000001
    r=list()
    while lam<100:
        
        cmd='cosmocalc --refr 89,-10,1024,80,%f,0,51,0.0065,1e-9 --script --dontSaveRun' % (lam)
        res=cpedsPythCommon.getStdOutValues(cmd,['ZDspace [deg]'],':')

        ref=res["ZDspace [deg]"].split(':')[-1]
#         print ref
        r.append(np.array([lam,float(ref)]))
        print 'lambda [cm]:',lam
        lam=lam*1.01
    np.savetxt('refrtest',r)
    sys.exit()

if option.encode:
#     cmd='mencoder mf://Moon*.png -mf w=800:h=600:fps=25:type=png -ovc raw     -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o MoonShape.avi'
    cmd='mencoder mf://Moon*.png -mf w=800:h=600:fps=10:type=png -ovc lavc     -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o MoonShape.avi'
    cpedsPythCommon.sayAndExecute("encoding", cmd,1)
    sys.exit()


ZD=option.ZDstart
frame=0

while ZD<option.ZDstop:
    
    cmd='cosmocalc --refrMoonAtDusk %f,%f,%f,%f,20,1013,100,5.5e-5,113,54,0.0065,1e-8 --dontSaveRun' % (ZD,option.ZDstop,option.dZD,option.moonSize)
    MoonShape=cpedsPythCommon.getStdOutValues(cmd,['horiz/vert obs','elev obs','elev true','dist to lower edge from centre', 'dist to upper edge from centre','minLon obs','maxLon obs','dist to upper edge from centre','dist to lower edge from centre'],':')
    frame+=1
    
    print MoonShape["elev obs"]
    print MoonShape["elev true"]
    print MoonShape["elev true"]
    print float(MoonShape["elev obs"])-float(MoonShape['dist to upper edge from centre'])
#     sys.exit()
#     sys.exit()
    cmd='plot_function.py  --figRows 1 --figCols 2 --plotDs=2,2  '
    cmd+=' --figSize 15,12 --grid  '
    cmd+=' MoonSpace MoonObs MoonSpaceShifted MoonObs --lc k --ls - --lc b --ls -   --ls -- --ls - '
    cmd+=' --xlabel "Relative azimuth [deg]" --ylabel "Elevation [deg]" --Hline 0 --Hline 0  '
    cmd+=' --save -o MoonShape_%03i.png' % frame
    cmd+=" --textxy '0,%s,%s,0,12,b;0,%s,%s,0,12,k'" % ( MoonShape["elev obs"],  MoonShape["elev obs"], MoonShape["elev true"],  MoonShape["elev true"])
    cmd+=" --textxy '0,%s,%s (%s),0,12,b;" % ( MoonShape["elev obs"],  MoonShape["elev obs"], MoonShape["horiz/vert obs"])
    cmd+=" %f,%f,%s,0,12,b; " % ( 0.0-option.moonSize/2,float(MoonShape["elev obs"])+float(MoonShape['dist to upper edge from centre']),MoonShape['dist to upper edge from centre'])
    cmd+=" %f,%f,%s,0,12,b' "  % ( 0.0+option.moonSize/2,float(MoonShape["elev obs"])-float(MoonShape['dist to lower edge from centre']),MoonShape['dist to lower edge from centre'])
#     cmd+=" --xmin %f --xmax %f " % (float(MoonShape["minLon obs"])-option.moonSize/2,float(MoonShape["maxLon obs"])+option.moonSize/2)
    cmd+=" --xmin %f --xmax %f " % (-0.5,0.5)
    cmd+=" --ymin %f --ymax %f " % (float(MoonShape["elev obs"])-option.moonSize*2.5,float(MoonShape["elev obs"])+option.moonSize)
    cmd+=" --equalAspect "
    cmd+=" --operX - --valsX 270 "
    cmd+=" --xticks 0.1 --yticks 0.1 --marginFactor 0 "
#     cmd+=" --Bleft 0.1 --Bright 0.1 --Btop 0.1 --Bbottom 0.1 --Bhspace 0.1 "
    print cmd
    cpedsPythCommon.sayAndExecute("plotting", cmd,1)

    ZD+=option.dZD






