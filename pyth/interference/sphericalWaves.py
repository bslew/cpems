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
parser.add_option("", "--sep", dest="separation", default=2.0, type="float", help='Sound sources separation [m]', metavar="VALUE")
parser.add_option("", "--freq", dest="freq", default=3000.0, type="float", help='Sound wave frequency [Hz]', metavar="VALUE")


# switches
parser.add_option("", "--save", action="store_true", dest="save", default=False, help="Save field source files")
# parser.add_option("", "--plot", action="store_true", dest="save", default=False, help="Save field source files")
parser.add_option("", "--wave4", action="store_true", dest="wave4", default=False, help="Generate pattern with 4 sources of sound wave")




(option, args) = parser.parse_args()

# check for required options
#if not option.rawDataFile:
#    parser.error("No --xxx option supplied. Please specify ....")

# cpedsPythCommon.saveHowWeWereCalled()



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# FUNCTIONS SPACE

def wave(A,x0,y0,l,x,y):
    X=x-x0
    Y=y-y0
    r=np.sqrt(X*X+Y*Y)
    return A*np.sin(2.0*np.pi/l*r)


def mkPlot(f,extent):
    import matplotlib.pyplot as plt
    cmap = plt.cm.seismic

    speakersX=[x1,x2]
    speakersY=[y1+D,y2+D]
    if option.wave4:
        speakersX=[x1,x2,x1,x2]
        speakersY=[y1+D,y2+D,y1,y2]
        
    fig=plt.figure(figsize=(10,8))
    ax=fig.gca()
    im=ax.imshow(f, interpolation='bilinear',cmap=cmap,extent=extent)
    ax.plot(speakersX, speakersY, 's', ms=20, mfc='y', mec='g')
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_title(r'Spherical sound wave interference: $\lambda = %.1f$ cm, $\nu =%.1f$ Hz' % (l*100,option.freq))
    # Use a proxy artist for the colorbar...
    ax.grid(True)
    cbar=fig.colorbar(im)
    cbar.ax.set_ylabel('amplitude')
    fig.savefig('sphericalWaves.png')

    plt.show()
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# GLOBAL VARIABLES SPACE



################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
# MAIN PROGRAM


N=1000 # 
D=10.0 # area of interest [m]
xmin=-D/2.0
xmax=D/2.0
ymin=0.0
ymax=D
dx=D/N
dy=D/N

sep=option.separation
A=1.0 # sound wave source 1 amplitude
x1=-sep/2 # sound wave source 1
y1=0.0 # sound wave source 1
cs=340.29 # speed of sound [m/s]
l=cs/option.freq # wavelength [m]

x2=sep/2 # sound wave source 2
y2=0.0 # sound wave source 2


print 'wavelength [cm]: ',l*100.0
print 'wave frequency [Hz]: ',option.freq
print 'pixel size [m]:', D/N

field1=np.arange(N*N, dtype='float').reshape([N,N])
field2=np.arange(N*N, dtype='float').reshape([N,N])
if option.wave4:
    field3=np.arange(N*N, dtype='float').reshape([N,N])
    field4=np.arange(N*N, dtype='float').reshape([N,N])

for i in range(N):
    x=xmin+i*dx
    for j in range(N):
        y=ymin+j*dy
        field1[j][i]=wave(A,x1,y1,l,x,y)
        field2[j][i]=wave(A,x2,y2,l,x,y)
        if option.wave4:
            field3[j][i]=wave(A,x1,y1+D,l,x,y)
            field4[j][i]=wave(A,x2,y2+D,l,x,y)

if option.wave4:
    mkPlot(field1+field2+field3+field4,extent=[xmin,xmax,ymin,ymax])
else:
    mkPlot(field1+field2,extent=[xmin,xmax,ymin,ymax])

# print field1
if option.save:
    np.savetxt('field1.txt',field1)
    np.savetxt('field2.txt',field2)
    np.savetxt('interference.txt',field1+field2)

