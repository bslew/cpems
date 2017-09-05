#!/usr/bin/env python2.7
import sys
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab
import os
from matplotlib.collections import PatchCollection
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from optparse import OptionParser

programDescription="""This is a simple visualisation tool
for plotting the outputs of the domain decoposition from the test-domain-decomposition
program.

Usage: plot_domain_decomposition pointSet_file domains_file.

If you don't have these files then run with option --generate optionally with -N option to specifiy how many points to generate and plot.

"""

parser = OptionParser(description=programDescription)
parser.add_option("-N", "--Npoints", dest="Npoints", default=0, type="int", help='number of points to decompose', metavar="VALUE")
#parser.add_option("", "--zmax", dest="zmax", default=0, type="float", help='maximal z value the color bar is mapped to', metavar="VALUE")
parser.add_option("", "--figSize", dest="figSize", default=10, type="int", help='size of the figure in inch', metavar="SIZE")
parser.add_option("", "--DPIgui", dest="DPIgui", default=70, type="int", help='DPI for plotting on the screen', metavar="NUM")
parser.add_option("", "--DPI", dest="DPI", default=70, type="int", help='DPI for saving to file', metavar="NUM")
parser.add_option("", "--ticks", dest="ticks", default=10, type="int", help='ticks to be plotted every this cells', metavar="NUM")
parser.add_option("", "--title", dest="title", default='', type="string", help='plot title', metavar="STRING")
parser.add_option("", "--xlabel", dest="xlabel", default='x', type="string", help='plot x label', metavar="STRING")
parser.add_option("", "--ylabel", dest="ylabel", default='y', type="string", help='plot y label', metavar="STRING")
#parser.add_option("", "--bad", dest="bad", default=-1, type="float", help='disable plotting pixels with these values ', metavar="VALUE")
parser.add_option("", "--xmin", dest="xmin", default=-1, type="float", help='minimal X value in the grid to show on X axis', metavar="VALUE")
parser.add_option("", "--xmax", dest="xmax", default=-1, type="float", help='maximal X value in the grid to show on X axis', metavar="VALUE")
parser.add_option("", "--ymin", dest="ymin", default=-1, type="float", help='minimal Y value in the grid to show on Y axis', metavar="VALUE")
parser.add_option("", "--ymax", dest="ymax", default=-1, type="float", help='maximal Y value in the grid to show on Y axis', metavar="VALUE")

parser.add_option("", "--rxmin", dest="rxmin", default=-1, type="float", help='rect def to overplot', metavar="VALUE")
parser.add_option("", "--rxmax", dest="rxmax", default=-1, type="float", help='rect def to overplot', metavar="VALUE")
parser.add_option("", "--rymin", dest="rymin", default=-1, type="float", help='rect def to overplot', metavar="VALUE")
parser.add_option("", "--rymax", dest="rymax", default=-1, type="float", help='rect def to overplot', metavar="VALUE")
#parser.add_option("", "--interpolation", dest="interpolation", default='nearest', type="string", help='type of interpolation to be used for plotting matrix: the allowed interpolation strings are those compatible with matplotlib eg: "nearest", "bilinear", "bicubic", "spline16", "spline36", "hanning", "hamming", "hermite", "kaiser", "quadric", "catrom", "gaussian", "bessel", "mitchell", "sinc", "lanczos" (default: nearest) ', metavar="STRING")
parser.add_option("-o", "", dest="outputFile", default='', type="string", help='name of the output file if --save option is used (default: "") ', metavar="STRING")

# switches
parser.add_option("", "--grid", action="store_true", dest="grid", default=False, help="triggers showing the grid on plot")
parser.add_option("", "--generate", action="store_true", dest="generate", default=False, help="triggers generating a random set of points for domain decomposition")
#parser.add_option("", "--CMbone", action="store_true", dest="CMbone", default=False, help="triggers showing the plot with bone colormap")
parser.add_option("", "--save", action="store_true", dest="save", default=False, help="triggers saving the plot a image")
#parser.add_option("", "--colorbar", action="store_true", dest="colorbar", default=False, help="triggers showing a color bar")
#parser.add_option("", "--zrangeSym", action="store_true", dest="zrangeSym", default=False, help="triggers showing z-range to be within the -max(abs(zmin),abs(zmax)) and +max(abs(zmin),abs(zmax)) ")
#parser.add_option("", "--enumSlices", action="store_true", dest="enumSlices", default=False, help="triggers enumerating slices in the plot title by the slice number")
#parser.add_option("", "--makeMovie", action="store_true", dest="makeMovie", default=False, help="triggers making movie from the slices")
#parser.add_option("", "--maskBad", action="store_true", dest="maskBad", default=False, help="triggers masking bad pixels")



(option, args) = parser.parse_args()



#@param boxSlice - array with slices
#@param sliceNo - index of the slice
#@param sliceName - file name prefix of the slice
def plotDomainDecomposition(ps,domains):
#    vmin=boxSlice.min()
#    vmax=boxSlice.max()
#    maxabs=max([abs(vmin),abs(vmax)])
#    print "min: %E" % vmin
#    print "max: %E" % vmax
#    print maxabs
#    print "max-min: %E" % (vmax-vmin)
#
#    Nx=size(boxSlice[:,0])
#    Ny=size(boxSlice[0])

    
    fig=figure(figsize=(option.figSize,option.figSize), dpi=option.DPIgui)
    ax = subplot(111)
    
    #ax1=fig.add_subplot(2,2,1)


#    xy = 0.3, 0.3,
#    width, height = 0.2, 0.5
#    p = mpatches.Rectangle(xy, width, height, facecolor="orange", edgecolor="red")
#    gca().add_patch(p)
#    show()

    print domains
    n=len(domains)
    for i in arange(n):
        xy=domains[i,0], domains[i,1]
        w=domains[i,3]-domains[i,0]
        h= domains[i,4]-domains[i,1]
        print xy,w,h
        if (domains[i,6]==1):
            p=mpatches.Rectangle(xy,w,h, color='r', lw=5*(1-i/n), fill=False, fc='b', hatch='/')
        else:
            p=mpatches.Rectangle(xy,w,h, color='r', lw=5*(1-i/n), fill=False)
        gca().add_patch(p)
#        hatch='/',
#        print i

    # overplot rectangle
#    xy=option.rxmin,option.rymin
#    w=option.rxmax-option.rxmin
#    h=option.rymax-option.rymin
#    p=mpatches.Rectangle(xy,w,h, color='k', lw=2, fill=False, fc='b')
#    gca().add_patch(p)
#    # overplot circles
#    xy=5,5
#    p=mpatches.Circle(xy,3, color='k', lw=2, fill=False, fc='b')
#    gca().add_patch(p)
#    xy=2.5,2.5
#    p=mpatches.Circle(xy,1, color='k', lw=2, fill=False, fc='b')
#    gca().add_patch(p)


        
    # draw point set
    scatter(ps[:,0],ps[:,1],s=20, c='b', marker='o')

    #
    # title and labels
    #

    tit="Domains count: %i. " % len(domains)
    if option.title!="":
        tit=tit+option.title
        
    title(tit)
    xlabel(option.xlabel)
    ylabel(option.ylabel)
        
    #
    # ticks business
    #
#    if option.xmin!=-1 and option.xmax!=-1:
#        ax.xaxis.set_major_formatter(formatterX)
#    else:
#        xticks(arange(0,Nx,option.ticks))
#        xlim([-0.5,Nx-0.5])
#    
#    if option.ymin!=-1 and option.ymax!=-1:
#        ax.yaxis.set_major_formatter(formatterY)
#    else:
#        yticks(arange(0,Ny,option.ticks))
#        ylim([-0.5,Ny-0.5])

    if option.grid:
        grid(True)
    
    if option.save:
        if option.outputFile!='':
            fname=option.outputFile
        else:
            fname=sliceName+".png"
        fig.savefig(fname, dpi=option.DPI)
    else:
        show()
        
    return fig

###########################################################################################
###########################################################################################
###########################################################################################
# MAIN PROGRAM
###########################################################################################
###########################################################################################
###########################################################################################


if len(args)!=2 and option.generate==False:
    print "Too few parameters given. Type --help for help"
    sys.exit(0)

if len(args)==2 and option.generate==False:
    #    load point set
    ps=np.loadtxt(args[0])
    #    load domains
    domains=np.loadtxt(args[1])
        
if len(args)==0 and option.generate==True:
    cmd='test-domain-decomposition %i' % option.Npoints
    os.system(cmd)
    #    load point set
    ps=np.loadtxt('pointSet')
    #    load domains
    domains=np.loadtxt('domains')



plotDomainDecomposition(ps,domains);

