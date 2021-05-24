#!/usr/bin/env python2.7
#!/usr/bin/env python
import sys
import numpy as np
import os
import re


#from matplotlib.collections import PatchCollection
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
import h5py
import pyfits
from optparse import OptionParser
from pyCPEDScommonFunctions import cpedsPythCommon



programDescription="""This is a simple matrix plotter.
The input text files should be stored as a matrix of numbers.
These will be plotted with matplotlib.\

The data stored in every row of the input file will be placed horizontally in the plot. The data stored in columns with larger index are placed right of the data stored
in the column with smaller index.
The data stored in every column of the input file are placed vertically in the plot. The data stored in rows with larger index are plotted below the data stored in the rows
with smaller index unless the flipY option is triggered.

Example of usage:\
To make a movie output a series of matrix data stored in files from density-plane_2.0000 to density-plane_2.0029 use eg. such command:
The resulting movie will have rendered 3 fremes per second.\
plot_matrix.py density-plane_2 --st 0 --en 29 --flipY --makeMovie --save --xmin 0 --xmax 10 --ymin 0 --ymax 10 --enumSlices --title "slice #: " --fps 3\\


To plot data stored in hdf5 file (test.hdf5) in dataset slice-012 save every slice from the 3d grid in that dataset starting from 0 and ending at 9 and 
to make movie with 3 frames per second use:
plot_matrix.py  test.hdf5  --hdf5 --hdf5dset slice-012 --colorbar  --st 0 --en 9 --fps 3 --makeMovie --save

"""

parser = OptionParser(formatter=None, description=programDescription)
#parser.add_option("", "--simPrefix", dest="simPrefix", default="", type="string", help='simulation prefix', metavar="PREFIX")

#parser.add_option("-f", "--file", dest="signal_file", default='277_ocraf_absron_010711.dat.totalPower.last500000', 
#                  help="name of the file with input signal. The file should consist of 16 columns with total power signal in order as in the native ocraf format: a1 a2 b1 b2 ...", metavar="FILE")
#parser.add_option("-D", "--boxSize", dest="boxSize", default=256, type="float", help='comoving size of the box [Mpc]', metavar="LENGTH")
parser.add_option("", "--zmin", dest="zmin", default=0, type="float", help='minimal z value the color bar is mapped to', metavar="VALUE")
parser.add_option("", "--zmax", dest="zmax", default=0, type="float", help='maximal z value the color bar is mapped to', metavar="VALUE")
parser.add_option("", "--fzmax", dest="fzmax", default=1.0, type="float", help='set the zmax value as this fraction of the matrix maximal value (default 1)', metavar="VALUE")
parser.add_option("", "--figSize", dest="figSize", default="12,12", type="string", help='size of the figure in inch eg. 19,8. To have A4 plot type A4.' , metavar="STRING")
parser.add_option("", "--extraTextFontSize", dest="extraTextFontSize", default=15, type="int",  help='extra texts font size (default: 15)', metavar="NUM")
parser.add_option("", "--fontSize", dest="fontSize", default=10, type="int", help='font size used for labels', metavar="VALUE")

#parser.add_option("", "--figSize", dest="figSize", default=10, type="int", help='size of the figure in inch', metavar="SIZE")
parser.add_option("", "--DPIgui", dest="DPIgui", default=70, type="int", help='DPI for plotting on the screen', metavar="NUM")
parser.add_option("", "--DPI", dest="DPI", default=70, type="int", help='DPI for saving to file', metavar="NUM")
parser.add_option("", "--ticks", dest="ticks", default=-1, type="int", help='ticks to be plotted every this cells', metavar="NUM")
parser.add_option("", "--xticks", dest="xticks", default=-1, type="float", help='xticks separation', metavar="VALUE")
parser.add_option("", "--yticks", dest="yticks", default=-1, type="float", help='yticks separation', metavar="VALUE")
parser.add_option("", "--ticksFmt", dest="ticksFmt", default="%.2f", type="string", help='default format for ticks (default: %.2f)' , metavar="STRING")
parser.add_option("", "--Rxlabels", dest="Rxlabels", default=0, type="float", help='rotate xticklabels by this angle [deg]', metavar="VALUE")
parser.add_option("", "--Rylabels", dest="Rylabels", default=0, type="float", help='rotate yticklabels by this angle [deg]', metavar="VALUE")
#parser.add_option("", "--ylabelsPos", dest="ylabelsPos", default="left", type="float", help='rotate yticklabels by this angle [deg]', metavar="VALUE")
parser.add_option("", "--ylabels", dest="ylabels", default="", type="string", help='comma separated list of yticklabels ', metavar="VALUE")
parser.add_option("", "--xlabels", dest="xlabels", default="", type="string", help='comma separated list of xticklabels ', metavar="VALUE")
parser.add_option("", "--title", dest="title", default='', type="string", help='plot title', metavar="STRING")
parser.add_option("", "--xlabel", dest="xlabel", default='x', type="string", help='plot x label', metavar="STRING")
parser.add_option("", "--ylabel", dest="ylabel", default='y', type="string", help='plot y label', metavar="STRING")
parser.add_option("", "--zlabel", dest="zlabel", default='z', type="string", help='plot z label', metavar="STRING")
parser.add_option("", "--maskBelow", dest="maskBelow", default=np.nan, type="float", help='disable plotting pixels with values below this value ', metavar="VALUE")
parser.add_option("", "--bad", dest="bad", default=-1, type="float", help='disable plotting pixels with these values ', metavar="VALUE")
parser.add_option("", "--xmin", dest="xmin", default=-1, type="float", help='minimal X value in the grid to show on X axis', metavar="VALUE")
parser.add_option("", "--xmax", dest="xmax", default=-1, type="float", help='maximal X value in the grid to show on X axis', metavar="VALUE")
parser.add_option("", "--ymin", dest="ymin", default=-1, type="float", help='minimal Y value in the grid to show on Y axis', metavar="VALUE")
parser.add_option("", "--ymax", dest="ymax", default=-1, type="float", help='maximal Y value in the grid to show on Y axis', metavar="VALUE")
parser.add_option("", "--interpolation", dest="interpolation", default='nearest', type="string", help='type of interpolation to be used for plotting matrix: the allowed interpolation strings are those compatible with matplotlib eg: "nearest", "bilinear", "bicubic", "spline16", "spline36", "hanning", "hamming", "hermite", "kaiser", "quadric", "catrom", "gaussian", "bessel", "mitchell", "sinc", "lanczos" (default: nearest) ', metavar="STRING")
parser.add_option("-o", "", dest="outputFile", default='', type="string", help='name of the output file if --save option is used (default: "") ', metavar="STRING")
parser.add_option("", "--contours", dest="contours", default='', type="string", help='''a comma separated list of contours to be overplotted. 
To use different contours for every plot separate sets of contours by a colon :. Eg. for imshow with no contours and second plot with contours only, use :0.1,0.2,0.5
(default: "") ''', metavar="STRING")
parser.add_option("", "--contours_w", dest="contours_w", default=0.5, type="float", help='contour line width', metavar="VALUE")
parser.add_option("", "--contours_fontSize", dest="contours_fontSize", default=10, type="float", help='contour font size', metavar="VALUE")

#parser.add_option("", "--hdf5dset", dest="hdf5dset", default='', type="string", help='''name of the dataset from the hdf5 file (default: "") 
#The hdf5 dataset if it is 3d dataset saved by mscsFunction3dregc then the x-coordinate increases rightwards in the plot, and y coordinate increases downwards.    
#''', metavar="STRING")
parser.add_option("", "--hdf5dset", dest="hdf5dset", default='', type="string", help='''comma separated list of names of the datasets from the hdf5 file 
The hdf5 dataset if it is 3d dataset saved by mscsFunction3dregc then the x-coordinate increases rightwards in the plot, and y coordinate increases downwards.
(default: "") ''', metavar="STRING")

parser.add_option("", "--textxy", dest="textxy",  type="string", help='print text on the plot. The text string should be: x,y,"text to be plotted"', metavar="NUM", action="append")
parser.add_option("", "--textfxfy", dest="textfxfy",  type="string", help='print text on the plot. The text string should be: fx,fy,"text to be plotted", where fx,fy are fractions of the field size where the text should be placed', metavar="NUM", action="append")

parser.add_option("", "--hdf5slice", dest="hdf5slice", default=0, type="int", help='for 3d grid data the slice to use along 3rd dimension (default: -1) ', metavar="INT")
parser.add_option("", "--binX", dest="binX", default=1, type="int", help='number of cells to bin along X (this operation is done before transpose). X indicates data that are to be plotted horizontally.', metavar="NUM")
parser.add_option("", "--binY", dest="binY", default=1, type="int", help='number of cells to bin along Y (this operation is done before transpose). Y indicates data that are to be plotted vertically.', metavar="NUM")
parser.add_option("", "--binMethod", dest="binMethod", default="mean", type="string", help='''the how binning is done.
By default it's just a mean value, but in case when small features are binned, their value will be much smaller. In order to preserve the 
original signal amplitudes it may be useful to use a different binning method. Use 'max' binning method to select the maximal value
out of the binned pixels to do that. (default: 'mean')''' , metavar="STRING")
parser.add_option("", "--hdu", dest="hdu", default=1, type="int", help='HDU index in the fits file of which array is to be plotted', metavar="NUM")
parser.add_option("", "--fitsSlice", dest="fitsSlice", default=-1, type="int", help='HDU multi-dimentional array (or picture) slice index of the array to be plotted', metavar="NUM")
parser.add_option("", "--rr", dest="rr", default="", type="string", help='coma separated list of rows to be removed after the data are loaded', metavar="NUM")
parser.add_option("", "--shiftX", dest="shiftX", default=0, type="int", help='shift periodically along X (this operation is done before after transposing). Positive values of shift will move rightwards.', metavar="NUM")
parser.add_option("", "--shiftY", dest="shiftY", default=0, type="int", help='shift periodically along Y (this operation is done before after transposing.). Positive values of shift will move downwards. ', metavar="NUM")
parser.add_option("", "--scaleZ", dest="scaleZ", default="1", type="string", help='Multiply all matrix values by this factor before plotting. Use comma separated list for scale by different values corresponding to different data files plotted', metavar="VALUE")
parser.add_option("", "--Bleft", dest="border_left", default=0.1, type="float", help='border size in the plot from the left', metavar="VALUE")
parser.add_option("", "--Bright", dest="border_right", default=0.1, type="float", help='border size in the plot from the right', metavar="VALUE")
parser.add_option("", "--Btop", dest="border_top", default=0.1, type="float", help='border size in the plot from the top', metavar="VALUE")
parser.add_option("", "--Bbottom", dest="border_bottom", default=0.1, type="float", help='border size in the plot from the bottom', metavar="VALUE")
parser.add_option("", "--Bhspace", dest="border_hspace", default=0.2, type="float", help='border horizontal space between the plots for multiplot mode', metavar="VALUE")
parser.add_option("", "--Bvspace", dest="border_vspace", default=0.2, type="float", help='border vertical space between the plots for multiplot mode', metavar="VALUE")
parser.add_option("", "--circ", dest="circ", type="string",  default='', help='(DEPRECIATED, use --circle now) circles to plot defined as x,y,r. eg. --circ 10,20,5 for circle at x=10 y=20 and radii=5 This creates a new patch on a plot.', metavar="STRING")
parser.add_option("", "--circle", dest="circle", type="string",  default='', help='A NEW extended way of plotting circles. Circles are defined by a semi-colon separated list. Each circle is a comma separated list of parameters: x,y,r[,width[,color,[linestyle]]]. eg. --circle 10,20,5,1,r,-- for circle at x=10 y=20 and radii=5, width=1 red, dashed. This creates a new patch on a plot.', metavar="STRING")
parser.add_option("", "--line", dest="line", type="string",  default='', help='line to plot defined as x,y,dx,dy where x,y is the center point and dx,dy offsets to the second point. eg. --line 10,20,5,6 for a line from x=10 y=20 to x=15 y=26.', metavar="STRING")
parser.add_option("", "--rect", dest="rect", type="string",  default='', help='rectangles to plot defined as x,y,dx,dy where x,y is the center point and dx,dy width,height. eg. --rect 10,20,5,6 for a rectangle at x=10 y=20 and width=5 and height=6. .', metavar="STRING")
parser.add_option("", "--ellipse", dest="ellipse", type="string",  default='', help='circles to plot defined as x,y,rx,ry,ang. eg. --circ 10,20,5,4,20 for ellipse at x=10 y=20 and radiix=4 and radiiy=4, tiled 20 deg to x axis, This creates a new patch on a plot.', metavar="STRING")
parser.add_option("", "--lc", dest="lc",  type="string", help='line color to be used for plotting dataset. (default order is: b,g,r,c,m,y,k,b,g,r,c,m,y,k,b,g,r,c,m,y,k)', metavar="STRING", action="append")
parser.add_option("", "--gc", dest="gc",  type="string", help='grid line color to be used for plotting dataset. (default order is: k)', metavar="STRING", action="append")
parser.add_option("", "--CM", dest="CM",  type="string", default='hot', help='''color palette to be used. Interesting choices are 
    hot, spectral, jet, afmhot, jet, none (default: "hot")''', metavar="STRING" )
parser.add_option("", "--cut", dest="cut",  default="", type="string", help='x1,x2,y1,y2 pixel numbers to cut a rectangle before plotting. This operation is done first. (default>: '')', metavar="STRING")
parser.add_option("", "--contcolor", dest="contcolor",  type="string", default='k', help='''comma separated list of contour colors to be usef for every plot (default: "k")''', metavar="STRING" )
parser.add_option("", "--colorbarPad", dest="colorbarPad",  type="float", default=0.05, help='''colorbar padding (default: 0.15)''', metavar="FLOAT" )
parser.add_option("", "--aspect", dest="aspect",  type="string", default="auto", help='''sets axes aspect ratio (default: "")
Possible values are: auto, equal, normal or num = number''', metavar="FLOAT" )
parser.add_option("", "--scaleExtent", dest="scaleExtent", default=1, type="float", help='rescale the extent of the axes range -- this affects tick labels calculation', metavar="VALUE")


# switches
parser.add_option("", "--centerXextent", action="store_true", dest="centerXextent", default=False, help='''
Changes X ticks of the plot subtracting a value so that the center of the image is at Xtick 0. Useful when --xmin and --xmax are read from hdf5 file.''')
parser.add_option("", "--centerYextent", action="store_true", dest="centerYextent", default=False, help='''
Changes Y ticks of the plot subtracting a value so that the center of the image is at Ytick 0. Useful when --ymin and --ymax are read from hdf5 file.''')
parser.add_option("", "--grid", action="store_true", dest="grid", default=False, help="triggers showing the grid on plot")
#parser.add_option("", "--CMcool", action="store_true", dest="CMcool", default=False, help="triggers showing the plot with cool colormap")
#parser.add_option("", "--CMbone", action="store_true", dest="CMbone", default=False, help="triggers showing the plot with bone colormap")
#parser.add_option("", "--CMnone", action="store_true", dest="CMnone", default=False, help="triggers showing the plot with without colormap")
parser.add_option("", "--save", action="store_true", dest="save", default=False, help="triggers saving the plot a image")
parser.add_option("", "--colorbar", action="store_true", dest="colorbar", default=False, help="triggers showing a color bar")
parser.add_option("", "--zrangeSym", action="store_true", dest="zrangeSym", default=False, help="triggers showing z-range to be within the -max(abs(zmin),abs(zmax)) and +max(abs(zmin),abs(zmax)) ")
parser.add_option("", "--enumSlices", action="store_true", dest="enumSlices", default=False, help="triggers enumerating slices in the plot title by the slice number")
parser.add_option("", "--makeMovie", action="store_true", dest="makeMovie", default=False, help="triggers making movie from the slices")
parser.add_option("", "--maskBad", action="store_true", dest="maskBad", default=False, help="triggers masking bad pixels")
parser.add_option("", "--plotSpecial1", action="store_true", dest="plotSpecial1", default=False, help="plot special block1")
parser.add_option("", "--flipY", action="store_true", dest="flipY", default=False, help="triggers showing the Y axis flipped along with the data")
parser.add_option("", "--flipHDFylabels", action="store_true", dest="flipHDFylabels", default=False, help="triggers loading ymin as ymax and ymax as ymin from hdf file. Only used along with option --hdf5readfndata")
parser.add_option("", "--flipX", action="store_true", dest="flipX", default=False, help="triggers showing the X axis flipped along with the data")
parser.add_option("", "--transposeData", action="store_true", dest="transposeData", default=False, help="transposes the data matrix before plotting it")
parser.add_option("", "--hdf5", action="store_true", dest="hdf5", default=False, help="indicates that the input file is in hdf5 format")
parser.add_option("", "--fits", action="store_true", dest="fits", default=False, help="indicates that the input file is in fits format")
parser.add_option("", "--hdf5dump", action="store_true", dest="hdf5dump", default=False, help="indicates to dump data to txt file as matrix")
parser.add_option("", "--hdf5zscramble", action="store_true", dest="hdf5zscramble", default=False, help="indicates to Zscramble data from all slices defined by --st and --en options")
parser.add_option("", "--hdf5zmaximize", action="store_true", dest="hdf5zmaximize", default=False, help="indicates to maximize the hdf5 data from all slices defined by --st and --en options")
parser.add_option("", "--hdf5zsum", action="store_true", dest="hdf5zsum", default=False, help="indicates to sum data from all slices defined by --st and --en options along z direction")
parser.add_option("", "--hdf5readfndata", action="store_true", dest="hdf5readfndata", default=False, help="indicates to load the mscsFunction3dregc data stored in the hdf5keys to setup the axes")
parser.add_option("", "--log", action="store_true", dest="log", default=False, help="triggers log plot of matrix values")
parser.add_option("", "--abs", action="store_true", dest="abs", default=False, help="triggers absolute value plot of matrix values")
parser.add_option("", "--useNewAxes", action="store_true", dest="useNewAxes", default=False, help='''for circles plotting, this allows to treat the 
given circle coordinates as the coordinates recalculated according to the provided xmin xmax etc and not as the pixel numbers.
The circle center will be placed where you asked according to the specified axes ranges. You must take care that the plotted signal,
coincides with these new coordinates orientation by using apropriate matrix transformations like --flip[XY] or --transpose options.
Typically inverting vertical axes may be needed (--flipY) when ymin > ymax to force y-values increase upwards.''')
parser.add_option("", "--joinPlots", action="store_true", dest="joinPlots", default=False, help="do not close figure (plot) after every plot. Plot data on the same plot (useful for eg image + contours on the same plot)")
parser.add_option("", "--scriptMode", action="store_true", dest="scriptMode", default=False, help="flag to indicate that we are working in script mode for generating pictures and without GUI.")


# movie options
parser.add_option("", "--st", dest="st", default=0, type="int", help='first slice', metavar="NUM")
parser.add_option("", "--en", dest="en", default=0, type="int", help='last slice', metavar="NUM")
parser.add_option("", "--fps", dest="fps", default=10, type="int", help='frames per second', metavar="NUM")


cpedsPythCommon.saveHowWeWereCalled()

(option, args) = parser.parse_args()

if option.scriptMode:
    import matplotlib
    matplotlib.use('Agg')

from pylab import *
import matplotlib.pyplot as plt
import matplotlib.mlab
from matplotlib.ticker import FuncFormatter
import matplotlib.path as mpath
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection



if option.enumSlices and option.title=="":
    print("WARNING: title was not given, enumeration will take no effect.\nUse space for title if you want to enforce single number enumeration.")

ylabels=list()
if option.ylabels!='':
    ylabels=option.ylabels.split(',')
    print(ylabels)

xlabels=list()
if option.xlabels!='':
    xlabels=option.xlabels.split(',')
    print(xlabels)

# the input parameters are integers
def getSymmetricTicks(Asize,step):
#     d=(vmax-vmin)/2.0
#     vmin=-d
#     vmax=d
    print('-------------')
#     if N/2
    N=int(Asize/step)
    N/=2
    print('Asize ',Asize)
    print('step ',step)
    t2=np.arange(Asize/2,Asize+step/10,step)
    t1=np.arange(Asize/2-step,0,-step)
    print(t1)
    print(t2)
    t=np.sort(np.hstack([t1,t2]))
#     if t[-1]!=-t[0]:
#         t=np.hstack([t,[-t[0]]])
#     print
#     print vmin,vmax,N
    print(t)
    print('len(t): ',len(t))
    print() 
    return np.asarray(t)

def convertNegativePixelCoordinates(coord,cMax):
    i=0
    for i in range(len(coord)):
        if coord[i]<0:
            coord[i]=cMax+coord[i]
    return coord

def getFileExtension(fname):
    ext=fname.split('.')
#    print ext
#    sys.exit()
    return ext[-1]


def newAxesConverterX(x):
    Nx=size(boxSlice[0])
    y1=0.0; y2=Nx+0.0;
    x1=option.xmin
    x2=option.xmax
    print(x1,x2,y1,y2,x)
    print((x-x1)*(y2-y1)/(x2-x1)+y1) 
    return ( (x-x1)*(y2-y1)/(x2-x1)+y1 )

def newAxesConverterY(x):
    Ny=size(boxSlice[:,0])
#    print "Ny",Ny
    y1=0.0; y2=Ny+0.0;
#    if option.ymin < option.ymax:
    x1=option.ymin
    x2=option.ymax
    return ( (x-x1)*(y2-y1)/(x2-x1)+y1 )

def newAxesConverterLength(x):
    if option.useNewAxes:
        Ny=size(boxSlice[:,0])
    #    print "Ny",Ny
        y1=0.0; y2=Ny+0.0;
    #    if option.ymin < option.ymax:
        x1=option.ymin
        x2=option.ymax
        dl=  ( (0-x1)*(y2-y1)/(x2-x1)+y1 ) - ( (x-x1)*(y2-y1)/(x2-x1)+y1 )
    else:
        dl=x
    return dl
    

def XTicksFormatter(x, pos):
    'The two args are the value and tick position'
#    Nx=size(boxSlice[0])-1
    Nx=size(boxSlice[0])
#    x1=0.5; x2=Nx+0.5;
    x1=0.0; x2=Nx+0.0;
    y1=option.xmin
    y2=option.xmax
    return option.ticksFmt % ( (x-x1)*(y2-y1)/(x2-x1)+y1 )

def YTicksFormatter(x, pos):
    'The two args are the value and tick position'
#    Ny=size(boxSlice[:,0])-1
    Ny=size(boxSlice[:,0])
#    x1=0.5; x2=Ny+0.5;
    x1=0.0; x2=Ny+0.0;
    y1=option.ymin
    y2=option.ymax
#    if option.flipY:
#        y2=option.ymin
#        y1=option.ymax
    return option.ticksFmt % ( (x-x1)*(y2-y1)/(x2-x1)+y1 )
#    return '%.2f' % ( (x-x1)*(y2-y1)/(x2-x1)+y1 )

def getMatInfo(boxSlice):
    vmin=nanmin(boxSlice)
    vmax=nanmax(boxSlice)
    maxabs=max([abs(vmin),abs(vmax)])
    print("matrix statistics:")
    print("min: %E" % vmin)
    print("max: %E" % vmax)
    print("max-min: %E" % (vmax-vmin))
    print("mean: %E" % np.mean(boxSlice))
    print("st.dev: %E" % np.std(boxSlice))
    print("max absolute: ",maxabs)
    print_mat_stats(boxSlice)
    print("------------------")
    return vmin,vmax,maxabs


#@param boxSlice - array with slices
#@param sliceNo - index of the slice
#@param sliceName - file name prefix of the slice
def makeMatrixPlot(boxSlice,plotNo,sliceName,fig,ax):
#    i=plotNo

    showPlot=False
    newFigure=False
    if option.joinPlots:
        if plotNo==len(boxSlices)-1:
            showPlot=True
        if plotNo==0:
            newFigure=True
    else:
        showPlot=True
        newFigure=True

    vmin=boxSlice.min()
    vmax=boxSlice.max()
    maxabs=max([abs(vmin),abs(vmax)])
    print("min: %E" % vmin)
    print("max: %E" % vmax)
    print(maxabs)
    print("max-min: %E" % (vmax-vmin))
    print("mean: %E" % np.mean(boxSlice))
    print("st.dev: %E" % np.std(boxSlice))
    print('minimal element: ',cpedsPythCommon.get_ij_min2Darray(boxSlice))
    print('maximal element: ',cpedsPythCommon.get_ij_max2Darray(boxSlice))

    Ny=size(boxSlice[:,0])
    Nx=size(boxSlice[0])
    print('Matrix size is: ')
    print('Nx: ',Nx)
    print('Ny: ',Ny)
    
    extent=(option.xmin,option.xmax,option.ymin,option.ymax)
#    extent=(0,Nx+1.0,Ny+1,0)
    extent=(0,Nx,Ny,0)
    
    if option.figSize=="A4":
        option.figSize='11.6929,8.2677'
    else:
        figSize=cpedsPythCommon.getFloatList(option.figSize)

    if newFigure:
#    fig=figure(figsize=(option.figSize,option.figSize), dpi=option.DPIgui)

        fig=figure(figsize=figSize, dpi=option.DPIgui)
        ax=subplot(111)
        fig.subplots_adjust(left=option.border_left, right=1-option.border_right, top=1-option.border_top, bottom=option.border_bottom,  hspace=option.border_vspace, wspace=option.border_hspace)
#    ax=axes()
#    else:
#        ax = fig.get_subplot()


    #
    # color map
    #
    
#    cmap=CM[i]
#    if option.CMcool:
#        cmap=cm.cool
#    if option.CMbone:
#        cmap=cm.bone

#     if CM!='':
#         if CM[plotNo % len(CM)]=='spectral':
#             cmap=cm.spectral
#         if CM[plotNo % len(CM)]=='spectral_r':
#             cmap=cm.spectral_r
#         if CM[plotNo % len(CM)]=='seismic':
#             cmap=cm.seismic
#         if CM[plotNo % len(CM)]=='jet':
#             cmap=cm.jet
#         if CM[plotNo % len(CM)]=='blues':
#             cmap=cm.blues
#         if CM[plotNo % len(CM)]=='afmhot':
#             cmap=cm.afmhot
#         if CM[plotNo % len(CM)]=='hot':
#             cmap=cm.hot
#         if CM[plotNo % len(CM)]=='none':
#             cmap=None

#     cmap=CM[plotNo % len(CM)]

    cmap = plt.get_cmap(CM[plotNo % len(CM)])
    #
    # mask business
    #     
    if option.maskBad:
        print('masking bad:')
        cmap.set_bad('w')
        nx=Nx
        ny=Ny
        for i in np.arange(nx):
            for j in np.arange(ny):
                if boxSlice[i,j]==option.bad:
                    boxSlice[i,j]=np.nan
        vmin,vmax,maxabs=getMatInfo(boxSlice)
        
#        vmin=nanmin(boxSlice)
#        vmax=nanmax(boxSlice)
#        maxabs=max([abs(vmin),abs(vmax)])
#        print "masked statistics:"
#        print "min: %E" % vmin
#        print "max: %E" % vmax
#        print "max-min: %E" % (vmax-vmin)
#        print maxabs

    if math.isnan(option.maskBelow)==False:
        print('masking Nan:')
        cmap.set_bad('w')
        ny=size(boxSlice[0])
        nx=size(boxSlice[:,0])
        print(nx,ny)
        for i in arange(nx):
            for j in arange(ny):
                if boxSlice[i,j]<option.maskBelow:
                    boxSlice[i,j]=np.nan
                    
        vmin,vmax,maxabs=getMatInfo(boxSlice)

    if plotNo==len(boxSlices)-1:
        #
        # plot circles
        #
        circX=circlePatchDef[0::3]
        circY=circlePatchDef[1::3]
        circR=circlePatchDef[2::3]
        
        if option.useNewAxes==False:
            circX=convertNegativePixelCoordinates(circX,Nx)
            circY=convertNegativePixelCoordinates(circY,Ny)
    
        if option.useNewAxes:
            tmpCirc=list()
            for v in circX:
                tmpCirc.append(newAxesConverterX(v))
            circX=tmpCirc
    
            tmpCirc=list()
            for v in circY:
                tmpCirc.append(newAxesConverterY(v))
            circY=tmpCirc
                
            tmpCirc=list()
            for v in circR:
                tmpCirc.append(newAxesConverterLength(v))
    #            tmpCirc.append(float(v))
            circR=tmpCirc
                
        
        
        if len(circlePatchDef)>0:
            axes(ax)
            for i in np.arange(len(circR)):
                print('Plotting circle (old way): ',circX[i],circY[i],circR[i])
                circ = Circle((circX[i],circY[i]), circR[i], fc=None, ec=option.lc[i % len(option.lc)], color=None, fill=False)
                ax.add_patch(circ)

        if option.circle!='':
            circles=option.circle.split(';')
            axes(ax)
            for circle in circles:
                circleSplit=re.split(r'(?<!\\),', circle) # negative lookbehind assertion
                print(circleSplit)
                if len(circleSplit)>=3:
                    cx=float(circleSplit[0])
                    cy=float(circleSplit[1])
                    if option.useNewAxes==False:
                        cx=convertNegativePixelCoordinates(cx, Nx)
                        cy=convertNegativePixelCoordinates(cy, Ny)
                    cx=newAxesConverterX(cx)
                    cy=newAxesConverterY(cy)
                    r =newAxesConverterLength( float(circleSplit[2]) )
                    cw=1
                    ccolor='k'
                    lstyle='solid'
#                     fill=False
                if len(circleSplit)>=4:
                    cw=float(circleSplit[3])
                if len(circleSplit)>=5:
                    ccolor=circleSplit[4]
                if len(circleSplit)>=6:
                    lstyle=circleSplit[5]

                print('Plotting circle: ',circleSplit)
#                 print cx,cy,r,cw,ccolor,lstyle
                circ = Circle((cx,cy), r, fc=None, ec=ccolor, linestyle=lstyle, linewidth=cw, color=None, fill=False)
                ax.add_patch(circ)

            
        #
        # plot ellipses
        #
        if len(ellipsePatchDef)>0:
            axes(ax)
            for i in np.arange(len(ellipsePatchDef)):
                ell = mpatches.Ellipse((ellipsePatchDef[1],ellipsePatchDef[0]), ellipsePatchDef[3],ellipsePatchDef[2], angle=90-ellipsePatchDef[4], fc=None, ec=option.lc[i % len(option.lc)], color=None, fill=False)
                ax.add_patch(ell)

        #
        # plot rectangles
        #
#        rectX=list(np.array(rectanglePatchDef[0::5])-np.array(rectanglePatchDef[2::5])/2)
#        rectY=list(np.array(rectanglePatchDef[1::5])-np.array(rectanglePatchDef[3::5])/2)
        rectX=rectanglePatchDef[0::5]
        rectY=rectanglePatchDef[1::5]
        rectW=rectanglePatchDef[2::5]
        rectH=rectanglePatchDef[3::5]
        rectA=rectanglePatchDef[4::5]
        
#        circX=convertNegativePixelCoordinates(circX,Nx)
#        circY=convertNegativePixelCoordinates(circY,Ny)
    
        if option.useNewAxes:
            tmpRect=list()
            for v in rectX:
                tmpRect.append(newAxesConverterX(v))
            rectX=tmpRect
    
            tmpRect=list()
            for v in rectY:
                tmpRect.append(newAxesConverterY(v))
            rectY=tmpRect
                
            tmpRect=list()
            for v in rectW:
                tmpRect.append(newAxesConverterLength(v))
            rectW=tmpRect

            tmpRect=list()
            for v in rectH:
                tmpRect.append(newAxesConverterLength(v))
            rectH=tmpRect

        if len(rectanglePatchDef)>0:
            axes(ax)
            for i in np.arange(len(rectX)):
#                rect = mpatches.Rectangle((rectanglePatchDef[0]-rectanglePatchDef[2]/2,rectanglePatchDef[1]-rectanglePatchDef[3]/2), rectanglePatchDef[2],rectanglePatchDef[3], angle=rectanglePatchDef[4], fc=None, ec=option.lc[i % len(option.lc)], color=None, fill=False)
                rect = mpatches.Rectangle((rectX[i]-rectW[i]/2,rectY[i]-rectH[i]/2), rectW[i],rectH[i], angle=rectA[i], fc=None, ec=option.lc[i % len(option.lc)], color=None, fill=False)
                ax.add_patch(rect)

        #
        # plot lines
        #
        lineX1=linePatchDef[0::5]
        lineY1=linePatchDef[1::5]
        lineX2=linePatchDef[2::5]
        lineY2=linePatchDef[3::5]
        lineWidth=linePatchDef[4::5]
        
#        circX=convertNegativePixelCoordinates(circX,Nx)
#        circY=convertNegativePixelCoordinates(circY,Ny)
    
        if option.useNewAxes:
            tmpLine=list()
            for v in lineX1:
                tmpLine.append(newAxesConverterX(v))
            lineX1=tmpLine
    
            tmpLine=list()
            for v in lineY1:
                tmpLine.append(newAxesConverterY(v))
            lineY1=tmpLine
                
            tmpLine=list()
            for v in lineX2:
                tmpLine.append(newAxesConverterLength(v))
            lineX2=tmpLine

            tmpLine=list()
            for v in lineY2:
                tmpLine.append(newAxesConverterLength(v))
            lineY2=tmpLine

        if len(linePatchDef)>0:
            axes(ax)
            for i in np.arange(len(lineX1)):
#                rect = mpatches.Rectangle((rectanglePatchDef[0]-rectanglePatchDef[2]/2,rectanglePatchDef[1]-rectanglePatchDef[3]/2), rectanglePatchDef[2],rectanglePatchDef[3], angle=rectanglePatchDef[4], fc=None, ec=option.lc[i % len(option.lc)], color=None, fill=False)
                st=[lineX1[i],lineX1[i]+lineX2[i]]
                en=[lineY1[i],lineY1[i]+lineY2[i]]
                l = plot(st,en,linewidth=lineWidth[i],color=option.lc[i % len(option.lc)])



    #    
    # z-range business
    #
    if option.zmin<option.zmax:
        vmin=option.zmin
        vmax=option.zmax
    else:
        if option.zrangeSym:
            vmin=-maxabs
            vmax=maxabs
    
    print(option.zmin,option.zmin)
    print(vmax)
    vmax=vmax*option.fzmax
#    imshow(boxSlice, interpolation='nearest', cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    print("CM",CM)
    if CM[plotNo % len(CM)]!='none':
#        imshow(boxSlice, interpolation=option.interpolation, cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
        imshow(boxSlice, interpolation=option.interpolation, cmap=cmap, vmin=vmin, vmax=vmax, aspect=option.aspect, extent=extent)
#        extent=(option.xmin,option.xmax,option.ymin,option.ymax)
    #ax1=fig.add_subplot(2,2,1)
    #imshow(boxSlice, interpolation='spline16')
    #imshow(boxSlice, interpolation='bicubic', vmin=-20000, vmax=20000)
    #imshow(boxSlice, interpolation='bicubic')


    #
    # extra texts
    #
    tstrIdx=0
    for tstr in option.textxy:
        if plotNo==len(boxSlices)-1:
            if tstr!='None':
                xyTextTuple=tstr.split(';')
                for tstr in xyTextTuple:
                    tstrSplit=tstr.split(',',2)
                    tx=float(tstrSplit[0])
                    ty=float(tstrSplit[1])
                    if option.useNewAxes:
                        tx=newAxesConverterX(tx)
                        ty=newAxesConverterX(ty) 
                    t=tstrSplit[2]
                    text(tx,ty,t,fontsize=option.extraTextFontSize, zorder=10000)
            
        tstrIdx+=1
    
    tstrIdx=0
    for tstr in option.textfxfy:
        if plotNo==len(boxSlices)-1:
            if tstr!='None':
                xyTextTuple=tstr.split(';')
                for tstr in xyTextTuple:
                    tstrSplit=tstr.split(',',2)
                    tx=float(tstrSplit[0])*Nx
                    ty=float(tstrSplit[1])*Ny
                    if option.useNewAxes:
                        tx=newAxesConverterX(tx)
                        ty=newAxesConverterX(ty) 
                    t=tstrSplit[2]
                    text(tx,ty,t,fontsize=option.extraTextFontSize, zorder=10000)
            
        tstrIdx+=1


    #
    # title and labels
    #
    if plotNo==len(boxSlices)-1:

        if option.title!="":
            tit=option.title
            if option.enumSlices:
                tit=tit+" %04i" % plotNo
            title(tit, fontsize=option.fontSize)
        xlabel(option.xlabel, fontsize=option.fontSize)
        ylabel(option.ylabel, fontsize=option.fontSize)
        
        #
        # ticks business
        #
        option.xmin*=option.scaleExtent
        option.ymin*=option.scaleExtent
        option.xmax*=option.scaleExtent
        option.ymax*=option.scaleExtent
        
        if option.centerXextent:
            dX=option.xmax-option.xmin
            option.xmin=-dX/2.0
            option.xmax=dX/2.0
        if option.centerYextent:
            dY=option.ymax-option.ymin
            option.ymin=-dY/2.0
            option.ymax=dY/2.0
        
        if option.xmin!=-1 or option.xmax!=-1:
            if option.xlabels!='':
                print('changing xlabels')
                xticks(arange(0,Nx,option.xticks), xlabels)
            else:
                if option.centerXextent:
    #                 xticks(getSymmetricTicks(option.xmin, option.xmax, option.xticks))
    #                 xticks(getSymmetricTicks(0, Nx, option.xticks))
                    xticks(getSymmetricTicks(Nx, option.xticks))
                else:
                    xticks(arange(0,Nx+1,option.xticks))
                ax.xaxis.set_major_formatter(formatterX)
                print('formatting x ticks:')
                xt= arange(option.xmin,option.xmax,option.xticks)
                print(xt)
                if len(xt)==0:
                    print('WARNING: something wrong with xticks ? if the ticks are now what you wanted consider changing/setting option.xticks')        
            setp( gca().get_xticklabels(), rotation=option.Rxlabels, horizontalalignment='center', verticalalignment='top')
        else:
            if option.ticks==-1:
                xtcs=Nx/10
                if xtcs==0:
                    xtcs=1
            else:
                xtcs=option.ticks
            xticks(arange(0,Nx,xtcs))
    #        xlim([-0.5,Nx-0.5])
        
        if option.ymin!=-1 or option.ymax!=-1:
            if option.ylabels!='':
                print('changing ylabels')
                yticks(arange(0,Ny,option.yticks), ylabels)
            else:
                if option.centerXextent:
                    yticks(getSymmetricTicks(Ny, option.yticks))
                else:
                    yticks(arange(0,Ny+1,option.yticks))
                ax.yaxis.set_major_formatter(formatterY)
                yt= arange(option.xmin,option.xmax,option.xticks)
                print(yt)
                if len(yt)==0:
                    print('WARNING: something wrong with yticks ? if the ticks are now what you wanted consider changing/setting option.yticks')        
    
            setp( gca().get_yticklabels(), rotation=option.Rylabels, horizontalalignment='right', verticalalignment='center')
        else:
            if option.ticks==-1:
                ytcs=Ny/10
                if ytcs==0:
                    ytcs=1
            else:
                ytcs=option.ticks
            yticks(arange(0,Ny,ytcs))
    #        ylim([-0.5,Ny-0.5])

        if option.fontSize>0:
            setp( ax.get_xticklabels(), fontsize=option.fontSize)
            setp( ax.get_yticklabels(), fontsize=option.fontSize)
        else:
            print('switching off xlabels')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
    
        if option.grid:
            grid(True, color=option.gc[0])
        
        if option.colorbar:
            cb=colorbar(pad=option.colorbarPad)
            cb.set_label(option.zlabel, fontsize=option.fontSize)
            setp( cb.ax.get_xticklabels(), fontsize=option.fontSize)
            setp( cb.ax.get_yticklabels(), fontsize=option.fontSize)
            axes(ax)

    if contours[plotNo % len(contours)]!='':
        axes(ax)
        cont=contours[plotNo % len(contours)].split(',')
        print('cont:',cont)
        if len(contours)>0:
            print('plotting contours')
            CS=contour(arange(Nx),arange(Ny),boxSlice,cont,linewidths=option.contours_w,colors=contcolor[plotNo % len(contcolor)])
#            CS=contour(boxSlice,cont,linewidths=0.5,colors=contcolor[plotNo % len(contcolor)], extent=extent)
            if option.contours_fontSize>0.0:
                clabel(CS, inline=1, fontsize=option.contours_fontSize, fmt=option.ticksFmt)
        
#    xlim([option.xmin,option.xmax])
#    ylim([option.ymin,option.ymax])

    if option.plotSpecial1:
        import specialBlock
        specialBlock.plotMatrixSpecialBlock1(ax)

        
            
    if showPlot:
        if option.save:
            if option.outputFile!='':
                fname=option.outputFile
            else:
                if option.hdf5:
                    fname=sliceName+"-"+option.hdf5dset+"-slice_%03i.png" % option.hdf5slice
                else:
                    fname=sliceName+".png"
            fig.savefig(fname, dpi=option.DPI)
        else:
            show()
        
        
    return fig,ax


def loadHDF5fnData(fname,dset):
#    x0=np.array([])
#    xMax=np.array([])
#    y0=np.array([])
#    yMax=np.array([])
    
    f = h5py.File(fname,'r')
    x0 = f[dset].attrs['x0']
    xMax = f[dset].attrs['xMax']
    y0 = f[dset].attrs['y0']
    yMax = f[dset].attrs['yMax']

    return x0,xMax,y0,yMax

def loadData(fname, sliceNo=0):
    print("* loading file: %s (slice: %li)" % (fname, sliceNo))
    
    if getFileExtension(fname)=='hdf5' or getFileExtension(fname)=='hdf':
        f = h5py.File(fname,'r')
        if option.st < option.en:
            slice = f[option.hdf5dset].value[:,:,sliceNo]
            option.hdf5slice=sliceNo
        else:
            slice = f[option.hdf5dset].value[:,:,option.hdf5slice]
        f.close()        
        
        if option.hdf5readfndata:
            xmin,xmax,ymin,ymax=loadHDF5fnData(fname,option.hdf5dset)
            option.xmin=xmin
            option.xmax=xmax
            option.ymin=ymin
            option.ymax=ymax
            if option.flipHDFylabels:
                option.ymin=ymax
                option.ymax=ymin
                
            print('Read hdf5 function ranges xmin,xmax,ymin,ymax: ',xmin,xmax,ymin,ymax)
            
        
        # convert the ordering of the slice to make first dimention (x) increase rightwards and second dimention (y) increase downwards
        slice=slice.T
#        print np.size(slice)
#        sys.exit()
    elif getFileExtension(fname)=='fits':
        hdulist = pyfits.open(fname)
        print(hdulist.info())
        if option.fitsSlice!=-1:
            tbdata = hdulist[sliceNo].data[option.fitsSlice] # assuming the first extension is a table
        else:
            tbdata = hdulist[sliceNo].data # assuming the first extension is a table
        if type(tbdata)!=type(list()): # assuming this is image
            return tbdata
        else:
            s=array([])
            for c in hdulist[sliceNo].columns.names:
                s=np.append(s, array(tbdata.field(c)))
            slice=transpose(s.reshape(len(hdulist[sliceNo].columns.names),-1))
    else:
        slice=np.loadtxt(fname)

    if option.rr!="":
        removeRows=option.rr.split(",")
#        print removeRows
        rr=0;
        for r in removeRows:
            rr=int(r)
            slice=delete(slice, s_[rr:rr+1],axis=0)
    print()
    print("    * data size after loading and trimming: %i rows %i cols" % (len(slice[:,0]),len(slice[0])))
    print()
    return slice
    

def print_mat_stats(slice):
#     slice=np.array(slice)
    print('size: {}'.format(slice.size))
    print('Sum of diagomal terms: {}'.format(np.sum(slice.diagonal())))
#     print('1^T x C x 1/n^2: {}'.format(np.mean(slice.diagonal())))
    print('Sum of off-diagomal terms: {}'.format(np.sum(slice)-np.sum(slice.diagonal())))
    print('square root of mean diagomal term: {}'.format(np.sqrt(np.mean(slice.diagonal()))))
    print('sqrt( 1^T x C x 1/n^2 ): {}'.format(
        np.sqrt(np.sum(slice))/slice.diagonal().size))
    print('sqrt( diag(C) x 1/size(diag(C))^2 ): {}'.format(
        np.sqrt(
            np.sum(
                slice.diagonal()
                ))/slice.diagonal().size))

    offdiag_sum=np.sum(slice)-np.sum(slice.diagonal())
    
    print('off-diag dev: {}'.format(
        np.sqrt(offdiag_sum)/slice.diagonal().size))

def binMatrix(slice,dim,binSize):
#    s=arange(20)
#    slice=array([])
#    for i in arange(10):
#        slice=np.append(slice,s,axis=0)
#    slice=transpose(slice.reshape(10,-1))
#    print slice
#    sys.exit()
    if dim==0:
        # perform binning along X direction
        n=len(slice[0])
        Nbin=int(float(n)/binSize)
        binslice=array([])
        for i in arange(Nbin):
            if option.binMethod=='mean':
                tmp=slice[:,binSize*i:binSize*i+binSize].sum(axis=1)/float(binSize)
            elif option.binMethod=='max':
                tmp=np.max(slice[:,binSize*i:binSize*i+binSize],axis=1)
            else:
                print('unknown bin method')
                sys.exit(-1)
#             binslice=np.append(binslice,tmp,axis=1)
            binslice=np.append(binslice,tmp)
        return transpose(binslice.reshape(Nbin,-1))
#        print slice
#        sys.exit()
    if dim==1:
        # perform binning along Y direction
        print("perform binning along Y direction")
#        slice=transpose(binMatrix(transpose(slice),0,binSize))
        n=len(slice[:,0])
        Nbin=int(float(n)/binSize)
        binslice=array([])
        for i in arange(Nbin):
            if option.binMethod=='mean':
                tmp=slice[binSize*i:binSize*i+binSize].sum(axis=0)/float(binSize)
            elif option.binMethod=='max':
                tmp=np.max(slice[binSize*i:binSize*i+binSize],axis=0)
            else:
                print('unknown bin method')
                sys.exit(-1)
#             binslice=np.append(binslice,tmp,axis=1)
            binslice=np.append(binslice,tmp)
        return binslice.reshape(Nbin,len(slice[0]))
#        print slice
    return slice


def processSlice(boxSlice, idx=0):
    print('caclulating matrix statistics before processing')
    getMatInfo(boxSlice)
    
    
    # perform slice cut before anything else
    if len(option.cut)==4:
        cut=option.cut
        boxSlice=boxSlice[cut[2]:cut[3],cut[0]:cut[1]]


    if option.scaleZ!=1:
        boxSlice=boxSlice*scaleZ[idx % len(scaleZ)]

    if option.binX>1:
        boxSlice=binMatrix(boxSlice,0,option.binX)
    if option.binY>1:
        boxSlice=binMatrix(boxSlice,1,option.binY)


    if option.log:
        boxSlice=log10(boxSlice)

    if option.abs:
        boxSlice=abs(boxSlice)

    if option.flipY:
        boxSlice=np.flipud(boxSlice)
    if option.flipX:
        boxSlice=np.transpose(boxSlice)
        boxSlice=np.flipud(boxSlice)
        boxSlice=np.transpose(boxSlice)
    if option.transposeData:
        boxSlice=np.transpose(boxSlice)
    if option.hdf5dump and option.hdf5zscramble==False:
#         np.savetxt(fname+'.'+option.hdf5dset+'.'+str(option.hdf5slice), boxSlice);
        np.savetxt(fname+".dump", boxSlice);

    if option.shiftX!=0:
        boxSlice = np.roll(boxSlice,option.shiftX,1)
    if option.shiftY!=0:
        boxSlice = np.roll(boxSlice,option.shiftY,0)

    return boxSlice

###########################################################################################
###########################################################################################
###########################################################################################
# MAIN PROGRAM
###########################################################################################
###########################################################################################
###########################################################################################

formatterX = FuncFormatter(XTicksFormatter)
formatterY = FuncFormatter(YTicksFormatter)

circlePatchDef=list()
if option.circ!='':
    circlePatchDef=cpedsPythCommon.getFloatList(option.circ)

ellipsePatchDef=list()
if option.ellipse!='':
    ellipsePatchDef=cpedsPythCommon.getFloatList(option.ellipse)

rectanglePatchDef=list()
if option.rect!='':
    rectanglePatchDef=cpedsPythCommon.getFloatList(option.rect)

linePatchDef=list()
if option.line!='':
    linePatchDef=cpedsPythCommon.getFloatList(option.line)

if type(option.textxy)!=type(list()):    option.textxy=list()
if type(option.textfxfy)!=type(list()):    option.textfxfy=list()
if type(option.lc)!=type(list()):    option.lc=list(['b','g','r','c','m','y','k','b','g','r','c','m','y','k','b','g','r','c','m','y','k'])
if type(option.gc)!=type(list()):    option.gc=list(['k'])
if type(option.cut)!=type(list()):    option.cut=cpedsPythCommon.getIntList(option.cut)
CM=option.CM.split(',')
contcolor=option.contcolor.split(",")
#print CM
#sys.exit(0)

contours=option.contours.split(':')
scaleZ=cpedsPythCommon.getFloatList(option.scaleZ)
hdf5dset=option.hdf5dset.split(',')

#print contours
#sys.exit()

if len(args)==0:
    if option.st==0 and option.en==0:
        print("Too few parameters given")
        sys.exit(0)
        
boxSliceScrambled=np.array([])
boxSlices=list()


    

matID=0
for fname in args:
    ext=getFileExtension(fname)
    option.hdf5dset=hdf5dset[matID % len(hdf5dset)]
    
    if ext=='hdf5' or ext=='hdf':

        if option.st < option.en:
            for i in arange(option.st,option.en+1):
                fname=args[0]
                boxSlice=loadData(fname,i)
        
                if option.hdf5zscramble:
                    if i==option.st:
                        boxSliceScrambled=boxSlice
                    else:
                        boxSliceScrambled=boxSliceScrambled+boxSlice
            
                    if i==option.en:
                        print()
                        print("Averaging slices")
                        print() 
                        boxSliceScrambled/=(option.en-option.st+1)
                        boxSlice=processSlice(boxSliceScrambled, len(boxSlices))
                        boxSlices.append(boxSlice)
                        if option.hdf5dump:
                            np.savetxt(fname+'.'+option.hdf5dset+'.'+str(option.st)+'-'+str(option.en), boxSlice);
                elif option.hdf5zmaximize:
                    if i==option.st:
                        boxSliceMax=boxSlice
                    else:
                        boxSliceMax=np.max([boxSliceMax,boxSlice],axis=0)
                    if i==option.en:
                        boxSlice=processSlice(boxSliceMax, len(boxSlices))
                        boxSlices.append(boxSlice)
                        if option.hdf5dump:
                            np.savetxt(fname+'.'+option.hdf5dset+'.'+str(option.st)+'-'+str(option.en), boxSlice);
    
        
                else:
                    boxSlice=processSlice(boxSlice, len(boxSlices))
                    boxSlices.append(boxSlice)
        else:
            if option.hdf5slice>=0:
                fname=args[0]
                boxSlice=loadData(fname,option.hdf5slice)
                boxSlice=processSlice(boxSlice, len(boxSlices))
                boxSlices.append(boxSlice)
            
    elif ext=='fits':
        fname=args[0]
        boxSlice=loadData(fname,option.hdu)
        boxSlice=processSlice(boxSlice, len(boxSlices))
        boxSlices.append(boxSlice)
    #    fname=args[0]+'.%04i' % i
    #    boxSlice=loadData(fname)
    #    boxSlices.append(boxSlice)
                    
    else:
        boxSlice=loadData(fname)
        boxSlice=processSlice(boxSlice, len(boxSlices))
        boxSlices.append(boxSlice)

    matID=matID+1
i=0
fig=None
ax=None

for boxSlice in boxSlices:


    print()
    print("Plotting slice: %i" % i)
    print()
    if i==0:
        fig=None
        ax=None
    fig,ax=makeMatrixPlot(boxSlice,i,fname,fig,ax)
    if option.joinPlots==False:
        close(fig)
    print("done")
    print("---------- ")
    print()
    i+=1
            
        

if option.makeMovie:
    print()
    print("Generating movie...")
    print()
    filePrefix=fname+"-"+option.hdf5dset+"-slice_%i" % option.hdf5slice
    print("Using files:")
    os.system('ls -1 %s*.png' % filePrefix)
    cmd="mencoder mf://%s*.png -mf type=png:fps=%i -ovc lavc -lavcopts  vcodec=mpeg4:mbd=1:vbitrate=2800  -oac copy -o %s.avi" % (filePrefix,option.fps, args[0]+"-"+option.hdf5dset)
    os.system(cmd)
    print("done")
    print("---------- ")
    print()
