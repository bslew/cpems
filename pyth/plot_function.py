#!/usr/bin/python2.7
# #!/usr/bin/env python
# -*- Encoding: utf-8 -*-
# coding: utf-8

import sys
import os
import time as tm
sys.path.append('./')
from optparse import OptionParser


##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#
# PROGRAM OPTIONS
#
##########################################################################################
##########################################################################################
##########################################################################################


programDescription = """This is a simple function plotter.
The input text files should be stored as a matrix of numbers.
The chosen columns will be plotted.
"""

parser = OptionParser(description=programDescription)
# parser.add_option("", "--simPrefix", dest="simPrefix", default="", type="string", help='simulation prefix', metavar="PREFIX")

# parser.add_option("-f", "--file", dest="signal_file", default='277_ocraf_absron_010711.dat.totalPower.last500000', 
#                  help="name of the file with input signal. The file should consist of 16 columns with total power signal in order as in the native ocraf format: a1 a2 b1 b2 ...", metavar="FILE")
# parser.add_option("-D", "--boxSize", dest="boxSize", default=256, type="float", help='comoving size of the box [Mpc]', metavar="LENGTH")
parser.add_option("", "--figSize", dest="figSize", default="12,12", type="string", help='size of the figure in inch eg. 19,8. To have A4 plot type A4.' , metavar="STRING")
parser.add_option("", "--figCols", dest="figCols", default='1', type="string", help='''Number of columns in the plot in multiplot mode.
    If a list of integers, the numbers indicate relative sub-plot sizes: e.g. '[1,2]' indicates two subplots with width ratio 1:2 respectively.
    If a single value is required use '[1]'. ''' , metavar="VALUE")
parser.add_option("", "--figRows", dest="figRows", default='1', type="string", help='Number of rows in the plot in multiplot mode' , metavar="VALUE")
parser.add_option("", "--dsPerPlot", dest="dsPerPlot", default=-1, type="int", help='Number of datasets per plot. By default all datasets are plotted in the same plot' , metavar="VALUE")
# parser.add_option("", "--subPlotRatios", dest="subPlotRatios", default='', type="string", help='Definition of relative subplot sizes. By default not used. Eg. for figRows=2 figCols=2 subPlotRations=[2,1],[4,1]' , metavar="VALUE")
parser.add_option("", "--plotDs", dest="plotDs", default='', type="string", help='Coma-separated list of datasets to be plotted on each plot figure. Must sum-up to the total numper of input files' , metavar="VALUE")
parser.add_option("", "--marginFactor", dest="marginFactor", default="0.05", type="float", help='size of the margin as a fraction of the horizontal and vertical extent of the data plotted (DEPRECIATED).' , metavar="STRING")
parser.add_option("", "--marginX", dest="marginX", default="0.0", type="float", help='plot margin X size' , metavar="VALUE")
parser.add_option("", "--marginY", dest="marginY", default="0.0", type="float", help='plot margin Y size' , metavar="VALUE")
parser.add_option("", "--DPIgui", dest="DPIgui", default=70, type="int", help='DPI for plotting on the screen', metavar="NUM")
parser.add_option("", "--DPI", dest="DPI", default=70, type="int", help='DPI for saving to file', metavar="NUM")
parser.add_option("", "--bgcolor", dest="bgcolor", default='white', type="string", help='color to use for the background of the figure', metavar="NUM")
parser.add_option("", "--xticks", dest="xticks", type="float", help='x ticks to be plotted every this cells', action="append")
parser.add_option("", "--xticksMinor", dest="xticksMinor", type="float", help='minor x ticks to be plotted every this cells', action="append")
parser.add_option("", "--yticksMinor", dest="yticksMinor", type="float", help='minor y ticks to be plotted every this cells', action="append")
parser.add_option("", "--yticks", dest="yticks", type="float", help='y ticks to be plotted every this cells', action="append")
parser.add_option("", "--title", dest="title", type="string", help='list of plot titles (for multiplot mode) and single title for single plot mode', metavar="STRING", action="append")
parser.add_option("", "--fontSizeTitle", dest="fontSizeTitle", type="int", help='title font size', metavar="STRING", default=15)
parser.add_option("", "--titleAlignH", dest="titleAlignH", default='center', type="string", help='title alignment (center, left, right)', metavar="STRING")
parser.add_option("", "--titleAlignV", dest="titleAlignV", default='center', type="string", help='title alignment (center, top, baseline)', metavar="STRING")
parser.add_option("", "--xlabel", dest="xlabel", type="string", help='plot x label', metavar="STRING", action="append")
parser.add_option("", "--x2label", dest="x2label", type="string", help='plot x2 label', metavar="STRING", action="append")
parser.add_option("", "--ylabel", dest="ylabel", type="string", help='plot y label', metavar="STRING", action="append")
parser.add_option("", "--fontSize", dest="fontSize", default='20', type="string", help='coma separated list of font size used for ticks per each plot', metavar="VALUE")
parser.add_option("", "--fontSizeCM", dest="fontSizeCM", default='20', type="string", help='coma separated list of font size used for ticks per each plot', metavar="VALUE")
parser.add_option("", "--fontSizeLegend", dest="fontSizeLegend", default='20', type="string", help='coma separated list of font size used for legends per each plot', metavar="VALUE")
parser.add_option("", "--fontSizeLabels", dest="fontSizeLabels", default='20', type="string", help='coma separated list of font size used for labels per each plot', metavar="VALUE")
# parser.add_option("", "--fontSizeXtickLabels", dest="fontSizeXtickLabels", default='18', type="string", help='coma separated list of font size used for XtickLabels per each plot', metavar="VALUE")
# parser.add_option("", "--fontSizeYtickLabels", dest="fontSizeYtickLabels", default='18', type="string", help='coma separated list of font size used for XtickLabels per each plot', metavar="VALUE")
parser.add_option("", "--logYaxes", dest="logYaxes", default='', type="string", help='coma separated list of log,lin combination each for eahc of the subplots to defined log/lin Y axes ', metavar="VALUE")
# parser.add_option("", "--maskBelow", dest="maskBelow", default=np.nan, type="float", help='disable plotting pixels with values below this value ', metavar="VALUE")
parser.add_option("", "--bad_val", dest="bad_val", default='None', type="string", help='disable plotting pixels with these values (map plots only)', metavar="VALUE")
# parser.add_option("", "--set_below", dest="set_below", default='None', type="string", help='set regions with values below VAL to a color eg. "-0.5,#0D0D0D" or "-0.5,0.5" (map plots only)', metavar="VALUE")
parser.add_option("", "--xmin", dest="xmin", type="float", help='minimal X value in the grid to show on X axis', action="append")
parser.add_option("", "--xmax", dest="xmax", type="float", help='maximal X value in the grid to show on X axis', action="append")
parser.add_option("", "--ymin", dest="ymin", type="float", help='minimal Y value in the grid to show on Y axis', action="append")
parser.add_option("", "--ymax", dest="ymax", type="float", help='maximal Y value in the grid to show on Y axis', action="append")
parser.add_option("", "--zmin", dest="zmin", default=-1, type="float", help='minimal Z value in the grid to show on Z axis (3d plots only)', metavar="VALUE")
parser.add_option("", "--zmax", dest="zmax", default=-1, type="float", help='maximal Z value in the grid to show on Z axis (3d plots only)', metavar="VALUE")
parser.add_option("", "--vmin", dest="vmin", type="float", help='Minimal value in the map. Below this value the setUnder option defines the color. If equal to vmax then will be calculated from the data.  (default: 0)', action="append")
parser.add_option("", "--vmax", dest="vmax", type="float", help='Maximal value in the map. Above this value the setAbove option defines the color. If equal to vmax then will be calculated from the data.  (default: 0)', action="append")

parser.add_option("", "--datexmin", dest="datexmin", type="string", help='minimal X value in the grid to show on X axis (only for ts type plots)', action="append")
parser.add_option("", "--datexmax", dest="datexmax", type="string", help='maximal X value in the grid to show on X axis (only for ts type plots)', action="append")
# parser.add_option("", "--xmin", dest="xmin", default=-1, type="float", help='minimal X value in the grid to show on X axis', action="append")
# parser.add_option("", "--xmax", dest="xmax", default=-1, type="float", help='maximal X value in the grid to show on X axis', action="append")
# parser.add_option("", "--ymin", dest="ymin", default=-1, type="float", help='minimal Y value in the grid to show on Y axis', action="append")
# parser.add_option("", "--ymax", dest="ymax", default=-1, type="float", help='maximal Y value in the grid to show on Y axis', action="append")
parser.add_option("", "--Qscale", dest="Qscale", default=1, type="float", help='Quiver plot arrows scale factor', metavar="VALUE")
# parser.add_option("", "--yerrType", dest="yerrType", default="", type="string", help='type of y-error bars to be plotted (default "" not plotted). Possible values are: shaded. Used with "fn" type plots only', metavar="VALUE")
parser.add_option("", "--Bleft", dest="border_left", default=0.1, type="float", help='border size in the plot from the left', metavar="VALUE")
parser.add_option("", "--Bright", dest="border_right", default=0.1, type="float", help='border size in the plot from the right', metavar="VALUE")
parser.add_option("", "--Btop", dest="border_top", default=0.1, type="float", help='border size in the plot from the top', metavar="VALUE")
parser.add_option("", "--Bbottom", dest="border_bottom", default=0.1, type="float", help='border size in the plot from the bottom', metavar="VALUE")
parser.add_option("", "--Bhspace", dest="border_hspace", default=0.2, type="float", help='border horizontal space between the plots for multiplot mode', metavar="VALUE")
parser.add_option("", "--Bvspace", dest="border_vspace", default=0.2, type="float", help='border vertical space between the plots for multiplot mode', metavar="VALUE")
# parser.add_option("", "--interpolation", dest="interpolation", default='nearest', type="string", help='type of interpolation to be used for plotting matrix: the allowed interpolation strings are those compatible with matplotlib eg: "nearest", "bilinear", "bicubic", "spline16", "spline36", "hanning", "hamming", "hermite", "kaiser", "quadric", "catrom", "gaussian", "bessel", "mitchell", "sinc", "lanczos" (default: nearest) ', metavar="STRING")
parser.add_option("-o", "--outputFile", dest="outputFile", default='', type="string", help='name of the output file if --save option is used (default: "") ', metavar="STRING")
# parser.add_option("", "--colx", dest="colx", default="0", type="string", help='comma separated string of columns to be plotted at x axis for each of the input files ', metavar="STRING", action="append")
# parser.add_option("", "--coly", dest="coly", default="1", type="string", help='comma separated string of columns to be plotted at y axis for each of the input files ', metavar="STRING", action="append")
parser.add_option("-x", "--colx", dest="colx", type="int", help='number of the column to be plotted at x axis for each of the input files', action="append")
parser.add_option("-y", "--coly", dest="coly", type="int", help='number of the column to be plotted at y axis for each of the input files', action="append")
parser.add_option("-z", "--colz", dest="colz", type="int", help='number of the column to be plotted at z axis for each of the input files', action="append")
parser.add_option("", "--badY", dest="badY", type="string", help='bad Y value that will be ignored. None for all data OK.', action="append")
parser.add_option("", "--zorder", dest="zorder", type="int", help='z-order for the plot. By default it is generated according to the plotting order', action="append")
parser.add_option("", "--x2", dest="x2", type="int", help='number of the column with data to be used with fillbtw type plots (eg -x 0 --x2 1 -y 2) ', action="append")
parser.add_option("", "--y2", dest="y2", type="int", help='number of the column with data to be used with fillbtw type plots (eg -x 0 -y 1 --y2 2)', action="append")
parser.add_option("", "--yerr", dest="colyerr", type="int", help='number of the column to be used for "fnshaded" type plots. This data is used as lower CL limits if yerr2 is also specified', action="append")
parser.add_option("", "--yerr2", dest="colyerr2", type="int", help='number of the column to be used for "fnshaded" type plots. This data is used as upper CL limits', action="append")
parser.add_option("", "--xerr", dest="colxerr", type="int", help='number of the column to be used for error type plots. This data is used as lower CL limits if xerr2 is also specified', action="append")
parser.add_option("", "--xerr2", dest="colxerr2", type="int", help='number of the column to be used for error type plots. This data is used as upper CL limits', action="append")
parser.add_option("-l", "--label", dest="label", type="string", help='label for the curve for each of the input files', action="append")
# parser.add_option("-w", "--width", dest="width",  help='line width for the curve for each of the input files', action="append_const", const=1.0)
parser.add_option("-w", "--width", dest="width", type="float", help='line width for the curve for each of the input files', action="append")
# parser.add_option("", "--pt", dest="pt", default="", type="string", help='marker for the data points', metavar="STRING", action="append")
# parser.add_option("", "--ls", dest="ls", default="-", type="string", help='line style: - -- -. :', metavar="STRING", action="append")
parser.add_option("", "--pt", dest="pt", type="string", help='marker for the data points', metavar="STRING", action="append")
parser.add_option("", "--ps", dest="ps", type="int", help='marker size for plots for each dataset', metavar="NUM", action="append")
parser.add_option("", "--mew", dest="markerEdgeWidth", type="int", help='marker edge width', action="append")
parser.add_option("", "--ls", dest="ls", type="string", help='line style to be used for plotting dataset. (default order is: -,-,-,-,-,-,-,--,--,--,--,--,--,--,-.,-.,-.,-.,-.,-.,-.) :', metavar="STRING", action="append")
parser.add_option("", "--lsdash", dest="lsdash", type="string", help='''line dashes styles to be used for plotting dataset. This applies only to the dataset that is plotted 
with - or -. line style (default order is: [10,20]).  To combine multiple styles use eg. [10,20],[10,10]''', metavar="STRING", action="append")
parser.add_option("", "--lc", dest="lc", type="string", help='line color to be used for plotting dataset. (default order is: b,g,r,c,m,y,k,b,g,r,c,m,y,k,b,g,r,c,m,y,k)', metavar="STRING", action="append")
parser.add_option("", "--pc", dest="pc", type="string", help='point color for each dataset', metavar="STRING", action="append")
parser.add_option("", "--fc", dest="fc", type="string", help='fill color to be used for filling rectangles and between lines. (default order is: b,g,r,c,m,y,k,b,g,r,c,m,y,k,b,g,r,c,m,y,k)', metavar="STRING", action="append")
parser.add_option("", "--ec", dest="ec", type="string", help='edge color to be used for plotting filled datasets. (default order is: b,g,r,c,m,y,k,b,g,r,c,m,y,k,b,g,r,c,m,y,k)', metavar="STRING", action="append")
parser.add_option("", "--fhs", dest="fhs", type="string", help='face hatch style. (possible values are:  /\|-+xoO.*)', metavar="STRING", action="append")
# parser.add_option("", "--olc", dest="olc",  type="string", help='overlie line color to be used for plotting horizontal or vertical lines in the plot. (default order is: b,g,r,c,m,y,k)', metavar="STRING", default='b,g,r,c,m,y,k')
parser.add_option("", "--olc", dest="olc", type="string", help='overlie line color to be used for plotting horizontal or vertical lines in the plot. (default order is: b,g,r,c,m,y,k)', metavar="STRING", action="append")
# parser.add_option("", "--ols", dest="ols",  type="string", help='overlie line style to be used for plotting horizontal or vertical lines in the plot. (default order is: -)', metavar="STRING", default='-')
parser.add_option("", "--ols", dest="ols", type="string", help='overlie line style to be used for plotting horizontal or vertical lines in the plot. A comma separated list. Use once per each subplot (default order is: -)', metavar="STRING", action="append")
parser.add_option("", "--ofc", dest="ofc", type="string", help='overlie fill color to be used for plotting on VHspan plots. (default order is: -)', metavar="STRING", action="append")
parser.add_option("", "--olw", dest="olw", type="string", help='overlie line width to be used for plotting on VHline plots. (default order is: -)', metavar="STRING", action="append")
parser.add_option("", "--oalpha", dest="oalpha", type="string", help='overlie alpha list. (default order is: 1)', metavar="STRING", default='1')
parser.add_option("", "--hdu", dest="hdu", type="int", help='HDU index in the fits file of which array is to be plotted', action="append")
parser.add_option("", "--wh", dest="where", type="string", help='dataset filter', metavar="NUM", action="append")
parser.add_option("", "--operY", dest="operY", type="string", default='', help='Comma separated list of operations to be performed on the input data Y before plotting. Use with valsY option.', metavar="STRING")
parser.add_option("", "--operX", dest="operX", type="string", default='', help='Comma separated list of operations to be performed on the input data X before plotting. Use with valsX option.', metavar="STRING")
parser.add_option("", "--valsY", dest="valueY", type="string", default='', help='Comma separated list of values to be used in Yoperations, for each intput file dataset.', metavar="STRING")
parser.add_option("", "--valsX", dest="valueX", type="string", default='', help='Comma separated list of values to be used in Xoperations, for each intput file dataset.', metavar="STRING")
parser.add_option("", "--Hline", dest="Hline", type="string", help='Comma separated list of values to mark horizontal lines on the plot', metavar="STRING", action="append")
parser.add_option("", "--Vline", dest="Vline", type="string", help='Comma separated list of values to mark vertical lines on the plot', metavar="STRING", action="append")
parser.add_option("", "--Vspan", dest="Vspan", type="string", help='Comma separated list of value-pairs to mark vertical bars on the plot. Each Vspan option corresponds to the number of plot in which it will be marked', metavar="STRING", action="append")
parser.add_option("", "--Hspan", dest="Hspan", type="string", help='Comma separated list of value-pairs to mark horizontal bars on the plot. Each Hspan option corresponds to the number of plot in which it will be marked', metavar="STRING", action="append")
parser.add_option("", "--plotVerticalLinesWithLabelsFromFile", dest="plotVerticalLinesWithLabelsFromFile", default="", type="string", help='indicates that you want to overplot the vertical lines using the data from file. In the first column should be the number at which the vertical line should be placed and the rest of the row is used as label for that line.', metavar="FILE")
parser.add_option("", "--plotLabelsFromFile", dest="plotLabelsFromFile", default="", type="string", help='indicates that you want to overplot points with labels loaded from a file. In the first two columns should be x,y coordinates and the third column should be a string used as label for that point.', metavar="FILE")
parser.add_option("", "--plotSpecial", dest="plotSpecial", default=-1, type="int", help='number of the special block to be called', metavar="NUM")
parser.add_option("", "--loadMask", dest="loadMask", default="", type="string", help='Loads mask as defined in this file. Each row specifies a range (from, to) to be marked on the plot. The mask data are loaded into global maskRanges array, so that it can be edit in the interactive mode', metavar="FILE")
parser.add_option("", "--startFrom", dest="startFrom", default=0, type="int", help='starts loading from indicated row', metavar="NUM")
parser.add_option("", "--rows", dest="rows", type="int", help='how many rows to read', metavar="NUM", action="append")
parser.add_option("", "--every", dest="every", type="int", help='plot every this data row', metavar="NUM", action="append")
parser.add_option("", "--bin", dest="bin", type="int", help='plot binned data', metavar="NUM", action="append")
# parser.add_option("-c", "--colColor", dest="colColor", default=-1, type="int", help='number of the column to be used to plot in color scale', metavar="NUM")
parser.add_option("-c", "--colColor", dest="colColor", type="int", help='number of the column to be used to plot in color scale', action="append")
parser.add_option("-s", "--colSize", dest="colSize", default=-1, type="int", help='number of the column to be used to plot data with markers size proportional to this column values', metavar="NUM")
parser.add_option("", "--legendLoc", dest="legendLoc", default='ur', type="string", help='default location of the plot legend: ur,ul, br,bl ("None" if not wanted in given panel)', metavar="STRING")
parser.add_option("", "--legendPatternLength", dest="legendPatternLength", default=3, type="float", help='default length of lines pattern that show up in legend (default: 3)', metavar="VALUE")
parser.add_option("", "--legendNcol", dest="legendNcol", default=1, type="int", help='default number of columns in plot legend', metavar="INT")
parser.add_option("", "--legendMode", dest="legendMode", default='default', type="string", help='default legend mode. Can be: default|expand', metavar="INT")
parser.add_option("", "--legendBorderAxesPad", dest="legendBorderAxesPad", default=0.5, type="float", help='default plot legend border padding from the plot frame', metavar="VALUE")
parser.add_option("", "--legendSplit", dest="legendSplit", default='', type="string", help='''Controls how/if the legend is split. 
    If empty, then the legend is in a single box. If eg. 2,3 then there will be two legends first containing two entries and the second three entries. 
    The --legendLoc should specify where the legends should be located using a coma-separated list eg. ul,br NOT IMPLEMENTED YET ''', metavar="INT")
# parser.add_option("", "--legendParams", dest="legendParams", default='', type="string",  help='''Extra legend parameters. 
#     Eg. 'bbox_to_anchor=(0., 1.02, 1., .102),ncol=2, mode="expand", borderaxespad=0.0'. ''', metavar="STRING")
parser.add_option("-F", "--format", dest="fileFormat", type="string", help='''list of input file formats according to which loading is done. 
It is overridden by global --fits modifier (default: txt) 
The available formats are: [txt,txtX,fits,fitsPL,fitsWMAP, hdf5]
txt - [commented] text file regular array
txtCol=X - text file/pipe/socket with white space separated row of numbers that is converted to a 2-d regular array using the X columns
fits - fits file
fitsPL - ?? what was that ?
fitsWMAP - WMAP maps fits files
hdf5 - hdf5 file reader for 2d and 3d grid data plotting''', action="append")
parser.add_option("", "--hdf5dset", dest="hdf5dset", default='', type="string", help='comma separated list of names of the datasets from the hdf5 file (default: "") ', metavar="STRING")
parser.add_option("", "--hdf5slice", dest="hdf5slice", default=0, type="int", help='for 3d grid data the slice to use along 3rd dimension (default: 0) ', metavar="INT")
parser.add_option("", "--histNbin", dest="histNbin", type="int", default=20, help="number of bins in histogram (only for -t hist) (default: 20)")
parser.add_option("", "--gridNbins", dest="gridNbins", type="string", default='25,25', help="number of bins when calculating points density using a histogram (only for -t scatDensCont and scatDensContFill plots) (default: 20,20)")
parser.add_option("", "--mapPts", dest="mapPts", type="int", default=600, help="number of points in the 'map' type plot. This number is passed to draw_maps program")
parser.add_option("", "--circ", dest="circ", type="string", default='', help='circles to plot defined as x,y,r. eg. --circ 10,20,5 for circle at x=10 y=20 and radii=5 for fn and scat type plots. This creates a new patch on a plot.', metavar="STRING")
parser.add_option("", "--rect", dest="rect", type="string", help='rectangles to plot defined as x,y,dx,dy,angle where x,y is the center point and dx,dy width,height. eg. --rect 10,20,5,6,0 for a rectangle at x=10 y=20 and width=5 and height=6. .', action="append", metavar="STRING")
parser.add_option("", "--gridColor", dest="gridColor", type="string", default='0.5', help='grid color', metavar="STRING")
parser.add_option("", "--dateFmt", dest="dateFmt", type="string", default='', help='date format for ts type plots: eg "%Y-%m-%d %H:%M:%S"', metavar="STRING")
parser.add_option("", "--dateFmtPlot", dest="dateFmtPlot", type="string", default='%Y-%m-%d %H:%M:%S', help='date format for ts type plots used for plotting: eg "%Y-%m-%d %H:%M:%S"', metavar="STRING")
parser.add_option("", "--Rxlabels", dest="Rxlabels", default=0, type="float", help='rotate xticklabels by this angle [deg]', metavar="VALUE")
parser.add_option("", "--Rylabels", dest="Rylabels", default=0, type="float", help='rotate yticklabels by this angle [deg]', metavar="VALUE")
parser.add_option("", "--removeXtickLabels", dest="removeXtickLabels", type="int", help='if 0 - not used; if single value - plot every this tick label; if coma-separated list - indexes of labels to be removed (NOT implemented yet). Useful for time series plots', action="append", metavar="VALUE")
parser.add_option("", "--removeYtickLabels", dest="removeYtickLabels", type="int", help='if 0 - not used; if single value - plot every this tick label; if coma-separated list - indexes of labels to be removed (NOT implemented yet). Useful for time series plots', action="append", metavar="VALUE")

# parser.add_option("", "--rect", dest="rect", type="string",  default='', help='rectangles to plot defined as x,y,dx,dy where x,y is the center point and dx,dy width,height. eg. --rect 10,20,5,6 for a rectangle at x=10 y=20 and width=5 and height=6. .', metavar="STRING")
# parser.add_option("", "--mk11line", dest="mk11line", type="string",  default='', help='generates line a=1 data to plot.   ', metavar="STRING")
parser.add_option("", "--textxy", dest="textxy", type="string", help='print text on the plot. The text string should be: x,y,"text to be plotted"[,rotation[,size]]', metavar="NUM", action="append")
parser.add_option("", "--textfxfy", dest="textfxfy", type="string", help='print text on the plot. The text string should be: x,y,"text to be plotted"[,rotation[,size]]', metavar="NUM", action="append")
parser.add_option("", "--annotate", dest="annotate", type="string", help='''annotate the plot with text and arrows. 
One option for For each sub-plot. For each sub-plot a string should be a colon-separated list of annotation definitions.
Each definition should consist of a comma-separated list with: x,y,xt,yt,"annotateion text"[,rotation[,size,[color,[arrows]]]]''', metavar="NUM", action="append")

# data generating options
parser.add_option("", "--plotCircle", dest="plotCircle", type="string", help='circle to plot defined as l,b,r. eg. --plotCircle 10,20,5 for circle at l=10 b=20 and radii=5 deg. This is generates a new dataset for plotting', metavar="STRING", action="append")


# spherical plots options
parser.add_option("-t", "--ptype", dest="plotType", type="string", help='''plot types for the dataset. The available plot types are: 
[fn,fnshaded, err, scat, scatContFill, scatContFillD, scatDensContFill, scatDensCont, hist, map,sphere, img, PR, circ,vect,ts,fillbtwX,fillbtwY] (default: fn)''', action="append")
parser.add_option("", "--proj", dest="proj", default="moll", type="string", help='projection type for --sphere mode plotting', metavar="STRING")
parser.add_option("", "--merid", dest="meridians", default="30", type="string", help='comma-separated list of meridians to plot. If single value is given then it is interpreted as a meridians separation [deg]', metavar="STRING")
parser.add_option("", "--meridLabels", dest="meridianLabels", default="", type="string", help='comma-separated list of meridian labels to plot.', metavar="STRING")
parser.add_option("", "--para", dest="parallels", default="15", type="string", help='comma-separated list of parallels to plot. If single value is given then it is interpreted as a parallels separation [deg]', metavar="STRING")
parser.add_option("", "--MPfontSize", dest="MPfontSize", default=15, type="int", help='font size of meridians and parallels labels (default: 15)', metavar="NUM")
parser.add_option("", "--extraTextFontSize", dest="extraTextFontSize", default=15, type="int", help='extra texts font size (default: 15)', metavar="NUM")
parser.add_option("", "--MPcolor", dest="MPcolor", default="k", type="string", help='Color the be used for plotting meridians and parallels ("map" type plots only) (default: k)', metavar="STRING")
parser.add_option("", "--lon0", dest="lon0", default=0, type="float", help='zero longitude for Mollweide/ortho projection plot. Use only for "sphere" plot tyle. (default: 0)', metavar="NUM")
parser.add_option("", "--lat0", dest="lat0", default=90, type="float", help='zero longitude for ortho projection plot. Use only for "sphere" plot tyle. (default: 90)', metavar="NUM")
parser.add_option("", "--lon1", dest="lon1", default=45, type="float", help='1st longitude for lcc projection plot. Also used in stere projs. to define the corners of the plot. Use only for "sphere" plot tyle. (default: 0)', metavar="NUM")
parser.add_option("", "--lat1", dest="lat1", default=10, type="float", help='1st longitude for lcc projection plot. Also used in stere projs. to define the corners of the plot. Use only for "sphere" plot tyle. (default: 0)', metavar="NUM")
parser.add_option("", "--lon2", dest="lon2", default=235, type="float", help='2nd longitude for lcc projection plot. Also used in stere projs. to define the corners of the plot.  Use only for "sphere" plot tyle. (default: 10)', metavar="NUM")
parser.add_option("", "--lat2", dest="lat2", default=10, type="float", help='2nd longitude for lcc projection plot. Also used in stere projs. to define the corners of the plot. Use only for "sphere" plot tyle. (default: 10)', metavar="NUM")
parser.add_option("", "--mapwidth", dest="mapwidth", default=500, type="float", help='Width in km of the projected map. Use only for "sphere" plot type. (default: 500)', metavar="NUM")
parser.add_option("", "--mapheight", dest="mapheight", default=500, type="float", help='Height in km of the projected map. Use only for "sphere" plot type. (default: 500)', metavar="NUM")
parser.add_option("", "--mapresolution", dest="mapresolution", default='c', type="string", help='Map country/ocean contours resolution. Use only for "sphere" plot type. (default: c)', metavar="CHAR")
parser.add_option("", "--levels", dest="levels", type="string", action="append", help='Number of levels for colorbar in the map plot, or list of comma-separated list of values to be used for lelvels ("map", scatDensCont* scatCont* type plots only) (default: 50)', metavar="LIST")
# parser.add_option("", "--vmin", dest="vmin", default=0, type="float", help='Minimal value in the map. Below this value the setUnder option defines the color. If equal to vmax then will be calculated from the data.  ("map" type plots only) (default: 0)', metavar="NUM")
# parser.add_option("", "--vmax", dest="vmax", default=0, type="float", help='Maximal value in the map. Above this value the setAbove option defines the color. If equal to vmax then will be calculated from the data.  ("map" type plots only) (default: 0)', metavar="NUM")
parser.add_option("", "--setAbove", dest="setAbove", default="r", type="string", help='Color the be used in the map type plot for values > vmax.  ("map" type plots only) (default: r)', metavar="STRING")
parser.add_option("", "--setBelow", dest="setBelow", default="b", type="string", help='Color the be used in the map type plot for values < vmin.  ("map" type plots only) (default: b)', metavar="STRING")
parser.add_option("", "--colorbar", action="store_true", dest="colorbar", default=False, help='triggers showing the colorbar on plot ("map" type plots only) ')
parser.add_option("", "--palette", dest="palette", type="string", help='name of the matplotlib palette to use. (jet, gray)', metavar="STRING", action="append")
parser.add_option("", "--IMGextent", dest="IMGextent", default=",", type="string", help='width,height of the picture to plot, given to math the axes on the other plots (need to be set by hand, for now).', metavar="STRING")
parser.add_option("", "--IMGalpha", dest="IMGalpha", default=0.5, type="float", help='alpha value for the background image for the img type plots. Only one value is accepted.', metavar="VALUE")
# parser.add_option("", "--IMGflipX", action="store_true", dest="IMGflipX", default=False, help='Controls image flips about X axis.')
parser.add_option("", "--IMGflipY", action="store_true", dest="IMGflipY", default=False, help='Controls image flips about Y axis.')
parser.add_option("", "--sphereContour", action="store_true", dest="sphereContour", default=False, help='use fontour plot for sphere type plots insetad of pcolormesh ("sphere" type plots only) ')
parser.add_option("", "--draw_maps_path", dest="draw_maps_path", default="/home/blew/programy/Mscs.devel/Mscs/bin/", type="string", help='directory where the draw_maps program should be searched for. (useful for php script executions)', metavar="STRING")
parser.add_option("", "--draw_maps_options", dest="draw_maps_options", default="-n 256", type="string", help='string with options for draw_maps program for plotting circles (default: -n 256)', metavar="STRING")
parser.add_option("", "--coastLines", action="store_true", dest="coastLines", default=False, help="triggers showing the coastLines on plot")
parser.add_option("", "--countryLines", action="store_true", dest="countryLines", default=False, help="triggers showing country on plot")
parser.add_option("", "--riverLines", action="store_true", dest="riverLines", default=False, help="triggers showing rivers on plot")
# parser.add_option("", "--zlabel", dest="zlabel", default='z', type="string", help='plot z label used for some spherical plots for colorbar descriptions', metavar="STRING")
parser.add_option("", "--zlabel", dest="zlabel", type="string", help='plot z label', metavar="STRING", action="append")
parser.add_option("", "--reverseLon", action="store_true", dest="reverseLon", default=False, help="triggers plotting longituges increasing leftwards")
parser.add_option("", "--cbLabelsNum", dest="cbLabelsNum", type='float', default=7, help='''This option specifies the number of colorbar labels (default: 7). ''')
parser.add_option("", "--ZticksFmt", dest="ZticksFmt", type='string', default='%.2f', help="format for the Zticklabels", metavar="STRING")

#parser.add_option("", "--CM", dest="CM", type="string", default='hot', help='''color palette to be used with 'scat' type plots. Interesting choices are 
#    hot, spectral, jet, afmhot, jet, none (default: "hot")''', metavar="STRING")
parser.add_option("", "--CM", dest="CM", type="string", help='''color palette to be used with 'scat' type plots. Interesting choices are 
    hot, spectral, jet, afmhot, jet, none (default: "hot"). 
    Currently one CMap is available per sub-plot.''', metavar="STRING", action="append")
parser.add_option("", "--CMrange", dest="CMrange", type='string', default='', help='''This option specifies 
    the range of the colormap that will be used for color mapping. E.g. '0.0,0.9 will clip the  colors in the colormap CM that 
    map to values between 0.9 and 1.0 (default: don't clip)''')


#
# switches
#
parser.add_option("", "--autolimits", action="store_true", dest="autolimits", default=False, help="recalculate limits automatically")
parser.add_option("", "--brokenYaxis", action="store_true", dest="brokenYaxis", default=False, help="triggers showing minor ticks")
parser.add_option("", "--minorTicks", action="store_true", dest="minorTicks", default=False, help="triggers showing minor ticks")
parser.add_option("", "--minutes", action="store_true", dest="minutes", default=False, help="triggers showing ticks every minute")
parser.add_option("", "--hours", action="store_true", dest="hours", default=False, help="triggers showing ticks every hour")
parser.add_option("", "--days", action="store_true", dest="days", default=False, help="triggers showing ticks every day")
parser.add_option("", "--months", action="store_true", dest="months", default=False, help="triggers showing ticks every month")
parser.add_option("", "--years", action="store_true", dest="years", default=False, help="triggers showing ticks every year")
parser.add_option("", "--grid", action="store_true", dest="grid", default=False, help="triggers showing the grid on plot")
parser.add_option("", "--save", action="store_true", dest="save", default=False, help="triggers saving the plot a image")
parser.add_option("", "--logX", action="store_true", dest="logX", default=False, help="triggers log X scale")
parser.add_option("", "--logY", action="store_true", dest="logY", default=False, help="triggers log Y scale")
parser.add_option("", "--logZ", action="store_true", dest="logZ", default=False, help="triggers log Z scale for colorbar")
parser.add_option("", "--linZtickLabels", action="store_true", dest="linZtickLabels", default=False, help='''to be used with --logZ option in colorbar plots. Converts z ticklabels back to the linear scale, 
    while preserving the linear color scale for logarithms of color-mapped values. The scale becomes non-linear. ''')

parser.add_option("", "--transpose", action="store_true", dest="transpose", default=False, help="transpose data array right after loading")
parser.add_option("", "--yerr2sided", action="store_true", dest="yerr2sided", default=False, help="if yerr column is loaded, this switch tells to plot shaded region assuming that the yerr data is already twosided error bar: y-yerr/2, y+yerr/2")
parser.add_option("", "--fits", action="store_true", dest="fits", default=False, help="indicates that the input file is in fits format")
parser.add_option("", "--PR", action="store_true", dest="PR", default=False, help="indicates that the input file a POVray DF3 file format array")
parser.add_option("", "--binDAQd", action="store_true", dest="binDAQd", default=False, help="indicates that the input file is in binDAQd format - works only with --big option. -y 2 indicates the first data channel. -y 10 - the switching signal")
parser.add_option("", "--sqrt", action="store_true", dest="sqrt", default=False, help="triggers plotting sqrt of the signal (useful for power spectra)")
parser.add_option("", "--absY", action="store_true", dest="absY", default=False, help="triggers plotting absolute value of the signal on y-axis. Performed after log and before sqrt and shift0")
parser.add_option("", "--absX", action="store_true", dest="absX", default=False, help="triggers plotting absolute value of the signal on x-axis. Performed after log and before sqrt and shift0")
parser.add_option("", "--shift0", action="store_true", dest="shift0", default=False, help="will shift signals to zero mean before plotting")
parser.add_option("", "--shift0X", action="store_true", dest="shift0X", default=False, help="will shift X values to start with zero before plotting")
parser.add_option("", "--shift0Xmax", action="store_true", dest="shift0Xmax", default=False, help="will shift X values so that maximal X will be zero")
parser.add_option("", "--noaxes", action="store_true", dest="noAxes", default=False, help="switches axes off.")
parser.add_option("", "--axesSwitch", action="append", dest="axesSwitch", type='int', help="switches axes on or off on each of the subplots.")
parser.add_option("", "--interactive", action="store_true", dest="interactive", default=False, help="enters interactive mode allowing data selection from the plot.")
# parser.add_option("", "--inter", action="store_true", dest="inter", default=False, help="enters a new interactive mode allowing data selection from the plot.")
parser.add_option("", "--big", action="store_true", dest="big", default=False, help="flag to indicate that we are loading a big file. This will optimize the loading and plotting")
parser.add_option("", "--asCol", action="store_true", dest="asCol", default=False, help="flag to indicate that we are loading a column data. Normally the load does not recognise verctor between column or row alignment, so this option helps to define which is it.")
parser.add_option("", "--scriptMode", action="store_true", dest="scriptMode", default=False, help="flag to indicate that we are working in script mode for generating pictures and without GUI.")
parser.add_option("", "--transparent", action="store_true", dest="transparent", default=False, help="flag to indicate that we want transparent background in the output figure .")
parser.add_option("", "--nolog", action="store_true", dest="nolog", default=False, help="blocks logging the run to a log file: .plot_function.log")
parser.add_option("", "--mkLabels", action="store_true", dest="mkLabels", default=False, help='make labels for each of the input files based on the file names')
parser.add_option("", "--polar", action="store_true", dest="polar", default=False, help="triggers plotting x data as phi and y data as R on polar plot")
parser.add_option("", "--proj3d", action="store_true", dest="proj3d", default=False, help="triggers plotting on 3d projection (scat3d plot only)")
parser.add_option("", "--chColorsLTypes", action="store_true", dest="chColorsLTypes", default=False, help="this will trigger changing both colors and line types for each dataset - useful for publications in gray scale to distinguish the datasets even if the colors are not displayed correctly")

parser.add_option("", "--histCuml", action="store_true", dest="histCuml", default=False, help='make histogram cumulative')
parser.add_option("", "--histNormed", action="store_true", dest="histNormed", default=False, help='make histogram normed')
parser.add_option("", "--scatDensUsePlotLimits", action="store_true", dest="scatDensUsePlotLimits", default=False, help='when making a density plot, use the --xmin --xmax and --ymin --ymax values instead of the data ranges to define the grid used for density calculations')


parser.add_option("", "--nestedGrid", dest="nestedGridLevel", default=-1, type="int", help='plots a nested grid over the plotted data with level of subdivision specified by this parameter (default: -1 = no nested grid). 0 = generates rectangle around the data, 1 - first square-tree scheme grid and so on ', metavar="NUM")
parser.add_option("", "--nestedGridOct", action="store_true", dest="nestedGridOct", default=False, help="the nested grid is an OCT tree type grid. If false, then any number of gric cells is possible")
parser.add_option("", "--nestedGridDset", dest="nestedGridDset", default=-1, type="int", help='plots the nested grid only for the specified dataset. By detault all datasets have their nested grid plotted. (default: -1). Value 0 makes nested grid plot only for the first dataset. Datasets are enumerated from 0.', metavar="NUM")
parser.add_option("", "--saveData", action="store_true", dest="saveData", default=False, help="triggers saving input data to plot_function.dumpFile. Usefult for single file plotting eg. saved in fits format ")
parser.add_option("", "--mk11line", action="store_true", dest="mk11line", default=False, help="triggers generating a=1 line data to plot. The dataset plot details will be taked as last... expplain this better ")
parser.add_option("", "--equalAspect", action="store_true", dest="equalAspect", default=False, help="triggers equal aspect ratio axes ")
parser.add_option("", "--shareX", action="store_true", dest="shareX", default=False, help="triggers sharing X axis in multi-subplot plots")
parser.add_option("", "--show", action="store_true", dest="show", default=False, help="triggers showing the plot even if it is saved to file. ")

parser.add_option("", "--cont", action="store_true", dest="continuousPotting", default=False, help="triggers continuous data load/plot sequence")
parser.add_option("", "--contAutoRescale", action="store_true", dest="contAutoRescale", default=False, help="continuousPotting option: rescales the plot axes according to the most recent data span")
parser.add_option("", "--clear", action="store_true", dest="clearFig", default=False, help="continuousPotting option: clears the plot after every update")
parser.add_option("", "--waitN", dest="waitN", default="1", type="float", help='continuousPotting option: wait N seconds before replotting', metavar="VALUE")

# parser.add_option("", "--zrangeSym", action="store_true", dest="zrangeSym", default=False, help="triggers showing z-range to be within the -max(abs(zmin),abs(zmax)) and +max(abs(zmin),abs(zmax)) ")
# parser.add_option("", "--enumSlices", action="store_true", dest="enumSlices", default=False, help="triggers enumerating slices in the plot title by the slice number")
# parser.add_option("", "--makeMovie", action="store_true", dest="makeMovie", default=False, help="triggers making movie from the slices")
# parser.add_option("", "--maskBad", action="store_true", dest="maskBad", default=False, help="triggers masking bad pixels")
# parser.add_option("", "--plotSpecial", action="store_true", dest="plotSpecial", default=False, help="plot special block")


# movie options
# parser.add_option("", "--st", dest="st", default=0, type="int", help='first slice', metavar="NUM")
# parser.add_option("", "--en", dest="en", default=0, type="int", help='last slice', metavar="NUM")
# parser.add_option("", "--fps", dest="fps", default=10, type="int", help='frames per second', metavar="NUM")


# POVRAY options
parser.add_option("-N", "--fieldPixNum", dest="fieldPixNum", default="200", type="string", help='comma separated list of numbers defining the field pixel number in each direction. eg. 200 for 200x200x200, or 100,200,300 for Nx=100, Ny=200 and Nz=300', metavar="STRING")
parser.add_option("", "--PRfov", dest="PRfov", default="80", type="float", help='field of view (80) [deg]', metavar="VALUE")
parser.add_option("", "--PRwidth", dest="PRwidth", default="600", type="float", help='picture width [pixels]', metavar="VALUE")
parser.add_option("", "--PRheight", dest="PRheight", default="600", type="float", help='picture height [pixels]', metavar="VALUE")
parser.add_option("", "--PRcamX", dest="PRcamX", default="0", type="float", help='camera location X coordinate [pixels]', metavar="VALUE")
parser.add_option("", "--PRcamY", dest="PRcamY", default="0", type="float", help='camera location Y coordinate [pixels]', metavar="VALUE")
parser.add_option("", "--PRcamZ", dest="PRcamZ", default="0", type="float", help='camera location Z coordinate [pixels]', metavar="VALUE")
parser.add_option("", "--PRgui", action="store_true", dest="PRgui", default=False, help="triggers displaying rendering progress in window")
parser.add_option("", "--PRhold", action="store_true", dest="PRhold", default=False, help="triggers pause when rendering is done")

# MAYAVI options
parser.add_option("", "--mayaPS", dest="mayaPS", default="1", type="float", help='maya point size', metavar="VAL")
parser.add_option("", "--Mvrange", dest="Mvrange", default="0,0", type="string", help='maya volume vmin,vmax parameters', metavar="STRING")
parser.add_option("", "--MdataScale", dest="MdataScale", default=1, type="float", help='scale factor applied to the data before plotting', metavar="STRING")
parser.add_option("", "--Mcontours", dest="Mcontours", default="", type="string", help='coma separated list of contours to be plotted with mayaContour type of plot', metavar="STRING")



(option, args) = parser.parse_args()

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#
# IMPORT REQUIRED PACKAGES
#
##########################################################################################
##########################################################################################
##########################################################################################

if option.scriptMode:
    import matplotlib
    matplotlib.use('Agg')

# if option.log:
import time
from datetime import datetime, date, time
# from datetime import *
import ast

from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab
from matplotlib.ticker import FuncFormatter, MultipleLocator, FormatStrFormatter
# from matplotlib.image import AxesImage

from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection

import matplotlib.path as mpath
import matplotlib.patches as mpatches
# if option.subPlotRations!='':
# if type(option.figRows)==type(list()) and type(option.figCols)==type(list()):
import matplotlib.gridspec as gridspec
import csv
import h5py
import pyfits
import re



# imports due for data selectors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.widgets import Button
from matplotlib.widgets import Cursor
from matplotlib.widgets import SpanSelector
import Tkinter as Tk
import tkMessageBox
import operator

# imports for img plot types
# try:
#    from PIL import Image
# except ImportError, exc:
#    raise SystemExit("PIL must be installed to run this example")
#
# import matplotlib.cbook as cbook
from matplotlib.cbook import get_sample_data
from pyCPEDScommonFunctions import cpedsPythCommon


















##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
#
# DEFINE GLOBAL VARIABLES
#
##########################################################################################
##########################################################################################
##########################################################################################

if type(option.width) != type(list()):    option.width = list([1])
if type(option.label) != type(list()):    option.label = list()
if option.mkLabels:
    for i in arange(len(args)):
        option.label.append(args[i])
plotLegendLabels = list()  # this is used to create the correct sub-list of labels for datasets in case some of the data are filtered out 

if type(option.textxy) != type(list()):    option.textxy = list()
if type(option.textfxfy) != type(list()):    option.textfxfy = list()
if type(option.annotate) != type(list()):    option.annotate = list()
if type(option.removeXtickLabels) != type(list()):    option.removeXtickLabels = list([0])
if type(option.removeYtickLabels) != type(list()):    option.removeYtickLabels = list([0])

if type(option.x2) != type(list()):    option.x2 = list()
if type(option.y2) != type(list()):    option.y2 = list()
if type(option.colx) != type(list()):    option.colx = list([0])
if type(option.coly) != type(list()):    option.coly = list([1])
if type(option.colz) != type(list()):    option.colz = list([2])
if type(option.badY) != type(list()):    option.badY = list(['None'])
if type(option.colColor) != type(list()):    option.colColor = list([-1])
if type(option.zorder) != type(list()):    option.zorder = list()
if type(option.xlabel) != type(list()):   option.xlabel = list(['x'])
if type(option.x2label) != type(list()):   option.x2label = list([None])
if type(option.ylabel) != type(list()):   option.ylabel = list(['y'])
if type(option.zlabel) != type(list()):   option.zlabel = list(['z'])
if type(option.CM) != type(list()):   option.CM = list(['hot'])
option.fontSize = cpedsPythCommon.getFloatList(option.fontSize)
option.fontSizeCM = cpedsPythCommon.getFloatList(option.fontSizeCM)
option.fontSizeLegend = cpedsPythCommon.getFloatList(option.fontSizeLegend)
option.fontSizeLabels = cpedsPythCommon.getFloatList(option.fontSizeLabels)

if option.bad_val=="None":
    option.bad_val=None
else:
    option.bad_val=float(option.bad_val)

# if option.set_below=="None":
#     option.set_below=None
# else:
#     set_below=option.set_below.split(',')
#     option.set_below=[float(set_below[0]),set_below[1]]


# option.fontSizeXtickLabels=cpedsPythCommon.getFloatList(option.fontSizeXtickLabels)
# option.fontSizeYtickLabels=cpedsPythCommon.getFloatList(option.fontSizeYtickLabels)
option.logYaxes = option.logYaxes.split(',')

print 'levels: ', option.levels
if type(option.levels) != type(list()):
    option.levels = list([50])
else:
    option.levels = [ (cpedsPythCommon.getFloatList(levels)) for levels in option.levels ]
print 'levels: ', option.levels
# for i in range(len(option.levels)):
#    print len(option.levels[i])
#    if len(option.levels[i])==1:
#        option.levels[i]=option.levels[i][0]
# print option.levels
# sys.exit()


colxNum = len(option.colx)
colyNum = len(option.coly)
colsxy = max(colxNum, colyNum)
inFilesNum = len(args)
if inFilesNum == 1:
	inFile = args[0]
	print colsxy
	print len(args)
	while colsxy > len(args):
		args.append(inFile)
		print args
if type(option.xmin) != type(list()):    option.xmin = list([None])
if type(option.xmax) != type(list()):    option.xmax = list([None])
if type(option.ymin) != type(list()):    option.ymin = list([None])
if type(option.ymax) != type(list()):    option.ymax = list([None])
if type(option.vmin) != type(list()):    option.vmin = list([None])
if type(option.vmax) != type(list()):    option.vmax = list([None])
if type(option.datexmin) != type(list()):    option.datexmin = list([-1])
if type(option.datexmax) != type(list()):    option.datexmax = list([-1])
if type(option.xticks) != type(list()):    option.xticks = list([0])
if type(option.xticksMinor) != type(list()):    option.xticksMinor = list([None])
if type(option.yticksMinor) != type(list()):    option.yticksMinor = list([None])
if type(option.yticks) != type(list()):    option.yticks = list([0])
if type(option.axesSwitch) != type(list()):    option.axesSwitch = list([1])
        
if type(option.colyerr) != type(list()):    option.colyerr = list([])
if type(option.colyerr2) != type(list()):    option.colyerr2 = list([])
if type(option.colxerr) != type(list()):    option.colxerr = list([])
if type(option.colxerr2) != type(list()):    option.colxerr2 = list([])
if type(option.hdu) != type(list()):    option.hdu = list([1])
if option.chColorsLTypes:
    if type(option.ls) != type(list()):    option.ls = list(['-', '--', '-.', '..'])
    if type(option.lc) != type(list()):    option.lc = list(['k', 'b', 'g', 'r'])
    if type(option.fc) != type(list()):    option.fc = list(['k', 'b', 'g', 'r'])
    if type(option.ec) != type(list()):    option.ec = list(['k', 'b', 'g', 'r'])
else:
    if type(option.ls) != type(list()):    option.ls = list(['-', '-', '-', '-', '-', '-', '-', '--', '--', '--', '--', '--', '--', '--', '-.', '-.', '-.', '-.', '-.', '-.', '-.'])
    if type(option.lc) != type(list()):    option.lc = list(['k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k'])
    if type(option.fc) != type(list()):    option.fc = list(['k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k'])
    if type(option.ec) != type(list()):    option.ec = list(['k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k', 'b', 'g', 'r', 'c', 'm', 'y', 'k'])

if type(option.lsdash) == type(list()):
    dstyleIdx = 0
    lsdashes = option.lsdash
    for lsd in lsdashes:
        lsdash = ast.literal_eval(lsd)
        option.lsdash[dstyleIdx] = lsdash
        dstyleIdx += 1
# print option.lsdash
# sys.exit(0)




# if type(option.fhs)!=type(list()):    option.fhs=list(['/','\\','|','-','+','x','o','O','.','*'])
if type(option.fhs) != type(list()):    option.fhs = list([None])
        
if type(option.plotType) != type(list()):    option.plotType = list(['fn'])

if option.plotType[0] == 'mayaVolCut':
    if type(option.palette) != type(list()):    option.palette = list(['spectral'])
else:
    if type(option.palette) != type(list()):    option.palette = list(['jet'])
    
if len(option.plotType) == 1 and option.plotType[0] == 'scat':
    if type(option.pt) != type(list()):    option.pt = list(['o'])
else:
    if type(option.pt) != type(list()):    option.pt = list([None])
if type(option.ps) != type(list()):    option.ps = list([5])
if type(option.markerEdgeWidth) != type(list()):    option.markerEdgeWidth = list([0])

if type(option.pc) != type(list()):    option.pc = list(['k', 'b', 'g', 'r', 'c', 'm', 'y'])
if type(option.title) != type(list()):    option.title = list([''])
if type(option.fileFormat) != type(list()):    option.fileFormat = list(['txt'])
if type(option.rows) != type(list()):    option.rows = list([0])
if type(option.every) != type(list()):    option.every = list([1])
if type(option.bin) != type(list()):    option.bin = list([1])
if type(option.where) != type(list()):    option.where = list()
# if type(option.val)!=type(list()):    option.val=list([])
# if type(option.oper)!=type(list()):    option.oper=list([])

if type(option.plotCircle) != type(list()):    option.plotCircle = list([])
# if type(option.circ)!=type(list()):    option.circ=list([])

option.oalpha = cpedsPythCommon.getFloatList(option.oalpha)

if option.Mcontours != "":
    mayaContours = cpedsPythCommon.getFloatList(option.Mcontours)

MINIMAL_VALUE_FOR_LOGPLOT = 1E-100
# rcParams["text.usetex"] = False
# import matplotlib.font_manager as fm
# fp1=fm.FontProperties(fname="/usr/share/fonts/liberation/LiberationSans-Regular.ttf") 

def toFloat(v):
    if v == '':
        return 0 
    return float(v)

def toInt(v):
    if v == '':
        return 0 
    return int(v)

Yoperations = option.operY.split(',')
YoperationValues = [ toFloat(v) for v in option.valueY.split(',')]
Xoperations = option.operX.split(',')
XoperationValues = [ toFloat(v) for v in option.valueX.split(',')]


globalWhichGridLines = 'major'
gridNbins = cpedsPythCommon.getIntList(option.gridNbins)
# print gridNbins
# sys.exit()

if type(option.Hline) != type(list()):    
    if option.Hline != None:
#        print option.Hline
#        sys.exit()
        option.Hline = list([cpedsPythCommon.getFloatList(option.Hline)])
    else:
        option.Hline = list()
else:
    for idx in range(len(option.Hline)):
        if option.Hline[idx] != '':
            option.Hline[idx] = cpedsPythCommon.getFloatList(option.Hline[idx])
        else:
            option.Hline[idx] = list()
    
if type(option.Vline) != type(list()):    
    if option.Vline != None:
        option.Vline = list([cpedsPythCommon.getFloatList(option.Vline)])
    else:
        option.Vline = list()
else:
    for idx in range(len(option.Vline)):
        if option.Vline[idx] != '':
            option.Vline[idx] = cpedsPythCommon.getFloatList(option.Vline[idx])
        else:
            option.Vline[idx] = list()


if type(option.Vspan) != type(list()):    option.Vspan = list()
else:
    for idx in range(len(option.Vspan)):
        if option.Vspan[idx] != '':
            option.Vspan[idx] = cpedsPythCommon.getFloatList(option.Vspan[idx])
        else:
            option.Vspan[idx] = list()
    
if type(option.Hspan) != type(list()):    option.Hspan = list()
else:
    for idx in range(len(option.Hspan)):
        if option.Hspan[idx] != '':
            option.Hspan[idx] = cpedsPythCommon.getFloatList(option.Hspan[idx])
        else:
            option.Hspan[idx] = list()

# option.olc=option.olc.split(',')
# option.ols=option.ols.split(',')

if type(option.ols) != type(list()):    option.ols = list([list('-')])
else:
    for idx in range(len(option.ols)):
        if option.ols[idx] != '':
            option.ols[idx] = option.ols[idx].split(',')
        else:
            option.ols[idx] = list()

if type(option.olc) != type(list()):    option.olc = list([list('k')])
else:
    for idx in range(len(option.olc)):
        if option.olc[idx] != '':
            option.olc[idx] = option.olc[idx].split(',')
        else:
            option.olc[idx] = list()

if type(option.ofc) != type(list()):    option.ofc = list([list('r')])
else:
    for idx in range(len(option.ofc)):
        if option.ofc[idx] != '':
            option.ofc[idx] = option.ofc[idx].split(',')
        else:
            option.ofc[idx] = list()

if type(option.olw) != type(list()):    option.olw = list([list([1])])
else:
    for idx in range(len(option.olw)):
        if option.olw[idx] != '':
#            option.olw[idx]=option.olw[idx].split(',')
            option.olw[idx] = cpedsPythCommon.getIntList(option.olw[idx])
        else:
            option.olw[idx] = list()

print 'olw: ', option.olw
# sys.exit()
print "Yoperations: "
print Yoperations
print "Yoperation values"
print YoperationValues
# sys.exit()

if len(option.meridians.split(',')) == 1:
    if option.meridians.isdigit():
        meridians = arange(0, 360, float(option.meridians))
    else:
        meridians = array([])
else:
    meridians = array([ toFloat(v) for v in option.meridians.split(',') ])
    
meridianLabels = None
if option.meridianLabels != '':
    meridianLabels = option.meridianLabels.split(',')
    

if len(option.parallels.split(',')) == 1:
    if option.parallels.isdigit():
        if option.proj == 'moll':
            parallels = arange(-90, 90, float(option.parallels))
        else:
            if option.proj == 'lcc':
                parallels = arange(-90, 90, float(option.parallels))
            else:
                if option.proj == 'ortho':
                    parallels = arange(0, 90, float(option.parallels))
                else:
                    parallels = array([])
    else:
        parallels = array([])
else:
    parallels = array([ toFloat(v) for v in option.parallels.split(',') ])

print meridians    
for pt in option.plotType:
    if pt == "map" or pt == "sphere":
        option.xlabel = list(['None'])
        option.ylabel = list(['None'])


mapBbox = None
mapPaths = None
projectedMapData = None
projectMap = None


if len(option.IMGextent.split(',')) == 1:
    option.IMGextent = array([100, 100])
else:
    option.IMGextent = array([ toFloat(v) for v in option.IMGextent.split(',') ])


if option.circ != '':
    circlePatchDef = [ toFloat(v) for v in option.circ.split(',')]
else:
    circlePatchDef = list()

if type(option.rect) != type(list()):    
    option.rect = list()
else:
    for idx in range(len(option.rect)):
        if option.rect[idx] != '':
            option.rect[idx] = option.rect[idx].split(';')  # for each sub-plot extract list of rectangles to plot
            for idx2 in range(len(option.rect[idx])): 
                option.rect[idx][idx2] = cpedsPythCommon.getFloatList(option.rect[idx][idx2])  # extract rectangle parameters for each rectangle
        else:
            option.rect[idx] = list()

# rectanglePatchDef=option.rect

legendLocation = list()
for lloc in option.legendLoc.split(','):
    if lloc == "ur":
        legendLocation.append('upper right')
    if lloc == "ul":
        legendLocation.append('upper left')
    if lloc == "br":
        legendLocation.append('lower right')
    if lloc == "bl":
        legendLocation.append('lower left')
    
    if lloc == "None":
        legendLocation.append('None')

scatterGlobalData = []
globalScatterColorbar = None
globalCurrentDatasetIdx = 0
globalCurrentDataset = None
globalCurrentDatasetX = None
globalCurrentDatasetY = None
globalStdInData = ''
globalPlotTypeIdx = -1
globalPlotDs = None
globalSubPlotGridSpec = None
globalAxes = None
globalLassoManager = None

fieldPixNuml = [ toInt(v) for v in option.fieldPixNum.split(',')]
if option.PRcamZ == 0:
    option.PRcamZ = -1.1 * fieldPixNuml[2 % len(fieldPixNuml)]


hdf5dset = option.hdf5dset.split(',')








###########################################################################################
###########################################################################################
###########################################################################################
# INTERACTIVE OPERATION MODE FUNCTIONS
###########################################################################################
###########################################################################################
###########################################################################################


maskLines = []
maskRanges = []
selectedLines = []
selectedLinesData = []
selectedLinesMode = False
cursorAsCross = False
cursor = None
spannerOn = False
span = None
spannerSelectLinesOn = False
spanSelectLines = None
binFromWidget = None
binWidthWidget = None
binWidthMultiplierWidget = None
binningParams = { 'st':1.0, 'bw': 10, 'gm':1.05 }
lastBinnedSpectra = []
saveWidgetMain = None
saveWidget = None
plotAxes = None
plotAxesList = []
shareX=None
plotFig = None
continuousPottingFirstTime = True
saveToFile = None
saveToFileWidget = None
periodicMaskSamplesNumer = 0
globalDeleteDataMode = False
globalLassoDataMode = False
globalMaskedPointsMap = np.array([])
globalMaskedPointsMarkers = []
globalNfiles = 0

def printInteractiveHelp():
    print "The following keys are active in the mask making mode:"
    print "space - sets and stores a new [hv]line and prints its info on terminal"
    print "d - deletes the last added 'selected line'"
    print "p - prints the information on the stored lines"
    print "w - saves the defined lines information to file called selectedMaskLines.txt and selectedMaskRanges.txt it also stores the binning information"
    print "e - do the spike extraction: save the list of ranges and levels to be removed, run the extract_spikes program and plot results in a new window"
    print "c - toggle cursor as cross"
    print "v - toggle spanner"
    print "f1 - prints this help"
    print "f2 - mark and add 'selected like'"
    print "f3 - toggle selectedLines edit mode and maskRanges edit mode"
    print "f4 - toggle spanner for electing lines"
    print "D - mark points for creating mask using picker "
    print "Q - mark points for creating mask using lasso "
    print "M - write data mask to file"
    print "m - load data mask from file"
    print "P - set period for the periodic operations"
    print ""
    print "1-mouse button - view info on selected point"
    print "2-mouse button - remove closest selectedLine/maskRange"


#--------------------------------------------------------------------------------------------------------
def loadMask(m):
    print 'Loading mask file: ', m
    tmp = np.loadtxt(m)
    tmp = np.asarray(tmp)
    tmp = tmp.reshape([-1, 2])
    return list(tmp)


def plotButtons(ax):
    bprev = Button(ax, 'Extract Spikes')
    bprev.on_clicked(extractSpikes)


def printCurrentMaskRanges():
    print "Current mask ranges are:"
    print maskRanges
def printCurrentMaskLines():
    print 'Current mask lines are:'
    print maskLines

def printSelectedLines():
    print 'Selected lines are:'
    print selectedLines


def plotSelectedLines():
    t = array(selectedLines)
    if len(t) > 0:
        plot(t[:, 0], t[:, 1], 'yo', ms=6, label='selected lines', picker=6)
        draw()

# def calculatePlotID(i,figRows,figCols):
#    global globalPlotDs
#    global globalNfiles
#    offset=0
# #    plotidx=0
#    subPlotIdx0=0
#    k=0
#    for ii in range(globalNfiles):
# #        for jj in range(figCols):
#        if i==k+offset:
#            return [figRows,figCols, subPlotIdx0 + 1],k
#        else:
#            k=k+1
#            if k==globalPlotDs[subPlotIdx0]+offset:
#                offset=offset+globalPlotDs[subPlotIdx0]
#                subPlotIdx0=subPlotIdx0+1
#                k=0
#    return [figRows,figCols, figRows*figCols],k

# def calculatePlotID(dsIdx,figRows,figCols):
#    global globalPlotDs
#    global globalNfiles
# #    offset=0
# #    plotidx=0
#    subPlotIdx0=0
#    subplotDsIdx=0
#    for dummyIdx in range(globalNfiles):
# #        for jj in range(figCols):
#        if dsIdx==subplotDsIdx: #+offset:
#            return [figRows,figCols, subPlotIdx0 + 1],subplotDsIdx
#        else:
#            subplotDsIdx=subplotDsIdx+1
#            if subplotDsIdx==globalPlotDs[subPlotIdx0]: #+offset:
# #                offset=offset+globalPlotDs[subPlotIdx0]
#                subPlotIdx0=subPlotIdx0+1
#                subplotDsIdx=0
#    return [figRows,figCols, figRows*figCols],subplotDsIdx

def calculatePlotID(dsIdx, figRows, figCols):
    global globalPlotDs
    global globalNfiles
#    offset=0
#    plotidx=0
    subPlotIdx0 = 0
    subplotDsIdx = 0
    for dummyIdx in range(globalNfiles):
#        for jj in range(figCols):
        if dsIdx == 0:  # +offset:
            return [figRows, figCols, subPlotIdx0 + 1], subplotDsIdx
        else:
            dsIdx = dsIdx - 1
            subplotDsIdx = subplotDsIdx + 1
            if subplotDsIdx == globalPlotDs[subPlotIdx0]:  # +offset:
#                offset=offset+globalPlotDs[subPlotIdx0]
                subPlotIdx0 = subPlotIdx0 + 1
                subplotDsIdx = 0
    return [figRows, figCols, figRows * figCols], subplotDsIdx


def isNewPlot(i, figRows, figCols):
    print 'DEBUG: ', i, figRows, figCols
    t, k = calculatePlotID(i, figRows, figCols)
    if k == 0:
        return True
    return False
    

def makeBinningInfoWidget(event):

    root = Tk.Tk()
    
    root.title("Power spectra binning parameters:")
    root["padx"] = 40
    root["pady"] = 20       

    # Create a text frame to hold the text Label and the Entry widget
    textFrame = Tk.Frame(root)
    
    # Create a Label in textFrame
    binFromLabel = Tk.Label(textFrame)
    binFromLabel["text"] = "Bin from: "
    binFromLabel.pack(side=Tk.LEFT)
    # Create an Entry Widget in textFrame
    global binFromWidget
    binFromWidget = Tk.Entry(textFrame)
    binFromWidget["width"] = 50
    binFromWidget.insert(0, str(binningParams['st']))
    binFromWidget.pack(side=Tk.LEFT)

    # Create a Label in textFrame
    binWidthLabel = Tk.Label(textFrame)
    binWidthLabel["text"] = "First bin width: "
    binWidthLabel.pack(side=Tk.LEFT)
    # Create an Entry Widget in textFrame
    global binWidthWidget
    binWidthWidget = Tk.Entry(textFrame)
    binWidthWidget["width"] = 50
    binWidthWidget.insert(0, str(binningParams['bw']))
    binWidthWidget.pack(side=Tk.LEFT)

    # Create a Label in textFrame
    binWidthMultiplierLabel = Tk.Label(textFrame)
    binWidthMultiplierLabel["text"] = "Bin width multiplier: "
    binWidthMultiplierLabel.pack(side=Tk.LEFT)
    # Create an Entry Widget in textFrame
    global binWidthMultiplierWidget
    binWidthMultiplierWidget = Tk.Entry(textFrame)
    binWidthMultiplierWidget["width"] = 50
    binWidthMultiplierWidget.insert(0, str(binningParams['gm']))
    binWidthMultiplierWidget.pack(side=Tk.LEFT)
    textFrame.pack()

    button = Tk.Button(root, text="Submit", command=setBinningParameters)
    button.pack()
    
    root.mainloop()

def deleteAllMaskRanges(event):
    global maskRanges
    axes(plotAxes)
    while len(maskRanges) > 0:
        del(gca().patches[-1]) 
        del maskRanges[-1]
    draw()
    printCurrentMaskRanges()

#    root = Tk.Tk()
#    
#    root.title("Periodic mask parameters:")
#    root["padx"] = 40
#    root["pady"] = 20       
#
#    # Create a text frame to hold the text Label and the Entry widget
#    textFrame = Tk.Frame(root)
#    
#    #Create a Label in textFrame
#    binFromLabel = Tk.Label(textFrame)
#    binFromLabel["text"] = "period [samples]: "
#    binFromLabel.pack(side=Tk.LEFT)
#    # Create an Entry Widget in textFrame
#    binFromWidget = Tk.Entry(textFrame)
#    binFromWidget["width"] = 50
#    binFromWidget.insert(0,str(10))
#    binFromWidget.pack(side=Tk.LEFT)
#
#
#    button = Tk.Button(root, text="Submit", command=makeMask(10))
#    button.pack()
#    
#    root.mainloop()
#    p,err=cpedsPythCommon.readUserValue('samples number', 'float', '10')
#    if err==0:
    makePeriodicMask(p)
    
def makePeriodicMask(samplesNum):
    print 'samples: ', samplesNum
    global maskRanges
    global periodicMaskSamplesNumer
    global globalCurrentDatasetX
    print 'len(maskRanges): ', len(maskRanges)
    
    if periodicMaskSamplesNumer == 1.0 and len(maskRanges) >= 2:
        xmin = maskRanges[-1][0]
        xmax = maskRanges[-1][1]
        delta = maskRanges[-1][0] - maskRanges[-2][0]
        xMax, err = cpedsPythCommon.readUserValue('enter xmax value to proceed with the mask', 'float', 0)
        x = xmin
        i = 1
        print 'xMax: ', xMax
        print 'x: ', x
        while x < xMax:
            r = [xmin + i * delta, xmax + i * delta]
            print 'adding range: ', r
            maskRanges.append(r)
            axvspan(r[0], r[1], facecolor='b', edgecolor=None, alpha=0.3, zorder=-100)
            x = xmax + i * delta
            i = i + 1
        
        periodicMaskSamplesNumer = 0
        
        plotFig.canvas.draw()
        printCurrentMaskRanges()
    

def saveLastBinnedSpectraWidget(event):
    global saveWidget
    global saveWidgetMain
    saveWidgetMain = Tk.Tk()
    
    saveWidgetMain.title("Save binned power spectra:")
    saveWidgetMain["padx"] = 40
    saveWidgetMain["pady"] = 20       

    # Create a text frame to hold the text Label and the Entry widget
    textFrame = Tk.Frame(saveWidgetMain)
    
    # Create a Label in textFrame
    saveLabel = Tk.Label(textFrame)
    saveLabel["text"] = "File name: "
    saveLabel.pack(side=Tk.LEFT)
    # Create an Entry Widget in textFrame
    saveWidget = Tk.Entry(textFrame)
    saveWidget["width"] = 50
    fname = dsName[0] + "-" + DA[0] + ".bin_st" + str(binningParams['st']) + "_bw" + str(binningParams['bw']) + "_gm" + str(binningParams['gm'])
    saveWidget.insert(0, fname)
    saveWidget.pack(side=Tk.LEFT)
    
    textFrame.pack()
    button = Tk.Button(saveWidgetMain, text="Save", command=saveLastBinnedSpectra)
#    button = Tk.Button(saveWidgetMain, text="Save", command=lambda arg=saveWidget.get(): saveLastBinnedSpectra(arg)) 

    button.pack()
    
    saveWidgetMain.mainloop()


def saveLastBinnedSpectra():
# def saveLastBinnedSpectra(fname):
#    fname=dsName[0]+"-"+DA[0]+".bin"
    global saveWidget
    global saveWidgetMain
    fname = saveWidget.get().strip()
    np.savetxt(fname, lastBinnedSpectra)
    print "saved spectra to file: " + fname
    saveWidgetMain.destroy()
    saveWidgetMain = None



def deleteAllMaskLines(event):
        axes(plotAxes)
        global maskLines
        print "mask lines" + str(len(maskLines))
        print "gca lines" + str(len(gca().lines))
        st = len(gca().lines) - 2 * len(maskLines)
        print "st" + str(st)
        del gca().lines[st:len(gca().lines)]
        maskLines = []
        draw()




def deleteMaskLine(event):
        axes(plotAxes)
        del(gca().lines[-1]) 
        del(gca().lines[-1]) 
#        del(linesData[-1])
        del(maskLines[-1])
        draw()

def deleteMaskRange(event):
        axes(plotAxes)
        maskRange = maskRanges[-1]
        del(gca().patches[-1]) 
#        del(linesData[-1])
        del(maskRanges[-1])
        draw()
        printCurrentMaskRanges()

        # update masked points map         
        maskMaskedPointsMapRange(maskRange[0], maskRange[1], 1)


def deleteLastSelectedLine(n=1):
        sl = axes(plotAxes).lines
        i = len(sl) - 1
        d = n
        print "selected lines count" + str(len(sl))
        while i >= 0:
            print i
            print sl[i].get_label()
            if sl[i].get_label() == 'selected lines':
#                sl.pop(i)
                del sl[i]
#                axes(ax1).lines.pop(i)
                d = d - 1
                if d == 0:
                    i = 0
#                i=i+1
            i = i - 1
        del(selectedLines[-1])
        draw()
        printSelectedLines()

def clearSelectedLines():
        
    sl = axes(plotAxes).lines
    i = len(sl) - 1
    print "selected lines count " + str(len(sl))
    while i >= 0:
        print i
        print sl[i].get_label()
        if sl[i].get_label() == 'selected lines':
            sl.pop(i)
        i = i - 1
    draw()

# n - line index to remove
def clearSelectedLine(n):
        
    sl = axes(plotAxes).lines
    i = 0
    k = 0
    print "selected lines count " + str(len(sl))
    while i < len(sl):
        print "i: " + str(i)
        print "k: " + str(k)
        print sl[i].get_label()
        if sl[i].get_label() == 'selected lines':
            if k == n:
                del sl[i]
                i = len(sl)
                print "removed"
            else:
                print "not removed"
            k = k + 1
        i = i + 1
        
    draw()
    
def deleteSelectedLine(event):
    sl = axes(plotAxes).lines
#        i=len(sl)-1
    print "selected lines count" + str(len(sl))
    x = event.xdata
    y = event.ydata
    slt = array(selectedLines)
    dx = (slt[:, 0] - x); dx2 = dx * dx
    dy = (slt[:, 1] - y); dy2 = dy * dy
    d = sqrt(dx2 + dy2)
    row = d.argmin()
#        print x
#        print y
#        print d
#    print "will remove line: "+ str(row)
#        j=0
#        while i>=0:
#            print i
#            print sl[i].get_label()
#            if sl[i].get_label()=='selected lines':
#                if j==row:
#                    del sl[i]
#                j=j+1
#            i=i-1

#    if (sqrt(event.x-*event.x+event.y*event.y)d[row]<10:
#    print "row:"+str(row)
    del(selectedLines[row])
#    print "row:"+str(row)
#    clearSelectedLine(row)
    clearSelectedLines()
    plotSelectedLines()
    draw()
    printSelectedLines()


def loadDataPointsMask():
    global globalCurrentDatasetIdx
    global globalMaskedPointsMarkers
    global globalMaskedPointsMap
    global globalCurrentDatasetX
    global globalCurrentDatasetY
    
    fname = args[globalCurrentDatasetIdx] + '.mask'
    print 'loading data points mask from file: ', fname
    globalMaskedPointsMap = np.loadtxt(fname, dtype=int)
    print globalMaskedPointsMap
    
    globalMaskedPointsMarkers = []

    for idx in range(len(globalMaskedPointsMap)):
        if globalMaskedPointsMap[idx] == 0:
            x = globalCurrentDatasetX[idx]
            y = globalCurrentDatasetY[idx]
            axes(plotAxes)
            PT, = plot(x, y, 'x', ms=option.ps[globalCurrentDatasetIdx % len(option.ps)] + 1, picker=True)
            globalMaskedPointsMarkers.append([PT, idx])

    
    draw()


def maskLastPickedDataPoint(event):
    global globalMaskedPointsMap
    global globalMaskedPointsMarkers
    global globalCurrentDataset
    global globalCurrentDatasetIdx
    
    remove = []
#    for i in range(len(globalMaskedPointsMarkers)):
#        if globalMaskedPointsMarkers[i]==event.artist:
#            remove.append(artist)
    print 'event.ind: ', event.ind
#    remove = [artist for artist in globalMaskedPointsMarkers if event.artist.contains(artist)]
    print globalMaskedPointsMarkers

    if all(globalMaskedPointsMap[event.ind]) == 0:
        globalMaskedPointsMap[event.ind] = 1
        remove = []
        for i in range(len(globalMaskedPointsMarkers)):
            if globalMaskedPointsMarkers[i][1] in event.ind:
#                globalMaskedPointsMarkers[i][0].remove()
                remove.append(i)
        print 'remove'
        print remove
        if len(remove) > 0:
            remove = sorted(remove, reverse=True)
        ax = gca()
        print remove
        print ax.lines
        for r in remove:
            del(ax.lines[r])
            del(globalMaskedPointsMarkers[r])
    else:
        globalMaskedPointsMap[event.ind] = 0
        for idx in event.ind:
            x = globalCurrentDatasetX[idx]
            y = globalCurrentDatasetY[idx]
            print 'this point will be masked'
            axes(plotAxes)
            PT, = plot(x, y, 'x', ms=option.ps[globalCurrentDatasetIdx % len(option.ps)] + 1, color='r', picker=True)
            globalMaskedPointsMarkers.append([PT, idx])
            print len(gca().lines)
        
#    event.artist.remove()
#    print 'event.ind: ',event.ind
    print globalMaskedPointsMap
     
#    if not remove:
#        x=globalCurrentDatasetX[event.ind]
#        y=globalCurrentDatasetY[event.ind]
# ##        x, y = ax.transData.inverted().transform_point([event.x, event.y])
#        print 'this point will be masked'
#        axes(plotAxes)
#        PT,=plot(x,y,'x',ms=option.ps[globalCurrentDatasetIdx % len(option.ps)]+1, picker=True)
# #        globalMaskedPointsMarkers.append(len(globalMaskedPointsMarkers))
#        globalMaskedPointsMarkers.append([PT,event.idx])
#    else:
#        for artist in remove:
#            artist.remove()

#        i=0
#        while i<len(gca().collections):
#            if len(gca().collections==
#            del(gca().patches[-1]) 
#            del(maskRanges[-1])
#            draw()
#            printCurrentMaskRanges()
#        ax=gca();


    draw()


def maskLastLassoDataPoints(ind):
    global globalMaskedPointsMap
    global globalMaskedPointsMarkers
    global globalCurrentDataset
    global globalCurrentDatasetIdx
    
    remove = []

    if all(globalMaskedPointsMap[ind]) == 0:
        globalMaskedPointsMap[ind] = 1
        remove = []
        for i in range(len(globalMaskedPointsMarkers)):
            if globalMaskedPointsMarkers[i][1] in ind:
                remove.append(i)
        print 'remove'
        print remove
        if len(remove) > 0:
            remove = sorted(remove, reverse=True)
        ax = gca()
        print remove
        print ax.lines
        for r in remove:
            del(ax.lines[r])
            del(globalMaskedPointsMarkers[r])
    else:
        globalMaskedPointsMap[ind] = 0
        for idx in ind:
            x = globalCurrentDatasetX[idx]
            y = globalCurrentDatasetY[idx]
            print 'this point will be masked'
            axes(plotAxes)
            PT, = plot(x, y, 'x', ms=option.ps[globalCurrentDatasetIdx % len(option.ps)] + 1, color='r', picker=True)
            globalMaskedPointsMarkers.append([PT, idx])
            print len(gca().lines)
        
    print globalMaskedPointsMap
    draw()





def toggleCursor(event):
    global cursor
    global cursorAsCross
    if cursorAsCross:
        cursorAsCross = False
        cursor = Cursor(plotAxes, useblit=True, color='green', linewidth=2)
    else:
        del cursor
        cursorAsCross = True
        cursor = None

def toggleSpanner(event):
    global span
    global spannerOn
#    global cursorAsCross
    if spannerOn:
        spannerOn = False
        if span != None:
            span = None
    else:
        spannerOn = True
        span = SpanSelector(plotAxes, onselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='g'))

def toggleSpannerSelectLines(event):
    global spanSelectLines
    global spannerSelectLinesOn
#    global cursorAsCross
    if spannerSelectLinesOn:
        spannerSelectLinesOn = False
        if spanSelectLines != None:
            spanSelectLines = None
    else:
        spannerSelectLinesOn = True
        spanSelectLines = SpanSelector(plotAxes, onselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='r'))

def extractSpikes(event):
        fn = '.selectedMaskLines.tmp';
        mask = array(maskLines)
        np.savetxt(fn, mask)
        fn = '.selectedMaskRanges.tmp';
        mask = array(maskRanges)
        np.savetxt(fn, mask)
        global binningParams
    
        s = dsName[i] + '.power/' + dsName[i] + '.' + DA[j] + '.p' + binnedSuffix
        rrastr = ''
        if len(maskLines) > 0:
            rrastr = ' --removeRangeAbove .selectedMaskLines.tmp '
        rrstr = ''
        if len(maskRanges) > 0:
            rrstr = ' --removeRange .selectedMaskRanges.tmp '

        cmd = 'extract_spikes ' + s + '  -o ' + s + '.nospikes.tmp' + ' --maskMode ' + rrstr + rrastr
        os.system(cmd)
        cmd = 'bin_function ' + s + '.nospikes.tmp --fromAsArg -f ' + str(binningParams['st']) + ' -d ' + str(binningParams['bw']) + ' --gm ' + str(binningParams['gm']) + '  --geo -o ' + s + '.nospikes.tmp.bin'
        os.system(cmd)

        tmpfig = figure(figsize=(15, 10), dpi=opt['DPIgui'])

        tmpax = tmpfig.add_subplot(1, 1, 1)
        s = dsName[i] + '.power/' + dsName[i] + '.' + DA[j] + '.p' + binnedSuffix + '.nospikes.tmp'
        spectra = loadPower(s)
        tmpax.plot(spectra[:, 0], spectra[:, 1], dsStyle[0], label=dsName[0] + "(" + DA[j] + ") no spikes", alpha=.5, picker=False)
        spectra2 = loadPower(s + '.bin')
        tmpax.plot(spectra2[:, 0], spectra2[:, 1], dsStyle[1] + 'v', label=dsName[0] + "(" + DA[j] + ") no spikes bin", alpha=.5, picker=False)
        print "Binned spectra consists of %li points" % len(spectra2[:, 1])
        print "Binned spectra consists of %li points" % len(spectra2[:, 1])

        setupAxes(spectra, 'loglinear', tmpax, False, False, False)
        global lastBinnedSpectra
        lastBinnedSpectra = spectra2
        printPowerSpectraInfo()
        show()


def simplifyMask(mask):
    print "* Removing redundant mask ranges"
    i = 0
    j = 0
    rr = 0
    while i < len(mask):
        j = 0
        while j < len(mask):
            if i != j and mask[i][0] <= mask[j][0] and mask[i][1] >= mask[j][1]:
                print "mask range: (%lf, %lf) includes (%lf, %lf). I will remove it." % (mask[i][0], mask[i][1], mask[j][0], mask[j][1])
                mask.pop(j)
                j = j - 1
                if i > j: i = i - 1
                rr = rr + 1
#            print "i: %li j: %li\n" % (i,j)
            j = j + 1
        i = i + 1
    print "* DONE. Removed %li redundant mask ranges" % rr

    return mask, rr


def on_key(event):
    global selectedLinesMode
    global globalCurrentDataset
    global globalCurrentDatasetIdx
    global periodicMaskSamplesNumer
# SELECTED LINES    
    if event.key == 'd':
#        deleteMaskLine(event)
        deleteLastSelectedLine(1)

    if event.key == 'f':
        fn = dsName[0] + "-" + DA[0] + '.selectedLines.txt'
        sl = array(selectedLines)
        np.savetxt(fn, sl)
        print "selected lines saved to file: " + fn

    if event.key == 'f2':
        xdata = event.xdata
        ydata = event.ydata
        print "adding new point data: %lf, %fE" % (xdata, ydata)
        selectedLines.append([xdata, ydata])
#        sl=axes(ax1).lines
#        sl.append

        printSelectedLines()
#        clearSelectedLines()
        plotSelectedLines()

    if event.key == 'f3':
        if selectedLinesMode == False:
            print "* SETTING SELECTED LINES MODE ON (2-MOUSE BUTTON WILL REMOVE SELECTED LIES)"
            selectedLinesMode = True
        else:
            print "* SETTING SELECTED LINES MODE OFF (2-MOUSE BUTTON WILL REMOVE SELECTED MASK RANGE)"
            selectedLinesMode = False


    if event.key == 'f4':
        toggleSpannerSelectLines(event)
        if (spannerSelectLinesOn):
            print "spanner SelectLines is ON"
        else:
            print "spanner SelectLines is OFF"



# SELECTED MASK RANGES AND LINES
    if event.key == 'i':
        printPowerSpectraInfo()

    if event.key == 'w':
        fn = 'lastPlotFunction.selectedMaskLines.txt'
#        fn='selectedMaskLines.txt';
        mask = array(maskLines)
        np.savetxt(fn, mask)
        print "mask lines saved to file: " + fn
        fn = 'lastPlotFunction.selectedMaskRanges.txt'
#        fn='selectedMaskRanges.txt';

        mask = sorted(maskRanges, key=operator.itemgetter(0))
        i = len(mask) - 1
        while i >= 0:
            if abs(mask[i][1] - mask[i][0]) < 1e-5:
                print "removing empty range: %lE %lE\n" % (mask[i][0], mask[i][1])
                mask.pop(i)
            i = i - 1
        mask, rr = simplifyMask(mask)
        np.savetxt(fn, array(mask))
        print "mask ranges saved to file: " + fn

#        fn=dsName[0]+"."+DA[0]+'.selectedBinning.txt'
#        selectedBinning=array([binningParams['st'],binningParams['bw'],binningParams['gm']])
#        np.savetxt(fn,selectedBinning)
#        print "selected binning parameters saved to file: "+fn
        
    if event.key == 'c':
        toggleCursor(event)

#    if event.key=='C':
    if event.key == 'v':
        toggleSpanner(event)
        if spannerOn:
            print "spanner is ON"
        else:
            print "spanner is OFF"

    if event.key == 'P':
        periodicMaskSamplesNumer = float(raw_input('Enter samples number: (1 to read from the last two spans)'))
#        periodicMaskSamplesNumer,err=cpedsPythCommon.readUserValue('Enter samples number', 'float', '10')
#        if err==0:
        makePeriodicMask(periodicMaskSamplesNumer)
            
    if event.key == ' ':
        print event.xdata, event.ydata
        line = axhline(y=event.ydata, xmax=1, color='r', lw=2)
        line = axvline(x=event.xdata, ymax=1, color='r', lw=2)
#        linesData.append([event.xdata, event.ydata])
        last = len(maskLines) - 1
        print last
        if (len(maskLines) == 0):
            maskLines.append([event.xdata, 138, event.ydata])
        else:
            maskLines.append([event.xdata, maskLines[last][0], event.ydata])
        draw()
    if event.key == 'p':
        printCurrentMaskLines()
        printCurrentMaskRanges()


    if event.key == 'e':
        extractSpikes()


    if event.key == 'D':
        global globalDeleteDataMode
        if globalDeleteDataMode:
            globalDeleteDataMode = False
            print 'globalDeleteDataMode: OFF'
        else:
            globalDeleteDataMode = True
            print 'globalDeleteDataMode: ON'
#        x = event.xdata
#        y = event.ydata
#        ind = event.ind
#        deleteSelectedDataPoint(event)
#        print 'deleted data point:', globalCurrentDataset[ind]
    if event.key == 'Q':
        global globalLassoDataMode
        if globalLassoDataMode:
            globalLassoDataMode = False
            print 'globalLassoDataMode: OFF'
        else:
            globalLassoDataMode = True
            print 'globalLassoDataMode: ON'
            global globalCurrentDatasetX, globalCurrentDatasetY, globalLassoManager
#             print 'globalCurrentDataset'
#             print globalCurrentDataset
#             print 'globalCurrentDataset stackXY'
#             print np.vstack([globalCurrentDatasetX,globalCurrentDatasetY]).T
            lassoData = [Datum(*xy) for xy in np.vstack([globalCurrentDatasetX, globalCurrentDatasetY]).T]
            ax = gca()
            globalLassoManager = LassoManager(ax, lassoData)
        
    if event.key == 'M':
        global globalMaskedPointsMap
        fname = args[globalCurrentDatasetIdx] + '.mask'
        print 'saving data mask to file: ', fname
        np.savetxt(fname, globalMaskedPointsMap.reshape((len(globalMaskedPointsMap), 1)), fmt='%i')
        fnameIn = args[globalCurrentDatasetIdx]
        fnameMasked = fnameIn + '.masked'
        fnameSelected = fnameIn + '.selected'
        if option.fits:
            print 'using read data: ', fnameIn
            data = globalCurrentDataset
            data = np.hstack([data, globalMaskedPointsMap.reshape((len(globalMaskedPointsMap), 1))])
            dataMasked = data[data[:, -1] == 1]
            dataMasked = dataMasked[:, 0:-1]
            np.savetxt(fnameMasked, dataMasked)
            dataSelected = data[data[:, 1] == 0]
            dataSelected = dataSlected[:, 0:-1]
            np.savetxt(fnameSelected, dataSelected)
        else:
            print 'reading data from file: ', fnameIn
            file = open(fnameIn, 'r')
            data = np.array(file.readlines())
            file.close()
#             data=globalCurrentDataset
            mask = globalMaskedPointsMap.reshape((len(globalMaskedPointsMap), 1))
            print mask
#             data=np.hstack([data,globalMaskedPointsMap])
#             data=np.hstack([
#                 data.reshape((len(globalMaskedPointsMap),1)),
#                 globalMaskedPointsMap.reshape((len(globalMaskedPointsMap),1))
#                 ])
            print 'saving masked data to file: ', fnameMasked
#            print data
            dataMasked = data[mask[:, -1] == 1]
            #        print dataMasked
            file = open(fnameMasked, "w")
            file.writelines(dataMasked)
            file.close()
            #        saveDataMaskToFile()

            dataSelected = data[mask[:, -1] == 0]
            file = open(fnameSelected, "w")
            file.writelines(dataSelected)
            file.close()

    if event.key == 'm':
        loadDataPointsMask()

# HELP    
    
    if event.key == 'f1':
        printInteractiveHelp()



# def deleteSelectedDataPoint(event):
#    x = event.xdata
#    y = event.ydata

#    ind = event.ind
#    thisline = event.artist
#    xdata = thisline.get_offsets()[:,0]
#    ydata = thisline.get_offsets()[:,1]
#    r=thisline.get_sizes()
#    print 'deleted data point:', globalCurrentDataset[ind]
    


def onselect(xmin, xmax):
#    indmin, indmax = np.searchsorted(x, (xmin, xmax))
#    indmax = min(len(x)-1, indmax)

#    thisx = x[indmin:indmax]
#    thisy = y[indmin:indmax]
#    line2.set_data(thisx, thisy)
#    ax2.set_xlim(thisx[0], thisx[-1])
#    ax2.set_ylim(thisy.min(), thisy.max())
    axes(plotAxes)
    if spannerOn:
        axvspan(xmin, xmax, facecolor='b', edgecolor=None, alpha=0.3)
        maskRanges.append([xmin, xmax])
    elif spannerSelectLinesOn:
        axvspan(xmin, xmax, facecolor='y', edgecolor=None, alpha=0.3)
        selectedLinesData.append([xmin, xmax])
        
    plotFig.canvas.draw()
    printCurrentMaskRanges()

    #
    # update data points mask
    #
    maskMaskedPointsMapRange(xmin, xmax)
    
def maskMaskedPointsMapRange(vmin, vmax, val=0):
    global globalCurrentDataset
    global globalCurrentDatasetX
    global globalMaskedPointsMap
#     print globalCurrentDataset
    
    idx = np.argwhere(np.logical_and(globalCurrentDatasetX >= vmin, globalCurrentDatasetX <= vmax))
    globalMaskedPointsMap[idx] = val
    
    print globalMaskedPointsMap

# def onclick(event):
#    print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
#        event.button, event.x, event.y, event.xdata, event.ydata)
    
# def onpick3(event):
#    ind = event.ind
#    print 'onpick2 line:', event.pickx, event.picky
#    print 'onpick3 scatter:', ind, np.take(x, ind), np.take(y, ind)

def deletePointedMaskRange(event):
#    axes(ax1)
#    del(gca().patches[-1]) 
# #        del(linesData[-1])
#    del(maskRanges[-1])
#    draw()
#    print "removing mask range: "
#    printCurrentMaskRanges()

    axes(plotAxes)
    x = event.xdata
    print "removing at x=%f " % x
    i = len(maskRanges) - 1
    while i >= 0:

#    for maskRange in maskRanges:
        maskRange = maskRanges[i]
        if x > maskRange[0] and x < maskRange[1]:
            maskRanges.pop(i)
            del(gca().patches[i]) 
            i = 0
        i = i - 1
    plotFig.canvas.draw()
    printCurrentMaskRanges()
    maskMaskedPointsMapRange(maskRange[0], maskRange[1], 1)

#    if isinstance(event.artist, Line2D):
#        thisline = event.artist
#        ydata = thisline.get_ydata()
#        ind = event.ind
#        print 'onpick1 line:', zip(npy.take(xdata, ind), npy.take(ydata, ind))


def on_pick(event):
    global globalCurrentDataset
    global globalDeleteDataMode
    global globalLassoDataMode

#    axes(ax1)
#    del(gca().patches[-1]) 
# #        del(linesData[-1])
#    del(maskRanges[-1])
#    draw()
#    print "removing mask range: "
#    printCurrentMaskRanges()
    mouseevent = event.mouseevent
    artist = event.artist

    if mouseevent.button == 1:
        if isinstance(event.artist, Line2D):  # if this is fn plot
            thisline = event.artist
            xdata = thisline.get_xdata()
            ydata = thisline.get_ydata()
            ind = event.ind
            print 'selected point data:', zip(np.take(xdata, ind), np.take(ydata, ind))
#        if isinstance(artist, AxesImage):
#            im = artist
#            A = im.get_array()
#            print 'onpick4 image', A.shape
#        if isinstance(event.artist, Rectangle):
#            patch = event.artist
#            print 'onpick1 patch:', patch.get_path()
        else:
            ind = event.ind
            thisline = event.artist
            xdata = thisline.get_offsets()[:, 0]
            ydata = thisline.get_offsets()[:, 1]
            r = thisline.get_sizes()
#            print r
#            print 'selected point data:', zip(np.take(xdata, ind), np.take(ydata, ind), ind , np.take(scatterGlobalData,ind))
#            print globalCurrentDataset
            print 'selected point data:', globalCurrentDataset[ind]
            if globalDeleteDataMode:
                 maskLastPickedDataPoint(event)


def on_mouse_button_down(event):
    
    global selectedLines
    global globalAxes, globalLassoManager
    
    
    if event.button == 2:
        if selectedLinesMode:
            deleteSelectedLine(event)
        else:
            deletePointedMaskRange(event)

    if event.button == 1:
        if globalLassoDataMode:
            globalLassoManager.onpress(event)
        
#    if event.button==3:
#        xdata = event.xdata
#        ydata = event.ydata
#        print "adding new point data: %lf, %fE" % (xdata,ydata)
#        selectedLines.append([xdata,ydata])
#        printSelectedLines()
#        t=array(selectedLines)
#        plot(t[:,0],t[:,1], 'yo', ms=4, label='selected lines', picker=True)
#        draw()

# def chooseFileNameWidget():
#    global saveToFile
#    global saveToFileWidget
#    fname=saveToFileWidget.get().strip()
#    np.savetxt(fname,lastBinnedSpectra)
#    print "saved spectra to file: "+fname
#    saveWidgetMain.destroy()
#    saveWidgetMain=None


def map_picker(qm, mouseevent):
    """
    find the points within a certain distance from the mouseclick in
    data coords and attach some extra attributes, pickx and picky
    which are the data points that were picked
    """
    global mapBbox
    global mapPaths
    global projectedMapData
    
    if projectedMapData == None:
        mapBbox = qm.get_datalim(gca().transData)
        mapPaths = qm.get_paths()
        projectedMapData = array([ np.mean(mapPaths[i].vertices[0:4], 0) for i in arange(len(mapPaths)) ])
#    m=array([ p[i].vertices[0] for i in arange(len(p)) ])
#    np.savetxt("dupa",m)
        
#    v2= np.mean(qm.get_paths()[1].vertices[0:4],0)
#    print v1
#    print v2
#    print v1-v2
    if mouseevent.xdata is None: return False, dict()
    xdata = projectedMapData[:, 0]
    ydata = projectedMapData[:, 1]
    maxd = mapBbox.width / 300
    d = np.sqrt((xdata - mouseevent.xdata) ** 2. + (ydata - mouseevent.ydata) ** 2.)

    ind = np.nonzero(np.less_equal(d, maxd))
    if len(ind):
        pickx = np.mean(np.take(xdata, ind))
        picky = np.mean(np.take(ydata, ind))
        props = dict(ind=ind, pickx=pickx, picky=picky)
        return True, props
    else:
        return False, dict()

def onMapPick(event):
    global projectMap
    global globalCurrentDatasetIdx
    global inFile
    x = event.pickx
    y = event.picky
    print 'onMapPick line:', x, y
    l, b = projectMap(x, y, inverse=True)
    l = 360 - l
    if l > 360:
        l = l - 360
    print 'converted coordinates: (l,b)=(%lf, %lf)' % (l, b)
    
    # find the indexes of the data in the lon and lat data for the current plot
    # find lon index
    lonIdx = min(range(len(inFile[globalCurrentDatasetIdx][1])), key=lambda i: abs(inFile[globalCurrentDatasetIdx][1][i] - l))
    print lonIdx, inFile[globalCurrentDatasetIdx][1][lonIdx]
    # find lat index
    latIdx = min(range(len(inFile[globalCurrentDatasetIdx][2])), key=lambda i: abs(inFile[globalCurrentDatasetIdx][2][i] - b))
    print latIdx, inFile[globalCurrentDatasetIdx][2][latIdx]
    print "data lon: %f, lat: %f" % (inFile[globalCurrentDatasetIdx][1][lonIdx], inFile[globalCurrentDatasetIdx][2][latIdx])
    lonIdx = -lonIdx
#    lonIdx=lonIdx+round(len(inFile[globalCurrentDatasetIdx][1])/2)
    if lonIdx > len(inFile[globalCurrentDatasetIdx][1]):
        lonIdx = lonIdx - len(inFile[globalCurrentDatasetIdx][1] - 1)
    print "lonIdx: %i" % lonIdx
    print np.shape(inFile[globalCurrentDatasetIdx][0])
    print "data value: %f" % (inFile[globalCurrentDatasetIdx][0][latIdx][lonIdx])



###### LASSO SELECT TOOL


class Datum(object):
#    colorin = mcolors.to_rgba("red")
#    colorout = mcolors.to_rgba("blue")
#     colorin='k'
#     colorout='b'

    def __init__(self, x, y, include=False):
        self.x = x
        self.y = y
#         if include:
#             self.color = self.colorin
#         else:
#             self.color = self.colorout


class LassoManager(object):
    def __init__(self, ax, data):
        self.axes = ax
        self.canvas = ax.figure.canvas
        self.data = data

        self.Nxy = len(data)

#         facecolors = [d.color for d in data]
        self.xys = [(d.x, d.y) for d in data]
#         self.collection = RegularPolyCollection(
#             6, sizes=(100,),
#             facecolors=facecolors,
#             offsets=self.xys,
#             transOffset=ax.transData)

#         ax.add_collection(self.collection)

#         self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)

    def callback(self, verts):
#         facecolors = self.collection.get_facecolors()
        p = path.Path(verts)
        ind = p.contains_points(self.xys)
        print ind
        maskLastLassoDataPoints(ind)
#         print self.xys
#        for i in range(len(self.xys)):
#            if ind[i]:
#                facecolors[i] = Datum.colorin
#            else:
#                facecolors[i] = Datum.colorout

        self.canvas.draw_idle()
        self.canvas.widgetlock.release(self.lasso)
        del self.lasso

    def onpress(self, event):
        if self.canvas.widgetlock.locked():
            return
        if event.inaxes is None:
            return
        self.lasso = Lasso(event.inaxes,
                           (event.xdata, event.ydata),
                           self.callback)
        # acquire a lock on the widget drawing
        self.canvas.widgetlock(self.lasso)



############################################################################################
############################################################################################
############################################################################################
def formatXaxisDates(figIdx):
    import matplotlib.dates as mdates
    seconds = mdates.SecondLocator()
    minutes = mdates.MinuteLocator()
    hours = mdates.HourLocator()
    days = mdates.DayLocator()
    years = mdates.YearLocator()  # every year
    months = mdates.MonthLocator()  # every month
#                yearsFmt = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')

    ax = gca()
    print(option.dateFmtPlot)
#     yearsFmt = mdates.DateFormatter("%Y-%m-%d\n%H:%M:%S")
    yearsFmt = mdates.DateFormatter(option.dateFmtPlot.decode("unicode_escape"))
    ax.xaxis.set_major_formatter(yearsFmt)
    # format the ticks
    ax.xaxis.set_major_locator(days)
    ax.xaxis.set_major_formatter(yearsFmt)
#                 ax.xaxis.set_minor_locator(hours)
#                 ax.xaxis.set_minor_locator(None)

    if option.years:
        ax.xaxis.set_major_locator(years)
        if option.minorTicks:
            ax.xaxis.set_minor_locator(months)
    if option.months:
        ax.xaxis.set_major_locator(months)
        if option.minorTicks:
            ax.xaxis.set_minor_locator(days)
    if option.days:
        ax.xaxis.set_major_locator(days)
        if option.minorTicks:
            ax.xaxis.set_minor_locator(hours)
    if option.hours:
        ax.xaxis.set_major_locator(hours)
        if option.minorTicks:
            ax.xaxis.set_minor_locator(minutes)
    if option.minutes:
        ax.xaxis.set_major_locator(minutes)
        if option.minorTicks:
            ax.xaxis.set_minor_locator(seconds)

#                datemin = datetime.date(2013, 1, 1)
#                datemax = datetime.date(2015, 1, 1)
    datemin = ''
    datemax = ''
    global globalCurrentDatasetX
    if globalCurrentDatasetX != None:
        datemin = globalCurrentDatasetX[0]
        datemax = globalCurrentDatasetX[-1]
        print datemin
        print datemax
    print option.datexmin[figIdx % len(option.datexmin)]
    print option.datexmax
    if option.datexmin[figIdx % len(option.datexmin)] != -1:
        datemin = datetime.datetime.strptime(option.datexmin[figIdx % len(option.datexmin)], option.dateFmt)
    if option.datexmax[figIdx % len(option.datexmax)] != -1:
        datemax = datetime.datetime.strptime(option.datexmax[figIdx % len(option.datexmax)], option.dateFmt)

    if datemin != '' and datemax != '':
        print 'datemin,datemax: ', datemin, datemax
        ax.set_xlim(datemin, datemax)
#                ax.format_xdata = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
#                fig.autofmt_xdate()

    ymin, ymax = ax.get_ylim()
    if option.ymin[figIdx % len(option.ymin)] != None:
        ymin = option.ymin[figIdx % len(option.ymin)]
    if option.ymax[figIdx % len(option.ymax)] != None:
        ymax = option.ymax[figIdx % len(option.ymax)]
    print 'setting y limits: %f, %f' % (ymin, ymax)
#                ylim([ymin,ymax])
    ax.set_ylim(ymin, ymax)
    print
    print 'Setting xticklabels rotation: ', option.Rxlabels    
    print
    setp(ax.get_xticklabels(), rotation=option.Rxlabels, horizontalalignment='center', verticalalignment='top') 
    setp(ax.get_yticklabels(), rotation=option.Rylabels) 
    xticks(rotation=option.Rxlabels)
    


def getPalette():
    palette = None

    if option.palette[i % len(option.palette)] == "jet":
        palette = matplotlib.cm.jet
    elif option.palette[i % len(option.palette)] == "gray":
        palette = matplotlib.cm.gray
    elif option.palette[i % len(option.palette)] == "hot":
        palette = matplotlib.cm.hot
    elif option.palette[i % len(option.palette)] == "afmhot":
        palette = matplotlib.cm.afmhot
    elif option.palette[i % len(option.palette)] == "spectral":
        palette = matplotlib.cm.Spectral
    elif option.palette[i % len(option.palette)] == "blues":
        palette = matplotlib.cm.Blues
    elif option.palette[i % len(option.palette)] == "seismic":
        palette = matplotlib.cm.seismic
    else:
        palette = matplotlib.cm.gist_yarg

    return palette

def filterRows(dataset, i):
#     if len(option.every) > 0:
#         dataset = dataset[::option.every[i % len(option.every)]]
    return dataset
    

def filterData(data, dataFilter, colNames=""):
    for f in dataFilter.split(','):
        print
        print "applying data filter: %s" % f
        filterType = ""
        filter = f.split('==')
        if len(filter) != 2:
            filter = f.split('>=')
            if len(filter) != 2:
                filter = f.split('<=')
                if len(filter) != 2:
                    filter = f.split('<')
                    if len(filter) != 2:
                        filter = f.split('>')
                        if len(filter) != 2:
                            print "unknown filter format, skipping"
                        else:
                            filterType = ">"
                    else:
                        filterType = "<"
                else:
                    filterType = "<="
            else:
                filterType = ">="
        else:
            filterType = "=="
        
        if filterType == "==" or filterType == "<=" or filterType == ">=" or filterType == ">" or filterType == "<" :
    
            filterName = filter[0]
            filterValue = float(filter[1])
            print "filterType: %s" % filterType
            print "filterName: %s" % filterName
            print "filterValue: %s" % filterValue
            print "searching for the requested column name"
            colIdx = -1
            i = 0
            if type(colNames) == type(list()):
                for colName in colNames:
                    if colName == filterName:
                        colIdx = i
                    i = i + 1
                if colIdx == -1:
                    print "error, filter name not found"
                else:
                    print "ok, filter name found in column %i" % colIdx
                    if filterType == "==":
                        data = data[data[:, colIdx] == filterValue]
                    if filterType == ">=":
                        data = data[data[:, colIdx] >= filterValue]
                    if filterType == "<=":
                        data = data[data[:, colIdx] <= filterValue]
                    if filterType == "<":
                        data = data[data[:, colIdx] < filterValue]
                    if filterType == ">":
                        data = data[data[:, colIdx] > filterValue]
                        
                            
    
        print "filtered data has now: %i rows" % len(data)
        
        
        
    return data

def performXdataOperations(datax, i):
    if len(Xoperations) > 0:
        if Xoperations[i % len(Xoperations)] == '+':
#            print type(datax[0])
#            print type(datetime.datetime(2013,1,1))
#            s='string'
#            print type(s)
            if type(datax[0]) == type(datetime.datetime(2013, 1, 1)):
                for ii in range(len(datax)):
                    datax[ii] = datax[ii] + datetime.timedelta(seconds=XoperationValues[i % len(XoperationValues)])
            else:                            
                datax = datax + XoperationValues[i % len(XoperationValues)]
        if Xoperations[i % len(Xoperations)] == '-':
            datax = datax - XoperationValues[i % len(XoperationValues)]
        if Xoperations[i % len(Xoperations)] == '*':
            datax = datax * XoperationValues[i % len(XoperationValues)]
        if Xoperations[i % len(Xoperations)] == '/':
            datax = datax / XoperationValues[i % len(XoperationValues)]
    return datax

def performYdataOperations(datay, i):
    if len(Yoperations) > 0:
        print "performing Yoperations for dataset " + str(i)
        if Yoperations[i % len(Yoperations)] == '+':
            datay = datay + YoperationValues[i % len(YoperationValues)]
        if Yoperations[i % len(Yoperations)] == '-':
            datay = datay - YoperationValues[i % len(YoperationValues)]
        if Yoperations[i % len(Yoperations)] == '*':
            datay = datay * YoperationValues[i % len(YoperationValues)]
        if Yoperations[i % len(Yoperations)] == '/':
            datay = datay / YoperationValues[i % len(YoperationValues)]
    return datay
    
def getHDF5dsetNames(fname):
    f = h5py.File(fname, 'r')
    k = h5py.File.keys(f)
    f.close()        
    return k

def plotPOVray_file(inFile):
    povfile = '''
        #include "metals.inc"
        
        #declare TWOPI = 6.283185307179586476925287;
        #declare RADIUS = 1;
        #declare NX = ;
        #declare NY = ;
        #declare NZ = ;
        #declare DD = <NX,NY,NZ>;
        #declare CC = DD / 2;
        #declare VP = <0,0,-NZ>;
        
        global_settings { 
            ambient_light <1,1,1> 
            assumed_gamma 1
        }
        
        camera {
           location VP
           up y
           right x
           angle FOV
           //sky <0,0,-1>
           look_at <0,0,0>
        }
        
        light_source {
           VP + <0,0,2*NZ>
           color rgb <1,1,1>
           media_interaction on
           media_attenuation on
           shadowless
        }
        
        #declare theinterior = interior {
           media {
                intervals 100
                ratio 0.5
                samples 2,2
                method 2
                emission <1,1,1> / 100
                absorption <1,1,1> / 1000
                scattering { 1, <0,0,0> }
                confidence 0.999
                variance 1/1000
                density {
                    density_file df3 
                    interpolate 1
                    color_map {
                        [0.00 rgb <0,0,0>]
                        [0.50 rgb <0,0,1>]
                        [1.00 rgb <1,1,1>]
                        [1.00 rgb <1,0,0>]
                    }
                }
            }
        }
        
        box {
            <0,0,0>, <1,1,1>
            pigment { rgbf 1 }
            interior { theinterior }
            hollow
            translate <-0.5,-0.5,-0.5>
            scale DD
            rotate <0,0,360*clock>
        }
        
        #declare bbox = texture {
           pigment { rgb <0.5,0.5,0.5> }
           finish { F_MetalB }
        }
        
        /* Corners of box */
        union {
            sphere { <0,0,0>, RADIUS texture {bbox} }
            sphere { <0,NY,0>, RADIUS texture {bbox} }
            sphere { <NX,NY,0>, RADIUS texture {bbox} }
            sphere { <NX,0,0>, RADIUS texture {bbox} }
            sphere { <0,0,NZ>, RADIUS texture {bbox} }
            sphere { <0,NY,NZ>, RADIUS texture {bbox} }
            sphere { <NX,NY,NZ>, RADIUS texture {bbox} }
            sphere { <NX,0,NZ>, RADIUS texture {bbox} }
            translate -CC
        }
        
        /* Main border of box */
        union {
            cylinder { <0,0,0>, <NX,0,0>, RADIUS texture {bbox} }
            cylinder { <0,0,0>, <0,NY,0>, RADIUS texture {bbox} }
            cylinder { <0,0,0>, <0,0,NZ>, RADIUS texture {bbox} }
            cylinder { <NX,NY,NZ>, <0,NY,NZ>, RADIUS texture {bbox} }
            cylinder { <NX,NY,NZ>, <NX,0,NZ>, RADIUS texture {bbox} }
            cylinder { <NX,NY,NZ>, <NX,NY,0>, RADIUS texture {bbox} }
            cylinder { <0,0,NZ>, <NX,0,NZ>, RADIUS texture {bbox} }
            cylinder { <0,0,NZ>, <0,NY,NZ>, RADIUS texture {bbox} }
            cylinder { <NX,NY,0>, <NX,0,0>, RADIUS texture {bbox} }
            cylinder { <NX,NY,0>, <0,NY,0>, RADIUS texture {bbox} }
            cylinder { <0,NY,0>, <0,NY,NZ>, RADIUS texture {bbox} }
            cylinder { <NX,0,0>, <NX,0,NZ>, RADIUS texture {bbox} }
            translate -CC
        }

    '''

    povfile = povfile.replace('NX = ', 'NX = %i' % fieldPixNuml[0 % len(fieldPixNuml)])
    povfile = povfile.replace('NY = ', 'NY = %i' % fieldPixNuml[1 % len(fieldPixNuml)])
    povfile = povfile.replace('NZ = ', 'NZ = %i' % fieldPixNuml[2 % len(fieldPixNuml)])
    povfile = povfile.replace('#declare VP = <0,0,-NZ>;', '#declare VP = <%f,%f,%f>;' % (option.PRcamX, option.PRcamY, option.PRcamZ,))
    
    
    povfile = povfile.replace('density_file df3', 'density_file df3 "%s"' % inFile)
    povfile = povfile.replace('angle FOV', 'angle %f' % option.PRfov)
    
#    print povfile
    outf = open("render.pov", 'w')
    outf.write(povfile)
    outf.close()
    guiflag = ' -d '
    guiholdflag = ''
    if option.PRgui:
        guiflag = ''
    if option.PRhold:
        guiholdflag = ' +P '
    cmd = 'povray %s %s -W%i -H%i render.pov -O%s.png' % (guiflag, guiholdflag , option.PRwidth, option.PRheight, option.outputFile)
    cpedsPythCommon.sayAndExecute("command", cmd, 1)


def plotCircles(ax, circlePatchDef):
    circ = Circle((circlePatchDef[0], circlePatchDef[1]), circlePatchDef[2], fc=None, ec='k', color=None, fill=False)
    ax.add_patch(circ)
    



def getDataMinMaxValues2(datax, datay):
    minx = min(datax)
    maxx = max(datax)
    miny = min(datay)
    maxy = max(datay)
    return minx, maxx, miny, maxy




def getDataMinMaxValues(ax, datax, datay, datasetNo):
    if option.polar:
        minx, maxx, miny, maxy = getDataMinMaxValues2(datax, datay)    
        return minx, maxx, miny, maxy
        
    xmin, xmax, ymin, ymax = ax.axis()
    print "getDataMinMaxValues: ax"
    print xmin, xmax, ymin, ymax

    minx, maxx, miny, maxy = getDataMinMaxValues2(datax, datay)
    if datasetNo == 0:
        print 'minx,maxx,miny,maxy: data'
        print minx, maxx, miny, maxy
        return minx, maxx, miny, maxy
    
    print "getDataMinMaxValues: data"
    print minx, maxx, miny, maxy
    if minx < xmin:
        xmin = minx;
    if miny < ymin:
        ymin = miny;
    if maxx > xmax:
        xmax = maxx;
    if maxy > ymax:
        ymax = maxy;
    
    
    return xmin, xmax, ymin, ymax;


def loadDataFromUDP(UDPparams):
    UDPparams = UDPparams.split("=")
    UDPaddr = UDPparams[1].split(":")
    if len(UDPaddr) != 4:
        print "WRONG UDP format. Should be UDP=host:port:multicast_flag:number_of_UDP_packets_to_read"
        print "where:"
        print "multicast_flag = 0 for non-multicast; 1 - for multi-cast"
        print "number_of_UDP_packets_to_read - number of UDP packages to read from the socket before plotting"
        sys.exit(-1)
    HOST = UDPaddr[0]
    PORT = UDPaddr[1]
    isMULTICAST = False 
    multi = UDPaddr[2]
    if multi == '1':
        isMULTICAST = True
    udpCount = int(UDPaddr[3])
    
    
    print 'listening to port: ', PORT
    print 'interface: ', HOST
    print 'multicast: ', isMULTICAST
    print 'UDP datagrams to read:', udpCount
#     while 1:
    bindata = cpedsPythCommon.getUDPdatagram(HOST, int(PORT), udpCount, isMULTICAST)
    adata = list()
    for dgram in bindata:
#         print dgram
        tmp = dgram.strip().split(' ')
        adata.append(asarray(tmp, dtype='float'))
    adata = asarray(adata, dtype='float').reshape((1, -1))
#     print adata
    return adata 
    
def loadDataFromFileStd(fname, startFrom=0, rowsCount=0, loadEvery=1, binSamples=-1, colx=0, coly=1, badX='None', badY='None'):
    d = np.loadtxt(args[i])
    print startFrom, rowsCount, loadEvery
    if rowsCount == 0:
        d = d[startFrom::loadEvery]
    else:
        d = d[startFrom:rowsCount:loadEvery]
        
    if badY!='None':
        d=d[d[:,coly]!=float(badY)]
    
    
    dbin = list()
    if binSamples > 1:
        for ci in range(len(d[0])):
            data = d[:, ci]
            slices = np.linspace(0, len(data), binSamples + 1, True).astype(np.int)
            counts = np.diff(slices)
            mean = np.add.reduceat(data, slices[:-1]) / counts
            dbin.append(mean)
        d = np.asanyarray(dbin).T
    print fname
    print d
    return d
    
def loadDataFromFile(fname, colx=0, coly=1, startFrom=0, rowsCount=-1, loadEvery=1, binSamples=-1, badX='None', badY='None'):
    
    l = list()
    lineNo = 0
    readRows = 0
#    rowsCount*=loadEvery
    global globalPlotTypeIdx

    if option.binDAQd:
        print 'loading data from file using y-column: %i and x column: %li' % (coly, colx)
        import struct
        infile = open (fname, 'rb')
        if rowsCount == -1:
            print 'loading 10000 records from input file'
            rowsCount = 10000

        for line in arange(rowsCount):
            line = infile.read(104)
            line = struct.unpack('q12d', line)
            if lineNo >= startFrom:
#                print line
                if colx == -1:
                    l.append([lineNo, line[coly]])
                else:
                    l.append([line[colx], line[coly]])
                readRows = readRows + 1
            lineNo = lineNo + 1
        
    else:    
        infile = open (fname, 'r')
        loadingDataTimeWithSpace = 0
        if option.plotType[globalPlotTypeIdx] == 'ts' or " " in option.dateFmt:
            loadingDataTimeWithSpace = option.dateFmt.count(' ')
            dateFmt=option.dateFmt.split()
            print 'Loading time sequence data with space[s] in format - will load %i columns (%i:%i and %i)' % (loadingDataTimeWithSpace + 1, colx, colx + loadingDataTimeWithSpace, coly)
            
        print 'loading data from file using y-column: %i and x column: %li' % (coly, colx)
        if rowsCount == -1:
            for line in infile:
                if lineNo >= startFrom:
                    if not re.match('#', line) and len(line) > 1:
                        if (lineNo - startFrom) % loadEvery == 0:
                            line = line.strip()
                            sline = line.split()
#                             print sline
                            if colx == -1:
                                l.append([lineNo, sline[coly]])
                            else:
                                if loadingDataTimeWithSpace > 0:
#                                     dt=[ str(int(float(tmp))) for tmp in sline[colx:colx+loadingDataTimeWithSpace+1] ]
                                    dt = [ tmp for tmp in sline[colx:colx + loadingDataTimeWithSpace + 1] ]
                                    sline[colx] = ' '.join(dt)
                                l.append([sline[colx], sline[coly]])
                            readRows = readRows + 1
                lineNo = lineNo + 1
        else:
            dtstr=lambda fmt,dt: str(float(dt)) if '.' in fmt else str(int(float(dt)))
            for line in infile:
                if lineNo >= startFrom:
                    if not re.match('#', line):
                        if (lineNo - startFrom) % loadEvery == 0:
                            line = line.strip()
                            sline = line.split()
                            if colx == -1:
                                l.append([lineNo, sline[coly]])
                            else:
                                if loadingDataTimeWithSpace > 0:
#                                     dt = [ str(int(float(tmp))) for tmp in sline[colx:colx + loadingDataTimeWithSpace + 1] ]
                                    dt = [ dtstr(dateFmt[f],tmp) for f,tmp in enumerate(sline[colx:colx + loadingDataTimeWithSpace + 1]) ]
                                    sline[colx] = ' '.join(dt)
                                l.append([sline[colx], sline[coly]])
                            readRows = readRows + 1
                lineNo = lineNo + 1
                if readRows == rowsCount: break
            
        infile.close()
        print "loaded %i lines" % readRows
    
#    print 'selecting every %i' % loadEvery
    loaded = array(l)
#    loaded=loaded[::loadEvery]
    print "selected %i lines" % len(loaded)
    
    if badY!='None':
        loaded=loaded[loaded[:,1]!=badY]
    
    
    if binSamples > 1:
        l1 = np.array(loaded[:, 0], dtype=float)
        l2 = np.array(loaded[:, 1], dtype=float)
        idx = 0
        loadedl = list()
        while idx < len(l1):
            binX = l1[idx:idx + binSamples].mean()
            binY = l2[idx:idx + binSamples].mean()
            loadedl.append([binX, binY])
            idx = idx + binSamples
        loaded = array(loadedl)
        
        
    
    
    return loaded



def ticksFormatter(x, pos):
    return '%.2E' % (x)


def XTicksFormatter(x, pos):
    'The two args are the value and tick position'
    Nx = size(boxSlice[:, 0])
#    x1=0.5; x2=Nx+0.5;
    x1 = 0.0; x2 = Nx + 0.0;
    y1 = option.xmin[0]
    y2 = option.xmax[0]
    return '%.2f' % ((x - x1) * (y2 - y1) / (x2 - x1) + y1)

def YTicksFormatter(x, pos):
    'The two args are the value and tick position'
    Ny = size(boxSlice[0])
#    x1=0.5; x2=Ny+0.5;
    x1 = 0.0; x2 = Ny + 0.0;
    y1 = option.ymin[0]
    y2 = option.ymax[0]
    return '%.2f' % ((x - x1) * (y2 - y1) / (x2 - x1) + y1)


def ZTicksFormatter(x, pos):
    'The two args are the value and tick position'
    Ny = size(boxSlice[0])
#    x1=0.5; x2=Ny+0.5;
    x1 = 0.0; x2 = Ny + 0.0;
    y1 = option.ymin[0]
    y2 = option.ymax[0]
    return '%.2f' % ((x - x1) * (y2 - y1) / (x2 - x1) + y1)

# def getMatInfo(boxSlice):
#    vmin=nanmin(boxSlice)
#    vmax=nanmax(boxSlice)
#    maxabs=max([abs(vmin),abs(vmax)])
#    print "masked statistics:"
#    print "min: %E" % vmin
#    print "max: %E" % vmax
#    print "max-min: %E" % (vmax-vmin)
#    print maxabs
#    return vmin,vmax,maxabs


# def getCol(i, s):
#    reader = csv.reader(s, skipinitialspace=True)
#    j=0;
#    print "getCOl"
#    for r in reader:
# #        print int(r[0])
#        if i==j:
#            return int(r[0])
#        j=j+1
#    return False

def removeNans3(datax, datay, datac):
    print "removing nans3"
    tmp = np.array(np.concatenate([[datax], [datay], [datac]], axis=0).T,dtype=np.float64)
    print "data length before nan removal: %li" % len(datax)
#    print tmp
    tmp = tmp[~np.isnan(tmp).any(axis=1)]
#    print tmp
    datax = tmp[:, 0]
    datay = tmp[:, 1]
    datac = tmp[:, 2]
    print "data length after nan removal: %li" % len(datax)
    return datax, datay, datac

def removeNans2(datax, datay):
    print "removing nans2"
    print "data length before nan removal: %li" % len(datax)
    tmp = np.array(np.concatenate([[datax], [datay]], axis=0).T,dtype=np.float64)
#    print tmp
    tmp = tmp[~np.isnan(tmp).any(axis=1)]
#    print tmp
    datax = tmp[:, 0]
    datay = tmp[:, 1]
    print "data length after nan removal: %li" % len(datax)
    return datax, datay

def removeNans1(datax):
    print "removing nans1"
    print "data length before nan removal: %li" % len(datax)
    tmp = np.asarray([datax], dtype='float').T
    tmp = tmp[~np.isnan(tmp).any(axis=1)]
    datax = tmp
    print "data length after nan removal: %li" % len(datax)
    return datax

def removeNegativesValues(datax, y1, y2):
    print "removing negative values"
    tmp = np.concatenate([[datax], [y1], [y2]], axis=0).T
#    print tmp
    tmp = tmp[tmp[:, 1] > 0]
    tmp = tmp[tmp[:, 2] > 0]
    return tmp[:, 0], tmp[:, 1], tmp[:, 2]







def makeMayaViPlot(dataList):
    print len(dataList)
    for i in np.arange(len(dataList)):
        print 'plotting dataset: ', i
        data = dataList[i]
        data *= option.MdataScale
        vvals = cpedsPythCommon.getFloatList(option.Mvrange)
        if vvals[0] == 0 and vvals[1] == 0:
            vvals[0] = np.amin(data)
            vvals[1] = np.amax(data)
            
        print 'vvals: ', vvals
    
        ###############################################################################################
        # mayaSurf type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaSurf':
    #        print data
    #        mlab.pipeline.volume(mlab.pipeline.scalar_scatter(data))
            mlab.surf(data, vmin=vvals[0], vmax=vvals[1])
    #        mlab.scalar_scatter(data)
            if option.colorbar:
                mlab.colorbar()
    
        ###############################################################################################
        # mayaScat type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaScat':
            print data
    #        mlab.pipeline.volume(mlab.pipeline.scalar_scatter(data))
            colx = option.colx[i % len(option.colx)]
            coly = option.coly[i % len(option.coly)]
            colz = option.colColor[i % len(option.colColor)]
            if size(data[0]) >= 4:
                colSize = option.colSize # [i % len(option.colSize)]

#                 mlab.points3d(data[:, 0], data[:, 1], data[:, 2], data[:, 3] / np.amax(data[:, 3]), scale_factor=option.mayaPS)
                mlab.points3d(data[:, colx], data[:, coly], data[:, colz], data[:, colSize] / np.amax(data[:, colSize]), scale_factor=option.mayaPS)
            else:
#                 mlab.points3d(data[:, 0], data[:, 1], data[:, 2], np.ones((len(data[:, 0]))), scale_factor=option.mayaPS)
                mlab.points3d(data[:, colx], data[:, coly], data[:, colz], np.ones((len(data[:, 0]))), scale_factor=option.mayaPS)
    #        mlab.scalar_scatter(data)
            if option.colorbar:
                mlab.colorbar()
        
        
        ###############################################################################################
        # mayaVolCut type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaVolCut':
            mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data), plane_orientation='x_axes', slice_index=len(data) / 2, vmin=vvals[0], vmax=vvals[1], colormap=option.palette[i % len(option.palette)])
            mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data), plane_orientation='y_axes', slice_index=len(data) / 2, vmin=vvals[0], vmax=vvals[1], colormap=option.palette[i % len(option.palette)])                        
            mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(data), plane_orientation='z_axes', slice_index=len(data) / 2, vmin=vvals[0], vmax=vvals[1], colormap=option.palette[i % len(option.palette)])
            if option.colorbar:
                mlab.colorbar()
    
        ###############################################################################################
        # mayaBarChart type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaBarChart':
    
            mlab.barchart(data, vmin=vvals[0], vmax=vvals[1])
            if option.colorbar:
                mlab.colorbar()

        ###############################################################################################
        # mayaVol type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaContour':
    #        data/=np.amax(data)
    #        print data
            contList = mayaContours
            mlab.contour3d(data, vmin=vvals[0], vmax=vvals[1], colormap=option.palette[i % len (option.palette)], contours=contList)
            if option.colorbar:
                mlab.colorbar()
    #            hold(True)
    
        ###############################################################################################
        # mayaVol type plot
        ###############################################################################################
        if option.plotType[i] == 'mayaVol':
    #        data/=np.amax(data)
    #        print data
            
            mlab.pipeline.volume(mlab.pipeline.scalar_field(data), vmin=vvals[0], vmax=vvals[1])
            if option.colorbar:
                mlab.colorbar()
    #            hold(True)

    if option.noAxes == False:
        if option.xmin[0] != None or option.xmax[0] != None or option.ymin[0] != None or option.ymax[0] != None or option.zmin != None or option.zmax != None:
            extentxyz = [option.xmin[0], option.xmax[0], option.ymin[0], option.ymax[0], option.zmin, option.zmax]
            mlab.axes(extent=extentxyz, xlabel=option.xlabel[i % len(option.xlabel)], ylabel=option.ylabel[i % len(option.ylabel)], zlabel=option.zlabel[i % len(option.zlabel)])
        else:
            extentxyz = None

    if option.save:
        if option.outputFile != '':
            fname = option.outputFile
            mlab.savefig(fname, size=cpedsPythCommon.getFloatList(option.figSize))
        else:
            print "no output file name given, will not save"
    else:
        mlab.show()
#                mlab.outline()

    

def mkGrid(datax, datay, gridNbins):
    Nbinx = gridNbins[0]
    print 'Nbinx', Nbinx
    xmin = np.amin(datax)
    xmax = np.amax(datax)
    ymin = np.amin(datay)
    ymax = np.amax(datay)
    binsx = np.linspace(xmin, xmax, Nbinx)
    if option.logX:
        binsx = np.logspace(np.log10(xmin), np.log10(xmax), Nbinx)
    Nbiny = gridNbins[1]
    print 'Nbiny', Nbiny
    binsy = np.linspace(ymin, ymax, Nbiny)
    if option.logY:
        binsy = np.logspace(np.log10(ymin), np.log10(ymax), Nbiny)

    return binsx, binsy


# usePlotLimits - if true then the grid will be defined based on the --xmin --xmax --ymin --ymax values and not
# based on the data ranges
def calculatePointsDensityDistribution(datax, datay, gridNbins, usePlotLimits=False):
    Nbinx = gridNbins[0]
    Nbiny = gridNbins[1]
    if usePlotLimits:
        binsx = Nbinx
        if option.logX:
            binsx = np.logspace(np.log10(np.amin(datax)), np.log10(np.amax(datax)), Nbinx)
        else:
            binsx = np.linspace(option.xmin, option.xmax, Nbinx)
        if option.logY:
            binsy = np.logspace(np.log10(np.amin(datay)), np.log10(np.amax(datay)), Nbiny)        
        else:
            binsy = np.linspace(option.ymin, option.ymax, Nbiny)
    else:
        binsx = Nbinx
        if option.logX:
            binsx = np.logspace(np.log10(np.amin(datax)), np.log10(np.amax(datax)), Nbinx)
        binsy = Nbiny
        if option.logY:
            binsy = np.logspace(np.log10(np.amin(datay)), np.log10(np.amax(datay)), Nbiny)
            
    H, xedges, yedges = np.histogram2d(datax, datay, bins=[binsx, binsy])
    dataxtmp = list()
    dataytmp = list()
    scatterGlobalData = list()
    for i in range(len(H[:, 0])):
        for j in range(len(H)):
            scatterGlobalData.append(H[i, j])
    
    datax = xedges[0:-1]
    datay = yedges[0:-1]
    H = H.T
    
    return datax, datay, H





# @param boxSlice - array with slices
# @param sliceNo - index of the slice
# @param sliceName - file name prefix of the slice
def makeFunctionPlot(inFile):
    global plotAxes
    global plotFig
    global continuousPottingFirstTime
    global globalCurrentDataset
    global globalCurrentDatasetX
    global globalCurrentDatasetY
    global globalNfiles
    Nfiles = len(inFile)
    globalNfiles = Nfiles
    
        
#    Nx=size(boxSlice[:,0])
#    Ny=size(boxSlice[0])
    if option.figSize == "A4":
        option.figSize = '11.6929,8.2677'

#    reader = csv.reader([option.figSize], skipinitialspace=True, delimiter=',')
#    figSize=list()
        
#    for r in reader:
#        figSize.append(r)
#    if len(figSize)==1:
#        figSize.append(figSize[0])
#    fig=figure(figsize=(float(figSize[0][0]),float(figSize[0][1])), dpi=option.DPIgui)
    
    figSize = cpedsPythCommon.getFloatList(option.figSize)
    print figSize

    if plotFig == None:
        fig = figure(figsize=(float(figSize[0]), float(figSize[1 % len(figSize)])), dpi=option.DPIgui)
        plotFig = fig
    else:
        fig = plotFig
        
    fig.set_facecolor(option.bgcolor)


    # plot all datasets in the same plot by default
    if option.dsPerPlot == -1:
        option.dsPerPlot = Nfiles
        
    if option.plotDs != '':
        globalPlotDs = cpedsPythCommon.getIntList(option.plotDs)
    
#    figCols=option.figCols
#    figRows=option.figRows
    
    figCols = ast.literal_eval(option.figCols)
    figRows = ast.literal_eval(option.figRows)
    print 'type(figCols): ', type(figCols)
    print 'type(figRows): ', type(figRows)
    if type(figCols) == type(list()) and type(figRows) == type(list()):
        print 'setting figCols and figRows'
        option.figCols = figCols
        option.figRows = figRows
        figCols = len(option.figCols)
        figRows = len(option.figRows)

    if type(figCols) == type(list()):
        print 'setting figCols'
        option.figCols = figCols
        figCols = len(option.figCols)
        print 'type(figCols): ', type(figCols)
        print 'type(figRows): ', type(figRows)
        if type(figRows) == type(int):
            print 'setting figRows2'
            option.figRows = list(figRows)

    if type(figRows) == type(list()):
        print 'setting figRows'
        option.figRows = figRows
        figRows = len(option.figRows)
        print 'type(figCols): ', type(figCols)
        print 'type(figRows): ', type(figRows)
        if type(figCols) == type(int):
            option.figCols = list([figCols])
            print 'setting fig cols2'
    
    
#    figRows=int( np.floor( float(np.round(float(Nfiles) / option.dsPerPlot))/float(figCols) )  )
    if figRows == 0:
        figRows = 1
#    print "np.round(float(Nfiles) / option.dsPerPlot): %lf" % np.round(float(Nfiles) / option.dsPerPlot)
#    print " %.15lf" % (float(np.round(float(Nfiles) / option.dsPerPlot))/float(figCols))
    print "figure plot cols: ", figCols
    print "figure plot rows: ", figRows
    print "figure plot cols: ", option.figCols
    print "figure plot rows: ", option.figRows
    print "figure datasets per plot: %i" % option.dsPerPlot
    print "figure datasets distribution: %s" % option.plotDs
    if type(option.figCols) == type(list()) or  type(option.figRows) == type(list()):
        print 'adjusting subplots'
        globalSubPlotGridSpec = gridspec.GridSpec(figRows, figCols, width_ratios=option.figCols, height_ratios=option.figRows)
#     sys.exit(1)

#    sys.exit()
#    figPlotNum=figRows*100+figCols*10+1
    
#    fig=figure(figsize=(12,12), dpi=option.DPIgui)
    if plotAxes == None:
        if option.polar:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
            print "Making polar axes"
        elif option.proj3d:
            ax = fig.add_subplot(111, projection='3d')
        else:
            ax = subplot(111)
    #        fig.subplots_adjust(left=option.border_left, right=1-option.border_right, top=1-option.border_top, bottom=option.border_bottom)
            fig.subplots_adjust(left=option.border_left, right=1 - option.border_right, top=1 - option.border_top, bottom=option.border_bottom, hspace=option.border_vspace, wspace=option.border_hspace)
        plotAxes = ax
    else:
        ax = plotAxes

    if continuousPottingFirstTime == False:
        ax.set_autoscale_on(True)

    if option.continuousPotting and continuousPottingFirstTime:
        continuousPottingFirstTime = False
        ion()
        show()
        return fig, ax

#     tight_layout(1.08, 1, 1, None)
#     ax.set_autoscale_on(False)
    # ax1=fig.add_subplot(2,2,1)

    #
    # define columns to plot
    #
#    reader = csv.reader(option.colx, skipinitialspace=True)
#    Xcols=list()
#    for r in reader:
#        if (r[0]!=''): Xcols.append(int(r[0]))
#    reader = csv.reader(option.coly, skipinitialspace=True)
#    Ycols=list()
#    for r in reader:
#        if (r[0]!=''): Ycols.append(int(r[0]))

    
#    for i in arange(Nfiles):
#        data=inFile[i]
#        colx=Xcols[i]
#        coly=Ycols[i]
#    plot(data[:,colx],data[:,coly])
    
#    print option.colx
    plotTypeIdx = 0
    global globalCurrentDatasetIdx
    global globalSubPlotGridSpec
    global globalPlotDs
    global globalAxes
    ax = None
    for i in arange(Nfiles):
        if option.saveData:
            np.savetxt("plot_function.dumpFile", inFile[i])
            
#         if option.transpose:
#             inFile[i]=inFile[i].T
            
        print "=================================================== PLOTTING DATASET %i ====================================================" % i
        plotTypeIdx = i % len(option.plotType)
        globalCurrentDatasetIdx = i
        figPlotNum = figRows * 100 + figCols * 10 + (i / option.dsPerPlot + 1)
        figIdx = int(i / option.dsPerPlot)
        newSubplot = True
        isLastDsInSubplot = False
            
        dsIdxInSubplot = i % option.dsPerPlot
        if i % option.dsPerPlot == 0:
            newSubplot = True
        else:
            newSubplot = False
            
        if i + 1 < Nfiles:
            if (i + 1) % option.dsPerPlot == 0:
                isLastDsInSubplot = True
            else:
                isLastDsInSubplot = False
        else:
            isLastDsInSubplot = True
            
            
        if globalPlotDs != None:
            figPlotNum, dsIdxInSubplot = calculatePlotID(i, figRows, figCols)
            figIdx = figPlotNum[2] - 1
            newSubplot = isNewPlot(i, figRows, figCols)
      
            if i + 1 < Nfiles:
                if isNewPlot(i + 1, figRows, figCols):
                    isLastDsInSubplot = True
                else:
                    isLastDsInSubplot = False

        global shareX
        if option.shareX and len(plotAxesList)>0 and figIdx>0:
            shareX=plotAxesList[0]


        print "figure plot number: ", figPlotNum
        print('shareX: ',shareX)
        print "first dataset in this subplot: ", newSubplot
        print "lsat dataset in this subplot: ", isLastDsInSubplot
        print "dataset idx in this subplot: ", dsIdxInSubplot
        if option.polar:
            ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
            print "Making polar axes"
        elif option.proj3d:
            ax = fig.add_subplot(111, projection='3d')
        else:
            if type(option.figCols) == type(list()) and type(option.figRows) == type(list()):
                if newSubplot:
#                     if ax==None:
                    ax = subplot(globalSubPlotGridSpec[figIdx], sharex=shareX)
#                     else:
#                         ax = subplot(globalSubPlotGridSpec[figIdx],sharex=ax)
            else:
                if type(figPlotNum) == type(list()):    
#                     if ax==None:
                    ax = subplot(figPlotNum[0], figPlotNum[1], figPlotNum[2], sharex=shareX)
#                     else:
#                         ax = subplot(figPlotNum[0],figPlotNum[1],figPlotNum[2],sharex=ax)
                else:
                    ax = subplot(figPlotNum, sharex=shareX)
                
        plotAxesList.append(ax)
        plotAxes = ax
        globalAxes = ax
        
        if len(inFile[i]) > 0:
            print "Plotting dataset %i of type: %s" % (i, option.plotType[plotTypeIdx])
            print "dataset has %i rows" % (len(inFile[i]))
            axes(ax)
            
            # select the right data and perform last-minute operations before plotting (only for plotTypes that are not maps)
            if option.plotType[plotTypeIdx] == 'map':
                data = inFile[i]
            if option.plotType[plotTypeIdx] == 'img':
                data = inFile[i]
    
            if (option.plotType[plotTypeIdx] != 'map' and option.plotType[plotTypeIdx] != 'img' and option.plotType[plotTypeIdx] != 'ts') or option.dateFmt == '%Y %m %d %H %M %S':
#             if (option.plotType[plotTypeIdx] != 'map' and option.plotType[plotTypeIdx] != 'img' and option.plotType[plotTypeIdx] != 'ts') or option.dateFmt == '':
                data = inFile[i]
                
                if option.plotType[plotTypeIdx] == 'barchartH':
                    barchartLabels = np.reshape(data[:, option.colx], (-1))
                    print barchartLabels
                    data = np.reshape(np.asarray(data[:, option.coly], dtype='float'), (-1))
#                     print data
                    data = np.vstack([np.arange(len(data)), data]).T
                    print data
                else:            
#                     data = array(inFile[i], dtype='float')
                    data=array(inFile[i])
#                 print data
#                 sys.exit()

#                 if option.big or len(option.rows)>0:
                if option.big > 0:
                    colx = 0
                    coly = 1
                    colColor = 2
                else:
                    colx = option.colx[i % len(option.colx)]
                    coly = option.coly[i % len(option.coly)]
                    colz = option.colz[i % len(option.colz)]
#                     colColor = option.colColor
                    colColor = option.colColor[i % len(option.colColor)]


                if option.logX:
                    data = data[data[:, colx] > 0]
                    if option.absX == False:
                        data = data[data[:, colx] > 0]
                if option.logY:
                    if option.absY == False:
                        data = data[data[:, coly] > 0]


                colxerr = ''
                if len(option.colxerr) > 0:
                    colxerr = option.colxerr[i % len(option.colxerr)]
                colxerr2 = ''
                if len(option.colxerr2) > 0:
                    colxerr2 = option.colxerr2[i % len(option.colxerr2)]
                
                colyerr = ''
                if len(option.colyerr) > 0:
                    colyerr = option.colyerr[i % len(option.colyerr)]
                colyerr2 = ''
                if len(option.colyerr2) > 0:
                    colyerr2 = option.colyerr2[i % len(option.colyerr2)]
                        
    #            print data[:,coly]
    
                if option.dateFmt != '':
                    dt_ncols = option.dateFmt.count(' ')
                    dt = list()
                    for sline in data:
#                         dt.append(' '.join([ str(int(float(tmp))) for tmp in sline[colx:colx + dt_ncols + 1] ]))
                        dt.append(sline[colx])
                    datex = list()
                    [ datex.append(datetime.datetime.strptime(str(tmpdate), option.dateFmt)) for tmpdate in dt ]
                    datax = datex
                else:
                    print "Using data columns: %i, %i" % (colx, coly)
                    print "Color column (if used): %i" % (colColor)
                    if colx == -1:
        #                if len(data.shape)==1:
        #                    datax=arange(len(data))
        #                else:
                        datax = arange(len(data[:, coly]))
                    else:
                        datax = data[:, colx]
        #                datax=array(data[:,colx], dtype='float')
                    
    #                if len(data.T)==1:
    #                    datay=data[:,0];
    #                else:
                datay = ''
                if option.plotType[plotTypeIdx] != 'hist':
                    datay = data[:, coly]
                if option.plotType[plotTypeIdx] == 'yerr' or option.plotType[plotTypeIdx] == 'fnshaded' or option.plotType[plotTypeIdx] == 'vect':
                    datayerr = data[:, colyerr]
                    if colyerr2 != '':
                        datayerr2 = data[:, colyerr2]
#            datay=array(data[:,coly], dtype='float')
                
                
                if option.plotType[plotTypeIdx] == 'barchartH':
                    print datax
                    print datay
                    ax.barh(datax, datay, align='center', color=option.fc[i % len(option.fc)], ecolor=option.ec[i % len(option.ec)])
                    ax.set_yticks(datax)
                    ax.set_yticklabels(barchartLabels)
                    ax.invert_yaxis()  # labels read top-to-bottom
                
                if option.plotType[plotTypeIdx] == 'hist':
                    datax = removeNans1(datax)
                    
                if option.logY:
                    for j in arange(len(datay)):
                        if datay[j] == 0:
                            datay[j] = MINIMAL_VALUE_FOR_LOGPLOT
    
                if option.absY:
                    for j in arange(len(datay)):
                        datay[j] = abs(datay[j])

                if option.absX:
                    for j in arange(len(datax)):
                        datax[j] = abs(datax[j])
                            
                if option.sqrt:
                    for j in arange(len(datay)):
                        datay[j] = sqrt(datay[j])
        
                if option.shift0Xmax:
                    datax = datax - datax[-1]
                if option.shift0X:
                    datax = datax - datax[0]
                if option.shift0:
                    datay=np.asanyarray(datay,dtype=float)
                    datay = datay - np.mean(datay)
                if option.polar:
                    print datax
                    datax = datax * np.pi / 180.0  # conversion to radians
                    print datax
                
                datax = performXdataOperations(datax, i)
                datay = performYdataOperations(datay, i)
        
        
                if option.colSize == -1:
                    sizeData = option.ps[i % len(option.ps)]
                    print "marker size to be used in next plot: %f" % sizeData
                else:
                    sizeData = data[:, option.colSize] * option.ps[i % len(option.ps)]
                    print "marker sizes to be used" 
                    print sizeData

                print "allocated marker sizes: "
                print option.ps
    
    
            ###############################################################################################
            ###############################################################################################
            ###############################################################################################
            #
            # make the various types of plots
            #
            ###############################################################################################
            ###############################################################################################
            ###############################################################################################
                globalCurrentDataset = data
#                print 'current global dataset'
#                print globalCurrentDataset
                globalCurrentDatasetX = datax
                globalCurrentDatasetY = datay
                if option.interactive:
                    global globalMaskedPointsMap
                    globalMaskedPointsMap = np.ones(len(globalCurrentDataset), dtype=int)
#                    print len(globalMaskedPointsMap)
#                    sys.exit()
        
            ###############################################################################################
            # img type plot
            ###############################################################################################
            if option.plotType[plotTypeIdx] == 'img':
    #            print ax
                axes(ax)
    #            sys.exit()
    #            img=Image.open(data)
                img = data
                orig = 'lower'
                if option.IMGflipY:
                    orig = 'upper'
                ax2 = fig.add_axes(ax)
                print img
                im = imshow(img, origin=orig, extent=(0, option.IMGextent[0], 0, option.IMGextent[1]), interpolation='bicubic', alpha=option.IMGalpha)
    #            hold(True)
    
            ###############################################################################################
            # map type plot
            ###############################################################################################
            
            if option.plotType[plotTypeIdx] == 'map':
    #            ax3 = fig.add_axes(ax)
                from mpl_toolkits.basemap import Basemap, shiftgrid
    
                axes(ax)
                global projectMap
                print data
                map = data[0]
                lons = data[1]
                lats = data[2]
#                 print 'map',map
                # this is a dirty hack to remove the data that cannot be plotted 
                # and if plotted the plot looks ugly.
                if option.proj == 'ortho' and option.lat0 == 90:
                    lats = lats[lats>0]
                    map=map[-len(lats):-1]
#                 print 'map',map
                
                
                map, lons = shiftgrid(180, map, lons, start=False)
                if option.proj == 'stere':
                    projectMap = Basemap(projection=option.proj, lon_0=0.5 * (lons[0] + lons[-1]), lat_0=option.lat0, lat_ts=option.lat0, llcrnrlat=option.lat1, urcrnrlat=option.lat2, llcrnrlon=option.lon1, urcrnrlon=option.lon2)
                else:
                    if option.sphereContour:
                        projectMap = Basemap(projection=option.proj, lon_0=0.5 * (lons[0] + lons[-1]))
        #                projectMap = Basemap(projection=option.proj,lon_0=option.lon0, lat_0=option.lat0)
                    else:
#                         projectMap = Basemap(projection=option.proj, lon_0=option.lon0, lat_0=option.lat0, llcrnrlat=option.lat1, urcrnrlat=option.lat2, llcrnrlon=option.lon1, urcrnrlon=option.lon2)
                        projectMap = Basemap(projection=option.proj, lon_0=option.lon0, lat_0=option.lat0)
                lons, lats = meshgrid(lons, lats)
                palette = getPalette()
                palette = cm.get_cmap(option.CM[figIdx % len(option.CM)], int(option.cbLabelsNum))

    #            palette = matplotlib.cm.jet
                palette.set_over(option.setAbove)  # , 1.0)
                palette.set_under(option.setBelow, 1)
                print 'map minimal value: %lE ' % np.amin(map.reshape(-1, 1))
                print 'map maximal value: %lE ' % np.amax(map.reshape(-1, 1))
                
    
                if option.vmin[figIdx % len(option.vmin)] == option.vmax[figIdx % len(option.vmin)]:
                    minv = np.amin(map.reshape(-1, 1))
                    maxv = np.amax(map.reshape(-1, 1))
                else:
                    minv = option.vmin[figIdx % len(option.vmin)]
                    maxv = option.vmax[figIdx % len(option.vmin)]
                    
                print "minv: " + str(minv)
                print "maxv: " + str(maxv)
                color_num = float(option.levels[i % len(option.levels)])
                delta = (maxv - minv) / color_num
                print 'delta: ', delta
                levels = arange(minv, maxv, delta)
                print 'levels: ', levels
                x, y = projectMap(lons, lats)
    
                
                
                # set bad values
#                 for j in arange(len(map)):
#                     for k in arange(len(map[j])):
#                         if map[j][k] == -1: 
#                             map[j][k] = nan
                if option.bad_val!=None:
                    palette.set_bad('#0d0d0d', float(option.bad_val))
                    
#                set masked values
#                 map = ma.array(map)
#                 for j in arange(len(map)):
#                     for k in arange(len(map[j])):
#                         if map[j][k]==option.bad_val: 
#                            map[j][k]=ma.masked
                
                
#                 map = np.ma.masked_invalid(map)
#                     map = np.ma.masked_where(map==option.bad_val,map)
                    map = np.ma.masked_values(map,option.bad_val)
                
#                 if option.set_below!=None:
#                     palette.set_under(option.set_below[1],option.set_below[0])
#                     levels[0]=option.set_below[0]
                
                # this works but is slow and leaves while lines along contours. the hack below doesn't seem to help
                if option.sphereContour:
                    cs = projectMap.contourf(x, y, map, levels, cmap=palette, norm=matplotlib.colors.normalize(vmin=levels[0], vmax=levels[-1], clip=True), extend='both', origin='upper')
                    #            for c in cs.collections:
                    #                c.set_linewidth( 0.0 )
                else:
                    # this is very fast and gives what I want
                    cs = pcolormesh(x, y, map, vmin=levels[0], vmax=levels[-1], norm=None, cmap=palette, picker=map_picker, edgecolor=None, lw=0)
    #            cs=projectMap.pcolormesh(x, y, mmap, vmin=levels[0], vmax=levels[-1], norm=None, cmap=palette, picker=map_picker, edgecolor=None, lw=0)
                # this works fast but doesn't project correctly
    #            cs=m1.imshow(map,cmap=palette)

                if option.plotLabelsFromFile != "":
                    plotLabelsFromFile(fig, ax, option.plotType[plotTypeIdx])
    
                axes(ax)

    
                # draw meridians and parallels
                if len(meridians) > 0:
#                     projectMap.drawmeridians(-meridians, labels=[0, 0, 0, 1], color=option.MPcolor)
    #                projectMap.drawmeridians(arange(180,180.2,0.1),labels=[0,0,0,1], color=option.MPcolor)
#                    if option.proj=='lcc':
#                        txtl=[]; txtb=[]; txt=[];
#                        for j in range(len(meridians)):
#                            txtl.append(meridians[j]);
#                            txtb.append(-5);
#                            txt.append('%1.1f' % (meridians[j]))
#    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
#                            xpt,ypt = projectMap(-txtl[j],txtb[j])
#                            text(xpt,ypt,txt[j],fontsize=option.MPfontSize)

                    if option.proj == 'moll':
                        projectMap.drawmeridians(-meridians, labels=[0, 0, 0, 1], color=option.MPcolor)
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(-5);
                            txt.append('%1.0f' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                    elif option.proj == 'ortho' and option.lat0 == 90:
                        projectMap.drawmeridians(-meridians, labels=[0, 0, 0, 0], color=option.MPcolor)
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(1);
                            txt.append(' %1.0f ' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            if txtl[j] > 90 and txtl[j] < 270:
                                if txtl[j] > 0 and txtl[j] <= 180:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, color=option.MPcolor, horizontalalignment='right', verticalalignment='bottom')
                                else:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, color=option.MPcolor, horizontalalignment='left', verticalalignment='bottom')
                            else:
                                if txtl[j] > 0 and txtl[j] <= 180:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, color=option.MPcolor, horizontalalignment='right', verticalalignment='top')
                                else:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, color=option.MPcolor, horizontalalignment='left', verticalalignment='top')
                    
                    else:
                        projectMap.drawmeridians(-meridians, labels=[0, 0, 0, 1], color=option.MPcolor)
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(-5);
                            if meridianLabels != None:
                                mltmp = meridianLabels[j].split(' ')
                                if len(mltmp) == 1: 
                                    mltmp.append('')
                                txt.append(u'%s\N{DEGREE SIGN}%s' % (mltmp[0], mltmp[1]))
                            else:
                                txt.append('%1.0f' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                        
                                
    
                if len(parallels) > 0:
                    if option.proj == 'moll':
                        projectMap.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=option.MPfontSize, color=option.MPcolor)
                    if (option.proj == 'ortho' or option.proj == 'stere') and option.lat0 == 90:
                        txtl = []; txtb = []; txt = [];
                        print 'parallels:'
                        print parallels
                        for j in range(len(parallels)):
                            txtl.append(0);
                            txtb.append(parallels[j]);
                            txt.append('%1.0f' % (parallels[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], color=option.MPcolor, fontsize=option.MPfontSize)
                        if len(parallels) > 0:
                            projectMap.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=option.MPfontSize, color=option.MPcolor)
#                     elif option.proj=='lcc' or option.proj=='stere':
#                         projectMap.drawparallels(parallels,labels=[True,True,True,True])
                    elif option.proj == 'stere':
                        projectMap.drawparallels(parallels, labels=[True, True, True, True])
                    else:
                        txtl = []; txtb = []; txt = [];
                        print 'parallels:'
                        print parallels
                        for j in range(len(parallels)):
                            txtl.append(0);
                            txtb.append(parallels[j]);
                            txt.append('%1.0f' % (parallels[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                        if len(parallels) > 0:
                            projectMap.drawparallels(parallels, labels=[1, 1, 1, 1], fontsize=option.MPfontSize, color=option.MPcolor)
                        
#                # draw meridians and parallels
#                if len(meridians)>0:
#                    projectMap.drawmeridians(-meridians,labels=[0,0,0,1], color=option.MPcolor)
#                    txtl=[]; txtb=[]; txt=[];
#                    for j in range(len(meridians)):
#                        txtl.append(meridians[j]);
#                        txtb.append(-5);
#                        txt.append('%1.1f' % (meridians[j]))
#                        xpt,ypt = projectMap(-txtl[j],txtb[j])
#                        text(xpt,ypt,txt[j],fontsize=option.MPfontSize)
#    
#                if len(parallels)>0:
#                    projectMap.drawparallels(parallels,labels=[1,0,0,0],fontsize=option.MPfontSize, color=option.MPcolor)
    
                projectMap.drawmapboundary(color='k', linewidth=1.0, ax=None)
                if option.colorbar:
                    cb = colorbar(cs)
#                     cb.set_label(option.zlabel)
#                     cb.set_label(option.zlabel[figIdx % len(option.zlabel)])
                    cb.set_label(option.zlabel[figIdx % len(option.zlabel)], size=option.fontSizeCM[figIdx % len(option.fontSizeCM)])

                    if option.fontSizeCM[figIdx % len(option.fontSizeCM)] > 0:
                        cb.ax.tick_params(labelsize=option.fontSizeCM[figIdx % len(option.fontSizeCM)])

                if option.interactive:
                    fig.canvas.mpl_connect('pick_event', onMapPick)


    
                
            ###############################################################################################
            # sphere type plot
            ###############################################################################################
            axes(ax)
            
            if option.plotType[plotTypeIdx] == 'sphere':
                from mpl_toolkits.basemap import Basemap, shiftgrid
                global projectMap
                if len(option.label) > 0:
                    plotLegendLabel = option.label[i % len(option.label)].decode('utf8')
                else:
                    plotLegendLabel = None
                
                lons = -datax
                lats = datay
#                 map = data[:, option.colColor]
                map = data[:, option.colColor[i % len(option.colColor)]]
    
                if projectMap == None:
                    if option.proj == 'moll' or option.proj == 'ortho':
                        projectMap = Basemap(projection=option.proj, lon_0=option.lon0, lat_0=option.lat0)
                    elif option.proj == 'stere':
                        projectMap = Basemap(projection=option.proj, lon_0=option.lon0, lat_0=option.lat0, lat_ts=option.lat0, llcrnrlat=option.lat1, urcrnrlat=option.lat2, llcrnrlon=option.lon1, urcrnrlon=option.lon2, resolution=option.mapresolution)
                    
                    elif option.proj == 'lcc':
                            projectMap = Basemap(projection=option.proj, lon_0=option.lon0, lat_0=option.lat0, lat_1=option.lat1, lat_2=option.lat2, resolution=option.mapresolution, width=option.mapwidth * 1000, height=option.mapheight * 1000)


    #                projectMap = Basemap(projection=option.proj,lon_0=0.5*(lons[0]+lons[-1]), lat_0=option.lat0)
                if option.reverseLon:
                    x, y = projectMap(lons, lats)
                else:
                    x, y = projectMap(-lons, lats)
                if option.pt[i % len(option.pt)] != 'o':
                    print '\n\nThe selected marker type for this plot is not suitable - ".", will change it to "o"\n\n'
                    option.pt[i % len(option.pt)] = 'o'
    
    
                palette = getPalette()
                palette.set_over(option.setAbove)  # , 1.0)
                palette.set_under(option.setBelow, 1)

                if option.vmin[figIdx % len(option.vmin)] == option.vmax[figIdx % len(option.vmin)]:
                    minv = np.amin(map.reshape(-1, 1))
                    maxv = np.amax(map.reshape(-1, 1))
                else:
                    minv = option.vmin[figIdx % len(option.vmin)]
                    maxv = option.vmax[figIdx % len(option.vmin)]
                    
                print "minv: " + str(minv)
                print "maxv: " + str(maxv)
                color_num = float(option.levels[i % len(option.levels)])
                delta = (maxv - minv) / color_num
                print 'delta: ', delta
                levels = arange(minv, maxv, delta)
                print 'levels: ', levels
                 
#                 if option.colColor == -1:
                if option.colColor[i % len(option.colColor)] == -1:
    #                global scatterGlobalData
    #                scatterGlobalData=option.pc[i % len(option.pc) ]
                    print 'will use points color: ', option.pc[i % len(option.pc) ]
                    scatter(x, y, c=option.pc[i % len(option.pc) ], s=sizeData, lw=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], marker=option.pt[i % len(option.pt)], label=plotLegendLabel)  # , marker=option.pt[i % len(option.pt)]
                else:

    #                global scatterGlobalData
    #                scatterGlobalData=map
                    scatter(x, y, c=map, s=sizeData, cmap=palette, norm=matplotlib.colors.normalize(vmin=levels[0], vmax=levels[-1], clip=True), lw=option.markerEdgeWidth[i % len(option.markerEdgeWidth)] , marker=option.pt[i % len(option.pt)], label=plotLegendLabel)
    #                scatter(x,y, c=map, s=40, lw=option.markerEdgeWidth , marker=option.pt[i % len(option.pt)])
    
                if option.plotLabelsFromFile != "":
                    plotLabelsFromFile(fig, ax, option.plotType[plotTypeIdx])
    
                axes(ax)

    #            circ=Circle((x[0],y[0]),radius=x[0]/5)
    #            ax.add_patch(circ)
    
                # draw meridians and parallels
                if len(meridians) > 0:
                    projectMap.drawmeridians(-meridians, labels=[0, 0, 0, 1], fontsize=option.MPfontSize, color=option.MPcolor)
    #                projectMap.drawmeridians(arange(180,180.2,0.1),labels=[0,0,0,1], color=option.MPcolor)
#                    if option.proj=='lcc':
#                        txtl=[]; txtb=[]; txt=[];
#                        for j in range(len(meridians)):
#                            txtl.append(meridians[j]);
#                            txtb.append(-5);
#                            txt.append('%1.1f' % (meridians[j]))
#    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
#                            xpt,ypt = projectMap(-txtl[j],txtb[j])
#                            text(xpt,ypt,txt[j],fontsize=option.MPfontSize)

                    if option.proj == 'moll':
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(-5);
                            if meridianLabels != None:
                                mltmp = meridianLabels[j].split(' ')
                                if len(mltmp) == 1: 
                                    mltmp.append('')
                                txt.append(u'%s\N{DEGREE SIGN}%s' % (mltmp[0], mltmp[1]))
                            else:
                                txt.append('%1.1f' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                    elif option.proj == 'ortho' and option.lat0 == 90:
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(1);
                            txt.append(' %1.1f ' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            if txtl[j] > 90 and txtl[j] < 270:
                                if txtl[j] > 0 and txtl[j] <= 180:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, horizontalalignment='right', verticalalignment='bottom')
                                else:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, horizontalalignment='left', verticalalignment='bottom')
                            else:
                                if txtl[j] > 0 and txtl[j] <= 180:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, horizontalalignment='right', verticalalignment='top')
                                else:
                                    text(xpt, ypt, txt[j], fontsize=option.MPfontSize, horizontalalignment='left', verticalalignment='top')
                    
                    else:
                        txtl = []; txtb = []; txt = [];
                        for j in range(len(meridians)):
                            txtl.append(meridians[j]);
                            txtb.append(-5);
                            if meridianLabels != None:
                                mltmp = meridianLabels[j].split(' ')
                                if len(mltmp) == 1: 
                                    mltmp.append('')
                                txt.append(u'%s\N{DEGREE SIGN}%s' % (mltmp[0], mltmp[1]))
                            else:
                                txt.append('%1.1f' % (meridians[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                                    
    
                if len(parallels) > 0:
                    if option.proj == 'moll':
                        projectMap.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=option.MPfontSize, color=option.MPcolor)
                    if (option.proj == 'ortho' or option.proj == 'stere') and option.lat0 == 90:
                        txtl = []; txtb = []; txt = [];
                        print 'parallels:'
                        print parallels
                        for j in range(len(parallels)):
                            txtl.append(0);
                            txtb.append(parallels[j]);
                            txt.append('%1.1f' % (parallels[j]))
    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
                            xpt, ypt = projectMap(-txtl[j], txtb[j])
                            text(xpt, ypt, txt[j], fontsize=option.MPfontSize)
                        if len(parallels) > 0:
                            projectMap.drawparallels(parallels, labels=[1, 0, 0, 0], fontsize=option.MPfontSize, color=option.MPcolor)
                    if option.proj == 'lcc':
                        projectMap.drawparallels(parallels, labels=[True, True, True, True])
#                        txtl=[]; txtb=[]; txt=[];
#                        print 'parallels:'
#                        print parallels
#                        for j in range(len(parallels)):
#                            txtl.append(0);
#                            txtb.append(parallels[j]);
#                            txt.append('%1.1f' % (parallels[j]))
#    #                        xpt,ypt = projectMap(txtl[j]+0.01,txtb[j])
#                            xpt,ypt = projectMap(-txtl[j],txtb[j])
#                            text(xpt,ypt,txt[j],fontsize=option.MPfontSize)
    
                if option.colorbar and option.CM[figIdx % len(option.CM)]!='None':
                    cb = colorbar()
#                     cb.set_label(option.zlabel)
                    cb.set_label(option.zlabel[figIdx % len(option.zlabel)], size=option.fontSizeCM[figIdx % len(option.fontSizeCM)])


                axes(ax)
                if option.xlabel[i % len(option.xlabel)] != "None":
                    xlabel(option.xlabel[i % len(option.xlabel)], fontsize=option.fontSizeLabels[i % len(option.fontSizeLabels)])
                if option.ylabel[i % len(option.ylabel)] != "None":
                    ylabel(option.ylabel[i % len(option.ylabel)], fontsize=option.fontSizeLabels[i % len(option.fontSizeLabels)])
#                if option.fontSize[i % len(option.fontSize)]>0:
#                 if option.fontSizeXtickLabels[i % len(option.fontSizeXtickLabels)]>0.0:
#                     setp( ax.get_yticklabels(), fontsize=option.fontSizeXtickLabels[i % len(option.fontSizeXtickLabels)])
#                 else:
#                     ax.get_xaxis().set_visible(False)
#                 if option.fontSizeYtickLabels[i % len(option.fontSizeYtickLabels)]>0.0:
#                     setp( ax.get_xticklabels(), fontsize=option.fontSizeYtickLabels[i % len(option.fontSizeYtickLabels)])
#                 else:
#                     ax.get_yaxis().set_visible(False)

                projectMap.drawmapboundary(color='k', linewidth=2.0, ax=None)
                if option.coastLines:
                    projectMap.drawcoastlines()
                if option.countryLines:
                    projectMap.drawcountries()
                if option.riverLines:
                    projectMap.drawrivers()
    
                if len(option.label) > 0 and i == Nfiles - 1:  # this isn't working as is should - FIX IT
                    legend(option.label.decode('utf8'), loc=legendLocation[figIdx % len(legendLocation)], prop={'size':option.fontSizeLegend[figIdx % len(option.fontSizeLegend)]}, ncol=option.legendNcol, mode=option.legendMode, borderaxespad=option.legendBorderAxesPad)



            if option.plotType[plotTypeIdx] == 'sphere' or option.plotType[plotTypeIdx] == 'map':

                #
                # circles            
                #
                if len(circlePatchDef) > 0:
                    axes(ax)
                    plotCircles(ax, circlePatchDef)



                #
                # title and labels
                #
                if isLastDsInSubplot:
                    if option.title[figIdx % len(option.title)] != "":
                        tit = option.title[figIdx % len(option.title)]
                        title(tit, fontsize=option.fontSizeTitle, horizontalalignment=option.titleAlignH, verticalalignment=option.titleAlignV)


                #
                # extra texts
                #
#                 tstrIdx=0
                if isLastDsInSubplot:
                    if len(option.textxy) > 0:
                        tstr = option.textxy[figIdx]
    #                    if tstrIdx==figIdx and i % option.dsPerPlot==0:
    #                     if tstrIdx==dsIdxInSubplot and newSubplot:
                            # do the first split for all texts that we want on this figure (if any)
                        if tstr != 'None':
                            xyTextTuple = tstr.split(';')
                            print 'xyTextTuple', xyTextTuple
                            for tstr in xyTextTuple:
                                tstrSplit = re.split(r'(?<!\\),', tstr)  # negative lookbehind assertion
                                if len(tstrSplit) >= 3:
                                    tx = float(tstrSplit[0])
                                    ty = float(tstrSplit[1])
                                    t = tstrSplit[2].replace('\,', ',')
                                    trot = 0
                                    tsiz = option.extraTextFontSize
                                    tcolor = 'k'
                                if len(tstrSplit) >= 4:
                                    trot = tstrSplit[3]
                                if len(tstrSplit) >= 5:
                                    tsiz = tstrSplit[4]
                                if len(tstrSplit) >= 6:
                                    tcolor = tstrSplit[5]
                                    
                                text(tx, ty, t, fontsize=tsiz, zorder=10000, rotation=trot, color=tcolor)
                        
#                     tstrIdx+=1
                    if len(option.textfxfy) > 0:
                        tstr = option.textfxfy[figIdx]
                        if tstr != 'None':
                            xyTextTuple = tstr.split(';')
                            print 'xyTextTuple', xyTextTuple
                            curXlimits = ax.get_xlim()
                            curYlimits = ax.get_ylim()
                            curXspan = curXlimits[1] - curXlimits[0]
                            curYspan = curYlimits[1] - curYlimits[0]
                            for tstr in xyTextTuple:
                                tstrSplit = re.split(r'(?<!\\),', tstr)  # negative lookbehind assertion
                                if len(tstrSplit) >= 3:
                                    tx = float(tstrSplit[0]) * (curXspan) + curXlimits[0]
                                    ty = float(tstrSplit[1]) * (curYspan) + curYlimits[0]
                                    print '*************************************'
                                    print '*************************************'
                                    print '*************************************'
                                    print '*************************************'
                                    print tx, ty
                                    t = tstrSplit[2].replace('\,', ',')
                                    trot = 0
                                    tsiz = option.extraTextFontSize
                                    tcolor = 'k'
                                if len(tstrSplit) >= 4:
                                    trot = tstrSplit[3]
                                if len(tstrSplit) >= 5:
                                    tsiz = tstrSplit[4]
                                if len(tstrSplit) >= 6:
                                    tcolor = tstrSplit[5]
                                    
                                text(tx, ty, t, fontsize=tsiz, zorder=10000, rotation=trot, color=tcolor)


            ###############################################################################################
            # time sequence type plot
            ###############################################################################################
            if option.plotType[plotTypeIdx] == 'ts':
                if len(option.label) > 0:
                    if option.label[i % len(option.label)] != "None":
                        plotLegendLabel = option.label[i % len(option.label)].decode('utf8')
                    else:
                        plotLegendLabel = None
                else:
                    plotLegendLabel = None
                
                datax = inFile[i][:, 0]
                datay = inFile[i][:, 1]
                import matplotlib.dates as mdates
                import time
                datex = list()
#                [ datex.append(datetime.datetime.strptime(tmpdate, '%M/%d/%Y')) for tmpdate in datax ]
                print(datax)
                try:
                    [ datex.append(datetime.datetime.strptime(str(tmpdate), option.dateFmt)) for tmpdate in datax ]
                except ValueError,msg:
                    print(tmpdate)
                    print(msg)
                    raise(ValueError)
#                print 'date X'
#                print datex
                print len(datex)
                datex = performXdataOperations(datex, i)
#                print 'date X after operations'
#                print datex
#                print len(datex)
#                print 'data Y'

                datay = np.asanyarray(datay, dtype='float')
                if option.shift0:
                    datay = datay - mean(datay)

                datay = performYdataOperations(datay, i)

#                print datay
#                plot(datex, datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)],mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)
#                plot(datex, datay, lw='1', marker='.')
                plot(datex, datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)], mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)

                seconds = mdates.SecondLocator()
                minutes = mdates.MinuteLocator()
                hours = mdates.HourLocator()
                days = mdates.DayLocator()
                years = mdates.YearLocator()  # every year
                months = mdates.MonthLocator()  # every month
#                yearsFmt = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
                yearsFmt = mdates.DateFormatter(option.dateFmtPlot)
                ax.xaxis.set_major_formatter(yearsFmt)
                # format the ticks
                ax.xaxis.set_major_locator(days)
                ax.xaxis.set_major_formatter(yearsFmt)
#                 ax.xaxis.set_minor_locator(hours)
#                 ax.xaxis.set_minor_locator(None)

                if option.years:
                    ax.xaxis.set_major_locator(years)
                    if option.minorTicks:
                        ax.xaxis.set_minor_locator(months)
                if option.months:
                    ax.xaxis.set_major_locator(months)
                    if option.minorTicks:
                        ax.xaxis.set_minor_locator(days)
                if option.days:
                    ax.xaxis.set_major_locator(days)
                    if option.minorTicks:
                        ax.xaxis.set_minor_locator(hours)
                if option.hours:
                    ax.xaxis.set_major_locator(hours)
                    if option.minorTicks:
                        ax.xaxis.set_minor_locator(minutes)
                if option.minutes:
                    ax.xaxis.set_major_locator(minutes)
                    if option.minorTicks:
                        ax.xaxis.set_minor_locator(seconds)

#                datemin = datetime.date(2013, 1, 1)
#                datemax = datetime.date(2015, 1, 1)
                datemin = datex[0]
                datemax = datex[-1]
                print datemin
                print datemax
                print option.datexmin[figIdx % len(option.datexmin)]
                print option.datexmax
                if option.datexmin[figIdx % len(option.datexmin)] != -1:
                    datemin = datetime.datetime.strptime(option.datexmin[figIdx % len(option.datexmin)], option.dateFmt)
                if option.datexmax[figIdx % len(option.datexmax)] != -1:
                    datemax = datetime.datetime.strptime(option.datexmax[figIdx % len(option.datexmax)], option.dateFmt)
                print 'datemin,datemax: ', datemin, datemax    
                ax.set_xlim(datemin, datemax)
#                ax.format_xdata = mdates.DateFormatter('%Y/%m/%d %H:%M:%S')
#                fig.autofmt_xdate()

                ymin, ymax = ax.get_ylim()
                if option.ymin[figIdx % len(option.ymin)] != None:
                    ymin = option.ymin[figIdx % len(option.ymin)]
                if option.ymax[figIdx % len(option.ymax)] != None:
                    ymax = option.ymax[figIdx % len(option.ymax)]
                print 'setting y limits: %f, %f' % (ymin, ymax)
#                ylim([ymin,ymax])
                ax.set_ylim(ymin, ymax)
                print
                print 'Setting xticklabels rotation: ', option.Rxlabels    
                print
                setp(ax.get_xticklabels(), rotation=option.Rxlabels, horizontalalignment='center', verticalalignment='top') 
                setp(ax.get_yticklabels(), rotation=option.Rylabels) 
                xticks(rotation=option.Rxlabels)
#                if option.xticks[figIdx % len(option.xticks)]!=0:
#                    tcs=ax.get_xticks()
#                    print tcs
#                    xticks(arange(min(tcs[0],option.xmin[figIdx % len(option.xmin)]),max(tcs[-1],option.xmax[figIdx % len(option.xmax)]),option.xticks[figIdx % len(option.xticks)]))
                
            ###############################################################################################
            # function type plot
            ###############################################################################################
    
            if option.plotType[plotTypeIdx] == 'vect' or option.plotType[plotTypeIdx] == 'ts' or option.plotType[plotTypeIdx] == 'fn' or option.plotType[plotTypeIdx] == 'step' or option.plotType[plotTypeIdx] == 'err' or option.plotType[plotTypeIdx] == 'fnshaded' or option.plotType[plotTypeIdx] == 'fillbtwX' or option.plotType[plotTypeIdx] == 'fillbtwY' or option.plotType[plotTypeIdx] == 'scat' or option.plotType[plotTypeIdx] == 'scat3d' or option.plotType[plotTypeIdx] == 'scatContFill' or option.plotType[plotTypeIdx] == 'scatContFillD' or option.plotType[plotTypeIdx] == 'scatDensContFill' or option.plotType[plotTypeIdx] == 'scatDensCont' or option.plotType[plotTypeIdx] == 'hist' or option.plotType[plotTypeIdx] == 'circ' or option.plotType[plotTypeIdx] == 'barchartH':
    #            if option.colColor!=-1 or option.colSize!=-1:
                if len(option.label) > 0:
                    if option.label[i % len(option.label)] != "None":
                        plotLegendLabel = option.label[i % len(option.label)].decode('utf8')
                    else:
                        plotLegendLabel = None
                else:
                    plotLegendLabel = None

                if option.plotType[plotTypeIdx] == 'scat'  or option.plotType[plotTypeIdx] == 'scat3d' or option.plotType[plotTypeIdx] == 'scatContFillD' or option.plotType[plotTypeIdx] == 'scatContFill' or option.plotType[plotTypeIdx] == 'scatDensContFill' or option.plotType[plotTypeIdx] == 'scatDensCont':  # or option.plotType[plotTypeIdx]=='scatcirc':
                    global scatterGlobalData
#                     if option.colColor == -1:
                    if option.colColor[i % len(option.colColor)] == -1:
                        data = arange(len(data))
                        scatterGlobalData = data
                    else:
#                        scatterGlobalData=data[:,option.colColor]
                        scatterGlobalData = data[:, colColor]
                        if option.logZ:
                            scatterGlobalData = np.log10(scatterGlobalData)
                    print "marker: %s" % option.pt[i % len(option.pt)]
                    print "marker size: %s" % option.ps[i % len(option.ps)]

                    if option.dateFmt == '':
                        datax, datay, scatterGlobalData = removeNans3(datax, datay, scatterGlobalData)
#                     palette = getPalette()

                    palette=None
                    if option.CM[figIdx % len(option.CM)]!='None':
                        palette = cm.get_cmap(option.CM[figIdx % len(option.CM)], int(option.cbLabelsNum))
                        if option.CMrange!='':
                            from matplotlib.colors import LinearSegmentedColormap
                            CMrange=cpedsPythCommon.getFloatList(option.CMrange)
                            # Remove the middle 40% of the RdBu_r colormap
    #                         interval = np.hstack([np.linspace(0, CMclip[0]), np.linspace(CMclip[1], 1)])
                            # remove top white part of the hot colormap
                            interval = np.linspace(CMrange[0], CMrange[1])
    #                         colors = plt.cm.hot(interval)
                            print 'palette: ',palette
                            colors=palette(interval)
                            cmap = LinearSegmentedColormap.from_list('name', colors)                        
                            palette = cm.get_cmap(cmap, int(option.cbLabelsNum))
                        
#                     palette = cm.get_cmap(option.CM, int(option.levels))

#                    scatter(datax,datay, c=scatterGlobalData, s=sizeData, lw=1, marker=option.pt[i % len(option.pt)], label=plotLegendLabel, cmap=palette,  zorder=option.zorder[i % len(option.zorder)], picker=True)
#                    print option.pt[i % len(option.pt)]
#                    cpedsPythCommon.waitEnter()
#                    if len(option.levels[i % len(option.levels)])==1:
#                        option.levels[i % len(option.levels)]=np.arange()
                    if option.pt[i % len(option.pt)] == '.' or option.pt[i % len(option.pt)] == None:
                        print 'changing the marker for the scatter type plot from "." to "o"'
                        option.pt[i % len(option.pt)] = "o"
#                        cpedsPythCommon.waitEnter()

                    if option.plotType[plotTypeIdx] == 'scat':
                        if option.colColor[i % len(option.colColor)] == -1:
                            scatter(datax, datay, c=option.pc[i % len(option.pc) ], s=sizeData, lw=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], marker=option.pt[i % len(option.pt)], label=plotLegendLabel, cmap=palette, zorder=option.zorder[i % len(option.zorder)], picker=True)
                        else:
                            cmax = max(scatterGlobalData)
                            cmin = min(scatterGlobalData)
                            if option.vmin[figIdx % len(option.vmin)] != option.vmax[figIdx % len(option.vmin)]:
#                             if option.vmin[figIdx % len(option.vmin)] != None and option.vmax[figIdx % len(option.vmin)] != None:
                                cmin = option.vmin[figIdx % len(option.vmin)]
                                cmax = option.vmax[figIdx % len(option.vmin)]

                            scatter(datax, datay, c=scatterGlobalData, s=sizeData, vmin=cmin, vmax=cmax, lw=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], marker=option.pt[i % len(option.pt)], label=plotLegendLabel, cmap=palette, zorder=option.zorder[i % len(option.zorder)], picker=True)
                    elif option.plotType[plotTypeIdx] == 'scat3d':
                            dataz = inFile[i][:,colz]
                            scatter(datax, datay, dataz, c=scatterGlobalData, marker=option.pt[i % len(option.pt)], zorder=option.zorder[i % len(option.zorder)])
                        
                    elif option.plotType[plotTypeIdx] == 'scatContFillD':
                        triang = Triangulation(datax, datay)
                        
                        # Mask off unwanted triangles.
                        xmid = datax[triang.triangles].mean(axis=1)
                        ymid = datay[triang.triangles].mean(axis=1)
#                        mask = np.where(xmid*xmid + ymid*ymid < 1e-3*1e-3, 1, 0)
                        mask = np.where(xmid * xmid + ymid * ymid < 1e-3 * 1e-3, 1, 0)
                        triang.set_mask(mask)
                        plt.triplot(triang, lw=0.5, color='grey')
                        # refine
                        refine = False
                        if refine:
                            refiner = UniformTriRefiner(triang)
                            tri_refi, z_test_refi = refiner.refine_field(scatterGlobalData, subdiv=3)
                            if len(option.levels[i % len(option.levels)]) == 1:
                                tricontourf(tri_refi, z_test_refi, option.levels[i % len(option.levels)][0], cmap=palette)
                            else:
                                tricontourf(tri_refi, z_test_refi, levels=option.levels[i % len(option.levels)], cmap=palette)
                        else:
                            if len(option.levels[i % len(option.levels)]) == 1:
                                tricontourf(triang, scatterGlobalData, option.levels[i % len(option.levels)][0], cmap=palette)
                            else:
                                tricontourf(triang, scatterGlobalData, levels=option.levels[i % len(option.levels)], cmap=palette)

                    elif option.plotType[plotTypeIdx] == 'scatContFill':
                        triang = Triangulation(datax, datay)
                        # make a regular grid for data interpolation
                        xi, yi = mkGrid(datax, datay, gridNbins)
                        xA, yA = np.meshgrid(xi, yi)
                        interp_lin = mtri.LinearTriInterpolator(triang, scatterGlobalData)
                        zA_lin = interp_lin(xA, yA)                        
                        
                        if len(option.levels[i % len(option.levels)]) == 1:
                            print option.levels[i % len(option.levels)]
                            contourf(xA, yA, zA_lin, int(option.levels[i % len(option.levels)][0]), cmap=palette, zorder=option.zorder[i % len(option.zorder)])
                        else:
#                             print option.levels[i % len(option.levels)]
#                             contourf(triang, zA_lin, levels=option.levels[i % len(option.levels)], cmap=palette, zorder=option.zorder[i % len(option.zorder)])
                            contourf(xA, yA, zA_lin, levels=option.levels[i % len(option.levels)], cmap=palette, zorder=option.zorder[i % len(option.zorder)])
                            
                    elif option.plotType[plotTypeIdx] == 'scatDensContFill':
                        datax, datay, H = calculatePointsDensityDistribution(datax, datay, gridNbins)
                        if len(option.levels[i % len(option.levels)]) == 1:
                            contourf(datax, datay, H, option.levels[i % len(option.levels)][0], cmap=palette, zorder=option.zorder[i % len(option.zorder)])
                            contour(datax, datay, H, option.levels[i % len(option.levels)][0], lw=option.width[i % len(option.width)], linestyle=option.ls[i % len(option.ls)], colors=option.lc[i % len(option.lc)], zorder=option.zorder[i % len(option.zorder)])
                        else:
                            contourf(datax, datay, H, levels=option.levels[i % len(option.levels)], cmap=palette, zorder=option.zorder[i % len(option.zorder)])
                            contour(datax, datay, H, levels=option.levels[i % len(option.levels)], lw=option.width[i % len(option.width)], linestyle=option.ls[i % len(option.ls)], colors=option.lc[i % len(option.lc)], zorder=option.zorder[i % len(option.zorder)])
                            
                    elif option.plotType[plotTypeIdx] == 'scatDensCont':
                        datax, datay, H = calculatePointsDensityDistribution(datax, datay, gridNbins)
                        CS = contour(datax, datay, H, levels=option.levels[i % len(option.levels)], lw=option.width[i % len(option.width)], linestyle=option.ls[i % len(option.ls)], colors=option.lc[i % len(option.lc)], zorder=option.zorder[i % len(option.zorder)])
                        clabel(CS, inline=1, fontsize=10, fmt=option.ZticksFmt)
                        
                        
#                        print H
#                    else:
#                        for idx in np.arange(len(datax)):
#                            print "generating circle: ",idx
#                            circ = Circle((datax[idx],datay[idx]), sizeData[idx], fc=None, ec='k', color=None, fill=False)
#                            ax.add_patch(circ)
                
                    
                        
                elif option.plotType[plotTypeIdx] == 'fn':
                    print
                    print "plot %i" % i
                    print "marker %s" % option.pt[i % len(option.pt)]
                    print "color %s" % option.pc[i % len(option.pc)]
                    print "all markers" 
                    print option.pt
                    print "all colors" 
                    print option.pc
                    print
                    if option.dateFmt == '':
                        datax, datay = removeNans2(datax, datay)

                    
                    plotLine = plot(datax, datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)], mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)


                    if option.ls[i % len(option.ls)] == '--' or option.ls[i % len(option.ls)] == '-.' or option.ls[i % len(option.ls)] == ':':
                        if type(option.lsdash) == type(list()):
                            plotLine[-1].set_dashes(option.lsdash[i % len(option.lsdash)])

                elif option.plotType[plotTypeIdx] == 'step':
                    print
                    print "plot %i" % i
                    print "marker %s" % option.pt[i % len(option.pt)]
                    print "color %s" % option.pc[i % len(option.pc)]
                    print "all markers" 
                    print option.pt
                    print "all colors" 
                    print option.pc
                    print
                    if option.dateFmt == '':
                        datax, datay = removeNans2(datax, datay)

                    
                    plotLine = step(datax, datay, lw=option.width[i % len(option.width)], where='mid', marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)], mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)


                    if option.ls[i % len(option.ls)] == '--' or option.ls[i % len(option.ls)] == '-.' or option.ls[i % len(option.ls)] == ':':
                        if type(option.lsdash) == type(list()):
                            plotLine[-1].set_dashes(option.lsdash[i % len(option.lsdash)])

                elif option.plotType[plotTypeIdx] == 'vect':
                    global scatterGlobalData
                    if option.colColor[i % len(option.colColor)] == -1:
                        data = arange(len(data))
                        scatterGlobalData = data
                    else:
#                        scatterGlobalData=data[:,option.colColor]
                        scatterGlobalData = data[:, colColor]
                        if option.logZ:
                            scatterGlobalData = np.log10(scatterGlobalData)
                    print "marker: %s" % option.pt[i % len(option.pt)]
                    print "marker size: %s" % option.ps[i % len(option.ps)]

                    palette = getPalette()
                    
                    
                    dataxerr = None
                    datayerr = None
                    print colxerr, colyerr
                    print data
                    if colyerr != '':
                        datayerr = data[:, colyerr]
                    if colyerr2 != '':
                        datayerr2 = data[:, colyerr2].T
#                        datayerr=np.vstack([datayerr,datayerr2]).T
                        datayerr = data[:, [colyerr, colyerr2]].T
#                        datayerr=np.concatenate([datayerr,datayerr2],axis=1)
                    if colxerr != '':
                        dataxerr = data[:, colxerr]
                    if colxerr2 != '':
                        dataxerr = data[:, [colxerr, colxerr2]].T
                    print dataxerr
                    dataxerr = dataxerr * 100
                    print dataxerr
                    datayerr = datayerr * 100
                    if option.colColor[i % len(option.colColor)] == -1:
                        Q = quiver(datax, datay, dataxerr, datayerr, scale=option.Qscale, lw=option.width[i % len(option.width)], linestyle=option.ls[i % len(option.ls)], color=option.lc[i % len(option.lc)])
                    else:
                        Q = quiver(datax, datay, dataxerr, datayerr, scatterGlobalData, scale=option.Qscale, lw=option.width[i % len(option.width)], linestyle=option.ls[i % len(option.ls)], cmap=palette)
#                    qk = quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})
#                        qk = quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W', fontproperties={'weight': 'bold'})

                elif option.plotType[plotTypeIdx] == 'err':
                    dataxerr = None
                    datayerr = None
                    print
                    print "plot %i" % i
                    print "marker %s" % option.pt[i % len(option.pt)]
                    print "color %s" % option.pc[i % len(option.pc)]
                    if colyerr != '':
                        datayerr = data[:, colyerr]
                    if colyerr2 != '':
                        datayerr2 = data[:, colyerr2].T
#                        datayerr=np.vstack([datayerr,datayerr2]).T
                        datayerr = data[:, [colyerr, colyerr2]].T
#                        datayerr=np.concatenate([datayerr,datayerr2],axis=1)
                    if colxerr == -1:
                        dataxerr = np.zeros(len(data))
                    elif colxerr != '':
                        dataxerr = data[:, colxerr]
                    if colxerr2 != '':
                        dataxerr = data[:, [colxerr, colxerr2]].T
#                        dataxerr=np.vstack([dataxerr,dataxerr2]).T
#                    print dataxerr[:,0]
#                    sys.exit()
#                        dataxerr=np.concatenate([dataxerr,dataxerr2],axis=1)
#                    datax, datay, datayerr = removeNans3(datax, datay, datayerr)
                    errorbar(datax, datay, xerr=dataxerr, yerr=datayerr, color=option.lc[i % len(option.lc)], ecolor=option.lc[i % len(option.lc)], elinewidth=option.width[i % len(option.width)], barsabove=False, capsize=3, fmt='-', lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)], mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)])

                elif option.plotType[plotTypeIdx] == 'fillbtwX':
                    x1 = datax
                    x2 = data[:, option.x2[i % len(option.x2)]]
                    y = datay
                    p = ax.fill_betweenx(datay, x1, x2, facecolor=option.fc[i % len(option.fc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)])
                        
#                    if option.fhs!='':
#                        from matplotlib.patches import PathPatch
#                        for path in p.get_paths():
#                            p1 = PathPatch(path, fc="none", hatch=option.fhs[i % len(option.fhs)])
#                            ax.add_patch(p1)
#                            p1.set_zorder(p.get_zorder()-0.1)
#                    p = plt.Rectangle((0, 0), 0, 0, facecolor=option.fc[i % len(option.fc)], hatch=option.fhs[i % len(option.fhs)], label=plotLegendLabel) # this is a hack to make labels display correcly
#                    ax.add_patch(p)
                elif option.plotType[plotTypeIdx] == 'fillbtwY':
                    y1 = datay
                    y2 = data[:, option.y2[i % len(option.y2)]]
                    print y1, y2
                    p = ax.fill_between(datax, y1, y2, facecolor=option.fc[i % len(option.fc)], linestyle=option.ls[i % len(option.ls)], edgecolor=option.lc[i % len(option.lc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)])
                        
                    if option.fhs != '':
                        from matplotlib.patches import PathPatch
                        for path in p.get_paths():
                            p1 = PathPatch(path, fc="none", hatch=option.fhs[i % len(option.fhs)])
                            ax.add_patch(p1)
                            p1.set_zorder(p.get_zorder() - 0.1)
                    p = plt.Rectangle((0, 0), 0, 0, facecolor=option.fc[i % len(option.fc)], edgecolor=option.lc[i % len(option.lc)], hatch=option.fhs[i % len(option.fhs)], label=plotLegendLabel)  # this is a hack to make labels display correcly
                    ax.add_patch(p)
#                    plot(datax,datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth, mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)
                    
                elif option.plotType[plotTypeIdx] == 'fnshaded':
                    print
                    print "plot %i" % i
                    print "marker %s" % option.pt[i % len(option.pt)]
                    print "color %s" % option.pc[i % len(option.pc)]
                    print "all markers" 
                    print option.pt
                    print "all colors" 
                    print option.pc
                    print
#                    datax, datay = removeNans2(datax, datay)
#                    print 'colyerr % i' % colyerr
#                    print 'colyerr2 % i' % colyerr2
                    if colyerr2 == '':
                        if option.yerr2sided:
                            y1 = datay - datayerr / 2
                            y2 = datay + datayerr / 2
                        else:
                            y1 = datay - datayerr
                            y2 = datay + datayerr
                    else:
                        y1 = datayerr
                        y2 = datayerr2
                        
                    if option.logY:
                        datax, y1, y2 = removeNegativesValues(datax, y1, y2)
#                    ax.fill_between(datax, y1, y2, facecolor=option.fc[i % len(option.fc)], hatch=option.fhs[i % len(option.fhs)], label=plotLegendLabel) # hatch is not supported in version 1.0 so hack below
                    p = ax.fill_between(datax, y1, y2, facecolor=option.fc[i % len(option.fc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)])
                    if option.fhs != '':
                        from matplotlib.patches import PathPatch
                        for path in p.get_paths():
                            p1 = PathPatch(path, fc="none", hatch=option.fhs[i % len(option.fhs)])
                            ax.add_patch(p1)
                            p1.set_zorder(p.get_zorder() - 0.1)
                    p = plt.Rectangle((0, 0), 0, 0, facecolor=option.fc[i % len(option.fc)], hatch=option.fhs[i % len(option.fhs)], label=plotLegendLabel)  # this is a hack to make labels display correcly
                    ax.add_patch(p)
#                    plot(datax,datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth, mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)
                elif option.plotType[plotTypeIdx] == 'circ':
                    axes(ax)
#                    if option.colColor==-1:
#                        data=arange(len(data))
#                        scatterGlobalData=data
#                    else:
#                        scatterGlobalData=data[:,option.colColor]
                    
#                    if option.colColor==-1:
                    for idx in np.arange(len(datax)):
                        print "generating circle: ", idx
                        circ = Circle((datax[idx], datay[idx]), sizeData[idx], fc=None, ec='k', color=None, fill=False)
                        ax.add_patch(circ)
#                    else:
#                        for idx in np.arange(len(datax)):
#                            print "generating circle: ",idx
#                            circ = Circle((datax[idx],datay[idx]), sizeData[idx], fc=None, ec='k', color=None, fill=False)
#                            ax.add_patch(circ)
                
                elif option.plotType[plotTypeIdx] == 'hist':
                    ###############################################################################################
                    # map type histogram
                    ###############################################################################################
#                    print datax
                    hist(datax, color=option.lc[i % len(option.lc)], bins=option.histNbin, label=plotLegendLabel, cumulative=option.histCuml, normed=option.histNormed)

#                elif option.plotType[plotTypeIdx]=='ts':
#                    datax=inFile[i][:,0]
#                    datay=inFile[i][:,1]
#                    import matplotlib.dates as mdates
#                    import time
#                    datex=list()
#                    [ datex.append(time.strptime(tmpdate, '%M/%d/%Y')) for tmpdate in datax ]
#                    print datax
#                    print datex
#                    plot(datex, datay, lw=option.width[i % len(option.width)], marker=option.pt[i % len(option.pt)], linestyle=option.ls[i % len(option.ls)], c=option.lc[i % len(option.lc)], markersize=option.ps[i % len(option.ps)], mew=option.markerEdgeWidth[i % len(option.markerEdgeWidth)], mec=option.ec[i % len(option.ec)],mfc=option.pc[i % len(option.pc)], label=plotLegendLabel, zorder=option.zorder[i % len(option.zorder)], picker=True)
#
#
#                    years    = mdates.YearLocator()   # every year
#                    months   = mdates.MonthLocator()  # every month
#                    yearsFmt = mdates.DateFormatter('%Y')
#                    # format the ticks
#                    ax.xaxis.set_major_locator(years)
#                    ax.xaxis.set_major_formatter(yearsFmt)
#                    ax.xaxis.set_minor_locator(months)
# #                    datemin = datetime.date(datex.min().year, 1, 1)
# #                    datemax = datetime.date(datex.max().year+1, 1, 1)
# #                    ax.set_xlim(datemin, datemax)
# #                    ax.format_xdata = mdates.DateFormatter('%Y-%m-%d')
# #                    fig.autofmt_xdate()


    
                if option.interactive:
                    cid = fig.canvas.mpl_connect('key_press_event', on_key)
                    fig.canvas.mpl_connect('pick_event', on_pick)
                    fig.canvas.mpl_connect('button_press_event', on_mouse_button_down)

    
            
                #
                # plot mask ranges
                #
                for maskRange in maskRanges:
                    axvspan(maskRange[0], maskRange[1], facecolor='b', edgecolor=None, alpha=0.3)
                    maskMaskedPointsMapRange(maskRange[0], maskRange[1])
        
#                figIdx=int(i/option.dsPerPlot)
#                if i % option.dsPerPlot==0:
#                    vi=0
#                    for v in HlinesValues:
#                        if int(vi/option.dsPerPlot) == figIdx:
#                            axhline(y=v, xmin=0, xmax=1, color=option.olc[vi % len(option.olc)], lw=option.width[vi % len(option.width)],label=None)
# #                        text(fx, -65, HLlabels[vi], withdash=False, rotation=90, verticalalignment='top', horizontalalignment='center')
#                        vi+=1
                        

#                figIdx=int(i/option.dsPerPlot)
#                if option.plotDs!=None:
#                    dummy,dsInFigIdxfigIdx=calculatePlotID(i,figRows,figCols)
                
                for Vspidx in range(len(option.Hline)):
                    if Vspidx == figIdx:
                        vi = 0
                        for vidx in range(len(option.Hline[Vspidx])):
                            axhline(y=option.Hline[Vspidx][vidx], xmin=0, xmax=1, color=option.olc[Vspidx % len(option.olc)][vidx % len(option.olc[Vspidx % len(option.olc)])], lw=option.olw[Vspidx % len(option.olw)][vidx % len(option.olw[Vspidx % len(option.olw)])], ls=option.ols[Vspidx % len(option.ols)][vidx % len(option.ols[Vspidx % len(option.ols)])], label=None)
                            vi += 1

                for Vspidx in range(len(option.Vline)):
                    if Vspidx == figIdx:
                        vi = 0
                        for vidx in range(len(option.Vline[Vspidx])):
#                            axvline(x=option.Vline[Vspidx][vidx], ymin=0, ymax=1, color=option.olc[vi % len(option.olc)], lw=1, ls=option.ols[vi% len(option.ols)], label=None)
                            axvline(x=option.Vline[Vspidx][vidx], ymin=0, ymax=1, color=option.olc[Vspidx % len(option.olc)][vidx % len(option.olc[Vspidx % len(option.olc)])], lw=option.olw[Vspidx % len(option.olw)][vidx % len(option.olw[Vspidx % len(option.olw)])], ls=option.ols[Vspidx % len(option.ols)][vidx % len(option.ols[Vspidx % len(option.ols)])], label=None)
                            vi += 1
                
                
                for Vspidx in range(len(option.Vspan)):
                    if Vspidx == figIdx:
                        
                        for vidx in range(len(option.Vspan[Vspidx]) / 2):
                            print 'printing vspan vi: ', vidx
                            x1 = option.Vspan[Vspidx][vidx * 2]
                            x2 = option.Vspan[Vspidx][vidx * 2 + 1]
                            axvspan(xmin=x1, xmax=x2, ymin=0, ymax=1, alpha=option.oalpha[vidx % len(option.oalpha)], color=option.ofc[Vspidx][vidx % len(option.ofc[Vspidx])], lw=option.olw[Vspidx % len(option.olw)][vidx % len(option.olw[Vspidx % len(option.olw)])], label=None)
#                        text(fx, -65, HLlabels[vi], withdash=False, rotation=90, verticalalignment='top', horizontalalignment='center')
    

                for Hspidx in range(len(option.Hspan)):
                    if Hspidx == figIdx:
                        
                        for vidx in range(len(option.Hspan[Hspidx]) / 2):
                            print 'printing vspan vi: ', vidx
                            x1 = option.Hspan[Hspidx][vidx * 2]
                            x2 = option.Hspan[Hspidx][vidx * 2 + 1]
                            axhspan(ymin=x1, ymax=x2, xmin=0, xmax=1, alpha=option.oalpha[vidx % len(option.oalpha)], color=option.ofc[Hspidx][vidx % len(option.ofc[Hspidx])], lw=option.olw[Hspidx % len(option.olw)][vidx % len(option.olw[Hspidx % len(option.olw)])], label=None, zorder=-100)
#                        text(fx, -65, HLlabels[vi], withdash=False, rotation=90, verticalalignment='top', horizontalalignment='center')

                #
                # extra texts
                #
#                 tstrIdx=0
#                 for tstr in option.textxy:
#                     if tstrIdx==dsIdxInSubplot and newSubplot:
                if isLastDsInSubplot:
                    if len(option.textxy) > 0:
                        tstr = option.textxy[figIdx]
                    else:
                        tstr = 'None'

                    # do the first split for all texts that we want on this figure (if any)
                    if tstr != 'None':
                        xyTextTuple = tstr.split(';')
                        print 'xyTextTuple', xyTextTuple
                        for tstr in xyTextTuple:
                            tstrSplit = re.split(r'(?<!\\),', tstr)  # negative lookbehind assertion
                            if len(tstrSplit) >= 3:
                                tx = float(tstrSplit[0])
                                ty = float(tstrSplit[1])
                                t = tstrSplit[2].replace('\,', ',')
                                trot = 0
                                tsiz = option.extraTextFontSize
                                tcolor = 'k'
                            if len(tstrSplit) >= 4:
                                trot = tstrSplit[3]
                            if len(tstrSplit) >= 5:
                                tsiz = tstrSplit[4]
                            if len(tstrSplit) >= 6:
                                tcolor = tstrSplit[5]
                                
                            text(tx, ty, t, fontsize=tsiz, zorder=10000, rotation=trot, color=tcolor)
                        
#                     tstrIdx+=1
                    if len(option.textfxfy) > 0:
                        tstr = option.textfxfy[figIdx % len(option.textfxfy)]
                        if tstr != 'None':
                            xyTextTuple = tstr.split(';')
                            print 'xyTextTuple', xyTextTuple
                            for tstr in xyTextTuple:
#                                 tstrSplit=tstr.split(',')
                                tstrSplit = re.split(r'(?<!\\),', tstr)  # negative lookbehind assertion
                                if len(tstrSplit) >= 3:
                                    tx = float(tstrSplit[0])
                                    ty = float(tstrSplit[1])
                                    t = tstrSplit[2].replace('\,', ',')
                                    trot = 0
                                    tsiz = option.extraTextFontSize
                                    tcolor = 'k'
                                if len(tstrSplit) >= 4:
                                    trot = tstrSplit[3]
                                if len(tstrSplit) >= 5:
                                    tsiz = tstrSplit[4]
                                if len(tstrSplit) >= 6:
                                    tcolor = tstrSplit[5]

                                ax.annotate(t, xy=(tx, ty), xycoords='axes fraction', fontsize=tsiz, horizontalalignment='left', verticalalignment='top', rotation=trot, color=tcolor)

                    if len(option.annotate) > 0:
                        tstr = option.annotate[figIdx]
                        if tstr != 'None':
                            xyTextTuple = tstr.split(';')
                            print 'xyTextTuple', xyTextTuple
                            for tstr in xyTextTuple:
                                tstrSplit = re.split(r'(?<!\\),', tstr)  # negative lookbehind assertion
                                
                                print tstrSplit
#                                 sys.exit()
                                if len(tstrSplit) >= 5:
                                    arrowx = float(tstrSplit[0])
                                    arrowy = float(tstrSplit[1])
                                    tx = float(tstrSplit[2])
                                    ty = float(tstrSplit[3])
                                    t = tstrSplit[4].replace('\,', ',')
                                    trot = 0
                                    tsiz = option.extraTextFontSize
                                    tcolor = 'k'
                                    
                                if len(tstrSplit) >= 6:
                                    trot = tstrSplit[5]
                                if len(tstrSplit) >= 7:
                                    tsiz = tstrSplit[6]
                                if len(tstrSplit) >= 8:
                                    tcolor = tstrSplit[7]

                                ax.annotate(t, xy=(arrowx, arrowy), xytext=(tx, ty), xycoords='data', textcoords='axes fraction', arrowprops=dict(facecolor=tcolor, shrink=0.01), fontsize=tsiz, horizontalalignment='left', verticalalignment='top', rotation=trot, color=tcolor)


                #
                # limits business
                #
                if option.plotType[plotTypeIdx] != 'hist' and option.plotType[plotTypeIdx] != 'ts' and option.plotType[plotTypeIdx] != 'barchartH' and option.dateFmt == '':
                    xmin, xmax, ymin, ymax = getDataMinMaxValues(ax, datax, datay, i)
                    print xmin, xmax, ymin, ymax
                    if option.nestedGridLevel >= 0:
                        if option.nestedGridDset == i or option.nestedGridDset == -1:
                            makeNestedGridPlot(getDataMinMaxValues2(datax, datay))

                    marginFactor = option.marginFactor
                    if (option.logX or option.logY) and marginFactor > 0:
                        print 'log scale requested and margin factor > 0, will change to 0 for safety'
                        marginFactor = 0
#                    if xmin <= 0 and option.logX:
#                        xmin=0;
#                    if ymin <= 0 and option.logY:
#                        xmin=0;
                    
                    if option.xmin[figIdx % len(option.xmin)] != None:
                        xmin = option.xmin[figIdx % len(option.xmin)]
                    if option.xmax[figIdx % len(option.xmax)] != None:
                        xmax = option.xmax[figIdx % len(option.xmax)]
                
                    if option.ymin[figIdx % len(option.ymin)] != None:
                        ymin = option.ymin[figIdx % len(option.ymin)]
                    if option.ymax[figIdx % len(option.ymax)] != None:
                        ymax = option.ymax[figIdx % len(option.ymax)]

                    xmin = xmin - option.marginX
                    xmax = xmax + option.marginX
                    ymin = ymin - option.marginY
                    ymax = ymax + option.marginY
                    print 'xmin,xmax,ymin,ymax: with margin'
                    print xmin, xmax, ymin, ymax
                    print 'setting x limits: %f, %f' % (xmin, xmax)
                    xlim([xmin, xmax])
                    print 'setting y limits: %f, %f' % (ymin, ymax)
                    ylim([ymin, ymax])
                    
                    
                if option.plotType[plotTypeIdx] == 'scat3d':
#                     if option.zmin!= -1:
                    zmin = option.zmin
#                     if option.zmax!= -1:
                    zmax = option.zmax
                    ax.set_xlim3d(xmin,xmax)
                    ax.set_ylim3d(ymin,ymax)
                    ax.set_zlim3d(zmin,zmax)

                if option.dateFmt != '':
                    formatXaxisDates(figIdx)

                if option.plotType[plotTypeIdx] == 'ts':
                    xmin, xmax, ymin, ymax = getDataMinMaxValues(ax, datax, datay, i)
#                     if option.xmin[figIdx % len(option.xmin)]!=-1:
#                         xmin=option.xmin[figIdx % len(option.xmin)]
#                     if option.xmax[figIdx % len(option.xmax)]!=-1:
#                         xmax=option.xmax[figIdx % len(option.xmax)]
                
                    if option.ymin[figIdx % len(option.ymin)] != None:
                        ymin = option.ymin[figIdx % len(option.ymin)]
                    if option.ymax[figIdx % len(option.ymax)] != None:
                        ymax = option.ymax[figIdx % len(option.ymax)]

#                     xmin=xmin-option.marginX
#                     xmax=xmax+option.marginX
                    ymin = ymin - option.marginY
                    ymax = ymax + option.marginY
#                     print 'xmin,xmax,ymin,ymax: with margin'
#                     print xmin,xmax,ymin,ymax
#                     print 'setting x limits: %f, %f' % (xmin,xmax)
#                     xlim([xmin,xmax])
                    print 'setting y limits: %f, %f' % (ymin, ymax)
                    ylim([ymin, ymax])
                #
                # log scales
                # 
#                if i % option.dsPerPlot == 0:
                if newSubplot:
                    if option.logX:
                        ax.set_xscale('log')
                    if option.logY:
                        ax.set_yscale('log')

                    global globalWhichGridLines
                    whichGridLines = globalWhichGridLines
                    if option.logX or option.logY:
                        whichGridLines = 'both'
                    if len(option.logYaxes) > 0:
                        if option.logYaxes[figIdx % len(option.logYaxes)] == 'log':
                            ax.set_yscale('log')
                            whichGridLines = 'both'
                        if option.logYaxes[figIdx % len(option.logYaxes)] == 'lin':
                            ax.set_yscale('linear')
                            whichGridLines = 'major'







                    
                #
                # title and labels
                #
                if isLastDsInSubplot:
                    print 'figIdx: ', figIdx
                    print 'title: ', option.title
                    print 'title: ', option.title[figIdx % len(option.title)]
                    if option.title[figIdx % len(option.title)] != "":
                        tit = option.title[figIdx % len(option.title)].decode('utf8')
                        title(tit, fontsize=option.fontSizeTitle, horizontalalignment=option.titleAlignH, verticalalignment=option.titleAlignV)

#                if i % option.dsPerPlot == 0:
#                if newSubplot:
                if isLastDsInSubplot:
#                print option.xlabel
#                cpedsPythCommon.waitEnter()
#                    labidx=i/option.dsPerPlot
                    labidx = figIdx
                    print 'labidx ', labidx
                    print 'i ', i
                    print 'option.fontSize[(2*labidx) % len(option.fontSize)]: ', option.fontSize[(2 * labidx) % len(option.fontSize)]
                    print 'option.fontSize[(2*labidx+1) % len(option.fontSize)]: ', option.fontSize[(2 * labidx + 1) % len(option.fontSize)]
                    if option.xlabel[labidx % len(option.xlabel)] != "None":
                        xlabel(option.xlabel[labidx % len(option.xlabel)].decode('utf8'), fontsize=option.fontSizeLabels[(2 * labidx) % len(option.fontSizeLabels)])
#                    if option.plotType[plotTypeIdx]!='ts':
                    if option.fontSize[(2 * labidx) % len(option.fontSize)] > 0:
                        setp(ax.get_xticklabels(), fontsize=option.fontSize[(2 * labidx) % len(option.fontSize)])
                    else:
                        print 'switching off xlabels'
                        ax.set_xticklabels([])
                        print 'done'
    #                        ax.set_xticks([])
                        
                    if option.ylabel[labidx % len(option.ylabel)] != "None":
                        ylabel(option.ylabel[labidx % len(option.ylabel)].decode('utf8'), fontsize=option.fontSizeLabels[(2 * labidx + 1) % len(option.fontSizeLabels)])
                    if option.fontSize[(2 * labidx + 1) % len(option.fontSize)] > 0:
                        setp(ax.get_yticklabels(), fontsize=option.fontSize[(2 * labidx + 1) % len(option.fontSize)])
                    else:
                        print 'switching off ylabels'
                        ax.set_yticklabels([])
#                        ax.set_yticks([])

                    print
                    print 'Setting xticklabels rotation'    
                    print
                    setp(ax.get_xticklabels(), rotation=option.Rxlabels)  # 
                    setp(ax.get_yticklabels(), rotation=option.Rylabels)  # , horizontalalignment='center', verticalalignment='top'
    
                        
    
                #
                # ticks business
                #
                if option.plotType[plotTypeIdx] != 'ts':
                    if option.xticks[figIdx % len(option.xticks)] != 0:
                        tcs = ax.get_xticks()
                        if option.xmin[figIdx % len(option.xmin)] != None:
                            x1 = option.xmin[figIdx % len(option.xmin)]
                        else:
                            x1 = tcs[0]
                        if option.xmax[figIdx % len(option.xmax)] != None:
                            x2 = option.xmax[figIdx % len(option.xmax)]
                        else:
                            x2 = tcs[-1]
                        
                        xticks(arange(x1, x2, option.xticks[figIdx % len(option.xticks)]))
#                         xticks(arange(min(tcs[0],option.xmin[figIdx % len(option.xmin)]),max(tcs[-1],option.xmax[figIdx % len(option.xmax)]),option.xticks[figIdx % len(option.xticks)]))
                    if option.yticks[figIdx % len(option.yticks)] != 0:
                        tcs = ax.get_yticks()
                        if option.ymin[figIdx % len(option.ymin)] != None:
                            y1 = option.ymin[figIdx % len(option.ymin)]
                        else:
                            y1 = tcs[0]
                        if option.ymax[figIdx % len(option.ymax)] != None:
                            y2 = option.ymax[figIdx % len(option.ymax)]
                        else:
                            y2 = tcs[-1]
                        yticks(arange(y1, y2, option.yticks[figIdx % len(option.yticks)]))
#                         yticks(arange(min(tcs[0], option.ymin[figIdx % len(option.ymin)]), max(tcs[-1], option.ymax[figIdx % len(option.ymax)]), option.yticks[figIdx % len(option.yticks)]))


                    if option.xticksMinor[figIdx % len(option.xticks)] != None:
                        minorLocatorX = MultipleLocator(option.xticksMinor[figIdx % len(option.xticksMinor)])
                        ax.xaxis.set_minor_locator(minorLocatorX)
                    if option.yticksMinor[figIdx % len(option.yticks)] != None:
                        minorLocatorY = MultipleLocator(option.yticksMinor[figIdx % len(option.yticksMinor)])
                        ax.yaxis.set_minor_locator(minorLocatorY)
#                     ax1.tick_params(axis='x',which='major',direction='out',length=4,width=4,color='b',pad=10,labelsize=20,labelcolor='g')
                    ax.tick_params(axis='x',which='major',direction='in',length=8,width=2,color='k',pad=10)
                    ax.tick_params(axis='y',which='major',direction='in',length=8,width=2,color='k',pad=10)
                    ax.tick_params(axis='x',which='minor',direction='in',length=4,width=1,color='k',pad=10)
                    ax.tick_params(axis='y',which='minor',direction='in',length=4,width=1,color='k',pad=10)

                    
                    
#                 setp( ax.get_xticklabels(), fontsize=option.fontSizeXtickLabels[i % len(option.fontSizeXtickLabels)])
#                 setp( ax.get_yticklabels(), fontsize=option.fontSizeYtickLabels[i % len(option.fontSizeYtickLabels)])
    
                if isLastDsInSubplot:
#                     if len(option.removeXtickLabels)>0:
                    if option.removeXtickLabels[figIdx % len(option.removeXtickLabels)] > 0:
                        setp(ax.get_xticklabels(), visible=False)
                        setp(ax.get_xticklabels()[::option.removeXtickLabels[figIdx % len(option.removeXtickLabels)]], visible=True)
#                     if len(option.removeYtickLabels)>0:
                    if option.removeYtickLabels[figIdx % len(option.removeYtickLabels)] > 0:
                        setp(ax.get_yticklabels(), visible=False)
                        setp(ax.get_yticklabels()[::option.removeYtickLabels[figIdx % len(option.removeYtickLabels)]], visible=True)
            
                    #
                    # broken axis stuff
                    #
                    if option.brokenYaxis:
                        d = .005  # how big to make the diagonal lines in axes coordinates
                        # arguments to pass to plot, just so we don't keep repeating them
                        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)

                        if figIdx == 0:
                            ax.spines['bottom'].set_visible(False)
                            ax.tick_params(labelbottom='off')
                            ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
                            ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
                        if figIdx == 1:
                            ax.spines['top'].set_visible(False)
                            ax.tick_params(labeltop='off')
                            kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
                            ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
                            ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
#                             if figIdx<option.figRows-1:
#                                 ax.spines['bottom'].set_visible(False)
#                                 ax.tick_params(labeltop='off')
#                                 kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
#                                 ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
#                                 ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
                                
    
    
                #
                # grid 
                #        
                if option.grid:
                    grid(True, whichGridLines, color=option.gridColor)
                
                #
                # vertical lines 
                #
                if option.plotVerticalLinesWithLabelsFromFile != "":
                    plotVerticalLinesWithLabelsFromFile(fig, ax)
                if option.plotLabelsFromFile != "":
                    plotLabelsFromFile(fig, ax, option.plotType[plotTypeIdx])
                
                
                #
                # legend
                #

                if option.plotType[plotTypeIdx] == 'scat':
                    if len(option.label) > 0:
                        # get handles
                        handles, labels = ax.get_legend_handles_labels()
                        # remove the errorbars
#                        handles[plotTypeIdx]=handles[plotTypeIdx][0]
                        if figRows == 1 and figCols == 1:
                            if i == Nfiles - 1:
#                                legend(handles, labels, loc=legendLocation, prop={'size':option.fontSizeLegend[int(i/option.dsPerPlot) % len(option.fontSizeLegend)]})
                                legend(handles, labels, loc=legendLocation[figIdx % len(legendLocation)], prop={'size':option.fontSizeLegend[figIdx % len(option.fontSizeLegend)]}, ncol=option.legendNcol, mode=option.legendMode, borderaxespad=option.legendBorderAxesPad)
                        else:
#                            legend(handles, labels, loc=legendLocation, prop={'size':option.fontSizeLegend[int(i/option.dsPerPlot) % len(option.fontSizeLegend)]})
                            legend(handles, labels, loc=legendLocation[figIdx % len(legendLocation)], prop={'size':option.fontSizeLegend[figIdx % len(option.fontSizeLegend)]}, ncol=option.legendNcol, mode=option.legendMode, borderaxespad=option.legendBorderAxesPad)
                                
                else:
                    if len(option.label) > 0:
                        # get handles
                        handles, labels = ax.get_legend_handles_labels()
                        # remove the errorbars
#                        handles[plotTypeIdx]=handles[plotTypeIdx][0]

                        if figRows == 1 and figCols == 1:
                            if i == Nfiles - 1:
#                                legend(handles, labels, loc=legendLocation, handlelength=option.legendPatternLength, prop={'size':option.fontSizeLegend[int(i/option.dsPerPlot) % len(option.fontSizeLegend)]})
                                legend(handles, labels, loc=legendLocation[figIdx % len(legendLocation)], handlelength=option.legendPatternLength, prop={'size':option.fontSizeLegend[figIdx % len(option.fontSizeLegend)]}, ncol=option.legendNcol, mode=option.legendMode, borderaxespad=option.legendBorderAxesPad)
                        else:
#                            legend(handles, labels, loc=legendLocation, handlelength=option.legendPatternLength, prop={'size':option.fontSizeLegend[int(i/option.dsPerPlot) % len(option.fontSizeLegend)]})
                            if legendLocation[figIdx % len(legendLocation)] != 'None':
                                legend(handles, labels, loc=legendLocation[figIdx % len(legendLocation)], handlelength=option.legendPatternLength, prop={'size':option.fontSizeLegend[figIdx % len(option.fontSizeLegend)]}, ncol=option.legendNcol, mode=option.legendMode, borderaxespad=option.legendBorderAxesPad)
#                                 legend.remove()

                #
                # colorbar
                # 
#                if option.plotType[plotTypeIdx]=='scat':
#                    global globalScatterColorbar
#                    if option.colorbar:
#                        if globalScatterColorbar==None:
#                            cmax=max(scatterGlobalData)
#                            cmin=min(scatterGlobalData)
#                            cdelta=(max(scatterGlobalData)-min(scatterGlobalData))/(option.cbLabelsNum-1)
#                            cbyticks=np.arange(cmin,cmax+cdelta,cdelta)
#                            print cbyticks
#                            cb=colorbar(ticks=cbyticks)
#                            globalScatterColorbar=cb
#                        else:
#                            xmin,xmax,ymin,ymax=globalScatterColorbar.ax.axis()
#                            print "getDataMinMaxValues: ax"
#                            print xmin,xmax,ymin,ymax
#                            sys.exit(0)
#                                     
#                            
#                            globalScatterColorbar=colorbar(ticks=cbyticks)
#                        
#                        cb.set_label(option.zlabel)
                    
                if option.plotType[plotTypeIdx] == 'vect' or option.plotType[plotTypeIdx] == 'scat' or option.plotType[plotTypeIdx] == 'scatContFill' or option.plotType[plotTypeIdx] == 'scatDensContFill':
                    if option.colorbar:
                        labidx = figIdx

                        cmax = max(scatterGlobalData)
                        cmin = min(scatterGlobalData)
                        if option.vmin[figIdx % len(option.vmin)] != option.vmax[figIdx % len(option.vmin)]:
#                         if option.vmin[figIdx % len(option.vmin)] != None and option.vmax[figIdx % len(option.vmin)] != None:
                            cmin = option.vmin[figIdx % len(option.vmin)]
                            cmax = option.vmax[figIdx % len(option.vmin)]
                        cdelta = (cmax - cmin) / (option.cbLabelsNum - 1)
                        cbyticks = np.arange(cmin, cmax + cdelta, cdelta)
                        print cbyticks
                        cb = colorbar(ticks=cbyticks)
#                         cb.ax.tick_params(labelsize=option.fontSize[i % len(option.fontSize)])
                        if option.fontSizeCM[figIdx % len(option.fontSizeCM)] > 0:
                            cb.ax.tick_params(labelsize=option.fontSizeCM[figIdx % len(option.fontSizeCM)])
                        else:
                            print 'switching off zlabels'
                            cb.ax.set_yticklabels([])



#                         cb.set_label(option.zlabel, size=option.fontSizeLabels[i % len(option.fontSizeLabels)])
                        cb.set_label(option.zlabel[figIdx % len(option.zlabel)], size=option.fontSizeCM[figIdx % len(option.fontSizeCM)])
                        
                        
                        
    #                    if option.colorbar:
                        if option.linZtickLabels:
                            print 'converting color domain labels from log space back to linear space'
                            cax = cb.ax
    #                        axes(cax)
    #                        print get_yticks()
    #                        cax.ylabels
    #                        cbyticks=cax.get_yticks()
    #                        print cbyticks
                            cbyticksLabelsNew = list()
                            for lab in cbyticks:
                                cbyticksLabelsNew.append(option.ZticksFmt % np.power(10, float(lab)))
    #                        yticks(cbyticks)
                            print 'the new colorbar labels will be'
                            print cbyticksLabelsNew
    #                        axes(ax)
                            
    #                        sys.exit()
    #                        cbyticksLabelsNew=list()
    #                        for lab in cbyticksLabels:
    #                            cbyticksLabelsNew.append('%.2f' % np.power(10,float(lab)))
                            cax.set_yticklabels(cbyticksLabelsNew)
                        else:
                            cax = cb.ax
                            cbyticksLabelsNew = list()
                            for lab in cbyticks:
                                cbyticksLabelsNew.append(option.ZticksFmt % float(lab))
                            cax.set_yticklabels(cbyticksLabelsNew)


                    if option.CM[figIdx % len(option.CM)]=='None':
                        fig.delaxes(cax)

                #
                # circles            
                #
                if len(circlePatchDef) > 0:
                    axes(ax)
                    plotCircles(ax, circlePatchDef)


                #
                # rectangles
                #
                print option.rect
                if len(option.rect) > 0:
                    rectanglePatchDef = option.rect[i % len(option.rect)]
                        
                    for rect_i in range(len(rectanglePatchDef)):
                        axes(ax)
                        rectX = rectanglePatchDef[rect_i][0]
                        rectY = rectanglePatchDef[rect_i][1]
                        rectW = rectanglePatchDef[rect_i][2]
                        rectH = rectanglePatchDef[rect_i][3]
                        rectA = rectanglePatchDef[rect_i][4]
                        print option.olc[i % len(option.olc)]
    #                         rect = mpatches.Rectangle((rectX[rect_i % len(rectX)]-rectW[rect_i % len(rectX)]/2,rectY[rect_i % len(rectX)]-rectH[rect_i % len(rectX)]/2), rectW[rect_i % len(rectX)],rectH[rect_i % len(rectX)], angle=rectA[rect_i % len(rectX)], fc=None, ec=option.olc[rect_i % len(option.olc)][0], color=None, fill=False)
                        rect = mpatches.Rectangle((rectX - rectW / 2, rectY - rectH / 2), rectW, rectH, angle=rectA, fc=None, ec=option.olc[rect_i % len(option.olc)][0], color=None, fill=False)
                        ax.add_patch(rect)
            


                if option.plotType[plotTypeIdx] != 'hist' and option.plotType[plotTypeIdx] != 'barchartH' and option.plotType[plotTypeIdx] != 'ts' and option.dateFmt == '':
                    print 'Reseting limits'
                    xlim([xmin, xmax])
                    ylim([ymin, ymax])


                if isLastDsInSubplot:
                    if len(option.x2)>0:
                        x2 = data[:, option.x2[figIdx % len(option.x2)]]
                        ax2=ax.twiny()
    
                        ax1_ticks = ax.get_xticks()
                        ax2.set_xticks(ax1_ticks)
                        print 'ax1_ticks: ',ax1_ticks
                        print 'ax2_ticks: ',ax2.get_xticks()
                        ax2.set_xticklabels(ax.get_xticklabels()) # this is a placeholder because the ticklabes are not defined until show()
                        print 'ax2_xticklabels: ',ax2.get_xticklabels()
                        ax2_tickLabels=list()
                        idx=range(len(datax))
                        for t,tick_idx in zip(ax1_ticks,range(len(ax1_ticks))):
                            # hide labels that have interpolated values from outside the defined range
                            # we use constant extrapolation which would give the labels wrong values, so we hide them
                            if (t<min(datax)):
                                print 'setting %i invisible (t: %f, min: %f, len: %i)' % (tick_idx,t,min(datax),len(ax2.get_xticklabels()))
                                setp(ax2.get_xticklabels()[tick_idx], visible=False)
                            elif (t>max(datax)):
                                print 'setting %i invisible (t: %f, max: %f, len: %i)' % (tick_idx,t,max(datax), len(ax2.get_xticklabels()))
                                setp(ax2.get_xticklabels()[tick_idx], visible=False)
                            ti=np.interp(t, datax,idx)
    # 
                            ax2_tickLabels.append('%.1f' % np.interp(ti,idx,x2))
    #                         ax2_tickLabels.append('%.1f' % x2[int(round(ti))])
    #                         ax2_tickLabels.append(x2[int(ti)])
                            print t,ti,np.interp(ti,idx,x2)
    #  
                        print 'ax ticks: ',ax.get_xticks()
                        print ax2_tickLabels
    #                    ax2.set_xticks(ax1_ticks) # set the locations of the xticks at the same places as in the initial x axis
                        ax2.set_xticks(ax1_ticks)
                        ax2.set_xlim(ax.get_xlim())
                        ax2.set_xticklabels(ax2_tickLabels)
                        setp( ax2.get_xticklabels(), fontsize=option.fontSize[figIdx % len(option.fontSize)])
                        
                        axes(ax)
                        if option.x2label[i % len(option.x2label)]!=None:
                            ax2.set_xlabel(option.x2label[figIdx % len(option.x2label)],fontsize=option.fontSizeLabels[figIdx % len(option.fontSizeLabels)])

            if option.plotSpecial != -1:
                import specialBlock
                specialBlock.plotSpecialBlock(ax, option.plotSpecial)

        if option.autolimits:
            ax.relim()
            # update ax.viewLim using the new dataLim
            ax.autoscale(tight=False)
        
        if option.axesSwitch[figIdx % len(option.axesSwitch)] == 0:
            axis('off')
        if option.noAxes:
            axis('off')
    
        if option.equalAspect:    
#             ax.set_aspect('equal', 'datalim')
            ax.set_aspect('equal', 'box')
#             axes().set_aspect('equal', adjustable='box', anchor='C')

    
####################################################################################
# INTERACTIVE STUFF    
    
    if option.interactive: 





#         from matplotlib.widgets import CheckButtons
#         check = CheckButtons(rax, ('2 Hz', '4 Hz', '6 Hz'), (False, True, True))


        
        
#        plotButtons(Bax)
#        if spannerOn==False:
#            span = SpanSelector(ax1, onselect, 'horizontal', useblit=True, rectprops=dict(alpha=0.5, facecolor='g') )
#        cursor = Cursor(ax1, useblit=True, color='green', linewidth=2 )
        Bax = plt.axes([0.9, 0.03, 0.1, 0.015])
        bExtract = Button(Bax, 'Extract Spikes')
        bExtract.on_clicked(extractSpikes)
        Bax = plt.axes([0.9, 0.015, 0.1, 0.015])
        bBinConfig = Button(Bax, 'Configure binning')
        bBinConfig.on_clicked(makeBinningInfoWidget)

        Bax = plt.axes([0.8, 0.03, 0.1, 0.015])
        bDeleteLines = Button(Bax, 'Delete all lines')
        bDeleteLines.on_clicked(deleteAllMaskLines)
        Bax = plt.axes([0.8, 0.015, 0.1, 0.015])
        bDeleteLine = Button(Bax, 'Delete line')
        bDeleteLine.on_clicked(deleteMaskLine)

        Bax = plt.axes([0.7, 0.015, 0.1, 0.03])
        bDeleteMaskRange = Button(Bax, 'Delete Mask Range')
        bDeleteMaskRange.on_clicked(deleteMaskRange)

        Bax = plt.axes([0.6, 0.015, 0.1, 0.03])
        bSaveLastBinnedSpectra = Button(Bax, 'Save last binned spectra')
        bSaveLastBinnedSpectra.on_clicked(saveLastBinnedSpectraWidget)

        Bax = plt.axes([0.5, 0.015, 0.1, 0.03])
        bSaveLastBinnedSpectra = Button(Bax, 'Delete All Mask ranges')
        bSaveLastBinnedSpectra.on_clicked(deleteAllMaskRanges)
    
    
#     ax.xaxis.set_label_coords(0.5, 0.5, transform = fig.transFigure)
    
#     ax.xlim(-0.4,0.4)
#     ax.axis([-1.4,1.4,0,10]) 


    if option.continuousPotting:
        if option.contAutoRescale:
            ax.relim()
            # update ax.viewLim using the new dataLim
            ax.autoscale_view()
        draw()
        if option.waitN:
            tm.sleep(option.waitN)
        if option.clearFig:
            clf()
#             continuousPottingFirstTime=True
        return fig, ax
    else:
    
        if option.save:
            if option.outputFile != '':
                fname = option.outputFile
                fig.savefig(fname, dpi=option.DPI, transparent=option.transparent)
            else:
                print "no output file name given, will not save"
                
            if option.show:
                show()
                
        else:
            show()
            
        
    return fig, ax
######################################################################################
    
    
    
    
    
    
    
    
    
    
    
    
def plotLabelsFromFile(fig, ax, plotType):
    print "-------------------------------"
    print "plotting labels from file"
    print "-------------------------------"
    # first read the data from the file
    f = open(option.plotLabelsFromFile, "r")
    c1 = list()
    c2 = list()
    Lon = list()
    Lat = list()
    txt = list()
    ps = list()
    pc = list()
    pt = list()
    ha = list()
    va = list()
    labelSize = list()
    for line in f:
        if line[0] != '#':
#            s=line.split([' ','\n'])
            s = line.decode("utf-8").split()
            print s
            print s[0]
            c1.append(s[0])
            c2.append(s[1])
            ps.append(int(s[2]))
            pc.append(s[3])
            pt.append(s[4])
            labelSize.append(int(s[5]))
            ha.append(s[6][1])
            va.append(s[6][0])
            s.pop(0)
            s.pop(0)
            s.pop(0)
            s.pop(0)
            s.pop(0)
            s.pop(0)
            s.pop(0)
            S = ''
            for part in s:
                S = S + ' ' + part
            txt.append(S)
    f.close()
    print ha
    print va

#    rcParams["text.usetex"] = False
#    import matplotlib.font_manager as fm
#    fp1=fm.FontProperties(fname="/usr/share/fonts/liberation/LiberationSans-Regular.ttf") 
    matplotlib.rc('font', **{'sans-serif' : 'Arial', 'family' : 'sans-serif'})
    hAlign = 'left'
    vAlign = 'top'
    if plotType == "sphere" or plotType == "map": 
        for iloc in np.arange(len(c1)):
#            Lon=list()
#            Lat=list()
            if option.reverseLon:
                xpt, ypt = projectMap(-c1[iloc], c2[iloc])
            else:
                xpt, ypt = projectMap(c1[iloc], c2[iloc])
            print ha[iloc]
            print va[iloc]
            if ha[iloc] == 'l':
                hAlign = 'left'
            if ha[iloc] == 'r':
                hAlign = 'right'
            if ha[iloc] == 'c':
                hAlign = 'center'
            if va[iloc] == 't':
                vAlign = 'top'
            if va[iloc] == 'b':
                vAlign = 'bottom'
            if va[iloc] == 'c':
                hAlign = 'center'
            print hAlign, vAlign
#            rc('text', fontsize=labelSize[iloc], color='#000000')
#            rc('xtick', labelsize=15, color='#000000')
#            rc('ytick', labelsize=15, color='#000000')
#    rc('text', labelsize=15, color='#000000')
#            text(xpt,ypt,txt[iloc],fontsize=labelSize[iloc], horizontalalignment=hAlign, verticalalignment=vAlign, fontproperties=fp1)
#            rcParams['text.usetex'] = True
            text(xpt, ypt, unicode(txt[iloc]), fontsize=labelSize[iloc], horizontalalignment=hAlign, verticalalignment=vAlign)
#            Lon.append(xpt)
#            Lat.append(ypt)

#            scatter(Lon,Lat, c=pc[iloc], s=ps[iloc], lw=option.markerEdgeWidth , marker=pt[iloc])
            scatter([xpt], [ypt], c=pc[iloc], s=ps[iloc], lw=option.markerEdgeWidth[i % len(option.markerEdgeWidth)] , marker=pt[iloc])
            
    else:
        print "plotLabelsFromFile(fig,ax,plotType): IS NOT SUPPORTED FOR OTHER PLOT TYPES THAN sphere YET, SORRY."    
        sys.exit()
    


def plotVerticalLinesWithLabelsFromFile(fig, ax):
    # first read the data from the file
    f = open(option.plotVerticalLinesWithLabelsFromFile, "r")
    freq = list()
    labels = list()
    for line in f:
        s = line.split(' \t')
        freq.append(s[0])
        s.pop(0)
        S = ''
        for part in s:
            S = S + ' ' + part
        labels.append(S)
    f.close()
    print freq
    print labels
    i = 0
#    dash_style = (
#    (0, 20, -15, 30, 10),
#    (1, 30, 0, 15, 10),
#    (0, 40, 15, 15, 10),
#    (1, 20, 30, 60, 10),
#    )
    for fx in freq:
        
        axvline(x=fx, color=option.olc[0])
        
#        (dd, dl, r, dr, dp) = dash_style[i]
        ax.text(fx, -65, labels[i], withdash=False, rotation=90, verticalalignment='top', horizontalalignment='center')
        i = i + 1
#    sys.exit()



def logThisRun():
    f = open(".plot_function.log", "a")
    a = ''
    for arg in sys.argv:
        a = a + r" %s" % arg
    dt = datetime.datetime.now()
    f.write(dt.ctime())
    f.write("\n")
    f.write(a)
    f.write("\n")
    f.close()


# def makeNestedGridPlot(xmin,xmax,ymin,ymax):
def makeNestedGridPlot(xminxmaxyminymax):
    rect = False
    ax = gca()
    if option.nestedGridLevel >= 0:
        rect = Rectangle((xminxmaxyminymax[0], xminxmaxyminymax[2]), xminxmaxyminymax[1] - xminxmaxyminymax[0], xminxmaxyminymax[3] - xminxmaxyminymax[2], edgecolor='r', facecolor="none")
        ax.add_patch(rect)
    
    r = xminxmaxyminymax
    nestedGridDivisions = option.nestedGridLevel
        
    if option.nestedGridLevel > 0:
        n = nestedGridDivisions - 1
        if option.nestedGridOct:
            nestedGridDivisions = np.power(2, option.nestedGridLevel)
            n = np.power(2, option.nestedGridLevel) - 1
        for l in arange(n):
            dx = [r[0], r[1]]
            y = r[2] + (l + 1.0) * (r[3] - r[2]) / nestedGridDivisions
            dy = [y, y]
            ax.plot(dx, dy, color='r', lw=1)
            x = r[0] + (l + 1.0) * (r[1] - r[0]) / nestedGridDivisions
            dx = [x, x]
            dy = [r[2], r[3]]
            ax.plot(dx, dy, color='r', lw=1)

        
#    if rect!=False:
    
###########################################################################################
###########################################################################################
###########################################################################################
# MAIN PROGRAM
###########################################################################################
###########################################################################################
###########################################################################################
if option.nolog == False:
    logThisRun()


if option.interactive:
    from matplotlib.widgets import Lasso
    from matplotlib.collections import RegularPolyCollection
    from matplotlib import colors as mcolors, path
    from numpy import nonzero
    from numpy.random import rand



formatterX = FuncFormatter(XTicksFormatter)
formatterY = FuncFormatter(YTicksFormatter)
formatterLog2LinZ = FuncFormatter(ZTicksFormatter)










while 1:

    
    
    # if len(args)==0:
    #    if option.st==0 and option.en==0:
    #        print "Too few parameters given"
    #        sys.exit(0)
    print
    print "Loading data"
    print
    inFile = list()
    for i in arange(len(args)):
        plotTypeIdx = i % len(option.plotType)
        globalPlotTypeIdx = plotTypeIdx
    
        print "=================================================== LOADING DATASET %i ====================================================" % i
        print "Loading file: %s" % args[i]
        if option.PR:
            print "this is not implemented yet"
            plotPOVray_file(args[i])
        else:
            if option.fits or option.fileFormat[i % len (option.fileFormat)] == 'fits':
                hdulist = pyfits.open(args[i])
                print hdulist.info()
                tbdata = hdulist[option.hdu[i % len(option.hdu)]].data  # assuming the first extension is a table
                s = array([])
                print "Selected HDU coluns names"
                colNames = hdulist[option.hdu[i % len(option.hdu)]].columns.names
                colFormats = hdulist[option.hdu[i % len(option.hdu)]].columns.formats
                print colNames
                print colFormats
                fidx = 0
                loadedColumns = list()
                numericColNo = -1
                for c in colNames:
                    if colFormats[fidx] == 'E' or colFormats[fidx] == '1E' or colFormats[fidx] == 'D' or colFormats[fidx] == '1D' or colFormats[fidx] == 'F' or colFormats[fidx] == 'F8.0':
                        numericColNo = numericColNo + 1
                        s = np.append(s, array(tbdata.field(c)))
                        loadedColumns.append(c)
                        print "%i) %s --> %i" % (fidx, c, numericColNo)
                    fidx = fidx + 1
                print 'len(loadedColumns)', len(loadedColumns)
                slice = transpose(s.reshape(len(loadedColumns), -1))  # make sure slice is an array, and not scalar in case of a single hdu in file
                if len(option.where) > 0:
                    slice = filterData(slice, option.where[i % len(option.where)], loadedColumns)
    #            if len(slice)==0:
    #                print
    #                print "the selected HDU contains not data, will not load"
    #                print
    #            else:
#                 if option.big or len(option.rows)>0:
                if option.big > 0:
    #                print slice
    #                slice=slice.reshape(-1,len(loadedColumns))
    #                print slice
    #                slice=slice[:,[option.colx,option.coly]]
    #                print slice
    #                slice=slice[0:option.rows[i % len(option.rows)]]
    #                print slice
    #                slice.reshape(10,-1)
    #                print slice
    #                sys.exit(0)
    #                inFile.append(slice)
                    inFile.append(slice[0:option.rows[i % len(option.rows)]][:, [option.colx[i % len(option.colx)], option.coly[i % len(option.coly)]]])  # this is actually stupid. the idea of "big" was to avoid loading all data from file, and not loading all and slicing
                else:
                    inFile.append(slice)
    #            if len(inFile[-1])==0:
    #                print "last filtered dataset has zero length. removing"
    #                inFile.pop()
                
            
            else:
                if option.fileFormat[i % len (option.fileFormat)] == 'hdf5':
                    f = h5py.File(args[i], 'r')
                    hdf5dsetName = hdf5dset[i % len(hdf5dset)]
                    if hdf5dsetName == 'all':
                        dsetnames = getHDF5dsetNames(args[i])
                        for dsn in dsetnames:
                            hdf5data = np.array(f[dsn].value)
                            if len(hdf5data) == 3:
                                slice = f[hdf5dset[i % len(hdf5dset)]].value[:, :, option.hdf5slice]
                            else:
                                slice = hdf5data
                            inFile.append(slice.T)
    
                    else:
                        hdf5data = np.array(f[hdf5dsetName].value)
                        if (len(hdf5data.shape) == 3) and hdf5data.shape[2] == 1:
    #                        print hdf5data.shape
    #                        sys.exit(0)
                            if hdf5data.shape[0] == 2:
                                hdf5data = hdf5data[:, :, 0].T
                            else:
                                hdf5data = hdf5data[:, :, 0]
#                             print hdf5data
                            
                            inFile.append(hdf5data)
                        elif option.plotType[plotTypeIdx] == 'mayaVol' and len(hdf5data) == 3:
                            inFile.append(np.array(hdf5data, dtype='float'))
                            print len(hdf5data)
                        else:
                            if len(hdf5data) == 3:
                                slice = f[hdf5dset[i % len(hdf5dset)]].value[:, :, option.hdf5slice]
                                print "slice"
    #                            print slice
                            else:
                                print "slice2"
                                slice = hdf5data
    #                            print slice
                            inFile.append(slice.T)
                    f.close()
                else:
    
                    if (option.fileFormat[i % len (option.fileFormat)] == 'txt') or ('txtCol' in option.fileFormat[i % len (option.fileFormat)]):
    #                    if option.plotType[plotTypeIdx]=='ts':
    #                        import datetime
    #                        import matplotlib.dates as mdates
    #                        import matplotlib.cbook as cbook
    #
    #                        years    = mdates.YearLocator()   # every year
    #                        months   = mdates.MonthLocator()  # every month
    #                        yearsFmt = mdates.DateFormatter('%Y')
    #                        
    #                        # load a numpy record array from yahoo csv data with fields date,
    #                        # open, close, volume, adj_close from the mpl-data/example directory.
    #                        # The record array stores python datetime.date as an object array in
    #                        # the date column
    #                        datafile = cbook.get_sample_data('goog.npy')
    #                        r = np.load(datafile).view(np.recarray)
                            
                            

                        if option.plotType[plotTypeIdx] == 'barchartH':
                            inFile.append(np.loadtxt(args[i], dtype="string"))
                            
                        if option.plotType[plotTypeIdx] == 'vect' or option.plotType[plotTypeIdx] == 'ts' or option.plotType[plotTypeIdx] == 'fn' or option.plotType[plotTypeIdx] == 'step' or option.plotType[plotTypeIdx] == 'err' or option.plotType[plotTypeIdx] == 'fnshaded' or option.plotType[plotTypeIdx] == 'fillbtwX'  or option.plotType[plotTypeIdx] == 'fillbtwY' or option.plotType[plotTypeIdx] == 'scat' or option.plotType[plotTypeIdx] == 'scat3d' or option.plotType[plotTypeIdx] == 'scatContFill' or option.plotType[plotTypeIdx] == 'scatContFillD' or option.plotType[plotTypeIdx] == 'scatDensContFill' or option.plotType[plotTypeIdx] == 'scatDensCont' or option.plotType[plotTypeIdx] == 'sphere' or option.plotType[plotTypeIdx] == 'hist' or option.plotType[plotTypeIdx] == 'circ' or option.plotType[plotTypeIdx] == 'mayaScat':
                            if option.big:  # or len(option.rows)>0:
                                if args[i] == "stdin":
                                    print "warning: option big cannot read from stdin yet. The result might not be what you wanted."
                                    sys.exit()
                                if 'txtCol' in option.fileFormat[i % len (option.fileFormat)]:
                                    print "warning: format txtCol=X cannot read from stdin yet."
                                    sys.exit()
                                inFile.append(loadDataFromFile(args[i], option.colx[i % len(option.colx)], option.coly[i % len(option.coly)], 
                                                               option.startFrom, option.rows[i % len(option.rows)], option.every[i % len(option.every)], 
                                                               option.bin[i % len(option.bin)], badY=option.badY[i % len(option.badY)]))
                            else:
                                if args[i] == "stdin":
                                    if globalStdInData == '':
                                        globalStdInData = np.loadtxt(sys.stdin)
        #                            inFile.append(np.loadtxt(sys.stdin))
                                    inFile.append(globalStdInData)
                                elif 'UDP=' in args[i]:
                                    from socket import *
                                    sys.path.append(os.environ['OCRA_TOOLKIT_DIR'] + '/scripts/fluxCalibration/')  # FOR EXTRA MODULES LOCATED IN LOCAL DIRECTORY
                                    from pyCPEDScommonFunctions import cpedsPythCommon                                
                                    inFile.append(loadDataFromUDP(args[i]))
                                    
                                else:
#                                     if len(option.rows)==0:
#                                         option.rows=list([1])
#                                     inFile.append(np.loadtxt(args[i]))

                                    inFile.append(loadDataFromFileStd(fname=args[i], startFrom=option.startFrom, rowsCount=option.rows[i % len(option.rows)], loadEvery=option.every[i % len(option.every)], binSamples=option.bin[i % len(option.bin)], 
                                                                      colx=option.colx[i % len(option.colx)],
                                                                      coly=option.coly[i % len(option.coly)],
                                                                      badY=option.badY[i % len(option.badY)]))
    
                                if 'txtCol' in option.fileFormat[i % len (option.fileFormat)]:
                                    txtNcols = int(option.fileFormat[i % len (option.fileFormat)].split('=')[1])
    #                                 print 'Ncols=',txtNcols
                                    inFile[-1] = np.reshape(inFile[-1], (-1, txtNcols))
    #                                 print inFile[-1]
    #                                 sys.exit()
                                    
                                inFile[-1] = filterRows(inFile[-1], i)
                        if option.plotType[plotTypeIdx] == 'scatContFill' or option.plotType[plotTypeIdx] == 'scatContFillD' or option.plotType[plotTypeIdx] == 'scatDensContFill':
                            from matplotlib.tri import Triangulation, UniformTriRefiner
                            import matplotlib.tri as mtri
    #                        import matplotlib.tri as tri
    #                        from matplotlib.tri import UniformTriRefiner
                        if option.plotType[plotTypeIdx] == 'img':
                #            from matplotlib._png import read_png
                            import matplotlib.image as mpimg
                            arr = mpimg.imread(args[i])
                #            print args[i]
                #            fn = get_sample_data(args[i])
                #            arr = read_png(fn)
                #            print arr
                            inFile.append(arr)
                #            sys.exit()
                #            inFile.append(cbook.get_sample_data(args[i]), asfileobj=True)
                #            sys.exit()
                
                #            dpi = rcParams['figure.dpi']
                #            figsize = lena.size[0]/dpi, lena.size[1]/dpi
                            
                #            figure(figsize=figsize)
                #            ax = axes([0,0,1,1], frameon=False)
                #            ax.set_axis_off()
        
        
                    if option.plotType[plotTypeIdx] == 'map':
                        cmd = ''
                        if args[i][-7:] == '-Tn-bin':
                            cmd = 'draw_maps %s -o lonlat --reverse_l --cyclic  --resX %i --resY %i' % (args[i], option.mapPts, option.mapPts)
                            cpedsPythCommon.sayAndExecute("Exporting data from binary file", cmd, 1)
                        if args[i][-5:] == '.fits':
                            if option.fileFormat[i % len (option.fileFormat)] == 'fitsWMAP':
                                cmd = 'draw_maps %s --ft fitsWMAP -o lonlat --reverse_l --cyclic  --resX %i --resY %i' % (args[i], option.mapPts, option.mapPts)
                            if option.fileFormat[i % len (option.fileFormat)] == 'fitsPL':
                                cmd = 'draw_maps %s --ft fitsPL -o lonlat --reverse_l --cyclic  --resX %i --resY %i' % (args[i], option.mapPts, option.mapPts)
                            cpedsPythCommon.sayAndExecute("Exporting data from binary file", cmd, 1)
            
                        inFile.append(array([np.loadtxt(args[i] + '.lonlatT'), np.loadtxt(args[i] + '.lon'), np.loadtxt(args[i] + '.lat')]))
        if option.transpose:
            inFile[-1] = inFile[-1].T
            print inFile[-1]

    
    for i in arange(len(inFile)):
        print i, len(inFile)
        plotTypeIdx = i
        if option.plotType[plotTypeIdx % len(option.plotType)] != 'img' and option.plotType[plotTypeIdx % len(option.plotType)] != 'map' and option.plotType[plotTypeIdx % len(option.plotType)] != 'PR':
            print 'shape is:'
            print inFile[i].shape
            print 'len shape is:'
            print len(inFile[i].shape)
            if len(inFile[i].shape) == 0:
                print 'single number'
                inFile[i] = array([[inFile[i]]])
                
            if len(inFile[i].shape) == 1 and inFile[i].shape[0] == 2:
                print 'single pair of numbers'
                inFile[i] = array([inFile[i]])
                
            if len(inFile[i].shape) == 1 and inFile[i].shape[0] > 2:
                print 'single vector'
                if option.asCol:
                    print "treating this as a column vector"
                    inFile[i] = array([inFile[i]]).transpose()
                else:
                    inFile[i] = array([inFile[i]])
        
    if option.PR:
    #    print inFile[i]
        sys.exit()
        
    print
    print "Generating data"
    print
    for i in arange(len(option.plotCircle)):
        print "Generating circles"
        lbr = array([ toFloat(v) for v in option.plotCircle[i].split(',') ])
        cmd = '%sdraw_maps nofile %s --tl %lf --tb %lf --tv 1 --ts %lf --plot_reg_type emptydot --dont_plot_zero -o lb' % (option.draw_maps_path, option.draw_maps_options, lbr[0], lbr[1], lbr[2])
        os.system(cmd)
        inFile.append(np.loadtxt('lb'))
        os.remove('lb')
        
    if option.mk11line:
        print "Generating 11line"
        xy = np.array([[option.xmin[0], option.xmin[0]], [option.xmax[0], option.xmax[0]]])
        print xy
        inFile.append(xy)
        option.colx.append(0)
        option.coly.append(1)
    #
    # load mask if requested
    #
    if option.loadMask != '':
        maskRanges = loadMask(option.loadMask)
        for r in maskRanges:
            maskMaskedPointsMapRange(r[0], r[1])
    
    
    #
    # generate z-order
    #
    if len(option.zorder) == 0:
        print 'GENERATING Z-ORDER'
        for i in np.arange(len(inFile)):
            option.zorder.append(i)
        
    
    print
    print "Plotting"
    print
    
    
    
    if option.plotType[0] == 'mayaVol' or option.plotType[0] == 'mayaVolCut' or option.plotType[0] == 'mayaScat' or option.plotType[0] == 'mayaSurf' or option.plotType[0] == 'mayaBarChart' or option.plotType[0] == 'mayaContour':
#         from enthought.mayavi import mlab
        from mayavi import mlab
        makeMayaViPlot(inFile)
        
    else:
        fig, ax = makeFunctionPlot(inFile)
    
    
    print "done"
    print "---------- "
    print
    
    if option.continuousPotting == False:
        sys.exit(0)
    
    # if option.makeMovie:
    #    print
    #    print "Generating movie..."
    #    print
    #    cmd="mencoder mf://*.png -mf type=png:fps=%i -ovc lavc -lavcopts  vcodec=mpeg4:mbd=1:vbitrate=2800  -oac copy -o %s.avi" % (option.fps, args[0])
    #    os.system(cmd)
    #    print "done"
    #    print "---------- "
    #    print
