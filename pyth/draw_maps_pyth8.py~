#!/usr/bin/env python
import sys # get the command line arguments
from pylab import *
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
#from matplotlib.basemap import Basemap, shiftgrid
#import matplotlib.numerix.ma as M
import matplotlib.colors as colors
#from matplotlib import *
import asciidata

#mapfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lonlatT'
#lonfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lon'
#latfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lat'

mapfile = sys.argv[1]+'.lonlatT'
lonfile = sys.argv[1]+'.lon'
latfile = sys.argv[1]+'.lat'
color_num = int(sys.argv[2])
minv = float(sys.argv[3])
maxv = float(sys.argv[4])
dpires = int(sys.argv[5])
bgcolor = sys.argv[6]
fgcolor = sys.argv[7]
mytitle = sys.argv[8]
mycolorbar = sys.argv[9]
mymeridians = sys.argv[10]
plotNS = sys.argv[11]
DeltaLatitude = float(sys.argv[12])
DeltaLongitude = float(sys.argv[13])
labels_format = sys.argv[14]
colorbar_label = sys.argv[15]
colorbar_fontsize = sys.argv[16]
meridians_fontsize = sys.argv[17]
plot_text = sys.argv[18]

#Cldata=sys.argv[19]
#Clsim=sys.argv[20]

# load data
map=array(mlab.load(mapfile))


#zdata=asciidata.open(mapfile);
#minv=zdata[0][0]
#maxv=zdata[0][0]
#for i in range(zdata.nrows):
#    for j in range(zdata.ncols):
#        if zdata[j][i] < minv:
#            minv=zdata[j][i]
#        if zdata[j][i] > maxv:
#            maxv=zdata[j][i]


#minv=min(tlist)
#maxv=max(tlist)
#minv=-1.791916E-04
#maxv=1.878548E-04
print(minv,maxv)


lons=array(mlab.load(lonfile))
lats=array(mlab.load(latfile))
map,lons = shiftgrid(180,map,lons,start=False)
#mapM,lons = shiftgrid(180,mapM,lons,start=False)


# Set up a colormap:
palette = matplotlib.cm.jet
palette.set_over('w', 1.0)
palette.set_under('k', 1.0)
#palette.set_bad('r', 1.0)
#palette.set_bad(alpha = 0.0)
#palette.set_under(alpha = 0.0)
# Alternatively, we could use
#palette.set_bad(alpha = 0.0)
# to make the bad region transparent.  This is the default.
# If you comment out all the palette.set* lines, you will see
# all the defaults; under and over will be colored with the
# first and last colors in the palette, respectively.
#Zm = ma.masked_where(Z > 1.2, Z)
# By setting vmin and vmax in the norm, we establish the
# range to which the regular palette color scale is applied.
# Anything above that range is colored based on palette.set_over, etc.
#im = imshow(Zm, interpolation='bilinear',
#    cmap=palette,
#    norm = colors.normalize(vmin = -1.0, vmax = 1.0, clip = False),
#    origin='lower', extent=[-3,3,-3,3])
#title('Green=low, Red=high, Blue=bad')
# colorbar(im) # colorbar does not yet show the under and over colors.


# make masked array
#map=M.array(map)
#threshold=-1e20
#mapM=M.masked_where(map < minv,map)
#mapM=M.masked_where(map < -4.1,map)
mapM=map
#mapM = M.masked_where(map > 1.2, map)
#map = where(map < minv,1.e10,map)
#mapM = M.masked_values(map, 1.e10)


# make topodat a masked array, masking values lower than sea level.
#topodat = where(topodat < 0.,1.e10,topodat)
#topodatm = ma.masked_values(topodat, 1.e10)
#palette = cm.YlOrRd
#palette.set_bad('aqua', 1.0)
# plot image over map with imshow.
#im = m.imshow(topodatm,palette,norm=colors.normalize(vmin=0.0,vmax=3000.0,clip=False))






figure(figsize=(15,10), dpi=dpires)
subplots_adjust(left=0.2,right=.98,top=0.97,bottom=0.11,wspace=0,hspace=0.08)
subplot(1,2,1)

#fig = m1.createfigure()
# Mollweide projection
#m1 = 
m1=Basemap(projection='moll',lon_0=0,lat_0=0) #.5*(lons[0]+lons[-1])
lons, lats = meshgrid(lons, lats)
delta=(maxv-minv)/color_num
levels = arange(minv-1*delta, maxv+1*delta,delta ); print(levels)
#levels = arange(minv, maxv, delta ); print(levels)

x, y = m1(lons, lats)

#fig=figure(figsize=(12,6), dpi=dpires)
#fig=figure(figsize=(4,2), dpi=dpires)
#print dpires





#cs = m1.contourf(x,y,mapM,color_num,cmap=cm.jet)
#cs = m1.contourf(x,y,mapM,color_num,cmap=palette, norm = matplotlib.colors.normalize(vmin=minv, vmax=maxv, clip=False) )
#cs = m1.contourf(x,y,mapM,colors=None, cmap=palette, norm = matplotlib.colors.normalize(vmin=minv, vmax=maxv, clip=False) )
#cs = m1.contourf(x,y,mapM, color_num, cmap=palette, norm = matplotlib.colors.normalize(vmin=minv, vmax=maxv, clip=False), extend='min' )

#cs = m1.contourf(x,y,mapM, levels, cmap=palette, norm = matplotlib.colors.normalize(minv, maxv, clip=False) )
contourf(x,y,mapM, levels, cmap=palette )

#cs = m1.imshow(mapM, cmap=palette, norm = matplotlib.colors.normalize(vmin=minv, vmax=maxv, clip=False) )

m1.drawmapboundary(color='k', linewidth=1.0, ax=None)
#cs = m1.contourf(x,y,mapM, cmap=palette, norm = colors.normalize(vmin=-2, vmax=2, clip=False) )
#im = imshow(mapM, interpolation='bilinear',    cmap=palette,    norm = colors.normalize(vmin = -1.0, vmax = 1.0, clip = False),    origin='lower', extent=[-3,3,-3,3])
#cs = imshow(mapM, cmap=palette, norm = colors.normalize(vmin=minv, vmax=maxv, clip=False) )






















































if mymeridians == '1' :
    m1.drawparallels(arange(-90.,90.,DeltaLatitude),labels=[1,0,0,0],fontsize=meridians_fontsize)
    m1.drawmeridians(arange(0.,420.,DeltaLongitude),labels=[0,0,0,1])
#   m1.drawparallels(arange(-90.,90.,20.),labels=[1,0,0,0])
#   m1.drawmeridians(arange(0.,420.,30.),labels=[0,0,0,1])

#title('Moll. proj:'+ sys.argv[1])
#title('>3$\sigma$ detections composite map: skewness')
#title('>3$\sigma$ detections composite map: kurtosis')
#rc('font',c=fgcolor)
#rc('lines', lw=2, c=fgcolor)
if mytitle != 'notitle' :
    title(mytitle,color=fgcolor, fontsize=meridians_fontsize)

if mycolorbar == 'moll' :
    #    colorbar(tickfmt='%.1lf',)
    #setp(a,pad=0.05)
    #clabel(cs, levels[1::2],  # label every second level
    #       inline=1,
    #       fmt='%1.1f',
    #       fontsize=14)
    if mymeridians == '1' :
#        cax = axes([0.15, 0.05, 0.7, 0.03])
        cax = axes([0.15, 0.07, 0.7, 0.05])
    else :
        cax = axes([0.15, 0.11, 0.7, 0.06])
        
    #cax = axes([0.0, 0.0, 0.02, 0.25])
    #setp(gca(), xticks=[-4.1e-5,0,4.1e-5])
    #setp(gca(), xticks=(minv,0,maxv))
    a=colorbar(mappable=None, cax=cax, orientation='horizontal', drawedges=False)
    #setp(a, xticks=(minv,0,maxv))
    #xlabel('K')
    #setp(a, xticks=[-4.1e-5,0,4.1e-5])
    #    setp(a, labels=['-4.1e-5','4.1e-5'])
    #xlabels = a.get_xticklabels()
    #setp(xlabels, fontsize=10, rotation=0)
    if mymeridians == '1' :
        subplots_adjust(left=0.05,right=1,top=1,bottom=0.15,wspace=0,hspace=0)

    else :
        #subplots_adjust(left=0,right=1,top=1,bottom=0.15,wspace=0,hspace=0) # normal
        subplots_adjust(left=0.05,right=.95,top=1,bottom=0.2,wspace=0,hspace=0) # experimental
    #zeros='%1.1E' % (0)
    if labels_format == 'E':
        minvs='%1.1E' % (minv)
        maxvs='%1.1E %s' % (maxv, colorbar_label)
        middlev='%1.1E' % ((minv+maxv)/2)

    if labels_format == 'f':
        minvs='%1.1f' % (minv)
        maxvs='$%1.1f [ %s ]$' % (maxv, 'n_\sigma')
        middlev='%1.1f' % ((minv+maxv)/2)
    
    
    text(minv, -0.5, minvs, horizontalalignment='center', verticalalignment='top', fontsize=colorbar_fontsize,rotation=0)
    text(maxv, -0.5, maxvs, horizontalalignment='center', verticalalignment='top', fontsize=colorbar_fontsize,rotation=0, fontname='Sans')
    #text(0, -0.5, '0', horizontalalignment='center', verticalalignment='top', fontsize=14,rotation=0)
    text((minv+maxv)/2, -0.5, middlev, horizontalalignment='center', verticalalignment='top', fontsize=colorbar_fontsize,rotation=0)
    #text(maxv+0.4*maxv, -0.5, '[K]', horizontalalignment='center', verticalalignment='top', fontsize=14,rotation=0)
    #xticks((minv,0,maxv), ('1','2','3'),fontsize=10)
    #zlabels = cs.clabels()
    #setp(xlabels, rotation=45)
    #setp(gca(),'xticklabels', [])

    #setp(a,aspect=10.0,rotation=45)
    

if plot_text == 'true' :
    axes(ax)  # make the original axes current again
    
    plot_text_data = asciidata.open(sys.argv[1]+'.from_map_taken_vals'); 
    txtl=[]; txtb=[]; txt=[];
    for i in range(plot_text_data.nrows):
        txtl.append(plot_text_data[0][i]);
        txtb.append(plot_text_data[1][i]);
        txt.append(plot_text_data[3][i]);
        xpt,ypt = m1(-txtl[i],txtb[i])         
        #m1.plot([xpt],[ypt],'bo')
        print(xpt, ypt,txt[i])
        text(xpt,ypt,txt[i],fontsize=meridians_fontsize)
        

#show()
#exit



show()
