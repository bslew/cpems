#!/usr/bin/env python
import sys # get the command line arguments
from pylab import *
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
import numpy as np
#from matplotlib.basemap import Basemap, shiftgrid
#import matplotlib.numerix.ma as M
import matplotlib.pyplot as plt
import matplotlib.mlab
#import matplotlib.colors as colors
#from matplotlib import *
#import asciidata
#import MA

mapfile = sys.argv[1]+'.lonlatT'
lonfile = sys.argv[1]+'.lon'
latfile = sys.argv[1]+'.lat'


# load data
map=np.array(mlab.load(mapfile))




#minv=MA.minimum(map)
#maxv=MA.maximum(map)
#print(minv,maxv)

lons=np.array(mlab.load(lonfile))
lats=np.array(mlab.load(latfile))
#tmp=lons
#tmp=lons[0]
#lons=lons-tmp
#lons[size(lons)-1]=360
#print lons
#lons[0]=0;
#print lons
map,lons = addcyclic(map,lons)
#lons[-1]=360;
#print lons
#shift=[]
#dl=(lons[1]-lons[0])/2.0
#for phi in arange(0,256):
#    shift.append(lons[0]-dl)
#print size(shift)
#print size(lons)
#sys.exit()
#lons = lons -shift
#print lons
#sys.exit()
map,lons = shiftgrid(180,map,lons,start=False)
#print lons
#sys.exit()

lons, lats = meshgrid(lons, lats)
print
#print lats
#lons[0]=-180
#color_num=50

# Set up a colormap:
palette = matplotlib.cm.jet
#palette.set_over('w', 1.0)
#palette.set_under('k', 1.0)


#mapM=map

figure(figsize=(16,8), dpi=70)

#lon_0 must be in the middle of the grid, so if the grid starts with
# -180.3 and ends at 179.7 it should be -0.3 or else...
lon0=lons[0][0]+180.
print lon0
m=Basemap(projection='moll',lon_0=lon0,lat_0=0) #.5*(lons[0]+lons[-1])
#delta=(maxv-minv)/color_num
#levels = arange(minv-1*delta, maxv+1*delta,delta ); print(levels)
levels=50
#levels = arange(minv, maxv, delta ); print(levels)

ax=gca()
x, y = m(lons, lats)

contourf(x,y,map, levels, cmap=matplotlib.pyplot.cm.jet )

parallels = np.arange(-60.,90,30.)
m.drawparallels(parallels,labels=[1,0,0,0])
meridians = np.arange(0.,360.,30.)
m.drawmeridians(meridians)
 
m.drawmapboundary(color='k', linewidth=1.0, ax=None)


#DeltaLatitude=45
#meridians_fontsize=12
#m.drawparallels(arange(-90.,90.,DeltaLatitude),labels=[1,0,0,0],fontsize=meridians_fontsize)

subplots_adjust(left=0.04,right=.8,top=0.99,bottom=0.01,wspace=0,hspace=0)
cax = axes([0.85, 0.15,  0.02, 0.7])
colorbar(cax=cax)
ylabel('T [K]')
#savefig(sys.argv[1]+'.png')
savefig(sys.argv[1]+'.eps')
#show()

#sys.exit()
