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
print lons
print lons

#map,lons = addcyclic(map,lons)
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
#map,lons = shiftgrid(180.5,map,lons,start=False)
lons, lats = meshgrid(lons, lats)
print lons
print
print lats
#lons[0]=-180
#color_num=50
#sys.exit()

# Set up a colormap:
palette = matplotlib.cm.jet
#palette.set_over('w', 1.0)
#palette.set_under('k', 1.0)


#mapM=map

figure(figsize=(16,8), dpi=70)

m=Basemap(projection='moll',lon_0=180,lat_0=0) #.5*(lons[0]+lons[-1])
#delta=(maxv-minv)/color_num
#levels = arange(minv-1*delta, maxv+1*delta,delta ); print(levels)
levels=50
#levels = arange(minv, maxv, delta ); print(levels)


x, y = m(lons, lats)

contourf(x,y,map, levels, cmap=matplotlib.pyplot.cm.jet )
 
m.drawmapboundary(color='k', linewidth=1.0, ax=None)

subplots_adjust(left=0.01,right=.8,top=0.99,bottom=0.01,wspace=0,hspace=0)
cax = axes([0.85, 0.15,  0.02, 0.7])
colorbar(cax=cax)
ylabel('T [K]')
show()

#sys.exit()

