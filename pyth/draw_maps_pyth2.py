#!/usr/bin/env python
import sys # get the command line arguments
from matplotlib.toolkits.basemap import Basemap, shiftgrid
from pylab import *

#mapfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lonlatT'
#lonfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lon'
#latfile = '/home/blew/programy/Mscs/WMAPstuff/SKregstat/sm1.0/3sigma_chisq_comb-K.lat'

mapfile = sys.argv[1]+'.lonlatT'
lonfile = sys.argv[1]+'.lon'
latfile = sys.argv[1]+'.lat'
color_num = int(sys.argv[2])

# load data
map=array(load(mapfile))
lons=array(load(lonfile))
lats=array(load(latfile))
map,lons = shiftgrid(180.,map,lons,start=False)

# Mollweide projection
m1 = Basemap(projection='moll',lon_0=0.5*(lons[0]+lons[-1]))
lons, lats = meshgrid(lons, lats)
fig1 = m1.createfigure()
x, y = m1(lons, lats)
cs = m1.contourf(x,y,map,color_num,cmap=cm.jet)
#m1.drawparallels(arange(-90.,90.,20.),labels=[1,0,0,0])
#m1.drawmeridians(arange(0.,420.,30.),labels=[0,0,0,1])
#title('Moll. proj:'+ sys.argv[1])
#title('>3$\sigma$ detections composite map: skewness')
#title('>3$\sigma$ detections composite map: kurtosis')
#title('Mollweide projection')
#colorbar(tickfmt='%.1lE');

savefig(sys.argv[1]+'-moll.png',dpi=300); #savefig(sys.argv[1]+'-moll.eps')

show()
