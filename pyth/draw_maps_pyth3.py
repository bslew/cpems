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
m1.drawparallels(arange(-90.,90.,20.),labels=[1,0,0,0])
m1.drawmeridians(arange(0.,420.,30.),labels=[0,0,0,1])
#title('Moll. proj:'+ sys.argv[1])
#title('>3$\sigma$ detections composite map: skewness')
#title('>3$\sigma$ detections composite map: kurtosis')
#title('Mollweide projection')
colorbar();

savefig(sys.argv[1]+'-moll.png',dpi=200); #savefig(sys.argv[1]+'-moll.eps')

# north polar projection.
m2 = Basemap(llcrnrlon=-45.,llcrnrlat=10.,urcrnrlon=135.,urcrnrlat=10., projection='stere', lat_0=90.,lon_0=0.,lat_ts=90.)
#            resolution='c',area_thresh=10000.,

fig2 = m2.createfigure()
x, y = m2(lons, lats)
cs = m2.contourf(x,y,map,color_num,cmap=cm.jet)

delat = 15.
circles = arange(0.,90.+delat,delat).tolist()+arange(-delat,-90.-delat,-delat).tolist()
m2.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
delon = 45.
meridians = arange(-180,180,delon)
m2.drawmeridians(meridians,labels=[1,1,1,1])
#title('North Stereographic Proj.:'+ sys.argv[1] )
#title('North Stereographic Proj.:' )
#colorbar();

savefig(sys.argv[1]+'-north_stere.png',dpi=200); #savefig(sys.argv[1]+'-north_stere.eps')

# south polar projection.
m3 = Basemap(llcrnrlon=-45.,llcrnrlat=-10.,urcrnrlon=135.,urcrnrlat=-10., projection='stere', lat_0=-90.,lon_0=0.,lat_ts=-90.)
#            resolution='c',area_thresh=10000.,

fig3 = m3.createfigure()
x, y = m3(lons, lats)
cs = m3.contourf(x,y,map,color_num,cmap=cm.jet)

delat = 15.
circles = arange(0.,90.+delat,delat).tolist()+arange(-delat,-90.-delat,-delat).tolist()
m3.drawparallels(circles,labels=[1,1,1,1])
# draw meridians
delon = 45.
meridians = arange(-180,180,delon)
m3.drawmeridians(meridians,labels=[1,1,1,1])
#title('South Stereographic Proj.:'+ sys.argv[1])
#title('South Stereographic Proj.:')

savefig(sys.argv[1]+'-south_stere.png',dpi=200); #savefig(sys.argv[1]+'-south_stere.eps')

#show()
