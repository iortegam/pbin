#! /usr/bin/python2.7
import srtm
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#from scipy.spatial import cKDTree
from geographiclib.geodesic import Geodesic
import sys
from dbfread import DBF
from collections import namedtuple
from matplotlib.backends.backend_pdf import PdfPages

#------------------------------------------------------------------
#             THE FOLLOWING INFORMATION IS TO CREEATE THE MAP
#------------------------------------------------------------------
origin = (40.0,-105.0)   #CENTER OF THE MAP
max_altitude = 4600               #MAXIMUM ALTITUDE TO CONSIDER IN THE TERRAIN

azimuth, distance = (270, 85000)  #DISTANCE ON THE LEFT (in meters) 
geolft = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

azimuth, distance = (90,128000)   #RIGHT
georgt = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

azimuth, distance = (0,88000)     #TOP
geotop = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

azimuth, distance = (180,88000)   #BOTTOM
geobot = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

domain_size = (500, 500)          #DOMAIN SIZE (roughly the RESOLUTION OF THE MAP)

#------------------------------------------------------------------
#                THE FOLLOWING IS ADDITIONAL INFORMATION 
#------------------------------------------------------------------
lonB, latB = -105.245, 40.035      #LOCATION OF NCAR FOOTHILLS LAB in Boulder
lonD, latD = -104.9903, 39.7292    #LOCATION OF Denver, CO

#------------------------------------------------------------------
#                      FILE WITH WELL LOCATIONS
#------------------------------------------------------------------
LatWell = []
LonWell = []
FacType = []
FacStat = []

#DATA FROM HERE: http://cogcc.state.co.us/data2.html#/downloads
table = DBF('/data/iortega/results/fl0/WELLS_SHP2/Wells.dbf', load=True)
#table = DBF('/data/iortega/results/fl0/OIL_AND_GAS_LOCATIONS_SHP/Oil_and_Gas_Locations.dbf', load=True)

for record in table:
	LatWell.append(record['Latitude'])
	LonWell.append(record['Longitude'])
	FacType.append(record['Facil_Type'])
	FacStat.append(record['Facil_Stat'])

LatWell = np.asarray(LatWell)
LonWell = np.asarray(LonWell)
FacType = np.asarray(FacType)
FacStat = np.asarray(FacStat)

print 'Number of Wells in CO: '+ str(len(LatWell))

inds = np.where( (LatWell > float(geobot['lat2'])) & (LatWell <  float(geotop['lat2'])) & (LonWell > float(geolft['lon2'])) & (LonWell < float(georgt['lon2']))    )[0]

print 'Number of Wells in the map: '+ str(len(inds))

LatWell = LatWell[inds]
LonWell = LonWell[inds]
FacType = FacType[inds]
FacStat = FacStat[inds]

indsAC = np.where( FacStat == 'AC'  )[0]   #Active 
indsPR = np.where( FacStat == 'PR'  )[0]   #Producing
indsDG = np.where( FacStat == 'DG'  )[0]   #Drilling
indsWO = np.where( FacStat == 'WO'  )[0]   #Drilling

LatWellAC = LatWell[indsAC]
LonWellAC = LonWell[indsAC]

LatWellPR = LatWell[indsPR]
LonWellPR = LonWell[indsPR]

LatWellDG = LatWell[indsDG]
LonWellDG = LonWell[indsDG]

LatWellWO = LatWell[indsWO]
LonWellWO = LonWell[indsWO]

#------------------------------------------------------------------
#                      START MAP
#------------------------------------------------------------------
domain=[domain_size,
        (geobot['lat2'],geotop['lat2']),
        (geolft['lon2'],georgt['lon2']),
        max_altitude]

geodata = srtm.get_data()

elev_array = geodata.get_image(*domain, mode='array')

lats=np.linspace(domain[1][0],domain[1][1],domain[0][1])
lons=np.linspace(domain[2][0],domain[2][1],domain[0][0])

latn=np.min(lats)
latx=np.max(lats)
lonn=np.min(lons)
lonx=np.max(lons)


fig,ax = plt.subplots(figsize=(10.5,7))

m= Basemap(llcrnrlon=geolft['lon2'], llcrnrlat=geobot['lat2'], urcrnrlon=georgt['lon2'], urcrnrlat=geotop['lat2'],
             resolution='i', projection='lcc', lat_0 = origin[0], lon_0 = origin[1], ax=ax) # projection='tmerc'

#cmap = plt.get_cmap('gist_ncar')
#cmap = plt.get_cmap('jet')
#cmap = plt.get_cmap('nipy_spectral')
#cmap = plt.get_cmap('terrain_r')
cmap = plt.get_cmap('terrain')

x, y = m(lons, lats)
p = m.pcolormesh(x,y,elev_array, cmap=cmap, vmin=500, vmax=domain[3])
p.cmap.set_over(cmap(.0))
p.cmap.set_under('w')
#ax.set_ylabel('Latitude [$^{\circ}$]', fontsize=14)
#ax.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)    
ax.tick_params(labelsize=14)
#ax.set_title('Location of FTS and O&NG', fontsize=16)
cbar = fig.colorbar(p, orientation='vertical', fraction=0.05, pad=0.035)
cbar.set_label('Altitude [m]', fontsize=14)
#cbar.set_ticklabels(labelsize=12)

xx, yy = m(LonWell, LatWell)
m.plot(xx, yy, c='k', marker='.', markersize=2, linestyle='None', alpha=0.35, label = 'Drilling and/or abandoned') 

xxx, yyy = m(LonWellPR, LatWellPR)
m.plot(xxx, yyy, c='m', marker='.', markersize=2, linestyle='None', alpha=0.5, label = 'Producing') 
#ax.plot(LonWellDG,LatWellDG, c='c', marker='.', markersize=2, alpha=0.5, label = 'Drilling') 
#ax.plot(LonWellWO,LonWellWO,'b.', markersize=2, alpha=0.75, label = 'Waiting on completion') 

x, y = m(lonB, latB)
x2, y2 = m(lonB-0.095, latB-0.085)
ax.plot(x,y,'r^', markersize=15)
plt.text(x2,y2,'NCAR', color='r', fontsize=16)

xD, yD = m(lonD, latD)
x2D, y2D = m(lonD-0.095, latD-0.085)
ax.plot(xD,yD,'k^', markersize=15)
plt.text(x2D,y2D,'Denver', color='k', fontsize=16)

ax.legend(loc=3, fontsize=14, markerscale=6)

#m.drawmapboundary(fill_color='aqua')
#m.fillcontinents(color='#cc9955', lake_color='aqua')
m.drawcounties()
# draw parallels.
parallels = np.arange(30.,50, 0.5)
m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14, linewidth=0.0, color='gray')
# draw meridians
meridians = np.arange(180.,360., 0.5)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14, linewidth=0.0, color='gray')
m.drawmapscale(-105.6, 40.7 , origin[1], origin[0], 50, barstyle='fancy')

plt.show(block=False)
#plt.savefig('/data/iortega/Manuscripts/Tropgases/fig//Map.pdf', bbox_inches='tight', dpi=80)
plt.savefig('/data/iortega/pbin/tropGases/fig/Map.pdf', bbox_inches='tight', dpi=80)

user_input = raw_input('Press any key to exit >>> ')
sys.exit()
