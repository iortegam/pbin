#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltFRAPPE.py
#
# Purpose:
#       The purpose of this program is to plt FRAPPE
#
# Notes:
#   
#
# Version History:
#       Created, August, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
import sys
import glob
import getopt

from scipy import interpolate

import matplotlib.dates as md
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
from itertools import izip
from numpy import *
import myfunctions as mf
from datetime import datetime, timedelta

import dataOutplts as dc
from mpl_toolkits.basemap import Basemap
import srtm
#from scipy.spatial import cKDTree
from geographiclib.geodesic import Geodesic
import sys
from dbfread import DBF
from scipy import stats

import PltClass as pc

                                    #-------------------------#
                                    # Define helper functions #
                                    #-------------------------#

def ckDir(dirName,logFlg=False,exit=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if logFlg: logFlg.error('Directory %s does not exist' % dirName)
        if exit: sys.exit()
        return False
    else:
        return True

def ckFile(fName,logFlg=False,exit=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if logFlg: logFlg.error('Unable to find file: %s' % fName)
        if exit: sys.exit()
        return False
    else:
        return True

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
       imshow(x, cmap=djet)
   """
   
    if type(cmap) == str:
        cmap = plt.get_cmap(cmap)
    colors_i = concatenate((linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)


def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    #setp(bp['fliers'][0], color='blue')
    #setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    #setp(bp['fliers'][2], color='red')
    #setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')


                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main(argv):

    #-------------------------------------------------
    # Initializations for C-130 insitu data
    #-------------------------------------------------
    dataDir     = '/ya4/Campaign/FL0/FRAPPE/C130/'
    
    #-------------------------------------------------
    # Date range to Read C-130 insitu data (ICART FILES)
    #-------------------------------------------------
    iyear      = 2014
    imnth      = 7
    iday       = 25 
    fyear      = 2014
    fmnth      = 8
    fday       = 20

    pltFile    = 'figFRAPPE/C130_FRAPPE.pdf'
    #-------------------------------------------------
    #                 FLAGS
    #-------------------------------------------------
    saveFlg         = True
    
    #--------------FLAGS TO READ DATA-----------------
    RFFlag          = True              #FLIGHT INFORMATION (ALWAYS TRUE TO INTERPOLATE LAT AND LON)
    
    C2H6Flag        = True              #C2H6
    COFlag          = True              #CO
    H2COFlag        = True              #H2CO
    NH3Flag         = True              #H2CO
    HCOOHFlag       = True              #H2CO

    #--------------FLAGS TO PLOT DATA-----------------
    pltC2H6toCO     = True              #Plt Ratio
    pltNH3toCO      = False              #Plt Ratio
    pltH2COtoCO     = False              #Plt Ratio

    pltMap          = True              #ALWAYS TRUE (IF AT LEAST ONE MAP IS CREATED)
   
    pltMapRF        = False             #Map of all RF
    pltMapC2H6      = True             #Map of all C2H6
    pltMapCO        = True             #Map of all CO
    pltMapH2CO      = True             #Map of all h2co
    pltMapNH3       = True             #Map of all NH3
    pltMapHCOOH     = True             #Map of all NH3
    
    pltMapC2H6toCO  = True             #Map of the Ratio
    pltMapNH3toCO   = False             #Map of the Ratio
    pltMapH2COtoCO  = False             #Map of the Ratio


    #---------------------------------------------------
    # LOCATIONS for AIR MASSES
    #---------------------------------------------------
    #locP = [ [-104.0, 40.5] , [-104.6, 39.2], [-106.75, 40.0]]
    #locP = [ [-104.6, 40.42] , [-104.6, 39.4], [-106.4, 40.0]]
    locP = [ [-104.6, 40.42] , [-105.245, 39.4], [-106.6, 39.5]]
    #maxD = 75.0          #Maximum Distance in Km
    maxD = 50.0
    clr  = ['red', 'blue',  'green' ]
    IDP  = ['O&NG/Feedlot', 'Urban', 'Background']
   
    lonB, latB = -105.245, 40.035      #LOCATION OF NCAR FOOTHILLS LAB in Boulder

    #------------------------------------------------------------------------------------------------------
    #
    #                                               START
    #
    #------------------------------------------------------------------------------------------------------

    #---------------------------------------------------
    # INITIALIZE A CLASS FOR MAPS
    #---------------------------------------------------
    if pltMap: mp = pc.MapClass(origin=(40.0,-105.0) , maxAlt=4600, minAlt=500., DistLeft=250000,  DistRight=350000,  DistTop=250000,  DistBottom=200000, saveFlg=saveFlg)

    #-------------------------------------------------
    #                 START
    #-------------------------------------------------
    if saveFlg: pdfsav = PdfPages(pltFile)
    else:       pdfsav = False

    i_date   = dt.date(iyear,imnth,iday)                                                     # Initial Day
    f_date   = dt.date(fyear,fmnth,fday)                                                     # Final Day
    ndays    = (f_date + dt.timedelta(days=1) - i_date).days
    dateList =[i_date + dt.timedelta(days=i) for i in range(0, ndays, 1)]  
    
    #--------------------------------------------
    # Walk through first level of directories and
    # collect directory names for processing
    #--------------------------------------------
    dirlist        = {}
    Data           = {}
    gases          = []
  

    for dd in dateList:
    	# Find year month and day strings
        yrstr   = "{0:02d}".format(dd.year)
        mnthstr = "{0:02d}".format(dd.month)
        daystr  = "{0:02d}".format(dd.day)


        filename = 'FRAPPE-NCAR-LRT-NAV_C130_'+yrstr+mnthstr+daystr+'*.ict'
        icarttfile = glob.glob( dataDir + filename )
        #if not icarttfile: continue
        for k in icarttfile:
            #dirlist.append(k)
            dirlist.setdefault('RF', []).append(k)

            
        if COFlag:

            filename = 'frappe-CO_C130_'+yrstr+mnthstr+daystr+'*.ict'
            icarttfile = glob.glob( dataDir + filename )
            #if not icarttfile: continue
            for k in icarttfile:
                dirlist.setdefault('CO', []).append(k)

            gases.append('CO')

        if C2H6Flag:

            filename = 'frappe-C2H6_C130_'+yrstr+mnthstr+daystr+'*.ict'
            icarttfile = glob.glob( dataDir + filename )
            #if not icarttfile: continue
            for k in icarttfile:
                dirlist.setdefault('C2H6', []).append(k)

            gases.append('C2H6')

        if H2COFlag:

            filename = 'frappe-CH2O_C130_'+yrstr+mnthstr+daystr+'*.ict'
            icarttfile = glob.glob( dataDir + filename )
            #if not icarttfile: continue
            for k in icarttfile:
                dirlist.setdefault('H2CO', []).append(k)

            gases.append('H2CO')

        if NH3Flag:

            filename = 'FRAPPE-NH3_C130_'+yrstr+mnthstr+daystr+'*.ict'
            icarttfile = glob.glob( dataDir + filename )
            #if not icarttfile: continue
            for k in icarttfile:
                dirlist.setdefault('NH3', []).append(k)

            gases.append('NH3')

        if HCOOHFlag:

            filename = 'FRAPPE-PCIMS_C130_'+yrstr+mnthstr+daystr+'*3sec.ict'
            icarttfile = glob.glob( dataDir + filename )
            #if not icarttfile: continue
            for k in icarttfile:
                dirlist.setdefault('HCOOH', []).append(k)

            gases.append('HCOOH')

    for i in dirlist:
        dirlist[i].sort()

    gases = list(set(gases))


    #--------------------------------------------
    # 
    #--------------------------------------------
    



    #--------------------------------------------
    #
    #--------------------------------------------

    print 'Reading C130:'
    altC130       = []
    LatC130       = []
    LonC130       = []
    DateC130      = []
    RF            = []

    for indMain,sngDir in enumerate(dirlist['RF']):
	    ckDir(sngDir)
    	    keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)
            DateC    = dateutc
    	    altC     = data[vnames[9]]   
            LatC     = data[vnames[11]]
            LonC     = data[vnames[14]]

            altC    = np.array(altC, dtype=np.float)
            LatC     = np.array(LatC, dtype=np.float)
            LonC     = np.array(LonC, dtype=np.float)
            DateC    = np.array(DateC)

            #index    = np.where( (altC >= 0.0) & (altC <= 3000.0) & (LatC >= -100.0)  & (LonC >= -3000.0) )[0]
            index    = np.where( (altC >= 0.0) & (LatC >= -100.0)  & (LonC >= -3000.0) )[0]

            altC     = altC[index]
            LatC     = LatC[index]
            LonC     = LonC[index]
            DateC    = DateC[index]

            DateC130.append(DateC)
            altC130.append(altC)
            LatC130.append(LatC)
            LonC130.append(LonC)
            RF.append(indMain + 1)

    altC130=np.array(altC130).flatten()
    LatC130=np.array(LatC130).flatten()
    LonC130=np.array(LonC130).flatten()
    RF=np.array(RF).flatten()
    DateC130=np.array(DateC130).flatten()
    doyC130 = [mf.toYearFraction(d) for d in DateC130]

    if 'CO' in gases:

        print '\nReading CO from C130 in-situ Files:'
        COC130        = []
        DateC130_CO   = []

        for indMain,sngDir in enumerate(dirlist['CO']):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            CO       = data[vnames[1]]   

            CO       = np.array(CO, dtype=np.float)
            DateC    = np.array(DateC)

            index    = np.where( (CO >= 0.0) )[0]

            CO       = CO[index]
            DateC    = DateC[index]

            DateC130_CO.append(DateC)
            COC130.append(CO)

        COC130=np.array(COC130).flatten()
        DateC130_CO=np.array(DateC130_CO).flatten()
        doyC130_CO = [mf.toYearFraction(d) for d in DateC130_CO]

        LatC130_CO_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_CO)]
        LonC130_CO_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_CO)]

        Data['CO'] = {'conc':COC130, 'DT':DateC130_CO, 'doy':doyC130_CO, 'lat':LatC130_CO_int, 'lon':LonC130_CO_int }
    

    if 'C2H6' in gases:
        print '\nReading C2H6 from C130 in-situ Files:'
        C2H6C130        = []
        C2H6C130_e      = []
        DateC130_C2H6   = []

        for indMain,sngDir in enumerate(dirlist['C2H6']):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            C2H6     = data[vnames[3]]
            C2H6_e   = data[vnames[8]]   

            C2H6     = np.array(C2H6, dtype=np.float)/1000.0  #ppt to ppb
            C2H6_e   = np.array(C2H6_e, dtype=np.float)/1000.0  #ppt to ppb
            DateC    = np.array(DateC)

            index    = np.where( (C2H6 >= 0.0) )[0]

            C2H6     = C2H6[index]
            C2H6_e   = C2H6_e[index]
            DateC    = DateC[index]

            DateC130_C2H6.append(DateC)
            C2H6C130.append(C2H6)
            C2H6C130_e.append(C2H6_e)

        C2H6C130      = np.array(C2H6C130).flatten()
        C2H6C130_e    = np.array(C2H6C130_e).flatten()
        DateC130_C2H6 =np.array(DateC130_C2H6).flatten()
        doyC130_C2H6  = [mf.toYearFraction(d) for d in DateC130_C2H6]

        LatC130_C2H6_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_C2H6)]
        LonC130_C2H6_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_C2H6)]

        Data['C2H6'] = {'conc':C2H6C130, 'DT':DateC130_C2H6, 'doy':doyC130_C2H6, 'lat':LatC130_C2H6_int, 'lon':LonC130_C2H6_int }

    if 'H2CO' in gases:
        print '\nReading H2CO from C130 in-situ Files:'
        H2COC130        = []
        DateC130_H2CO   = []

        for indMain,sngDir in enumerate(dirlist['H2CO']):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            H2CO     = data[vnames[3]]   

            H2CO     = np.array(H2CO, dtype=np.float)/1000.0  #ppt to ppb
            DateC    = np.array(DateC)

            index    = np.where( (H2CO >= 0.0) )[0]

            H2CO     = H2CO[index]
            DateC    = DateC[index]

            DateC130_H2CO.append(DateC)
            H2COC130.append(H2CO)

        H2COC130=np.array(H2COC130).flatten()
        DateC130_H2CO=np.array(DateC130_H2CO).flatten()
        doyC130_H2CO = [mf.toYearFraction(d) for d in DateC130_H2CO]

        LatC130_H2CO_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_H2CO)]
        LonC130_H2CO_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_H2CO)]

        Data['H2CO'] = {'conc':H2COC130, 'DT':DateC130_H2CO, 'doy':doyC130_H2CO, 'lat':LatC130_H2CO_int, 'lon':LonC130_H2CO_int }

    if 'NH3' in gases:
        print '\nReading NH3 from C130 in-situ Files:'
        NH3C130        = []
        DateC130_NH3   = []

        for indMain,sngDir in enumerate(dirlist['NH3']):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

            DateC    = dateutc
            NH3     = data[vnames[3]]   

            NH3     = np.array(NH3, dtype=np.float) 
            DateC    = np.array(DateC)

            index    = np.where( (NH3 >= 0.0) )[0]

            NH3     = NH3[index]
            DateC    = DateC[index]

            DateC130_NH3.append(DateC)
            NH3C130.append(NH3)

        NH3C130=np.array(NH3C130).flatten()
        DateC130_NH3=np.array(DateC130_NH3).flatten()
        doyC130_NH3 = [mf.toYearFraction(d) for d in DateC130_NH3]

        LatC130_NH3_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_NH3)]
        LonC130_NH3_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_NH3)]

        Data['NH3'] = {'conc':NH3C130, 'DT':DateC130_NH3, 'doy':doyC130_NH3, 'lat':LatC130_NH3_int, 'lon':LonC130_NH3_int }

    if 'HCOOH' in gases:
        print '\nReading HCOOH from C130 in-situ Files:'
        HCOOHC130        = []
        DateC130_HCOOH   = []

        for indMain,sngDir in enumerate(dirlist['HCOOH']):
            ckDir(sngDir)
            keys, data,dateutc, vnames = mf.read_ICARTT(sngDir)

        
            DateC     = dateutc
            HCOOH     = data[vnames[5]]   

            HCOOH     = np.array(HCOOH, dtype=np.float)/1000.0  #ppt to ppb
            DateC    = np.array(DateC)

            index    = np.where( (HCOOH >= 0.0) )[0]

            HCOOH     = HCOOH[index]
            DateC    = DateC[index]

            DateC130_HCOOH.append(DateC)
            HCOOHC130.append(HCOOH)

        HCOOHC130=np.array(HCOOHC130).flatten()
        DateC130_HCOOH=np.array(DateC130_HCOOH).flatten()
        doyC130_HCOOH = [mf.toYearFraction(d) for d in DateC130_HCOOH]

        LatC130_HCOOH_int  = [np.interp(d, doyC130[i], LatC130[i]) for i, d in enumerate(doyC130_HCOOH)]
        LonC130_HCOOH_int  = [np.interp(d, doyC130[i], LonC130[i]) for i, d in enumerate(doyC130_HCOOH)]

        Data['HCOOH'] = {'conc':HCOOHC130, 'DT':DateC130_HCOOH, 'doy':doyC130_HCOOH, 'lat':LatC130_HCOOH_int, 'lon':LonC130_HCOOH_int }

    #---------------------------------------------------------------------------------------
    #
    #                                    PLOTS
    #
    #---------------------------------------------------------------------------------------

    ngases = len(gases)

    #-------------------------------------------------------
    # Box plot of Concentrations
    #-------------------------------------------------------

    locations    = [ [1,2,3], [5,6,7], [9,10,11], [13,14,15], [17,18, 19]]
    lticks       = [2, 6, 10, 14, 18]

    gases        = ['C2H6', 'NH3', 'CO', 'HCOOH', 'H2CO']
    

    fig,  ax = plt.subplots(figsize=(10, 6))
    gasname_w = []

    for k, g in enumerate(gases):

        y    = []

        lat      = np.asarray(Data[g]['lat'])
        lon      = np.asarray(Data[g]['lon'])
        conc     = np.asarray(Data[g]['conc'])

        lat      = np.concatenate(lat)
        lon      = np.concatenate(lon)
        conc     = np.concatenate(conc)

        for p, loc in enumerate(locP):

            print IDP[p]

            Dist     = []
        
            for i, la in enumerate(lat):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lon[i]), abs(lat[i]) )
                Dist.append(float(c))

            Dist     = np.asarray(Dist)
            inds     = np.where(Dist <= maxD)[0]
            conc_i   = np.asarray(conc[inds]) 

            y.append(conc_i)        

            IQR = stats.iqr(conc_i, rng=(25,75))
            print 'IQR = ', str(IQR)

            PCy = np.percentile(conc_i, [25, 50, 75])
            print 'PCy = ', str(PCy)
           
            Qy = stats.mstats.hdquantiles(conc_i, prob=(0.75,0.9))
            print 'Qy = ', str(Qy)

            me   = np.median(conc_i)

            print np.median(conc_i) + PCy[2] + (IQR)
            print np.std(conc_i)*2.698
            print np.median(conc_i)

        y = np.asarray(y)

        if g.upper() == 'CO':    
            y = y/5.

            maxyi    = [np.amax(yi) for yi in y] 
            ax.text(lticks[k], np.amax(maxyi) +( np.amax(maxyi)*0.05), r'x5', fontsize=18)

        if ( (g.upper() == 'H2CO') or (g.upper() == 'HCOOH') ):  

            y = y*10.

            maxyi    = [np.amax(yi) for yi in y] 
            ax.text(lticks[k], np.amax(maxyi) +( np.amax(maxyi)*0.05), r'/10', fontsize=18)

        meanpointprops = dict(marker='o', markeredgecolor='black', markerfacecolor='white', markersize=6)

        bp = ax.boxplot(y, whis=[5, 95], positions = [locations[k][0], locations[k][1], locations[k][2]], widths = 0.85,  showfliers=True, showmeans=True, meanprops=meanpointprops, patch_artist=True)
        
        maxloc = locations[k][2] + 1.0

        gasname = mf.getgasname(g)
        gasname_w.append(gasname)

        for pi, patch in enumerate(bp['boxes']):
            patch.set_facecolor(clr[pi])
            patch.set( linewidth=2)

        for whisker in bp['whiskers']:
            whisker.set(color='k', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='k', linewidth=2)

        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set(color='k', linewidth=2)

        for fly in bp['fliers']:
            fly.set(markersize=2)
        
    ax.set_xlim((0,maxloc))
    ax.set_ylim((-1,100))
    ax.grid(True,which='y', alpha=0.5)
    #ax.set_yscale("log", nonposy='clip')
    ax.set_xticklabels(gasname_w)
    ax.set_xticks(lticks[0:ngases])
    ax.xaxis.set_tick_params(which='major',labelsize=18)
    ax.yaxis.set_tick_params(which='major',labelsize=18)
    ax.set_ylabel('VMR [ppb]', fontsize = 18)
    #ax.set_ylabel(idtitle, fontsize = 18)

    # draw temporary red and blue lines and use them to create a legend
    hB, = ax.plot([1,1], linewidth=2.0, linestyle='-', color=clr[0])
    hR, = ax.plot([1,1], linewidth=2.0, linestyle='-', color=clr[1])
    hC, = ax.plot([1,1], linewidth=2.0, linestyle='-', color=clr[2])
    #ax.legend((hB, hR, hC),('O&NG', 'Urban', 'Background'), prop={'size':16})
    ax.legend((hB, hR, hC),('O&NG/Feedlot', 'Urban', 'Background'), prop={'size':16}, loc='upper center', bbox_to_anchor=(0.5, 1.12),  fancybox=True, ncol=3) 
    #legend((hB, hR),(seasons))
    hB.set_visible(False)
    hR.set_visible(False)
    hC.set_visible(False)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)      
    else:           
        plt.show(block=False)

    fig.savefig('figFRAPPE/box_conc.pdf', bbox_inches='tight')

    #-------------------------------------------------------
    # Box plot of Ratios
    #-------------------------------------------------------

    locations    = [ [1,2,3], [5,6,7], [9,10,11], [13,14,15]]
    lticks       = [2, 6, 10, 14]

    fig,  ax = plt.subplots(figsize=(10, 6))
    gasname_w = []

    #-------------------------------------------------------
    #
    #-------------------------------------------------------
    lat_x     = np.asarray(Data['CO']['lat'])
    lon_x     = np.asarray(Data['CO']['lon'])
    dt_x      = np.asarray(Data['CO']['DT'])
    doy_x     = np.asarray(Data['CO']['doy'])
    conc_x    = np.asarray(Data['CO']['conc'])

    lat_x      = np.concatenate(lat_x)
    lon_x      = np.concatenate(lon_x)
    conc_x     = np.concatenate(conc_x)
    dt_x       = np.concatenate(dt_x)
    doy_x      = np.concatenate(doy_x)

    #-------------------------------------------------------
    #
    #-------------------------------------------------------
    lat_y     = np.asarray(Data['C2H6']['lat'])
    lon_y     = np.asarray(Data['C2H6']['lon'])
    dt_y      = np.asarray(Data['C2H6']['DT'])
    doy_y     = np.asarray(Data['C2H6']['doy'])
    conc_y    = np.asarray(Data['C2H6']['conc'])

    lat_y      = np.concatenate(lat_y)
    lon_y      = np.concatenate(lon_y)
    conc_y     = np.concatenate(conc_y)
    dt_y       = np.concatenate(dt_y)
    doy_y      = np.concatenate(doy_y)

    gases2     = [g for g in gases if g != 'CO']

    for k, g in enumerate(gases2):

        Ra    = []

        #-------------------------------------------------------
        #
        #-------------------------------------------------------
        lat      = np.asarray(Data[g]['lat'])
        lon      = np.asarray(Data[g]['lon'])
        conc     = np.asarray(Data[g]['conc'])
        doy     = np.asarray(Data[g]['doy'])

        lat      = np.concatenate(lat)
        lon      = np.concatenate(lon)
        conc     = np.concatenate(conc)
        doy      = np.concatenate(doy)

        #-------------------------------------------------------

        for p, loc in enumerate(locP):


            Dist     = []

            # if g.upper() != 'NH3':

            #     print 'HERE - {}'.format(g)
        
            #     for i, la in enumerate(lat_y):
            #         c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lon_y[i]), abs(lat_y[i]) )
            #         Dist.append(float(c))

            #     Dist      = np.asarray(Dist)
            #     inds      = np.where(Dist <= maxD)[0]

            #     doy_y2    = np.asarray(doy_y[inds])
            #     conc_y2   = np.asarray(conc_y[inds])
                
            #     #conc_i    = np.asarray(conc[inds])
            #     #doy_i     = np.asarray(doy[inds])

            #     conc_x_i = interpolate.interp1d(doy_x, conc_x, kind='nearest',   fill_value=-999, bounds_error=False)(doy_y2)
            #     conc_i = interpolate.interp1d(doy, conc, kind='nearest',   fill_value=-999, bounds_error=False)(doy_y2)

            #     #IQR = stats.iqr(conc_i, rng=(25,75))
            #     #print 'IQR = ', str(IQR)

            #     #PCy = np.percentile(conc_i, [25, 50, 75])
            #     #print 'PCy = ', str(PCy)
               
            #     #Qy = stats.mstats.hdquantiles(conc_i, prob=(0.75,0.9))
            #     #print 'Qy = ', str(Qy)

            #     if p == 0: 
            #         max_i  = 0. #np.median(conc_y2) + PCy[2]# + (IQR)
                    
            #     else:      
            #         max_i = 0.

            #     inds   = np.where( (conc_y2 >= max_i) & (conc_x_i > 0.) & (conc_i > 0.) )[0]
            #     Ra.append(conc_i[inds]/conc_x_i[inds])


            # else:

            print 'HERE 2 - {}'.format(g)

            for i, la in enumerate(lat):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lon[i]), abs(lat[i]) )
                Dist.append(float(c))

            Dist      = np.asarray(Dist)
            inds      = np.where(Dist <= maxD)[0]
            
            conc_i    = np.asarray(conc[inds])
            doy_i     = np.asarray(doy[inds])

            conc_x_i  = interpolate.interp1d(doy_x, conc_x, kind='nearest',   fill_value=-999, bounds_error=False)(doy_i)
        
            #IQR = stats.iqr(conc_i, rng=(25,75))
            #print 'IQR = ', str(IQR)

            #PCy = np.percentile(conc_i, [25, 50, 75])
            #print 'PCy = ', str(PCy)
           
            #Qy = stats.mstats.hdquantiles(conc_i, prob=(0.75,0.9))
            #print 'Qy = ', str(Qy)

            if p == 0: 
                max_i  = 0. #np.median(conc_i)# + PCy[2]# + (IQR)
                
            else:      
                max_i = 0.

            inds   = np.where( (conc_i >= max_i) & (conc_x_i > 0.) )[0]
            Ra.append(conc_i[inds]/conc_x_i[inds])


            odr  = mf.orthoregress(conc_x_i[inds], conc_i[inds], xerr=conc_x_i[inds]*0.05, yerr=conc_i[inds]*0.1, InError=False)
            slope       = float(odr[0])
            intercept   = float(odr[1])

            print 'gas = ', str(g)
            print 'loc = ', IDP[p]
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'Median Ra = ', str(np.median(conc_i[inds]/conc_x_i[inds]))
            print 'Mean Ra = ', str(np.mean(conc_i[inds]/conc_x_i[inds]))
            print 'std Ra = ', str(np.std(conc_i[inds]/conc_x_i[inds]))



        Ra = np.asarray(Ra)

        if ( (g.upper() == 'H2CO') or (g.upper() == 'HCOOH') ):  
            Ra      = Ra*10.
            maxRa    = [np.amax(r) for r in Ra] 
            ax.text(lticks[k], np.amax(maxRa) +( np.amax(maxRa)*0.05), r'x10', fontsize=18)


        meanpointprops = dict(marker='o', markeredgecolor='black', markerfacecolor='white', markersize=6)

        bp = ax.boxplot(Ra, whis=1.5, positions = [locations[k][0], locations[k][1], locations[k][2]], widths = 0.85,  showfliers=True, showmeans=True, meanprops=meanpointprops, patch_artist=True)
        
        maxloc = locations[k][2] + 1.0

        gasname = mf.getgasname(g)
        gasname_w.append(gasname)

        for pi, patch in enumerate(bp['boxes']):
            patch.set_facecolor(clr[pi])
            patch.set( linewidth=2)

        for whisker in bp['whiskers']:
            whisker.set(color='k', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='k', linewidth=2)

        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set(color='k', linewidth=2)

        for fly in bp['fliers']:
            fly.set(markersize=2)
            
    ax.set_xlim((0,maxloc))
    ax.set_ylim((-0.05,1.0))
    #ax.set_yscale("log", nonposy='clip')
    ax.grid(True,which='y', alpha=0.35)
    ax.set_xticklabels(gasname_w)
    ax.set_xticks(lticks[0:ngases])
    ax.xaxis.set_tick_params(which='major',labelsize=18)
    ax.yaxis.set_tick_params(which='major',labelsize=18)
    ax.set_ylabel('ER [ppb/ppb]', fontsize = 18)
    #ax.set_ylabel(idtitle, fontsize = 18)

    # draw temporary red and blue lines and use them to create a legend
    hB, = ax.plot([1,1], linewidth=2.0, linestyle='--', color=clr[0])
    hR, = ax.plot([1,1], linewidth=2.0, linestyle='--', color=clr[1])
    hC, = ax.plot([1,1], linewidth=2.0, linestyle='--', color=clr[2])
    ax.legend((hB, hR, hC),('O&NG/Feedlot', 'Urban', 'Background'), prop={'size':16}, loc='upper center', bbox_to_anchor=(0.5, 1.12),  fancybox=True, ncol=3)  
    #legend((hB, hR),(seasons))
    hB.set_visible(False)
    hR.set_visible(False)
    hC.set_visible(False)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)      
    else:           
        plt.show(block=False)

    fig.savefig('figFRAPPE/box_Ra.pdf', bbox_inches='tight')

    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()  

    # -----------

    if pltC2H6toCO:
        
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(C2H6C130)
        y_e      = np.concatenate(C2H6C130_e)
        yDate    = np.concatenate(DateC130_C2H6)
        ydoy     = mf.toYearFraction(yDate)


        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)
        
        y_int     = interpolate.interp1d(ydoy, y, kind='nearest',   fill_value=-999, bounds_error=False)(xdoy)
        y_e_int   = interpolate.interp1d(ydoy, y_e, kind='nearest',  fill_value=-999, bounds_error=False)(xdoy)

        x_int     = interpolate.interp1d(xdoy, x,  kind='nearest',  fill_value=-999, bounds_error=False)(ydoy)

        # inds = []

        # for xd in xdoy:

        #     ind_i = mf.nearestind(xd, ydoy)
        #     inds.append(ind_i)
         
        # x_int   = x[inds]


        la_int    = interpolate.interp1d(ddoy, la, fill_value=-999,  bounds_error=False)(ydoy)
        lo_int    = interpolate.interp1d(ddoy, lo, fill_value=-999,  bounds_error=False)(ydoy)

        #inds       = np.where(y_int > 0. )[0]
        inds      = np.where( (y > 0.) & (x_int > 0.) )[0]

        x_int   = x_int[inds]
        y       = y[inds]
        y_e     = y_e[inds]
        lo_int  = lo_int[inds]
        la_int  = la_int[inds]

        # Dist     = []

        # for i, k in enumerate(x):
        #     c = mf.haversine(abs(locP[-1][0]), abs(locP[-1][1]), abs(lo_int[i]), abs(la_int[i]) )
        #     Dist.append(float(c))

        # Dist = np.asarray(Dist)
        # inds = np.where(Dist <= maxD)[0]

        # xbkg = np.nanmean(np.asarray(x[inds]))
        # ybkg = np.nanmean(np.asarray(y_int[inds]))

        # print ybkg
        # print xbkg

        Ra        = np.true_divide(y, x_int)


        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        for p, loc in enumerate(locP):

            Dist     = []

            for i, k in enumerate(y):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1   = np.asarray(x_int[inds]) 
            y1   = np.asarray(y[inds]) 
            y1_e = np.asarray(y_e[inds]) 

            # x1bkg =  np.nanmin(x1)
            # y1bkg =  np.nanmin(y1)

            # x1 = x1 - x1bkg
            # y1 = y1 - y1bkg
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' C2H6 to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)
            print 'std_err =', str(std_err)
            print 'std_err =', str(std_err*np.sqrt(len(x1)))

            print 'mean C2H6 = ', str(np.mean(y1))
            print 'median C2H6 = ', str(np.median(y1))

            print 'mean CO = ', str(np.mean(x1))
            print 'median CO = ', str(np.median(x1))




            print 'Ratio = {} +/- {}'.format(np.mean(Ra[inds]), np.std(Ra[inds]) )
            print 'Ratio = {} +/- {}'.format(np.median(Ra[inds]), np.std(Ra[inds]) )

            odr  = mf.orthoregress(x1, y1, xerr=x1*0.05, yerr=y1_e, InError=False)
            slope       = float(odr[0])
            intercept   = float(odr[1])

            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)

            #slope_e     = float(odrErr[0])
            #intercept_e = float(odrErr[1])

            #print '\nSlope: {0:.2f} +/- {1:.2f}'.format(float(odr[0]), float(odrErr[0]))
            #print 'Intercept = {0:.3f} +/- {1:.3f}'.format(float(odr[1]), float(odrErr[1]))


            xx   = range(0, 500, 2)
            xx   = np.asarray(xx)
            Fit  = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
        
        ax.set_ylabel('C$_2$H$_6$ VMR [ppbv]', fontsize=18)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=18)
        ax.set_ylim(0, 62.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.legend(prop={'size':18})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapC2H6toCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.01, zmax=0.15, ztitle='Ratio C$_2$H$_6$ to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_C2H6toCO.pdf')

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()     

    if pltNH3toCO:
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(NH3C130)
        yDate    = np.concatenate(DateC130_NH3)
        ydoy     = mf.toYearFraction(yDate)

        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)

        y_int     = interpolate.interp1d(ydoy, y,  bounds_error=False)(xdoy)
        la_int    = interpolate.interp1d(ddoy, la,  bounds_error=False)(xdoy)
        lo_int    = interpolate.interp1d(ddoy, lo,   bounds_error=False)(xdoy)

        Ra        = np.true_divide(y_int, x)

        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        #ax.plot(x, y_int, '.k', markersize=1, linewidth=0, alpha=0.15)
        #ax.scatter(x, y_int, facecolors='white', s=50, color='gray', alpha=0.15)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y_int)
        print '\nCorrelation of all data' + ' NH3 to CO'
        print 'Slope = ', str(slope)
        print 'intercept = ', str(intercept)
        print 'r_value = ', str(r_value)
        print 'std_err =', str(std_err)

        for p, loc in enumerate(locP):
            Dist     = []

            for i, k in enumerate(x):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1 = np.asarray(x[inds])
            y1 = np.asarray(y_int[inds])
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' NH3 to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)

            xx = range(0, 500, 2)
            xx = np.asarray(xx)
            Fit = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
            
        
        ax.set_ylabel('NH$_3$ VMR [ppbv]', fontsize=18)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=18)
        ax.set_ylim(0, 40.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=18)
        ax.legend(prop={'size':18})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapNH3toCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.01, zmax=0.075, ztitle='Ratio NH$_3$ to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_NH3toCO.pdf')
    
    if pltH2COtoCO:
        x        = np.concatenate(COC130)
        xDate    = np.concatenate(DateC130_CO)
        xdoy     = mf.toYearFraction(xDate)

        y        = np.concatenate(H2COC130)
        yDate    = np.concatenate(DateC130_H2CO)
        ydoy     = mf.toYearFraction(yDate)

        la       = np.concatenate(LatC130)
        lo       = np.concatenate(LonC130)
        d        = np.concatenate(DateC130)
        ddoy     = mf.toYearFraction(d)

        y_int     = interpolate.interp1d(ydoy, y,  bounds_error=False)(xdoy)
        la_int    = interpolate.interp1d(ddoy, la,  bounds_error=False)(xdoy)
        lo_int    = interpolate.interp1d(ddoy, lo,   bounds_error=False)(xdoy)

        Ra        = np.true_divide(y_int, x)

        fig, ax = plt.subplots(figsize=(10,6), sharex=True)

        #ax.plot(x, y_int, '.k', markersize=1, linewidth=0, alpha=0.15)
        #ax.scatter(x, y_int, facecolors='white', s=50, color='gray', alpha=0.15)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y_int)
        print '\nCorrelation of all data' + ' H2CO to CO'
        print 'Slope = ', str(slope)
        print 'intercept = ', str(intercept)
        print 'r_value = ', str(r_value)
        print 'std_err =', str(std_err)

        for p, loc in enumerate(locP):
            Dist     = []

            for i, k in enumerate(x):
                c = mf.haversine(abs(loc[0]), abs(loc[1]), abs(lo_int[i]), abs(la_int[i]) )
                Dist.append(float(c))

            Dist = np.asarray(Dist)
            inds = np.where(Dist <= maxD)[0]

            x1 = np.asarray(x[inds])
            y1 = np.asarray(y_int[inds])
            
            ax.scatter(x1, y1, facecolors='white', s=70, color=clr[p], label=IDP[p])

            slope, intercept, r_value, p_value, std_err = stats.linregress(x1, y1)
            print '\nCorrelation of '+ IDP[p] + ' H2CO to CO'
            print 'Slope = ', str(slope)
            print 'intercept = ', str(intercept)
            print 'r_value = ', str(r_value)
            print 'std_err =', str(std_err)

            xx = range(0, 500, 2)
            xx = np.asarray(xx)
            Fit = slope*np.asarray(xx) + (intercept)#*math.exp(-1.6/7.4)
            ax.plot(xx, Fit, color=clr[p],  linewidth=2.0, linestyle='--')
            ax.fill_between(xx,Fit - Fit*0.2 ,Fit + Fit*0.2 ,alpha=0.25, color=clr[p])
            
        
        ax.set_ylabel('H$_2$CO VMR [ppbv]', fontsize=16)
        ax.set_xlabel('CO VMR [ppbv]', fontsize=16)
        ax.set_ylim(0, 7.0)
        ax.set_xlim(50, 400.0)
        ax.grid(True)
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.legend(prop={'size':12})
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.95)

        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if pltMapH2COtoCO:
            mp.pltMapZ( LonData=lo_int, LatData=la_int, zData=Ra, zmin=0.005, zmax=0.04, ztitle='Ratio H$_2$CO to CO [ppb]', 
                        LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_H2COtoCO.pdf')
   
    #------------------------------------------------------------------
    #                     PLOT MAP OF C2H6  
    #------------------------------------------------------------------
    if pltMapC2H6:

        C2H6C130    = np.concatenate(C2H6C130)
        LonC130_int = np.concatenate(LonC130_C2H6_int)
        LatC130_int = np.concatenate(LatC130_C2H6_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=C2H6C130, zmin=0.5, zmax=10, ztitle='C$_2$H$_6$ [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_C2H6.pdf')
    
    #------------------------------------------------------------------
    #                     PLOT MAP OF CO  
    #------------------------------------------------------------------
    if pltMapCO:

        COC130      = np.concatenate(COC130)
        LonC130_int = np.concatenate(LonC130_CO_int)
        LatC130_int = np.concatenate(LatC130_CO_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=COC130, zmin=50, zmax=150, ztitle='CO [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_CO.pdf')

    #------------------------------------------------------------------
    #                     PLOT MAP OF H2CO  
    #------------------------------------------------------------------
    if pltMapH2CO:

        H2COC130      = np.concatenate(H2COC130)
        LonC130_int = np.concatenate(LonC130_H2CO_int)
        LatC130_int = np.concatenate(LatC130_H2CO_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=H2COC130, zmin=0.1, zmax=4.0, ztitle='H$_2$CO [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_H2CO.pdf')

    #------------------------------------------------------------------
    #                     PLOT MAP OF NH3  
    #------------------------------------------------------------------
    if pltMapNH3:

        NH3C130      = np.concatenate(NH3C130)
        LonC130_int = np.concatenate(LonC130_NH3_int)
        LatC130_int = np.concatenate(LatC130_NH3_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=NH3C130, zmin=0.05, zmax=15.0, ztitle='NH$_3$ [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_NH3.pdf')

    #------------------------------------------------------------------
    #                     PLOT MAP OF NH3  
    #------------------------------------------------------------------
    if pltMapHCOOH:

        HCOOHC130   = np.concatenate(HCOOHC130)
        LonC130_int = np.concatenate(LonC130_HCOOH_int)
        LatC130_int = np.concatenate(LatC130_HCOOH_int)

        mp.pltMapZ( LonData=LonC130_int, LatData=LatC130_int, zData=HCOOHC130, zmin=0.05, zmax=3.0, ztitle='HCOOH [ppb]', 
                    LatID=latB, LonID=lonB, SaveFile='figFRAPPE/C130_FRAPPE_MAP_HCOOH.pdf')
    
    #------------------------------------------------------------------
    #                     PLOT MAP WITH RF TRACKS  
    #------------------------------------------------------------------
    if pltMapRF:
        mp.pltMapRF(LonData=LonC130, LatData=LatC130, zData=RF, ztitle='Research Flight', 
                    LatID=latB, LonID=lonB, locP=locP, SaveFile='figFRAPPE/C130_FRAPPE_MAP_RF.pdf')

      #--------------------------------
      # Pause so user can look at plots
      #--------------------------------
    if saveFlg: 
        pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])