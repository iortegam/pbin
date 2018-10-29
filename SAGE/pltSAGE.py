#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltSAGE-FTS.py
#
# Purpose:
#  
#
# Notes:
#   
#
# Version History:
#       Created, Sep, 2017  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#

import sys
import os
import getopt
import Class_SAGE as cs
import time
from matplotlib.backends.backend_pdf import PdfPages
import myfunctions as mf
import ClassFTS as cfts
from itertools import izip

import datetime as dt
import numpy as np
from numpy import *

import matplotlib


import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
from cycler import cycler

from scipy import linspace, polyval, polyfit, sqrt, stats, randn


from shapely.geometry.polygon import Polygon as Polygon2
from matplotlib.patches import Polygon as PatchPolygon
from shapely.geometry.polygon import Polygon

from mpl_toolkits.basemap import Basemap

from scipy import interpolate


#from descartes import PolygonPatch

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

def main():


    #-----------------------------------------------------------------------------------------
    #                             Initializations for FTIR
    #-----------------------------------------------------------------------------------------
    readftsFlg         = True

    saveFlg           = True 
    
    loc                = 'fl0'
    gasName            = ['o3']             
    
    if loc.lower() == 'mlo':
        ver                = ['Current']          # Name of retrieval version to process
        ctlF               = ['sfit4.ctl'] 
    elif loc.lower() == 'fl0':
        ver                = ['Current_WP']          # Name of retrieval version to process
        ctlF               = ['sfit4.ctl'] 
    else:
        print 'Error!!'
        exit()
    #----------------------
    # First Level Retrieval Directory
    #----------------------
    retDir             = '/data1/ebaumer/'+loc.lower()

    pltFile           = '/data/iortega/pbin/SAGE/SAGEIII-FTS-'+loc.upper()+'.pdf'

    

    #------
    # Flags
    #------
    errorFlg           = True                   # Flag to process error data
    fltrFlg            = True                   # Flag to filter the data
    byYrFlg            = False                  # Flag to create plots for each individual year in date range
    szaFlg             = True                   # Flag to filter based on min and max SZA
    dofFlg             = True                   # Flag to filter based on min DOFs
    pcNegFlg           = True                  # Flag to filter profiles with negative partial columns
    tcNegFlg           = True                  # Flagsag to filter profiles with negative total columns
    tcMMFlg            = False                  # Flag to filter based on min and max total column amount
    cnvrgFlg           = True                   # Flag to filter profiles that did not converge
    rmsFlg             = True                   # Flag to filter based on max RMS
    chiFlg             = False                  # Flag to filter based on max CHI_2_Y
    mnthFlg            = False                  # Flag to filter based on 

    mnths              = [6, 7, 8]                # Months to filter on (these are the months to include data)
    maxRMS             = [1.5]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0]                      # Min DOFs for filtering

    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    pColsFlg           = True
    pCols              = [ [16., 24.], [24., 34.], [16., 34.] ]         #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2017
    imnth              = 6
    iday               = 1
    fyear              = 2018
    fmnth              = 12
    

    fday               = 30
    
    if loc.lower() == 'mlo':
        latfts             = 19.4
        lonfts             = -155.57   # 204.43

        thick    = [14.,     11.5,     9.3,     8.2,     7.23,    6.37 ,   5.62   , 4.93   , 4.3,
                    3.74,    3.24 ,   2.8 ,    2.43 ,   2.11 ,   1.86 ,   1.68  ,  1.55,    1.48,
                    1.44,    1.41 ,   1.37 ,   1.34 ,   1.3 ,    1.26 ,   1.23  ,  1.19 ,   1.16,
                    1.12,    1.08 ,   1.05 ,   1.01 ,   0.98 ,   0.942 ,  0.9061 , 0.8801,  0.8442,
                    0.8182 , 0.7823  ,0.7463 , 0.7204 , 0.6844]

        mid      = [113.,     100.25,    89.85,    81.1,     73.385,   66.585,   60.59,    55.315,
                    50.7,     46.68 ,   43.19  ,  40.17 ,   37.555 ,  35.285 ,  33.3 ,    31.53,
                   29.915 ,  28.4 ,    26.94,    25.515,   24.125 ,  22.77,    21.45 ,   20.17,
                   18.925 ,  17.715 ,  16.54,    15.4  ,   14.3 ,    13.235  , 12.205,   11.21,
                   10.249,    9.3249 ,  8.4318 ,  7.5697,   6.7385,   5.9382 ,  5.174  ,  4.4406,
                   3.7382 ] 

    elif loc.lower() == 'fl0':
        latfts             = 40.4
        lonfts             = -105.24   # 204.43
    
        thick = [14. ,    11.5,     9.3,     8.2,     7.23,    6.37,    5.62,    4.93,    4.3,
                 3.74,    3.24 ,   2.8 ,    2.43,    2.11 ,   1.86  ,  1.68,    1.55,    1.48,
                1.44,    1.41,    1.37,    1.34 ,   1.3  ,   1.26 ,   1.23 ,   1.19 ,   1.16,
                1.12,    1.08 ,   1.05 ,   1.01,    0.98 ,   0.9414,  0.9043,  0.8772 , 0.8401,
                0.813 ,  0.7759 , 0.7388,  0.7117 , 0.6746 , 0.6474,  0.6103 , 0.5732 ]

        mid   = [113.,     100.25 ,   89.85  ,  81.1,     73.385,   66.585 ,  60.59  ,  55.315,
                50.7  ,   46.68  ,  43.19  ,  40.17  ,  37.555 ,  35.285  , 33.3  ,   31.53,
                29.915 ,  28.4  ,   26.94 ,   25.515 ,  24.125 ,  22.77 ,   21.45  ,  20.17,
                18.925 ,  17.715 ,  16.54  ,  15.4  ,   14.3 ,    13.235 ,  12.205  , 11.21,
                10.2493,   9.3264,   8.4356 ,  7.5769 ,  6.7504 ,  5.9559 ,  5.1986 ,  4.4734,
                3.7803 ,  3.1193 ,  2.4904 ,  1.8986]

    else:
        print 'Error:'
        exit()

    thick    = [float(t)*1000.*100. for t in thick] 
    thick    = np.asarray(thick)

    mid      = [float(t) for t in mid] 
    mid      = np.asarray(mid)


    #-----------------------------------------------------------------------------------------
    #                             Initializations for SAGE III/ISS
    #-----------------------------------------------------------------------------------------

    readSAGEFlg         = True

    DirSAGE3_ISS        = '/data1/Campaign/Satellite/SAGEIII-ISS'    
    
    saveSAGEFlg         = False
    pltFilESAGE3        = '/data/iortega/pbin/SAGE/SAGE3_ISSS.pdf'

    pltsage             = False


                                    #----------------------------#
                                    #                            #
                                    #        --- START ---       #
                                    #                            #
                                    #----------------------------#

    if saveFlg: pdfsav = PdfPages(pltFile)
    

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------
    if readftsFlg:
   
        fts = cfts.FTSClass( gasName, retDir,  ctlF, ver, iyear,imnth, iday, fyear, fmnth, fday)

        fts.ReadFTS(fltrFlg=fltrFlg, sclfct=sclfct, sclname=sclfctName, mnthFltr= mnths, mnthFltFlg=mnthFlg,
          errFlg=errorFlg, minSZA=minSZA, maxSZA=maxSZA, maxRMS=maxRMS, minTC=minTC, maxTC= maxTC,
          minDOF=minDOF, maxCHI=maxCHI, dofFlg=dofFlg, rmsFlg=rmsFlg, tcFlg=tcNegFlg,
          pcFlg=pcNegFlg, szaFlg=szaFlg,cnvrgFlg=cnvrgFlg, chiFlg=chiFlg,tcMMflg=tcMMFlg, pColsFlg=False, pCols=pCols)
   
    
    #-------------------------
    #    -- SAGE III --
    #-------------------------
    if pltsage:  
        ckDir(os.path.dirname(os.path.realpath(pltFilESAGE3)),exit=True)

    if readSAGEFlg:
        ckDir(DirSAGE3_ISS, exit=True)

        sage = cs.SAGE3Cls(DirSAGE3_ISS, saveFlg=saveSAGEFlg, outFname=pltFilESAGE3)
        sage.ReadVariables()

        if pltsage: sage.pltSAGE3()
        if saveSAGEFlg: sage.closeFig()

    #-------------------------

    #-------------------------
    GasVer = []

    for (g,v) in izip(gasName, ver):
        GasVer.append(g+'_'+v)

    for gv in GasVer:

        #-------------------------
        # Find Distances close to FTS
        #-------------------------
        latsage       = sage.Lat
        lonsage       = sage.Lon

        #distance      = np.asarray([mf.haversine(lonfts, latfts, lonsage[i],  latsage[i]) for i in range(len(lonsage))])
        #indsDist      = np.where(distance < 1000.)[0]

        if loc.lower() == 'mlo' : indsDist          = np.where((latsage >= 14.) & (latsage <= 26.) & (lonsage <= -140.) & (lonsage >= -166.) )[0]
        #elif loc.lower() == 'fl0' : indsDist        = np.where((latsage >= 35.) & (latsage <= 45.) & (lonsage <= -85.) & (lonsage >= -125.) )[0]
        elif loc.lower() == 'fl0' : indsDist        = np.where((latsage >= 35.) & (latsage <= 45.) & (lonsage <= -95.) & (lonsage >= -115.) )[0]

        else: 
            print 'Error!!'
            exit()
        

        indsDist          = np.asarray(indsDist)


        prfsage       = sage.O3_prf_conc[indsDist]
        dtsage        = sage.dateTime[indsDist]
        altsage       = sage.altitude[indsDist]
        latsage2      = sage.Lat[indsDist]
        lonsage2      = sage.Lon[indsDist]
        datesage      = [ dt.date(singDate.year, singDate.month, singDate.day) for singDate in dtsage  ]
        datesage      = np.asarray(datesage)

        print 'Total Number of SAGE profiles < max distance = ' +str(len(dtsage))+'\n'


        #-------------------------
        # FIND Same days
        #-------------------------
        listdsage     = [ dt.date(singDate.year, singDate.month, singDate.day) for singDate in dtsage  ]
        doysage       = mf.toYearFraction(listdsage)


        listdfts      = [ dt.date(singDate.year, singDate.month, singDate.day) for singDate in fts.dates[gv]  ]
        doyfts        = mf.toYearFraction(listdfts)
  
        intrsctVals   = np.intersect1d(doysage, doyfts, assume_unique=False)

        indssage      = np.nonzero( np.in1d( doysage, intrsctVals, assume_unique=False ) )[0]
        indsfts       = np.nonzero( np.in1d( doyfts, intrsctVals, assume_unique=False ) )[0]

        indssage      = np.asarray(indssage)
        indsfts       = np.asarray(indsfts)

        #-------------------------
        #
        #-------------------------
        prfsage        = prfsage[indssage]
        dtsage         = dtsage[indssage]
        altsage        = altsage[indssage]
        latsage3       = latsage2[indssage]
        lonsage3       = lonsage2[indssage]


        prftfts       = fts.rPrfMol[gv][indsfts, :] /thick#; prftfts  = prftfts[:,indsalt ] /thick[indsalt]
        aprftfts      = fts.aPrfMol[gv][indsfts, :] /thick#; aprftfts = aprftfts[:,indsalt ] /thick[indsalt]
        datesfts      = fts.dates[gv][indsfts]
        altfts        = mid #fts.alt[gv]
        avkSCFfts     = fts.avkSCF[gv][indsfts, :, :]#; avkSCFfts = avkSCFfts[:, indsalt, :]; avkSCFfts = avkSCFfts[:, :, indsalt]

        Airmassfts    = fts.Airmass[gv][indsfts, :]

        avkSCFftsMean = np.mean(avkSCFfts, axis=0)

        indsalt       = np.where((fts.alt[gv] > 10.) & (fts.alt[gv] < 50.) )[0]
        
        prftfts2      = prftfts[:,indsalt ] 
        aprftfts2     = aprftfts[:,indsalt ]
        avkSCFfts2    = avkSCFfts[:, indsalt, :]; avkSCFfts2 = avkSCFfts2[:, :, indsalt]
        altfts2       = altfts[indsalt]

        avkSCFfts2     = np.mean(avkSCFfts2, axis=0)

        datesftslist   = np.asarray([dt.date(d.year,d.month,d.day) for d in datesfts])

      
        #------
        # Daily
        #------
        dailyVals         = mf.dailyAvg(prftfts, datesfts, dateAxis=0, meanAxis=0)
        prftfts_daily     = dailyVals['dailyAvg']
        dates_daily       = dailyVals['dates']
        prftfts_std_daily = dailyVals['std']

        prftfts_daily3    = prftfts_daily

        dailyVals         = mf.dailyAvg(aprftfts, datesfts, dateAxis=0, meanAxis=0)
        aprftfts_daily     = dailyVals['dailyAvg']

        dailyVals         = mf.dailyAvg(prftfts2, datesfts, dateAxis=0, meanAxis=0)
        prftfts_daily2     = dailyVals['dailyAvg']
        dates_daily2       = dailyVals['dates']
        prftfts_std_daily2 = dailyVals['std']

        dailyVals          = mf.dailyAvg(aprftfts2, datesfts, dateAxis=0, meanAxis=0)
        aprftfts_daily2     = dailyVals['dailyAvg']

        print 'Total Number of coincident dates between FTS and SAGE = ' +str(len(dates_daily))+'\n'

        
        prfsage_inte   = np.zeros( (prftfts_daily.shape) )
        prfsage_res    = np.zeros( (prftfts_daily.shape) )
        prfsage_resApr = np.zeros( (prftfts_daily.shape) )
        prfsage_res_rel    = np.zeros( (prftfts_daily.shape) )
        prfsage_resApr_rel = np.zeros( (prftfts_daily.shape) )

        #-------------------------
        #
        #-------------------------
        prfsage_inte2   = np.zeros( (prftfts_daily.shape) )
        prfsage_smth   = np.zeros( (prftfts_daily.shape) )

        prfsage_res2    = np.zeros( (prftfts_daily.shape) )
        prfsage_res_rel2 = np.zeros( (prftfts_daily.shape) )


        for i in range(len(dates_daily)):  

            indsday = np.where(datesftslist == dates_daily[i])[0]
            AKdaily = np.mean(avkSCFfts[indsday, :, :], axis=0)

            #print 'shape', AKdaily.shape

            prfsage_i = prfsage[i,:]
            altsage_i = altsage[i,:]

            indsnan     = np.argwhere(np.isnan(prfsage_i))
            if len(indsnan) >=1:
                prfsage_i  = np.delete(prfsage_i, indsnan)
                altsage_i  = np.delete(altsage_i, indsnan)

            
            #prfsage_interp     = interpolate.interp1d(altsage_i, prfsage_i,fill_value='extrapolate', bounds_error=False, kind='cubic')(altfts2 )
            prfsage_interp        = interpolate.interp1d(altsage_i, prfsage_i, bounds_error=False, kind='linear')(altfts )
            prfsage_inte[i,:]     = prfsage_interp

            prfsage_res[i, :]     = (prftfts_daily[i, :] - prfsage_interp)
            prfsage_resApr[i, :]  = (aprftfts_daily[i, :] - prfsage_interp)

            prfsage_res_rel[i, :]     = (prftfts_daily[i, :] - prfsage_interp)/ prftfts_daily[i, :] *100.
            prfsage_resApr_rel[i, :]  = (aprftfts_daily[i, :] - prfsage_interp)/ prftfts_daily[i, :]*100.

            #-------------------------
            # Partial Columns
            #-------------------------
            prfsage_interp2       = interpolate.interp1d(altsage_i, prfsage_i,fill_value =( prftfts_daily2[i, -1],prftfts_daily2[i, 0]) , bounds_error=False, kind='cubic')(altfts)
            prfsage_inte2[i, :]  = prfsage_interp2

            indsnan          = np.argwhere(np.isnan(prfsage_interp))

            prfsage_smth[i, :]   = aprftfts_daily[i, :] + np.dot(AKdaily, (prfsage_interp2 -  aprftfts_daily[i, :]))

            prfsage_res2[i, :]     = (prftfts_daily[i, :] - prfsage_smth[i, :])
            prfsage_res_rel2[i, :]     = (prftfts_daily[i, :] - prfsage_smth[i, :])/ prftfts_daily[i, :] *100.

            prfsage_smth[i, indsnan]  = np.nan
            prfsage_res2[i, indsnan]  = np.nan
            prfsage_res_rel2[i, indsnan]  = np.nan
            prftfts_daily3[i, indsnan] = np.nan


        
        #-------------------------
        # Partial Columns
        #-------------------------


        #-------------------------
        #
        #-------------------------

        # fig, ax  = plt.subplots(figsize=(6,6))

        # ax.plot(lonsage, latsage, 'k.' , markersize=1)
        # ax.plot(lonsage2, latsage2, 'k.', color='blue', markersize=5)
        # ax.plot(lonfts, latfts, 'r^', markersize=5)
        #     #ax.plot(O3_prf_conc[i,:]/1e13,altitude[i,:], 'k.')
        
        # #ax.fill_betweenx(alt,Prfmean-prfSTD,Prfmean+prfSTD,alpha=0.5,color='0.75')  
        # ax.set_title('Track - SAGE III/ISS - {} - {}'.format(dtsage[0], dtsage[-1]), fontsize=14)
        # ax.set_ylabel('Latitude', fontsize=14)
        # ax.set_xlabel('Longitude', fontsize=14) 
        # ax.set_xlim(xmin=-200, xmax=200) 
        # ax.set_ylim(ymin=-80, ymax=80)    
        # ax.grid(True,which='both')
        # ax.tick_params(labelsize=14)

        
        # fig.subplots_adjust(left=0.12, bottom=0.1,top=0.94, right=0.95)

        # plt.show(block=False)

        #-------------------------
        #
        #-------------------------

        fig, ax = plt.subplots(figsize=(10,6))

        m = Basemap(projection='cyl', resolution='c',
            llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180,  lat_0=45, lon_0=-20)

        xfts, yfts = m(lonfts, latfts)
        ax.plot(xfts, yfts, 'r^', markersize=10)

        m.drawcoastlines(color='black')
        m.drawcountries(color='lightgrey')

        m.drawmeridians(np.arange(0,360,60), linewidth=.2, labels=[1,0,0,1], labelstyle='+/-', color='grey' ) 
        m.drawparallels(np.arange(-90,90,30), linewidth=.2, labels=[1,0,0,1], labelstyle='+/-', color='grey')

        x, y = m(lonsage, latsage)
        ax.plot(x,y,'.k', markersize=5)

        ax.set_title('Track - SAGE III/ISS - {} - {}'.format(datesage[0], datesage[-1]), fontsize=14)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)

        

        #-------------------------
        #
        #-------------------------

        fig, ax = plt.subplots(figsize=(10,6))

        if loc.lower() == 'mlo':
  
            m = Basemap(llcrnrlat=14,urcrnrlat=26, llcrnrlon=-166,urcrnrlon=-147, rsphere=(6378137.00,6356752.3142),
                resolution='l',area_thresh=1000.,projection='lcc', lat_1=20,lon_0=-156)

        elif loc.lower() == 'fl0':

            m = Basemap(llcrnrlat=35,urcrnrlat=45, llcrnrlon=-115,urcrnrlon=-95, rsphere=(6378137.00,6356752.3142),
                resolution='l',area_thresh=1000.,projection='lcc', lat_1=40,lon_0=-105)

        xfts, yfts = m(lonfts, latfts)
        
        m.drawcoastlines(color='black')
        m.drawcountries(color='lightgrey')

        m.drawstates(color='gray')

        m.drawparallels(np.arange(10., 40., 2), labels=[1,0,0,0], alpha=0.25)
        m.drawmeridians(np.arange(-180., 181., 2), labels=[0,0,0,1],alpha=0.25)

        ax.plot(xfts, yfts, 'r^', markersize=10, label = 'MLO')

        x, y = m(lonsage, latsage)
        ax.plot(x,y,'.k', markersize=10)

        #x, y = m(lonsage2, latsage2)
        #ax.plot(x,y,'.b', markersize=10, label='< Distance')

        x, y = m(lonsage3, latsage3)
        ax.plot(x,y,'.r', markersize=10. , label='Coincident Dates')

        ax.set_title('Track - Zoom In - SAGE III/ISS - {} - {}'.format(datesage[0], datesage[-1]), fontsize=14)
    
        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)

        #-------------------------
        #
        #-------------------------
        clmap         = 'jet'
        npanels      = int(math.ceil(len(dates_daily)/3.0))

        #-------------------------
        fig = plt.figure(figsize=(14,14))
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.15, hspace=0.125)
    
        for i in range(len(dates_daily)):
            #----------------------------
            #
            #----------------------------
            ax = plt.Subplot(fig, outer_grid[i])

            p1, =  ax.plot(prftfts_daily[i, :], altfts, linewidth=2.5, label='ret-FTIR', color='b', zorder=1)
            ax.scatter(prftfts_daily[i, :], altfts,facecolors='white', s=10, color='b', zorder=2)
            ax.text(0.8, 0.92, dates_daily[i], va='center',transform=ax.transAxes,fontsize=9)
            #ax.errorbar(dates_mnth[g], vmr_mnth[g], fmt='o', yerr=vmr_std_mnth[g] ,markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)

            p2, = ax.plot(aprftfts_daily[i, :], altfts, linewidth=2.5, label='apr-FTIR', color='gray', zorder=1)

            
            p3,  = ax.plot(prfsage[i, :], altsage[i,:], color='green', linewidth=2.5, zorder=3, label='SAGEIII-ISS')
            ax.scatter(prfsage[i, :], altsage[i,:], color='green', s=10, facecolors='white', zorder=4)


            p4,  = ax.plot(prfsage_inte[i, :], altfts, color='red', linewidth=2.5, zorder=3, label='SAGEIII-ISS - interpolated')
            ax.scatter(prfsage_inte[i, :], altfts, color='red', s=10, facecolors='white', zorder=4)

            p5,  = ax.plot(prfsage_smth[i, :], altfts, color='orange', linewidth=2.5, zorder=3, label='SAGEIII-ISS - smoothed')
            ax.scatter(prfsage_smth[i, :], altfts, color='orange', s=10, facecolors='white', zorder=4)
            
            # Create a legend for the first line.
            if i == 0: 
                first_legend  = fig.legend(handles=[p1, p2, p3, p4, p5], prop={'size':11},loc='center left',bbox_to_anchor=(0.7,0.15))
               
            
            ax.grid(True, alpha=0.35)
            ax.tick_params(labelsize=14)
            ax.grid(True,which='both', alpha=0.35)  
            #ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)
            ax.set_xlim(xmin=0)
            ax.set_ylim(ymin=10, ymax=50)

            fig.add_subplot(ax)

        all_axes  = fig.get_axes()

        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)


        fig.text(0.02, 0.5, 'Altitude [km]', fontsize=18, va='center', rotation='vertical')
        fig.text(0.5, 0.02, 'Concentration [molec/cm$^{3}$]', fontsize=18, ha='center', rotation='horizontal')

        fig.subplots_adjust(bottom=0.075,top = 0.95, left=0.095, right=0.95)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)


        #-------------------------
        fig = plt.figure(figsize=(14,14))
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.15, hspace=0.125)
    
        for i in range(len(dates_daily)):
            #----------------------------
            #
            #----------------------------
            ax = plt.Subplot(fig, outer_grid[i])

            p1, =  ax.plot(prfsage_res[i, :], altfts, linewidth=2.5, label='FTIR - SAGEIII interpolated', color='red', zorder=1)
            ax.scatter(prfsage_res[i, :], altfts,facecolors='white', s=10, color='red', zorder=2)

            ax.text(0.8, 0.92, dates_daily[i], va='center',transform=ax.transAxes,fontsize=9)
            #ax.errorbar(dates_mnth[g], vmr_mnth[g], fmt='o', yerr=vmr_std_mnth[g] ,markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)

            p2, =  ax.plot(prfsage_res2[i, :], altfts, linewidth=2.5, label='FTIR - SAGEIII smoothed', color='orange', zorder=1)
            ax.scatter(prfsage_res2[i, :], altfts,facecolors='white', s=10, color='orange', zorder=2)

            #ax.axvline(x = 0., color='k', linestyle='--')

            # Create a legend for the first line.
            if i == 0: 
                first_legend  = fig.legend(handles=[p1, p2], prop={'size':11},loc='center left',bbox_to_anchor=(0.7,0.2))
               
            
            ax.grid(True, alpha=0.35)
            ax.tick_params(labelsize=14)
            ax.grid(True,which='both', alpha=0.35)  
            #ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)
            ax.set_xlim(-1.5e12, 1.5e12)
            ax.set_ylim(ymin=10, ymax=50)

            fig.add_subplot(ax)

        all_axes  = fig.get_axes()

        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)


        fig.text(0.02, 0.5, 'Altitude [km]', fontsize=18, va='center', rotation='vertical')
        fig.text(0.5, 0.02, 'Difference (FTIR - SAGE III) [molec/cm$^{3}$]', fontsize=18, ha='center', rotation='horizontal')

        fig.subplots_adjust(bottom=0.075,top = 0.95, left=0.095, right=0.95)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
        #-------------------------
        #-------------------------
       
        fig = plt.figure(figsize=(14,14))
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.15, hspace=0.125)
    
        for i in range(len(dates_daily)):
            #----------------------------
            #
            #----------------------------
            ax = plt.Subplot(fig, outer_grid[i])

            p1, =  ax.plot(prfsage_res_rel[i, :], altfts, linewidth=2.5, label='FTIR - SAGEIII interpolated', color='red', zorder=1)
            ax.scatter(prfsage_res_rel[i, :], altfts,facecolors='white', s=10, color='red', zorder=2)
            ax.text(0.8, 0.92, dates_daily[i], va='center',transform=ax.transAxes,fontsize=9)
            #ax.errorbar(dates_mnth[g], vmr_mnth[g], fmt='o', yerr=vmr_std_mnth[g] ,markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)

            p2, =  ax.plot(prfsage_res_rel2[i, :], altfts, linewidth=2.5, label='FTIR - SAGEIII smoothed', color='orange', zorder=1)
            ax.scatter(prfsage_res_rel2[i, :], altfts,facecolors='white', s=10, color='orange', zorder=2)

            #p2, = ax.plot(prfsage_resApr_rel[i, :], altfts, linewidth=2.5, label='aFTIR', color='gray', zorder=1)

            ax.axvline(x = 0., color='k', linestyle='--')

            # Create a legend for the first line.
            if i == 0: 
                first_legend  = fig.legend(handles=[p1, p2], prop={'size':12},loc='center left',bbox_to_anchor=(0.7,0.2))
               
            
            ax.grid(True, alpha=0.35)
            ax.tick_params(labelsize=14)
            ax.grid(True,which='both', alpha=0.35)  
            #ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)
            ax.set_xlim(-150, 150)
            ax.set_ylim(ymin=10, ymax=50)

            fig.add_subplot(ax)

        all_axes  = fig.get_axes()

        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)


        fig.text(0.02, 0.5, 'Altitude [km]', fontsize=18, va='center', rotation='vertical')
        fig.text(0.5, 0.02, 'Relative Difference (FTIR - SAGEIII / FTIR) [%]', fontsize=18, ha='center', rotation='horizontal')

        fig.subplots_adjust(bottom=0.075,top = 0.95, left=0.095, right=0.95)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
        #-------------------------

        if len(datesfts) > 1:

            fig,(ax1,ax2) = plt.subplots(1,2,sharey=True, figsize=(10,7))
            
            prf_relMean = np.nanmean(prfsage_res_rel, axis=0)
            prf_relStd  = np.nanmean(prfsage_res_rel, axis=0)

            prf_relMean2 = np.nanmean(prfsage_res_rel2, axis=0)
            prf_relStd2  = np.nanmean(prfsage_res_rel2, axis=0)

            #ax1.plot(prf_relMean,altfts,color='k')
            #ax1.fill_betweenx(altfts,prf_relMean-prf_relStd, prf_relMean+prf_relStd,alpha=0.5,color='0.75')

            ax1.plot(prf_relMean,altfts,color='k')
            ax1.errorbar(prf_relMean, altfts,yerr=altfts*0, xerr=prf_relStd,fmt='o',markersize=6,ecolor='grey', color='k', elinewidth=1.5, zorder=1)
            ax1.scatter(prf_relMean, altfts, color='k', s=45, facecolors='white', linewidths=2, zorder=2)
           
            #ax2.plot(prf_relMean2,altfts,color='k')
            #ax2.fill_betweenx(altfts,prf_relMean2-prf_relStd2, prf_relMean2+prf_relStd2,alpha=0.5,color='0.75')

            ax2.plot(prf_relMean2,altfts,color='k')
            ax2.errorbar(prf_relMean2, altfts,yerr=altfts*0, xerr=prf_relStd2,fmt='o',markersize=6,ecolor='grey', color='k', elinewidth=1.5, zorder=1)
            ax2.scatter(prf_relMean2, altfts, color='k', s=45, facecolors='white', linewidths=2, zorder=2)
            
            ax1.grid(True,which='both')
            ax2.grid(True,which='both')
            
            ax1.set_title('Relative difference (Mean +/- std)\nSAGE III interpolated') 
            ax2.set_title('Relative difference (Mean +/- std)\nSAGE III smoothed') 

            
            ax1.set_ylabel('Altitude [km]', fontsize=14)
            ax1.set_xlabel('Rel Difference (FTIR - SAGEIII / FTIR) [%]', fontsize=13)
            ax2.set_xlabel('Rel Difference (FTIR - SAGEIII / FTIR) [%]', fontsize=13)
            
            ax1.tick_params(labelsize=14)
            ax2.tick_params(labelsize=14)

            ax1.set_xlim(-150, 150)
            ax2.set_xlim(-150, 150)

            ax1.set_ylim(ymin=10, ymax=50)
            ax2.set_ylim(ymin=10, ymax=50)


            if saveFlg: 
                pdfsav.savefig(fig,dpi=200)
                        #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
            else:
                plt.show(block=False)
       


        #---------------------------------
        # Plot : Averaging Kernel Smoothing Function (row of avk)
        #---------------------------------

        avkSCFftsMean

        clmap      = 'jet'
        cm         = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(10,7))
        gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        axc       = plt.subplot(gs[2])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(altfts), vmax=np.max(altfts))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(altfts)

        #---------------------------------
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in altfts])
        
        for i in range(len(altfts)):
            ax.plot(avkSCFftsMean[i,:],altfts)
            
        ax.set_ylabel('Altitude [km]', fontsize=14)
        ax.set_xlabel('AK', fontsize=14)
        #ax.grid(True, alpha=0.5)
        #ax.set_title('(a)', loc='left', fontsize=14)
        #ax.text(0.025, 0.95,'(a)', fontsize=16,transform=ax.transAxes)   

        cbaxes = fig.add_axes([0.4, 0.55, 0.02, 0.4]) 
        cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
        cbar.set_label('Altitude [km]', fontsize=14)
        #ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
        ax.tick_params(labelsize=14)
        ax.axvline(x=0,color='k', linestyle='--')
        ax.set_ylim((0,60))
      

        #----------------------------------------
        # Calculate total column averaging kernel
        #----------------------------------------
        nobs            = len(datesfts)
        nlyrs           = len(altfts)
        avkTC           = np.zeros((nobs,nlyrs))

        for i in range(0,nobs):
            AirMtemp  = np.squeeze(Airmassfts[i,:])
            akTemp    = np.squeeze(avkSCFfts[i,:,:])
            AirMinv   = np.diag(1.0/AirMtemp)
            avkTC[i,:] = np.dot(np.dot(AirMtemp,akTemp),AirMinv)

        avkTCAv        = np.mean(avkTC, axis=0)

        #---------------------------------        
        #axb.plot(np.sum(avkSCFav[gasVer],axis=0), alt[gasVer],color='k')
        axb.plot(avkTCAv, altfts,color='k')

        #axb.grid(True,  alpha=0.5)
        axb.set_xlabel('Total Column AK', fontsize=14)
        axb.tick_params(labelsize=14)
        #axb.set_title('(b)', loc='left', fontsize=14)
        #axb.text(0.725, 0.95,'(b)', fontsize=16,transform=axb.transAxes)  

        major_ticks = np.arange(0,3, 1)
        axb.set_xticks(major_ticks) 
        axb.set_ylim((0,60))

        #---------------------------------
        dofs_cs = np.cumsum(np.diag(avkSCFftsMean)[::-1])[::-1]
        axc.plot(dofs_cs,altfts,color='k',label='Cumulative Sum of DOFS (starting at surface)')
        #axc.plot(np.sum(avkSCFav[gasVer],axis=0), alt[gasVer],color='k')
        xval = range(0,int(np.ceil(max(dofs_cs)))+2)

        if pColsFlg:

            for pcol in pCols: 
                ind1 = mf.nearestind(pcol[0], altfts)
                ind2 = mf.nearestind(pcol[1], altfts)                
                axc.fill_between(xval,altfts[ind1],altfts[ind2],alpha=0.5,color='0.75')  
                axc.axhline(altfts[ind2],color='k',linestyle='--')
                dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
                axc.text(0.15,(altfts[ind1]+altfts[ind2])/2.0, 
                         'DOFs = {2:.3f}'.format(altfts[ind1],altfts[ind2],dofsPcol),
                         fontsize=9)

        #ind1         = mf.nearestind(Pcol[0], alt[gasVer])
        #ind2         = mf.nearestind(Pcol[1], alt[gasVer])

        #axc.fill_between(xval,alt[gasVer][ind1],alt[gasVer][ind2],alpha=0.5,color='0.75')  
        #axc.axhline(alt[gasVer][ind2],color='k',linestyle='--')
        #dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
        #axc.text(0.15,(alt[idhdf][ind1]+alt[idhdf][ind2])/2.0, 
        #         'DOFs for layer {0:.2f}-{1:.2f}[km] = {2:.3f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol),
        #         fontsize=9)
        #axc.set_title('DOFs Profile - ' +str(pltID[i]))
        #axc.set_ylabel('Altitude [km]')
        
        axc.set_xlabel('Cumulative\nSum of DOFS', fontsize=14)  
        axc.tick_params(labelsize=14)
        #axc.set_title('(c)', loc='left', fontsize=14)
        #axc.text(0.025, 0.95,'(c)', fontsize=16,transform=axc.transAxes)  
        #axc.set_title('DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol), fontsize=9)    
        
        major_ticks = np.arange(0,6, 1)
        axc.set_xticks(major_ticks)

        axc.set_ylim((0,60))
        #axc.grid(True,which='both',  alpha=0.5)
      
        plt.tight_layout()

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
        
        
        #---------------------------------
        #
        #---------------------------------
        if pColsFlg:

            one2one = np.arange(0, 1e14, 1e12)

            fig, ax = plt.subplots(len(pCols), figsize=(7, 12))

            for i, pcol in enumerate(pCols): 

                inds = np.where( (altfts >= pcol[0]) & (altfts <=pcol[1])  )[0]

                TCpfts           = np.nansum(prftfts_daily3[:,inds], axis=1)#*thick[indsalt][inds]
                TCpsage_smth     = np.nansum(prfsage_smth[:,inds], axis=1)#*thick[indsalt][inds]
                TCpsage_inte     = np.nansum(prfsage_inte2[:,inds], axis=1)#*thick[indsalt][inds]

                #ax[0].errorbar(mnthVals['dates'],mnthVals['mnthlyAvg'], linestyle='None', yerr=mnthVals['std'],ecolor=clr[i-1])
                ax[i].scatter(TCpfts,TCpsage_inte, color='red', s=35, label = 'SAGEIII-ISS')
                ax[i].scatter(TCpfts,TCpsage_smth, color='green', s=35, label='SAGEIII-ISS - smoothed')
                    
                if i == 0: ax[i].legend(prop={'size':11})
                ax[i].grid(True)
                ax[i].set_ylabel('Partial Column - SAGE III \n[molec/cm$^2$]', fontsize=14)
                if i ==len(pCols)-1: ax[i].set_xlabel('Partial Column - FTIR [molec/cm$^2$]', fontsize=14)
                ax[i].set_xlim(np.min(TCpsage_smth) - np.min(TCpsage_smth)*0.2, np.max(TCpsage_smth) + np.max(TCpsage_smth)*0.2)
                ax[i].set_ylim(np.min(TCpsage_smth) - np.min(TCpsage_smth)*0.2, np.max(TCpsage_smth) + np.max(TCpsage_smth)*0.2)

                ax[i].set_title('Altitude Layer '+str(altfts[inds][-1])+' - '+str(altfts[inds][0])+' km',
                                  multialignment='center',fontsize=14)

                ax[i].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
                ax[i].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
                ax[i].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2)

                ax[i].tick_params(labelsize=14)

                print 'Altitude Layer '+str(altfts[inds][-1])+' - '+str(altfts[inds][0])+' km'
                print 'Mean +/- std - SAGEIII: {} +/- {}'.format(np.nanmean(TCpsage_smth), np.nanstd(TCpsage_smth))
                print 'Mean +/- std - FTS:     {} +/- {}'.format(np.nanmean(TCpfts), np.nanstd(TCpfts))
                

                slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(TCpfts,TCpsage_inte)
                print slopelr, interceptlr, r_valueln, p_valuelr

                ax[i].text(0.06, 0.9, 'Coeff Correlation r = {:.3f}'.format(float(r_valueln)), color='red', va='center',transform=ax[i].transAxes,fontsize=16)

                slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(TCpfts,TCpsage_smth)
                print slopelr, interceptlr, r_valueln, p_valuelr

                ax[i].text(0.06, 0.8, 'Coeff Correlation r = {:.3f}'.format(float(r_valueln)), color='green', va='center',transform=ax[i].transAxes,fontsize=16)

        fig.subplots_adjust(bottom=0.08,top=0.95, left=0.15, right=0.95)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            #plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)



        
        
    





        if saveFlg: pdfsav.close()
        else:
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()  


    #print('\nFinished Plots.......\n')

    #--------------------------------
    # Pause so user can look at plots
    #--------------------------------
    if  pltsage:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program     


if __name__ == "__main__":
    main()