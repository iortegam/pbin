import datetime as dt
from datetime import datetime, timedelta
import time
import math
import sys
import numpy as np
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
import os
import csv
import itertools
from collections import OrderedDict
import os
from os import listdir
from os.path import isfile, join
import re
#import statsmodels.api as sm
from scipy.integrate import simps
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import glob
import myfunctions as mf
from scipy import interpolate
from pylab import setp
from scipy.io import netcdf
import getopt
from itertools import izip
from numpy import *
import dataOutplts as dc
from mpl_toolkits.basemap import Basemap
import srtm
from geographiclib.geodesic import Geodesic
from dbfread import DBF

#from osgeo import ogr, osr
import matplotlib.cm as cmx
import matplotlib.colors as colors


def getseason(dates):
    
    #seas2get = {'Summer':(dt.date(2014,6,21), dt.date(2014,9,22)),
    #           'Autumn':(dt.date(2014,9,23), dt.date(2014,12,20)),
    #           'Spring':(dt.date(2014,3,21), dt.date(2014,6,20))}

    seas2get = {'Summer':(dt.date(2014,5,1), dt.date(2014,11,1))}
               #'Autumn':(dt.date(2014,9,23), dt.date(2014,12,20)),
               #'Spring':(dt.date(2014,3,21), dt.date(2014,6,20))}
    
    dd = np.asarray([dt.date(2014,d.month,d.day) for d in dates])
    
    season = []

    for indDay in dd: 

        for seas,(season_start, season_end) in seas2get.items():
            if indDay>=season_start and indDay<= season_end:
                season.append(seas)
            else:
                season.append('Winter')

    return season

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

def getcolor():
    clr  = ['red', 'blue',  'green' ]
    return clr

def createCircleAroundWithRadius(lat, lon, radiusMiles):
    #ring = ogr.Geometry(ogr.wkbLinearRing)
    latArray = []
    lonArray = []

    for brng in range(0,360, 1):
        lat2, lon2 = getLocation(lat,lon,brng,radiusMiles)
        latArray.append(lat2)
        lonArray.append(lon2)

    return lonArray,latArray


def getLocation(lat1, lon1, brng, distanceMiles):
    lat1 = lat1 * math.pi/ 180.0
    lon1 = lon1 * math.pi / 180.0
    #earth radius
    R = 6378.1 #Km
    #R = ~ 3959 MilesR = 3959

    distanceMiles = distanceMiles/R

    #brng = (brng / 90)* math.pi / 2

    brng = (float(brng) / 90)* math.pi / 2

    lat2 = math.asin(math.sin(lat1) * math.cos(distanceMiles) 
    + math.cos(lat1) * math.sin(distanceMiles) * math.cos(brng))

    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(distanceMiles)
    * math.cos(lat1),math.cos(distanceMiles)-math.sin(lat1)*math.sin(lat2))

    lon2 = 180.0 * lon2/ math.pi
    lat2 = 180.0 * lat2/ math.pi

    return lat2, lon2


class myPltClass():

    def __init__(self, gases, versions, iyear, imnth, iday, fyear, fmnth, fday, outFname='', saveFlg= False,  yrsFlg=True):

        self.ngases     = len(gases)
        self.gases      = gases
        self.versions   = versions
        self.gv         = []

        self.iyear      = iyear
        self.imnth      = imnth
        self.iday       = iday
        self.fyear      = fyear
        self.fmnth      = fmnth
        self.fday       = fday

        self.idate      = dt.date(iyear,imnth, iday)
        self.fdate      = dt.date(fyear,fmnth, fday)

        for i, gas in enumerate(gases):
            self.gv.append(gas+'_'+versions[i])

        if saveFlg:
            self.pdfsav = PdfPages(outFname)
        else: self.pdfsav = False

        if yrsFlg:
            self.yrsFlg=True
        else:
            self.yrsFlg=False

    def closeFig(self):
        self.pdfsav.close()


    #---------------------------TOTAL COLUMN PLOTS OF SEVERAL TRACE GASE----------------------
    def pltts(self, DT, TC, vmrP, TCp, pCols):
        ''' Plot Time Series of Total Column for several gases'''

        #-------------------------------------------------------
        #                          PLOTS
        #-------------------------------------------------------
      
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)    
        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        #DateFmt      = DateFormatter('%m\n%Y')
        if self.yrsFlg: DateFmt = DateFormatter('%Y')
        else: DateFmt      = DateFormatter('%m\n%Y')

        print '\nPrinting Plots:\n'

        gvStr  = self.gv
        ngases = self.ngases
        gases  = self.gases

        if ngases >= 3: 
            fsize         = (8, 13)
            bottomadjust  = 0.075
            topadjust     = 0.93
        else:
            fsize         = (8, 9)
            bottomadjust  = 0.09
            topadjust     = 0.91


        Rate_d      = []
        Rate_e_d    = []
        Rate_m      = []
        Rate_e_m    = []

        Rate_vmr_d      = []
        Rate_vmr_e_d    = []
        Rate_vmr_m      = []
        Rate_vmr_e_m    = []
        
        gasname_w   = []
        
        fig,  ax = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig2, ax2 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig3, ax3 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig4, ax4 = plt.subplots(ngases,figsize=fsize, sharex=True)
        fig5, ax5 = plt.subplots(ngases,figsize=fsize, sharex=True)


        for k, gv in enumerate(gvStr):

            #----------------------------------------------------------------------
            #analysis of time series of DAILY averages
            #----------------------------------------------------------------------
            dailyVals        = mf.dailyAvg(TC[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            DT_d             = dailyVals['dates']
            TC_d             = dailyVals['dailyAvg']
            TCsd_d           = dailyVals['std']
            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)

            #-----------------------------------
            # bootstrap resampling information of DAILY averages
            #-----------------------------------
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot

            Rate_d.append(np.divide(res[1], np.mean(dailyVals['dailyAvg']))  *100.0 )
            Rate_e_d.append(np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0 )
           
            #----------------------------------------------------------------------
            #analysis of time series of Monthly averages
            #----------------------------------------------------------------------
            mnthlyVals       = mf.mnthlyAvg(TC[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            DT_m             = mnthlyVals['dates']
            TC_m             = mnthlyVals['mnthlyAvg']
            TCsd_m           = mnthlyVals['std']
            int_fourier_m    = intercept
            sl_fourier_m     = slope
            p_fourier_m      = pfourier
            FAT_m            = f_drift(dateYearFrac)
            FATAV_m          = f_driftfourier(dateYearFrac)

            #-----------------------------------
            # bootstrap resampling information of Monthly averages
            #-----------------------------------
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

            int_boot_m       = intercept_boot
            sl_boot_m        = slope_boot
            p_boot_m         = pfourier_boot
            Rate_m.append( np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))*100.0)
            Rate_e_m.append( np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0)

            gasname = mf.getgasname(gases[k])
            gasname_w.append(gasname)

            #-----------------------------------
            ax[k].plot(DT[gv],TC[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
            ax[k].scatter(DT_d,TC_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
            #-------fitted trend
            ax[k].plot(DT_d,FAT_d,linewidth=2.5)
            ax[k].plot(DT_d,FATAV_d, linewidth=2.5)
                
            ax[k].grid(True)
            ax[k].set_ylim([np.min(TC[gv])-0.05*np.min(TC[gv]), np.max(TC[gv])+0.1*np.max(TC[gv])])
           
            if gases[k].lower() == 'ch4': ax[k].set_ylim(np.min(TC[gv])-0.02*np.min(TC[gv]), np.max(TC[gv])+0.02*np.max(TC[gv]))
     
            if self.yrsFlg:
                ax[k].xaxis.set_major_locator(yearsLc)
                ax[k].xaxis.set_minor_locator(months)
                ax[k].xaxis.set_major_formatter(DateFmt) 
                ax[k].xaxis.set_tick_params(which='major',labelsize=12)
                ax[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax[k].set_xlim(self.idate, self.fdate)
                if k == ngases-1: ax[k].set_xlabel('Year', fontsize = 16)
            else:
                ax[k].xaxis.set_major_locator(monthsAll)
                ax[k].xaxis.set_major_formatter(DateFmt)
                ax[k].set_xlim((dt.date(self.iyear,1,1), dt.date(self.iyear,12,31)))
                ax[k].xaxis.set_minor_locator(AutoMinorLocator())
                if k == ngases-1: ax[k].set_xlabel('Month/Year', fontsize = 16)

            if ngases <=2: 
                ax[k].xaxis.set_tick_params(which='major',labelsize=14)
                ax[k].yaxis.set_tick_params(which='major',labelsize=14)
                ax[k].xaxis.set_tick_params(which='minor',labelbottom='off')

              
            ax[k].set_title(gasname, fontsize=14)
            fig.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
                          horizontalalignment='right',
                          verticalalignment='center',
                          rotation='vertical',
                          transform=ax[k].transAxes)
            #fig.subplots_adjust(bottom=0.05, top=0.95, left = 0.12, right = 0.98)

            #fig.autofmt_xdate()

            fig.suptitle('All data and daily averages', fontsize=18)
            fig.subplots_adjust(bottom=bottomadjust,top=topadjust)
                   
            #-----------------------------------
            ax2[k].plot(DT[gv],TC[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
            ax2[k].scatter(DT_m,TC_m,facecolors='red', edgecolors='black', s=35, label= 'Monthly averages')
            
            #-------fitted trend 
            ax2[k].plot(DT_m,FAT_m,label='Fitted Anual Trend',linewidth=2.5)
            ax2[k].plot(DT_m,FATAV_m,label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)
            #ax2[k].text(0.02,0.85,"Fitted trend -- slope: {0:.3E} ({1:.2f}%)".format(sl_fourier_m[k], Rate_m[k]),transform=ax2[k].transAxes, fontsize = 14) 
            
            #ax2[k].errorbar(DT_m[k],TC_m[k],yerr=TCsd_m[k],fmt='k.',markersize=0,ecolor='red', capthick=2)
            ax2[k].set_title(gasname, fontsize=14)
            ax2[k].grid(True)
            ax2[k].set_ylim([np.min(TC_m)-0.1*np.min(TC_m), np.max(TC_m)+0.13*np.max(TC_m)])
            if gases[k].lower() == 'ch4': ax2[k].set_ylim(np.min(TC[gv])-0.02*np.min(TC[gv]), np.max(TC[gv])+0.02*np.max(TC[gv]))
            ##ax2[k].set_ylabel('Total Column$\cdot$[molecules cm$^{-2}$]')

           #ax2[k].text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax2[k].transAxes)
           # ax2[k].text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(ds['totClmn'])*100.0),transform=ax2[k].transAxes) 

                
            if self.yrsFlg:
                #plt.xticks(rotation=45)
                ax2[k].xaxis.set_major_locator(yearsLc)
                ax2[k].xaxis.set_minor_locator(months)
                ax2[k].xaxis.set_major_formatter(DateFmt) 
                ax2[k].xaxis.set_tick_params(which='major',labelsize=12)
                ax2[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax2[k].set_xlim((self.idate, self.fdate))
                if k == ngases-1: ax2[k].set_xlabel('Year', fontsize = 16)
            else:
                ax2[k].xaxis.set_major_locator(monthsAll)
                ax2[k].xaxis.set_major_formatter(DateFmt)
                ax2[k].set_xlim((dt.date(self.iyear,1,1), dt.date(self.iyear,12,31)))
                ax2[k].xaxis.set_minor_locator(AutoMinorLocator())
                if k == ngases-1: ax2[k].set_xlabel('Month/Year', fontsize = 16)

            if ngases <=2: 
                ax2[k].xaxis.set_tick_params(which='major',labelsize=14)
                ax2[k].yaxis.set_tick_params(which='major',labelsize=14)
                ax2[k].xaxis.set_tick_params(which='minor',labelbottom='off')

                #fig2.autofmt_xdate()

            fig2.text(0.075, 0.5, 'Total Column [molecules$\cdot$cm$^{-2}$]', fontsize = 16,
            horizontalalignment='right',
            verticalalignment='center',
            rotation='vertical',
            transform=ax2[k].transAxes)

            fig2.suptitle('All data and monthly averages $\pm$ standard deviation', fontsize=18)
            fig2.subplots_adjust(bottom=bottomadjust,top=topadjust)

            #------------------------------------
            # Plot time series of partial columns with daily averages
            #------------------------------------
            
            for pcol in pCols:

                #----------------------------------------------------------------------
                #analysis of time series of DAILY averages
                #----------------------------------------------------------------------
                dailyVals        = mf.dailyAvg(vmrP[gv], DT[gv],dateAxis=1, meanAxis=0)
                dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
                weights          = np.ones_like(dateYearFrac)
                res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
                intercept, slope, pfourier = res[0:3]
                f_drift, f_fourier, f_driftfourier = res[3:6]

                vmr_d             = dailyVals['dailyAvg']
                vmrsd_d           = dailyVals['std']
                int_fourier_d    = intercept
                sl_fourier_d     = slope
                p_fourier_d      = pfourier
                FAT_d            = f_drift(dateYearFrac)
                FATAV_d          = f_driftfourier(dateYearFrac)

                #-----------------------------------
                # bootstrap resampling information of DAILY averages
                #-----------------------------------
                perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)

                int_boot_d       = intercept_boot
                sl_boot_d        = slope_boot
                p_boot_d         = pfourier_boot
                Rate_vmr_d.append(np.divide(res[1], np.mean(dailyVals['dailyAvg']))*100.0)
                Rate_vmr_e_d.append(np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0)

                #----------------------------------------------------------------------
                #analysis of time series of Monthly averages
                #----------------------------------------------------------------------
                mnthlyVals       = mf.mnthlyAvg(vmrP[gv], DT[gv],dateAxis=1, meanAxis=0)
                dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
                weights          = np.ones_like(dateYearFrac)
                res              = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
                intercept, slope, pfourier = res[0:3]
                f_drift, f_fourier, f_driftfourier = res[3:6]

                vmr_m             = mnthlyVals['dates']
                vmr_m             = mnthlyVals['mnthlyAvg']
                vmrsd_m           = mnthlyVals['std']
                int_fourier_m    = intercept
                sl_fourier_m     = slope
                p_fourier_m      = pfourier
                FAT_m            = f_drift(dateYearFrac)
                FATAV_m          = f_driftfourier(dateYearFrac)

                #-----------------------------------
                # bootstrap resampling information of Monthly averages
                #-----------------------------------
                perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

                int_boot_m       = intercept_boot
                sl_boot_m        = slope_boot
                p_boot_m         = pfourier_boot
                Rate_vmr_m .append( np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))*100.0)
                Rate_vmr_e_m.append( np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0)

                ax3[k].plot(DT[gv],vmrP[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                ax3[k].scatter(DT_d,vmr_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')


                #-----Fitted trend
                ax3[k].plot(DT_d,FAT_d,linewidth=2.5)
                ax3[k].plot(DT_d,FATAV_d, linewidth=2.5) 
                
                ax3[k].grid(True)
                ax3[k].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                if gases[k].lower() == 'ch4': ax3[k].set_ylim(np.min(vmrP[gv])-0.02*np.min(vmrP[gv]), np.max(vmrP[gv])+0.02*np.max(vmrP[gv]))
                
                
                if self.yrsFlg:
                
                    ax3[k].xaxis.set_major_locator(yearsLc)
                    ax3[k].xaxis.set_minor_locator(months)
                    ax3[k].xaxis.set_major_formatter(DateFmt) 
                    ax3[k].xaxis.set_tick_params(which='major',labelsize=12)
                    ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                    ax3[k].set_xlim((self.idate, self.fdate))
                    if k == ngases-1: ax3[k].set_xlabel('Year', fontsize = 16)
                else:
                    ax3[k].xaxis.set_major_locator(monthsAll)
                    ax3[k].xaxis.set_major_formatter(DateFmt)
                    ax3[k].set_xlim((dt.date(self.iyear,1,1), dt.date(self.iyear,12,31)))
                    ax3[k].xaxis.set_minor_locator(AutoMinorLocator())
                    if k == ngases-1: ax3[k].set_xlabel('Month/Year', fontsize = 16)

                if ngases <=2: 
                    ax3[k].xaxis.set_tick_params(which='major',labelsize=14)
                    ax3[k].yaxis.set_tick_params(which='major',labelsize=14)
                    #ax3[k].set_ylabel(gasname + ' [ppb$_v$]', fontsize = 16)
                    ax3[k].xaxis.set_tick_params(which='minor',labelbottom='off')

                  
                ax3[k].set_title(gasname, fontsize=14)
                fig3.text(0.075, 0.5, 'VMR [ppb$_{v}$] Dry Air', fontsize = 16,
                    horizontalalignment='right',
                    verticalalignment='center',
                    rotation='vertical',
                    transform=ax3[k].transAxes)

             
                
                fig3.suptitle('All data and daily averages - weighted VMR ('+str(pCols[0][0]) + '-' + str(pCols[0][1])+ 'km)', fontsize=18)
                fig3.subplots_adjust(bottom=bottomadjust,top = topadjust)

            #------------------------------------
            # Plot time series of partial columns with monthly averages
            #------------------------------------

                ax4[k].plot(DT[gv],vmrP[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                ax4[k].scatter(DT_m,vmr_m,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                
                #----fitted trend
                ax4[k].plot(DT_m,FAT_m,label='Fitted Anual Trend',linewidth=2.5)
                ax4[k].plot(DT_m,FATAV_m,label='Fitted Anual Trend + intra-annual variability', linewidth=2.5)

                ax4[k].set_title(gasname, fontsize=14)
                ax4[k].grid(True)
                ax4[k].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                if gases[k].lower() == 'ch4': ax4[k].set_ylim(np.min(vmrP[gv])-0.02*np.min(vmrP[gv]), np.max(vmrP[gv])+0.02*np.max(vmrP[gv]))
                
                
                if self.yrsFlg:
                
                    ax4[k].xaxis.set_major_locator(yearsLc)
                    ax4[k].xaxis.set_minor_locator(months)
                    ax4[k].xaxis.set_major_formatter(DateFmt) 
                    ax4[k].xaxis.set_tick_params(which='major',labelsize=12)
                    ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')
                    ax4[k].set_xlim((self.idate, self.fdate))
                    if k == ngases-1: ax4[k].set_xlabel('Year', fontsize = 16)
                else:
                    ax4[k].xaxis.set_major_locator(monthsAll)
                    ax4[k].xaxis.set_major_formatter(DateFmt)
                    ax4[k].set_xlim((dt.date(self.iyear,1,1), dt.date(self.iyear,12,31)))
                    ax4[k].xaxis.set_minor_locator(AutoMinorLocator())
                    if k == ngases-1: ax4[k].set_xlabel('Month/Year', fontsize = 16)

                if ngases <=2: 
                    ax4[k].xaxis.set_tick_params(which='major',labelsize=14)
                    ax4[k].yaxis.set_tick_params(which='major',labelsize=14)
                    #ax4[k].set_ylabel(gasname + ' [ppb$_v$]', fontsize = 16)
                    ax4[k].xaxis.set_tick_params(which='minor',labelbottom='off')

                 
                ax4[k].set_title(gasname)
                fig4.text(0.075, 0.5, 'VMR [ppb$_{v}$] Dry Air', fontsize = 16,
                    horizontalalignment='right',
                    verticalalignment='center',
                    rotation='vertical',
                    transform=ax4[k].transAxes)
                    
                
                fig4.suptitle('All data and monthly averages - weighted VMR ('+str(pCols[0][0]) + '-' + str(pCols[0][1])+ 'km)', fontsize=18)
                fig4.subplots_adjust(bottom=bottomadjust,top=topadjust)

            #----------------------------
            # Plot total columns by month
            #----------------------------
            month    = np.array([d.month for d in DT[gv]])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))
        
            for i,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[i] = np.mean(vmrP[gv][inds])
                mnthSTD[i]  = np.std(vmrP[gv][inds])   
            
            #ax5[k].plot(mnthSort,mnthMean,'k.',markersize=6)
           
            #ax5[k].errorbar(mnthSort,mnthMean,yerr=mnthSTD, fmt='o', color='white', markersize=7, ecolor='red', capthick=2)
            ax5[k].plot(mnthSort,mnthMean, c='red', marker='o', markersize=10)   # edgecolors='red', facecolors='white',
            ax5[k].fill_between(mnthSort,mnthMean-mnthSTD,mnthMean+mnthSTD,alpha=0.5, facecolor='red', color='0.75')   

            #ax5[k].scatter(mnthSort,mnthMean, facecolors='white', edgecolors='red', alpha=0.75, s=7, label='Dayly averages')     
            ax5[k].grid(True,which='both')
            #ax5[k].set_ylabel('Retrieved Total Column\n[molecules cm$^{-2}$]',multialignment='center')
            if k == ngases-1: ax5[k].set_xlabel('Month', fontsize = 16)
           
            ax5[k].set_xlim((0,13))
            ax5[k].set_xticks(range(1,13))

            if ngases <=2: 
                ax5[k].xaxis.set_tick_params(which='major',labelsize=14)
                ax5[k].yaxis.set_tick_params(which='major',labelsize=14)
                ax5[k].set_ylabel(gasname + ' [ppb$_v$]', fontsize = 16)
                ax5[k].xaxis.set_tick_params(which='minor',labelbottom='off')

            else:  
                ax5[k].set_title(gasname)
                fig5.text(0.065, 0.5, r'VMR [ppb$_{v}$] Dry Air (monthly mean with Standard Deviation)', fontsize = 16,
                horizontalalignment='right',
                verticalalignment='center',
                rotation='vertical',
                transform=ax5[k].transAxes)

            fig5.subplots_adjust(bottom=bottomadjust,top=topadjust)



        ind = np.arange(ngases)
        fign , ax = plt.subplots()
        ax.bar(ind, Rate_d, width = 0.27, align='center', color = 'r', yerr=Rate_e_d, ecolor = 'k', label = 'Total Column')
        ax.bar(ind+0.27,Rate_vmr_d, width = 0.27, align='center', color = 'b', yerr=Rate_vmr_e_d, ecolor = 'k', label = 'weighted VMR')
        ax.yaxis.grid(True)
        ax.set_xticks(ind)
        ax.set_ylabel('Annual rate of change (%)', fontsize = 16)
        ax.set_xticklabels(gasname_w)
        ax.set_xticks(ind+0.27)
        ax.set_xlabel('Gas', fontsize = 16)
        ax.axhline(0, color='black', lw=1)
        ax.legend(prop={'size':12})
        ax.xaxis.set_tick_params(which='major',labelsize=14)
        ax.yaxis.set_tick_params(which='major',labelsize=14)
        ax.set_title('Annual rate of change (%) - Daily values')

        fign2 , ax = plt.subplots()
        ax.bar(ind, Rate_m, width = 0.27, align='center', color = 'r', yerr=Rate_e_m, ecolor = 'k', label = 'Total Column')
        ax.bar(ind+0.27,Rate_vmr_m, width = 0.27, align='center', color = 'b', yerr=Rate_vmr_e_m, ecolor = 'k', label = 'weighted VMR')
        ax.yaxis.grid(True)
        ax.set_xticks(ind)
        ax.set_ylabel('Annual rate of change (%)', fontsize = 16)
        ax.set_xticklabels(gasname_w)
        ax.set_xticks(ind+0.27)
        ax.set_xlabel('Gas', fontsize = 16)
        ax.axhline(0, color='black', lw=1)
        ax.legend(prop={'size':12})
        ax.xaxis.set_tick_params(which='major',labelsize=14)
        ax.yaxis.set_tick_params(which='major',labelsize=14)
        ax.set_title('Annual rate of change (%) - Monthly values')


        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            self.pdfsav.savefig(fig2,dpi=200)
            self.pdfsav.savefig(fig3,dpi=200)
            self.pdfsav.savefig(fig4,dpi=200)
            self.pdfsav.savefig(fig5,dpi=200)
            self.pdfsav.savefig(fign,dpi=200)
            self.pdfsav.savefig(fign2,dpi=200)
        else:
            plt.show(block=False)

    #---------------------------TOTAL COLUMN PLOTS OF SEVERAL TRACE GASE----------------------
    def pltts_2(self, DT, TC, vmrP, TCp, pCols):
        ''' Plot Time Series of Total Column for several gases'''

        #-------------------------------------------------------
        #                          PLOTS
        #-------------------------------------------------------
      
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)    
        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        #DateFmt      = DateFormatter('%m\n%Y')
        if self.yrsFlg: DateFmt = DateFormatter('%Y')
        else: DateFmt      = DateFormatter('%m\n%Y')

        print '\nPrinting Plots:\n'

        gvStr  = self.gv
        ngases = self.ngases
        gases  = self.gases


        Rate_d          = []
        Rate_e_d        = []
        Rate_m          = []
        Rate_e_m        = []

        Rate_vmr_d      = []
        Rate_vmr_e_d    = []
        Rate_vmr_m      = []
        Rate_vmr_e_m    = []
        
        gasname_w   = []
        
        fig3,   ax3 = plt.subplots(ngases/2, 2, figsize=(15, 13), sharex=True)
        fig4,   ax4 = plt.subplots(ngases/2, 2, figsize=(15, 13), sharex=True)

        for k, gv in enumerate(gvStr):

            #----------------------------------------------------------------------
            #analysis of time series of DAILY averages
            #----------------------------------------------------------------------
            dailyVals        = mf.dailyAvg(TC[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            DT_d             = dailyVals['dates']
            TC_d             = dailyVals['dailyAvg']
            TCsd_d           = dailyVals['std']
            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)


            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot

            Rate_d.append(np.divide(res[1], np.mean(dailyVals['dailyAvg']))  *100.0 )
            Rate_e_d.append(np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0 )
           
            #----------------------------------------------------------------------
            #analysis of time series of Monthly averages
            #----------------------------------------------------------------------
            mnthlyVals       = mf.mnthlyAvg(TC[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            DT_m             = mnthlyVals['dates']
            TC_m             = mnthlyVals['mnthlyAvg']
            TCsd_m           = mnthlyVals['std']
            int_fourier_m    = intercept
            sl_fourier_m     = slope
            p_fourier_m      = pfourier
            FAT_m            = f_drift(dateYearFrac)
            FATAV_m          = f_driftfourier(dateYearFrac)

            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

            int_boot_m       = intercept_boot
            sl_boot_m        = slope_boot
            p_boot_m         = pfourier_boot
            Rate_m.append( np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))*100.0)
            Rate_e_m.append( np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0)

            gasname = mf.getgasname(gases[k])
            gasname_w.append(gasname)

            #print 'hello' 
            #print 'hello'
            dailyVals        = mf.dailyAvg(vmrP[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            vmr_d             = dailyVals['dailyAvg']
            vmrsd_d           = dailyVals['std']
            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)

            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)

            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot
            Rate_vmr_d.append(np.divide(res[1], np.mean(dailyVals['dailyAvg']))*100.0)
            Rate_vmr_e_d.append(np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0)

            mnthlyVals       = mf.mnthlyAvg(vmrP[gv], DT[gv],dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            vmr_m             = mnthlyVals['dates']
            vmr_m             = mnthlyVals['mnthlyAvg']
            vmrsd_m           = mnthlyVals['std']
            int_fourier_m    = intercept
            sl_fourier_m     = slope
            p_fourier_m      = pfourier
            FAT_m            = f_drift(dateYearFrac)
            FATAV_m          = f_driftfourier(dateYearFrac)

            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)

            int_boot_m       = intercept_boot
            sl_boot_m        = slope_boot
            p_boot_m         = pfourier_boot
            Rate_vmr_m .append( np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))*100.0)
            Rate_vmr_e_m.append( np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0)

            
            if k <= 3:
                
                ax3[k, 0].plot(DT[gv],vmrP[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                ax3[k, 0].scatter(DT_d,vmr_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                ax3[k, 0].plot(DT_d,FAT_d, color='blue', linewidth=2.5)
                ax3[k, 0].plot(DT_d,FATAV_d, color='green', linewidth=2.5) 
                ax3[k, 0].grid(True)
                ax3[k, 0].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                ax3[k, 0].tick_params(labelsize=18)
                ax3[k, 0].set_title(gasname, fontsize=18)

                ax3[k, 0].xaxis.set_major_locator(yearsLc)
                ax3[k, 0].xaxis.set_minor_locator(months)
                ax3[k, 0].xaxis.set_major_formatter(DateFmt)
                ax3[k, 0].xaxis.set_tick_params(which='major',labelsize=18)
                ax3[k, 0].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax3[k, 0].set_xlim((self.idate, self.fdate))
                if k == (ngases/2)-1: ax3[k, 0].set_xlabel('Year', fontsize = 18)
                if gases[k].lower() == 'ch4': ax3[k, 0].set_ylim(np.min(vmrP[gv])-0.02*np.min(vmrP[gv]), np.max(vmrP[gv])+0.02*np.max(vmrP[gv]))

                ax3[k, 0].set_ylabel('VMR [ppb$_v$]', fontsize=18)
                

            else:
                
                ax3[k-(ngases/2), 1].plot(DT[gv],vmrP[gv], color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                ax3[k-(ngases/2), 1].scatter(DT_d,vmr_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                ax3[k-(ngases/2), 1].plot(DT_d,FAT_d, color='blue', linewidth=2.5)
                ax3[k-(ngases/2), 1].plot(DT_d,FATAV_d, color='green', linewidth=2.5) 
                ax3[k-(ngases/2), 1].grid(True)
                ax3[k-(ngases/2), 1].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                ax3[k-(ngases/2), 1].set_title(gasname, fontsize=18)

                ax3[k-(ngases/2), 1].xaxis.set_major_locator(yearsLc)
                ax3[k-(ngases/2), 1].xaxis.set_minor_locator(months)
                ax3[k-(ngases/2), 1].xaxis.set_major_formatter(DateFmt)
                ax3[k-(ngases/2), 1].xaxis.set_tick_params(which='major',labelsize=18)
                ax3[k-(ngases/2), 1].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax3[k-(ngases/2), 1].set_xlim((self.idate, self.fdate))
                if k == ngases-1: ax3[k-(ngases/2), 1].set_xlabel('Year', fontsize = 18)

                ax3[k-(ngases/2), 1].set_ylabel('VMR [ppb$_v$]', fontsize=18)
                
                ax3[k-(ngases/2), 1].tick_params(labelsize=18)

            fig3.subplots_adjust(bottom=0.065,top = 0.98, left=0.085, right=0.98)

        #if self.pdfsav:
        #   self.pdfsav.savefig(fig3, dpi=200)
        #else:
        #   plt.show(block=False)


            #----------------------------
            # Plot total columns by month
            #----------------------------
            month    = np.array([d.month for d in DT[gv]])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))

            for i,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[i] = np.mean(vmrP[gv][inds])
                mnthSTD[i]  = np.std(vmrP[gv][inds])

            year     = np.array([d.year for d in DT[gv]])
            yearSort = list(set(year))
            deltaVMR = []

            for i,y in enumerate(yearSort):
                dtMean_i  = [dt.date(y, int(m), 1) for m in mnthSort]
                dtMean_i  = np.asarray(dtMean_i)

                indsY  = np.where(year == y)[0]

                doyMean     = mf.toYearFraction(dtMean_i)
                doyYear     = mf.toYearFraction(DT[gv][indsY])

                mnthYear     = interpolate.interp1d(doyMean, mnthMean, fill_value='extrapolate', bounds_error=False, kind='linear')(doyYear)
                deltaVMR.append( (vmrP[gv][indsY] - mnthYear)/mnthYear * 100.)


            deltaVMR = np.asarray(deltaVMR)
            dVMR = [y for x in deltaVMR for y in x]
            dVMR = np.asarray(dVMR)

            #####HERE
            daysall = [dt.date(da.year, da.month, da.day) for da in DT[gv]]

            dailyVals        = mf.dailyAvg(vmrP[gv], DT[gv], dateAxis=1, meanAxis=0)
            DT_d             = dailyVals['dates']
            y_d              = dailyVals['dailyAvg']
            deltaVMR = []

            for i, item in enumerate(DT_d):
                diff = np.asarray(daysall) - item
                inds = np.where( diff == dt.timedelta(0) )[0]

                deltaVMR.append(vmrP[gv][inds] - y_d[i])
            deltaVMR = np.asarray(deltaVMR)
            dVMR = [y for x in deltaVMR for y in x]
            dVMR = np.asarray(dVMR)


            if k <= 3:
                
                ax4[k, 0].plot(DT[gv], dVMR, color='k', linestyle='None', marker ='.', markersize=4)
                #ax4[k, 0].scatter(DT_d,vmr_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                ax4[k, 0].grid(True)
                #ax4[k, 0].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                ax4[k, 0].tick_params(labelsize=18)
                ax4[k, 0].set_title(gasname, fontsize=18)

                ax4[k, 0].xaxis.set_major_locator(yearsLc)
                ax4[k, 0].xaxis.set_minor_locator(months)
                ax4[k, 0].xaxis.set_major_formatter(DateFmt)
                ax4[k, 0].xaxis.set_tick_params(which='major',labelsize=18)
                ax4[k, 0].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax4[k, 0].set_xlim((self.idate, self.fdate))
                if k == (ngases/2)-1: ax4[k, 0].set_xlabel('Year', fontsize = 18)
                #if gases[k].lower() == 'ch4': ax4[k, 0].set_ylim(np.min(vmrP[gv])-0.02*np.min(vmrP[gv]), np.max(vmrP[gv])+0.02*np.max(vmrP[gv]))

                ax4[k, 0].set_ylabel('VMR [ppb$_v$]', fontsize=18)

            else:
                
                ax4[k-(ngases/2), 1].plot(DT[gv],dVMR, color='k', linestyle='None', marker ='.', markersize=4, label='All data')
                #ax4[k-(ngases/2), 1].scatter(DT_d,vmr_d,facecolors='red', edgecolors='black', s=35, label='Dayly averages')
                ax4[k-(ngases/2), 1].grid(True)
                #ax4[k-(ngases/2), 1].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
                ax4[k-(ngases/2), 1].set_title(gasname, fontsize=18)

                ax4[k-(ngases/2), 1].xaxis.set_major_locator(yearsLc)
                ax4[k-(ngases/2), 1].xaxis.set_minor_locator(months)
                ax4[k-(ngases/2), 1].xaxis.set_major_formatter(DateFmt)
                ax3[k-(ngases/2), 1].xaxis.set_tick_params(which='major',labelsize=18)
                ax4[k-(ngases/2), 1].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax4[k-(ngases/2), 1].set_xlim((self.idate, self.fdate))
                if k == ngases-1: ax4[k-(ngases/2), 1].set_xlabel('Year', fontsize = 18)

                ax4[k-(ngases/2), 1].set_ylabel('VMR [ppb$_v$]', fontsize=18)
                
                ax4[k-(ngases/2), 1].tick_params(labelsize=18)


            fig4.subplots_adjust(bottom=0.065,top = 0.98, left=0.085, right=0.98)

        if self.pdfsav:
            self.pdfsav.savefig(fig4, dpi=200)
        else:
            plt.show(block=False)




    def pltbox(self, DT, yData, idtitle=''):
        ''' Plot box plots for different seasons'''

        gvStr  = self.gv
        ngases = self.ngases
        gases  = self.gases


        seasons   = ('Winter', 'Summer')

        #-------------------------------------------------------
        #Define parameters for plots
        #-------------------------------------------------------
        clmap        = 'jet'
        locations    = [ [1,2], [4,5], [7,8], [10,11], [13,14], [16,17], [19,20], [22,23], [25,26] ]
        lticks       = [1.5, 4.5, 7.5, 10.5, 13.5, 16.5, 19.5, 22.5, 25.5]
        clr          = ['red', 'blue', 'green']
        l = 0

        fig,  ax = plt.subplots(figsize=(10, 6))
        gasname_w = []

        for k, gv in enumerate(gvStr):
            s = getseason(DT[gv])
            s = np.asarray(s)

            Data = []
            
            for season in seasons:
                inds   = np.where( s == season )[0]
                Data.append(yData[gv][inds])
            
            if gases[k] == 'co':
                Data = np.asarray(Data)/100.0
                maxy = [np.amax(d) for d in Data] 
                ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.05), r'x100', fontsize=18)
            
            if gases[k] == 'ch4': 
                Data = np.asarray(Data)/5000.0
                maxy = [np.amax(d) for d in Data] 
                ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.1), r'x5000', fontsize=18)
                
            if gases[k] == 'o3': 
                Data = np.asarray(Data)/100.0
                maxy = [np.amax(d) for d in Data] 
                ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.05), r'x100', fontsize=18)
            


            Data = np.asarray(Data)
            meanpointprops = dict(marker='o', markeredgecolor='black',
                                  markerfacecolor='white', markersize=4)

            bp = ax.boxplot(Data, positions = [locations[k][0], locations[k][1]], widths = 0.6,  showmeans=True, meanprops=meanpointprops)
            setBoxColors(bp)
            maxloc = locations[k][1] + 1.0

            gasname = mf.getgasname(gases[k])
            gasname_w.append(gasname)
            
        ax.set_xlim((0,maxloc))
        #ax.set_yscale("log", nonposy='clip')
        ax.set_xticklabels(gasname_w)
        ax.set_xticks(lticks[0:ngases])
        ax.xaxis.set_tick_params(which='major',labelsize=18)
        ax.yaxis.set_tick_params(which='major',labelsize=18)
        ax.set_xlabel('Gas', fontsize = 18)
        ax.set_ylabel(idtitle, fontsize = 18)


        # draw temporary red and blue lines and use them to create a legend
        hB, = ax.plot([1,1],'b-')
        hR, = ax.plot([1,1],'r-')
        ax.legend((hB, hR),('Cold season', 'Warm season'), prop={'size':16})
        #legend((hB, hR),(seasons))
        hB.set_visible(False)
        hR.set_visible(False)

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)


    def pltcorr(self, DT, yData, gasx= '', idtitle='', 
               weatherFlag=False, wdir=0.0, wspeed=0.0, temp=0.0, dt_weather=dt.datetime(2010,1,1),
               qqFlag=False):

        ''' Plot correlation plot of gases versus co using wind data'''

        gvStr  = self.gv
        ngases = self.ngases
        gases  = self.gases

        #-------------------------------------------------------
        #Find version in x-axis
        #-------------------------------------------------------
        for k, g in enumerate(gases):
            if g.lower() == gasx: index_co = k

        verxplt = gases[index_co]+'_'+self.versions[index_co]

        print '\nCorrelation Analysis '+ '(' + idtitle + ')'

        #-------------------------------------------------------
        #Remove x in the list of gases
        #-------------------------------------------------------
        listgv = list(gvStr)
        listgv.remove(verxplt)

        listgas = list(gases)
        listgas.remove(gasx)
        gasname_x = mf.getgasname(gasx)

        #-------------------------------------------------------
        #Define parameters for plots
        #-------------------------------------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)    
        yearsLc      = YearLocator()
        monthsAll    = MonthLocator()
        months       = MonthLocator()
        #DateFmt      = DateFormatter('%m\n%Y')
        if self.yrsFlg: DateFmt = DateFormatter('%Y')
        else: DateFmt      = DateFormatter('%m\n%Y')
 
        fig0, ax0 = plt.subplots(ngases-1,figsize=(8, 13), sharex=True)
        fig, ax   = plt.subplots(1,ngases-1, figsize=(24,6), sharex=True)
        #fig, ax   = plt.subplots(1,ngases-1, figsize=(14,6), sharex=True)
        fig1, ax1   = plt.subplots(1,ngases-1, figsize=(24,6))
        
        if qqFlag: 
            fig2, ax2 = plt.subplots(1,ngases-1, figsize=(24,6), sharex=True)
            fig3, ax3 = plt.subplots(1,ngases-1, figsize=(24,6))
            fig4, ax4 = plt.subplots(1,ngases-1, figsize=(24,6), sharex=True)
            fig5, ax5 = plt.subplots(1,ngases-1, figsize=(24,6), sharex=True)

        #-------------------------------------------------------
        #Define the sectors for plot wind direction
        #-------------------------------------------------------
        if weatherFlag: wdir_comp = [ [0.0, 120.0],[120.0, 240.0], [240.0, 360.0]  ]

        #-------------------------------------------------------
        #                   START ANALYSIS
        #-------------------------------------------------------
        for p, g in enumerate(listgv):

            gasname = mf.getgasname(listgas[p])
            print '\nGas = ', listgas[p]

            Wdir_fix       = []
            y_wdirc        = []
            x_wdirc        = []

            #-------------------------------------------------------
            #Interpolate CO to gas
            #-------------------------------------------------------
            y         = yData[g]
            ydate     = DT[g]

            x         = yData[verxplt]
            xdate     = DT[verxplt]

            doy_w     = mf.toYearFraction(ydate)
            doy_co_w  = mf.toYearFraction(xdate)
            x_coi     = np.interp(doy_w, doy_co_w, x)
            Ratio     = np.divide(y, x_coi)

            #-------------------------------------------------------
            #WIND
            #-------------------------------------------------------
            if weatherFlag:
                doy_weather     = mf.toYearFraction(dt_weather)
                wspeed_interp   = np.interp(doy_w, doy_weather, wspeed)
                wdir_interp     = np.interp(doy_w, doy_weather, wdir)
                
                #some elements of wind are negative (delete them. they are few points)
                inds            = np.where(wdir_interp < 0.0)[0]
                x_coi           = np.delete(x_coi, inds)
                y               = np.delete(y, inds)
                ydate           = np.delete(ydate, inds)
                Ratio           = np.delete(Ratio, inds)
                
                wdir_interp2    = np.delete(wdir_interp, inds)
                #print 'elements deleted: %s' %len(inds)

                for wc in range(len(wdir_comp)): 
                    ind_wdir = np.where((wdir_interp2 >= wdir_comp[wc][0]) & (wdir_interp2 <= wdir_comp[wc][1]))[0]
                    y_wdirc.append(y[ind_wdir])
                    x_wdirc.append(x[ind_wdir])
                    wdir_interp2[ind_wdir] = (wdir_comp[wc][0] + wdir_comp[wc][1])/2.0

            #-------------------------------------------------------
            #CALCULATE ANOMALIES
            #-------------------------------------------------------
            daysall = [dt.date(da.year, da.month, da.day) for da in ydate]

            dailyVals        = mf.dailyAvg(y, ydate, dateAxis=1, meanAxis=0)
            DT_d             = dailyVals['dates']
            y_d              = dailyVals['dailyAvg']

            dailyVals        = mf.dailyAvg(x_coi, ydate, dateAxis=1, meanAxis=0)
            x_d              = dailyVals['dailyAvg']

            deltaVMRy      = []
            deltaVMRx      = []

            for i, item in enumerate(DT_d):
                diff = np.asarray(daysall) - item
                inds = np.where( diff == dt.timedelta(0) )[0]

                deltaVMRy.append(y[inds] - y_d[i])
                deltaVMRx.append(x_coi[inds] - x_d[i])


            deltaVMRy = np.asarray(deltaVMRy)
            dVMRy = [yy for xx in deltaVMRy for yy in xx]
            dVMRy = np.asarray(dVMRy)

            deltaVMRx = np.asarray(deltaVMRx)
            dVMRx = [yyy for xxx in deltaVMRx for yyy in xxx]
            dVMRx = np.asarray(dVMRx)

            #-------------------------------------------------------
            #QUANTILE - QUANTILE ANALYSIS
            #-------------------------------------------------------
            if qqFlag:
                y_q    = []
                x_q    = []

                for i, item in enumerate(DT_d):

                    diff = np.asarray(daysall) - item
                    inds = np.where( diff == dt.timedelta(0) )[0]
                         
                    y_q.append( np.divide( (y[inds] - y_d[i]), y_d[i]) *100.0 )
                    x_q.append( np.divide( (x_coi[inds] - x_d[i]), x_d[i])*100.0  )

                y_q_all    = []
                x_q_all    = []
                
                for r, tt in enumerate(y_q):
                    y_q_all = np.concatenate((tt, y_q_all))
                    x_q_all = np.concatenate((x_q[r], x_q_all))

                x_q_all      = np.asarray(x_q_all)
                y_q_all      = np.asarray(y_q_all)

                #-------------------------------------------------------------------------------
                #SORT
                #-------------------------------------------------------------------------------

                ypair = zip( *(y_q_all, y, Ratio, ydate, x_coi))
                ypair.sort()

                y_q_all_sorted  = [i1 for i1, i2, i3, i4, i5 in ypair]
                y_sorted        = [i2 for i1, i2, i3, i4, i5 in ypair]
                Ra_sorted       = [i3 for i1, i2, i3, i4, i5 in ypair]
                date_sorted     = [i4 for i1, i2, i3, i4, i5 in ypair]
                x2_sorted       = [i5 for i1, i2, i3, i4, i5 in ypair]

                xpair           = zip(x_q_all, x_coi)
                xpair.sort()

                x_q_all_sorted  = [i1 for i1, i2 in xpair]
                x_sorted        = [i2 for i1, i2 in xpair]
                

                Ra_sorted       = np.asarray(Ra_sorted)
                date_sorted     = np.asarray(date_sorted)
                x2_sorted       = np.asarray(x2_sorted)

                x_q_all_sorted   = np.asarray(x_q_all_sorted)
                y_q_all_sorted   = np.asarray(y_q_all_sorted)

                x_sorted        = np.asarray(x_sorted)
                y_sorted        = np.asarray(y_sorted)

                #-------------------------------------------------------------------------------
                #STATISTICS: QUANTILES, PERCENTILES, etc
                #------------------------------------------------------------------------------
                IQR = stats.iqr(y_q_all_sorted, rng=(25,75))
                print 'IQR = ', str(IQR)

                PCy = np.percentile(y_q_all_sorted, [75, 90])
                print 'PCy = ', str(PCy)
                PCx = np.percentile(x_q_all_sorted, [75, 90])
                print 'PCx = ', str(PCx)

                Qy = stats.mstats.hdquantiles(y_q_all_sorted, prob=(0.1,0.9))
                print 'Qy = ', str(Qy)
                Qx = stats.mstats.hdquantiles(x_q_all_sorted, prob=(0.1,0.9))
                print 'Qx = ', str(Qx)

                indsfit = np.where( (y_q_all_sorted > Qy[0]) & (y_q_all_sorted < Qy[1]) & (x_q_all_sorted > Qx[0]) & (x_q_all_sorted < Qx[1]) )[0]
                inds = np.where( (y_q_all_sorted > Qy[1]) )[0]#| (y_q_all_sorted < Qy[0]) ) [0]

                print 'Ratio in the fit region = ', str(np.mean(Ra_sorted[indsfit]))
                print 'Ratio out the fit region = ', str(np.mean(Ra_sorted[inds]))
                #-------------------------------------------------------------------------------
                #UNSORT X AND Y BASED AND SORT BASED ON DATES
                #-------------------------------------------------------------------------------

                ypair = zip( *(date_sorted, y_sorted, x2_sorted, Ra_sorted))
                ypair.sort()

                ydate2        = [i1 for i1, i2, i3, i4 in ypair]    
                y2            = [i2 for i1, i2, i3, i4 in ypair]
                x2            = [i3 for i1, i2, i3, i4 in ypair]
                Ratio2        = [i4 for i1, i2, i3, i4 in ypair]

                ypair = zip( *(date_sorted[inds], y_sorted[inds], x2_sorted[inds], Ra_sorted[inds]))
                ypair.sort()

                ydate3        = [i1 for i1, i2, i3, i4 in ypair]    
                y3            = [i2 for i1, i2, i3, i4 in ypair]
                x3            = [i3 for i1, i2, i3, i4d in ypair]
                Ratio3        = [i4 for i1, i2, i3, i4 in ypair]

            #-------------------------------------------------------
            #TIME SERIES OF THE RATIOS
            #-------------------------------------------------------
            dailyVals        = mf.dailyAvg(Ratio, ydate, dateAxis=1, meanAxis=0)
            dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            weights          = np.ones_like(dateYearFrac)
            res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            ydate_d          = dailyVals['dates']
            Ratio_d          = dailyVals['dailyAvg']
            ysd_d            = dailyVals['std']
            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)



            #-----------------------------------
            # bootstrap resampling information of DAILY averages
            #-----------------------------------
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)

            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot
            #Rate_vmr_d.append(np.divide(res[1], np.mean(dailyVals['dailyAvg']))*100.0)
            #Rate_vmr_e_d.append(np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0)

            #print 'Data sets must have same size in corrcoef: xData = {}, yData = {}'.format(xData.shape[0],yData.shape[0])

            print 'Mean (+/-) of Ratio between {0:} and  {1:}:  {2:.3f} +/-  {3:.4f}'.format(str(listgas[p]), str(gasname_x), np.mean(Ratio), np.std(Ratio))
            mu     = np.mean(Ratio)
            sigma  = np.std(Ratio)


            #-------------------------------------------------------
            #             PLOT: TIME SERIES OF THE RATIOS
            #-------------------------------------------------------
            ax0[p].plot(ydate, Ratio, color='k', linestyle='None', marker ='.', markersize=4, label='All data')
            ax0[p].scatter(ydate_d,Ratio_d,facecolors='red', edgecolors='black', s=35, label='Dayly Averages')

            if listgas[p].upper() == 'C2H6': ax0[p].fill_between(ydate_d, mu+sigma, mu-sigma, color = 'yellow', facecolor='yellow', alpha=0.5)
        
            if qqFlag: ax0[p].scatter(ydate3,Ratio3,facecolors='green', edgecolors='black', s=35)
            #-------fitted trend
            ax0[p].plot(DT_d,FAT_d,linewidth=2.5)
            ax0[p].plot(DT_d,FATAV_d, linewidth=2.5)
             
            ax0[p].grid(True)
            ax0[p].set_ylim([np.min(Ratio)-0.05*np.min(Ratio), np.max(Ratio)+0.1*np.max(Ratio)])
            
            if self.yrsFlg:
                ax0[p].xaxis.set_major_locator(yearsLc)
                ax0[p].xaxis.set_minor_locator(months)
                ax0[p].xaxis.set_major_formatter(DateFmt) 
                ax0[p].xaxis.set_tick_params(which='major',labelsize=12)
                ax0[p].xaxis.set_tick_params(which='minor',labelbottom='off')
                ax0[p].set_xlim(self.idate, self.fdate)
                if p == ngases-2: ax0[p].set_xlabel('Year', fontsize = 16)
            else:
                ax0[p].xaxis.set_major_locator(monthsAll)
                ax0[p].xaxis.set_major_formatter(DateFmt)
                ax0[p].set_xlim((dt.date(self.iyear,1,1), dt.date(self.iyear,12,31)))
                ax0[p].xaxis.set_minor_locator(AutoMinorLocator())
                if p == ngases-2: ax0[p].set_xlabel('Month/Year', fontsize = 16)

            if ngases <=2: 
                ax0[p].xaxis.set_tick_params(which='major',labelsize=14)
                ax0[p].yaxis.set_tick_params(which='major',labelsize=14)
                ax0[p].xaxis.set_tick_params(which='minor',labelbottom='off')

              
            ax0[p].set_title(gasname, fontsize=14)
            fig0.text(0.075, 0.5, 'Ratio of gas with '+gasname_x +' - '+idtitle , fontsize = 16,
                          horizontalalignment='right',
                          verticalalignment='center',
                          rotation='vertical',
                          transform=ax0[p].transAxes)

            #fig.autofmt_xdate()
            #fig.suptitle('All data and daily averages', fontsize=18)
            fig0.subplots_adjust(bottom=0.07,top=0.95, left=0.15)

            
            #-------------------------------------------------------
            #             PLOT: CORRELATION PLOT (it can be color coded by wind)
            #-------------------------------------------------------
            ###if weatherFlag: norm = colors.Normalize(vmin=np.nanmin(wdir_interp2), vmax=np.nanmax(wdir_interp2) )
            bounds = [0.0]
            if weatherFlag:
                bounds.extend([int(wdir_comp[ii][1]) for ii in range(len(wdir_comp))])   #   [wdir_comp[0][0], 90, 180, 270, 360]
                norm = colors.BoundaryNorm(bounds, cm.N)
                sc = ax[p].scatter(x_coi, y, s =25, c=wdir_interp2, cmap=cm,norm=norm) 
    
                #cax  = fig.add_axes([0.25, 0.92, 0.5, 0.05])
                #cbar = fig.colorbar(sc, cax=cax, format='%2i', orientation='horizontal', boundaries= bounds )
                #cbar.set_label('Wind direction [relative to North]')
                
            else:
                ax[p].scatter(x_coi, y, s =25, c='red', alpha=0.75, edgecolors='k')

            slope, intercept, r_value, p_value, std_err = stats.linregress(x_coi,y)
            print 'Slope = ', str(slope)
            print 'Intercept = ', str(intercept)
            print 'r-value = ', str(r_value)

            if listgas[p].upper() == 'C2H6':

                xx = range(50, 250, 1)
                xx = np.asarray(xx)
                
                FitNG = 0.18*np.asarray(xx) + (-13.4)#*math.exp(-1.6/7.4)
                FitNGp = FitNG + FitNG*0.2   #(0.18+0.18*0.15)*np.asarray(xx) + (-13.4)
                FitNGm = FitNG - FitNG*0.2   #(0.18-0.18*0.15)*np.asarray(xx) + (-13.4)

                FitAP = 0.021*np.asarray(xx) + (-0.71)#*math.exp(-1.6/7.4)
                FitAPp = FitAP + FitAP*0.2 #(0.021+0.021*0.15)*np.asarray(xx) + (-0.71)
                FitAPm = FitAP - FitAP*0.2 #(0.021-0.021*0.15)*np.asarray(xx) + (-0.71)
                

                FitBG = 0.013*np.asarray(xx) + (-0.029)#*math.exp(-1.6/7.4)
                FitBGp = FitBG + FitBG*0.2 #(0.013+0.013*0.15)*np.asarray(xx) + (-0.029)
                FitBGm = FitBG - FitBG*0.2 #(0.013-0.013*0.15)*np.asarray(xx) + (-0.029)

                ####FitNG = 0.18*math.exp(-1.6/7.4)*xx + (-13.4)*math.exp(-1.6/7.4)
                ####FitAP = 0.021*math.exp(-1.6/7.4)*xx + (-0.71)*math.exp(-1.6/7.4)
                ####FitBG = 0.013*math.exp(-1.6/7.4)*xx + (-0.029)*math.exp(-1.6/7.4)

                #ax[p].plot(xx, FitNG, color='red',  linewidth=2.0, linestyle='--', label='O&NG')
                #####ax[p].plot(xx, FitNGp, color='red',  linewidth=2.0, linestyle='--', label='O&NG')
                #####ax[p].plot(xx, FitNGm, color='red',  linewidth=2.0, linestyle='--', label='O&NG')

                #ax[p].fill_between(xx,FitNGm,FitNGp,alpha=0.25, color='r')
                
                #ax[p].plot(np.asarray(xx), FitAP, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #####ax[p].plot(np.asarray(xx), FitAPp, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #####ax[p].plot(np.asarray(xx), FitAPm, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #ax[p].fill_between(xx,FitAPm,FitAPp,alpha=0.25, color='blue')
                
                #ax[p].plot(np.asarray(xx), FitBG, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #####ax[p].plot(np.asarray(xx), FitBGp, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #####ax[p].plot(np.asarray(xx), FitBGm, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #ax[p].fill_between(xx,FitBGm,FitBGp,alpha=0.25, color='green')
                #ax[p].set_ylim(0, 3.0)
                #ax[p].set_xlim(50, 220)
    

            if listgas[p].upper() == 'NH3':

                xx = range(50, 250, 1)
                xx = np.asarray(xx)
                
                FitNG = 0.31*np.asarray(xx) + (-25.8)#*math.exp(-1.6/7.4)
                FitNGp = FitNG + FitNG*0.2   #(0.18+0.18*0.15)*np.asarray(xx) + (-13.4)
                FitNGm = FitNG - FitNG*0.2   #(0.18-0.18*0.15)*np.asarray(xx) + (-13.4)

                FitAP = 0.014*np.asarray(xx) + (-0.49)#*math.exp(-1.6/7.4)
                FitAPp = FitAP + FitAP*0.2 #(0.021+0.021*0.15)*np.asarray(xx) + (-0.71)
                FitAPm = FitAP - FitAP*0.2 #(0.021-0.021*0.15)*np.asarray(xx) + (-0.71)
                

                FitBG = 0.0073*np.asarray(xx) + (-0.25)#*math.exp(-1.6/7.4)
                FitBGp = FitBG + FitBG*0.2 #(0.013+0.013*0.15)*np.asarray(xx) + (-0.029)
                FitBGm = FitBG - FitBG*0.2 #(0.013-0.013*0.15)*np.asarray(xx) + (-0.029)

                #####FitNG = 0.18*math.exp(-1.6/7.4)*xx + (-13.4)*math.exp(-1.6/7.4)
                #####FitAP = 0.021*math.exp(-1.6/7.4)*xx + (-0.71)*math.exp(-1.6/7.4)
                #####FitBG = 0.013*math.exp(-1.6/7.4)*xx + (-0.029)*math.exp(-1.6/7.4)

                #ax[p].plot(xx, FitNG, color='red',  linewidth=2.0, linestyle='--', label='O&NG')
                #####ax[p].plot(xx, FitNGp, color='red',  linewidth=2.0, linestyle='--', label='O&NG')
                #####ax[p].plot(xx, FitNGm, color='red',  linewidth=2.0, linestyle='--', label='O&NG')

                #ax[p].fill_between(xx,FitNGm,FitNGp,alpha=0.25, color='r')
                
                #ax[p].plot(np.asarray(xx), FitAP, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #####ax[p].plot(np.asarray(xx), FitAPp, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #####ax[p].plot(np.asarray(xx), FitAPm, color='blue',  linewidth=2.0, linestyle='--', label='Anthro')
                #ax[p].fill_between(xx,FitAPm,FitAPp,alpha=0.25, color='blue')
                
                #ax[p].plot(np.asarray(xx), FitBG, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #####ax[p].plot(np.asarray(xx), FitBGp, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #####ax[p].plot(np.asarray(xx), FitBGm, color='green',  linewidth=2.0, linestyle='--', label='BKG')
                #ax[p].fill_between(xx,FitBGm,FitBGp,alpha=0.25, color='green')

                #ax[p].set_ylim(0, 1.5)

            ax[p].tick_params(axis='both', which='major', labelsize=12)
            ax[p].grid(True)

           
            
            if p == 0: ax[p].set_ylabel(idtitle,multialignment='center',fontsize = 14)
            #if p == (ngases/2.0): ax[p].set_xlabel(r'CO Total Column [molecules$\cdot$cm$^{-2}$]',multialignment='center',fontsize = 14)
            fig.text(0.5, 0.05, gasname_x + ' '+idtitle, fontsize = 14,
            horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax[p].transAxes)
            #if p == len(listgv)/2.0: ax[p].set_xlabel(r'CO '+idtitle,multialignment='center',fontsize = 14)
            ax[p].set_title(gasname)
             
            #ax[p].plot(x_coi, fit, 'k', label='Fitted line')

            fig.subplots_adjust(bottom=0.15,top=0.9, left=0.075, right=0.97)
            

            if weatherFlag:
                print bounds
                cax  = fig.add_axes([0.25, 0.92, 0.5, 0.05])
                cbar = fig.colorbar(sc, cax=cax, format='%2i', orientation='horizontal', boundaries= bounds )
                cbar.set_label('Wind direction [relative to North]', fontsize=14)


            #-------------------------------------------------------
            #   PLOT: CORRELATION PLOT OF ANOMALIES
            #-------------------------------------------------------

            ax1[p].scatter(dVMRx, dVMRy, s =25, c='red', alpha=0.75, edgecolors='k')
            ax1[p].grid(True)
            ax1[p].tick_params(axis='both', which='major', labelsize=14)
            ax1[p].set_title(gasname)
            ax1[p].set_ylabel('$\Delta$'+gasname + ' [ppb]',multialignment='center',fontsize = 14)
            ax1[p].set_xlabel('$\Delta$'+gasname_x + ' [ppb]',multialignment='center',fontsize = 14)

            #print '\nMean (+/-) of Ratio between {0:} and  {1:}:  {2:.3f} +/-  {3:.4f} (Anomalies)'.format(str(listgas[p]), str(gasname_x), np.mean(dVMRy/dVMRx), np.std(dVMRy/dVMRx))
            #mu     = np.mean(Ratio)
            #sigma  = np.std(Ratio)
            
            slope, intercept, r_value, p_value, std_err = stats.linregress(dVMRx,dVMRy)
            print 'Slope (Anomaly) = ', str(slope)
            print 'Intercept  (Anomaly) = ', str(intercept)
            print 'r-value  (Anomaly) = ', str(r_value)

            fig1.subplots_adjust(bottom=0.1,top=0.9, left=0.075, right=0.97)

            #-------------------------------------------------------
            #   PLOT: QQ
            #-------------------------------------------------------
                
            if qqFlag:

                slope, intercept, r_value, p_value, std_err = stats.linregress(x_q_all_sorted[indsfit], y_q_all_sorted[indsfit])
                fit = slope*x_q_all_sorted + intercept

                indy1 = mf.nearestind(Qy[0], y_q_all_sorted)
                indy2 = mf.nearestind(Qy[1], y_q_all_sorted)    

                indx1 = mf.nearestind(Qx[0], x_q_all_sorted)
                indx2 = mf.nearestind(Qx[1], x_q_all_sorted)

                print 'Fraction of O&NG:', str( float(len(inds))/float(len(y_q_all_sorted)) )
                #ydate2 = np.asarray(ydate2)
                #print 'Dates in green:', ydate3

                #-------------------------------------------------------
                #             PLOT: Q-Q 
                #-------------------------------------------------------        
                ax2[p].axhline(Qy[0],color='gray',linestyle='--', linewidth=2.0)
                ax2[p].axhline(Qy[1],color='gray',linestyle='--', linewidth=2.0)

                ax2[p].axvline(Qx[0],color='gray',linestyle='--', linewidth=2.0)
                ax2[p].axvline(Qx[1],color='gray',linestyle='--', linewidth=2.0)

                ax2[p].scatter(x_q_all_sorted,y_q_all_sorted, s =25, c='red', alpha=0.75, edgecolors='k')
                ax2[p].scatter(x_q_all_sorted[inds],y_q_all_sorted[inds], s =25, c='green', alpha=0.75, edgecolors='k')
                
                ax2[p].plot(x_q_all_sorted, fit, 'k',  linewidth=2.0, linestyle='--', label='Fitted line')
                ax2[p].grid(True)

                #ax2[p].broken_barh([(110, 30)], (10, 9), facecolors='blue')

                if p == 0: ax2[p].set_ylabel('$\Delta$ [$\%$]',multialignment='center',fontsize = 14)
                ax2[p].set_title(gasname)

                fig2.suptitle(idtitle, fontsize=14)
                fig2.text(0.5, 0.075, '$\Delta_{'+ gasx.upper()+'}$ [$\%$]', fontsize = 14,
                horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax2[p].transAxes)
                fig2.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)

                #-------------------------------------------------------
                #             PLOT: FREQUENCY
                #-------------------------------------------------------

                res = stats.relfreq(y_q_all_sorted, numbins=25)
                x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)

                ax3[p].bar(x, res.frequency, width=res.binsize)
                ax3[p].set_ylim(bottom=0)
                ax3[p].set_xlim([x.min(), x.max()]) 
                ax3[p].tick_params(axis='both',which='both')   #labelsize=12
                ax3[p].axvline(Qy[0],color='gray',linestyle='--', linewidth=2.0)
                ax3[p].axvline(Qy[1],color='gray',linestyle='--', linewidth=2.0)
                if p == 0: ax3[p].set_ylabel('Frequency',multialignment='center',fontsize = 14)
                ax3[p].set_title(gasname)
                #ax3[p].set_ylim(gasname)
                #ax3[p].grid(True)

                fig3.text(0.5, 0.075, '$\Delta$ [$\%$]', fontsize = 14,
                horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax3[p].transAxes)
                fig3.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)

                #-------------------------------------------------------
                #             PLOT: CORRELATION PLOT (highlighting "local plumes") 
                #-------------------------------------------------------
                ax4[p].scatter(x2, y2, s =25, c='red', alpha=0.75, edgecolors='k')
                ax4[p].scatter(x3, y3, s =25, c='green', alpha=0.75, edgecolors='k')
                ax4[p].grid(True)
                #ax4[p].set_ylim(bottom=0)
                if p == 0: ax[p].set_ylabel(idtitle,multialignment='center',fontsize = 14)
                fig4.text(0.5, 0.05, gasname_x+' '+idtitle, fontsize = 14,
                horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax[p].transAxes)
                ax4[p].set_title(gasname)   
                fig4.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)

                #-------------------------------------------------------
                #             PLOT: CORRELATION PLOT (highlighting "local plumes") 
                #-------------------------------------------------------
                ax5[p].scatter(x_q_all, y_q_all, s =25, c='red', alpha=0.75, edgecolors='k')
                ax5[p].scatter(x_q_all[inds],y_q_all_sorted[inds], s =25, c='green', alpha=0.75, edgecolors='k')
                ax5[p].grid(True)
                if p == 0: ax[p].set_ylabel(idtitle,multialignment='center',fontsize = 14)
                fig5.text(0.5, 0.05, gasname_x+' '+idtitle, fontsize = 14,
                horizontalalignment='center', verticalalignment='center', rotation='horizontal', transform=ax[p].transAxes)
                ax5[p].set_title(gasname)   
                fig5.subplots_adjust(bottom=0.15,top=0.78, left=0.05, right=0.97)


        if self.pdfsav:
            self.pdfsav.savefig(fig,dpi=200)
            self.pdfsav.savefig(fig0,dpi=200)
            self.pdfsav.savefig(fig2,dpi=200)
            self.pdfsav.savefig(fig3,dpi=200)
            self.pdfsav.savefig(fig4,dpi=200)
        else: plt.show(block=False)
        

class MapClass():

    def __init__(self, origin=(0,0), maxAlt=5000.0, minAlt=500., DistLeft=25000,  DistRight=25000,  DistTop=25000,  DistBottom=25000, saveFlg=False):

        #------------------------------------------------------------------
        #             THE FOLLOWING INFORMATION IS TO CREEATE THE MAP
        #------------------------------------------------------------------
        self.origin = origin   #CENTER OF THE MAP
        max_altitude = maxAlt               #MAXIMUM ALTITUDE TO CONSIDER IN THE TERRAIN
        min_altitude = minAlt

        azimuth, distance = (270, DistLeft)  #DISTANCE ON THE LEFT (in meters) 
        self.geolft = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

        azimuth, distance = (90, DistRight)   #RIGHT
        self.georgt = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

        azimuth, distance = (0,DistTop)     #TOP
        self.geotop = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

        azimuth, distance = (180,DistBottom)   #BOTTOM
        self.geobot = Geodesic.WGS84.Direct(origin[0],origin[1],azimuth,distance)

        domain_size = (500, 500)          #DOMAIN SIZE (roughly the RESOLUTION OF THE MAP)
        #domain_size = (100, 100)

        self.domain=[domain_size,
                (self.geobot['lat2'],self.geotop['lat2']),
                (self.geolft['lon2'],self.georgt['lon2']),
                max_altitude, min_altitude]

        self.geodata = srtm.get_data()

        self.elev_array = self.geodata.get_image(*self.domain, mode='array')

        self.lats=np.linspace(self.domain[1][0],self.domain[1][1],self.domain[0][1])
        self.lons=np.linspace(self.domain[2][0],self.domain[2][1],self.domain[0][0])

        self.SaveFlg=saveFlg 



    def pltMapZ(self, LonData=0., LatData=0., zData=0., zmin=0., zmax=0.0, ztitle='', LatID=0., LonID=0., SaveFile=''):

        '''
        Initially Developed to plot FRAPPE (C-130 Data)
        This plot a map with Lat and Lon color coded by zData:
        LonData and Lat Data = Vectors with Longitude
        zData                = Vector with zData
        zmin and zmax        = Optional min and max values for the colo code plot
        ztitle               = String for title of the color coded plot
        latID and LonID      = Optional coordinates to plot a triangle
        SaveFile             = Name and path of the map to save

        '''

        if self.SaveFlg: pdfsav = PdfPages(SaveFile)
        else:        pdfsav = False

        fig,ax = plt.subplots(figsize=(10.5,7))

        m= Basemap(llcrnrlon=self.geolft['lon2'], llcrnrlat=self.geobot['lat2'], urcrnrlon=self.georgt['lon2'], urcrnrlat=self.geotop['lat2'],
                     resolution='i', projection='lcc', lat_0 = self.origin[0], lon_0 = self.origin[1], ax=ax) # projection='tmerc'

        cmap = plt.get_cmap('terrain')

        # x, y = m(self.lons, self.lats)
        # p = m.pcolormesh(x,y,self.elev_array, cmap=cmap, vmin=self.domain[4], vmax=self.domain[3])
        # p.cmap.set_over(cmap(.0))
        # p.cmap.set_under('w')  
        
        ax.tick_params(labelsize=14)

        clmap   = 'jet'
        cm      = plt.get_cmap(clmap)

        xx, yy = m(LonData,LatData)

        if (zmin and zmax): 
            norm = colors.Normalize( vmin=zmin, vmax=zmax)
            sc = m.scatter(xx, yy, c=zData, cmap=cm, s=4,  edgecolors='none', norm=norm)
        else:
            sc = m.scatter(xx, yy, c=zData, cmap=cm, s=4,  edgecolors='none')

        cax  = fig.add_axes([0.88, 0.14, 0.03, 0.8])
        cbar = fig.colorbar(sc, cax=cax)
        
        cbar.set_label(ztitle, fontsize=14)

        ax.tick_params(axis='both', which='major', labelsize=14)


        m.drawcounties()
        # draw parallels.
        parallels = np.arange(30.,50, 1.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14, linewidth=0.0, color='gray')
        # draw meridians
        meridians = np.arange(180.,360., 1.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14, linewidth=0.0, color='gray')
        m.drawmapscale(-107.1, 42.05 , self.origin[1], self.origin[0], 100, barstyle='fancy')

        x, y = m(LonID, LatID)
        x2, y2 = m(LonID-0.15, LatID-0.25)
        ax.plot(x,y,'k^', markersize=15)
        plt.text(x2,y2,'NCAR', color='k', fontsize=16)


        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.86)

        # -------
        radii = 40

        locP = [ [-104.6, 40.42] , [-105.245, 39.6], [-106.6, 39.7]]

        for l in locP:

            centerlon = l[0]
            centerlat = l[1]

            X,Y = createCircleAroundWithRadius(centerlat,centerlon, radii)

            X,Y = m(X,Y)
            ax.plot(X,Y,marker=None,color='green',linewidth=2)

            x,y = m(centerlon,centerlat)
            ax.plot(x,y ,marker='D',color='green',markersize=4)
        
        if self.SaveFlg: 
           pdfsav.savefig(fig,dpi=200)
           
        else:           
            plt.show(block=False)

        if self.SaveFlg: pdfsav.close()


    def pltMapRF(self, LonData=0., LatData=0., zData=0., ztitle='', LatID=0., LonID=0., locP=(0.,0.), SaveFile=''):

        '''
        Initially Developed to plot All RF FRAPPE (C-130 Data) - Note the Lat and Lon Data is differennt than pltMapZ

        This plot a map with Lat and Lon color coded by zData:
        LonData and Lat Data = Vectors with Longitude
        zData                = Vector with zData
        ztitle               = String for title of the color coded plot
        latID and LonID      = Optional coordinates to plot a triangle
        locP                 = e.g.,  [ [], [], [] ]   with locations of Long and Lat to plot
        SaveFile             = Name and path of the map to save

        '''

        if self.SaveFlg: pdfsav = PdfPages(SaveFile)
        else:        pdfsav = False

        fig,ax = plt.subplots(figsize=(10.5,7))

        m= Basemap(llcrnrlon=self.geolft['lon2'], llcrnrlat=self.geobot['lat2'], urcrnrlon=self.georgt['lon2'], urcrnrlat=self.geotop['lat2'],
                     resolution='i', projection='lcc', lat_0 = self.origin[0], lon_0 = self.origin[1], ax=ax) # projection='tmerc'

        cmap = plt.get_cmap('terrain')

        x, y = m(self.lons, self.lats)
        p = m.pcolormesh(x,y,self.elev_array, cmap=cmap, vmin=self.domain[4], vmax=self.domain[3])
        p.cmap.set_over(cmap(.0))
        p.cmap.set_under('w')  
        ax.tick_params(labelsize=14)

        c = cmap_discretize('jet', len(LonData))
        cNorm          = colors.Normalize( vmin=np.nanmin(zData), vmax=np.nanmax(zData) )
        scalarMap      = mplcm.ScalarMappable(norm=cNorm,  cmap=c )
        scalarMap.set_array(c)

        ax.set_color_cycle( [scalarMap.to_rgba(x) for x in zData] )

        for i in range(len(LonData)):
            xx, yy = m(LonData[i],LatData[i])
            m.plot(xx, yy,  '.', markersize = 2,linestyle='None', linewidth=2, alpha=0.5)
        
        ax.tick_params(axis='both', which='major', labelsize=14)

        cax  = fig.add_axes([0.88, 0.14, 0.03, 0.8])
        cbar = fig.colorbar(scalarMap, cax=cax)
        cbar.set_label('Research Flight', fontsize=14)
        
        lab = np.arange(1, len(LonData)+1, 1)
        loc = lab# + 0.5

        cbar.ax.set_yticklabels(lab)  # vertically oriented colorbar
        cbar.set_ticks(loc)
        cbar.set_ticklabels(lab)

        m.drawcounties()
        # draw parallels.
        parallels = np.arange(30.,50, 1.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=14, linewidth=0.0, color='gray')
        # draw meridians
        meridians = np.arange(180.,360., 1.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=14, linewidth=0.0, color='gray')
        m.drawmapscale(-107.1, 42.05 , self.origin[1], self.origin[0], 100, barstyle='fancy')

        x, y = m(LonID, LatID)
        x2, y2 = m(LonID-0.15, LatID-0.25)
        ax.plot(x,y,'k^', markersize=15)
        #ax.text(x2,y2,'NCAR', color='k', fontsize=16)

        clr = getcolor()

        for j, p in enumerate(locP):
            x3, y3 = m(p[0], p[1])
            ax.plot(x3,y3,'r^', color=clr[j], markersize=15)
            #ax.scatter(x3,y3, facecolors='white', s=60, color=clr[j])
        
            Delta_deg = 0.5
            m.tissot(x3, y3, Delta_deg, 100, facecolor='green', edgecolor='black',linewidth=1, alpha=0.5)
            #circle = plt.Circle((x, y), radius=150.0, fc='y')
            #plt.gca().add_patch(circle)

            # define the epicentral distance Delta and plot this as a circle on the map:
            
        fig.subplots_adjust(left = 0.12, bottom=0.12, top=0.96, right = 0.86)
        
        if self.SaveFlg: 
           pdfsav.savefig(fig,dpi=100)
           
        else:           
            plt.show(block=False)

        if self.SaveFlg: pdfsav.close()

