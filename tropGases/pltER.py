#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltts.py
#
# Purpose:
#         Plot time series of multiple species retrieved with FTIR columns/VMRs
#         Note: See below for inputs
#
# Notes:
#   
#
# Version History:
#       Created, May, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
import numpy.ma as ma
import sys
import glob

from scipy import interpolate

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
from itertools import izip
from numpy import *
import myfunctions as mf
import dataOutClass as dc
from collections                     import OrderedDict
import PltClass as mp

from scipy import linspace, polyval, polyfit, sqrt, stats, randn


import ClassFTS as cfts
import myfunctions as mf

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

def closeFig(self):
    self.pdfsav.close()

def sumzip(*items):
    return [sum(values) for values in zip(*items)]

                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #-----------------------------------------------------------------------------------------
    #                             Initializations for FTIR
    #-----------------------------------------------------------------------------------------
    loc                = 'fl0'

    gasName            = ['co',  'c2h6' , 'ch4', 'nh3', 'h2co', 'hcooh', 'hcn']
    ver                = ['Current_v3', 'Current_v2', 'Current_WP', 'Current_v2', 'Current_v8', 'Current_v3', 'Current_v1']            

    saveFlg            = False 

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2017
    fmnth              = 12
    fday               = 31

    #----------------------
    # Date range for FRAPPE
    #----------------------
    frappe_i           = dt.datetime(2014, 7, 14)
    frappe_f           = dt.datetime(2014, 8, 20)
    
    #----------------------
    # First Level Retrieval Directory
    #----------------------
    pltDir            =  '/data/iortega/pbin/tropGases/fig/'
    DataDir           =  '/data/iortega/results/'+loc.lower()+'/data/'

    pltFile           =  pltDir + 'pltER.pdf'

    #----------------------
    # FLEXPART
    #----------------------
    FlexDir = '/ya4/Campaign/'+loc.upper()+'/FRAPPE//flexpart2/FLX_SCF_12h.txt'

    cols, indexToName = mf.getColumns(FlexDir, headerrow=1, delim=',', header=True)

    SCsd = np.asarray(cols['Date'])  
    SCdate = [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]), int(d[8:10])) for d in SCsd]
    SCdate = np.asarray(SCdate)
    SCNE   = np.asarray(cols['SC_NE'])
    SCNW   = np.asarray(cols['SC_NW'])
    SCSE   = np.asarray(cols['SC_SE'])
    SCSW   = np.asarray(cols['SC_SW'])

    #----------------------
    # Wind - EOL
    #----------------------
    weatherDir      = '/data1/ancillary_data/fl0/eol/'
    weatherFileTag  = 'v2'

    wdir, wspeed, temp, rh, dt_weather = mf.weatherout(loc, weatherDir, weatherFileTag, iyear, fyear )

     

                                    #----------------------------#
                                    #                            #
                                    #        --- START ---       #
                                    #                            #
                                    #----------------------------#

    gasFile   = [gasName[i]+'_'+ver[i] for i in range(len(gasName))]  

    if saveFlg: pdfsav = PdfPages(pltFile)

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------
    fts = {}

    for gi, gf in enumerate(gasFile):
        cols, indexToName = mf.getColumns(DataDir + gf + '.dat', headerrow=0, delim=',', header=True)

        for index in indexToName:
            fts[gasName[gi]+'_'+indexToName[index]]  = cols[indexToName[index]]

    #---------------------------------
    #Add Time Stamp
    #---------------------------------
    for g in gasName:
        
        date = fts[g+'_YYYY-MM-DD']
        time    = fts[g+'_HH:MM:SS']

        yyyy = [int(d[0:4]) for d in date]
        mm   = [int(d[5:7]) for d in date]
        dd   = [int(d[8:10]) for d in date]

        hh   = [int(d[0:2]) for d in time]
        mi   = [int(d[3:5]) for d in time]
        ss   = [int(d[6:8]) for d in time]

        dateTime =   [dt.datetime(yyyy[i], mm[i], dd[i], hh[i], mi[i], ss[i]) for i, d in enumerate(dd)]

        fts[g+'_dt'] = dateTime
    
    #---------------------------------
    #  -- Define Global Parameters --
    #---------------------------------
    ngases       = len(gasFile)
    npanels      = int(math.ceil(ngases/2.0))

    if iyear == fyear: yrsFlg = False
    else: yrsFlg = True

    yearsLc            = YearLocator()
    months             = MonthLocator()
    dayLc              = DayLocator()
    mondays            = WeekdayLocator(MONDAY)

    if yrsFlg: DateFmt = DateFormatter('%Y')
    else: DateFmt      = DateFormatter('%m\n%Y')

    #---------------------------------
    #     Analysis
    #---------------------------------
    vmr           = {}
    dates         = {}

    residuals     = {}
    residuals_bkg  = {}

    dvmr          = {}
    
    vmr_Enhan     = {}
    res_Enhan     = {}
    dates_Enhan   = {}
    dates_Enhan2  = {}
    dvmr_Enhan    = {}

    resMean       = {}
    resStd        = {}

    gas_corr      = {}

    vmr_NE         = {}
    dates_NE       = {}

    vmr_SE         = {}
    dates_SE       = {}

    vmr_WE         = {}
    dates_WE       = {}
    dates_WE2       = {}

    #-----------
    #
    #----------
    vmr_EnhanCO    = {}
    dvmr_EnhanCO   = {}
    vmrCO_WE       = {}


    for i, g in enumerate(gasName):

        dates[g]           = np.asarray(fts[g+'_dt'])
        vmr[g]             = np.asarray(fts[g+'_wVMR'], dtype=np.float32)

        # yearsall  = [ singDate.year for singDate in dates[g]]
        # yearsall  = np.asarray(yearsall)
        # yearsList = list(set(yearsall))
        # yearsList.sort()

        # for y in yearsList:
        #     indsy = np.where(yearsall == y)
        #     print indsy

        # exit()
        #----------------------------
        # Drift & Fourier
        #----------------------------
        dateYearFrac        = mf.toYearFraction(dates[g])
        weights             = np.ones_like(dateYearFrac)
        res                 = mf.fit_driftfourier(dateYearFrac, vmr[g], weights, 2, half_period=1)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        #----------------------------
        # Get Residuals
        #----------------------------
        res                 = vmr[g] - f_driftfourier(dateYearFrac)
        residuals[g]        = res

        res_Mean            = np.mean(res)
        res_Std             = np.std(res)

        resMean[g]          = res_Mean
        resStd[g]           = res_Std
        #----------------------------
        # Get Residuals above (possible enhancements due to O&NG)
        #----------------------------
        #vmr_above           = []
        #residual_above      = []
        #dt_above            = []

        inds                = np.where(res >= (res_Mean + res_Std*2))[0]

        vmr_Enhan[g]        = np.asarray(vmr[g][inds])
        dates_Enhan[g]      = np.asarray(dates[g][inds])
        dates_Enhan2[g]     = np.asarray(dates[g][inds])
        res_Enhan[g]        = np.asarray(res[inds])

        #------------------------------
        # Interpolate CO to enhancements
        #------------------------------
        datesCO           = np.asarray(fts['co_dt'])
        dateYearFracCO    = mf.toYearFraction(datesCO)
        vmrCO             = np.asarray(fts['co_wVMR'], dtype=np.float32)
        vmr_EnhanCO[g]    = interpolate.interp1d(dateYearFracCO, vmrCO, kind='nearest', bounds_error=False)(dateYearFrac[inds])

        #for ii, jj in enumerate(res):

        #     if jj > (res_Mean + res_Std*1.5):
        #         residual_above.append(jj)
        #         dt_above.append(dates[g][ii])
        #         vmr_above.append(vmr[g][ii])

        # vmr_Enhan[g]   = np.asarray(vmr_above)
        # dates_Enhan[g]  = np.asarray(dt_above)
        # res_Enhan[g]   = np.asarray(residual_above)

        #----------------------------
        # Get NE dates and possible Enhancements (possible enhancements due to O&NG)
        #----------------------------
        ydoywind         = mf.toYearFraction(dt_weather)
        wdir_interp      = interpolate.interp1d(ydoywind, wdir, axis=0, bounds_error=False)(dateYearFrac)

        inds_wd          = np.logical_and(wdir_interp >= 0., wdir_interp <= 90.)
        vmr_NE[g]        = vmr[g][inds_wd]
        dates_NE[g]      = dates[g][inds_wd]

        inds_wd          = np.logical_and(wdir_interp >= 90., wdir_interp <= 180.)
        vmr_SE[g]        = vmr[g][inds_wd]
        dates_SE[g]      = dates[g][inds_wd]

        inds_wd          = np.logical_and(wdir_interp >= 180., wdir_interp <= 360.)
        vmr_WE[g]        = vmr[g][inds_wd]
        dates_WE[g]      = dates[g][inds_wd]
        dates_WE2[g]     = dates[g][inds_wd]

        residuals_bkg[g] = residuals[g][inds_wd]

        doy_WE           = mf.toYearFraction(dates_WE[g])
        doy_Enhan        = mf.toYearFraction(dates_Enhan[g])

        #----------------------------
        #
        #----------------------------
        if g.upper() != 'HCN':
            intrsctVals      = np.intersect1d(doy_WE, doy_Enhan, assume_unique=False)
                        
            inds1            = np.nonzero( np.in1d( doy_WE, intrsctVals, assume_unique=False ) )[0]
            inds2            = np.nonzero( np.in1d( doy_Enhan, intrsctVals, assume_unique=False ) )[0]

            dates_Enhan[g]   = np.delete(dates_Enhan[g] ,inds2,axis=0)
            vmr_Enhan[g]     = np.delete(vmr_Enhan[g] ,inds2,axis=0)

            dates_WE[g]       = np.delete(dates_WE[g] ,inds1,axis=0)
            vmr_WE[g]         = np.delete(vmr_WE[g] ,inds1,axis=0)

        #----------------------------
        #
        #----------------------------
            vmr_EnhanCO[g]      = np.delete(vmr_EnhanCO[g] ,inds2,axis=0)

        vmrCO_WE[g]         = interpolate.interp1d(dateYearFracCO, vmrCO, kind='nearest', bounds_error=False)(dateYearFrac[inds_wd])
        if g.upper() != 'HCN': vmrCO_WE[g]         = np.delete(vmrCO_WE[g] ,inds1,axis=0)

        #----------------------------
        # Substracting background
        #----------------------------
        month               = np.array([d.month for d in dates[g]])
        mnthSort            = list(set(month))
        mnthMean            = np.zeros(len(mnthSort))
        mnthSTD             = np.zeros(len(mnthSort))

        month_bkg           = np.array([d.month for d in dates_WE[g]])
       
        mnthSort_bkg        = list(set(month_bkg))
        mnthMean_bkg        = np.zeros(len(mnthSort))
        mnthSTD_bkg         = np.zeros(len(mnthSort))
        #dmonth_bkg          = np.zeros(len(mnthSort))

        month_Enhan         = np.array([d.month for d in dates_Enhan[g]])

        dvmr_Enhan[g]       = np.zeros(len(vmr_Enhan[g]))
        dvmr[g]             = np.zeros(len(dates[g]))

        #----------------------------
        # 
        #----------------------------
        dvmr_EnhanCO[g]       = np.zeros(len(vmr_Enhan[g]))
        mnthMean_bkgCO        = np.zeros(len(mnthSort))

  
        for i,m in enumerate(mnthSort):

            inds           = np.where(month == m)[0]
            mnthMean[i]    = np.mean(vmr[g][inds])
            mnthSTD[i]     = np.std(vmr[g][inds])

            #dmonth_bkg[i]   = dt.date(2015, int(m), 15)

            inds_bkg        = np.where(month_bkg == m)[0]
            mnthMean_bkg[i] = np.median(vmr_WE[g][inds_bkg])
            mnthSTD_bkg[i]  = np.std(vmr_WE[g][inds_bkg])   

            dvmr[g][inds]       = vmr[g][inds] - mnthMean_bkg[i]

            inds                = np.where(month_Enhan == m)[0]
            dvmr_Enhan[g][inds] = vmr_Enhan[g][inds] - mnthMean_bkg[i]

            #----------------------------
            #
            #---------------------------- 
            mnthMean_bkgCO[i]     = np.min(vmrCO_WE[g][inds_bkg]) #- np.std(vmrCO_WE[g][inds_bkg]) 
            dvmr_EnhanCO[g][inds] = vmr_EnhanCO[g][inds] - mnthMean_bkgCO[i]

        print 'Gas: {}'.format(g)
        print 'Number of days with Enahncements  ({}): {}'.format(g, float(len(dvmr_Enhan[g])))
        print 'Percent of days with Enahncements ({}): {}'.format(g, float(len(dvmr_Enhan[g]))/float(len(dates[g])) * 100.)

        #print len(dvmr_EnhanCO[g])
        indspos = np.where(dvmr_EnhanCO[g] > 0)[0]
        indsneg = np.where(dvmr_EnhanCO[g] < 0)[0]
        #print len(indspos)
        #print len(indsneg)


        #odr, odrErr  = mf.orthoregress(vmr_EnhanCO[g], vmr_Enhan[g], xerr=vmr_EnhanCO[g]*0.05, yerr=vmr_Enhan[g]*0.075, InError=True)
        #slope = float(odr[0])
        #intercept = float(odr[1])
        #print slope 
        
        #odr, odrErr  = mf.orthoregress(dvmr_EnhanCO[g], dvmr_Enhan[g], xerr=dvmr_EnhanCO[g]*0.05, yerr=dvmr_Enhan[g]*0.075, InError=True)
        #slope = float(odr[0])
        #intercept = float(odr[1])
        #print slope 


        #fig,ax  = plt.subplots(figsize=(7,9))
        #ax.scatter(dvmr_EnhanCO[g], dvmr_Enhan[g], color='gray', s=25, facecolors='white', zorder=1)
        #ax.scatter(dates_Enhan[g], dvmr_Enhan[g], color='red', s=25, facecolors='white', zorder=1)  
        #plt.show(block=False)
    
    #user_input = raw_input('Press any key to exit >>> ')

    #---------------------------------
    #     PLOT: -- Time Series --
    #---------------------------------
    ngases       = len(gasFile)
    npanels      = int(math.ceil(ngases/2.0))

    xmin      = dt.date(iyear, imnth, iday)
    xmax      = dt.date(fyear, fmnth, fday)
    clmap     = 'jet'

    #fig0 = plt.figure(figsize=(8,12))
    #fig  = plt.figure(figsize=(8,12))
    #fig2 = plt.figure(figsize=(8,12))
    #fig3 = plt.figure(figsize=(8,12))

    fig0 = plt.figure(figsize=(12,10))
    fig  = plt.figure(figsize=(12,10))
    fig2 = plt.figure(figsize=(12,10))
    fig3 = plt.figure(figsize=(12,10))


    #outer_grid = gridspec.GridSpec(ngases, 1, wspace=0.15, hspace=0.125)
    outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName):

        gasStr = mf.getgasname(g)

        #----------------------------
        # PLOT: Time Series of VMR (all Data and drift & f_driftfourier)
        #----------------------------
        #----------------------------
        # Drift & Fourier
        #----------------------------
        dateYearFrac        = mf.toYearFraction(dates[g])
        weights             = np.ones_like(dateYearFrac)
        res                 = mf.fit_driftfourier(dateYearFrac, vmr[g], weights, 2)
        intercept, slope, pfourier = res[0:3]
        f_drift, f_fourier, f_driftfourier = res[3:6]

        #---------------------------------------------------
        #To make a continuous fit
        #---------------------------------------------------
        numdays = (dates[g].max() + dt.timedelta(days=1) - dates[g].min()).days
        dates2  = [dates[g].min() + dt.timedelta(days=x) for x in range(0, numdays)]
        dates2  = np.asarray(dates2)
        dateYearFrac2 = mf.toYearFraction(dates2)

        ax0 = plt.Subplot(fig0, outer_grid[i])

        ax0.plot(dates[g],vmr[g], color='gray', linestyle='None', marker ='.', markersize=4, zorder=1)
        ax0.scatter(dates[g],vmr[g], color='r', s=25, facecolors='white', zorder=2)
        ax0.plot(dates2,f_driftfourier(dateYearFrac2), color='green', linewidth=3.5, zorder=3)
        ax0.plot(dates2,f_drift(dateYearFrac2), color='blue', linewidth=3.5, zorder=3)
        
        ax0.grid(True, alpha=0.35)
        ax0.tick_params(labelsize=16)
        ax0.grid(True,which='both', alpha=0.35)  
        ax0.text(0.03, 0.9, gasStr, va='center',transform=ax0.transAxes,fontsize=24)

        ax0.set_xlim(xmin, xmax)

        ax0.xaxis.set_major_locator(yearsLc)
        ax0.xaxis.set_minor_locator(months)
        ax0.xaxis.set_major_formatter(DateFmt)

        if i == 0: 
            ax0.set_title('Weigthed VMR [1.6 - 8.0 km]', fontsize=18)
            #ax0.legend(prop={'size':12}, loc = 1)

        #ax.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
        #ax.set_xlim(frappe_i, frappe_f)

        fig0.add_subplot(ax0)

        #----------------------------
        # PLOT: Residuals in VMR (red points enhancements, dotted line STD*2)
        #----------------------------

        ax = plt.Subplot(fig, outer_grid[i])
 
        ax.scatter(dates[g],residuals[g], color='gray', s=25, facecolors='white', zorder=1, label='All')
        ax.scatter(dates_Enhan2[g],res_Enhan[g], color='red', s=25, facecolors='white', zorder=2, label='Enhancements')
        ax.scatter(dates_WE2[g],residuals_bkg[g], color='blue', s=25, facecolors='white', zorder=3,  label='Background')

        ax.axhline(resMean[g]+resStd[g]*2,color='k',linestyle='--', linewidth=2)
        ax.axhline(resMean[g]-resStd[g]*2,color='k',linestyle='--', linewidth=2)

        ax.grid(True, alpha=0.35)
        ax.tick_params(labelsize=16)
        ax.grid(True,which='both', alpha=0.35)  
        ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

        #ax2_2 = ax2.twinx()
        #ax2_2.plot(dt_weather, wdir, color='green', linewidth=2)

        ax.set_xlim(xmin, xmax)

        ax.xaxis.set_major_locator(yearsLc)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(DateFmt)


        #if i == 0: 
           # ax.set_title('Detrended Weigthed VMR & Enhancements\nDotted line represents natural variability', fontsize=18)
        #    ax.legend(prop={'size':11}, loc='upper center', bbox_to_anchor=(0.5, 1.2), fancybox=True, ncol=3) 


        #ax2.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
        #ax2.set_xlim(frappe_i, frappe_f)

        fig.add_subplot(ax)

        #----------------------------
        # PLOT: Time Series of VMR (all Data in gray, red points enhancements and blue point background vmr)
        #----------------------------
        ax2 = plt.Subplot(fig2, outer_grid[i])
  
        ax2.scatter(dates[g],vmr[g], color='gray', s=25, facecolors='white', zorder=1, label='All')
        ax2.scatter(dates_Enhan[g],vmr_Enhan[g], color='red', s=25, facecolors='white', zorder=2, label='Enhancements')
        ax2.scatter(dates_WE[g],vmr_WE[g], color='blue', s=25, facecolors='white', zorder=3,  label='Background')

        if g == 'co': 
            ax2.scatter(dates_Enhan[g],vmr_EnhanCO[g], color='orange', s=25, facecolors='white', zorder=2, label='Enhancements')
            ax2.scatter(dates_WE[g],vmrCO_WE[g], color='purple', s=25, facecolors='white', zorder=4,  label='Background')
        
        ax2.grid(True, alpha=0.35)
        ax2.tick_params(labelsize=16)
        ax2.grid(True,which='both', alpha=0.35)  
        ax2.text(0.03, 0.9, gasStr, va='center',transform=ax2.transAxes,fontsize=24)

        ax2.set_xlim(xmin, xmax)
        ax2.set_ylim(ymin = 0)

        ax2.xaxis.set_major_locator(yearsLc)
        ax2.xaxis.set_minor_locator(months)
        ax2.xaxis.set_major_formatter(DateFmt)
        
        if i == 0: 
            #ax2.set_title('Identifying background VMR', fontsize=18)        
            ax2.legend(prop={'size':12}, loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, ncol=3) 

        fig2.add_subplot(ax2)

        #----------------------------
        # PLOT:delta VMR
        #----------------------------

        ax3 = plt.Subplot(fig3, outer_grid[i])

        #ax3.scatter(dates_Enhan[g],vmr_Enhan[g], color='red', s=25, facecolors='white', zorder=1)
        ax3.scatter(dates_Enhan[g],dvmr_Enhan[g], color='red', s=35, facecolors='white', zorder=2)
        if g == 'co': ax3.scatter(dates_Enhan[g],dvmr_EnhanCO[g], color='orange', s=35, facecolors='white', zorder=2)

        ax3.grid(True, alpha=0.35)
        ax3.tick_params(labelsize=16)
        ax3.grid(True,which='both', alpha=0.35)  
        ax3.text(0.03, 0.9, gasStr, va='center',transform=ax3.transAxes,fontsize=24)

        #ax2_2 = ax2.twinx()
        #ax2_2.plot(dt_weather, wdir, color='green', linewidth=2)

        ax3.set_xlim(xmin, xmax)

        ax3.xaxis.set_major_locator(yearsLc)
        ax3.xaxis.set_minor_locator(months)
        ax3.xaxis.set_major_formatter(DateFmt)

        if i == 0: 
            ax3.set_title('Substracting background in enhancements', fontsize=18)        
            #ax3.legend(prop={'size':12}, loc = 1)

        #ax3.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
        #ax2.set_xlim(frappe_i, frappe_f)

        fig3.add_subplot(ax3)

    all_axes0 = fig0.get_axes()
    all_axes  = fig.get_axes()
    all_axes2 = fig2.get_axes()
    all_axes3 = fig3.get_axes()
    #show only the outside spines
    for ax in all_axes0:
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

    for ax in all_axes2:
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

    for ax in all_axes3:
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


    all_axes0[-2].set_xlabel('Year', fontsize=16)
    all_axes[-2].set_xlabel('Year', fontsize=16)
    all_axes2[-2].set_xlabel('Year', fontsize=16)
    all_axes3[-2].set_xlabel('Year', fontsize=16)

    all_axes0[-1].set_xlabel('Year', fontsize=16)
    all_axes[-1].set_xlabel('Year', fontsize=16)
    all_axes2[-1].set_xlabel('Year', fontsize=16)
    all_axes3[-1].set_xlabel('Year', fontsize=16)
    
    fig0.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')
    fig.text(0.02, 0.5, 'Residuals  [ppb]', fontsize=18, va='center', rotation='vertical')
    fig2.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')
    fig3.text(0.02, 0.5, 'Delta -VMR [ppb]', fontsize=18, va='center', rotation='vertical')

    fig0.autofmt_xdate()
    fig.autofmt_xdate()
    fig2.autofmt_xdate()
    fig3.autofmt_xdate()

    fig0.subplots_adjust(bottom=0.095,top = 0.94, left=0.11, right=0.95)
    fig.subplots_adjust(bottom=0.096,top = 0.94, left=0.11, right=0.95)
    fig2.subplots_adjust(bottom=0.095,top = 0.94, left=0.11, right=0.95)
    fig3.subplots_adjust(bottom=0.095,top = 0.94, left=0.11, right=0.95)

    if saveFlg:
        pdfsav.savefig(fig0,dpi=200)
        pdfsav.savefig(fig,dpi=200)
        pdfsav.savefig(fig2,dpi=200)
        pdfsav.savefig(fig3,dpi=200)
        #fig.savefig(pltDir+'Time_Series.pdf', bbox_inches='tight')
        fig.savefig(pltDir+'Time_Series_residuals.pdf', bbox_inches='tight')
        fig2.savefig(pltDir+'Time_Series_enh.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)
    
    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()
 

    #--------------------------------
    # Plot: Show by month
    #--------------------------------

    fig = plt.figure(figsize=(8,12))

    outer_grid = gridspec.GridSpec(ngases, 1, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName):

        gasStr = mf.getgasname(g)

        #----------------------------
        # 
        #----------------------------
        dates_Mnth        = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates[g]]
        dates_Enhan_Mnth  = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates_Enhan[g]]
        dates_WE_Mnth     = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates_WE[g]]
    

        ax = plt.Subplot(fig, outer_grid[i])

        #ax.scatter(dates_Mnth,vmr[g], color='gray', s=10, facecolors='white', zorder=1, label='All')
        ax.plot(dates_Mnth,vmr[g], color='gray', linestyle='None', marker ='.', markersize=10, label='All', zorder=1)
        ax.scatter(dates_Enhan_Mnth,vmr_Enhan[g], color='red', s=35, facecolors='white', zorder=2, label='Enhancements')
        ax.scatter(dates_WE_Mnth,vmr_WE[g], color='blue', s=35, facecolors='white', zorder=3,  label='Background')
        #ax.plot(dmonth_bkg[g],mnthMean_bkg[g], color='blue', linestyle='-', zorder=4, linewidths=2)

    
        
        ax.grid(True, alpha=0.35)
        ax.tick_params(labelsize=16)
        ax.grid(True,which='both', alpha=0.35)  
        ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

        ax.set_xlim(dt.date(2015,1,1), dt.date(2015,12,31))
        ax.set_ylim(ymin=0)

        ax.xaxis.set_major_locator(months)
        #ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(DateFormatter('%b'))

        if i == 0: 
            ax.legend(prop={'size':12}, loc='upper center', bbox_to_anchor=(0.5, 1.2), fancybox=True, ncol=3) 
            #ax0.set_title('Weigthed VMR [1.6 - 8.0 km]', fontsize=18)
            #ax0.legend(prop={'size':12}, loc = 1)

        #ax.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
        #ax.set_xlim(frappe_i, frappe_f)

        fig.add_subplot(ax)
    
    fig.autofmt_xdate()
    #fig.subplots_adjust(bottom=0.085,top = 0.94, left=0.11, right=0.95)

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

    all_axes[-1].set_xlabel('Month', fontsize=16)
    
    fig.text(0.02, 0.5, 'wVMR [ppb]', fontsize=18, va='center', rotation='vertical')

    #fig.autofmt_xdate()

    fig.subplots_adjust(bottom=0.085,top = 0.96, left=0.11, right=0.95)

    if saveFlg:
        pdfsav.savefig(fig,dpi=200)
    
    else:
        plt.show(block=False)

    plt.savefig(pltDir+'Enhancement_mnth.pdf', bbox_inches='tight')    


    #---------------------------------
    #     PLOT: -- Time Series --
    #---------------------------------
    dates_x      = np.asarray(dates['co'])
    vmr_x        = np.asarray(vmr['co'])
    doy_x        = mf.toYearFraction(dates_x)

    datesEnh_x   = np.asarray(dates_Enhan['co'])
    doyEnh_x     = mf.toYearFraction(datesEnh_x)   
    vmrEnh_x     = np.asarray(vmr_Enhan['co'])
    dvmrEnh_x    = np.asarray(dvmr_Enhan['co'])

    dates_x_WE   = np.asarray(dates_WE['co'])

    #gasName2     = [g for g in gasName if g.lower() != 'co']

    gasName2     = ['c2h6', 'nh3', 'h2co', 'hcooh', 'hcn']
    MX           = [ 30., 17., 30., 46., 27.]
    MCO          = 28. 
 
    #yearsall  = [ singDate.year for singDate in dates_Enhan['c2h6']]
    #yearsall  = np.asarray(yearsall)
    #yearsList = list(set(yearsall))
    #yearsList.sort()

    #yearsList     = [ [2010, 2014], [2015, 2017] ]
    #yearsListStr  = ['2010 - 2014', '2015 - 2017']

    yearsList     = [ [2010, 2017] ]
    yearsListStr  = ['2010 - 2017']

    dER1          = {}
    dERstd1       = {}
    moi1          = [5, 6, 7, 8, 9, 10]

    dER2         = {}
    dERstd2      = {}
    moi2         = [1, 2, 3, 4, 11, 12]

    #---------------------------------
    #
    #---------------------------------
    fig  = plt.figure(figsize=(12,8))
    fig2  = plt.figure(figsize=(12,8))
    
    #outer_grid   = gridspec.GridSpec(len(gasName2), 1, wspace=0.15, hspace=0.125)
    #outer_grid = gridspec.GridSpec(int(math.ceil(len(gasName2)/2.0))  , 2, wspace=0.15, hspace=0.125)

    outer_grid = gridspec.GridSpec(2  , 2, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName2):

        gasStr = mf.getgasname(g)

        # da = []
        # for dd in dates[g]:
        #     da.append(dt.date( int(dd.year), int(dd.month), int(dd.day) ) )

        # self.datefts = np.asarray(self.datefts)
        # self.sonde['date'] = np.asarray(self.sonde['date'])


        # doy_fts   = dc.toYearFraction(self.datefts)
        # doy_sonde = dc.toYearFraction(self.sonde['date'])

        # #-------------------------------------------------------
        # #FINDING COINCIDENT DATES
        # #-------------------------------------------------------
        # intrsctVals = np.intersect1d(doy_fts, doy_sonde, assume_unique=False)
        
        # inds1       = np.nonzero( np.in1d( doy_sonde, intrsctVals, assume_unique=False ) )[0]
        # inds2       = np.nonzero( np.in1d( doy_fts, intrsctVals, assume_unique=False ) )[0]

        datesEnh         = np.asarray(dates_Enhan[g])
        vmrEnh           = np.asarray(vmr_Enhan[g])
        dvmrEnh          = np.asarray(dvmr_Enhan[g])
        vmr_y             = np.asarray(vmr[g])
        doyEnh           = mf.toYearFraction(datesEnh) 

        doy              = mf.toYearFraction(dates[g])   

        vmr_x_interpol      = interpolate.interp1d(doy_x, vmr_x, axis=0,fill_value='extrapolate', bounds_error=False)(doy)
        Ratio               = vmr_y/vmr_x_interpol
        
        vmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, vmrEnh_x, axis=0, fill_value='extrapolate', bounds_error=False)(doyEnh)
        RatioEnh            = vmrEnh/vmrEnh_x_interpol

        dvmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, dvmrEnh_x, axis=0, kind='nearest', fill_value='extrapolate',bounds_error=False)(doyEnh)

        dRatioEnh            = dvmrEnh/dvmrEnh_x_interpol

        #yearG1  = [ singDate.year for singDate in dates_Enhan[g] if singDate.month in moi1]
        #yearG1  = np.asarray(yearG1)

        yearG1  = [ singDate.year for singDate in dates_Enhan[g]]
        yearG1  = np.asarray(yearG1)
        
        #dERstd1[g] = np.asarray(dERstd1[g], dtype=np.float32)

        dRatioEnh2   = dvmrEnh/np.asarray(dvmr_EnhanCO[g])

        dates_Enhan2 = dates_Enhan[g]

        #yearG1  = [ singDate.year for singDate in dates_Enhan[g]]
        #yearG1  = np.asarray(yearG1)
        

        if gasStr.upper() == 'HCOOH':
            indsok = np.where(dRatioEnh2 < 0.4)[0]

            dRatioEnh2  = dRatioEnh2[indsok]
            dates_Enhan2 = dates_Enhan[g][indsok]

        if gasStr.upper() == 'HCN':
            indsok = np.where(dRatioEnh2 < 0.2)[0]

            dRatioEnh2  = dRatioEnh2[indsok]
            dates_Enhan2 = dates_Enhan[g][indsok]

        for y in yearsList:

            yearG1  = [ singDate.year for singDate in dates_Enhan2]
            yearG1  = np.asarray(yearG1)
            #indsy = np.where(yearG1 == y)[0]

            indsy = np.where( (yearG1 >= y[0]) &  (yearG1 <= y[1]))[0]

            #print y
            #print yearG1[indsy]
     
            if len(indsy) >=1 :
                #dER1.setdefault(g,[]).append(np.nanmean(dRatioEnh[indsy]))
                #dERstd1.setdefault(g,[]).append(np.nanstd(dRatioEnh[indsy]))

                dER1.setdefault(g,[]).append(dRatioEnh2[indsy])
                #dERstd1.setdefault(g,[]).append(np.nanstd(dRatioEnh[indsy]))

            else:
                dER1.setdefault(g,[]).append(float('nan'))
                dERstd1.setdefault(g,[]).append(float('nan'))

        dER1[g]    = np.asarray(dER1[g], dtype=np.float32)


        indspos = np.where((dRatioEnh2 > 0) & (dRatioEnh2 < 0.4) )[0]

        print 'Negative percent dCO = {}'.format( (len(dRatioEnh2) - len(indspos)) /float(len(dRatioEnh2)) )

        odr, odrErr  = mf.orthoregress(vmrEnh_x_interpol, vmrEnh, xerr=vmrEnh_x_interpol*0.05, yerr=vmrEnh*0.075, InError=True)
        slope = float(odr[0])
        intercept = float(odr[1])
        print slope 

        slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(vmrEnh_x_interpol, vmrEnh)
        rvalue_i = float(r_valueln)
        print rvalue_i
        
        odr, odrErr  = mf.orthoregress(dvmrEnh_x_interpol, dvmrEnh, xerr=dvmrEnh_x_interpol*0.05, yerr=dvmrEnh*0.075, InError=True)
        slope = float(odr[0])
        intercept = float(odr[1])
        print slope 
       
        #----------------------------
        # 
        #----------------------------
        if gasStr.upper() == 'HCN': 
            EFCO = 89.0
           
        else: EFCO = 59.91

        
        print 'Mean Enhancement ratio of {} : {} +/- {}'.format(gasStr, np.nanmean(dRatioEnh), np.nanstd(dRatioEnh) )
        print 'Mean Enhancement - 2 ratio of {} : {} +/- {}'.format(gasStr, np.nanmean(dRatioEnh2[indspos]), 2.*np.nanstd(dRatioEnh2[indspos])/np.sqrt(len(dRatioEnh2[indspos])) )

        print 'Median Enhancement ratio of {} : {} +/- {}'.format(gasStr, np.nanmedian(dRatioEnh), np.nanstd(dRatioEnh) )
        print 'Median Enhancement - 2 ratio of {} : {} +/- {}'.format(gasStr, np.nanmedian(dRatioEnh2[indspos]), 2.*np.nanstd(dRatioEnh2[indspos])/np.sqrt(len(dRatioEnh2[indspos])) )

        print 'Median Emmision factor ratio of {} : {} +/- {}'.format(gasStr, EFCO *np.nanmedian(dRatioEnh)* (MX[i]/MCO), EFCO * np.nanstd(dRatioEnh)* (MX[i]/MCO) )
        print 'Median Emmision factor ratio of {} : {} +/- {}'.format(gasStr, EFCO *np.nanmedian(dRatioEnh2[indspos])* (MX[i]/MCO), EFCO * 2.*np.nanstd(dRatioEnh2[indspos])/np.sqrt(len(dRatioEnh2[indspos]))*(MX[i]/MCO) )
        
        if gasStr.upper() == 'HCN':
            print '\nFIRE HCN EMISSION'
            d1 = dt.datetime(2015, 8, 20)
            d2 = dt.datetime(2015, 8, 22)

            indsDates = np.where( (datesEnh[indspos] > d1) & (datesEnh[indspos] < d2) )[0]
            
            Factor = np.exp(1.5/75.)/np.exp(1.5/30)

            print Factor
            print datesEnh[indspos][indsDates]

            #print 'Median Enhancement ratio of {} : {} +/- {}'.format(gasStr, np.nanmedian(dRatioEnh[indsDates]), np.nanstd(dRatioEnh[indsDates]) )
            print 'Median Enhancement - 2 ratio of {} : {} +/- {}'.format(gasStr, np.nanmean(dRatioEnh2[indspos][indsDates]), 2.*np.nanstd(dRatioEnh2[indspos][indsDates])/np.sqrt(len(dRatioEnh2[indspos][indsDates])) )
            print 'Median Enhancement w Factor - 2 ratio of {} : {} +/- {}'.format(gasStr, Factor * np.nanmean(dRatioEnh2[indspos][indsDates]), 2.*np.nanstd(dRatioEnh2[indspos][indsDates])/np.sqrt(len(dRatioEnh2[indspos][indsDates])) )

            #print 'Median Emmision factor ratio of {} : {} +/- {}'.format(gasStr, 59.91 *np.nanmedian(dRatioEnh)* (MX[i]/MCO), 59.91 * np.nanstd(dRatioEnh)* (MX[i]/MCO) )
            print 'Median Emmision factor ratio of {} : {} +/- {}'.format(gasStr, EFCO * Factor*np.nanmean(dRatioEnh2[indspos][indsDates])* (MX[i]/MCO), EFCO * 2.*np.nanstd(dRatioEnh2[indspos][indsDates])/np.sqrt(len(dRatioEnh2[indspos][indsDates]))*(MX[i]/MCO) )
            
        if gasStr.upper() == 'HCN': continue
        
        ax = plt.Subplot(fig, outer_grid[i])
  
        #ax.scatter(dates[g],Ratio, color='gray', s=30, facecolors='white', zorder=1, label='All')
        ax.plot(dates[g],Ratio, color='gray', linestyle='None', marker ='.', markersize=10, label='All', zorder=1, alpha=0.75)
        ax.scatter(dates[g],Ratio, color='gray', s=30, facecolors='white', zorder=2, label='All')
        #ax.scatter(dates_Enhan[g],RatioEnh, color='red', s=25, facecolors='white', zorder=1)
        #ax.scatter(dates_Enhan[g],dRatioEnh, color='red', s=35, facecolors='white', zorder=2, label='Enhancements')

        ax.scatter(dates_Enhan2[indspos],dRatioEnh2[indspos], color='red', s=35, facecolors='white', zorder=3, label='Enhancements-2')
        
        ax.grid(True, alpha=0.35)
        ax.tick_params(labelsize=16)
        #ax.grid(True,which='both', alpha=0.35)  
        ax.text(0.03, 0.9, gasStr,  va='center',transform=ax.transAxes,fontsize=24)

        ax.xaxis.set_major_locator(yearsLc)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(DateFmt)

        ax.set_ylim(ymin=0)
        ax.set_xlim(xmin, xmax)

        #if i == 0: 
            #ax.set_title('ER (ppb/ppb)', fontsize=18)        
            #ax.legend(prop={'size':12}, loc = 1)

        #ax.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
        #ax.set_xlim(frappe_i, frappe_f)

        fig.add_subplot(ax)

        #----------------------------
        # 
        #----------------------------
        years = [ singDate.year for singDate in dates_Enhan2 ]

        tcks = range(np.min(years),np.max(years)+2)
        cm   = plt.get_cmap(clmap)
        norm = colors.BoundaryNorm(tcks,cm.N)                       
        #sc1 = ax1.scatter(sza,rms,c=years,cmap=cm,norm=norm)
        #ax2.scatter(sza,dofs,c=years,cmap=cm,norm=norm)   

        ax1 = plt.Subplot(fig2, outer_grid[i])

        dates_Mnth         = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates[g]]
        dates_Enhan_Mnth   = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates_Enhan[g]]
        dates_Enhan_Mnth2  = [dt.datetime(2015, d.month, d.day, d.hour, d.minute, d.second) for d in dates_Enhan2[indspos]]

        ax1 = plt.Subplot(fig2, outer_grid[i])
        
        ax1.plot(dates_Mnth,Ratio, color='gray', linestyle='None', marker ='.', markersize=10, label='All', zorder=1, alpha=0.75)
        #ax1.scatter(dates_Mnth,Ratio, color='gray', s=30, facecolors='white', zorder=1, label='All')
        #ax.scatter(dates_Enhan[g],RatioEnh, color='red', s=25, facecolors='white', zorder=1)
        #ax1.scatter(dates_Enhan_Mnth,dRatioEnh, color='red', s=35, facecolors='white', zorder=2, label='Enhancements')

        #sc1 = ax1.scatter(dates_Enhan_Mnth,dRatioEnh, c=years,cmap=cm,norm=norm, s=35, zorder=2, label='Enhancements')
        #sc1  = ax1.scatter(dates_Enhan_Mnth2,dRatioEnh2, c=years,cmap=cm,norm=norm, s=45,  zorder=2, label='Enhancements-2')
   
        #ax1.scatter(dates_Enhan_Mnth,dRatioEnh, color='red', s=35, facecolors='white', zorder=2, label='Enhancements')
        ax1.scatter(dates_Enhan_Mnth2,dRatioEnh2[indspos], color='blue', s=45, facecolors='white', zorder=2, label='Enhancements-2')
        #ax1.text(0.03, 0.9, gasStr+'/CO', va='center',transform=ax1.transAxes,fontsize=24)
        #ax1.text(0.03, 0.9, '$\Delta $ '+gasStr+'to \Delta $ CO', va='center',transform=ax1.transAxes,fontsize=24)
        #ax1.text(0.03, 0.9, '$\Delta $ '+gasStr+'to \Delta $ CO', va='center',transform=ax1.transAxes,fontsize=24)
        ax1.text(0.03, 0.9, gasStr, va='center',transform=ax1.transAxes,fontsize=24)
        
        ax1.grid(True, alpha=0.35)
        ax1.tick_params(labelsize=16)
        #ax1.grid(True,which='both', alpha=0.35)  
        #ax1.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

        ax1.set_xlim(dt.date(2015,1,1), dt.date(2015,12,31))
        ax1.set_ylim(ymin=0)

        ax1.xaxis.set_major_locator(months)
        #ax.xaxis.set_minor_locator(months)
        ax1.xaxis.set_major_formatter(DateFormatter('%b'))

        #if i == 0: 
        #    ax1.legend(prop={'size':12}, loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, ncol=3) 

        fig2.add_subplot(ax1)

    all_axes   = fig.get_axes()
    all_axes2  = fig2.get_axes()
   
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

    #show only the outside spines
    for ax in all_axes2:
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

            if i == 3:
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)

    
    

    all_axes[-1].set_xlabel('Year', fontsize=16)
    all_axes[-2].set_xlabel('Year', fontsize=16)
    all_axes[-3].set_xlabel('Year', fontsize=16)
    all_axes[-4].set_xlabel('Year', fontsize=16)
    
    all_axes2[-1].set_xlabel('Month', fontsize=16)
    all_axes2[-2].set_xlabel('Month', fontsize=16)

    fig.text(0.02, 0.5,  'ER [ppb/ppb]', fontsize=18, va='center', rotation='vertical')
    fig2.text(0.02, 0.5, 'ER [ppb/ppb]', fontsize=18, va='center', rotation='vertical')

    fig.autofmt_xdate()
    fig2.autofmt_xdate()

    fig.subplots_adjust(bottom=0.1,top = 0.96, left=0.11, right=0.95)
    fig2.subplots_adjust(bottom=0.085,top = 0.96, left=0.11, right=0.95)


    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
    else:
        plt.show(block=False)
    
    fig.savefig(pltDir+'EnhancementRatios_series.pdf', bbox_inches='tight')
    fig2.savefig(pltDir+'EnhancementRatios_mnth.pdf', bbox_inches='tight')


    #----------------------------
    # Whiskers
    #----------------------------

    #-------------------------------------------------------
    #Define parameters for plots
    #-------------------------------------------------------
    clmap        = 'jet'
    #locations    = [ [1,2], [4,5], [7,8], [10,11], [13,14]]
    lticks       = [1, 2, 3, 4, 5, 6]

    locations    = np.array(range(len(gasName2))) + 1
    
    l = 0

    fig,  ax = plt.subplots(figsize=(10, 6))
    gasname_w = []

    for i, g in enumerate(gasName2):

        gasStr = mf.getgasname(g)
    
        Data = np.asarray(dER1[g][0])
        meanpointprops = dict(marker='o', markeredgecolor='red',
                              markerfacecolor='white', markersize=10)


        bp = ax.boxplot(Data,   positions = [locations[i]], widths = 0.65, showmeans=True, meanprops=meanpointprops, patch_artist=True)
        #setBoxColors(bp)

        maxloc = locations[i] + 1.0

        gasname_w.append(gasStr)

        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
            patch.set( linewidth=2)

        for whisker in bp['whiskers']:
            whisker.set(color='blue', linewidth=2)

        ## change color and linewidth of the caps
        for cap in bp['caps']:
            cap.set(color='blue', linewidth=2)

        ## change color and linewidth of the medians
        for median in bp['medians']:
            median.set(linewidth=2)

        maxy = [np.amax(d) for d in Data] 

        #ax.text(lticks[i], np.amax(maxy)+( np.amax(maxy)*0.05), r'N={}'.format(len(Data)), fontsize=18)


        
    ax.set_xlim((0,maxloc))
    ax.set_ylim(ymin=0, ymax=0.2)
    ax.yaxis.grid(True, alpha = 0.25)
    #ax.set_yscale("log", nonposy='clip')
    ax.set_xticklabels(gasname_w)
    ax.set_xticks(locations[0:len(gasName2)])
    ax.xaxis.set_tick_params(which='major',labelsize=18)
    ax.yaxis.set_tick_params(which='major',labelsize=18)
    #ax.set_xlabel('Gas', fontsize = 18)
    ax.set_ylabel('ER [ppb/ppb]', fontsize = 18)


    # draw temporary red and blue lines and use them to create a legend
    # hB, = ax.plot([1,1],'b-')
    # hR, = ax.plot([1,1],'r-')
    # ax.legend((hB, hR),('Cold season', 'Warm season'), prop={'size':16})
    # #legend((hB, hR),(seasons))
    # hB.set_visible(False)
    # hR.set_visible(False)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
    else:
        plt.show(block=False) 

    fig.savefig(pltDir+'Box_ER.pdf', bbox_inches='tight')
    
    user_input = raw_input('Press any key to exit >>> ')
    sys.exit()

    

    #----------------------------
    # Bar: Enhancement ratio (Warm season)
    #----------------------------
    fig,ax = plt.subplots(figsize=(10,7))
    ind    = np.arange(len(yearsList))
    #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))

    #print dER1['nh3']
    #print dERstd1['nh3']

    #print dER1['c2h6']
    #print dERstd1['c2h6']

    #print dER1['h2co']
    #print dERstd1['h2co']
    
    ax.bar(ind-0.27, dER1['c2h6'], 0.27, yerr=dERstd1['c2h6'], align='center', color = 'r', ecolor = 'k', label = 'C$_2$H$_6$')
    ax.bar(ind, dER1['nh3'], 0.27, yerr=dERstd1['nh3'], align='center', color = 'b', ecolor = 'k', label = 'NH$_3$')
    #ax.bar(ind+0.27, dER1['ch4']/10., 0.27, yerr=dERstd1['ch4']/10., align='center', color = 'g', ecolor = 'k', label = 'H$_2$CO')  #yerr=slope_TC*0
    #ax.xaxis.grid(True)
    ax.yaxis.grid(True, alpha=0.35)
    #ax.set_xlabel('Year')
    ax.set_ylabel('Enhancement Ratios (ER)', fontsize=16)
    ax.set_xticks(ind)
    #ax.set_title('Median ER - Warm Season (May to October)', multialignment='center', fontsize=18)

    ax.tick_params(labelsize=16)
    #ax.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_title(subtitle+'\n(1996 - 2016)', multialignment='center')
    #ax.set_xlabel('Site')
    ax.set_xlim(-0.5,len(yearsList)-0.5)
    ax.axvline(0, color='black', lw=1)
    ax.legend(prop={'size':16}, loc = 1)
    #ax.invert_yaxis()
    xtickNames = ax.set_xticklabels(yearsListStr)
    plt.setp(xtickNames, rotation=0, fontsize=16)

    fig.subplots_adjust(bottom=0.075,top=0.95, left=0.11, right=0.95)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
    else:
        plt.show(block=False) 

    fig.savefig(pltDir+'EnhancementRatios_bars.pdf', bbox_inches='tight')

    user_input = raw_input('Press any key to exit >>> ')
    sys.exit()

    #----------------------------
    # Bar: Enhancement ratio (Cold season)
    #----------------------------

    # fig,ax = plt.subplots(figsize=(10,7))
    # ind    = np.arange(len(yearsList))
    # #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    # ax.bar(ind-0.27, dER2['c2h6'], 0.27, yerr=dERstd2['c2h6'], align='center', color = 'r', ecolor = 'k', label = 'C$_2$H$_6$')
    # ax.bar(ind, dER2['nh3'], 0.27, yerr=dERstd2['nh3'], align='center', color = 'b', ecolor = 'k', label = 'NH$_3$')
    # ax.bar(ind+0.27, dER2['h2co'], 0.27, yerr=dERstd2['h2co'], align='center', color = 'g', ecolor = 'k', label = 'H$_2$CO')  #yerr=slope_TC*0
    # ax.xaxis.grid(True)
    # ax.yaxis.grid(True)
    # ax.set_xlabel('Year')
    # ax.set_ylabel('Enhancement Ratio')
    # ax.set_xticks(ind)
    # #ax.set_yticklabels(np.transpose(pltID), rotation=0)
    # ax.set_title('Median ER - Cold Season (November to April)', multialignment='center', fontsize=18)
    # #ax.set_xlabel('Site')
    # ax.set_xlim(-0.5,len(yearsList))
    # ax.axvline(0, color='black', lw=1)
    # ax.legend(prop={'size':12}, loc = 1)
    # #ax.invert_yaxis()
    # xtickNames = ax.set_xticklabels(yearsList)
    # plt.setp(xtickNames, rotation=0, fontsize=11)

    # fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)

    # if saveFlg: 
    #     pdfsav.savefig(fig,dpi=200)
    # else:
    #     plt.show(block=False) 

    
    #----------------------------
    # Bar: Emmission factors (Cold season)
    #----------------------------

    print 'Mean Enhancement ratio of {} : {} +/- {}'.format(gasStr, np.nanmean(dRatioEnh), np.nanstd(dRatioEnh) )

    fig,ax = plt.subplots(figsize=(10,7))
    ind    = np.arange(len(yearsList))
    #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    #ax.bar(ind-0.27, 59.91 * dER['c2h6'] * (28./28.), 0.27, yerr=dERstd['c2h6'], align='center', color = 'r', ecolor = 'k', label = 'C2H6')
    #ax.bar(ind, dER['nh3']*59.91 * (17./23.), 0.27, yerr=dERstd['nh3'], align='center', color = 'b', ecolor = 'k', label = 'NH3')
    #ax.bar(ind+0.27, dER['h2co']*59.91 * (17./23.), 0.27, yerr=dERstd['h2co'], align='center', color = 'g', ecolor = 'k', label = 'H2CO')  #yerr=slope_TC*0

    ax.bar(ind, 59.91 * dER1['c2h6'] * (30./28.), 0.27, yerr=dERstd1['c2h6']* 59.91 *(30./28.), align='center', color = 'r', ecolor = 'k', label = 'C$_2$H$_6$')

    
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlabel('Year')
    ax.set_ylabel('Emmision Factor (EF) (Tons/day)')
    ax.set_xticks(ind)
    ax.set_title('Warm Season', multialignment='center')
    ax.set_title('Median Emmision Factor - Warm Season (May to October)', multialignment='center', fontsize=18)
    #ax.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_title(subtitle+'\n(1996 - 2016)', multialignment='center')
    #ax.set_xlabel('Site')
    ax.set_xlim(-0.5,len(yearsList)-0.5)
    ax.axvline(0, color='black', lw=1)
    ax.legend(prop={'size':12}, loc = 1)
    #ax.invert_yaxis()
    xtickNames = ax.set_xticklabels(yearsList)
    plt.setp(xtickNames, rotation=0, fontsize=11)

    fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
    else:
        plt.show(block=False)  

    #----------------------------
    # Bar: Emmission ratio (Warm season)
    #----------------------------

    # fig,ax = plt.subplots(figsize=(10,7))
    # ind    = np.arange(len(yearsList))
    # #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    # #ax.bar(ind-0.27, 59.91 * dER['c2h6'] * (28./28.), 0.27, yerr=dERstd['c2h6'], align='center', color = 'r', ecolor = 'k', label = 'C2H6')
    # #ax.bar(ind, dER['nh3']*59.91 * (17./23.), 0.27, yerr=dERstd['nh3'], align='center', color = 'b', ecolor = 'k', label = 'NH3')
    # #ax.bar(ind+0.27, dER['h2co']*59.91 * (17./23.), 0.27, yerr=dERstd['h2co'], align='center', color = 'g', ecolor = 'k', label = 'H2CO')  #yerr=slope_TC*0

    # ax.bar(ind, 59.91 * dER2['c2h6'] * (28./28.), 0.27, yerr=dERstd2['c2h6'] * 59.91 *(28./28.) , align='center', color = 'r', ecolor = 'k', label = 'C2H6')
    
    # ax.xaxis.grid(True)
    # ax.yaxis.grid(True)
    # ax.set_xlabel('Year')
    # ax.set_ylabel('Emmision Factor (Tons/day)')
    # ax.set_xticks(ind)
    # ax.set_title('Median Emmision Factor - Cold Season (November to April)', multialignment='center', fontsize=18)
    # #ax.set_yticklabels(np.transpose(pltID), rotation=0)
    # #ax.set_title(subtitle+'\n(1996 - 2016)', multialignment='center')
    # #ax.set_xlabel('Site')
    # ax.set_xlim(-0.5,len(yearsList))
    # ax.axvline(0, color='black', lw=1)
    # ax.legend(prop={'size':12}, loc = 1)
    # #ax.invert_yaxis()
    # xtickNames = ax.set_xticklabels(yearsList)
    # plt.setp(xtickNames, rotation=0, fontsize=11)

    # fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)

    # if saveFlg: 
    #     pdfsav.savefig(fig,dpi=200)
    # else:
    #     plt.show(block=False)  

  

    # #---------------------------------
    # #     PLOT: -- Time Series FRAPPE --
    # #---------------------------------

    
    # for i, g in enumerate(gasName2):

    #     gasStr = mf.getgasname(g)

    #     if g.lower() != 'co':

    #         doy              = mf.toYearFraction(dates[g])   

    #         vmr_x_interpol   = interpolate.interp1d(doy_x, vmr_x, axis=0, fill_value='extrapolate', bounds_error=False)(doy)
            
    #         Ratio            = vmr[g]/vmr_x_interpol

    #         datesEnh         = np.asarray(dates_Enhan[g])
    #         vmrEnh           = np.asarray(vmr_Enhan[g])
    #         dvmrEnh           = np.asarray(dvmr_Enhan[g])
    #         doyEnh           = mf.toYearFraction(datesEnh) 

    #         vmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, vmrEnh_x, axis=0, bounds_error=False)(doyEnh)

    #         RatioEnh            = vmrEnh/vmrEnh_x_interpol

    #         dvmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, dvmrEnh_x, axis=0, bounds_error=False)(doyEnh)

    #         dRatioEnh            = dvmrEnh/dvmrEnh_x_interpol


    #         #fig , (ax, ax1, ax2) = plt.subplots(3, figsize=(10,8), gridspec_kw = {'height_ratios':[3,1,1]}, sharex=True)
            
    #         ydoysc          = mf.toYearFraction(SCdate)
    #         SCNE_interp     = interpolate.interp1d(ydoysc, SCNE, axis=0, bounds_error=False)(doy)
    #         SCNW_interp     = interpolate.interp1d(ydoysc, SCNW, axis=0, bounds_error=False)(doy)
    #         SCSE_interp     = interpolate.interp1d(ydoysc, SCSE, axis=0, bounds_error=False)(doy)
    #         SCSW_interp     = interpolate.interp1d(ydoysc, SCSW, axis=0, bounds_error=False)(doy)
    #         SCwEST          =  SCSW_interp# + SCNW_interp
    #         SCNE_interp     = SCNE_interp + SCNW_interp

    #         ydoysc          = mf.toYearFraction(dt_weather)
    #         wdir_interp     = interpolate.interp1d(ydoysc, wdir, axis=0, bounds_error=False)(doy)


    #         inds_wd         = np.logical_and(wdir_interp >= 0, wdir_interp <= 90.)
    #         wdirNE          = wdir_interp[inds_wd]
    #         wVMRNE          = vmr[g][inds_wd]
    #         datesNE         = dates[g][inds_wd]
           
    #         #----------------------------
    #         #
    #         #----------------------------    
    #         fig , (ax, ax2) = plt.subplots(2, figsize=(10,8), gridspec_kw = {'height_ratios':[3,1]}, sharex=True)

    #         ax.scatter(dates[g],Ratio, color='gray', s=25, facecolors='white', zorder=1)
    #         ax.scatter(dates_Enhan[g],dRatioEnh, color='green', s=25, facecolors='white', zorder=1)
            
    #         ax.grid(True, alpha=0.35)
    #         ax.tick_params(labelsize=14)
    #         ax.grid(True,which='both', alpha=0.35)  
    #         ax.text(0.03, 0.9, gasStr+'/CO', va='center',transform=ax.transAxes,fontsize=24)
    #         ax.tick_params(labelsize=16)
    #         ax.set_ylim(ymin=0)

        
    #         ax2.bar(dates[g], SCNE_interp, width=0.2,  align='center', color = 'r', label='North-East', alpha=1 )
    #         ax2.bar(dates[g], SCSE_interp,width=0.2, align='center', color = 'blue',  label='South-East', alpha=1,  bottom =sumzip(SCNE_interp))
    #         ax2.bar(dates[g], SCwEST,width=0.2, align='center', color = 'green', label='West', alpha=1, bottom =sumzip(SCNE_interp,SCSE_interp))

    #         ax2.grid(True, alpha=0.35)
    #         ax2.tick_params(labelsize=14)
    #         ax2.xaxis.set_major_locator(mondays)
    #         ax2.xaxis.set_minor_locator(dayLc)
    #         ax2.xaxis.set_major_formatter(DateFormatter('%b %d'))

    #         #ax.axvspan(frappe_i, frappe_f,  alpha=0.5, color='yellow')
    #         ax2.set_xlim(frappe_i, frappe_f)
           
    #         fig.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')

    #         fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.095, right=0.96)

    #     if saveFlg: 
    #         pdfsav.savefig(fig,dpi=200)
    #         #plt.savefig(pltDir+'Time_Series_FRAPPE.pdf', bbox_inches='tight')
    #     else:
    #         plt.show(block=False)

    #---------------------------------
    #     PLOT: -- Correlation --
    #---------------------------------

    dates_x      = np.asarray(dates['co'])
    vmr_x        = np.asarray(vmr['co'])
    doy_x        = mf.toYearFraction(dates_x)

    datesEnh_x   = np.asarray(dates_Enhan['co'])
    vmrEnh_x     = np.asarray(vmr_Enhan['co'])
    doyEnh_x     = mf.toYearFraction(datesEnh_x) 

    dvmrEnh_x     = np.asarray(dvmr_Enhan['co'])

    fig  = plt.figure(figsize=(7,10))
    #fig2  = plt.figure(figsize=(7,10))

    outer_grid = gridspec.GridSpec((ngases-1), 1, wspace=0.15, hspace=0.125)

    gasName2 = [g for g in gasName if g.lower() != 'co']

    for i, g in enumerate(gasName2):

        gasStr = mf.getgasname(g)

        if g.lower() != 'co':

            doy              = mf.toYearFraction(dates[g])   

            vmr_x_interpol   = interpolate.interp1d(doy_x, vmr_x, axis=0, fill_value='extrapolate', bounds_error=False)(doy)
          
            Ratio            = vmr[g]/vmr_x_interpol

            datesEnh         = np.asarray(dates_Enhan[g])
            vmrEnh           = np.asarray(vmr_Enhan[g])
            dvmrEnh           = np.asarray(dvmr_Enhan[g])
            doyEnh           = mf.toYearFraction(datesEnh) 

            vmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, vmrEnh_x, axis=0, bounds_error=False)(doyEnh)

            RatioEnh            = vmrEnh/vmrEnh_x_interpol

            dvmrEnh_x_interpol   = interpolate.interp1d(doyEnh_x, dvmrEnh_x, axis=0, bounds_error=False)(doyEnh)

            dRatioEnh            = dvmrEnh/dvmrEnh_x_interpol


            #----------------------------
            #
            #----------------------------
            ax = plt.Subplot(fig, outer_grid[i])

            yearsall  = [ singDate.year for singDate in dates_Enhan[g]]
            yearsall  = np.asarray(yearsall)
            yearsList = list(set(yearsall))
            yearsList.sort()

            clr = mf.clrplt()

            for i, y in enumerate(yearsList):
                indsy = np.where(yearsall == y)[0]

            #     odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx_e, yerr=yy_e, InError=True)
            #     slope = float(odr[0])
            #     intercept = float(odr[1])

                ax.scatter(dvmrEnh_x_interpol[indsy],dvmrEnh[indsy], s=25, color=clr[i], zorder=1, label = y)



            #ax.scatter(vmr_x_interpol,vmr[g], color='gray', s=25, facecolors='white', zorder=1)
            #ax.scatter(vmrEnh_x_interpol,vmrEnh, color='red', s=25, facecolors='white', zorder=1)
            #ax.scatter(dvmrEnh_x_interpol,dvmrEnh, color='green', s=25, facecolors='white', zorder=1)
            if i == 0:  ax.set_title('Yearly correlation of X to CO', fontsize=18)        

            ax.legend(prop={'size':9})

            ax.grid(True, alpha=0.35)
            ax.tick_params(labelsize=16)
            ax.grid(True,which='both', alpha=0.35)  
            ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

            fig.add_subplot(ax)

            #----------------------------
            #
            #----------------------------



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


    all_axes[-1].set_xlabel('CO VMR [ppb]', fontsize=16)

    fig.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')
    fig.subplots_adjust(bottom=0.065,top = 0.94, left=0.15, right=0.95)


    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        #plt.savefig(pltDir+'Time_Series.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)
    
    if saveFlg:     
        pdfsav.close()
    else:           
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()    
   
     # Exit program       


 
if __name__ == "__main__":
    main()