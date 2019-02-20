#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltFTS.py
#
# Purpose:
#         Plot time series of multiple species retrieved with FTIR columns/VMRs
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
from datetime import timedelta
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
    
    gasName            = ['co',           'c2h2',         'c2h6',         'ch4',         'nh3',          'hcn',         'h2co',         'hcooh']         # 'o3',
    ver                = ['Current_v3',   'Current_v2',   'Current_v2',   'Current_WP',  'Current_v2',   'Current_v1', 'Current_v8',   'Current_v3']    # 'Current_WP',  

    #gasName            = ['co',        'h2co' ,        'c2h2',       'hcooh',         'c2h6',         'nh3',         'ch4',                 'hcn'        ,        ]         # 'o3',
    #ver                = ['Current_v3', 'Current_v8',  'Current_v2',  'Current_v3',   'Current_v2',   'Current_v2',  'Current_WP',   'Current_v1'  ]    # 'Current_WP',     

    #loc                = 'tab'
    #gasName            = ['co',            'c2h6',        'hcn']   
    #ver                = ['Current_B3',    'Current_v2',   'Current_WP']       

    saveFlg           = True 
    camFlg            = True

    DirData           = '/data/iortega/results/'+loc.lower()+'/data/'
    #DirData           = '/data/iortega/results/tab/data/'
    pltDir            =  '/data/iortega/pbin/tropGases/fig/'

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2017
    fmnth              = 12        # FRANCO ET AL MIGHT BE 2014 07 31
    fday               = 31


    if camFlg: pltFile  =  pltDir + 'pltFTSCAM-FINN.pdf'
    else: pltFile  =  pltDir + 'pltFTS.pdf'
    


                                    #----------------------------#
                                    #                            #
                                    #        --- START ---       #
                                    #                            #
                                    #----------------------------#

    if saveFlg: pdfsav = PdfPages(pltFile)

    gasFile   = [gasName[i]+'_'+ver[i] for i in range(len(gasName))]   


    #---------------------------------
    #     -- Read CAM-CHEM --
    #---------------------------------
    if camFlg:
        cam = {}

        for gi, gf in enumerate(gasName):
            #try:
                cols, indexToName = mf.getColumns(DirData + gf + '_CAM-CHEM_v2.dat', headerrow=0, delim=',', header=True)

                for index in indexToName:
                    cam[gasName[gi]+'_'+indexToName[index]]  = cols[indexToName[index]]
            #except:
            #    print gf + ' is not in CAM-CHEM'
            #    continue

        #---------------------------------
        #Add Time Stamp - CAM
        #---------------------------------
        for g in gasName:

            try:
            
                date    = cam[g+'_YYYY-MM-DD']
            
                yyyy = [int(d[0:4]) for d in date]
                mm   = [int(d[5:7]) for d in date]
                dd   = [int(d[8:10]) for d in date]

                dateTime =   [dt.date(yyyy[i], mm[i], dd[i]) for i, d in enumerate(dd)]

                cam[g+'_dt'] = dateTime

            except: pass

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------
    fts = {}

    for gi, gf in enumerate(gasFile):
        cols, indexToName = mf.getColumns(DirData + gf + '.dat', headerrow=0, delim=',', header=True)

        for index in indexToName:
            fts[gasName[gi]+'_'+indexToName[index]]  = cols[indexToName[index]]

    #---------------------------------
    #Add Time Stamp - FTS
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

    clr = mf.clrplt()

    xmin      = dt.datetime(iyear, imnth, iday)
    xmax      = dt.datetime(fyear, fmnth, fday)

    #---------------------------------
    #     Analysis - Variables FTS
    #---------------------------------
    dates                 = {}
    vmr                   = {}

    driftFourier          = {}
    drift                 = {}

    slope                 = []
    slope_e               = []

    vmr_mnth              = {}
    vmr_std_mnth          = {}
    dates_mnth            = {}

    driftFourier_mnth     = {}
    drift_mnth            = {}

    driftFourier_c        = {}
    dates_c               = {}

    driftFourier_c_mnth   = {}
    dates_c_mnth          = {}
    
    slope_mnth            = []
    slope_mnth_e          = []

    GasStrall             = []



    #---------------------------------
    #     Analysis - Variables CAM-CHEM
    #---------------------------------

    if camFlg:
        vmr_CAM              = {}
        vmr_std_CAM          = {}
        dates_CAM            = {}

        vmr_CAM_int          = {}
        vmr_std_CAM_int      = {}

        driftFourier_CAM     = {}
        drift_CAM            = {}

        slope_CAM            = []
        slope_e_CAM          = []

        rd                   = []
        rd_std               = []

    #---------------------------------

    
    for i, g in enumerate(gasName):

        idates        = np.asarray(fts[g+'_dt'])

        indsDates       = np.where((idates >= xmin) & (idates < xmax))[0]

        dates[g]        = np.asarray(fts[g+'_dt'])[indsDates]
        #vmr[g]          = np.asarray(fts[g+'_TC'], dtype=np.float32)[indsDates]
        vmr[g]          = np.asarray(fts[g+'_wVMR'], dtype=np.float32)[indsDates]

        gasStr = mf.getgasname(g)
        GasStrall.append(gasStr)

        #----------------------------
        # Drift & Fourier - ALL
        #----------------------------
        dateYearFrac        = mf.toYearFraction(dates[g])
        weights             = np.ones_like(dateYearFrac)
        res                 = mf.fit_driftfourier(dateYearFrac, vmr[g], weights, 2, half_period=1)
        f_drift, f_fourier, f_driftfourier, res_std, A, df_drift = res[3:9]
 
        driftFourier[g]     = f_driftfourier(dateYearFrac)
        drift[g]            = f_drift(dateYearFrac)
        
        res_b               = mf.cf_driftfourier(dateYearFrac, vmr[g], weights, 2, half_period=1)
        perc, intercept_b, slope_b, pfourier_b = res_b

        print "\nRate of Change ({}) = {:.2f} +/- {:.3f}%) - ALL".format(g, np.mean(slope_b)/np.mean(vmr[g])*100.0, np.std(slope_b)/np.mean(vmr[g])*100.0)

        slope.append(np.mean(slope_b)/np.mean(vmr[g])*100.0)
        slope_e.append(np.std(slope_b)/np.mean(vmr[g])*100.0)

        #---------------------------------------------------
        #To make a continuous fit
        #---------------------------------------------------
        numdays = (dates[g].max() + dt.timedelta(days=1) - dates[g].min()).days
        dates2  = [dates[g].min() + dt.timedelta(days=x) for x in range(0, numdays)]
        dates2  = np.asarray(dates2)
        dateYearFrac2 = mf.toYearFraction(dates2)

        dates_c[g]    = dates2
        driftFourier_c[g] = f_driftfourier(dateYearFrac2)

        #----------------------------
        # Drift & Fourier - MONTHLY
        #----------------------------

        #mnthlyVals          = mf.ndaysAvg(vmr[g], dates[g], dateAxis=1, meanAxis=0, deltaDays = dt.timedelta(days=1))
        mnthlyVals          = mf.mnthlyAvg(vmr[g], dates[g], dateAxis=1, meanAxis=0)
        vmr_mnth[g]         = mnthlyVals['mnthlyAvg']
        dates_mnth[g]       = mnthlyVals['dates']
        vmr_std_mnth[g]     = mnthlyVals['std']
        
        years               = np.asarray([dt.date(d.year,1,1) for d in mnthlyVals['dates']])

        dateYearFrac        = mf.toYearFraction(dates_mnth[g])
        weights             = np.ones_like(dateYearFrac)
        res                 = mf.fit_driftfourier(dateYearFrac, vmr_mnth[g], weights, 2, half_period=1)
        f_drift, f_fourier, f_driftfourier = res[3:6]
 
        driftFourier_mnth[g] = f_driftfourier(dateYearFrac)
        drift_mnth[g]        = f_drift(dateYearFrac)

        res_b               = mf.cf_driftfourier(dateYearFrac, vmr_mnth[g], weights, 2, half_period=1)
        perc, intercept_b, slope_b, pfourier_b = res_b
        
        print "Rate of Change ({}) = {:.2f} +/- {:.3f}%) - MONTHLY".format(g, np.mean(slope_b)/np.mean(vmr_mnth[g])*100.0, np.std(slope_b)/np.mean(vmr_mnth[g])*100.0)

        slope_mnth.append(np.mean(slope_b)/np.mean(vmr_mnth[g])*100.0)
        slope_mnth_e.append(np.std(slope_b)/np.mean(vmr_mnth[g])*100.0)

        #---------------------------------------------------
        #To make a continuous fit
        #---------------------------------------------------

        numdays = (dates_mnth[g].max() + dt.timedelta(days=1) - dates_mnth[g].min()).days
        dates2  = [dates_mnth[g].min() + dt.timedelta(days=x) for x in range(0, numdays)]
        dates2  = np.asarray(dates2)
        dateYearFrac2 = mf.toYearFraction(dates2)

        dates_c_mnth[g]    = dates2
        driftFourier_c_mnth[g] = f_driftfourier(dateYearFrac2)


        Amp   = np.sum(res[2]**2)
        Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0

        print 'Amplitude of {}: {}'.format(g, Amp)

        
        #---------------------------------
        # CAM-CHEM
        #---------------------------------
        if camFlg:
            #try:
            dates_CAM[g]        = np.asarray(cam[g+'_dt'])
            ##vmr[g]          = np.asarray(fts[g+'_TC'], dtype=np.float32)
            
            #vmr_CAM[g]          = np.asarray(cam[g+'_wVMR'], dtype=np.float32)
            #vmr_std_CAM[g]      = np.asarray(cam[g+'_wVMR'], dtype=np.float32)

            vmr_CAM[g]          = np.asarray(cam[g+'_wVMRsmth'], dtype=np.float32)
            vmr_std_CAM[g]      = np.asarray(cam[g+'_wVMRsmth_std'], dtype=np.float32)

            #----------------------------
            # Drift & Fourier - ALL
            #----------------------------
            dateYearFrac        = mf.toYearFraction(dates_CAM[g])
            weights             = np.ones_like(dateYearFrac)
            res                 = mf.fit_driftfourier(dateYearFrac, vmr_CAM[g], weights, 2, half_period=1)
            f_drift, f_fourier, f_driftfourier, res_std, A, df_drift = res[3:9]
            
            driftFourier_CAM[g]     = f_driftfourier(dateYearFrac)
            drift_CAM[g]            = f_drift(dateYearFrac)
            
            res_b               = mf.cf_driftfourier(dateYearFrac, vmr_CAM[g], weights, 2, half_period=1)
            perc, intercept_b, slope_b, pfourier_b = res_b

            print "Rate of Change CAM-CHEM ({}) = {:.2f} +/- {:.3f}%) - ALL".format(g, res[1]/np.mean(vmr_CAM[g])*100.0, np.std(slope_b)/np.mean(vmr_CAM[g])*100.0)

            slope_CAM.append(res[1]/np.mean(vmr_CAM[g])*100.0)
            slope_e_CAM.append(np.std(slope_b)/np.mean(vmr_CAM[g])*100.0)
            #except:
            #    slope_CAM.append(np.nan)
            #    slope_e_CAM.append(np.nan)

            doy_ftir             = mf.toYearFraction(dates_mnth[g]) 
            doy_cam              = mf.toYearFraction(dates_CAM[g])   

            vmr_CAM_int[g]       = interpolate.interp1d(doy_cam, vmr_CAM[g], axis=0, bounds_error=False)(doy_ftir)
            vmr_std_CAM_int[g]   = interpolate.interp1d(doy_cam, vmr_std_CAM[g], axis=0, bounds_error=False)(doy_ftir)

            diff                = (vmr_mnth[g] - vmr_CAM_int[g])/vmr_mnth[g]  * 100.
            rd.append(np.nanmean(diff))
            rd_std.append(np.nanstd(diff))

            print "Relative difference in percent of ({}) = {:.2f} +/- {:.3f}%) - ALL".format(g, np.nanmean(diff), np.nanstd(diff))

          
    #---------------------------------
    #     PLOT: -- Time Series - ALL --
    #---------------------------------

    # clmap     = 'jet'

    # fig = plt.figure(figsize=(12,12))
    
    # if ngases > 4: outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)
    # else: outer_grid = gridspec.GridSpec(ngases, 1, wspace=0.15, hspace=0.125)

    # for i, g in enumerate(gasName):

    #     gasStr = mf.getgasname(g)

    #     #----------------------------
    #     #
    #     #----------------------------
    #     ax = plt.Subplot(fig, outer_grid[i])

    #     ax.plot(dates[g],vmr[g], color='gray', linestyle='None', marker ='.', markersize=4, label='All data', zorder=1)
    #     ax.scatter(dates[g],vmr[g], color='r', s=25, facecolors='white', zorder=2)
    #     ax.plot(dates_c[g],driftFourier_c[g], color='green', linewidth=2.5, zorder=3)
    #     ax.plot(dates[g],drift[g], color='blue', linewidth=2.5, zorder=3)

    #     try:
    #         ax.plot(dates_CAM[g],vmr_CAM[g], linewidth = 4.0, color='orange', label='CAM-Chem')
    #     except: pass               
        
    #     ax.grid(True, alpha=0.35)
    #     ax.tick_params(labelsize=16)
    #     ax.grid(True,which='both', alpha=0.35)  
    #     ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

    #     #if g.lower() != 'ch4' or g.lower() != 'co': ax.set_ylim(ymin=0) 
    #     if g.lower() != 'ch4': ax.set_ylim(ymin=0) 

    #     ax.xaxis.set_major_locator(yearsLc)
    #     ax.xaxis.set_minor_locator(months)
    #     ax.xaxis.set_major_formatter(DateFmt)

    #     ax.set_xlim(xmin, xmax)

    #     fig.add_subplot(ax)

    # all_axes  = fig.get_axes()

    # #show only the outside spines
    # for ax in all_axes:
    #     for sp in ax.spines.values():
    #         sp.set_visible(False)
    #         plt.setp(ax.get_xticklabels(), visible=False)
    #         ax.spines['top'].set_visible(True)
    #         ax.spines['left'].set_visible(True)
    #         ax.spines['right'].set_visible(True)
    #         ax.spines['bottom'].set_visible(True)
    #         if ax.is_last_row():
    #             ax.spines['bottom'].set_visible(True)
    #             plt.setp(ax.get_xticklabels(), visible=True)


    # all_axes[-1].set_xlabel('Year', fontsize=16)
    # all_axes[-2].set_xlabel('Year', fontsize=16)

    # fig.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')
    # #fig.text(0.02, 0.5, 'Total Column [molec/cm$^2$]', fontsize=18, va='center', rotation='vertical')

    # fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.095, right=0.96)
    # #fig.autofmt_xdate()

    # if saveFlg: 
    #     pdfsav.savefig(fig,dpi=200)
    #     #plt.savefig(pltDir+'Time_Series_all.pdf', bbox_inches='tight')
    # else:
    #     plt.show(block=False)

    #---------------------------------
    #     PLOT: -- Time Series MONTHLY --
    #---------------------------------
    clmap     = 'jet'

    fig = plt.figure(figsize=(14,12))

    if ngases > 4: outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)
    else: outer_grid = gridspec.GridSpec(ngases, 1, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName):

        gasStr = mf.getgasname(g)

        #----------------------------
        #
        #----------------------------
        ax = plt.Subplot(fig, outer_grid[i])

        p1,  = ax.plot(dates[g], vmr[g], color='gray', linestyle='None', marker ='.', markersize=6, label='FTIR - all', zorder=1)
        ax.errorbar(dates_mnth[g], vmr_mnth[g], fmt='o', yerr=vmr_std_mnth[g] ,markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)
        p3  = ax.scatter(dates_mnth[g],vmr_mnth[g], color='red', s=45, facecolors='white', linewidths=2, zorder=2, label='FTIR - Monthly')
        p4,  = ax.plot(dates_c_mnth[g],driftFourier_c_mnth[g], color='green', linewidth=3, zorder=3, label='FTIR - Seasonal + trend fit')
        p5,  = ax.plot(dates_mnth[g],drift_mnth[g], color='blue', linewidth=2.5, zorder=3, label = 'FTIR - Trend')

        try:
            #if g.lower() == 'c2h6': f = 2.5
            #else: f = 1.
            #ax.plot(dates_CAM[g],vmr_CAM[g]*f, linewidth = 3.0, color='orange',zorder=2)
            ax.errorbar(dates_CAM[g],vmr_CAM[g], fmt='o', yerr=vmr_std_CAM[g] ,markersize=6, color='orange', ecolor='orange', elinewidth=1.5, zorder=1)
            p6  = ax.scatter(dates_CAM[g],vmr_CAM[g], color='orange', s=30, facecolors='white', linewidths=2, zorder=3, label='CAM-Chem - Monthly')
            
        except: pass  

        # Create a legend for the first line.
        if i == 0: 
            first_legend  = fig.legend(handles=[p1, p3], prop={'size':12},loc='center left',bbox_to_anchor=(0.03,0.12), bbox_transform=ax.transAxes)
            #plt.gca().add_artist(first_legend)

            if camFlg: second_legend = fig.legend(handles=[p6], prop={'size':12},loc='center left',bbox_to_anchor=(0.45,0.1), bbox_transform=ax.transAxes)

            #fig.legend(loc='upper right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes, title='fig.legend\nax.transAxes')

        if i ==2:    #second_legend = ax.legend(handles=[p6], prop={'size':10},loc='center left',bbox_to_anchor=(0.1,0.1))
            third_legend = fig.legend(handles=[p4, p5], prop={'size':12},loc='center left',bbox_to_anchor=(0.41,1), bbox_transform=ax.transAxes)
        #if i ==0: ax.legend(loc='center left',bbox_to_anchor=(0.8,0.7),prop={'size':10})
        
        ax.grid(True, alpha=0.35)
        ax.tick_params(labelsize=20)
        ax.grid(True,which='both', alpha=0.35)  
        ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=30)

        if g.lower() != 'ch4': ax.set_ylim(ymin=0) 

        ax.xaxis.set_major_locator(yearsLc)
        ax.xaxis.set_minor_locator(months)
        ax.xaxis.set_major_formatter(DateFmt)

        ax.set_xlim(xmin, xmax)

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


    all_axes[-1].set_xlabel('Year', fontsize=20)
    all_axes[-2].set_xlabel('Year', fontsize=20)

    fig.text(0.02, 0.5, 'wVMR [ppb]', fontsize=18, va='center', rotation='vertical')

    fig.autofmt_xdate()
    fig.subplots_adjust(bottom=0.085,top = 0.98, left=0.095, right=0.96)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        plt.savefig(pltDir+'Time_Series_mnth.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)
    


    #---------------------------------------------------
    # PLOT: --  CORRELATION--
    #---------------------------------------------------

    if camFlg:
        fig = plt.figure(figsize=(12,12))

        outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.2, hspace=0.2)

        for i, g in enumerate(gasName):

            xx   = np.asarray(vmr_mnth[g][~np.isnan(vmr_CAM_int[g])], dtype=np.float32)
            yy   = np.asarray(vmr_CAM_int[g][~np.isnan(vmr_CAM_int[g])], dtype=np.float32)

            xx_e = np.asarray(vmr_std_mnth[g][~np.isnan(vmr_CAM_int[g])], dtype=np.float32)
            yy_e = np.asarray(vmr_std_CAM_int[g][~np.isnan(vmr_CAM_int[g])], dtype=np.float32)


            odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx*0.1, yerr=yy_e, InError=True)
            slope_i      = float(odr[0])
            intercept_i = float(odr[1])

            slope_i_e     = float(odrErr[0])
            intercept_i_e = float(odrErr[1])

            slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(xx, yy)
            rvalue_i = float(r_valueln)

            bias = mf.bias(xx, yy)

            print '\n', g
            print 'Slope: {0:.3f} +/- {1:.3f}'.format(float(odr[0]), float(odrErr[0]))
            print 'Intercept = {0:.3f} +/- {1:.3f}'.format(float(odr[1]), float(odrErr[1]))
            print 'R value = {0:.2f}'.format(float(r_valueln))
            print 'Bias [%] = {0:.3f}'.format(float(bias)/np.mean(xx) * 100.)



            gasname = mf.getgasname(g)

            ax = plt.Subplot(fig, outer_grid[i])

            ax.errorbar(vmr_mnth[g], vmr_CAM_int[g], xerr=vmr_std_mnth[g], yerr=vmr_std_CAM_int[g],fmt='o', markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)
               
            ax.grid(True, alpha=0.35)
            
            #if i ==0: ax.legend(loc='center left',bbox_to_anchor=(0.8,0.75),prop={'size':10})
            ax.tick_params(labelsize=16)
            #ax.set_xlim((0.8, 12.2))

            ax.text(0.03, 0.9, gasname, va='center',transform=ax.transAxes,fontsize=24)

            ax.text(0.25,0.9,"Slope: {0:.3f}".format(slope_i),transform=ax.transAxes,  fontsize=11, color='k')
            ax.text(0.25,0.82,"Intercept: {:.3f}".format(intercept_i),transform=ax.transAxes,  fontsize=11, color='k')
            ax.text(0.25,0.74,"r-value: {0:.2f}".format(rvalue_i),transform=ax.transAxes,  fontsize=11, color='k')
            ax.text(0.25,0.68,"Bias [%]: {0:.3f}".format(float(bias)/np.mean(xx) * 100.),transform=ax.transAxes,  fontsize=11, color='k')


            if g.lower() != 'ch4':
                ax.set_xlim(xmin=0, xmax=np.max(vmr_mnth[g]) + np.max(vmr_mnth[g]) *0.3)
                ax.set_ylim(ymin=0, ymax=np.max(vmr_mnth[g]) + np.max(vmr_mnth[g]) *0.3)
            else:
                ax.set_xlim(xmin=np.min(vmr_mnth[g]) - np.min(vmr_mnth[g]) *0.05, xmax=np.max(vmr_mnth[g]) + np.max(vmr_mnth[g]) *0.05)
                ax.set_ylim(ymin=np.min(vmr_mnth[g]) - np.min(vmr_mnth[g]) *0.05, ymax=np.max(vmr_mnth[g]) + np.max(vmr_mnth[g]) *0.05)
            
            ymin, ymax = ax.get_ylim()
     
            one2one = np.arange(ymin, ymax, 0.1)

            ax.plot(one2one, one2one, color='black', linestyle='--', linewidth=2.0)

            fig.add_subplot(ax)

        all_axes  = fig.get_axes()

        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=True)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)


        all_axes[-1].set_xlabel('wVMR - FTIR [ppb]', fontsize=16)
        all_axes[-2].set_xlabel('wVMR - FTIR [ppb]', fontsize=16)

        fig.text(0.02, 0.5, 'wVMR - CAM-Chem [ppb]', fontsize=18, va='center', rotation='vertical')

        fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.085, right=0.98)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(pltDir+'Correlation_FTSvsCAM.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()           # Exit program



    # #---------------------------------
    # #     PLOT: -- Time Series MONTHLY -- FOR TAB (REBECCA'S PROPOSAL)
    # #---------------------------------
    # xmin      = dt.date(iyear, imnth, iday)
    # xmax      = dt.date(fyear, fmnth, fday)
    # clmap     = 'jet'

    # fig = plt.figure(figsize=(12,12))

    # if ngases > 4: outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)
    # else: outer_grid = gridspec.GridSpec(ngases, 1, wspace=0.15, hspace=0.125)

    # for i, g in enumerate(gasName):

    #     gasStr = mf.getgasname(g)

    #     #----------------------------
    #     #
    #     #----------------------------
    #     ax = plt.Subplot(fig, outer_grid[i])

    #     ax.plot(dates[g], vmr[g], color='gray', linestyle='None', marker ='.', markersize=6, label='All data', zorder=1)
    #     ax.errorbar(dates_mnth[g], vmr_mnth[g], fmt='o', yerr=vmr_std_mnth[g] ,markersize=6, color='red', ecolor='red', elinewidth=1.5, zorder=1)
    #     #ax.plot(dates_mnth[g],vmr_mnth[g], color='red', s=45, facecolors='white',  linestyle='-', linewidth = 2.0, label='All data', zorder=2) 

    #     ax.scatter(dates_mnth[g],vmr_mnth[g], color='red', s=55, facecolors='white', linewidths=2, zorder=2)
    #     ax.plot(dates_c[g],driftFourier_c[g], color='green', linewidth=2.5, zorder=3)
    #     #ax.plot(dates_mnth[g],drift_mnth[g], color='blue', linewidth=2.5, zorder=3)
        
    #     ax.grid(True, alpha=0.35)
    #     ax.tick_params(labelsize=16)
    #     ax.grid(True,which='both', alpha=0.35)  
    #     ax.text(0.03, 0.9, gasStr, va='center',transform=ax.transAxes,fontsize=24)

    #     #if g.lower() != 'ch4': ax.set_ylim(ymin=0) 

    #     ax.xaxis.set_major_locator(yearsLc)
    #     ax.xaxis.set_minor_locator(months)
    #     ax.xaxis.set_major_formatter(DateFmt)

    #     ax.set_xlim(xmin, xmax)

    #     fig.add_subplot(ax)

    # all_axes  = fig.get_axes()

    # #show only the outside spines
    # for ax in all_axes:
    #     for sp in ax.spines.values():
    #         sp.set_visible(False)
    #         plt.setp(ax.get_xticklabels(), visible=False)
    #         ax.spines['top'].set_visible(True)
    #         ax.spines['left'].set_visible(True)
    #         ax.spines['right'].set_visible(True)
    #         ax.spines['bottom'].set_visible(True)
    #         if ax.is_last_row():
    #             ax.spines['bottom'].set_visible(True)
    #             plt.setp(ax.get_xticklabels(), visible=True)


    # all_axes[-1].set_xlabel('Year', fontsize=16)
    # #all_axes[-2].set_xlabel('Year', fontsize=16)

    # #fig.text(0.02, 0.5, 'VMR [ppb]', fontsize=18, va='center', rotation='vertical')
    # fig.text(0.02, 0.5, 'Total Column [molec/cm$^2$]', fontsize=18, va='center', rotation='vertical')

    # fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.095, right=0.96)
    # fig.autofmt_xdate()

    # if saveFlg: 
    #     pdfsav.savefig(fig,dpi=200)
    #     plt.savefig(DirData+'Time_Series_mnth.pdf', bbox_inches='tight')
    # else:
    #     plt.show(block=False)
    # plt.savefig(DirData+'Time_Series_mnth.pdf', bbox_inches='tight')

    # user_input = raw_input('Press any key to exit >>> ')
    # sys.exit() 


    #---------------------------------------------------
    # PLOT: --  MONTHLY BY YEAR--
    #---------------------------------------------------
    fig = plt.figure(figsize=(12,12))

    outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName):

        years     = np.asarray([d.year for d in dates_mnth[g]])

        #----------------------------------------
        # Create a list of months to loop through
        #----------------------------------------
        uniqueyears     = list(set(years))          # Find a list of unique months
        uniqueyears.sort()
        uniqueyears     = np.asarray(uniqueyears)

        if camFlg:
            try:
                yearsCAM          = np.asarray([d.year for d in dates_CAM[g]])
                uniqueyearsCAM    = list(set(yearsCAM))          # Find a list of unique months
                uniqueyearsCAM.sort()
                uniqueyearsCAM    = np.asarray(uniqueyearsCAM)
            except: pass

        gasname = mf.getgasname(g)

        ax = plt.Subplot(fig, outer_grid[i])

        for ii, yyyy in enumerate(uniqueyears):

            indsY = where(years == yyyy)[0]
            month = np.array([d.month for d in dates_mnth[g][indsY]])

            ax.errorbar(month, vmr_mnth[g][indsY], yerr=vmr_std_mnth[g][indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy)
            #ax.errorbar(month, vmr_mnth[g][indsY]/np.nanmean(vmr_mnth[g][indsY]), yerr=vmr_std_mnth[g][indsY]/np.nanmean(vmr_mnth[g][indsY]),fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy.year)
            #ax.plot(month, vmr_mnth[g][indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)

            try:
                indsY = where(yearsCAM == yyyy)[0]
                if len(indsY) >= 1:
                    month = np.array([d.month for d in dates_CAM[g][indsY]])
                    ax.plot(month, vmr_CAM[g][indsY], color=clr[ii], linewidth = 2.0)
                    #ax.errorbar(month, vmr_CAM[g][indsY], yerr=vmr_std_CAM[g][indsY],fmt='none', markersize=1, color=clr[ii], ecolor=clr[ii], label=yyyy, elinewidth=2.5)
                    #ax.plot(month, vmr_CAM[g][indsY]/np.nanmean(vmr_CAM[g][indsY]), color=clr[ii], alpha=0.5, linewidth = 2.0)
            except: pass


        ax.grid(True, alpha=0.35)
        
        if i ==0: ax.legend(loc='center left',bbox_to_anchor=(0.8,0.75),prop={'size':10})
        ax.tick_params(labelsize=16)
        ax.set_xlim((0.8, 12.2))

        ax.text(0.03, 0.9, gasname, va='center',transform=ax.transAxes,fontsize=24)

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


    all_axes[-1].set_xlabel('Month', fontsize=16)
    all_axes[-2].set_xlabel('Month', fontsize=16)

    fig.text(0.02, 0.5, 'wVMR [ppb]', fontsize=18, va='center', rotation='vertical')

    fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.085, right=0.98)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        plt.savefig(pltDir+'MonthlyAvg.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)

    #---------------------------------------------------
    # PLOT: --  MONTHLY for selected years YEAR--
    #---------------------------------------------------
    clr2 = ['blue', 'red']
   
    fig = plt.figure(figsize=(12,12))

    outer_grid = gridspec.GridSpec(npanels, 2, wspace=0.15, hspace=0.125)

    for i, g in enumerate(gasName):
        ci = 0

        #years     = np.asarray([dt.date(d.year,1,1) for d in dates_mnth[g]])
        years     = np.asarray([d.year for d in dates_mnth[g]])

        #----------------------------------------
        # Create a list of months to loop through
        #----------------------------------------
        uniqueyears     = list(set(years))          # Find a list of unique months
        uniqueyears.sort()
        uniqueyears     = np.asarray(uniqueyears)

        if camFlg:
            try:
                yearsCAM          = np.asarray([d.year for d in dates_CAM[g]])
                uniqueyearsCAM    = list(set(yearsCAM))          # Find a list of unique months
                uniqueyearsCAM.sort()
                uniqueyearsCAM    = np.asarray(uniqueyearsCAM)
            except: pass

        gasname = mf.getgasname(g)

        ax = plt.Subplot(fig, outer_grid[i])

        clr = mf.clrplt()

        for ii, yyyy in enumerate(uniqueyears):

            indsY = where(years == uniqueyears[ii] )[0]
            month = np.array([d.month for d in dates_mnth[g][indsY]])

            if int(yyyy) == 2014 or int(yyyy) == 2017:

                #ax.errorbar(month, vmr_mnth[g][indsY], yerr=vmr_std_mnth[g][indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy)
                #ax.plot(month, vmr_mnth[g][indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)

                ax.fill_between(month, vmr_mnth[g][indsY] - vmr_std_mnth[g][indsY], vmr_mnth[g][indsY] + vmr_std_mnth[g][indsY], color=clr2[ci], alpha=0.25, zorder=2)
                #ax.plot(month, vmr_mnth[g][indsY] - vmr_std_mnth[g][indsY], color=clr[ii], alpha=0.25, linewidth = 2.0)
                #ax.plot(month, vmr_mnth[g][indsY] + vmr_std_mnth[g][indsY], color=clr[ii], alpha=0.25, linewidth = 2.0)
                
                ax.plot(month, vmr_mnth[g][indsY], color=clr2[ci], alpha=0.5, linewidth = 2.0, zorder=3)
                ax.scatter(month, vmr_mnth[g][indsY],  color=clr2[ci], s=55, facecolors='white', linewidths=2, zorder=4,  label=yyyy)

                try:
                    indsY = where(yearsCAM == yyyy)[0]
                    if len(indsY) >= 1:
                        month = np.array([d.month for d in dates_CAM[g][indsY]])
                        ax.plot(month, vmr_CAM[g][indsY], color=clr2[ci], alpha=0.5, linewidth = 2.0)
                        #ax.plot(month, vmr_CAM[g][indsY]/np.nanmean(vmr_CAM[g][indsY]), color=clr[ii], alpha=0.5, linewidth = 2.0)
                except: pass

                ci+=1

            #else: ax.plot(month, vmr_mnth[g][indsY], color='k', alpha=0.5, linewidth = 2.0, zorder=1)


        ax.grid(True, alpha=0.35)
        
        if i ==0: ax.legend(loc='center left',bbox_to_anchor=(0.8,0.7),prop={'size':10})
        ax.tick_params(labelsize=16)
        ax.set_xlim((0.8, 12.2))

        ax.text(0.03, 0.9, gasname, va='center',transform=ax.transAxes,fontsize=24)

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


    all_axes[-1].set_xlabel('Month', fontsize=16)
    all_axes[-2].set_xlabel('Month', fontsize=16)

    fig.text(0.02, 0.5, 'wVMR [ppb]', fontsize=18, va='center', rotation='vertical')

    fig.subplots_adjust(bottom=0.065,top = 0.98, left=0.085, right=0.98)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        #plt.savefig(pltDir+'MonthlyAvg_years.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)

    #---------------------------------------------------
    # PLOT: --  MONTHLY for selected years and selected gases--
    #---------------------------------------------------
    fig = plt.figure(figsize=(8,8))

    outer_grid = gridspec.GridSpec(2, 1, wspace=0.2, hspace=0.2)

    gasNameOI            = ['c2h6', 'ch4']     

    for i, g in enumerate(gasNameOI):
        ci = 0

        #years     = np.asarray([dt.date(d.year,1,1) for d in dates_mnth[g]])
        years     = np.asarray([d.year for d in dates_mnth[g]])

        #----------------------------------------
        # Create a list of months to loop through
        #----------------------------------------
        uniqueyears     = list(set(years))          # Find a list of unique months
        uniqueyears.sort()
        uniqueyears     = np.asarray(uniqueyears)

        if camFlg:
            try:
                yearsCAM          = np.asarray([d.year for d in dates_CAM[g]])
                uniqueyearsCAM    = list(set(yearsCAM))          # Find a list of unique months
                uniqueyearsCAM.sort()
                uniqueyearsCAM    = np.asarray(uniqueyearsCAM)
            except: pass

        gasname = mf.getgasname(g)

        ax = plt.Subplot(fig, outer_grid[i])

        clr = mf.clrplt()

        for ii, yyyy in enumerate(uniqueyears):

            indsY = where(years == uniqueyears[ii] )[0]
            month = np.array([d.month for d in dates_mnth[g][indsY]])

            if int(yyyy) == 2014 or int(yyyy) == 2017:

                #ax.errorbar(month, vmr_mnth[g][indsY], yerr=vmr_std_mnth[g][indsY],fmt='o', markersize=7.5, color=clr[ii], ecolor=clr[ii], label=yyyy)
                #ax.plot(month, vmr_mnth[g][indsY], color=clr[ii], alpha=0.5, linewidth = 2.0)

                ax.fill_between(month, vmr_mnth[g][indsY] - vmr_std_mnth[g][indsY], vmr_mnth[g][indsY] + vmr_std_mnth[g][indsY], color=clr2[ci], alpha=0.25, zorder=2)
                #ax.plot(month, vmr_mnth[g][indsY] - vmr_std_mnth[g][indsY], color=clr[ii], alpha=0.25, linewidth = 2.0)
                #ax.plot(month, vmr_mnth[g][indsY] + vmr_std_mnth[g][indsY], color=clr[ii], alpha=0.25, linewidth = 2.0)
                
                ax.plot(month, vmr_mnth[g][indsY], color=clr2[ci], alpha=0.5, linewidth = 2.0, zorder=3)
                ax.scatter(month, vmr_mnth[g][indsY],  color=clr2[ci], s=55, facecolors='white', linewidths=2, zorder=4,  label=yyyy)

                try:
                    indsY = where(yearsCAM == yyyy)[0]
                    if len(indsY) >= 1:
                        month = np.array([d.month for d in dates_CAM[g][indsY]])
                        ax.plot(month, vmr_CAM[g][indsY], color=clr2[ci], alpha=0.5, linewidth = 2.0)
                        #ax.plot(month, vmr_CAM[g][indsY]/np.nanmean(vmr_CAM[g][indsY]), color=clr[ii], alpha=0.5, linewidth = 2.0)
                except: pass

                ci += 1

            #else: ax.plot(month, vmr_mnth[g][indsY], color='k', alpha=0.5, linewidth = 2.0, zorder=1)


        ax.grid(True, alpha=0.35)
        
        #if i ==0: ax.legend(loc='center left',bbox_to_anchor=(0.8,0.7),prop={'size':16})
        if i ==0: ax.legend(loc=1,prop={'size':14})
        ax.tick_params(labelsize=16)
        ax.set_xlim((0.8, 12.2))

        ax.text(0.03, 0.9, gasname, va='center',transform=ax.transAxes,fontsize=24)

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


    all_axes[-1].set_xlabel('Month', fontsize=16)
    #all_axes[-2].set_xlabel('Month', fontsize=16)

    fig.text(0.02, 0.5, 'wVMR [ppb]', fontsize=18, va='center', rotation='vertical')

    fig.subplots_adjust(bottom=0.1,top = 0.95, left=0.15, right=0.95)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        plt.savefig(pltDir+'MonthlyAvg_C2H6_CH4.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)



    #---------------------------------
    #     PLOT: -- Box Plot --
    #---------------------------------
    seasons   = ('Winter', 'Summer')

    #-------------------------------------------------------
    #Define parameters for plots
    #-------------------------------------------------------
    
    locations    = [[i, i+1] for i in range(1, ngases*3, 3)]
    lticks       = [float(i)+0.5 for i in range(1, ngases*3, 3)]

    fig,  ax = plt.subplots(figsize=(10, 6))
    gasname_w = []

    for k, g in enumerate(gasName):
        s = mf.getseason(dates[g])
        s = np.asarray(s)

        Data = []
        
        for season in seasons:
            inds   = np.where( s == season )[0]
            Data.append(vmr[g][inds])
        
        if gasName[k].lower() == 'co':
            Data = np.asarray(Data)/100.0
            maxy = [np.amax(d) for d in Data] 
            ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.05), r'x100', fontsize=18)
        
        if gasName[k].lower() == 'ch4': 
            Data = np.asarray(Data)/5000.0
            maxy = [np.amax(d) for d in Data] 
            ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.1), r'x5000', fontsize=18)
            
        if gasName[k].lower() == 'o3': 
            Data = np.asarray(Data)/100.0
            maxy = [np.amax(d) for d in Data] 
            ax.text(lticks[k]-1.0, np.amax(maxy)+( np.amax(maxy)*0.05), r'x100', fontsize=18)
    
        Data = np.asarray(Data)
        meanpointprops = dict(marker='o', markeredgecolor='black',
                              markerfacecolor='white', markersize=4)

        bp = ax.boxplot(Data, positions = [locations[k][0], locations[k][1]], widths = 0.6,  showmeans=True, meanprops=meanpointprops)
        mf.setBoxColors(bp)
        maxloc = locations[k][1] + 1.0

        gasname = mf.getgasname(gasName[k])
        gasname_w.append(gasname)
        
    ax.set_xlim((0,maxloc))
    #ax.set_yscale("log", nonposy='clip')
    ax.set_xticklabels(gasname_w)
    ax.set_xticks(lticks[0:ngases])
    ax.xaxis.set_tick_params(which='major',labelsize=18)
    ax.yaxis.set_tick_params(which='major',labelsize=18)
    ax.set_xlabel('Gas', fontsize = 18)
    ax.set_ylabel('Weighted VMR [ppbv]', fontsize = 18)
    ax.tick_params(labelsize=18)


    # draw temporary red and blue lines and use them to create a legend
    hB, = ax.plot([1,1],'b-')
    hR, = ax.plot([1,1],'r-')
    ax.legend((hB, hR),('Cold season', 'Warm season'), prop={'size':16})
    #legend((hB, hR),(seasons))
    hB.set_visible(False)
    hR.set_visible(False)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        #plt.savefig(pltDir+'Box.pdf', bbox_inches='tight')
    else:
        plt.show(block=False)
    
    #---------------------------------
    #     PLOT: -- RATE OF CHANGE BAR PLOT --
    #---------------------------------

    N   = len(slope)
    ind = np.arange(N)  

    if camFlg:

        fig, (ax, ax2) = plt.subplots(2, figsize=(9, 8), sharex=True)

        ax.bar(ind+(0.27/2.), slope, width = 0.27, color = 'blue', yerr=slope_e, ecolor = 'k', label = 'FTIR - all', align='center')
        ax.bar(ind+(0.27 + 0.27/2.), slope_mnth, width = 0.27, color = 'r', yerr=slope_mnth_e, ecolor = 'k',label = 'FTIR - Monthly', align='center')
        ax.bar(ind+(0.27*2+0.27/2.), slope_CAM, width = 0.27, color = 'orange', yerr=slope_e_CAM, ecolor = 'k', label = 'CAM-Chem - Monthly', align='center')
        #ax.bar(ind+0.4,slope_mnth, width = 0.4, align='center', color = 'b', yerr=slope_mnth_e, ecolor = 'k', label = 'Monthly')
        ax.set_ylabel('Rate of change [%$\cdot$yr$^{-1}$]', fontsize=20)
        #ax.set_xlabel('Gas', fontsize=18)
        #ax.set_xticks(ind+0..27*3/2.)
        #ax.set_xticklabels(GasStrall,  rotation=45)
        ax.tick_params(labelsize=20)
        ax.legend(prop={'size':14}, loc=2)
        ax.axhline(0, color='black', lw=1)
        ax.yaxis.grid(True)
        ax.xaxis.set_tick_params(which='major',labelsize=20)
        ax.yaxis.set_tick_params(which='major',labelsize=20)
        ax.set_xlim(-0.2, N)


        ax2.bar(ind, rd, width = 0.27*3, color = 'r', yerr=rd_std, ecolor = 'k')
        
        ax2.set_ylabel('Relative Difference [%]', fontsize=20)
        #ax.set_xlabel('Gas', fontsize=18)
        ax2.set_xticks(ind+0.27*3/2.)
        ax2.set_xticklabels(GasStrall,  rotation=45)
        ax2.tick_params(labelsize=20)
        ax2.axhline(0, color='black', lw=1)
        ax2.yaxis.grid(True)
        ax2.xaxis.set_tick_params(which='major',labelsize=20)
        ax2.yaxis.set_tick_params(which='major',labelsize=20)
        ax2.set_xlim(-0.2, N)


        fig.subplots_adjust(left=0.135, bottom=0.15, right=0.95, top=0.95)

        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(pltDir+'Trends.pdf', bbox_inches='tight')
        else:           
            plt.show(block=False)

    else:
        fig, ax = plt.subplots(figsize=(10, 6))

        #ax.bar(ind, slope_mnth, width = 0.8, color = 'r', yerr=slope_mnth_e, ecolor = 'k')
        ax.bar(ind, slope, width = 0.8, color = 'r', yerr=slope_e, ecolor = 'k')
        #ax.bar(ind+0.4,slope_mnth, width = 0.4, align='center', color = 'b', yerr=slope_mnth_e, ecolor = 'k', label = 'Monthly')
        ax.set_ylabel('Annual rate of change [%]', fontsize=18)
        #ax.set_xlabel('Gas', fontsize=18)
        ax.set_xticks(ind+0.4)
        ax.set_xticklabels(GasStrall,  rotation=45)
        ax.tick_params(labelsize=18)
        #ax.legend(prop={'size':12}, loc=4)
        ax.axhline(0, color='black', lw=1)
        ax.yaxis.grid(True)
        ax.xaxis.set_tick_params(which='major',labelsize=18)
        ax.yaxis.set_tick_params(which='major',labelsize=18)
        ax.set_xlim(-0.2, N)
        #ax.set_ylim(-15, 7)
        #ax.set_title('Annual rate of change (%)', fontsize=14)
        fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95)

        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
        else:           
            plt.show(block=False)



    # fig, ax = plt.subplots()
    # Ra1 = ax.bar(ind, slope, width = 0.4, align='center', color = 'r', yerr=slope_e, ecolor = 'k', label = 'All')
    # ax.bar(ind+0.4,slope_mnth, width = 0.4, align='center', color = 'b', yerr=slope_mnth_e, ecolor = 'k', label = 'Monthly')
    # # add some text for labels, title and axes ticks
    # ax.set_ylabel('Annual rate of change (%)', fontsize=18)
    # #ax.xticks(ind + 0.35/2., gasName)
    # ax.set_xlabel('Gas', fontsize=18)
    # #ax.set_xticks(ind+0.4)
    # ax.set_xticks(ind+0.35/2.0)
    # ax.set_xticklabels(GasStrall,  rotation=45)
    # ax.tick_params(labelsize=18)
    # ax.legend(prop={'size':12}, loc=4)
    # ax.axhline(0, color='black', lw=1)
    # ax.yaxis.grid(True)
    # ax.xaxis.set_tick_params(which='major',labelsize=18)
    # ax.yaxis.set_tick_params(which='major',labelsize=18)
    # #ax.set_xlim(-0.3, N)
    # #ax.set_ylim(-15, 7)
    # #ax.set_title('Annual rate of change (%)', fontsize=14)
    # fig.subplots_adjust(left=0.125, bottom=0.15, right=0.95, top=0.95)

    # if saveFlg:     pdfsav.savefig(fig,dpi=200)
    # else:           plt.show(block=False)

 

    if saveFlg: pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program  

 
if __name__ == "__main__":
    main()