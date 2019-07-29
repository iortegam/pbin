#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltHDFGases.py
#
# Purpose:
#         Plot time series of multiple species using HDF GEOMS
#
# Notes:  
#         Initially created for OCS (using trop height) and see trend in trop/strat OCS
#   
#
# Version History:
#       Created, Feb, 2017  Ivan Ortega (iortega@ucar.edu)

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
import math

from scipy import interpolate
from scipy.integrate import simps

import matplotlib.dates as md
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
from collections                     import OrderedDict
import PltClass as mp
#import sondeclass as rs
import classOCS as dc
import scipy.stats.mstats as stats
from mpl_toolkits.basemap import Basemap

try:
    from itertools import product
except ImportError:
    # product is new in v 2.6
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)

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



    
                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #-----------------------------------------------------------------------------------------
    #                             Initializations
    #-----------------------------------------------------------------------------------------
    gasName        = 'ocs'

    tcolFlg        = True        # if True it will show/plot only total columns

    #-------------------------------------
    #Global Data Directory
    #-------------------------------------
    GDataDir     = '/data1/projects/ocs/'
    #-------------------------------------
    #three letter ID
    #-------------------------------------
    locs         = ['bld', 'wlg', 'jfj', 'mlo', 'tab', 'tor', 'eur', 'stp', 'ldr', 'rkb', 'tsk', 'ahs', 'mai', 'std', 'zgp', 'kir', 'iza', 'par', 'bre', 'nya', 'pmb', 'alz']

    #-------------------------------------
    #ID in the HDF Files
    #-------------------------------------
    locID        = ['boulder', 'wollongong', 'jungfraujoch', 'mauna.loa.h', 'thule', '_toronto_', '_eureka_', 'st.petersburg', '_laud_120hr', 'rikubetsu', 
                    'tsukuba', 'ahts',  'reunion.maido' , 'stdenis', 'zugspitze', 'kiruna', 'izana', 'paris', 'bremen', 'ny.alesund', 'paramaribo', 'altzomoni']                       
    #-------------------------------------
    #Names in Plots
    #-------------------------------------   
    pltID        = ['Boulder', 'Wollongong', 'Jungfraujoch', 'Mauna Loa', 'Thule', 'Toronto', 'Eureka', 'St Petersburg', 'Lauder', 'Rikubetsu', 
                   'Tsukuba', 'AHTS', 'Maido', 'St Denis', 'Zugspitze', 'Kiruna', 'Izana', 'Paris', 'Bremen', 'Ny Alesund', 'Paramaribo', 'Altzomoni'] 

    
    #------
    # Flags
    #------
    AvgType       = 'Monthly'               # Average type: 'Monthly' or  'Daily'
    smthFlg       = False                   # smooth data?
    period        = 1.0                     # Period for the Fourier Series
    fitFlg        = True                    # Caulcualte Fits?

    #------
    # Flags
    #------
    saveFlg       = False                   # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
    errorFlg      = False                  # Flag to process error data
    fltrFlg       = True                   # Flag to filter data based on dates and negative columns/profiles

    dateFlg       = True                   # Flag to filter based on min and max dates
   
    iyear         = 1985   
    imonth        = 1
    iday          = 1
    fyear         = 2020
    fmonth        = 12
    fday          = 31
    
    sclfct        = 1.0E12                 # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName    = 'ppt'                  # Name of scale factor for labeling plots
    TCsclfct      = 1.0e15
    TCsclfctName  = 'x10$^{15}$'

    pColsFlg      = True                   # Calculate tropospheric and stratospheric columns?

    pltPcol       = False                  # plot the time series in partial columns
    pltWvmr       = True                   # plot the time series in weighted VMR
    
    Adth          = 16.0                   # Altitude in km of tropopause in case NCEP or DTH is not available
    offH          = 5.0                    # Additional altitude above the tropopause height

    #-------------------------------------
    # Flag for Plots
    #-------------------------------------
    mapPltFlg     = True
    prfPltFlg     = True
    aksPltFlg     = True
    errPltFlg     = True 
    mntPltFlg     = True


    yrsROC        = [ [1996, 2018], [1996, 2002], [2002, 2008], [2008, 2015 ], [2015, 2019] ]



                                            #----------------------------#
                                            #        --- START ---       #
                                            #----------------------------#



    pltFile  =  '/data1/projects/ocs/figures/HDF_'+gasName.upper()+'_wVMR'+'_'+AvgType+'_2019.pdf'

    #-------------------------------------
    # Check file and directories
    #-------------------------------------
    dataDir    = [GDataDir+l+'/'  for l in locs]

    for d in dataDir:  ckDir(d,exit=True)

    ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------------------
    # For total columns in Luader subsitute
    #-------------------------------------
    if tcolFlg: 
        locID = ['_laud_' if x == '_laud_120hr' else x for x in locID]

    #-------------------------
    # Create Instance of Class
    #-------------------------
    hdf = dc.pltOCS(gasName, dataDir, locID, pltID, locs, iyear=iyear, imonth=imonth, iday=iday, fyear=fyear, fmonth=fmonth, fday=fday,  saveFlg= saveFlg, outFname=pltFile, errorFlg=errorFlg, pColsFlg=pColsFlg, fltrFlg=fltrFlg)

    #-------------------------
    # Total Columns
    #-------------------------
    if tcolFlg:
        resTC = dc.AnalTS2(hdf.npanels2, hdf.dates, hdf.totClmn, hdf.pltID, hdf.Lat, hdf.locsID, fits=fitFlg, AvgType=AvgType, pltFig=fitFlg, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Total Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=4.5 ,ymax=12)
        resTC = np.asarray(resTC)

        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()
    
    #-------------------------
    # Map
    #-------------------------
    if mapPltFlg: hdf.pltMap()

    #-------------------------
    # Profiles
    #-------------------------
    if prfPltFlg:  hdf.pltPrf()

    #-------------------------
    # Error profiles 
    #-------------------------
    if errPltFlg: 
        if errorFlg: hdf.pltPrfErr()

    #-------------------------
    # AKs
    #-------------------------
    if aksPltFlg: hdf.pltAK(pltWvmr=pltWvmr)


    

    
    
    #---------------------------------------------------
    # Tropospheric Weighted VMR ==> Retrieval
    #---------------------------------------------------
    print '\nPlot: Tropospheric Weighted VMR:\n'

    resvmrLC = dc.AnalTS(hdf.npanels3-1, hdf.dates2, hdf.WvmrTrop12, hdf.pltID2, hdf.Lat2, hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Low Tropospheric wVMR', unitsStr=sclfctName, period=period, ymin=250 ,ymax=650, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrLC = np.asarray(resvmrLC)

    print '\nPlot: Tropospheric Weighted VMR Anomalies:\n' 
    resvmrLCAnom = dc.AnalTSAnom(hdf.npanels3-1, hdf.dates2, hdf.WvmrTrop12, hdf.pltID2, hdf.Lat2,hdf.locID2,  fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Low Tropospheric wVMR Anomalies', unitsStr=sclfctName, period=period, qboFlg=False, ymin=-150 ,ymax=150, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrLCAnom = np.asarray(resvmrLCAnom)

    
    # #---------------------------------------------------
    # # Free Tropospheric Weighted VMR ==> Retrieval
    # #---------------------------------------------------
    print '\nPlot: Free Tropospheric Weighted VMR:\n' 
    resvmrLC2 = dc.AnalTS(hdf.npanels3-1, hdf.dates2, hdf.WvmrTrop22, hdf.pltID2, hdf.Lat2, hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Free Tropospheric wVMR', unitsStr=sclfctName, period=period, ymin=350, ymax=600, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrLC2 = np.asarray(resvmrLC2)

    print '\nPlot: Free Tropospheric Weighted VMR Anomalies:\n' 
    resvmrLCAnom2 = dc.AnalTSAnom(hdf.npanels3-1, hdf.dates2, hdf.WvmrTrop22, hdf.pltID2, hdf.Lat2,hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Free Tropospheric wVMR Anomalies', unitsStr=sclfctName, period=period, qboFlg=False, ymin=-100, ymax=120, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrLCAnom2 = np.asarray(resvmrLCAnom2)

    
    #---------------------------------------------------
    # Stratospheric Weighted VMR ==> Retrieval
    #---------------------------------------------------
    print '\nPlot: Stratospheric Weighted VMR:\n' 
    resvmrSC = dc.AnalTS(hdf.npanels3-1, hdf.dates2, hdf.WvmrStrat2, hdf.pltID2, hdf.Lat2, hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Stratospheric wVMR', unitsStr=sclfctName, period=period, ymin=100 ,ymax=450, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrSC = np.asarray(resvmrSC)

    print '\nPlot: Stratospheric Weighted VMR Anomalies:\n' 
    resvmrSCAnom = dc.AnalTSAnom(hdf.npanels3-1, hdf.dates2, hdf.WvmrStrat2, hdf.pltID2, hdf.Lat2, hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Stratospheric wVMR Anomalies', unitsStr=sclfctName, period=period, ymin=-100 ,ymax=100, qboFlg=False, xmin=hdf.xmin, xmax=hdf.xmax)
    resvmrSCAnom = np.asarray(resvmrSCAnom)


    #print '\nPlot: Boundary Layer Column Anomalies:\n' 
    #resLCAnom = dc.AnalTSAnom(hdf.npanels3-1, hdf.dates2, hdf.WvmrTrop12, hdf.pltID2, hdf.Lat2, hdf.locID2, fits = fitFlg, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=hdf.pdfsav, ytypeStr='Low Tropospheric wVMR Anomalies', unitsStr=sclfctName,period=period, xmin=hdf.xmin, xmax=hdf.xmax)
    #resLCAnom = np.asarray(resLCAnom)

    #-------------------------
    # Monthly Prf
    #-------------------------
    if mntPltFlg: hdf.pltMnthlyVMR()

    #---------------------------------------------------
    # Bar plot (ROC)
    #---------------------------------------------------
    if fitFlg: dc.hbarplt3(resvmrLCAnom, resvmrLCAnom2, resvmrSCAnom, hdf.locID2, b1_label='Low Tropospheric', b2_label='Free Tropospheric', b3_label='Stratospheric', subtitle='ROC - Anomalies wVMR', saveFlg=saveFlg, pdfsav=hdf.pdfsav)
    
    #---------------------------------------------------
    # Hemispheric 
    #---------------------------------------------------
    if fitFlg: dc.latPlt(resvmrLC, resvmrLC2, resvmrSC, hdf.Lat2, saveFlg=saveFlg, pdfsav=hdf.pdfsav)

    #---------------------------------------------------
    # Save/Show figures
    #---------------------------------------------------
    if saveFlg:
        hdf.closeFig()   
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()


 
if __name__ == "__main__":
    main()
