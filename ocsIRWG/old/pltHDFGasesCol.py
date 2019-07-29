#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltHDFGasesCol.py
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
import HDFClassRead as dc
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


    
                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #-----------------------------------------------------------------------------------------
    #                             Initializations
    #-----------------------------------------------------------------------------------------

    #-------------------------------------
    #Global Data Directory
    #-------------------------------------
    GDataDir     = '/data1/projects/ocs/'
    #-------------------------------------
    #three letter ID
    #-------------------------------------
    #locs         = ['ahs', 'ldr', 'bld', 'kir']
    #locs         = ['kir', 'iza', 'bld', 'stp', 'jfj', 'wlo']
    #locs         = ['bld', 'wlo', 'jfj', 'mlo', 'tab', 'tor', 'eur', 'stp', 'ldr', 'rkb', 'tsk', 'ahs', 'mai', 'std', 'zgp', 'kir', 'iza', 'par', 'bre', 'nya', 'pmb', 'alt']
    locs         = ['bld', 'wlg', 'jfj', 'mlo', 'tab', 'tor', 'eur', 'spu', 'ldr', 'rkb', 'tsk', 'ahs', 'mai', 'std', 'zgp', 'kir', 'iza', 'par', 'bre', 'nya', 'pmb', 'alz']

    #-------------------------------------
    #ID in the HDF Files
    #-------------------------------------
    #locID        = ['arrival.heights', 'laud_120hr', 'boulder', 'kiruna']#, 'kiruna', 'maido' , 'stdenis']
    #locID        = ['kiruna', 'izana', 'boulder', 'st.petersburg', 'jungfraujoch', 'wollongong']
    locID        = ['boulder', 'wollongong', 'jungfraujoch', 'mauna.loa.h', 'thule', '_toronto_', '_eureka_', 'st.petersburg', '_laud_', 'rikubetsu', 
                   'tsukuba', 'ahts',  'reunion.maido' , 'stdenis', 'zugspitze', 'kiruna', 'izana', 'paris', 'bremen', 'ny.alesund', 'paramaribo', 'altzomoni']                       
    #-------------------------------------
    #Names in Plots
    #-------------------------------------   
    #pltID        = [ 'AHTS', 'Lauder', 'Boulder', 'Kiruna']
    #pltID        = ['Kiruna', 'Izana', 'Boulder', 'St Petersburg', 'Jungfraujoch', 'Wollongong']
    pltID        = ['Boulder', 'Wollongong', 'Jungfraujoch', 'Mauna Loa', 'Thule', 'Toronto', 'Eureka', 'St Petersburg', 'Lauder', 'Rikubetsu', 
                   'Tsukuba', 'AHTS', 'Maido', 'St Denis', 'Zugspitze', 'Kiruna', 'Izana', 'Paris', 'Bremen', 'Ny Alesund', 'Paramaribo', 'Altzomoni'] 
     
    #-------------------------------------
    #Inputs
    #-------------------------------------
    gasName        = 'ocs'

    AvgType        = 'Monthly'   #'Monthly'  'Daily'
    smthFlg        = False
    period         = 1.0
    fitFlg         = True

    #------
    # Flags
    #------
    saveFlg       = False                  # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
    errorFlg      = False                  # Flag to process error data
    fltrFlg       = True                   # Flag to filter the data

    dateFlg       = True                  # Flag to filter based on min and max dates
    tcFlg         = True                   # Flag to filter total column amount < 0
    tcMMFlg       = True                   # Flag to filter based on min and max total column amount
    pcFlg         = True                     # Flag to filter profiles with negative partial columns
    szaFlg        = True                   # Flag to filter based on min and max SZA    

    minSZA        = 0.0                    # Min SZA for filtering
    maxSZA        = 90.0                   # Max SZA for filtering
    maxTC         = 1.0e25                 # Max Total column amount for filtering
    minTC         = 0.0                    # Min Total column amount for filtering

    iyear         = 1984   
    imonth        = 1
    iday          = 1
    fyear         = 2020
    fmonth        = 12
    fday          = 31
    
    sclfct        = 1.0E12                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName    = 'ppt'                 # Name of scale factor for labeling plots
    TCsclfct      = 1.0e15
    TCsclfctName  = 'x10$^{15}$'


    pltPcol       = False                  #plot the time series in partial columns
    pltWvmr       = True                   #plot the time series in weighted VMR
    
    Adth          = 16.0                   #Altitude in km of tropopause in case NCEP or DTH is not available
    offH          = 5.0                    #Additional altitude above the tropopause height

    alt_maido = [97.25,        89.84999847,  81.09999847,  73.38500214,  66.58499908,
  60.59000015,  55.31499863,  50.70000076,  46.68000031,  43.18999863,
  40.16999817,  37.55500031,  35.28499985,  33.29999924,  31.53000069,
  29.91500092,  28.39999962,  26.94000053,  25.51499939,  24.125,       22.77000046,
  21.45000076,  20.17000008,  18.92499924,  17.71500015,  16.54000092,
  15.39999962,  14.30000019,  13.23499966,  12.20499992,  11.21000004,  10.25,
   9.32499981,   8.43500042,   7.57499981,   6.74499989,   5.94999981,
   5.19000006,   4.46000004,   3.7650001,    3.0999999 ,   2.4649999]

    #-------------------------------------
    # Flag for Plots
    #-------------------------------------
    prfFlg   = False
    mapFlg   = False
    aksFlg   = False





                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    #-------------------------------------
    #Name of PDF with Figures
    #-------------------------------------
    if (pltPcol) and not (pltWvmr):      pltFile =  '/data1/projects/ocs/figures/HDF_'+gasName.upper()+'_tpp_'+str('2std')+'_pCol'+'_'+AvgType+'.pdf'
    elif (pltWvmr) and not (pltPcol):    pltFile  =  '/data1/projects/ocs/figures/HDF_'+gasName.upper()+'_tpp_'+str('std')+'_wVMR'+'_'+AvgType+'.pdf'
    elif (pltPcol) & (pltWvmr):          pltFile  =  '/data1/projects/ocs/figures/HDF_'+gasName.upper()+'_tpp_'+str('std')+'_pCol_wVMR'+'_'+AvgType+'.pdf'
    else: pltFile = 'test.pdf'

    #-------------------------------------
    # Check file and directories
    #-------------------------------------
    dataDir    = [GDataDir+l+'/'  for l in locs]

    for d in dataDir:  ckDir(d,exit=True)
    ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------------------
    # Create instance of output data class   
    #-------------------------------------
    statDataCl = OrderedDict()
    
    Group = zip(dataDir,locID, pltID, locs)
    Group.sort(key=lambda Group: Group[2])
    pltID.sort()
    locs = [l for dd, id, pl, l in Group]


    for dd, id, pl, l in Group:
        #-------------------------------------
        # Some HDF files are in specific folder: change here accordingly
        #-------------------------------------
        if pl == 'Wollongong':      dd = dd + 'ocs_hippov2/'
        elif pl == 'Jungfraujoch' : dd = dd + 'hdf_OCS/'
        elif pl == 'Toronto' :      dd = dd + 'OCS/'
        #elif pl == 'Eureka' :       dd = dd + 'OCS/'
        elif pl == 'Rikubetsu':     dd = dd + 'HDF_Fil4/'
        elif pl == 'Tsukuba' :      dd = dd + 'HDFfiles/'
        elif pl == 'Zugspitze':     dd = dd + 'OCS_Zugspitze/'
        elif pl == 'Kiruna':        dd = dd + 'OCS_Kiruna/'
        elif pl == 'Izana':         dd = dd + 'OCS_Izana/'
        elif pl == 'St Petersburg': dd = dd + 'HDF_OCS_SPb_O3_atm16/'
        elif pl == 'Paris':         dd = dd + '2019_Paris/'
        else: dd = dd

        statDataCl[pl] = dc.ReadHDFData(dd, id, gasName)

    #-------------------------------------
    # Variables from HDF files 
    #-------------------------------------
    datesJD2K    = OrderedDict()
    totClmn      = OrderedDict()   #retrieved total column (molec/cm2)
    atotClmn     = OrderedDict()   #apriori total column (molec/cm2)
    sza          = OrderedDict()   #Solar Zenith Angle

    #-------------------------------------
    # Variables calculated 
    #-------------------------------------
    #alt_orig     = OrderedDict()
    dates        = OrderedDict()

    totWvmr      = OrderedDict()    #Weightet VMR A priori
    atotWvmr     = OrderedDict()


    Lat          = []
    Lon          = []

    if errorFlg:
        tot_rnd       = OrderedDict()
        tot_sys       = OrderedDict()
        tot_std       = OrderedDict()
        vmr_rnd_err   = OrderedDict()
        vmr_sys_err   = OrderedDict()
        vmr_tot_err   = OrderedDict()


    for ii, idhdf in enumerate(pltID):

        print idhdf

        datesJD2K[idhdf]    = statDataCl[idhdf].HDF[statDataCl[idhdf].getDatetimeName()]
        dates[idhdf]        = dc.jdf_2_datetime(datesJD2K[idhdf])
        sza[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAngleSolarZenithAstronomicalName()]
        

        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()+'VAR_SI_CONVERSION']
        totClmn[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()]*float(conv[0][1]) * (6.02e23) /100./100. / TCsclfct
        atotClmn[idhdf]     = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarAprioriName()]*float(conv[0][1]) * (6.02e23) /100./100. / TCsclfct

        #----------------------------------------
        #READ LAT/LON/HEIGHT OF INSTRUMENT
        #----------------------------------------
        Lat_i           = statDataCl[idhdf].HDF[statDataCl[idhdf].getLatitudeInstrumentName()]
        Lon_i           = statDataCl[idhdf].HDF[statDataCl[idhdf].getLongitudeInstrumentName()]
        alt_instru      = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeInstrumentName()]

        Lat.append(float(Lat_i[0]))
        Lon.append(float(Lon_i[0]))

        print '\n'
        print idhdf
        print 'Latitude          = {0:.2f}'.format(Lat_i[0])
        print 'Longitude         = {0:.2f}'.format(Lon_i[0])
        print 'Altitude of Instr = {0:.2f}'.format(alt_instru[0])


        #----------------------------------------
        #OBTAIN ERROR VARIABLES
        #---------------------------------------- 
        if errorFlg:                               
            tot_rnd[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintyRandomName()]
            tot_sys[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintySystematicName()]
            tot_std[idhdf]      = np.sqrt(tot_rnd[idhdf]**2 + tot_sys[idhdf]**2)

            npnts               = np.shape(tot_std[idhdf])[0]
            

            for i in range(npnts):
                conv    = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()+'VAR_SI_CONVERSION']  
                cov_rnd = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()]
                cov_sys = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintySystematicName()]

              
               

        #----------------------------------------
        # FILTER DATA
        #----------------------------------------
        if fltrFlg: statDataCl[idhdf].fltrData(statDataCl[idhdf].PrimaryGas,iyear=iyear, imonth=imonth, iday=iday, fyear=fyear, fmonth=fmonth, fday=fday, minsza=minSZA,
                                               mxsza=maxSZA,minTC=minTC,maxTC=maxTC, tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,tcMMFlg=tcMMFlg, dateFlg=dateFlg)
        else:    statDataCl[idhdf].inds = np.array([]) 
        
        try:
            dates[idhdf]    = np.delete(dates[idhdf], statDataCl[idhdf].inds)
            sza[idhdf]      = np.delete(sza[idhdf], statDataCl[idhdf].inds)
            totClmn[idhdf]  = np.delete(totClmn[idhdf], statDataCl[idhdf].inds)
            atotClmn[idhdf] = np.delete(atotClmn[idhdf], statDataCl[idhdf].inds)

        except Exception as errmsg:
            print '\nError: ', errmsg

        if errorFlg:
            try:
            
                tot_rnd[idhdf]      = np.delete(tot_rnd[idhdf],statDataCl[idhdf].inds)
                tot_sys[idhdf]      = np.delete(tot_sys[idhdf],statDataCl[idhdf].inds)
                tot_std[idhdf]      = np.delete(tot_std[idhdf],statDataCl[idhdf].inds)
            except Exception as errmsg:
                print '\nError: ', errmsg

    

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                           PLOTS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------
    print '\nPrinting Plots.......\n'
    if saveFlg: pdfsav = PdfPages(pltFile)
    else: pdfsav = ''
    clr = mf.clrplt()

    #---------------------------------------------------
    # Order data based on +Lat to -Lat
    #---------------------------------------------------
    pltID  = [y for (x, y) in sorted(zip(Lat,pltID), reverse=True)]
    locsID = [y for (x, y) in sorted(zip(Lat,locs), reverse=True)]
    Lon    = [y for (x, y) in sorted(zip(Lat,Lon), reverse=True)]
    
    Lat    = sorted(Lat, reverse=True)
    pltID  = np.asarray(pltID)
    locsID = np.asarray(locsID)



    #---------------------------------------------------
    # Defining variable for plots
    #---------------------------------------------------
    npanels   = len(pltID)

    xmin      = dt.date(iyear, imonth, iday)
    xmax      = dt.date(fyear, fmonth, fday)
    clmap     = 'jet'
    
    npanels2  = int(math.ceil(npanels/4.0))

    print '\nPlot: Averaged Total Columns:\n' 
    resTC = AnalTS2(npanels2, dates, totClmn, pltID, Lat, locsID, fits=True, AvgType=AvgType, pltFig=fitFlg, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Total Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=4.5 ,ymax=12)
    resTC = np.asarray(resTC)
    
    user_input = raw_input('Press any key to exit >>> ')
    sys.exit()

    

#------------------------------------------------------------------------------------------------------------------------------------------------------
#                                                                END
#------------------------------------------------------------------------------------------------------------------------------------------------------
                                    
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


def squiggle_xy(a, b, c, d, i=np.arange(0.0, 2*np.pi, 0.05)):
    return np.sin(i*a)*np.cos(i*b), np.sin(i*c)*np.cos(i*d)

def AnalTSAnom(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1, qboFlg=False):
    
    #--------------------
    #Slope and Amplitude of time series
    #--------------------
    slope       = []   #slope
    slope_e     = []   #slope error

    slope_p1    = []   #slope 1995 - 2002
    slope_p1_e  = []   #slope error 

    slope_p2    = []   #slope 2002 - 2008
    slope_p2_e  = []   #slope error

    slope_p3    = []   #slope  2008 - 2016
    slope_p3_e  = []   #slope error

    amp         = []   #Amplitude 

    avgD        = []
    stdD        = []

    avgD_p1     = []
    stdD_p1     = []

    avgD_p2     = []
    stdD_p2     = []

    avgD_p3     = []
    stdD_p3     = []

    xmin      = dt.date(1985, 1, 1)
    xmax      = dt.date(2019, 12, 31)


    if pltFig:    
        
        fig  = plt.figure(figsize=(18,13))
        fig2 = plt.figure(figsize=(18,13))  
        
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.05, hspace=0.075)

    fi = 0

    for i, idhdf in enumerate(pltID):

        ax   = plt.Subplot(fig, outer_grid[i])
        ax2  = plt.Subplot(fig2, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_facecolor('lightcyan')
            ax2.set_facecolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_facecolor('lightgreen')
            ax2.set_facecolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_facecolor('mistyrose')
            ax2.set_facecolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_facecolor('cornsilk')
            ax2.set_facecolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_facecolor('lightgrey')
            ax2.set_facecolor('lightgrey')
        
        else:
            ax.set_facecolor('lightgrey')
            ax2.set_facecolor('lightgrey')
        
        fig.add_subplot(ax)

        fig2.add_subplot(ax2)

        #----------------------------
        if AvgType == 'Daily':
            Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      = Avg['dailyAvg']

            
        elif AvgType == 'Monthly':
            Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      =  Avg['mnthlyAvg']

        elif AvgType == 'none':
            Dates        = xDates[idhdf]
            dateYearFrac = mf.toYearFraction(Dates)
            AvgData      = yData[idhdf]

        else:
            print 'Error: Define average type: Daily, Monthly, or none'
            exit()

        indsnan = [ni for ni, n in enumerate(AvgData) if np.isnan(n)]

        if len(indsnan) == len(AvgData):

            slope.append(float('nan'))
            slope_e.append(float('nan'))

            slope_p1.append(float('nan'))
            slope_p1_e.append(float('nan'))

            slope_p2.append(float('nan'))
            slope_p2_e.append(float('nan'))

            slope_p3.append(float('nan'))
            slope_p3_e.append(float('nan'))

            amp.append(float('nan'))

            avgD.append(float('nan'))
            stdD.append(float('nan'))

            avgD_p1.append(float('nan'))
            stdD_p1.append(float('nan'))

            avgD_p2.append(float('nan'))
            stdD_p2.append(float('nan'))

            avgD_p3.append(float('nan'))
            stdD_p3.append(float('nan'))


        else:
            #continue

    
            #---------------------------------------------------
            #To make a continuous fit
            #---------------------------------------------------
            numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
            dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
            dates2  = np.asarray(dates2)
            dateYearFrac2 = mf.toYearFraction(dates2)
            #---------------------------------------------------

            yyyy = [sngdate.year for sngdate in Dates]
            yyyy = np.asarray(yyyy)

            months = [sngdate.month for sngdate in Dates]
            months = np.asarray(months)


            #--------------------
            # Apply savitzky golay filter (Comment out if not wanted)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            #--------------------
            # Anomalies
            #--------------------
            if AvgType == 'Monthly':

                weights      = np.ones_like(dateYearFrac)

                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                AnnualCy        = f_fourier(dateYearFrac)
                Anomaly         = AvgData - AnnualCy


                # Anomaly  = np.zeros(len(Dates))

                # month    = np.array([d.month for d in Dates])
                # mnthSort = list(set(month))
                
                # mnthMean = np.zeros(len(mnthSort))
                # mnthSTD  = np.zeros(len(mnthSort))
                
                # for j,m in enumerate(mnthSort):
                #     inds        = np.where(month == m)[0]
                #     mnthMean[j] = np.mean(AvgData[inds])
                #     mnthSTD[j]  = np.std(AvgData[inds])

                #     Anomaly[inds]   = AvgData[inds]  - mnthMean[j]

            if AvgType == 'Daily':
                weights      = np.ones_like(dateYearFrac)

                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                AnnualCy        = f_fourier(dateYearFrac)
                Anomaly         = AvgData - AnnualCy

            #--------------------
            # Removing Additional Modulation (QBO, solar cycle, etc) --> TO BE IMPROVED
            #--------------------
            if qboFlg:
                weights      = np.ones_like(dateYearFrac)

                #res          = mf.fit_driftfourier_qbo(dateYearFrac, Anomaly, weights, 2, half_period=2.3)
                res          = mf.fit_qbo(dateYearFrac, Anomaly+300., weights, 2, half_period=2.33)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                Mdl         = f_fourier(dateYearFrac)
                Anomaly2    = Anomaly - Mdl + 300.

            else:

                Anomaly2    = Anomaly

            #--------------------
            # Plot Initial Anomalies (removing seasonal cycle)
            #--------------------
            if pltFig:
                
                ax.plot(Dates, Anomaly,'k.',markersize=0)
                ax.scatter(Dates, Anomaly, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)
                #print len(Dates)
                #print len(f_fourier(dateYeardateYearFracFrac))
                if qboFlg: ax.plot(dates2,f_fourier(dateYearFrac2)-300.,label='Fitted Anual Trend',linewidth=2.0)

                ax.axhline(y=0, linestyle='--', linewidth=1.5, color='gray', alpha=0.5)


                #ax2 = plt.Subplot(fig2, outer_grid[fi])
                fi += 1
                ax2.plot(Dates, Anomaly2,'k.',markersize=0)
                ax2.scatter(Dates, Anomaly2, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                ax2.axhline(y=0, linestyle='--', linewidth=1.5, color='gray', alpha=0.5)

                #--------------------
                #
                #--------------------

                if qboFlg:
                    qboDir        = '/data1/projects/ocs/'                 # Name of station location
                    qbofile       = 'QBO.dat' 

                    cols, indexToName = mf.getColumns(qboDir + qbofile, headerrow=8, delim=' ', header=True)

                    YYMM  = cols[indexToName[1]]
                    qbo50 = cols[indexToName[3]]
                    qbo30 = cols[indexToName[5]]

                    qboDate = []

                    for q in YYMM:
                        YY = q[0:2]
                        MM = q[2:4]

                        if int(YY) <= 20 :
                            YYYY = '20'+YY
                        else:
                            YYYY = '19'+YY

                        qboDate.append(dt.date(int(YYYY), int(MM), 15))

                    qboDate    = np.asarray(qboDate)
                    qbo50   = np.asarray(qbo50, dtype=float)
                    qbo30   = np.asarray(qbo30, dtype=float)

                    axx = ax.twinx()
                    axx.plot(qboDate, qbo30, color= 'gray')
                #fig.add_subplot(axx)

            #--------------------

            if fits:
                
                #---------------------------------------------------
                #LINEAR FIT AND DRIFT FOR ALL YEARS
                #---------------------------------------------------
                #if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Paris'):
                yoi    = [1996, 2018]

                if idhdf == 'AHTS':
                    

                    indx1  = np.where( (yyyy >= yoi[0]) &  (yyyy <= yoi[1]) & (months >= 1) & (months <= 12) )[0]

                else:
                    indx1  = np.where( (yyyy >= yoi[0]) &  (yyyy <= yoi[1]))[0]

                weights      = np.ones_like(dateYearFrac)

                res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A, df_drift = res[3:9]

                res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yoi[0], yoi[-1])

                if (idhdf == 'Zugspitze') or (idhdf == 'Jungfraujoch') or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Thule'):

                    #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                    

                    slope.append(res[1]/np.mean(AvgData[indx1])*100.0)
                    slope_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                    #print res[1]/np.mean(AvgData)*100.0

                    Amp   = np.sum(res[2]**2)
                    Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                    amp.append(Amp)

                    avgD.append(np.mean(AvgData[indx1]))
                    stdD.append(np.std(AvgData[indx1]))

                    ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color='k')

                else:

                    slope.append(float('nan'))
                    slope_e.append(float('nan'))

                    avgD.append(float('nan'))
                    stdD.append(float('nan'))

                

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, Anomaly, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, Anomaly, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                #---------------------------------------------------
                # POLY FIT FOR STATIONS WITH LONG TIME SERIES (TO KNOW INFLECTION POINTS)
                #---------------------------------------------------
                if (idhdf != 'Eureka') & (idhdf != 'St Petersburg') & (idhdf != 'Boulder') & (idhdf != 'Maido') &  (idhdf != 'St Denis') & (idhdf != 'Bremen') & (idhdf != 'Paris') & (idhdf != 'Altzomoni'):

                    res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly2, weights, 2, half_period=period)
                    f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]

            #---------------------------------------------------
            # Finding local min & Max (decrease or increase rate of change)
            #---------------------------------------------------
                    if fits: 
                        #roc  = df_drift(dateYearFrac2)/np.mean(AvgData)*100.0
                        #roc  = np.asarray(roc)
                        #zero_crossings = np.where(np.diff(np.sign(roc)))[0]

                        a = diff(sign(diff(f_driftfourier(dateYearFrac2)))).nonzero()[0] + 1 # local min+max
                        b = (diff(sign(diff(f_driftfourier(dateYearFrac2)))) > 0).nonzero()[0] + 1 # local min
                        c = (diff(sign(diff(f_driftfourier(dateYearFrac2)))) < 0).nonzero()[0] + 1 # local max
                        
                        #ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                        #ax.plot(dates2[a], f_driftfourier(dateYearFrac2)[a], "o", color='green')
                        ##ax.plot(dates2[c], f_driftfourier(dateYearFrac2)[c],  "o", color='red')

                
                        print idhdf
                        print 'Number of crossings = {}'.format(len(dates2[a]))
                        print 'Dates of crossings (minimums)  = {}'.format(dates2[b])
                        print 'Dates of crossings (maximums)  = {}'.format(dates2[c])

                # if pltFig:
                #     #ax.plot(dates2,f_fourier(dateYearFrac2),label='Fitted Anual Trend', linewidth=2.0)
                #     ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #     ax2.plot(dates2,f_drift(dateYearFrac2),label='Fitted Anual Trend',linewidth=2.0)

                #     if idhdf.lower() == 'jungfraujoch':

                #         fig, ax = plt.subplots(2, figsize=(10,6), sharex=True)
                #         ax[0].plot(dates2, f_driftfourier(dateYearFrac2), linewidth=2.0)
                #         ax[0].plot(dates2[b], f_driftfourier(dateYearFrac2)[b], "o", color='green')
                #         ax[0].plot(dates2[c], f_driftfourier(dateYearFrac2)[c],  "o", color='red')

                #         ax[1].plot(dates2, df_drift(dateYearFrac2), linewidth=2.0)

                #         plt.show(block= False)

                #         user_input = raw_input('Press any key to exit >>> ')
                #         sys.exit()   
                
            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    
                # ax = plt.Subplot(fig, outer_grid[i])
                # ax.plot(Dates, AvgData,'k.',markersize=0)
                # ax.scatter(Dates,AvgData, s=10, facecolor='lightgray', edgecolor='k',alpha=0.85)

                # if yData2:
                #     ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                # if yData3:
                #     ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                # if yData4:
                #     dtpStd      = np.nanstd(yData4[idhdf])
                #     dtpMean     = np.nanmean(yData4[idhdf])

                #     ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                #     ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  \

                

                


                #---------------------------------------------------


                # if (idhdf == 'Paris') or (idhdf == 'Altzomoni'):


                #     slope_p1.append(float('nan'))
                #     slope_p1_e.append(float('nan'))

                #     slope_p2.append(float('nan'))
                #     slope_p2_e.append(float('nan'))

                #     slope_p3.append(float('nan'))
                #     slope_p3_e.append(float('nan'))

                #     avgD_p1.append(float('nan'))
                #     stdD_p1.append(float('nan'))

                #     avgD_p2.append(float('nan'))
                #     stdD_p2.append(float('nan'))

                #     avgD_p3.append(float('nan'))
                #     stdD_p3.append(float('nan'))



                if (idhdf == 'Eureka') or (idhdf == 'Kiruna') or (idhdf == 'St Petersburg') or (idhdf == 'Bremen') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Tsukuba')  or (idhdf == 'Izana') or (idhdf == 'Paramaribo') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                    yoi  = [[2009, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color='red')


                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    clr_yoi = ['green', 'red']

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
        
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                      
                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])
                  
                        ##ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Rikubetsu') :

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))

                    clr_yoi = ['blue', 'green']
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])



                elif ( (idhdf == 'Jungfraujoch')  or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze') or (idhdf == 'Mauna Loa') ):

                    yoi  = [[1996, 2002], [2002, 2008], [2009, 2018]]

                    clr_yoi = ['blue', 'green', 'red']
                
                    for ii, y in enumerate(yoi):

                        
                        if idhdf == 'AHTS':
                        
                            indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]) & (months >= 1) & (months <= 12) )[0]

                            #ax.plot(Dates[indx1], Anomaly[indx1],'k.',markersize=0)
                            #ax.scatter(Dates[indx1], Anomaly[indx1], s=40, facecolor='red', edgecolor='k',alpha=0.85)

                        else:
    
                            indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])

                        
                    ##ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)


        if pltFig:
            ax.set_xlim(xmin, xmax)

            ax.grid(True, color='gray', alpha=0.25)
            ax.tick_params(which='both',labelsize=10)
            #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
            ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
            #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
            #start, end = ax1[i].get_xlim()
            #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
            #ax1[i].set_ylim(bottom=0)

            yearsLc1      = YearLocator(2)
            yearsLc2      = YearLocator(1)
            months        = MonthLocator()
            DateFmt       = DateFormatter('%Y')

            #plt.xticks(rotation=45)
            ax.xaxis.set_major_locator(yearsLc1)
            ax.xaxis.set_minor_locator(yearsLc2)
            #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax.xaxis.set_major_formatter(DateFmt)
            #ax.set_xlabel('Year')
            #ax1.xaxis.set_tick_params(which='major', pad=15)  
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')
            ax.tick_params(which='both')

            if (ymin and ymax):
                ax.set_ylim(ymin, ymax)
          
            #ax.set_xticks([])
            #ax.set_yticks([])
            #fig.add_subplot(ax)


            ax2.set_xlim(xmin, xmax)
            #ax2.set_ylim(-7, 7)
            #ax2.axhline(y=0, linestyle='--', color='k')
            ax2.grid(True, color='gray', alpha=0.25)
            ax2.tick_params(which='both',labelsize=10)
            ax2.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')

            ax2.xaxis.set_major_locator(yearsLc1)
            ax2.xaxis.set_minor_locator(yearsLc2)
            #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax2.xaxis.set_major_formatter(DateFmt)
            #ax.set_xlabel('Year')
            #ax1.xaxis.set_tick_params(which='major', pad=15)  
            ax2.xaxis.set_tick_params(which='minor',labelbottom='off')
            ax2.tick_params(which='both')

        

            if (ymin and ymax):
                ax2.set_ylim(ymin, ymax)
    
            #fig2.add_subplot(ax2)


            #---------------------------------------------------

    if pltFig:
    
        # all_axes = fig.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')
        
        # fig.autofmt_xdate()
        # fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        # fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        # fig.suptitle(ytypeStr, fontsize=16  )
        
        #-----
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
            
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)

        
        for i in range(-1, -5, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        #fig.text(0.02, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

        if saveFlg: pdfsav.savefig(fig,dpi=200)
        else: 
            plt.show(block=False)

        ytypeStr = ytypeStr.replace(' ', '')

        if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
        else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')

    #---------------------------------------------------
    #
    #---------------------------------------------------
    if pltFig:
    
        # all_axes = fig2.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')
        
        # fig2.autofmt_xdate()
        # fig2.text(0.03, 0.5, 'Rate of change [%/y]', fontsize=16, va='center', rotation='vertical')
        # fig2.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        # fig2.suptitle('Rate of change of VMR Anomalies', fontsize=16  )

        #-----
        all_axes = fig2.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
            
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)

        for i in range(-1, -5, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        #fig2.autofmt_xdate()

        #fig2.text(0.02, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #fig2.text(0.0075, 0.5, 'Rate of change [%/y]', fontsize=16, va='center', rotation='vertical')
        fig2.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig2.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

        
        if saveFlg: pdfsav.savefig(fig2,dpi=200)
        else: 
            plt.show(block=False)

        ytypeStr = ytypeStr.replace(' ', '')

        #if fits: fig2.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'ROC_wFit.pdf', bbox_inches='tight')
        #else: fig2.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'ROC.pdf', bbox_inches='tight')


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)


def AnalTS(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1):
    
    #--------------------
    #Slope and Amplitude of time series
    #--------------------
    slope       = []   #slope
    slope_e     = []   #slope error

    slope_p1    = []   #slope 1995 - 2002
    slope_p1_e  = []   #slope error 

    slope_p2    = []   #slope 2002 - 2008
    slope_p2_e  = []   #slope error

    slope_p3    = []   #slope  2008 - 2016
    slope_p3_e  = []   #slope error

    amp         = []   #Amplitude 

    avgD        = []
    stdD        = []

    avgD_p1     = []
    stdD_p1     = []

    avgD_p2     = []
    stdD_p2     = []

    avgD_p3     = []
    stdD_p3     = []

    xmin      = dt.date(1985, 1, 1)
    xmax      = dt.date(2019, 12, 31)

    if pltFig:    
        fig = plt.figure(figsize=(18,13))  

        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.05, hspace=0.075)
    #--------------
    # with open('/data1/projects/ocs/figures/fig/'+ytypeStr+'.txt','w') as fopen:
    #     fopen.write('#Site\tLatitude\tMean TPH\tstd TPH\tMax TPH\tMin TPH\tAmplitude\n')

    
    #     strformat = '{:<15}\t'+'\t'.join('{:>5.1f}' for i in range(6)) + '\n'
        
    fi = 0
    for i, idhdf in enumerate(pltID):

        ax   = plt.Subplot(fig, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_facecolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_facecolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_facecolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_facecolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_facecolor('lightgrey')
        
        else:
            ax.set_facecolor('lightgrey')
        
        fig.add_subplot(ax)

        
        #----------------------------
        if AvgType == 'Daily':
            Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      = Avg['dailyAvg']

            
        elif AvgType == 'Monthly':
            Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      =  Avg['mnthlyAvg']

        elif AvgType == 'none':
            Dates        = xDates[idhdf]
            dateYearFrac = mf.toYearFraction(Dates)
            AvgData      = yData[idhdf]

        else:
            print 'Error: Define average type: Daily, Monthly, or none'
            exit()

        indsnan = [ni for ni, n in enumerate(AvgData) if np.isnan(n)]



        if len(indsnan) == len(AvgData):

    
        #if (idhdf == 'Mauna Loa') or (idhdf == 'Altzomoni') or (idhdf == 'Jungfraujoch') :

            slope.append(float('nan'))
            slope_e.append(float('nan'))

            slope_p1.append(float('nan'))
            slope_p1_e.append(float('nan'))

            slope_p2.append(float('nan'))
            slope_p2_e.append(float('nan'))

            slope_p3.append(float('nan'))
            slope_p3_e.append(float('nan'))

            amp.append(float('nan'))

            avgD.append(float('nan'))
            stdD.append(float('nan'))

            avgD_p1.append(float('nan'))
            stdD_p1.append(float('nan'))

            avgD_p2.append(float('nan'))
            stdD_p2.append(float('nan'))

            avgD_p3.append(float('nan'))
            stdD_p3.append(float('nan'))

        else:
            #continue



            #--------------------
            #Apply savitzky golay filter (Comment out if not wated)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            if fits:

                weights      = np.ones_like(dateYearFrac)

                #---------------------------------------------------
                # To make a continuous fit
                #---------------------------------------------------
                numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
                dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
                dates2  = np.asarray(dates2)
                dateYearFrac2 = mf.toYearFraction(dates2)
                #---------------------------------------------------

                yyyy = [sngdate.year for sngdate in Dates]
                yyyy = np.asarray(yyyy)
                
                #---------------------------------------------------
                # Calculate the linear regression of all years and estimate Fourier seris
                #---------------------------------------------------
                
                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                slope.append(res[1]/np.mean(AvgData)*100.0)
                slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                Amp   = np.sum(res[2]**2)
                Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                amp.append(Amp)

                avgD.append(np.mean(AvgData))
                stdD.append(np.std(AvgData))


                #fopen.write(strformat.format(idhdf, float(Lat[i]), np.mean(AvgData), np.std(AvgData), np.max(AvgData), np.min(AvgData), float(Amp)))


                # if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Bremen') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                    #res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                    #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]


            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    

                #ax = plt.Subplot(fig, outer_grid[fi])
                fi +=1
                ax.plot(Dates, AvgData,'k.',markersize=0)
                ax.scatter(Dates,AvgData, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                if yData2:
                    ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                if yData3:
                    ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                if yData4:
                    dtpStd      = np.nanstd(yData4[idhdf])
                    dtpMean     = np.nanmean(yData4[idhdf])

                    ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                    ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  

                # Lat[i] = float(Lat[i])

                # if Lat[i] >= 50.:     
                #     ax.set_facecolor('lightcyan')
                    
                # elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                #     ax.set_facecolor('lightgreen')
                   
                # elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                #     ax.set_facecolor('mistyrose')
                    
                # elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                #     ax.set_facecolor('cornsilk')
                    
                # elif (Lat[i] < -50.):
                #     ax.set_facecolor('lightgrey')
                
                # else:
                #     ax.set_facecolor('lightgrey')
     
            if fits:

                if pltFig:
                    #ax.plot(dates2,f_fourier(dateYearFrac2),label='Intra-annual variability', linewidth=1.0)
                    ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #---------------------------------------------------

                if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Bremen') or  (idhdf == 'Paris') or (idhdf == 'Altzomoni') or (idhdf == 'Tsukuba') or (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paramaribo'):

                    yoi  = [[2010, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))



                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                        roc    =  np.mean(roc)

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])


                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                elif (idhdf == 'Rikubetsu'):

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                  
                        #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Jungfraujoch')  or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze'):

                    yoi  = [[1995, 2002], [2002, 2008], [2009, 2018]]
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                    #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

        if pltFig:
                ax.set_xlim(xmin, xmax)

                ax.grid(True, color='gray', alpha=0.25)
                ax.tick_params(which='both',labelsize=10)
                #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
                ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
                #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
                #start, end = ax1[i].get_xlim()
                #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
                #ax1[i].set_ylim(bottom=0)

                yearsLc1      = YearLocator(2)
                yearsLc2      = YearLocator(1)
                months        = MonthLocator()
                DateFmt       = DateFormatter('%Y')


                #plt.xticks(rotation=45)
                ax.xaxis.set_major_locator(yearsLc1)
                ax.xaxis.set_minor_locator(yearsLc2)
                #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)
                #ax.set_xlabel('Year')
                #ax1.xaxis.set_tick_params(which='major', pad=15)  
                ax.xaxis.set_tick_params(which='minor',labelbottom='off')
                ax.tick_params(which='both')

                if (ymin and ymax):
                    ax.set_ylim(ymin, ymax)
              
                #ax.set_xticks([])
                #ax.set_yticks([])
                #fig.add_subplot(ax)

    if pltFig:
    
        # all_axes = fig.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        
        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')

        
        # fig.autofmt_xdate()

        # fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        # plt.suptitle(ytypeStr, fontsize=16  )

        # fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

        #-----
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
                # if ax.is_last_row():
                #     ax.spines['bottom'].set_visible(True)
                #     plt.setp(ax.get_xticklabels(), visible=True)

                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)


        for i in range(-1, -5, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        #fig.autofmt_xdate()

        fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

        
        if saveFlg: pdfsav.savefig(fig,dpi=200)
        else: 
            plt.show(block=False)

        ytypeStr = ytypeStr.replace(' ', '')

        if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
        else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)

def AnalTS2(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1):
        
        #--------------------
        #Slope and Amplitude of time series
        #--------------------
        slope       = []   #slope
        slope_e     = []   #slope error

        slope_p1    = []   #slope 1995 - 2002
        slope_p1_e  = []   #slope error 

        slope_p2    = []   #slope 2002 - 2008
        slope_p2_e  = []   #slope error

        slope_p3    = []   #slope  2008 - 2016
        slope_p3_e  = []   #slope error

        amp         = []   #Amplitude 

        avgD        = []
        stdD        = []

        avgD_p1     = []
        stdD_p1     = []

        avgD_p2     = []
        stdD_p2     = []

        avgD_p3     = []
        stdD_p3     = []

        xmin      = dt.date(1985, 1, 1)
        xmax      = dt.date(2019, 12, 31)

        if pltFig:    
            #fig = plt.figure(figsize=(18,13))  
            fig = plt.figure(figsize=(18,13))

            outer_grid = gridspec.GridSpec(npanels, 4, wspace=0.075, hspace=0.075)
            #outer_grid = gridspec.GridSpec(npanels2, 4, wspace=0.11, hspace=0.085) 

        for i, idhdf in enumerate(pltID):
            
            #----------------------------
            if AvgType == 'Daily':
                Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
                Dates        = Avg['dates']
                dateYearFrac = mf.toYearFraction(Avg['dates'])
                AvgData      = Avg['dailyAvg']

                
            elif AvgType == 'Monthly':
                Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
                Dates        = Avg['dates']
                dateYearFrac = mf.toYearFraction(Avg['dates'])
                AvgData      =  Avg['mnthlyAvg']

            elif AvgType == 'none':
                Dates        = xDates[idhdf]
                dateYearFrac = mf.toYearFraction(Dates)
                AvgData      = yData[idhdf]

            else:
                print 'Error: Define average type: Daily, Monthly, or none'
                exit()

            #--------------------
            #Apply savitzky golay filter (Comment out if not wated)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            if fits:

                weights      = np.ones_like(dateYearFrac)

                #---------------------------------------------------
                # To make a continuous fit
                #---------------------------------------------------
                numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
                dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
                dates2  = np.asarray(dates2)
                dateYearFrac2 = mf.toYearFraction(dates2)
                #---------------------------------------------------

                yyyy = [sngdate.year for sngdate in Dates]
                yyyy = np.asarray(yyyy)
                
                #---------------------------------------------------
                # Calculate the linear regression of all years and estimate Fourier seris
                #---------------------------------------------------
                
                #res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                #f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 3, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                slope.append(res[1]/np.mean(AvgData)*100.0)
                slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                Amp   = np.sum(res[2]**2)
                Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                amp.append(Amp)

                avgD.append(np.mean(AvgData))
                stdD.append(np.std(AvgData))


                # if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Bremen') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                    #res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                    #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]


            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    
                ax = plt.Subplot(fig, outer_grid[i])
                ax.plot(Dates, AvgData,'k.',markersize=0)
                ax.scatter(Dates,AvgData, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                if yData2:
                    ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                if yData3:
                    ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                if yData4:
                    dtpStd      = np.nanstd(yData4[idhdf])
                    dtpMean     = np.nanmean(yData4[idhdf])

                    ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                    ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  

                Lat[i] = float(Lat[i])

                if Lat[i] >= 50.:     
                    ax.set_facecolor('lightcyan')
                    
                elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                    ax.set_facecolor('lightgreen')
                   
                elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                    ax.set_facecolor('mistyrose')
                    
                elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                    ax.set_facecolor('cornsilk')
                    
                elif (Lat[i] < -50.):
                    ax.set_facecolor('lightgrey')
                
                else:
                    ax.set_facecolor('lightgrey')
     
            if fits:

                if pltFig:
                    #ax.plot(dates2,f_fourier(dateYearFrac2),label='Intra-annual variability', linewidth=1.0)
                    ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #---------------------------------------------------

                if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Bremen') or  (idhdf == 'Paris') or (idhdf == 'Altzomoni') or (idhdf == 'Tsukuba') or (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paramaribo'):

                    yoi  = [[2010, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))



                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                        roc    =  np.mean(roc)

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])


                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                elif (idhdf == 'Rikubetsu'):

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                  
                        #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Jungfraujoch')  or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze'):

                    yoi  = [[1995, 2002], [2002, 2008], [2009, 2018]]
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                    #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

            if pltFig:
                ax.set_xlim(xmin, xmax)

                ax.grid(True, color='gray', alpha=0.25)
                ax.tick_params(which='both',labelsize=10)
                #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
                ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
                #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
                #start, end = ax1[i].get_xlim()
                #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
                #ax1[i].set_ylim(bottom=0)

                yearsLc1      = YearLocator(2)
                yearsLc2      = YearLocator(1)
                months        = MonthLocator()
                DateFmt       = DateFormatter('%Y')


                #plt.xticks(rotation=45)
                ax.xaxis.set_major_locator(yearsLc1)
                ax.xaxis.set_minor_locator(yearsLc2)
                #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)
                #ax.set_xlabel('Year')
                #ax1.xaxis.set_tick_params(which='major', pad=15)  
                ax.xaxis.set_tick_params(which='minor',labelbottom='off')
                ax.tick_params(which='both')

                if (ymin and ymax):
                    ax.set_ylim(ymin, ymax)
              
                #ax.set_xticks([])
                #ax.set_yticks([])
                fig.add_subplot(ax)

        if pltFig:
        
            #-----
            all_axes = fig.get_axes()
            #show only the outside spines
            for i, ax in enumerate(all_axes):

                for sp in ax.spines.values():
                    sp.set_visible(False)
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.setp(ax.get_yticklabels(), visible=False)
                    ax.spines['top'].set_visible(True)
                    ax.spines['left'].set_visible(True)
                    ax.spines['right'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                
                    if ax.is_first_col():
                        ax.spines['left'].set_visible(True)
                        plt.setp(ax.get_yticklabels(), visible=True)
                        ax.tick_params(labelsize = 14)

            #fig.autofmt_xdate()

            for i in range(-1, -5, -1):

                all_axes[i].spines['bottom'].set_visible(True)
                plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
                all_axes[i].set_zorder(1)

                all_axes[i].set_xlabel('Year', fontsize=14)
                all_axes[i].tick_params(labelsize = 14)

            #fig.autofmt_xdate()

            fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
            #plt.suptitle(ytypeStr, fontsize=16  )

            fig.subplots_adjust(left=0.05, bottom=0.075, right=0.99, top=0.99)

            
            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else: 
                plt.show(block=False)

            ytypeStr = ytypeStr.replace(' ', '')


            if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
            else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')


        return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
                avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)


def hbarplt3(b1, b2, b3, pltID, b1_label='', b2_label='', b3_label='', subtitle='', saveFlg=False, pdfsav=' '):

    pltID        = np.asarray(pltID)

    pltID        = [p.upper() for p in pltID]


    #-------------------------------------------------
    # Calculate partial columns and weighted VMR
    #-------------------------------------------------

    with open('/data1/projects/ocs/figures/fig/'+'ROC_lowTrop.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b1[0][pi], b1[1][pi], b1[2][pi], b1[3][pi], b1[4][pi], b1[5][pi], b1[6][pi], b1[7][pi]     ) )

    with open('/data1/projects/ocs/figures/fig/'+'ROC_freeTrop.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b2[0][pi], b2[1][pi], b2[2][pi], b2[3][pi], b2[4][pi], b2[5][pi], b2[6][pi], b2[7][pi]     ) )

    with open('/data1/projects/ocs/figures/fig/'+'ROC_Strat.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b3[0][pi], b3[1][pi], b3[2][pi], b3[3][pi], b3[4][pi], b3[5][pi], b3[6][pi], b3[7][pi]     ) )


    
    #---------------------------------------------------
    # Bar plot: Three different periods ==> Retrieval
    #---------------------------------------------------
    

    ind = np.arange(len(b1[0]))
    
    fig, (ax, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=(15, 8.5))
    #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    ax.barh(ind-0.27, b1[0], 0.27, xerr=b1[1], align='center', color = 'r', ecolor = 'k', label = b1_label)
    ax.barh(ind, b2[0], 0.27, xerr=b2[1], align='center', color = 'b', ecolor = 'k', label = b2_label)
    ax.barh(ind+0.27, b3[0], 0.27, xerr=b3[1], align='center', color = 'g', ecolor = 'k', label = b3_label)  #yerr=slope_TC*0
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlabel('Rate of change [%/y]')
    ax.set_yticks(ind)
    ax.set_yticklabels(np.transpose(pltID), rotation=0)
    ax.set_title('1996 - 2016', multialignment='center')
    #ax.set_xlabel('Site')
    ax.set_xlim(-2.0, 2.0)
    ax.axvline(0, color='black', lw=1)
    ax.legend(prop={'size':8}, loc = 1)
    #ax.invert_yaxis()

    ax2.barh(ind-0.27, b1[2], 0.27, xerr=b1[3], align='center', color = 'r', ecolor = 'k', label = b1_label)
    ax2.barh(ind, b2[2], 0.27, xerr=b2[3], align='center', color = 'b', ecolor = 'k', label = b2_label)
    ax2.barh(ind+0.27, b3[2], 0.27, xerr=b3[3], align='center', color = 'g', ecolor = 'k', label = b3_label)  
    ax2.xaxis.grid(True)
    ax2.yaxis.grid(True)
    ax2.set_xlabel('Rate of change [%/y]')
    ax2.set_title('1996 - 2002', multialignment='center')
    ax2.set_yticks(ind)
    ax2.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    #ax2.legend(prop={'size':8}, loc = 1)
    ax2.set_xlim(-2.0, 2.0)
    ax2.axvline(0, color='black', lw=1)
    #ax2.invert_yaxis()

    #ax2.invert_yaxis()
    #ax2.yticks([])

    ax3.barh(ind-0.27, b1[4], 0.27, xerr=b1[5], align='center', color = 'r', ecolor = 'k')
    ax3.barh(ind, b2[4], 0.27, xerr=b2[5], align='center', color = 'b', ecolor = 'k')
    ax3.barh(ind+0.27, b3[4], 0.27, xerr=b3[5], align='center', color = 'g', ecolor = 'k')  
    ax3.xaxis.grid(True)
    ax3.yaxis.grid(True)
    ax3.set_xlabel('Rate of change [%/y]')
    ax3.set_title('2002 - 2008', multialignment='center')
    ax3.set_yticks(ind)
    #ax3.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    ax3.set_xlim(-2.0, 2.0)
    ax3.axvline(0, color='black', lw=1)

    #ax3.invert_yaxis()

    ax4.barh(ind-0.27, b1[6], 0.27, xerr=b1[7], align='center', color = 'r', ecolor = 'k')
    ax4.barh(ind, b2[6], 0.27, xerr=b2[7], align='center', color = 'b', ecolor = 'k')
    ax4.barh(ind+0.27, b3[6], 0.27, xerr=b3[7], align='center', color = 'g', ecolor = 'k')  
    ax4.xaxis.grid(True)
    ax4.yaxis.grid(True)
    ax4.set_xlabel('Rate of change [%/y]')
    ax4.set_title('2009 - 2018', multialignment='center')
    ax4.set_yticks(ind)
    #ax4.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    ax4.set_xlim(-2.0, 2.0)
    ax4.axvline(0, color='black', lw=1)

    #ax4.invert_yaxis()

    plt.gca().invert_yaxis()
    #fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.97)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        fig.savefig('/data1/projects/ocs/figures/fig/'+subtitle+'.pdf', bbox_inches='tight')
    else:       plt.show(block=False)

    subtitle = subtitle.replace(' ', '')

    


 
if __name__ == "__main__":
    main()
