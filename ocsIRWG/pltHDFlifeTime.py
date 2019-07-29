#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltHDFlifeTime.py
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

    region       = 'poleS' 

    #-------------------------------------
    #Global Data Directory
    #-------------------------------------
    GDataDir     = '/data1/projects/'
    #-------------------------------------
    #three letter ID ; ID in the HDF Files
    #-------------------------------------
        #-------------------------------------
    #
    #-------------------------------------
    if region    == 'poleN':
        locs         = ['kir', 'tab', 'spu', 'bre', 'eur', 'nya']
        locID        = ['kiruna', 'thule', 'st.petersburg', 'bremen', '_eureka_', 'ny.alesund'] 
        pltID        = ['Kiruna', 'Thule', 'St Petersburg', 'Bremen', 'Eureka', 'Ny Alesund'] 
    
    elif region == 'middleN':
        locs         = ['zgp', 'rkb',  'iza']  #'tor',  'jfj'
        locID        = ['zugspitze',  'rikubetsu', 'izana']  #'_toronto_', 'jungfraujoch', 
        pltID        = ['Zugspitze',  'Rikubetsu', 'Izana'] #'Toronto',  'Jungfraujoch',
    
    elif region == 'tropics':
        locs         = ['mlo', 'alz']#, 'pmb']
        locID        = ['mauna.loa.h', 'altzomoni']#, 'paramaribo'] 
        pltID        = ['Mauna Loa', 'Altzomoni']#, 'Paramaribo']
    
    elif region == 'middleS':
        locs         = ['std', 'mai', 'wlg', 'ldr']
        locID        = ['stdenis', 'maido', 'wollongong', 'niwa001' ] 
        pltID        = ['St Denis', 'Maido', 'Wollongong', 'Lauder']

    elif region == 'poleS':
        locs         = ['ahs']
        locID        = ['arrival.heights'] 
        pltID        = ['AHTS']
    else:
        print 'An error ocurred: region is not defined'
        exit()
    #-------------------------------------
    #Inputs
    #-------------------------------------
    gasName1      = 'ocs'
    gasName2      = 'n2o'

    AvgType        = 'Monthly'   #'Monthly'  'Daily'
    smthFlg        = False
    period         = 1.0
    fitFlg         = False

    ColFlg         = False


    #------
    # Flags
    #------
    saveFlg       = True                  # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
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

    iyear         = 2009   
    imonth        = 1
    iday          = 1
    fyear         = 2016
    fmonth        = 12
    fday          = 31
    
    sclfct        = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName    = 'ppb'                 # Name of scale factor for labeling plots
    
    TCsclfct      = 1.0e16
    TCsclfctName  = 'x10$^{16}$'

    pColsFlg      = True                   #Calculate tropospheric and stratospheric columns?

    pltPcol       = False                  #plot the time series in partial columns
    pltWvmr       = True                   #plot the time series in weighted VMR
    
    Adth          = 16.0                   #Altitude in km of tropopause in case NCEP or DTH is not available
    offH          = 5.0                    #Additional altitude above the tropopause height



    #-------------------------------------
    # Flag for Plots
    #-------------------------------------


                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    #-------------------------------------
    #Name of PDF with Figures
    #-------------------------------------
    if ColFlg: pltFile = GDataDir+'/ocs/figures/LifeTime_OCS_Column_'+region+'.pdf'
    else: pltFile = GDataDir+'/ocs/figures/LifeTime_OCS_pCol_'+region+'.pdf'

    if saveFlg: pdfsav = PdfPages(pltFile)
    else: pdfsav = ''

    #-------------------------------------
    # Check file and directories
    #-------------------------------------
    dataDir1    = [GDataDir+gasName1+'/'+l+'/'  for l in locs]
    dataDir2    = [GDataDir+gasName2+'/'+l+'/'  for l in locs]

    for d in dataDir1:  ckDir(d,exit=True)
    for d in dataDir2:  ckDir(d,exit=True)
    ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------------------
    # Create instance of output data class   
    #-------------------------------------
    statDataCl  = OrderedDict()
    statDataCl2 = OrderedDict()

    Group1 = zip(dataDir1,locID, pltID, locs)
    Group1.sort(key=lambda Group1: Group1[2])

    Group2 = zip(dataDir2,locID, pltID, locs)
    Group2.sort(key=lambda Group2: Group2[2])

    locs = [l for dd, id, pl, l in Group1]
    pltID.sort()

    
    for dd, id, pl, l in Group1:

        #-------------------------------------
        # Some HDF files are in specific folder: change here accordingly
        #-------------------------------------
        if pl == 'Wollongong':      dd = dd + 'ocs_hippov2/'
        elif pl == 'Jungfraujoch' : dd = dd + 'OCS.39_1b3144b4fe4a58f29f1f_/'
        elif pl == 'Toronto' :      dd = dd + 'OCS/'
        elif pl == 'Eureka' :       dd = dd + 'OCS/'
        elif pl == 'Rikubetsu':     dd = dd + 'HDF_Fil4/'
        elif pl == 'Tsukuba' :      dd = dd + 'HDFfiles/'
        elif pl == 'Zugspitze':     dd = dd + 'OCS_Zugspitze/'
        elif pl == 'Kiruna':        dd = dd + 'OCS_Kiruna/'
        elif pl == 'Izana':         dd = dd + 'OCS_Izana/'
        elif pl == 'St Petersburg': dd = dd + 'HDF_OCS_SPb_O3_atm16/'
        elif pl == 'Paris':         dd = dd + '2019_Paris/'
        else: dd = dd

        statDataCl[pl]  = dc.ReadHDFData(dd, id, gasName1)
    
    for dd, id, pl, l in Group2:

        statDataCl2[pl] = dc.ReadHDFData(dd, id, gasName2)


    #-------------------------------------
    # Variables from HDF files 
    #-------------------------------------
    datesJD2K    = OrderedDict()
    rPrf         = OrderedDict();  rPrf_2         = OrderedDict() #retrieved Prf in mixing ratio
    aPrf         = OrderedDict()   #apriori Prf in mixing ratio
    rPrfMol      = OrderedDict();  rPrfMol_2      = OrderedDict()   #retrieved Prf partial Column (molec/cm2)
    aPrfMol      = OrderedDict()   #apriori Prf partial Column (molec/cm2)
    totClmn      = OrderedDict();  totClmn_2      = OrderedDict() #retrieved total column (molec/cm2)
    atotClmn     = OrderedDict()   #apriori total column (molec/cm2)
    avkVMR       = OrderedDict()   #Averaging kernel (VMR)
    avkTC        = OrderedDict()   #Averaging kernel total column
    alt          = OrderedDict()   #Altitude 
    sza          = OrderedDict()   #Solar Zenith Angle
    TempPrf      = OrderedDict()   #Temperature Profile
    PresPrf      = OrderedDict()   #Pressure Profile

    #-------------------------------------
    # Variables calculated 
    #-------------------------------------
    #alt_orig     = OrderedDict()
    dates        = OrderedDict()
    dates_2       = OrderedDict()
    
    avkSCF       = OrderedDict()   #Averaging kernel (scale factor)
    dofs         = OrderedDict()   #degrees of freedom
    AirMPrf      = OrderedDict(); AirMPrf_2    = OrderedDict()   #Airmass
    rPrfMol      = OrderedDict()   #retrieved Prf in molec/cm2
    aPrfMol      = OrderedDict()   #apriori Prf in molec/cm2

    totWvmr      = OrderedDict()    #Weightet VMR A priori
    atotWvmr     = OrderedDict()

    alttpp       = OrderedDict()
    alttpp2      = OrderedDict()

    altbl1       = OrderedDict()
    altbl2       = OrderedDict()

    altft1       = OrderedDict()
    altft2       = OrderedDict()

    altst1       = OrderedDict()
    altst2       = OrderedDict()

    Lat          = []
    Lon          = []

    if errorFlg:
        tot_rnd       = OrderedDict()
        tot_sys       = OrderedDict()
        tot_std       = OrderedDict()
        vmr_rnd_err   = OrderedDict()
        vmr_sys_err   = OrderedDict()
        vmr_tot_err   = OrderedDict()

    if pColsFlg:
        dtp           = OrderedDict()
        datesdtp      = OrderedDict()
        
        PcolStrat     = OrderedDict()   #partial columns
        PcolTrop1     = OrderedDict()
        PcolTrop2     = OrderedDict()

        PcolStratapr  = OrderedDict()   #partial columns A priori
        PcolTropapr1  = OrderedDict()
        PcolTropapr2  = OrderedDict()

        WvmrStrat     = OrderedDict(); WvmrStrat_2     = OrderedDict()   #Weighted VMR
        WvmrTrop1     = OrderedDict()
        WvmrTrop1_2     = OrderedDict()
        WvmrTrop2     = OrderedDict(); WvmrTrop2_2     = OrderedDict()

        WvmrStratapr  = OrderedDict()    #Weighted VMR A priori
        WvmrTropapr1  = OrderedDict()
        WvmrTropapr2  = OrderedDict()

        rPcol         = OrderedDict(); rPcol_2         = OrderedDict() 
        aPcol         = OrderedDict()

        rPvmr         = OrderedDict(); rPvmr_2         = OrderedDict()
        aPvmr         = OrderedDict()


    for ii, idhdf in enumerate(pltID):

        print idhdf

        datesJD2K[idhdf]    = statDataCl[idhdf].HDF[statDataCl[idhdf].getDatetimeName()]
        dates[idhdf]        = dc.jdf_2_datetime(datesJD2K[idhdf])

        datesJD2K_2         = statDataCl2[idhdf].HDF[statDataCl2[idhdf].getDatetimeName()]
        dates_2[idhdf]      = dc.jdf_2_datetime(datesJD2K_2)
        

        alt[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeName()]
        sza[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAngleSolarZenithAstronomicalName()]
        
        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()+'VAR_SI_CONVERSION']            
        rPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()]*float(conv[0][1])*sclfct

        conv_2              = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl2[idhdf].getMixingRatioAbsorptionSolarName()+'VAR_SI_CONVERSION']    
        rPrf_2[idhdf]       = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()]*float(conv_2[0][1])*sclfct

        aPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAprioriName()]*float(conv[0][1])*sclfct

        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()+'VAR_SI_CONVERSION']
        rPrfMol[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()]*float(conv[0][1])*(6.02e23/100./100.)
        
        conv_2                = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl2[idhdf].getColumnPartialAbsorptionSolarName()+'VAR_SI_CONVERSION']
        rPrfMol_2[idhdf]    = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl2[idhdf].getColumnPartialAbsorptionSolarName()]*float(conv_2[0][1])*(6.02e23/100./100.)

        aPrfMol[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarAprioriName()]*float(conv[0][1])*(6.02e23/100./100.)

        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()+'VAR_SI_CONVERSION']
        totClmn[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()]*float(conv[0][1]) * (6.02e23) /100./100. / TCsclfct

        conv_2                = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl2[idhdf].getColumnAbsorptionSolarName()+'VAR_SI_CONVERSION']
        totClmn_2[idhdf]     = statDataCl2[idhdf].HDF[statDataCl2[idhdf].PrimaryGas.upper()+'.'+statDataCl2[idhdf].getColumnAbsorptionSolarName()]*float(conv_2[0][1]) * (6.02e23) /100./100. / TCsclfct

        atotClmn[idhdf]     = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarAprioriName()]*float(conv[0][1]) * (6.02e23) /100./100. / TCsclfct
        
        PresPrf[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].getPressureIndependentName()]
        TempPrf[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].getTemperatureIndependentName()]

        AltBo               = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeBoundariesName()]
             
        nobs                = rPrf[idhdf].shape[0]
        n_layer             = rPrf[idhdf].shape[1]

        if statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName() in statDataCl[idhdf].HDF.keys():
            avkVMR[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName()]
            avkTC[idhdf]        = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarAvkName()]
        else:
            avkVMR[idhdf]  = np.empty([nobs,n_layer,n_layer])
            avkTC[idhdf]   = np.empty([nobs,n_layer,n_layer])
            avkVMR[idhdf].fill('nan')
            avkTC[idhdf].fill('nan')

        #----------------------------------------
        #CALCULATED AIR MASS
        #----------------------------------------
        AirMPrf[idhdf]     =  np.divide(rPrfMol[idhdf], rPrf[idhdf])*sclfct

        AirMPrf_2[idhdf]   =  np.divide(rPrfMol_2[idhdf], rPrf_2[idhdf])*sclfct

        #----------------------------------------
        #EXTRACT SINGLE ALTITUDE VECTOR
        #----------------------------------------
        if (idhdf == 'Kiruna') or (idhdf == 'Zugspitze') or (idhdf == 'Izana') or (idhdf == 'Paris'):
            alt[idhdf]          = alt[idhdf][0, :]
        else:
            alt[idhdf]          = alt[idhdf][0:n_layer]

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
        #CALCULATE SCALING FACTOR AK
        #----------------------------------------
        if statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName() in statDataCl[idhdf].HDF.keys():
            avkSCF[idhdf]  = np.zeros((nobs,n_layer,n_layer))

            for obs in range(0,nobs):
                Iapriori        = np.zeros((n_layer,n_layer))
                IaprioriInv     = np.zeros((n_layer,n_layer))
                np.fill_diagonal(Iapriori, aPrf[idhdf][obs])
                np.fill_diagonal(IaprioriInv, 1.0 / (aPrf[idhdf][obs]))
                avkSCF[idhdf][obs,:,:] = np.dot(np.dot(IaprioriInv,np.squeeze(avkVMR[idhdf][obs,:,:])),Iapriori)

            dofs[idhdf]         = np.asarray([np.trace(aki) for aki in avkSCF[idhdf]])
        else:
            avkSCF[idhdf]  = np.zeros((nobs,n_layer,n_layer))
            avkSCF[idhdf].fill('nan')


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
            rPrf[idhdf]     = np.delete(rPrf[idhdf], statDataCl[idhdf].inds, axis=0)
            rPrfMol[idhdf]  = np.delete(rPrfMol[idhdf], statDataCl[idhdf].inds, axis=0)
            aPrf[idhdf]     = np.delete(aPrf[idhdf], statDataCl[idhdf].inds, axis=0)
            aPrfMol[idhdf]  = np.delete(aPrfMol[idhdf], statDataCl[idhdf].inds, axis=0)
            avkVMR[idhdf]   = np.delete(avkVMR[idhdf], statDataCl[idhdf].inds, axis=0)
            avkSCF[idhdf]   = np.delete(avkSCF[idhdf], statDataCl[idhdf].inds, axis=0)
            avkTC[idhdf]    = np.delete(avkTC[idhdf], statDataCl[idhdf].inds, axis=0)
            AirMPrf[idhdf]  = np.delete(AirMPrf[idhdf], statDataCl[idhdf].inds, axis=0)

        except Exception as errmsg:
            print '\nError: ', errmsg


        if pColsFlg:

            #---------------------------------------------------
            #STATISTICS OF TROPOPAUSE HEIGHT BASED ON DAILY AVERAGES
            #---------------------------------------------------
            # AvgTpp       = mf.dailyAvg(dtp[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
            # AvgTpp       = AvgTpp['dailyAvg']

            # maxTpp       = np.max(AvgTpp)
            # minTpp       = np.min(AvgTpp)
            # meanTpp      = np.mean(AvgTpp)
            # stdTpp       = np.std(AvgTpp)

            #print '\nMean TPH: {0:.2f} +/- {1:.2f}'.format(meanTpp, stdTpp)

            #----------------------------------------------------
            #
            #----------------------------------------------------
            if float(Lat_i[0]) >=70.: 
                meanTpp = 8.8
                stdTpp  = 1.2
            
            elif (float(Lat_i[0]) >= 60.0) & (float(Lat_i[0]) < 70.0):
                meanTpp = 9.8
                stdTpp  = 1.3
            
            elif (float(Lat_i[0]) >= 50.0) & (float(Lat_i[0]) < 60.0):
                meanTpp = 10.9
                stdTpp  = 1.2

            elif (float(Lat_i[0]) >= 40.0) & (float(Lat_i[0]) < 50.0):
                meanTpp = 11.6
                stdTpp  = 1.6

            elif (float(Lat_i[0]) >= 30.0) & (float(Lat_i[0]) < 40.0):
                meanTpp = 12.9 #12.58
                stdTpp  = 2.4  #2.72

            elif (float(Lat_i[0]) >= 20.0) & (float(Lat_i[0]) < 30.0):
                meanTpp = 15.0
                stdTpp  = 1.3

            elif (float(Lat_i[0]) >= -25.0) & (float(Lat_i[0]) < 20.0):
                meanTpp = 16.5
                stdTpp  = 0.4

            elif (float(Lat_i[0]) >= -40.0) & (float(Lat_i[0]) < -25.0):
                meanTpp = 12.3
                stdTpp  = 2.2

            elif (float(Lat_i[0]) >= -50.0) & (float(Lat_i[0]) < -40.0):
                meanTpp = 11.1
                stdTpp  = 1.3

            elif float(Lat_i[0]) < -50:
                meanTpp = 8.8
                stdTpp  = 1.7


            partialCols  = [ [0.0, 4.0], [4.0, (meanTpp - stdTpp*2.)], [(meanTpp+stdTpp*2.), 40.] ]


            for ii, pc in enumerate(partialCols):

                inds = np.where( (alt[idhdf] >= pc[0]) & (alt[idhdf] <= pc[1])  )[0]

                #---------------------------------------------------
                #THESE SITES REPORT INCREASING ALTITUDE
                #---------------------------------------------------
                if (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):       
                

                    rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,inds], axis=1)
                    aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,inds], axis=1)

                    rPcol_2[idhdf+str(pc)]  = np.sum(rPrfMol_2[idhdf][:,inds], axis=1)

                    try:

                        rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,inds], weights=AirMPrf[idhdf][:,inds],axis=1)
                        aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,inds], weights=AirMPrf[idhdf][:,inds],axis=1)

                        rPvmr_2[idhdf+str(pc)]  = np.average(rPrf_2[idhdf][:,inds], weights=AirMPrf_2[idhdf][:,inds],axis=1)
                    
                    except Exception as errmsg:
                        rPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        rPvmr[idhdf+str(pc)][:] = float('nan')

                        rPvmr_2[idhdf+str(pc)]    = np.zeros(len(rPrfMol_2[idhdf][:,0]))
                        rPvmr_2[idhdf+str(pc)][:] = float('nan')

                        aPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        aPvmr[idhdf+str(pc)][:] = float('nan')

                else: 


                    rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,inds], axis=1)
                    aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,inds], axis=1)

                    rPcol_2[idhdf+str(pc)]  = np.sum(rPrfMol_2[idhdf][:,inds], axis=1)

                    try:

                        rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,inds], weights=AirMPrf[idhdf][:,inds],axis=1)
                        aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,inds], weights=AirMPrf[idhdf][:,inds],axis=1)

                        rPvmr_2[idhdf+str(pc)]  = np.average(rPrf_2[idhdf][:,inds], weights=AirMPrf_2[idhdf][:,inds],axis=1)

                    except Exception as errmsg:
                        rPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        rPvmr[idhdf+str(pc)][:] = float('nan')

                        aPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        aPvmr[idhdf+str(pc)][:] = float('nan')

                        rPvmr_2[idhdf+str(pc)]    = np.zeros(len(rPrfMol_2[idhdf][:,0]))
                        rPvmr_2[idhdf+str(pc)][:] = float('nan')

                if ii == 0:
                    PcolTrop1[idhdf]     = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolTropapr1[idhdf]  = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrTrop1[idhdf]     = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrTropapr1[idhdf]  = np.asarray(aPvmr[idhdf+str(pc)])

                    altbl1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altbl1[idhdf][:]    = np.asarray(alt[idhdf][inds[-1]])

                    altbl2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altbl2[idhdf][:]    = np.asarray(alt[idhdf][inds[0]])

                    WvmrTrop1_2[idhdf]     = np.asarray(rPvmr_2[idhdf+str(pc)])

                elif ii == 1:
                    PcolTrop2[idhdf]     = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolTropapr2[idhdf]  = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrTrop2[idhdf]     = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrTropapr2[idhdf]  = np.asarray(aPvmr[idhdf+str(pc)])

                    altft1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altft1[idhdf][:]    = np.asarray(alt[idhdf][inds[-1]])

                    altft2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altft2[idhdf][:]    = np.asarray(alt[idhdf][inds[0]])


                    WvmrTrop2_2[idhdf]     = np.asarray(rPvmr_2[idhdf+str(pc)])

                elif ii == 2:
                    PcolStrat[idhdf]    = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolStratapr[idhdf] = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrStrat[idhdf]    = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrStratapr[idhdf] = np.asarray(aPvmr[idhdf+str(pc)])

                    altst1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altst1[idhdf][:]    = np.asarray(alt[idhdf][inds[-1]])

                    altst2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altst2[idhdf][:]    = np.asarray(alt[idhdf][inds[0]])

                    WvmrStrat_2[idhdf]   = np.asarray(rPvmr_2[idhdf+str(pc)])


        totWvmr[idhdf]  = np.average(rPrf[idhdf], axis=1, weights=AirMPrf[idhdf])
        atotWvmr[idhdf] = np.average(aPrf[idhdf], axis=1, weights=AirMPrf[idhdf])
    
    
    clmap = 'jet'
    cm           = plt.get_cmap(clmap)
    yearsLc      = YearLocator()
    daysLc       = DayLocator()
    months       = MonthLocator()
    DateFmt      = DateFormatter('%m')
    
    
    fig, ax   = plt.subplots(2, figsize=(8, 9), sharex=True)
    fig2, ax2 = plt.subplots(figsize=(7, 6))

    OCS_all   = []
    OCS_e_all = []
    
    N2O_all   = []
    N2O_e_all = []

    OCS_trop_all = []
    OCS_trop_e_all = []

    N2O_trop_all = []
    N2O_trop_e_all = []

    lifetime_all = []

    for i, idhdf in enumerate(pltID):

        if ColFlg:

            Avg            = mf.mnthlyAvg(totClmn[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
            Dates          = Avg['dates']
            dateYearFrac   = mf.toYearFraction(Avg['dates'])
            AvgData        =  Avg['mnthlyAvg']
            std            =  Avg['std']

            Avg_2          = mf.mnthlyAvg(totClmn_2[idhdf], dates_2[idhdf], dateAxis=1, meanAxis=0)
            Dates_2        = Avg_2['dates']
            dateYearFrac_2 = mf.toYearFraction(Avg_2['dates'])
            AvgData_2      =  Avg_2['mnthlyAvg']
            std_2          =  Avg_2['std']

        else:

            Avg            = mf.mnthlyAvg(WvmrStrat[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
            Dates          = Avg['dates']
            dateYearFrac   = mf.toYearFraction(Avg['dates'])
            AvgData        =  Avg['mnthlyAvg']
            std            =  Avg['std']

            Avg_2          = mf.mnthlyAvg(WvmrStrat_2[idhdf], dates_2[idhdf], dateAxis=1, meanAxis=0)
            Dates_2        = Avg_2['dates']
            dateYearFrac_2 = mf.toYearFraction(Avg_2['dates'])
            AvgData_2      =  Avg_2['mnthlyAvg']
            std_2          =  Avg_2['std']


        AvgTrop        = mf.mnthlyAvg(WvmrTrop2[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
        OCStrop        =  AvgTrop['mnthlyAvg']
        OCStrop_e      = AvgTrop['std']
        
        AvgTrop2       = mf.mnthlyAvg(WvmrTrop2_2[idhdf], dates_2[idhdf], dateAxis=1, meanAxis=0)
        N2Otrop        =  AvgTrop2['mnthlyAvg']
        N2Otrop_e      = AvgTrop2['std']


        intrsctVals = np.intersect1d(dateYearFrac, dateYearFrac_2, assume_unique=False)
        
        inds1       = np.nonzero( np.in1d( dateYearFrac, intrsctVals, assume_unique=False ) )[0]
        inds2       = np.nonzero( np.in1d( dateYearFrac_2, intrsctVals, assume_unique=False ) )[0]

        print '\n'
        print idhdf
        #print 'Total Number of Monthly OCS = ' +str(len(dateYearFrac))
        #print 'Total Number of Monthly N2O = ' +str(len(dateYearFrac_2))
        #print 'Total Number of coincident dates between OCS and N2O = ' +str(len(intrsctVals))


        AvgData   = AvgData[inds1]
        AvgData_2 = AvgData_2[inds2]

        std       = std[inds1]
        std_2     = std_2[inds2]

        indsZero   = np.where(std <= 0.)[0]
        indsZero_2 = np.where(std_2 <= 0.)[0]

        std[indsZero] =  AvgData[indsZero]*0.05 
        std_2[indsZero_2] =  AvgData_2[indsZero_2]*0.05    
        
        Dates     = Dates[inds1]
        Dates_2   = Dates_2[inds2]

        OCS_all.extend(AvgData)
        N2O_all.extend(AvgData_2)

        OCS_e_all.extend(std)
        N2O_e_all.extend(std_2)

        OCS_trop_all.extend(OCStrop[inds1])

        N2O_trop_all.extend(N2Otrop[inds2])


        meanTropOCS   = np.nanmean(OCStrop[inds1])
        meanTropN2O   = np.nanmean(N2Otrop[inds2])


        ax[i].plot(Dates, AvgData,   linestyle='-', marker ='', color='b', label='OCS')
        ax[i].scatter(Dates, AvgData, s=35, edgecolor='k', color='b')
        ax[i].set_title(idhdf)

        axr = ax[i].twinx()

        axr.plot(Dates_2, AvgData_2,   linestyle='-', marker ='', color='r', label='N2O')
        axr.scatter(Dates_2, AvgData_2, s=35, edgecolor='k', color='r')

        if i == 0:  
            ax[i].legend(prop={'size':12}, loc=2)
            axr.legend(prop={'size':12}, loc=3)

        
        ax2.plot(AvgData, AvgData_2, linestyle='none', marker ='')
        ax2.scatter(AvgData, AvgData_2, s=35, edgecolor='k', label=idhdf)


        odr, odrErr  = mf.orthoregress(AvgData, AvgData_2, xerr= std, yerr=std_2,  InError=True)
        slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(AvgData, AvgData_2)

        slope      = float(odr[0])
        slope_e    = float(odrErr[0])
        
        intercept  = float(odr[1])
        intercept_e  = float(odrErr[1])

        ax[i].xaxis.set_tick_params(which='major',labelsize=12)
        ax[i].xaxis.set_tick_params(which='minor',labelbottom='off')
        if ColFlg: ax[i].set_ylabel('OCS [{}]'.format(TCsclfctName), fontsize=16)
        else: ax[i].set_ylabel('OCS [ppb]', fontsize=16)
        #ax[i].set_xlabel('OCS [ppt]', fontsize=16)
        ax[i].tick_params(axis='both', which='major', labelsize=14)
        ax[i].grid(True)

        if ColFlg: axr.set_ylabel('N$_2$O  [{}]'.format(TCsclfctName), fontsize=16)
        else: axr.set_ylabel('N$_2$O [ppb]', fontsize=16)

        #ax[i].set_xlabel('OCS [ppt]', fontsize=16)
        axr.tick_params(axis='both', which='major', labelsize=14)
        axr.grid(True)

        if i == 0:  
            ax[i].legend(prop={'size':12}, loc=2)
            axr.legend(prop={'size':12}, loc=3)

        lifetime = slope * (meanTropOCS / meanTropN2O) * 117.

        lifetime_e  = np.sqrt( (slope_e/slope)**2 + (20./117.)**2  +  (np.std(OCStrop[inds1])/np.mean(OCStrop[inds1]))**2 +  (np.std(N2Otrop[inds2])/np.mean(N2Otrop[inds2]))**2   ) * lifetime

        print '\nSlope: {0:.2f} +/- {1:.2f}'.format(slope, slope_e)
        print 'Intercept = {0:.3f} +/- {1:.3f}'.format(intercept, intercept_e)
        print 'R value = {0:.2f}'.format(float(r_valueln))
        print 'Trop OCS [ppb] = {0:.3f} +/- {1:.3f}'.format(np.mean(OCStrop[inds1]), np.std(OCStrop[inds1]))
        print 'Trop N2O [ppb] = {0:.3f} +/- {1:.3f}'.format(np.mean(N2Otrop[inds2]), np.std(N2Otrop[inds2]))

        print 'Lifetime = {0:.2f} +/- {1:.2f}'.format(float(lifetime), float(lifetime_e))

        lifetime_all.append(lifetime)


    
    OCS_all        = np.asarray(OCS_all)
    OCS_e_all      = np.asarray(OCS_e_all)
    N2O_e_all      = np.asarray(N2O_e_all)
    N2O_all        = np.asarray(N2O_all)

    OCS_trop_all   = np.asarray(OCS_trop_all)
    N2O_trop_all   = np.asarray(N2O_trop_all)

    odr, odrErr  = mf.orthoregress(OCS_all, N2O_all, xerr=OCS_e_all, yerr=N2O_e_all, InError=True)
    slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(OCS_all, N2O_all)

    slope      = float(odr[0])
    slope_e    = float(odrErr[0])
    
    intercept  = float(odr[1])
    intercept_e  = float(odrErr[1])

    lifetime    = slope * (np.mean(OCS_trop_all) / np.mean(N2O_trop_all)) * 117.
    lifetime_e  = np.sqrt( (slope_e/slope)**2 + (20./117.)**2  +  (np.std(OCS_trop_all)/np.mean(OCS_trop_all))**2 +  (np.std(N2O_trop_all)/np.mean(N2O_trop_all))**2   ) * lifetime
   
    
    print '\nAll'
    print '\nSlope        = {0:.2f} +/- {1:.2f}'.format(slope, slope_e)
    print 'Intercept      = {0:.3f} +/- {1:.3f}'.format(intercept, intercept_e)
    print 'R value        = {0:.2f}'.format(float(r_valueln))
    print 'Trop OCS [ppb] = {0:.3f} +/- {1:.3f}'.format(np.mean(OCS_trop_all), np.std(OCS_trop_all))
    print 'Trop N2O [ppb] = {0:.3f} +/- {1:.3f}'.format(np.mean(N2O_trop_all), np.std(N2O_trop_all))
    print 'Lifetime       = {0:.2f} +/- {1:.2f}'.format(float(lifetime), float(lifetime_e))

    lifetime_all = np.asarray(lifetime_all)

    print 'Lifetime all (Mean)   = {0:.2f} +/- {1:.2f}'.format(np.mean(lifetime_all), np.std(lifetime_all))
    print 'Lifetime all (Median) = {0:.2f} +/- {1:.2f}'.format(np.median(lifetime_all), np.std(lifetime_all))




    ax2.legend(prop={'size':10})


    ax2.xaxis.set_tick_params(which='major',labelsize=12)
    ax2.xaxis.set_tick_params(which='minor',labelbottom='off')

    if ColFlg:
        ax2.set_xlabel('OCS [{} mole/cm$^2$]'.format(TCsclfctName), fontsize=16)
        ax2.set_ylabel('N$_2$O [{} mole/cm$^2$]'.format(TCsclfctName), fontsize=16)

    else:
        ax2.set_ylabel('N$_2$O [ppb]', fontsize=16)
        ax2.set_xlabel('OCS [ppb]', fontsize=16)
    ax2.tick_params(axis='both', which='major', labelsize=14)
    #ax2.set_ylim(300, 350)
    #ax2.set_xlim(0.3, 0.5)
    ax2.grid(True)

    fig.subplots_adjust(left = 0.12, bottom=0.075, top=0.95, right = 0.9)
    fig2.subplots_adjust(left = 0.12, bottom=0.12, top=0.95, right = 0.95)


    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        pdfsav.savefig(fig2,dpi=200)
        pdfsav.close() 
    else:           
        plt.show(block= False)
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


if __name__ == "__main__":
    main()