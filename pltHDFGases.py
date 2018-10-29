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
    #locs         = ['ldr', 'ldr', 'bld']#, 'kir', 'mai', 'std']
    #locs         = ['kir', 'iza', 'bld', 'stp', 'jfj', 'wlo']
    locs         = ['bld', 'wlo', 'jfj', 'mlo', 'tab', 'tor', 'eur', 'stp', 'ldr', 'rkb', 'tsk', 'ldr', 'mai', 'std', 'zgp', 'kir', 'iza', 'par', 'bre', 'nya', 'pmb', 'alt']

    #-------------------------------------
    #ID in the HDF Files
    #-------------------------------------
    #locID        = ['ahts', 'laud_120hr', 'boulder']#, 'kiruna', 'maido' , 'stdenis']
    #locID        = ['kiruna', 'izana', 'boulder', 'st.petersburg', 'jungfraujoch', 'wollongong']
    locID        = ['boulder', 'wollongong', 'jungfraujoch', 'mauna.loa.h', 'thule', '_toronto_', '_eureka_', 'st.petersburg', 'laud_120hr', 'rikubetsu', 
                   'tsukuba', 'ahts',  'maido' , 'stdenis', 'zugspitze', 'kiruna', 'izana', 'paris', 'bremen', 'ny.alesund', 'paramaribo', 'altzomoni']                       
    #-------------------------------------
    #Names in Plots
    #-------------------------------------   
    #pltID        = [ 'AHTS', 'Lauder', 'Boulder']#, 'Kiruna', 'Maido', 'StD-Maido']
    #pltID        = ['Kiruna', 'Izana', 'Boulder', 'St Petersburg', 'Jungfraujoch', 'Wollongong']
    pltID        = ['Boulder', 'Wollongong', 'Jungfraujoch', 'Mauna Loa', 'Thule', 'Toronto', 'Eureka', 'St Petersburg', 'Lauder', 'Rikubetsu', 
                   'Tsukuba', 'AHTS', 'Maido', 'StD-Maido', 'Zugspitze', 'Kiruna', 'Izana', 'Paris', 'Bremen', 'Ny Alesund', 'Paramaribo', 'Altzomoni'] 
     
    #-------------------------------------
    #Inputs
    #-------------------------------------
    gasName        = 'ocs'

    AvgType        = 'Monthly'   #'Monthly'  'Daily'
    smthFlg        = True
    period         = 1.0



    #------
    # Flags
    #------
    saveFlg       = False                  # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
    errorFlg      = True                  # Flag to process error data
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

    iyear         = 1993   
    imonth        = 1
    iday          = 1
    fyear         = 2016
    fmonth        = 12
    fday          = 31
    
    sclfct        = 1.0E12                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName    = 'pptv'                 # Name of scale factor for labeling plots
    TCsclfct      = 1.0e15
    TCsclfctName  = 'x10$^{15}$'

    pColsFlg      = True                   #Calculate tropospheric and stratospheric columns?

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
    Group = zip(dataDir,locID, pltID)
    Group.sort(key=lambda Group: Group[2])
    pltID.sort()

    for dd, id, pl in Group:
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
        else: dd = dd

        statDataCl[pl] = dc.ReadHDFData(dd, id, gasName)

    #-------------------------------------
    # Variables from HDF files 
    #-------------------------------------
    datesJD2K    = OrderedDict()
    rPrf         = OrderedDict()   #retrieved Prf in mixing ratio
    aPrf         = OrderedDict()   #apriori Prf in mixing ratio
    rPrfMol      = OrderedDict()   #retrieved Prf partial Column (molec/cm2)
    aPrfMol      = OrderedDict()   #apriori Prf partial Column (molec/cm2)
    totClmn      = OrderedDict()   #retrieved total column (molec/cm2)
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
    avkSCF       = OrderedDict()   #Averaging kernel (scale factor)
    dofs         = OrderedDict()   #degrees of freedom
    AirMPrf      = OrderedDict()   #Airmass
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
        PcolTrop1      = OrderedDict()
        PcolTrop2     = OrderedDict()

        PcolStratapr  = OrderedDict()   #partial columns A priori
        PcolTropapr1   = OrderedDict()
        PcolTropapr2  = OrderedDict()

        WvmrStrat     = OrderedDict()   #Weighted VMR
        WvmrTrop1      = OrderedDict()
        WvmrTrop2     = OrderedDict()

        WvmrStratapr  = OrderedDict()    #Weighted VMR A priori
        WvmrTropapr1   = OrderedDict()
        WvmrTropapr2  = OrderedDict()

        rPcol         = OrderedDict() 
        aPcol         = OrderedDict()

        rPvmr         = OrderedDict()
        aPvmr         = OrderedDict()


    for ii, idhdf in enumerate(pltID):

        datesJD2K[idhdf]    = statDataCl[idhdf].HDF[statDataCl[idhdf].getDatetimeName()]
        dates[idhdf]        = dc.jdf_2_datetime(datesJD2K[idhdf])
        alt[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeName()]
        sza[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAngleSolarZenithAstronomicalName()]
        
        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()+'VAR_SI_CONVERSION']            
        rPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()]*float(conv[0][1])*sclfct
        aPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAprioriName()]*float(conv[0][1])*sclfct

        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()+'VAR_SI_CONVERSION']
        rPrfMol[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()]*float(conv[0][1])*(6.02e23/100./100.)
        aPrfMol[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarAprioriName()]*float(conv[0][1])*(6.02e23/100./100.)

        conv                = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()+'VAR_SI_CONVERSION']
        totClmn[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()]*float(conv[0][1]) * (6.02e23) /100./100. / TCsclfct
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
        #PROFFIT IS A LEVEL-BASED & SFIT IS LAYER-BASED: THEY DEFINE THE ALTITUDE AND BOUNDARY VARIABLES DIFFERENT
        #----------------------------------------
        # if (idhdf == 'Kiruna') or (idhdf == 'Zugspitze') or (idhdf == 'Izana') or (idhdf == 'Paris'):

        #     alt[idhdf]          = alt[idhdf][0, :]
        #     AltBo2              = list(AltBo[0,0,:])
        #     AltBo2              = np.asarray(AltBo2)
  
        #     alt_orig[idhdf]     = alt[idhdf]
        #     dz                  = np.absolute((AltBo2[1:] - AltBo2[:-1]))*1000.*100.
        #     alt2                = np.absolute((AltBo2[1:]  + AltBo2[:-1])*0.5)
        #     rPrf[idhdf]         = np.absolute((rPrf[idhdf][:,1:]  + rPrf[idhdf][:,:-1])*0.5)
        #     TempPrf             = np.absolute((TempPrf[:,1:]  + TempPrf[:,:-1])*0.5)
        #     PresPrf             = np.absolute((PresPrf[:,1:]  + PresPrf[:,:-1])*0.5)
        #     aPrf[idhdf]         = np.absolute((aPrf[idhdf][:,1:]  + aPrf[idhdf][:,:-1])*0.5)

        #     alt[idhdf]          = alt2
            
        #     #rho                 = (PresPrf*100.00) / ((TempPrf)*8.314) * 6.022e23 /100. /100. /100.
        #     #AirMPrf[idhdf]      = rho*dz

        # else:
        #     alt[idhdf]          = alt[idhdf][0:n_layer] 
        #     AltBo2              = list(AltBo[0,:])
        #     AltBo2.extend([AltBo[1,:][-1]])
        #     AltBo2              = np.asarray(AltBo2)

        #     alt_orig[idhdf]     = alt[idhdf]
        #     dz                  = np.absolute((AltBo2[1:] - AltBo2[:-1]))*1000.*100.
        #     alt2                = np.absolute((AltBo2[1:]  + AltBo2[:-1])*0.5)     
        #     #rho                 = (PresPrf*100.00) / ((TempPrf)*8.314) * 6.022e23 /100. /100. /100.
        #     #AirMPrf[idhdf]      = rho*dz

        #rho                 = (PresPrf*100.00) / ((TempPrf)*8.314) * 6.022e23 /100. /100. /100.
        #AirMPrf[idhdf]      = rho*dz

        #rPrfMol[idhdf]      = rPrf[idhdf]*AirMPrf[idhdf] / sclfct
        #aPrfMol[idhdf]      = aPrf[idhdf]*AirMPrf[idhdf] / sclfct


        #----------------------------------------
        #INTERPOLATE ST DENIS TO MAIDO ALTITUDE
        #----------------------------------------
        # if idhdf == 'St Denis':
        #     rPrf[idhdf]    = np.asarray(interpolate.interp1d(alt[idhdf], rPrf[idhdf], axis=1)(alt_maido))
        #     aPrf[idhdf]    = np.asarray(interpolate.interp1d(alt[idhdf], aPrf[idhdf], axis=1)(alt_maido))
        #     AirMPrf[idhdf] = np.asarray(interpolate.interp1d(alt[idhdf], AirMPrf[idhdf], axis=1)(alt_maido))
        #     rPrfMol[idhdf] = np.asarray(interpolate.interp1d(alt[idhdf], rPrfMol[idhdf], axis=1)(alt_maido))
        #     aPrfMol[idhdf] = np.asarray(interpolate.interp1d(alt[idhdf], aPrfMol[idhdf], axis=1)(alt_maido))

        #     alt[idhdf]     = np.asarray(alt_maido)

        
        
        #----------------------------------------
        #CREATE A FILE WITH DATE AND TIME (TO SEND AND GET BACK THE TROPOPAUSE HEIGHT)
        #----------------------------------------
        # dirc     = Group[ii]
        # Fileout  = dirc[0] + 'Dates_'+idhdf+'.ascii'
        # with open(Fileout, 'wb') as fopen:
        #     fopen.write('#Location:   {0:15}\n'.format(idhdf))
        #     fopen.write('#Latitude:   {0:8.4f} [deg, positive North]\n'.format(float(Lat_i[0])))
        #     fopen.write('#Longitude:  {0:8.4f} [deg, positive East]\n'.format(float(Lon_i[0])))
        #     fopen.write('#Altitude:   {0:6.4f} [km]\n'.format(float(alt_instru[0])))
        #     fopen.write('YYYYMMDD     hhmmss [UT]\n')
        #     for dd in dates[idhdf]:
        #         YYYYMMDD = '{0:4d}{1:02d}{2:02d}'.format(dd.year, dd.month, dd.day)
        #         hhmmss   = '{0:02d}:{1:02d}:{2:02d}'.format(dd.hour, dd.minute, dd.second)
        #         fopen.write('{0:13}{1:13}\n'.format(YYYYMMDD, hhmmss))
        # exit()
        #----------------------------------------
        if pColsFlg:

            try:
                dirc     = Group[ii]
                #Filein   = dirc[0] + 'Dates_'+idhdf+'_dtp.ascii'
                Filein   = dirc[0] + 'Dates_'+idhdf+'_dNCEP.ascii'

                
                with open(Filein,'r') as fopen:
                    lines       = fopen.readlines()
                    YYYY        = np.array([int(x[0:4]) for x in lines[5:]])
                    MM          = np.array([int(x[4:6]) for x in lines[5:]])
                    DD          = np.array([int(x[6:8]) for x in lines[5:]])
                    HH          = np.array([int(x[13:15]) for x in lines[5:]])
                    MI          = np.array([int(x[16:18]) for x in lines[5:]])
                    SS          = np.array([int(x[19:21]) for x in lines[5:]])
                    datesdtp[idhdf]     = [dt.datetime(ye, mo, da, ho, mi, se) for ye, mo, da, ho, mi, se in zip(YYYY, MM, DD, HH, MI, SS)  ] 
                    dtp[idhdf]  = np.array([float(x[26:33]) for x in lines[5:]])
                    #dtp[idhdf]  = np.array([float(x[27:33]) for x in lines[5:]])


            except:
                print 'Missing Dynamical TH file for '+idhdf+', a dtp will be assumed = '+str(Adth)
                datesdtp[idhdf] = dates[idhdf]
                datesdtp[idhdf] = np.asarray(datesdtp[idhdf])
                dtp[idhdf]      = np.zeros(len(datesdtp[idhdf]))
                dtp[idhdf][:]   = Adth

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
        #OBTAIN ERROR VARIABLES
        #---------------------------------------- 
        if errorFlg:                               
            tot_rnd[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintyRandomName()]
            tot_sys[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintySystematicName()]
            tot_std[idhdf]      = np.sqrt(tot_rnd[idhdf]**2 + tot_sys[idhdf]**2)

            npnts               = np.shape(tot_std[idhdf])[0]
            nlvls               = np.shape(alt[idhdf])[0]
            vmr_rnd_err[idhdf]  = np.zeros((npnts,nlvls))
            vmr_sys_err[idhdf]  = np.zeros((npnts,nlvls))

            for i in range(npnts):
                conv    = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()+'VAR_SI_CONVERSION']  
                cov_rnd = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()]
                cov_sys = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintySystematicName()]

                #vmr_rnd_err[idhdf][i,:] = np.diag(cov_rnd[i][:,:])*float(conv[0][1])*sclfct**2
                #vmr_sys_err[idhdf][i,:] = np.diag(cov_sys[i][:,:])*float(conv[0][1])*sclfct**2

                vmr_rnd_err[idhdf][i,:] = np.diag(cov_rnd[i][:,:])*float(conv[0][1])
                vmr_sys_err[idhdf][i,:] = np.diag(cov_sys[i][:,:])*float(conv[0][1])

            vmr_tot_err[idhdf]  = np.sqrt(vmr_rnd_err[idhdf]+ vmr_sys_err[idhdf]) *sclfct
            vmr_rnd_err[idhdf]  = np.sqrt(vmr_rnd_err[idhdf])*sclfct
            vmr_sys_err[idhdf]  = np.sqrt(vmr_sys_err[idhdf])*sclfct

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
            dtp[idhdf]      = np.delete(dtp[idhdf], statDataCl[idhdf].inds)
            datesdtp[idhdf] = np.delete(datesdtp[idhdf], statDataCl[idhdf].inds)

        if errorFlg:
            try:
                vmr_rnd_err[idhdf]  = np.delete(vmr_rnd_err[idhdf],statDataCl[idhdf].inds,axis=0)
                vmr_sys_err[idhdf]  = np.delete(vmr_sys_err[idhdf],statDataCl[idhdf].inds,axis=0)  
                vmr_tot_err[idhdf]  = np.delete(vmr_tot_err[idhdf],statDataCl[idhdf].inds,axis=0)

                tot_rnd[idhdf]      = np.delete(tot_rnd[idhdf],statDataCl[idhdf].inds)
                tot_sys[idhdf]      = np.delete(tot_sys[idhdf],statDataCl[idhdf].inds)
                tot_std[idhdf]      = np.delete(tot_std[idhdf],statDataCl[idhdf].inds)
            except Exception as errmsg:
                print '\nError: ', errmsg


        if pColsFlg:

            #---------------------------------------------------
            #STATISTICS OF TROPOPAUSE HEIGHT BASED ON DAILY AVERAGES
            #---------------------------------------------------
            AvgTpp       = mf.dailyAvg(dtp[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
            AvgTpp       = AvgTpp['dailyAvg']

            maxTpp       = np.max(AvgTpp)
            minTpp       = np.min(AvgTpp)
            meanTpp      = np.mean(AvgTpp)
            stdTpp       = np.std(AvgTpp)

            print '\nMean TPH: {0:.2f} +/- {1:.2f}'.format(meanTpp, stdTpp)

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

            alttpp[idhdf]       = np.zeros((len(rPrfMol[idhdf][:,0])))
            alttpp2[idhdf]      = np.zeros((len(rPrfMol[idhdf][:,0])))

            ind1         = mf.nearestind(partialCols[1][1], alt[idhdf])
            ind2         = mf.nearestind(partialCols[2][0], alt[idhdf])

            alttpp[idhdf][:]      = alt[idhdf][ind1]
            alttpp2[idhdf][:]     = alt[idhdf][ind2]

            for ii, pc in enumerate(partialCols):

                ind1         = mf.nearestind(pc[0], alt[idhdf])
                ind2         = mf.nearestind(pc[1], alt[idhdf])

                #---------------------------------------------------
                #THESE SITES REPORT INCREASING ALTITUDE
                #---------------------------------------------------
                if (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):       
                    
                    rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,ind1:ind2], axis=1)
                    aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,ind1:ind2], axis=1)

                    try:
                        rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,ind1:ind2], weights=AirMPrf[idhdf][:,ind1:ind2],axis=1)
                        aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,ind1:ind2], weights=AirMPrf[idhdf][:,ind1:ind2],axis=1)
                    
                    except Exception as errmsg:
                        rPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        rPvmr[idhdf+str(pc)][:] = float('nan')

                        aPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        aPvmr[idhdf+str(pc)][:] = float('nan')

                else: 
                    rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,ind2:ind1], axis=1)
                    aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,ind2:ind1], axis=1)

                    try:
                        rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,ind2:ind1], weights=AirMPrf[idhdf][:,ind2:ind1],axis=1)
                        aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,ind2:ind1], weights=AirMPrf[idhdf][:,ind2:ind1],axis=1)

                    except Exception as errmsg:
                        rPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        rPvmr[idhdf+str(pc)][:] = float('nan')

                        aPvmr[idhdf+str(pc)]    = np.zeros(len(rPrfMol[idhdf][:,0]))
                        aPvmr[idhdf+str(pc)][:] = float('nan')

                if ii == 0:
                    PcolTrop1[idhdf]     = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolTropapr1[idhdf]  = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrTrop1[idhdf]     = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrTropapr1[idhdf]  = np.asarray(aPvmr[idhdf+str(pc)])

                    altbl1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altbl1[idhdf][:]    = np.asarray(alt[idhdf][ind1])

                    altbl2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altbl2[idhdf][:]    = np.asarray(alt[idhdf][ind2])

                elif ii == 1:
                    PcolTrop2[idhdf]     = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolTropapr2[idhdf]  = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrTrop2[idhdf]     = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrTropapr2[idhdf]  = np.asarray(aPvmr[idhdf+str(pc)])

                    altft1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altft1[idhdf][:]    = np.asarray(alt[idhdf][ind1])

                    altft2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altft2[idhdf][:]    = np.asarray(alt[idhdf][ind2])

                elif ii == 2:
                    PcolStrat[idhdf]    = np.asarray(rPcol[idhdf+str(pc)])/TCsclfct
                    PcolStratapr[idhdf] = np.asarray(aPcol[idhdf+str(pc)])/TCsclfct

                    WvmrStrat[idhdf]    = np.asarray(rPvmr[idhdf+str(pc)])
                    WvmrStratapr[idhdf] = np.asarray(aPvmr[idhdf+str(pc)])

                    altst1[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altst1[idhdf][:]    = np.asarray(alt[idhdf][ind1])

                    altst2[idhdf]       = np.zeros(len(rPrfMol[idhdf][:,0]))
                    altst2[idhdf][:]    = np.asarray(alt[idhdf][ind2])


        totWvmr[idhdf]  = np.average(rPrf[idhdf], axis=1, weights=AirMPrf[idhdf])
        atotWvmr[idhdf] = np.average(aPrf[idhdf], axis=1, weights=AirMPrf[idhdf])

    #----------------------------
    #CONCATENATE st denis and maido
    #---------------------------- 
    dates2        = {}
    
    PcolTrop12     = {}
    PcolTrop22     = {}
    PcolStrat2    = {}
    
    PcolTropapr12  = {}
    PcolTropapr22  = {}
    PcolStratapr2 = {}

    WvmrTrop12     = {}
    WvmrTrop22     = {}
    WvmrStrat2    = {}
    
    WvmrTropapr12  = {}
    WvmrTropapr22  = {}
    WvmrStratapr2 = {}

    totClmn2        = {}

    dtp2           = {}
    alttpp12        = {}
    alttpp22       = {}

    pltID2        = []
    Lat2          = []
    Lon2          = []

    for i, idhdf in enumerate(pltID):

        if idhdf == 'Maido': continue
        
        if idhdf  == 'StD-Maido':
            dates2[idhdf]          = np.concatenate( (dates[idhdf], dates['Maido']))

            totClmn2[idhdf]       = np.concatenate( (totClmn[idhdf], totClmn['Maido']))


            PcolTrop12[idhdf]       = np.concatenate( (PcolTrop1[idhdf], PcolTrop1['Maido']))
            PcolTrop22[idhdf]       = np.concatenate( (PcolTrop2[idhdf], PcolTrop2['Maido']))
            PcolStrat2[idhdf]      = np.concatenate( (PcolStrat[idhdf], PcolStrat['Maido']))
            
            PcolTropapr12[idhdf]    = np.concatenate( (PcolTropapr1[idhdf], PcolTropapr1['Maido']))
            PcolTropapr22[idhdf]    = np.concatenate( (PcolTropapr2[idhdf], PcolTropapr2['Maido']))
            PcolStratapr2[idhdf]   = np.concatenate( (PcolStratapr[idhdf], PcolStratapr['Maido']))

            WvmrTrop12[idhdf]       = np.concatenate( (WvmrTrop1[idhdf], WvmrTrop1['Maido']))
            WvmrTrop22[idhdf]       = np.concatenate( (WvmrTrop2[idhdf], WvmrTrop2['Maido']))
            WvmrStrat2[idhdf]      = np.concatenate( (WvmrStrat[idhdf], WvmrStrat['Maido']))
            
            WvmrTropapr12[idhdf]    = np.concatenate( (WvmrTropapr1[idhdf], WvmrTropapr1['Maido']))
            WvmrTropapr22[idhdf]    = np.concatenate( (WvmrTropapr2[idhdf], WvmrTropapr2['Maido']))
            WvmrStratapr2[idhdf]   = np.concatenate( (WvmrStratapr[idhdf], WvmrStratapr['Maido']))

            dtp2[idhdf]              = np.concatenate( (dtp[idhdf], dtp['Maido']))
            alttpp12[idhdf]           = np.concatenate( (alttpp[idhdf], alttpp['Maido']))
            alttpp22[idhdf]           = np.concatenate( (alttpp2[idhdf], alttpp2['Maido']))

            pltID2.append(pltID[i])
            Lat2.append(Lat[i])
            Lon2.append(Lon[i])
        else:
            dates2[idhdf]         = dates[idhdf]
            totClmn2[idhdf]       = totClmn[idhdf]
            PcolTrop12[idhdf]     = PcolTrop1[idhdf] 
            PcolTrop22[idhdf]     = PcolTrop2[idhdf] 
            PcolStrat2[idhdf]     = PcolStrat[idhdf]
            PcolTropapr12[idhdf]   = PcolTropapr1[idhdf]
            PcolTropapr22[idhdf]   = PcolTropapr2[idhdf]
            PcolStratapr2[idhdf]  = PcolStratapr[idhdf]

            WvmrTrop12[idhdf]      = WvmrTrop1[idhdf] 
            WvmrTrop22[idhdf]      = WvmrTrop2[idhdf] 
            WvmrStrat2[idhdf]     = WvmrStrat[idhdf]
            WvmrTropapr12[idhdf]   = WvmrTropapr1[idhdf]
            WvmrTropapr22[idhdf]   = WvmrTropapr2[idhdf]
            WvmrStratapr2[idhdf]  = WvmrStratapr[idhdf]

            dtp2[idhdf]            = dtp[idhdf]
            alttpp12[idhdf]        = alttpp[idhdf]
            alttpp22[idhdf]        = alttpp2[idhdf]

            pltID2.append(pltID[i])
            Lat2.append(Lat[i])
            Lon2.append(Lon[i])

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                           PLOTS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------
    print '\nPrinting Plots.......\n'
    if saveFlg: pdfsav = PdfPages(pltFile)
    else: pdfsav = ''
    clr = mf.clrplt()

    #---------------------------------------------------
    # Determine if multiple years
    #---------------------------------------------------
    years  = []
    Nyears = []

    for i, idhdf in enumerate(pltID):
        years.append([ singDate.year for singDate in dates[idhdf]] )
    years = np.asarray(years)

    Nyears = [len(list(set(y))) for y in years]
    indY = np.max(Nyears)
    
    if indY > 1: yrsFlg = True         
    else:        yrsFlg = False

    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)    
    yearsLc      = YearLocator()
    monthsAll    = MonthLocator()
    #months       = MonthLocator(bymonth=1,bymonthday=1)
    months       = MonthLocator()

    if yrsFlg: DateFmt      = DateFormatter('%Y')
    else:      DateFmt      = DateFormatter('%m\n%Y')


    #---------------------------------------------------
    # Order data based on +Lat to -Lat
    #---------------------------------------------------
    pltID = [y for (x,y) in sorted(zip(Lat,pltID), reverse=True)]
    Lon = [y for (x,y) in sorted(zip(Lat,Lon), reverse=True)]
    Lat   = sorted(Lat, reverse=True)
    pltID = np.asarray(pltID)

    pltID2 = [y for (x,y) in sorted(zip(Lat2,pltID2), reverse=True)]
    Lon2 = [y for (x,y) in sorted(zip(Lat2,Lon2), reverse=True)]
    Lat2   = sorted(Lat2, reverse=True)
    pltID2 = np.asarray(pltID2)

    
    #---------------------------------------------------
    # Map with Sites
    #---------------------------------------------------
    fig = plt.figure(figsize=(11,7))

    eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0, lat_0=0, lon_0=-130)
    eq_map.drawcoastlines()
    eq_map.drawcountries()
    eq_map.fillcontinents(color = 'gray', alpha=0.5)
    eq_map.drawmapboundary()
    eq_map.drawmeridians(np.arange(0, 360, 30))
    eq_map.drawparallels(np.arange(-90, 90, 30))
     
    for lo, la, idhdf in zip(Lon, Lat, pltID):
        x,y = eq_map(lo, la)
        eq_map.plot(x, y, 'yo', markersize=10.0)
        plt.text(x+10000,y-800000, idhdf, fontsize=12, color='red')

    plt.title("IRWG-NDACC Sites participating in the OCS global study")

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else: 
        plt.show(block=False)
       #user_input = raw_input('Press any key to exit >>> ')
       #sys.exit()

    #---------------------------------------------------
    # Defining variable for plots
    #---------------------------------------------------
    npanels   = len(pltID)
    npanels2  = int(math.ceil(npanels/4.0))
    xmin      = dt.date(iyear, imonth, iday)
    xmax      = dt.date(fyear, fmonth, fday)
    clmap     = 'jet'

    #---------------------------------------------------
    # Mean vertical profiles
    #---------------------------------------------------
    fig = plt.figure(figsize=(12,12))

    outer_grid = gridspec.GridSpec(npanels2, 4, wspace=0.11, hspace=0.085)

    for i, idhdf in enumerate(pltID):

        ax = plt.Subplot(fig, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_axis_bgcolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_axis_bgcolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_axis_bgcolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_axis_bgcolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_axis_bgcolor('lightgrey')
        
        else:
            ax.set_axis_bgcolor('lightgrey')

        if int(nobs) > 1:

            prfMean    = np.nanmean(rPrf[idhdf],axis=0)
            prfSTD     = np.nanstd(rPrf[idhdf],axis=0)
            avkVMRAv   = np.nanmean(avkVMR[idhdf], axis=0)
            avkSCFAv   = np.nanmean(avkSCF[idhdf], axis=0)
            ax.plot(prfMean,alt[idhdf],label='Retrieved')
            ax.fill_betweenx(alt[idhdf],prfMean-prfSTD,prfMean+prfSTD,alpha=0.25)

            dtpMean    = np.nanmean(dtp[idhdf])
            dtpStd     = np.nanstd(dtp[idhdf])

            ax.axhline(y=dtpMean, linewidth=1.5, color='r', alpha=0.5)  
            ax.axhline(y=dtpMean - dtpStd*2., linewidth=1.5, color='r', alpha=0.5)  
            ax.axhline(y=dtpMean + dtpStd*2, linewidth=1.5, color='r', alpha=0.5)

            altftmean  =   np.nanmean(alttpp[idhdf])
            altstmean  =   np.nanmean(alttpp2[idhdf])

            ax.axhline(y=altftmean, linewidth=1.5, color='green', alpha=0.5) 
            ax.axhline(y=altstmean, linewidth=1.5, color='green', alpha=0.5)

        else:
            ax.plot(rPrf[0],alt,label='Retrieved')


        ax.plot(aPrf[idhdf][0],alt[idhdf],linestyle='--', label='-Apriori', color='r')

        ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

        ax.grid(True,which='both')    
        #ax.legend(prop={'size':10})    
        #ax.set_ylabel('Altitude [km]')
        ax.tick_params(which='both',labelsize=10)
        ax.set_ylim(0,50)
        ax.set_xlim(0,600) 
        #ax.set_title('Mean profile of'+ gasName.upper())
      
        fig.add_subplot(ax)

    all_axes = fig.get_axes()
    #show only the outside spines
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
        if ax.is_last_row():
            ax.spines['bottom'].set_visible(True)
            plt.setp(ax.get_xticklabels(), visible=True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
        if ax.is_last_col():
            ax.spines['right'].set_visible(True)

    if (npanels % 2 == 1): #even

        all_axes[-2].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        all_axes[-2].set_zorder(1)

    all_axes[-1].set_xlabel('VMR ['+sclfctName+']')
    all_axes[-2].set_xlabel('VMR ['+sclfctName+']')

    fig.text(0.03, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
    plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
    #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
    fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

    
    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else: 
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()

    #---------------------------------------------------
    # Mean Error vertical profiles
    #---------------------------------------------------
    fig = plt.figure(figsize=(12,12))

    outer_grid = gridspec.GridSpec(npanels2, 4, wspace=0.11, hspace=0.085)

    for i, idhdf in enumerate(pltID):

        ax = plt.Subplot(fig, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_axis_bgcolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_axis_bgcolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_axis_bgcolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_axis_bgcolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_axis_bgcolor('lightgrey')
        
        else:
            ax.set_axis_bgcolor('lightgrey')

        if int(nobs) > 1:

            vmr_rnd_err_Mean  = np.nanmean(vmr_rnd_err[idhdf],axis=0)   
            vmr_sys_err_Mean  = np.nanmean(vmr_sys_err[idhdf],axis=0)
            vmr_tot_err_Mean  = np.nanmean(vmr_tot_err[idhdf],axis=0)

            vmr_tot_err_std   = np.nanstd(vmr_tot_err[idhdf],axis=0)

           
            ax.plot(vmr_rnd_err_Mean,alt[idhdf],label='Random Error', color='r')
            ax.plot(vmr_sys_err_Mean,alt[idhdf],label='Systematic Error', color='b')
            ax.plot(vmr_tot_err_Mean,alt[idhdf],label='Total Error', color='k')
            #ax.fill_betweenx(alt[idhdf],vmr_tot_err_Mean-vmr_tot_err_std,vmr_tot_err_Mean+vmr_tot_err_std,alpha=0.25, color='k')
          
        else:
            ax.plot(rPrf[0],alt,label='Retrieved')

        ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

        ax.grid(True,which='both')    
        #ax.legend(prop={'size':10})    
        #ax.set_ylabel('Altitude [km]')
        ax.tick_params(which='both',labelsize=10)
        ax.set_ylim(0,50)
        ax.set_xlim(0,60) 
        #ax.set_title('Mean profile of'+ gasName.upper())
      
        fig.add_subplot(ax)

    all_axes = fig.get_axes()
    #show only the outside spines
    for ax in all_axes:
        for sp in ax.spines.values():
            sp.set_visible(False)
            plt.setp(ax.get_xticklabels(), visible=False)
        if ax.is_first_row():
            ax.spines['top'].set_visible(True)
        if ax.is_last_row():
            ax.spines['bottom'].set_visible(True)
            plt.setp(ax.get_xticklabels(), visible=True)
        if ax.is_first_col():
            ax.spines['left'].set_visible(True)
        if ax.is_last_col():
            ax.spines['right'].set_visible(True)

    if (npanels % 2 == 1): #even

        all_axes[-2].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        all_axes[-2].set_zorder(1)

    all_axes[-1].set_xlabel('VMR ['+sclfctName+']')
    all_axes[-2].set_xlabel('VMR ['+sclfctName+']')

    fig.text(0.03, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
    plt.suptitle('{} Uncertainty Vertical Profiles'.format(gasName.upper()), fontsize=16  )
    #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
    fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

    
    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else: 
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()

    #----------------------------
    #CONCATENATE st denis and maido
    #---------------------------- 
    # Lat = [li for i, li in enumerate(Lat) if (pltID[i] != 'Maido') & (pltID[i] != 'St Denis')]
    # Lon = [li for i, li in enumerate(Lon) if (pltID[i] != 'Maido') & (pltID[i] != 'St Denis')]

    # pltID = [i for i in pltID if (i != 'Maido') & (i != 'St Denis')]

    # pltID.append('MaSDenis')
    # Lat.append(-21.08)
    # Lon.append(55.38)

    # dates['MaSDenis']         = np.concatenate( (dates['Maido'], dates['St Denis']))

    # PcolTrop1['MaSDenis']     = np.concatenate( (PcolTrop1['St Denis'], PcolTrop1['Maido']))
    # PcolTrop2['MaSDenis']     = np.concatenate( (PcolTrop2['St Denis'], PcolTrop2['Maido']))
    # PcolStrat['MaSDenis']     = np.concatenate( (PcolStrat['St Denis'], PcolStrat['Maido']))
   
    # PcolTropapr1['MaSDenis']  = np.concatenate( (PcolTropapr1['St Denis'], PcolTropapr1['Maido']))
    # PcolTropapr2['MaSDenis']  = np.concatenate( (PcolTropapr2['St Denis'], PcolTropapr2['Maido']))
    # PcolStratapr['MaSDenis']  = np.concatenate( (PcolStratapr['St Denis'], PcolStratapr['Maido']))

    # WvmrTrop1['MaSDenis']     = np.concatenate( (WvmrTrop1['St Denis'], WvmrTrop1['Maido']))
    # WvmrTrop2['MaSDenis']     = np.concatenate( (WvmrTrop2['St Denis'], WvmrTrop2['Maido']))
    # WvmrStrat['MaSDenis']     = np.concatenate( (WvmrStrat['St Denis'], WvmrStrat['Maido']))
    
    # WvmrTropapr1['MaSDenis']  = np.concatenate( (WvmrTropapr1['St Denis'], WvmrTropapr1['Maido']))
    # WvmrTropapr2['MaSDenis']  = np.concatenate( (WvmrTropapr2['St Denis'], WvmrTropapr2['Maido']))
    # WvmrStratapr['MaSDenis']  = np.concatenate( (WvmrStratapr['St Denis'], WvmrStratapr['Maido']))

    # dtp['MaSDenis']  = np.concatenate( (dtp['St Denis'], dtp['Maido']))

    # alttpp['MaSDenis']  = np.concatenate( (alttpp['St Denis'], alttpp['Maido']))
    # alttpp2['MaSDenis']  = np.concatenate( (alttpp2['St Denis'], alttpp2['Maido']))

    npanels2  = int(math.ceil(npanels/3.0))
    print npanels2

    #---------------------------------------------------
    # Tropopause Height
    #---------------------------------------------------
    print '\nPlot: Tropopause height:\n' 
    ##resTH = AnalTS(npanels2-1, dates2, dtp2, pltID2, Lat2, fits = True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Tropopause height', unitsStr=' km', yData2=alttpp12, yData3=alttpp22, yData4=dtp2, ymin=4 ,ymax=19)
    resTH = AnalTS(npanels2-1, dates2, dtp2, pltID2, Lat2, fits = True, AvgType=AvgType, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Tropopause height', unitsStr=' km', ymin=4 ,ymax=19)
    resTH = np.asarray(resTH)
    
    if pltPcol:
        
        #---------------------------------------------------
        # Time Series of Total Columns (All points)
        # Res ==> slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp 
        #---------------------------------------------------
        #print '\nPlot: All Total Columns:\n'
        #res    = AnalTS(npanels2, dates, totClmn, pltID, Lat, fits=False, AvgType='none', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Total Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$')

        #---------------------------------------------------
        # Time Series of Total Columns (multiple panels) -- Averaged values
        #---------------------------------------------------
        print '\nPlot: Averaged Total Columns:\n' 
        resTC = AnalTS(npanels2-1, dates2, totClmn2, pltID2, Lat2, fits=True, AvgType=AvgType, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Total Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=5.0 ,ymax=11)
        resTC = np.asarray(resTC)

        #---------------------------------------------------
        # Time Series of Total Columns Apriori (multiple panels) -- Averaged values
        #---------------------------------------------------

        #print '\nPlot: Averaged Total Columns Apriori:\n' 
        #resTCa = AnalTS(npanels2, dates, atotClmn, pltID, Lat, fits=True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Total Column (Apriori)', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=5.0 ,ymax=11)
        #resTCa = np.asarray(resTCa)
            
        if pColsFlg:

            #---------------------------------------------------
            # Boundary Layer Columns ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Boundary Layer Column:\n' 
            resLC = AnalTS(npanels2-1, dates2, PcolTrop12, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Low Tropospheric Partial Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', period=period)#, ymin=0.9 ,ymax=6.0)
            resLC = np.asarray(resLC)

            print '\nPlot: Boundary Layer Column Anomalies:\n' 
            resLCAnom = AnalTSAnom(npanels2-1, dates2, PcolTrop12, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Low Tropospheric Partial Column Anomalies', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$',period=period)
            resLCAnom = np.asarray(resLCAnom)

            #---------------------------------------------------
            # Free Tropospheric Columns ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Free Tropospheric Column:\n' 
            resLC2 = AnalTS(npanels2-1, dates2, PcolTrop22, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Free Tropospheric Partial Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', period=period)#, ymin=0.9, ymax=6.0)
            resLC2 = np.asarray(resLC2)

            print '\nPlot: Boundary Layer Column Anomalies:\n' 
            resLC2Anom = AnalTSAnom(npanels2-1, dates2, PcolTrop22, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Free Tropospheric Partial Column Anomalies', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', period=period)
            resLC2Anom = np.asarray(resLC2Anom)

            #---------------------------------------------------
            # Stratospheric Columns ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Stratospheric Column:\n' 
            resSC = AnalTS(npanels2-1, dates2, PcolStrat2, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Partial Column', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', period=period)#, ymin=0.001 ,ymax=2.5)
            resSC = np.asarray(resSC)

            print '\nPlot: Stratospheric Column Anomalies:\n' 
            resSCAnom = AnalTSAnom(npanels2-1, dates2, PcolStrat2, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Partial Column Anomalies', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', period=period)
            resSCAnom = np.asarray(resSCAnom)

            #---------------------------------------------------
            # Tropospheric Columns ==> Apriori
            #---------------------------------------------------
            #print '\nPlot: Tropospheric Column Apriori:\n' 
            #resLCa = AnalTS(npanels2, dates, PcolTropapr, pltID, Lat, fits = True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Tropospheric Partial Column (Apriori)', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=2.0 ,ymax=9.5)
            #resLCa = np.asarray(resLCa)

            #---------------------------------------------------
            # Stratospheric Columns ==> Apriori
            #---------------------------------------------------
            #print '\nPlot: Stratospheric Column Apriori:\n' 
            #resSCa = AnalTS(npanels2, dates, PcolStratapr, pltID, Lat, fits = True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Partial Column (Apriori)', unitsStr=TCsclfctName+' molecules$\cdot$cm$^{-2}$', ymin=0.001 ,ymax=2.5)
            #resSCa = np.asarray(resSCa)

            #---------------------------------------------------
            # Plot total columns by month
            #---------------------------------------------------
            # print '\nPlot: Monthly mean Partial Columns:\n' 
            # fig, (ax1, ax2, ax3)  = plt.subplots(3, figsize=(10, 10))

            # colormap = plt.cm.gist_ncar
            # ax1.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])
            # ax2.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])
            # ax3.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])
            # #colors = [colormap(i) for i in np.linspace(0, 0.9,len(pltID))]
            # #clrs = plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])

            # for i, idhdf in enumerate(pltID):

            #     month    = np.array([d.month for d in dates[idhdf]])
            #     mnthSort = list(set(month))
                
            #     mnthMean = np.zeros(len(mnthSort))
            #     mnthSTD  = np.zeros(len(mnthSort))

            #     mnthMean2 = np.zeros(len(mnthSort))
            #     mnthSTD2  = np.zeros(len(mnthSort))

            #     mnthMean3 = np.zeros(len(mnthSort))
            #     mnthSTD3  = np.zeros(len(mnthSort))
                
            #     for i,m in enumerate(mnthSort):
            #         inds        = np.where(month == m)[0]
            #         mnthMean[i] = np.mean(PcolTrop1[idhdf][inds])
            #         mnthSTD[i]  = np.std(PcolTrop1[idhdf][inds])

            #         mnthMean2[i] = np.mean(PcolTrop2[idhdf][inds])
            #         mnthSTD2[i]  = np.std(PcolTrop2[idhdf][inds])

            #         mnthMean3[i] = np.mean(PcolStrat[idhdf][inds])
            #         mnthSTD3[i]  = np.std(PcolStrat[idhdf][inds])
                    
            #     ax1.errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6, label=idhdf)
            #     ax2.errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6)
            #     ax3.errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)
            
            # ax1.grid(True,which='both')
            # ax1.set_ylabel('Low Tropospheric\n['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')
            # #ax1.set_xlabel('Month')
            # #ax1.set_title('Retrieved Monthly Mean with Standard Deviation')
            # ax1.set_xlim((0,13))
            # ax1.set_xticks(range(1,13))
            # ax1.legend(ncol=4, prop={'size':12}, loc = 'upper center', bbox_to_anchor=[0.5, 1.3])

            # ax2.grid(True,which='both')
            # ax2.set_ylabel('Free Tropospheric\n['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')
            # ax2.set_xlabel('Month')
            # #ax1.set_title('Retrieved Monthly Mean with Standard Deviation')
            # ax2.set_xlim((0,13))
            # ax2.set_xticks(range(1,13))

            # ax3.grid(True,which='both')
            # ax3.set_ylabel('Stratospheric\n['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')
            # ax3.set_xlabel('Month')
            # #ax1.set_title('Retrieved Monthly Mean with Standard Deviation')
            # ax3.set_xlim((0,13))
            # ax3.set_xticks(range(1,13))
            
            # if saveFlg: pdfsav.savefig(fig,dpi=200)
            # else: 
            #     plt.show(block=False)
            # #user_input = raw_input('Press any key to exit >>> ')
            # #sys.exit()

            #---------------------------------------------------
            # Plot total Columns by month
            #---------------------------------------------------
            print '\nPlot: Monthly mean Partial Columns:\n' 
            fig, ax  = plt.subplots(3, 5, figsize=(17, 10), sharey=False, sharex=True)
            
            for i, idhdf in enumerate(pltID2):

                month    = np.array([d.month for d in dates2[idhdf]])
                mnthSort = list(set(month))
                
                mnthMean = np.zeros(len(mnthSort))
                mnthSTD  = np.zeros(len(mnthSort))

                mnthMean2 = np.zeros(len(mnthSort))
                mnthSTD2  = np.zeros(len(mnthSort))

                mnthMean3 = np.zeros(len(mnthSort))
                mnthSTD3  = np.zeros(len(mnthSort))
                
                for j,m in enumerate(mnthSort):
                    inds        = np.where(month == m)[0]
                    mnthMean[j] = np.mean(PcolTrop12[idhdf][inds])
                    mnthSTD[j]  = np.std(PcolTrop12[idhdf][inds])

                    mnthMean2[j] = np.mean(PcolTrop22[idhdf][inds])
                    mnthSTD2[j]  = np.std(PcolTrop22[idhdf][inds])

                    mnthMean3[j] = np.mean(PcolStrat2[idhdf][inds])
                    mnthSTD3[j]  = np.std(PcolStrat2[idhdf][inds])

                Lat[i] = float(Lat[i])

                if Lat[i] >= 50.:     
                    ax[0,0].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,0].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,0].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,0].grid(True,which='both')
                    ax[1,0].grid(True,which='both')
                    ax[2,0].grid(True,which='both')

                    ax[0,0].set_ylabel('Low Tropospheric Columns ['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')
                    ax[1,0].set_ylabel('Free Tropospheric Columns ['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')
                    ax[2,0].set_ylabel('Stratospheric Columns ['+ TCsclfctName+ ' molecules cm$^{-2}$]',multialignment='center')

                    ax[0,0].set_xlim((0,13))
                    ax[0,0].set_xticks(range(1,13))
                    ax[1,0].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,0].set_xlabel('Month')
                    ax[0, 0].set_title('50$^{\circ}$N - 90$^{\circ}$N', fontsize=14)

                    #ax[0, 0].set_ylim(350, 580)
                    #ax[1, 0].set_ylim(350, 580)
                    #ax[2, 0].set_ylim(100, 350)
                    
                elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                    ax[0,1].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,1].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,1].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,1].grid(True,which='both')
                    ax[1,1].grid(True,which='both')
                    ax[2,1].grid(True,which='both')

                    ax[1,1].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,1].set_xlabel('Month')

                    ax[0,1].set_title('20$^{\circ}$N - 50$^{\circ}$N', fontsize=14)

                    #ax[0, 1].set_ylim(350, 580)
                    #ax[1, 1].set_ylim(350, 580)
                    #ax[2, 1].set_ylim(100, 350)
                   
                elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                    ax[0,2].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,2].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,2].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,2].grid(True,which='both')
                    ax[1,2].grid(True,which='both')
                    ax[2,2].grid(True,which='both')

                    ax[1,2].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,2].set_xlabel('Month')
                    ax[0,2].set_title('20$^{\circ}$S - 20$^{\circ}$N', fontsize=14)

                    #ax[0, 2].set_ylim(350, 580)
                    #ax[1, 2].set_ylim(350, 580)
                    #ax[2, 2].set_ylim(100, 350)
                    
                elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                    ax[0,3].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,3].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,3].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,3].grid(True,which='both')
                    ax[1,3].grid(True,which='both')
                    ax[2,3].grid(True,which='both')

                    ax[1,3].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,3].set_xlabel('Month')

                    ax[0,3].set_title('50$^{\circ}$S - 20$^{\circ}$S', fontsize=14)

                    #ax[0, 3].set_ylim(350, 580)
                    #ax[1, 3].set_ylim(350, 580)
                    #ax[2, 3].set_ylim(100, 350)

                    
                elif (Lat[i] < -50.):
                    ax[0,4].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,4].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,4].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,4].grid(True,which='both')
                    ax[1,4].grid(True,which='both')
                    ax[2,4].grid(True,which='both')

                    ax[1,4].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,4].set_xlabel('Month')
                    ax[0,4].set_title('90$^{\circ}$S - 50$^{\circ}$S', fontsize=14)

                    #ax[0, 4].set_ylim(350, 580)
                    #ax[1, 4].set_ylim(350, 580)
                    #ax[2, 4].set_ylim(100, 350)

            fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.98, wspace=0.15, hspace=0.15)
            
            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else: 
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

            #---------------------------------------------------
            # Bar plot: Three different periods ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Trends Columns:\n'
            hbarplt3(resLC, resLC2, resSC, pltID2, b1_label='Low Tropospheric', b2_label='Free Tropospheric', b3_label='Stratospheric', subtitle='Retrieval', saveFlg=saveFlg, pdfsav=pdfsav)
            hbarplt3(resLCAnom, resLC2Anom, resSCAnom, pltID2, b1_label='Low Tropospheric', b2_label='Free Tropospheric', b3_label='Stratospheric', subtitle='Retrieval - Anomalies', saveFlg=saveFlg, pdfsav=pdfsav)

            #---------------------------------------------------
            # Bar plot: Three different periods ==> Apriori
            #---------------------------------------------------
            #hbarplt3(resTCa, resLCa, resSCa, pltID, b1_label='Tot Column', b2_label='Trop Column', b3_label='Strat Column', subtitle='A priori', saveFlg=saveFlg, pdfsav=pdfsav)

            # #---------------------------------------------------
            # # Bar plot: Three different periods ==> Retrieval - Apriori
            # #---------------------------------------------------
                    
            ind = np.arange(len(resSC[0]))
            
            # fig, (ax, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=(15, 8.5))
            
            # ax.barh(ind-0.27, resTC[0]-resTCa[0], 0.27, xerr=np.sqrt(resTC[1]**2 + resTCa[1]**2) , align='center', color = 'r', ecolor = 'k', label = 'Tot Column')
            # ax.barh(ind, resLC[0]-resLCa[0], 0.27, xerr=np.sqrt(resLC[1]**2 + resLCa[1]**2), align='center', color = 'b', ecolor = 'k', label = 'Trop Column')
            # ax.barh(ind+0.27, resSC[0]-resSCa[0], 0.27, xerr=np.sqrt(resSC[1]**2 + resSCa[1]**2), align='center', color = 'g', ecolor = 'k', label = 'Strat Column')  #yerr=slope_TC*0
            # ax.xaxis.grid(True)
            # ax.yaxis.grid(True)
            # ax.set_xlabel('Rate of change (%/y)')
            # ax.set_yticks(ind)
            # ax.set_yticklabels(np.transpose(pltID), rotation=0)
            # ax.set_title('Retrieval - Apriori\n(years submitted)', multialignment='center')
            # #ax.set_xlabel('Site')
            # ax.set_xlim(-2.0, 2.0)
            # ax.axvline(0, color='black', lw=1)
            # ax.legend(prop={'size':8}, loc = 1)
            # #ax.invert_yaxis()

            # ax2.barh(ind-0.27, resTC[2]-resTCa[2], 0.27, xerr=np.sqrt(resTC[3]**2 + resTCa[3]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax2.barh(ind, resLC[2]-resLCa[2], 0.27, xerr=np.sqrt(resLC[3]**2 + resLCa[3]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax2.barh(ind+0.27, resSC[2]-resSCa[2], 0.27, xerr=np.sqrt(resSC[3]**2 + resSCa[3]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax2.xaxis.grid(True)
            # ax2.yaxis.grid(True)
            # ax2.set_xlabel('Rate of change (%/y))')
            # ax2.set_title('Retrieval - Apriori\n(1995 - 2002)', multialignment='center')
            # ax2.set_yticks(ind)
            # #ax2.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax2.set_xlim(-2.0, 2.0)
            # ax2.axvline(0, color='black', lw=1)

            # #ax2.invert_yaxis()
            # #ax2.yticks([])

            # ax3.barh(ind-0.27, resTC[4]-resTCa[4], 0.27, xerr=np.sqrt(resTC[5]**2 + resTCa[5]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax3.barh(ind, resLC[4]-resLCa[4], 0.27, xerr=np.sqrt(resLC[5]**2 + resLCa[5]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax3.barh(ind+0.27, resSC[4]-resSCa[4], 0.27, xerr=np.sqrt(resSC[5]**2 + resSCa[5]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax3.xaxis.grid(True)
            # ax3.yaxis.grid(True)
            # ax3.set_xlabel('Rate of change (%/y))')
            # ax3.set_title('Retrieval - Apriori\n(2002 - 2008)', multialignment='center')
            # ax3.set_yticks(ind)
            # #ax3.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax3.set_xlim(-2.0, 2.0)
            # ax3.axvline(0, color='black', lw=1)

            # #ax3.invert_yaxis()

            # ax4.barh(ind-0.27, resTC[6]-resTCa[6], 0.27, xerr=np.sqrt(resTC[7]**2 + resTCa[7]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax4.barh(ind, resLC[6]-resLCa[6], 0.27, xerr=np.sqrt(resLC[7]**2 + resLCa[7]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax4.barh(ind+0.27, resSC[6]-resSCa[6], 0.27, xerr=np.sqrt(resSC[7]**2 + resSCa[7]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax4.xaxis.grid(True)
            # ax4.yaxis.grid(True)
            # ax4.set_xlabel('Rate of change (%/y))')
            # ax4.set_title('Retrieval - Apriori\n(2008 - 2016)', multialignment='center')
            # ax4.set_yticks(ind)
            # #ax4.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax4.set_xlim(-2.0, 2.0)
            # ax4.axvline(0, color='black', lw=1)

            # #ax4.invert_yaxis()

            # plt.gca().invert_yaxis()
            # #fig.tight_layout()
            # fig.subplots_adjust(left=0.1, right=0.95)

            # if saveFlg: pdfsav.savefig(fig,dpi=200)
            # else:       plt.show(block=False)

            # for i, idhdf in enumerate(pltID):

            #     print '\nRate of change (Retrieval  - Apriori) - Columns: {}\n'.format(idhdf)

            #     print "Total Columns - all years   = {:.2f} +/- {:.2f}% ".format(resTC[0][i]-resTCa[0][i], np.sqrt(resTC[1][i]**2 + resTCa[1][i]**2))
            #     print "Total Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resTC[2][i]-resTCa[2][i], np.sqrt(resTC[3][i]**2 + resTCa[3][i]**2))
            #     print "Total Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resTC[4][i]-resTCa[4][i], np.sqrt(resTC[5][i]**2 + resTCa[5][i]**2))
            #     print "Total Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resTC[6][i]-resTCa[6][i], np.sqrt(resTC[7][i]**2 + resTCa[7][i]**2))

            #     print "Trop Columns - all years    = {:.2f} +/- {:.2f}% ".format(resLC[0][i]-resLCa[0][i], np.sqrt(resLC[1][i]**2 + resLCa[1][i]**2))
            #     print "Trop Columns - 1995 - 2002  = {:.2f} +/- {:.2f}% ".format(resLC[2][i]-resLCa[2][i], np.sqrt(resLC[3][i]**2 + resLCa[3][i]**2))
            #     print "Trop Columns - 2002 - 2008  = {:.2f} +/- {:.2f}% ".format(resLC[4][i]-resLCa[4][i], np.sqrt(resLC[5][i]**2 + resLCa[5][i]**2))
            #     print "Trop Columns - 2008 - 2016  = {:.2f} +/- {:.2f}% ".format(resLC[6][i]-resLCa[6][i], np.sqrt(resLC[7][i]**2 + resLCa[7][i]**2))

            #     print "Strat Columns - all years   = {:.2f} +/- {:.2f}% ".format(resSC[0][i]-resSCa[0][i], np.sqrt(resSC[1][i]**2 + resSCa[1][i]**2))
            #     print "Strat Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resSC[2][i]-resSCa[2][i], np.sqrt(resSC[3][i]**2 + resSCa[3][i]**2))
            #     print "Strat Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resSC[4][i]-resSCa[4][i], np.sqrt(resSC[5][i]**2 + resSCa[5][i]**2))
            #     print "Strat Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resSC[6][i]-resSCa[6][i], np.sqrt(resSC[7][i]**2 + resSCa[7][i]**2))

            #---------------------------------------------------
            #Hemispheric Columns
            #---------------------------------------------------
            print '\nPlot: Hemispheric Differences:\n'
            Lat2 = np.asarray(Lat2)
            
            LatBin = range(-90, 100, 10)
            LatMid = [(i+i+10)*0.5 for i in LatBin[:-1]]

            PCLCBin   = []
            PCFTBin   = []
            PCSCBin   = []

            PCLCBin_e = []
            PCFTBin_e = []
            PCSCBin_e = []

            for i in range(len(LatMid)):
                inds = np.where( (Lat2 >= LatBin[i]) & (Lat2 < LatBin[i+1]))[0]

                if len(inds) > 1:
                    PCLCBin.append(np.mean(resLC[9][inds]))
                    PCFTBin.append(np.mean(resLC2[9][inds]))
                    PCSCBin.append(np.mean(resSC[9][inds]))

                    PCLCBin_e.append(np.sqrt (np.sum(resLC[10][inds]**2))/np.sqrt(len(inds)))
                    PCFTBin_e.append(np.sqrt (np.sum(resLC2[10][inds]**2))/np.sqrt(len(inds)))
                    PCSCBin_e.append(np.sqrt (np.sum(resSC[10][inds]**2))/np.sqrt(len(inds)))

                elif len(inds) == 1:

                    PCLCBin.append(resLC[9][inds])
                    PCFTBin.append(resLC2[9][inds])
                    PCSCBin.append(resSC[9][inds])

                    PCLCBin_e.append(resLC[10][inds])
                    PCFTBin_e.append(resLC2[10][inds])
                    PCSCBin_e.append(resSC[10][inds])

                else:
                    PCLCBin.append(np.nan)
                    PCFTBin.append(np.nan)
                    PCSCBin.append(np.nan)

                    PCLCBin_e.append(np.nan)
                    PCFTBin_e.append(np.nan)
                    PCSCBin_e.append(np.nan)

            PCLCBin = np.asarray(PCLCBin)
            PCFTBin = np.asarray(PCFTBin)
            PCSCBin = np.asarray(PCSCBin)

            PCLCBin_e = np.asarray(PCLCBin_e)
            PCFTBin_e = np.asarray(PCFTBin_e)
            PCSCBin_e = np.asarray(PCSCBin_e)

            fig, ax = plt.subplots()

            #ax.errorbar(Lat, resvmrTC[9], yerr=resvmrTC[10], fmt='o', color='red', ecolor='red', label ='Total')
            ax.errorbar(LatMid, PCLCBin, yerr=PCLCBin_e, fmt='o', color='red', ecolor='red', label ='Low Tropospheric')
            ax.errorbar(LatMid, PCFTBin, yerr=PCFTBin_e, fmt='o', color='blue', ecolor='blue', label= 'Free Tropospheric')
            ax.errorbar(LatMid, PCSCBin, yerr=PCSCBin_e, fmt='o', color='green', ecolor='green', label = 'Stratospheric')
            ax.grid(True)
            ax.set_title('Years submitted', multialignment='center')
            ax.legend(prop={'size':8}, loc=4)

            ax.set_ylabel('Weighted VMR ['+sclfctName+']')
            ax.set_xlabel('Latitude')
            
            fig.subplots_adjust(left=0.1, right=0.95)

            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else:       plt.show(block=False)

            #----------------------------------------
            #Amplitude (Bar plot)
            #----------------------------------------
            fig, ax = plt.subplots()
            ax.barh(ind-0.27, resLC[8], 0.27, align='center', color = 'r', ecolor = 'k', label = 'Low Tropospheric')
            ax.barh(ind, resLC2[8], 0.27, align='center', color = 'b', ecolor = 'k', label = 'Free Tropospheric')
            ax.barh(ind+0.27, resSC[8],0.27, align='center', color = 'g', ecolor = 'k', label = 'Stratospheric')  #yerr=slope_TC*0
            
            ax.xaxis.grid(True, alpha=0.5)
            ax.set_xlabel('Amplitude ['+TCsclfctName+' molecules$\cdot$cm$^{-2}$]')
            ax.set_yticks(ind)
            ax.set_yticklabels(pltID, rotation=0)
            ax.set_title('Amplitude (Partial Columns)', multialignment='center')
            #ax.set_xlabel('Site')
            ax.axvline(0, color='black', lw=1)
            ax.legend(prop={'size':10}, loc = 4)
            ax.set_xlim(0, 2.5)

            plt.gca().invert_yaxis()
            fig.subplots_adjust(left=0.2)

            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else:       plt.show(block=False)

    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #SIMILAR ANALYSIS BUT USING WEIGHTED VMR
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    if pltWvmr:

        #---------------------------------------------------
        # Time Series of Total Weighted VMR  (multiple panels) -- Averaged values
        #---------------------------------------------------
        #print '\nPlot: Averaged WEIGHTED VMR:\n' 
        #resvmrTC = AnalTS(npanels2-1, dates2, totWvmr2, pltID2, Lat2, fits=True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Weighted VMR', unitsStr=sclfctName, ymin=340 ,ymax=600)
        #resvmrTC = np.asarray(resvmrTC)

        #---------------------------------------------------
        # Time Series of Total Weighted VMR Apriori (multiple panels) -- Averaged values
        #---------------------------------------------------
        #print '\nPlot: Averaged WEIGHTED VMR Apriori:\n' 
        #resvmrTCa = AnalTS(npanels2, dates, atotWvmr, pltID, Lat, fits=True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Weighted VMR (Apriori)', unitsStr=sclfctName, ymin=340 ,ymax=600)
        #resvmrTCa = np.asarray(resvmrTCa)

        if pColsFlg:

            #---------------------------------------------------
            # Tropospheric Weighted VMR ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Tropospheric Weighted VMR:\n' 
            resvmrLC = AnalTS(npanels2-1, dates2, WvmrTrop12, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Low Tropospheric Weighted VMR', unitsStr=sclfctName, period=period)#, ymin=340 ,ymax=600)
            resvmrLC = np.asarray(resvmrLC)

            print '\nPlot: Tropospheric Weighted VMR Anomalies:\n' 
            resvmrLCAnom = AnalTSAnom(npanels2-1, dates2, WvmrTrop12, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Low Tropospheric Weighted VMR Anomalies', unitsStr=sclfctName, period=period, qboFlg=False)
            resvmrLCAnom = np.asarray(resvmrLCAnom)

            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

            
            #---------------------------------------------------
            # Free Tropospheric Weighted VMR ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Free Tropospheric Weighted VMR:\n' 
            resvmrLC2 = AnalTS(npanels2-1, dates2, WvmrTrop22, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Free Tropospheric Partial Column', unitsStr=sclfctName, period=period)#, ymin=0.9, ymax=6.0)
            resvmrLC2 = np.asarray(resvmrLC2)

            print '\nPlot: Free Tropospheric Weighted VMR Anomalies:\n' 
            resvmrLCAnom2 = AnalTSAnom(npanels2-1, dates2, WvmrTrop22, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Free Tropospheric Partial Column Anomalies', unitsStr=sclfctName, period=period, qboFlg=False)
            resvmrLCAnom2 = np.asarray(resvmrLCAnom2)

            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
            #---------------------------------------------------
            # Stratospheric Weighted VMR ==> Retrieval
            #---------------------------------------------------
            print '\nPlot: Stratospheric Weighted VMR:\n' 
            resvmrSC = AnalTS(npanels2-1, dates2, WvmrStrat2, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Weighted VMR', unitsStr=sclfctName, period=period)#, ymin=50 ,ymax=450)
            resvmrSC = np.asarray(resvmrSC)

            print '\nPlot: Stratospheric Weighted VMR Anomalies:\n' 
            resvmrSCAnom = AnalTSAnom(npanels2-1, dates2, WvmrStrat2, pltID2, Lat2, fits = True, AvgType=AvgType, smthFlg=smthFlg, pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Weighted VMR Anomalies', unitsStr=sclfctName, period=period, qboFlg=False)
            resvmrSCAnom = np.asarray(resvmrSCAnom)

            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

            #---------------------------------------------------
            # Tropospheric Weighted VMR ==> Apriori
            #---------------------------------------------------
            #print '\nPlot: Tropospheric Weighted VMR Apriori:\n' 
            #resvmrLCa = AnalTS(npanels2, dates, WvmrTropapr, pltID, Lat, fits = True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Tropospheric Weighted VMR (Apriori)', unitsStr=sclfctName, ymin=340 ,ymax=550)
            #resvmrLCa = np.asarray(resvmrLCa)

            #---------------------------------------------------
            # Stratospheric Weighted VMR ==> Apriori
            #---------------------------------------------------
            #print '\nPlot: Stratospheric Weighted VMR Apriori:\n' 
            #resvmrSCa = AnalTS(npanels2, dates, WvmrStratapr, pltID, Lat, fits = True, AvgType='Daily', pltFig=True, saveFlg=saveFlg, pdfsav=pdfsav, ytypeStr='Stratospheric Weighted VMR (Apriori)', unitsStr=sclfctName, ymin=50 ,ymax=450)
            #resvmrSCa = np.asarray(resvmrSCa)

            # #---------------------------------------------------
            # # Plot total Weighted VMR by month
            # #---------------------------------------------------
            # print '\nPlot: Monthly mean Partial Columns:\n' 
            # fig, (ax1, ax2)  = plt.subplots(2, figsize=(10, 10))

            # colormap = plt.cm.gist_ncar
            # ax1.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])
            # ax2.set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, len(pltID))])

            # for i, idhdf in enumerate(pltID):

            #     month    = np.array([d.month for d in dates[idhdf]])
            #     mnthSort = list(set(month))
                
            #     mnthMean = np.zeros(len(mnthSort))
            #     mnthSTD  = np.zeros(len(mnthSort))

            #     mnthMean2 = np.zeros(len(mnthSort))
            #     mnthSTD2  = np.zeros(len(mnthSort))
                
            #     for i,m in enumerate(mnthSort):
            #         inds        = np.where(month == m)[0]
            #         mnthMean[i] = np.mean(WvmrTrop1[idhdf][inds])
            #         mnthSTD[i]  = np.std(WvmrTrop1[idhdf][inds])

            #         mnthMean2[i] = np.mean(WvmrStrat[idhdf][inds])
            #         mnthSTD2[i]  = np.std(WvmrStrat[idhdf][inds])
                    
            #     ax1.errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6, label=idhdf)

            #     ax2.errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6)
            
            # ax1.grid(True,which='both')
            # ax1.set_ylabel('Tropospheric weighted VMR ['+sclfctName+']',multialignment='center')
            # #ax1.set_xlabel('Month')
            # #ax1.set_title('Retrieved Monthly Mean with Standard Deviation')
            # ax1.set_xlim((0,13))
            # ax1.set_xticks(range(1,13))
            # ax1.legend(ncol=4, prop={'size':12}, loc = 'upper center', bbox_to_anchor=[0.5, 1.3])

            # ax2.grid(True,which='both')
            # ax2.set_ylabel('Stratospheric weighted VMR '+sclfctName+']',multialignment='center')
            # ax2.set_xlabel('Month')
            # #ax1.set_title('Retrieved Monthly Mean with Standard Deviation')
            # ax2.set_xlim((0,13))
            # ax2.set_xticks(range(1,13))
           
            # if saveFlg: pdfsav.savefig(fig,dpi=200)
            # else: 
            #     plt.show(block=False)
            #     #user_input = raw_input('Press any key to exit >>> ')
            #     #sys.exit()

            #---------------------------------------------------
            # Plot total Weighted VMR by month
            #---------------------------------------------------
            print '\nPlot: Monthly mean Partial Columns:\n' 
            fig, ax  = plt.subplots(3, 5, figsize=(17, 10), sharey=False, sharex=True)
           
            for i, idhdf in enumerate(pltID2):

                month    = np.array([d.month for d in dates2[idhdf]])
                mnthSort = list(set(month))
                
                mnthMean = np.zeros(len(mnthSort))
                mnthSTD  = np.zeros(len(mnthSort))

                mnthMean2 = np.zeros(len(mnthSort))
                mnthSTD2  = np.zeros(len(mnthSort))

                mnthMean3 = np.zeros(len(mnthSort))
                mnthSTD3  = np.zeros(len(mnthSort))
                
                for j,m in enumerate(mnthSort):
                    inds        = np.where(month == m)[0]
                    mnthMean[j] = np.mean(WvmrTrop12[idhdf][inds])
                    mnthSTD[j]  = np.std(WvmrTrop12[idhdf][inds])

                    mnthMean2[j] = np.mean(WvmrTrop22[idhdf][inds])
                    mnthSTD2[j]  = np.std(WvmrTrop22[idhdf][inds])

                    mnthMean3[j] = np.mean(WvmrStrat2[idhdf][inds])
                    mnthSTD3[j]  = np.std(WvmrStrat2[idhdf][inds])

                Lat2[i] = float(Lat2[i])

                if Lat2[i] >= 50.:     
                    ax[0,0].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,0].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,0].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,0].grid(True,which='both')
                    ax[1,0].grid(True,which='both')
                    ax[2,0].grid(True,which='both')

                    ax[0,0].set_ylabel('Low Tropospheric wVMR ['+sclfctName+']',multialignment='center')
                    ax[1,0].set_ylabel('Free Tropospheric wVMR ['+sclfctName+']',multialignment='center')
                    ax[2,0].set_ylabel('Stratospheric wVMR ['+sclfctName+']',multialignment='center')

                    ax[0,0].set_xlim((0,13))
                    ax[0,0].set_xticks(range(1,13))
                    ax[1,0].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,0].set_xlabel('Month')
                    ax[0, 0].set_title('50$^{\circ}$N - 90$^{\circ}$N', fontsize=14)

                    ax[0, 0].set_ylim(350, 580)
                    ax[1, 0].set_ylim(350, 580)
                    ax[2, 0].set_ylim(100, 450)
                    
                elif (Lat2[i] >= 20.) & (Lat2[i] < 50.):
                    ax[0,1].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,1].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,1].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,1].grid(True,which='both')
                    ax[1,1].grid(True,which='both')
                    ax[2,1].grid(True,which='both')

                    ax[1,1].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,1].set_xlabel('Month')

                    ax[0,1].set_title('20$^{\circ}$N - 50$^{\circ}$N', fontsize=14)

                    ax[0, 1].set_ylim(350, 580)
                    ax[1, 1].set_ylim(350, 580)
                    ax[2, 1].set_ylim(100, 450)
                   
                elif (Lat2[i] >= -20.)  & (Lat2[i] < 20.):
                    ax[0,2].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,2].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,2].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,2].grid(True,which='both')
                    ax[1,2].grid(True,which='both')
                    ax[2,2].grid(True,which='both')

                    ax[1,2].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,2].set_xlabel('Month')
                    ax[0,2].set_title('20$^{\circ}$S - 20$^{\circ}$N', fontsize=14)

                    ax[0, 2].set_ylim(350, 580)
                    ax[1, 2].set_ylim(350, 580)
                    ax[2, 2].set_ylim(100, 450)
                    
                elif (Lat2[i] >= -50.)  & (Lat2[i] < -20.):
                    ax[0,3].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,3].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,3].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,3].grid(True,which='both')
                    ax[1,3].grid(True,which='both')
                    ax[2,3].grid(True,which='both')

                    ax[1,3].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,3].set_xlabel('Month')

                    ax[0,3].set_title('50$^{\circ}$S - 20$^{\circ}$S', fontsize=14)

                    ax[0, 3].set_ylim(350, 580)
                    ax[1, 3].set_ylim(350, 580)
                    ax[2, 3].set_ylim(100, 450)

                    
                elif (Lat2[i] < -50.):
                    ax[0,4].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                    ax[1,4].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                    ax[2,4].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                    ax[0,4].grid(True,which='both')
                    ax[1,4].grid(True,which='both')
                    ax[2,4].grid(True,which='both')

                    ax[1,4].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                    ax[2,4].set_xlabel('Month')
                    ax[0,4].set_title('90$^{\circ}$S - 50$^{\circ}$S', fontsize=14)

                    ax[0, 4].set_ylim(350, 580)
                    ax[1, 4].set_ylim(350, 580)
                    ax[2, 4].set_ylim(100, 450)

            fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.98, wspace=0.15, hspace=0.15)
           
            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else: 
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

            #---------------------------------------------------
            # Bar plot: Three different periods ==> Retrieval
            #---------------------------------------------------
            #hbarplt3(resvmrLC, resvmrLC2, resvmrSC, pltID2, b1_label='Low Tropospheric', b2_label='Free Tropospheric', b3_label='Stratospheric', subtitle='Retrieval', saveFlg=saveFlg, pdfsav=pdfsav)
            hbarplt3(resvmrLCAnom, resvmrLCAnom2, resvmrSCAnom, pltID2, b1_label='Low Tropospheric', b2_label='Free Tropospheric', b3_label='Stratospheric', subtitle='Retrieval - Anomalies', saveFlg=saveFlg, pdfsav=pdfsav)


            #---------------------------------------------------
            # Bar plot: Three different periods ==> Apriori
            #---------------------------------------------------
            #hbarplt3(resvmrTCa, resvmrLCa, resvmrSCa, pltID, b1_label='Tot Column', b2_label='Trop Column', b3_label='Strat Column', subtitle='A priori', saveFlg=saveFlg, pdfsav=pdfsav)

            #---------------------------------------------------
            # Bar plot: Three different periods ==> Retrieval - Apriori
            #---------------------------------------------------
            # pltID        = np.asarray(pltID)

            ind = np.arange(len(resvmrLC[0]))
            
            # fig, (ax, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=(15, 8.5))
            
            # ax.barh(ind-0.27, resvmrTC[0]-resvmrTCa[0], 0.27, xerr=np.sqrt(resvmrTC[1]**2 + resvmrTCa[1]**2) , align='center', color = 'r', ecolor = 'k', label = 'Tot Column')
            # ax.barh(ind, resvmrLC[0]-resvmrLCa[0], 0.27, xerr=np.sqrt(resvmrLC[1]**2 + resvmrLCa[1]**2), align='center', color = 'b', ecolor = 'k', label = 'Trop Column')
            # ax.barh(ind+0.27, resvmrSC[0]-resvmrSCa[0], 0.27, xerr=np.sqrt(resvmrSC[1]**2 + resvmrSCa[1]**2), align='center', color = 'g', ecolor = 'k', label = 'Strat Column')  #yerr=slope_TC*0
            # ax.xaxis.grid(True)
            # ax.yaxis.grid(True)
            # ax.set_xlabel('Rate of change (%/y)')
            # ax.set_yticks(ind)
            # ax.set_yticklabels(np.transpose(pltID), rotation=0)
            # ax.set_title('Retrieval - Apriori\n(years submitted)', multialignment='center')
            # #ax.set_xlabel('Site')
            # ax.set_xlim(-2.0, 2.0)
            # ax.axvline(0, color='black', lw=1)
            # ax.legend(prop={'size':8}, loc = 1)
            # #ax.invert_yaxis()

            # ax2.barh(ind-0.27, resvmrTC[2]-resvmrTCa[2], 0.27, xerr=np.sqrt(resvmrTC[3]**2 + resvmrTCa[3]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax2.barh(ind, resvmrLC[2]-resvmrLCa[2], 0.27, xerr=np.sqrt(resvmrLC[3]**2 + resvmrLCa[3]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax2.barh(ind+0.27, resvmrSC[2]-resvmrSCa[2], 0.27, xerr=np.sqrt(resvmrSC[3]**2 + resvmrSCa[3]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax2.xaxis.grid(True)
            # ax2.yaxis.grid(True)
            # ax2.set_xlabel('Rate of change (%/y))')
            # ax2.set_title('Retrieval - Apriori\n(1995 - 2002)', multialignment='center')
            # ax2.set_yticks(ind)
            # #ax2.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax2.set_xlim(-2.0, 2.0)
            # ax2.axvline(0, color='black', lw=1)

            # #ax2.invert_yaxis()
            # #ax2.yticks([])

            # ax3.barh(ind-0.27, resvmrTC[4]-resvmrTCa[4], 0.27, xerr=np.sqrt(resvmrTC[5]**2 + resvmrTCa[5]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax3.barh(ind, resvmrLC[4]-resvmrLCa[4], 0.27, xerr=np.sqrt(resvmrLC[5]**2 + resvmrLCa[5]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax3.barh(ind+0.27, resvmrSC[4]-resvmrSCa[4], 0.27, xerr=np.sqrt(resvmrSC[5]**2 + resvmrSCa[5]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax3.xaxis.grid(True)
            # ax2.yaxis.grid(True)
            # ax3.set_xlabel('Rate of change (%/y))')
            # ax3.set_title('Retrieval - Apriori\n(2002 - 2008)', multialignment='center')
            # ax3.set_yticks(ind)
            # #ax3.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax3.set_xlim(-2.0, 2.0)
            # ax3.axvline(0, color='black', lw=1)

            # #ax3.invert_yaxis()

            # ax4.barh(ind-0.27, resvmrTC[6]-resvmrTCa[6], 0.27, xerr=np.sqrt(resvmrTC[7]**2 + resvmrTCa[7]**2), align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax4.barh(ind, resvmrLC[6]-resvmrLCa[6], 0.27, xerr=np.sqrt(resvmrLC[7]**2 + resvmrLCa[7]**2), align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax4.barh(ind+0.27, resvmrSC[6]-resvmrSCa[6], 0.27, xerr=np.sqrt(resvmrSC[7]**2 + resvmrSCa[7]**2), align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            # ax4.xaxis.grid(True)
            # ax4.yaxis.grid(True)
            # ax4.set_xlabel('Rate of change (%/y))')
            # ax4.set_title('Retrieval - Apriori\n(2008 - 2016)', multialignment='center')
            # ax4.set_yticks(ind)
            # #ax4.set_yticklabels(np.transpose(pltID), rotation=0)
            # #ax.set_xlabel('Site')
            # ax4.set_xlim(-2.0, 2.0)
            # ax4.axvline(0, color='black', lw=1)

            # #ax4.invert_yaxis()

            # plt.gca().invert_yaxis()
            # #fig.tight_layout()
            # fig.subplots_adjust(left=0.1, right=0.95)

            # if saveFlg: pdfsav.savefig(fig,dpi=200)
            # else:       plt.show(block=False)

            # for i, idhdf in enumerate(pltID):

            #     print '\nRate of change (Retrieval) - VMR: {}\n'.format(idhdf)

            #     print "Total Columns - all years   = {:.2f} +/- {:.2f}% ".format(resvmrTC[0][i], resvmrTC[1][i])
            #     print "Total Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resvmrTC[2][i], resvmrTC[3][i])
            #     print "Total Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resvmrTC[4][i], resvmrTC[5][i])
            #     print "Total Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resvmrTC[6][i], resvmrTC[7][i])

            #     print "Trop Columns - all years    = {:.2f} +/- {:.2f}% ".format(resvmrLC[0][i], resvmrLC[1][i])
            #     print "Trop Columns - 1995 - 2002  = {:.2f} +/- {:.2f}% ".format(resvmrLC[2][i], resvmrLC[3][i])
            #     print "Trop Columns - 2002 - 2008  = {:.2f} +/- {:.2f}% ".format(resvmrLC[4][i], resvmrLC[5][i])
            #     print "Trop Columns - 2008 - 2016  = {:.2f} +/- {:.2f}% ".format(resvmrLC[6][i], resvmrLC[7][i])

            #     print "Strat Columns - all years   = {:.2f} +/- {:.2f}% ".format(resvmrSC[0][i], resvmrSC[1][i])
            #     print "Strat Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resvmrSC[2][i], resvmrSC[3][i])
            #     print "Strat Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resvmrSC[4][i], resvmrSC[5][i])
            #     print "Strat Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resvmrSC[6][i], resvmrSC[7][i])

            # for i, idhdf in enumerate(pltID):

            #     print '\nRate of change (Retrieval  - Apriori) - VMR: {}\n'.format(idhdf)

            #     print "Total Columns - all years   = {:.2f} +/- {:.2f}% ".format(resvmrTC[0][i]-resvmrTCa[0][i], np.sqrt(resvmrTC[1][i]**2 + resvmrTCa[1][i]**2))
            #     print "Total Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resvmrTC[2][i]-resvmrTCa[2][i], np.sqrt(resvmrTC[3][i]**2 + resvmrTCa[3][i]**2))
            #     print "Total Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resvmrTC[4][i]-resvmrTCa[4][i], np.sqrt(resvmrTC[5][i]**2 + resvmrTCa[5][i]**2))
            #     print "Total Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resvmrTC[6][i]-resvmrTCa[6][i], np.sqrt(resvmrTC[7][i]**2 + resvmrTCa[7][i]**2))

            #     print "Trop Columns - all years    = {:.2f} +/- {:.2f}% ".format(resvmrLC[0][i]-resvmrLCa[0][i], np.sqrt(resvmrLC[1][i]**2 + resvmrLCa[1][i]**2))
            #     print "Trop Columns - 1995 - 2002  = {:.2f} +/- {:.2f}% ".format(resvmrLC[2][i]-resvmrLCa[2][i], np.sqrt(resvmrLC[3][i]**2 + resvmrLCa[3][i]**2))
            #     print "Trop Columns - 2002 - 2008  = {:.2f} +/- {:.2f}% ".format(resvmrLC[4][i]-resvmrLCa[4][i], np.sqrt(resvmrLC[5][i]**2 + resvmrLCa[5][i]**2))
            #     print "Trop Columns - 2008 - 2016  = {:.2f} +/- {:.2f}% ".format(resvmrLC[6][i]-resvmrLCa[6][i], np.sqrt(resvmrLC[7][i]**2 + resvmrLCa[7][i]**2))

            #     print "Strat Columns - all years   = {:.2f} +/- {:.2f}% ".format(resvmrSC[0][i]-resvmrSCa[0][i], np.sqrt(resvmrSC[1][i]**2 + resvmrSCa[1][i]**2))
            #     print "Strat Columns - 1995 - 2002 = {:.2f} +/- {:.2f}% ".format(resvmrSC[2][i]-resvmrSCa[2][i], np.sqrt(resvmrSC[3][i]**2 + resvmrSCa[3][i]**2))
            #     print "Strat Columns - 2002 - 2008 = {:.2f} +/- {:.2f}% ".format(resvmrSC[4][i]-resvmrSCa[4][i], np.sqrt(resvmrSC[5][i]**2 + resvmrSCa[5][i]**2))
            #     print "Strat Columns - 2008 - 2016 = {:.2f} +/- {:.2f}% ".format(resvmrSC[6][i]-resvmrSCa[6][i], np.sqrt(resvmrSC[7][i]**2 + resvmrSCa[7][i]**2))



            #---------------------------------------------------
            #Hemispheric Columns
            #---------------------------------------------------
            print '\nPlot: Hemispheric Differences:\n'
            Lat2 = np.asarray(Lat2)

            fig, (ax, ax2, ax3, ax4) = plt.subplots(4, 1, sharey=True, sharex=True, figsize=(7, 12))

            #ax.errorbar(Lat, resvmrTC[9], yerr=resvmrTC[10], fmt='o', color='red', ecolor='red', label ='Total')
            ax.errorbar(Lat2, resvmrLC[9], yerr=resvmrLC[10], fmt='o', color='red', ecolor='red', label ='Low Tropospheric')
            ax.errorbar(Lat2, resvmrLC2[9], yerr=resvmrLC2[10], fmt='o', color='blue', ecolor='blue', label= 'Free Tropospheric')
            ax.errorbar(Lat2, resvmrSC[9], yerr=resvmrSC[10], fmt='o', color='green', ecolor='green', label = 'Stratospheric')
            ax.grid(True)
            ax.set_title('Years submitted', multialignment='center')
            ax.legend(prop={'size':8}, loc=4)

            ax2.errorbar(Lat2, resvmrLC[11], yerr=resvmrLC[12], fmt='o', color='red', ecolor='red')
            ax2.errorbar(Lat2, resvmrLC2[11], yerr=resvmrLC2[12], fmt='o', color='blue', ecolor='blue')
            ax2.errorbar(Lat2, resvmrSC[11], yerr=resvmrSC[12], fmt='o', color='green', ecolor='green')
            ax2.grid(True)
            ax2.set_title('1998 - 2002', multialignment='center')

            ax3.errorbar(Lat2, resvmrLC[13], yerr=resvmrLC[14], fmt='o', color='red', ecolor='red')
            ax3.errorbar(Lat2, resvmrLC2[13], yerr=resvmrLC2[14], fmt='o', color='blue', ecolor='blue')
            ax3.errorbar(Lat2, resvmrSC[13], yerr=resvmrSC[14], fmt='o', color='green', ecolor='green')
            ax3.grid(True)
            ax3.set_title('2002 - 2008', multialignment='center')

            ax4.errorbar(Lat2, resvmrLC[15], yerr=resvmrLC[16], fmt='o', color='red', ecolor='red')
            ax4.errorbar(Lat2, resvmrLC2[15], yerr=resvmrLC2[16], fmt='o', color='blue', ecolor='blue')
            ax4.errorbar(Lat2, resvmrSC[15], yerr=resvmrSC[16], fmt='o', color='green', ecolor='green')
            ax4.grid(True)
            ax4.set_title('2008 - 2016', multialignment='center')
            ax4.set_xlabel('Latitude')
            ax4.set_xlim(-90, 90)

            fig.text(0.015, 0.5, 'Weighted VMR ['+sclfctName+']', fontsize=16, va='center', rotation='vertical')
            
            fig.subplots_adjust(left=0.1, right=0.95)

            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else:       plt.show(block=False)

            #---------------------------------------------------
            #Hemispheric Columns
            #---------------------------------------------------
            print '\nPlot: Hemispheric Differences:\n'
            Lat2 = np.asarray(Lat2)
            
            LatBin = range(-90, 100, 10)
            LatMid = [(i+i+10)*0.5 for i in LatBin[:-1]]

            vmrLCBin   = []
            vmrFTBin   = []
            vmrSCBin   = []

            vmrLCBin_e = []
            vmrFTBin_e = []
            vmrSCBin_e = []

            for i in range(len(LatMid)):
                inds = np.where( (Lat2 >= LatBin[i]) & (Lat2 < LatBin[i+1]))[0]

                if len(inds) > 1:
                    vmrLCBin.append(np.mean(resvmrLC[9][inds]))
                    vmrFTBin.append(np.mean(resvmrLC2[9][inds]))
                    vmrSCBin.append(np.mean(resvmrSC[9][inds]))

                    vmrLCBin_e.append(np.sqrt (np.sum(resvmrLC[10][inds]**2))/np.sqrt(len(inds)))
                    vmrFTBin_e.append(np.sqrt (np.sum(resvmrLC2[10][inds]**2))/np.sqrt(len(inds)))
                    vmrSCBin_e.append(np.sqrt (np.sum(resvmrSC[10][inds]**2))/np.sqrt(len(inds)))

                elif len(inds) == 1:

                    vmrLCBin.append(resvmrLC[9][inds])
                    vmrFTBin.append(resvmrLC2[9][inds])
                    vmrSCBin.append(resvmrSC[9][inds])

                    vmrLCBin_e.append(resvmrLC[10][inds])
                    vmrFTBin_e.append(resvmrLC2[10][inds])
                    vmrSCBin_e.append(resvmrSC[10][inds])

                else:
                    vmrLCBin.append(np.nan)
                    vmrFTBin.append(np.nan)
                    vmrSCBin.append(np.nan)

                    vmrLCBin_e.append(np.nan)
                    vmrFTBin_e.append(np.nan)
                    vmrSCBin_e.append(np.nan)

            vmrLCBin = np.asarray(vmrLCBin)
            vmrFTBin = np.asarray(vmrFTBin)
            vmrSCBin = np.asarray(vmrSCBin)

            vmrLCBin_e = np.asarray(vmrLCBin_e)
            vmrFTBin_e = np.asarray(vmrFTBin_e)
            vmrSCBin_e = np.asarray(vmrSCBin_e)

            fig, ax = plt.subplots()

            #ax.errorbar(Lat, resvmrTC[9], yerr=resvmrTC[10], fmt='o', color='red', ecolor='red', label ='Total')
            ax.errorbar(LatMid, vmrLCBin, yerr=vmrLCBin_e, fmt='o', color='red', ecolor='red', label ='Low Tropospheric')
            ax.errorbar(LatMid, vmrFTBin, yerr=vmrFTBin_e, fmt='o', color='blue', ecolor='blue', label= 'Free Tropospheric')
            ax.errorbar(LatMid, vmrSCBin, yerr=vmrSCBin_e, fmt='o', color='green', ecolor='green', label = 'Stratospheric')
            ax.grid(True)
            ax.set_title('Years submitted', multialignment='center')
            ax.legend(prop={'size':8}, loc=4)

            ax.set_ylabel('Weighted VMR ['+sclfctName+']')
            ax.set_xlabel('Latitude')
            
            fig.subplots_adjust(left=0.1, right=0.95)

            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else:       plt.show(block=False)



            #----------------------------------------
            #Amplitude (Bar plot)
            #----------------------------------------
            # fig, ax = plt.subplots()
            # ax.barh(ind-0.27, resvmrTC[8], 0.27, align='center', color = 'r', ecolor = 'k', label = 'Total Column')
            # ax.barh(ind, resvmrLC[8], 0.27, align='center', color = 'b', ecolor = 'k', label = 'Tropospheric Column')
            # ax.barh(ind+0.27, resvmrSC[8],0.27, align='center', color = 'g', ecolor = 'k', label = 'Stratospheric Column')  #yerr=slope_TC*0
            
            # ax.xaxis.grid(True, alpha=0.5)
            # ax.set_xlabel('Amplitude - Weighted VMR ['+sclfctName+']')
            # ax.set_yticks(ind)
            # ax.set_yticklabels(pltID, rotation=0)
            # ax.set_title('Amplitude', multialignment='center')
            # #ax.set_xlabel('Site')
            # ax.axvline(0, color='black', lw=1)
            # ax.legend(prop={'size':10}, loc = 4)
            # #ax.set_xlim(0, 2.5)

            # plt.gca().invert_yaxis()
            # fig.subplots_adjust(left=0.2)

            # if saveFlg: pdfsav.savefig(fig,dpi=200)
            # else:       plt.show(block=False)

            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        
    
    #----------------------------------------
    #Annual Rate of Change of TH (Bar plot)
    #----------------------------------------
    # fig, ax = plt.subplots()

    # ax.barh(ind, resTH[0], 0.27, xerr=resTH[1], align='center', color = 'r', ecolor = 'k')
   
    # ax.xaxis.grid(True, alpha=0.5)
    # ax.set_xlabel('Annual rate of change (%)')
    # ax.set_yticks(ind)
    # ax.set_yticklabels(pltID, rotation=0)
    # ax.set_title('Annual Rate of change (Tropospheric Height)', multialignment='center')
    # #ax.set_xlabel('Site')
    # ax.axvline(0, color='black', lw=1)

    # plt.gca().invert_yaxis()

    # fig.subplots_adjust(left=0.2)

    # if saveFlg: pdfsav.savefig(fig,dpi=200)
    # else:       plt.show(block=False)

        
    # #----------------------------------------
    # #Amplitude TH (Bar plot)
    # #----------------------------------------
    # fig, ax = plt.subplots()

    # ax.barh(ind, resTH[8], 0.27, align='center', color = 'r', ecolor = 'k')
    
    # ax.xaxis.grid(True, alpha=0.5)
    # ax.set_xlabel('Amplitude [km]')
    # ax.set_yticks(ind)
    # ax.set_yticklabels(pltID, rotation=0)
    # ax.set_title('Amplitude (Tropospheric Height)', multialignment='center')
    # #ax.set_xlabel('Site')
    # ax.axvline(0, color='black', lw=1)

    # plt.gca().invert_yaxis()
    # fig.subplots_adjust(left=0.2)

    # if saveFlg: pdfsav.savefig(fig,dpi=200)
    # else:       plt.show(block=False)

    #---------------------------------------------------
    # Averaging Kernel Smoothing Function (row of avk)
    #---------------------------------------------------
    npanels2  = int(math.ceil(npanels/4.0))

    fig = plt.figure(figsize=(15,15))

    outer_grid = gridspec.GridSpec(npanels2, 4, wspace=0.17, hspace=0.135)

    for i, idhdf in enumerate(pltID):

        if pltWvmr: pltAKAv   = np.nanmean(avkVMR[idhdf], axis=0)
        else: pltAKAv   = np.nanmean(avkSCF[idhdf], axis=0)
        
        avkSCFAv   = np.nanmean(avkSCF[idhdf], axis=0)

        #---------
        #Total Column AK
        #---------
        avkTCAv    = np.nanmean(avkTC[idhdf], axis=0)
        
        dtpMean    = np.nanmean(dtp[idhdf])
        dtpStd     = np.nanstd(dtp[idhdf])
        Pcol       = [dtpMean, 120.]

        #ax = plt.Subplot(fig, outer_grid[i])
        #gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
        gs        = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer_grid[i], width_ratios=[3,1,1])
        axa       = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        axc       = plt.subplot(gs[2])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt[idhdf]), vmax=np.max(alt[idhdf]))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt[idhdf])
        axa.set_color_cycle([scalarMap.to_rgba(x) for x in alt[idhdf]])
        
        for j in range(len(alt[idhdf])):
            axa.plot(pltAKAv[j,:],alt[idhdf])

        axa.axhline(y=dtpMean, linewidth=1.5, color='r')  
        axa.axhline(y=dtpMean - dtpStd, linewidth=1.5, color='r', alpha=0.5)  
        axa.axhline(y=dtpMean + dtpStd, linewidth=1.5, color='r', alpha=0.5)  
            
        #axa.set_ylabel('Altitude [km]')
        #axa.set_xlabel('Averaging Kernels')
        axa.grid(True, alpha=0.5)
        axa.set_ylim(0, 40)
        axa.set_xlim(-0.15, 0.25)
        axa.xaxis.set_ticks(np.arange(-0.15, 0.25, 0.1))
        #axa.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
        #axa.xticks(np.arange(-0.15, 0.25, 0.1))

        axa.tick_params(axis = 'both', which = 'major', labelsize = 8)
        axa.tick_params(axis = 'both', which = 'minor', labelsize = 0)

        if pltWvmr:
            if i == npanels-1: axa.set_xlabel('VMR AK')
            if i == npanels-2: axa.set_xlabel('VMR AK')
            if i == npanels-3: axa.set_xlabel('VMR AK')
        else:
            if i == npanels-1: axa.set_xlabel('AK')
            if i == npanels-2: axa.set_xlabel('AK')
            if i == npanels-3: axa.set_xlabel('AK')

        axa.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.875), xycoords='axes fraction', fontsize=10, ha='left')
        #all_axes[-2].set_xlabel('Averaging Kernels')
        #axa.set_xticks([])
        #axa.set_yticks([])

        #cbaxes = fig.add_axes([0.45, 0.47, 0.02, 0.35])
        #cbar = fig.colorbar(scalarMap,orientation='vertical', cax = cbaxes)
        #cbar.set_label('Altitude [km]')
        #ax.set_title('Averaging Kernels Scale Factor - ' +str(pltID[i]))
        #plt.suptitle('Averaging Kernels Scale Factor - ' +str(pltID[i]), fontsize=16)
        
        #axb.plot(np.sum(avkSCFAv,axis=0),alt[idhdf],color='k')
        axb.plot(avkTCAv,alt[idhdf],color='k')
        
        axb.grid(True)
        axb.axhline(y=dtpMean, linewidth=1.5, color='r')  
        axb.axhline(y=dtpMean - dtpStd, linewidth=1.5, color='r', alpha=0.5)  
        axb.axhline(y=dtpMean + dtpStd, linewidth=1.5, color='r', alpha=0.5)
        axb.set_ylim(0, 40)
        axb.set_xlim(0, 2)
        axb.xaxis.set_ticks(np.arange(0,2,1))
        #axb.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
        #axb.xticks(np.arange(0, 2, 1.0))

        axb.tick_params(axis = 'both', which = 'major', labelsize = 8)
        axb.tick_params(axis = 'both', which = 'minor', labelsize = 0)
        if i == npanels-1: axb.set_xlabel('TC AK')
        if i == npanels-2: axb.set_xlabel('TC AK')
        if i == npanels-3: axb.set_xlabel('TC AK')

        if alt[idhdf][0] > alt[idhdf][-1]:
            dofs_cs = np.cumsum(np.diag(avkSCFAv)[::-1])[::-1]
        else:
            dofs_cs = np.cumsum(np.diag(avkSCFAv))

        
        axc.plot(dofs_cs,alt[idhdf],color='k',label='Cumulative Sum of DOFS (starting at surface)')
        

        try:
            xval = range(0,int(np.ceil(max(dofs_cs)))+3)

        except Exception as errmsg:
            print '\nError: ', errmsg
            xval = range(0, 2)

        ind1         = mf.nearestind(Pcol[0], alt[idhdf])
        ind2         = mf.nearestind(Pcol[1], alt[idhdf])

        axc.fill_between(xval,alt[idhdf][ind1],alt[idhdf][ind2],alpha=0.5,color='0.75')  
        axc.axhline(alt[idhdf][ind2],color='k',linestyle='--')
        dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
        #axc.text(0.15,(alt[idhdf][ind1]+alt[idhdf][ind2])/2.0, 
        #         'DOFs for layer {0:.2f}-{1:.2f}[km] = {2:.3f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol),
        #         fontsize=9)
        #axc.set_title('DOFs Profile - ' +str(pltID[i]))
        #axc.set_ylabel('Altitude [km]')
        
        #axc.set_xlabel('Cumulative\nSum of DOFS')  
        #axc.set_title('sDOF = {:.2f}'.format(dofsPcol), fontsize=9)    
        axc.set_ylim((0,40))
        axc.grid(True,which='both')
        axc.set_ylim(0, 40)
        axc.set_xlim(0, 3)
        axc.xaxis.set_ticks(np.arange(0,3,1))
        #axc.xaxis.set_major_formatter(ticker.FormatStrFormatter('%1f'))
        #axb.xticks(np.arange(0, 3, 1.0))
        axc.tick_params(axis='x',which='both',labelsize=8)
        axc.tick_params(axis = 'both', which = 'major', labelsize = 8)
        axc.tick_params(axis = 'both', which = 'minor', labelsize = 0)
        if i == npanels-1: axc.set_xlabel('Cumulative\nSum of DOFS')
        if i == npanels-2: axc.set_xlabel('Cumulative\nSum of DOFS')
        if i == npanels-3: axc.set_xlabel('Cumulative\nSum of DOFS')

        print idhdf
        print 'DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol)
        print 'DOFs all = {0:.2f}'.format(np.trace(avkSCFAv))
      
        fig.add_subplot(axa)
        fig.add_subplot(axb)
        fig.add_subplot(axc)  

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

    #if (npanels % 2 == 1): #even

    #    all_axes[-2].spines['bottom'].set_visible(True)
    #    plt.setp(all_axes[-2].get_xticklabels(), visible=True)
    #    all_axes[-2].set_zorder(1)

    

    fig.text(0.03, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
    plt.suptitle('{} Averaging Kernels'.format(gasName.upper()), fontsize=16  )
    #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
    fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

    
    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else: 
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()
        

    if saveFlg: pdfsav.close()
    else:
        plt.show(block=False)
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

def AnalTSAnom(npanels, xDates, yData, pltID, Lat, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1, qboFlg=False):
    
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

    xmin      = dt.date(1993, 1, 1)
    xmax      = dt.date(2016, 12, 31)


    if pltFig:    
        
        fig  = plt.figure(figsize=(18,13))
        fig2 = plt.figure(figsize=(18,13))  
        
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.2, hspace=0.075)

    for i, idhdf in enumerate(pltID):

        print idhdf

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
            ax  = plt.Subplot(fig, outer_grid[i])
            fig.add_subplot(ax)
            ax.plot(Dates, Anomaly,'k.',markersize=0)
            ax.scatter(Dates, Anomaly, s=20, facecolor='lightgray', edgecolor='k',alpha=0.85)
            #print len(Dates)
            #print len(f_fourier(dateYeardateYearFracFrac))
            if qboFlg: ax.plot(dates2,f_fourier(dateYearFrac2)-300.,label='Fitted Anual Trend',linewidth=2.0)


            ax2 = plt.Subplot(fig2, outer_grid[i])
            ax2.plot(Dates, Anomaly2,'k.',markersize=0)
            ax2.scatter(Dates, Anomaly2, s=20, facecolor='lightgray', edgecolor='k',alpha=0.85)

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
            yoi    = [1996, 2016]
            indx1  = np.where( (yyyy >= yoi[0]) &  (yyyy <= yoi[1]))[0]

            weights      = np.ones_like(dateYearFrac)

            res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
            f_drift, f_fourier, f_driftfourier,  res_std, A, df_drift = res[3:9]

            res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
            perc, intercept_b, slope_b, pfourier_b = res_b

            print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yoi[0], yoi[-1])

            if (idhdf == 'Zugspitze') or (idhdf == 'Jungfraujoch') or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS'):

                #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                

                slope.append(res[1]/np.mean(AvgData[indx1])*100.0)
                slope_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                #print res[1]/np.mean(AvgData)*100.0

                Amp   = np.sum(res[2]**2)
                Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                amp.append(Amp)

                avgD.append(np.mean(AvgData[indx1]))
                stdD.append(np.std(AvgData[indx1]))

                ax2.plot(Dates,f_drift(dateYearFrac),label='Fitted Anual Trend',linewidth=2.0, color='k')

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
            if (idhdf != 'Eureka') & (idhdf != 'St Petersburg') & (idhdf != 'Boulder') & (idhdf != 'Maido') &  (idhdf != 'StD-Maido') & (idhdf != 'Bremen') & (idhdf != 'Paris') & (idhdf != 'Altzomoni'):

                res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly2, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]

        #---------------------------------------------------
        #---------------------------------------------------
        roc  = df_drift(dateYearFrac2)/np.mean(AvgData)*100.0
        roc  = np.asarray(roc)

        zero_crossings = np.where(np.diff(np.sign(roc)))[0]
    
        print idhdf
        print 'Number of crossings = {}'.format(len(dates2[zero_crossings]))
        print 'Dates of crossings  = {}'.format(dates2[zero_crossings])

        # if len(zero_crossings) <= 3:

        #     for iz, z in zero_crossings:
        #         if iz == 0: print np.mean(roc[:z]) 
        #         if iz == 1: print np.mean(roc[:z]) 
        
            
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
            #     ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  

            Lat[i] = float(Lat[i])

            if Lat[i] >= 50.:     
                ax.set_axis_bgcolor('lightcyan')
                ax2.set_axis_bgcolor('lightcyan')
                
            elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                ax.set_axis_bgcolor('lightgreen')
                ax2.set_axis_bgcolor('lightgreen')
               
            elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                ax.set_axis_bgcolor('mistyrose')
                ax2.set_axis_bgcolor('mistyrose')
                
            elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                ax.set_axis_bgcolor('cornsilk')
                ax2.set_axis_bgcolor('cornsilk')
                
            elif (Lat[i] < -50.):
                ax.set_axis_bgcolor('lightgrey')
                ax2.set_axis_bgcolor('lightgrey')
            
            else:
                ax.set_axis_bgcolor('lightgrey')
                ax2.set_axis_bgcolor('lightgrey')
 
        if fits:

            #if pltFig:
                #ax.plot(dates2,f_fourier(dateYearFrac2),label='Fitted Anual Trend', linewidth=2.0)
                #ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #ax2.plot(dates2,f_drift(dateYearFrac2),label='Fitted Anual Trend',linewidth=2.0)

                #ax2.plot(dates2, roc, linewidth=2.0)

            indsmin = np.where(f_drift(dateYearFrac2) == np.min(f_drift(dateYearFrac2)))[0]
            indsmax = np.where(f_drift(dateYearFrac2) == np.max(f_drift(dateYearFrac2)))[0]

            #---------------------------------------------------


            if (idhdf == 'Paris') or (idhdf == 'Altzomoni'):


                slope_p1.append(float('nan'))
                slope_p1_e.append(float('nan'))

                slope_p2.append(float('nan'))
                slope_p2_e.append(float('nan'))

                slope_p3.append(float('nan'))
                slope_p3_e.append(float('nan'))

                avgD_p1.append(float('nan'))
                stdD_p1.append(float('nan'))

                avgD_p2.append(float('nan'))
                stdD_p2.append(float('nan'))

                avgD_p3.append(float('nan'))
                stdD_p3.append(float('nan'))



            elif (idhdf == 'Eureka') or (idhdf == 'Kiruna') or (idhdf == 'St Petersburg') or (idhdf == 'Bremen') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Tsukuba')  or (idhdf == 'Izana') or (idhdf == 'Paramaribo'):

                yoi  = [[2009, 2016]]

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

                    ax2.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color='red')


            elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund'):

                yoi  = [[2002, 2008], [2009, 2016]]

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

                    ax2.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])
              
                    #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

            elif (idhdf == 'Rikubetsu') or (idhdf == 'Mauna Loa'):

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

                    ax2.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])



            elif (idhdf == 'Jungfraujoch')  or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze') :

                yoi  = [[1996, 2002], [2002, 2008], [2009, 2016]]

                clr_yoi = ['blue', 'green', 'red']
            
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

                    if ii == 2:
                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))

                    ax2.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])

                #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

        if pltFig:
            yearsLc      = YearLocator()
            months       = MonthLocator()
            DateFmt      = DateFormatter('%Y')

            ax.set_xlim(xmin, xmax)
            ax.grid(True, color='gray', alpha=0.5)
            ax.tick_params(which='both',labelsize=10)
            ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

            ax.xaxis.set_major_locator(yearsLc)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')

            if (ymin and ymax):
                ax.set_ylim(ymin, ymax)

            if qboFlg:
                axx.set_xlim(xmin, xmax) 
                axx.tick_params(which='both',labelsize=10)
                axx.tick_params(axis = 'both', which = 'minor', labelsize = 0)
                axx.set_ylim(-400, 300) 






            #axx.set_xlim(xmin, xmax)

            #fig.add_subplot(axx, ax)
            #fig.add_subplot(axx)

            #---------------------------------------------------
            #
            #---------------------------------------------------
            ax2.set_xlim(xmin, xmax)
            #ax2.set_ylim(-7, 7)
            #ax2.axhline(y=0, linestyle='--', color='k')
            ax2.grid(True, color='gray', alpha=0.5)
            ax2.tick_params(which='both',labelsize=10)
            ax2.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

            ax2.xaxis.set_major_locator(yearsLc)
            ax2.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            ax2.xaxis.set_major_formatter(DateFmt)
            ax2.xaxis.set_tick_params(which='minor',labelbottom='off')

            if (ymin and ymax):
                ax2.set_ylim(ymin, ymax)
    
            fig2.add_subplot(ax2)

            #---------------------------------------------------

    if pltFig:
    
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (npanels % 2 == 1): #even
            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        all_axes[-1].set_xlabel('Year')
        all_axes[-2].set_xlabel('Year')
        all_axes[-3].set_xlabel('Year')
        
        fig.autofmt_xdate()
        fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        fig.suptitle(ytypeStr, fontsize=16  )
        
        if saveFlg: pdfsav.savefig(fig,dpi=200)
        else: 
            plt.show(block=False)

    #---------------------------------------------------
    #
    #---------------------------------------------------
    if pltFig:
    
        all_axes = fig2.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (npanels % 2 == 1): #even
            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        all_axes[-1].set_xlabel('Year')
        all_axes[-2].set_xlabel('Year')
        all_axes[-3].set_xlabel('Year')
        
        fig2.autofmt_xdate()
        fig2.text(0.03, 0.5, 'Rate of change [%/y]', fontsize=16, va='center', rotation='vertical')
        fig2.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        fig2.suptitle('Rate of change of VMR Anomalies', fontsize=16  )
        
        if saveFlg: pdfsav.savefig(fig2,dpi=200)
        else: 
            plt.show(block=False)


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)


def AnalTS(npanels, xDates, yData, pltID, Lat, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1):
    
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

    xmin      = dt.date(1993, 1, 1)
    xmax      = dt.date(2016, 12, 31)

    if pltFig:    
        fig = plt.figure(figsize=(18,13))  

        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.1, hspace=0.075)

    for i, idhdf in enumerate(pltID):

        #----------------------------
        #CONCATENATE st denis and maido
        #----------------------------
        # if idhdf == 'St Denis':
        #     yData[idhdf]    = np.asarray(yData[idhdf])
        #     yData['Maido']  = np.asarray(yData['Maido'])

        #     xDates[idhdf]    = np.asarray(xDates[idhdf])
        #     xDates['Maido']  = np.asarray(xDates['Maido'])

        #     #if (len(yData[idhdf]) != len(xDates[idhdf])): 
        #     yData[idhdf]  = np.concatenate( (yData[idhdf], yData['Maido']))


        #     if (len(xDates[idhdf]) != len(yData[idhdf])): 
        #         xDates[idhdf] = np.concatenate( (xDates[idhdf], xDates['Maido']))
        
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
            ax.scatter(Dates,AvgData, s=20, facecolor='lightgray', edgecolor='k',alpha=0.85)

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
                ax.set_axis_bgcolor('lightcyan')
                
            elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                ax.set_axis_bgcolor('lightgreen')
               
            elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                ax.set_axis_bgcolor('mistyrose')
                
            elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                ax.set_axis_bgcolor('cornsilk')
                
            elif (Lat[i] < -50.):
                ax.set_axis_bgcolor('lightgrey')
            
            else:
                ax.set_axis_bgcolor('lightgrey')
 
        if fits:

            if pltFig:
                #ax.plot(dates2,f_fourier(dateYearFrac2),label='Intra-annual variability', linewidth=1.0)
                ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
            #---------------------------------------------------

            if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Bremen') or  (idhdf == 'Paris') or (idhdf == 'Altzomoni') or (idhdf == 'Tsukuba') or (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paramaribo'):

                yoi  = [[2010, 2016]]

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

                yoi  = [[2002, 2008], [2010, 2016]]

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

                yoi  = [[1995, 2002], [2002, 2008], [2010, 2016]]
            
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

            ax.grid(True, color='gray', alpha=0.5)
            ax.tick_params(which='both',labelsize=10)
            ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
            #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
            #start, end = ax1[i].get_xlim()
            #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
            #ax1[i].set_ylim(bottom=0)

            yearsLc      = YearLocator()
            months       = MonthLocator()
            DateFmt      = DateFormatter('%Y')

            #plt.xticks(rotation=45)
            ax.xaxis.set_major_locator(yearsLc)
            #ax.xaxis.set_minor_locator(months)
            ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax.xaxis.set_major_formatter(DateFmt)
            #ax.set_xlabel('Year')
            #ax1.xaxis.set_tick_params(which='major', pad=15)  
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')

            if (ymin and ymax):
                ax.set_ylim(ymin, ymax)
          
            #ax.set_xticks([])
            #ax.set_yticks([])
            fig.add_subplot(ax)

    if pltFig:
    
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
            if ax.is_first_row():
                ax.spines['top'].set_visible(True)
            if ax.is_last_row():
                ax.spines['bottom'].set_visible(True)
                plt.setp(ax.get_xticklabels(), visible=True)
            if ax.is_first_col():
                ax.spines['left'].set_visible(True)
            if ax.is_last_col():
                ax.spines['right'].set_visible(True)

        if (npanels % 2 == 1): #even
            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)

        
        all_axes[-1].set_xlabel('Year')
        all_axes[-2].set_xlabel('Year')
        all_axes[-3].set_xlabel('Year')

        
        fig.autofmt_xdate()

        fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        plt.suptitle(ytypeStr, fontsize=16  )

        fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        
        if saveFlg: pdfsav.savefig(fig,dpi=200)
        else: 
            plt.show(block=False)


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)


def hbarplt3(b1, b2, b3, pltID, b1_label='', b2_label='', b3_label='', subtitle='', saveFlg=False, pdfsav=' '):
    
    #---------------------------------------------------
    # Bar plot: Three different periods ==> Retrieval
    #---------------------------------------------------
    pltID        = np.asarray(pltID)

    ind = np.arange(len(b1[0]))
    
    fig, (ax, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=(15, 8.5))
    #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    ax.barh(ind-0.27, b1[0], 0.27, xerr=b1[1], align='center', color = 'r', ecolor = 'k', label = b1_label)
    ax.barh(ind, b2[0], 0.27, xerr=b2[1], align='center', color = 'b', ecolor = 'k', label = b2_label)
    ax.barh(ind+0.27, b3[0], 0.27, xerr=b3[1], align='center', color = 'g', ecolor = 'k', label = b3_label)  #yerr=slope_TC*0
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlabel('Rate of change (%/y)')
    ax.set_yticks(ind)
    ax.set_yticklabels(np.transpose(pltID), rotation=0)
    ax.set_title(subtitle+'\n(1996 - 2016)', multialignment='center')
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
    ax2.set_xlabel('Rate of change (%/y))')
    ax2.set_title(subtitle+'\n(1996 - 2002)', multialignment='center')
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
    ax3.set_xlabel('Rate of change (%/y))')
    ax3.set_title(subtitle+'\n(2002 - 2008)', multialignment='center')
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
    ax4.set_xlabel('Rate of change (%/y))')
    ax4.set_title(subtitle+'\n(2009 - 2016)', multialignment='center')
    ax4.set_yticks(ind)
    #ax4.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    ax4.set_xlim(-2.0, 2.0)
    ax4.axvline(0, color='black', lw=1)

    #ax4.invert_yaxis()

    plt.gca().invert_yaxis()
    #fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else:       plt.show(block=False)


 
if __name__ == "__main__":
    main()
