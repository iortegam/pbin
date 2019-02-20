#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltHIPPO_ACE.py
#
# Purpose:
#       1) Read and plot Outputs from HIPPO
#       2) Read and plot Outputs from ACE-FTS
#       3) Initially used to know OCS variablity and create Sa matrix binned by latitude
#
# Notes:
#   
#
# Version History:
#       Created, Sep, 2016  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io               import netcdf
import os
import datetime             as dt
import numpy                as np
import numpy.ma             as ma
import sys
import glob

from scipy                   import interpolate

import matplotlib.dates      as md
import matplotlib.dates      as md
from matplotlib.dates        import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY

import matplotlib.pyplot     as plt
from matplotlib.ticker       import MaxNLocator
from matplotlib.ticker       import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm         as mplcm
import matplotlib.colors     as colors
import matplotlib.gridspec   as gridspec
import matplotlib.gridspec   as gridspec
from itertools               import izip
from numpy                   import *
import myfunctions           as mf
from netCDF4                 import Dataset
import HIPPO_ACE_outClass    as hc
from collections             import OrderedDict
from scipy                   import linspace, polyval, polyfit, sqrt, stats, randn
import shutil
import subprocess as sp
from scipy.optimize import curve_fit

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

def chMod(PrntPath):
    for dirpath, dirnames, filenames in os.walk(PrntPath):
        try:    os.chmod(dirpath,0o777)
        except: pass        
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            try:    os.chmod(path, 0o777)   
            except: pass 

def closeFig(self):
    self.pdfsav.close()

def subProcRun( sysCall, logF=False, shellFlg=False ):
    '''This runs a system command and directs the stdout and stderr'''
    rtn = sp.Popen( sysCall, stdout=sp.PIPE, stderr=sp.PIPE, shell=shellFlg )
   
    stdoutInfo, stderrInfo = rtn.communicate()

    print stdoutInfo
    print stderrInfo


    if logF:
        if stdoutInfo: logF.info( stdoutInfo )
        if stderrInfo: logF.error( stderrInfo )
               
    return (stdoutInfo,stderrInfo)

def donwloadACE(dataDir, iyear, fyear):
    #-------------------------------------
    #INITIALLY USED TO RETRIEVE SINGLE FILES FROM ACE-FTS v.2.
    #HOWEVER, THE SINGLE NETCDF FILES FROM v3.5 ARE BETTER:  https://databace.scisat.ca/level2/ace_v3.5/
    #-------------------------------------
    i_date   = dt.date(iyear,1,1)                                                     
    f_date   = dt.date(fyear,12,31)

    ndays = (f_date + dt.timedelta(days=1) - i_date).days
    ListD =[i_date + dt.timedelta(days=i) for i in range(0, ndays, 1)]
    years = [ singDate.year for singDate in ListD]               
    years = list(set(years))                                            

    months = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']

    chMod(dataDir)


    for y in years:
        for i in months:
            #print 'http://www.ace.uwaterloo.ca/ACE-FTS_v2.2/'+str(y)+'-'+str(i)+'.tgz'
            #cmnd = ['wget','-r','-l1','-nd','-N','--no-parent','-P'+dataDir, 'http://www.ace.uwaterloo.ca/ACE-FTS_v2.2/'+str(y)+'-'+str(i)+'.tgz'  ]
            #(stdoutInfo,stderrInfo) = subProcRun( cmnd )

            cmnd = ['wget','-l1','-nc', '--user=iortega','--password=grindcore','--no-check-certificate', '-P'+dataDir, 'https://databace.scisat.ca/level2/ace_v3.5/ASC/'+str(y)+'-'+str(i)+'.zip']
            (stdoutInfo,stderrInfo) = subProcRun( cmnd )

            print stdoutInfo
            print stderrInfo

            cmnd = ['unzip', dataDir+'/'+str(y)+'-'+str(i)+'.zip', '-d '+dataDir]
          
            (stdoutInfo,stderrInfo) = subProcRun( cmnd )
            chMod(dataDir+str(y)+'-'+str(i))

            print stdoutInfo
            print stderrInfo

            #exit()

def model_func(t, A, K, C):
    return A * np.exp(K * t) + C

def segmnt(seq,n):
    '''Yeilds successive n-sized segments from seq'''
    for i in xrange(0,len(seq),n):
        yield seq[i:i+n] 


                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():

    #-------------------------------------
    #INITIALIZATIONS FOR HIPPO
    #-------------------------------------
    dataDir       = '/data1/Campaign/HIPPO/'
    file          = 'HIPPO_discrete_continuous_merge_20121129.tbl'   #MERGE FILE: http://data.eol.ucar.edu/master_list/?project=HIPPO-2
    
    HIPPOFlg       = False

    pltPr2Flg      = True    #PLOT PROFILES BINNED BY LATITUDE
    pltRFFlg       = False   #PLOT OF MAP AWITH RESEARCH FLIGHT
    pltPrFlg       = False   #PLOT ALL PROFILES (A LOT OF PROFILES)

    #-------------------------------------
    #INITIALIZATIONS FOR ACE-FTS
    #-------------------------------------
    dataDirACE    = '/data1/iortega/ACEv3p5/'
    gasACE        = 'ocs'
    iyearACE      = 2004
    imonthACE     = 1
    fyearACE      = 2013
    fmonthACE     = 12
    
    ACEFlg        = True
    
    pltPrACEFlg   = True
    DownFlg       = False

    #-------------------------------------
    #INDEPENDENT FLAGS
    #-------------------------------------
    saveFlg       = False
    gas           = 'ocs'
    #BinLat        = [ [-90.0, -30.0], [-30.0, 30.0], [30.0, 90.0] ]
    #BinID         = ['South', 'Tropics', 'North']
    #ClrID         = ['r', 'b', 'g']

    BinLat        = [ [-90.0, -50.0], [-50.0, -20.0], [-20.0, 20.0], [20.0, 50.0], [50.0, 90.0] ]
    BinID         = ['90-50S', '50-20S', '20S-20N', '20-50N', '50-90N']
    ClrID         = ['r', 'b', 'g', 'c', 'm', 'y']

    #BinLat        = [ [-90.0, -70.0], [-70.0, -50.0], [-50.0, -30.0], [-30.0, -10.0], [-10.0, 10.0], [10.0, 30.0], [30.0, 50.0], [50.0, 70.0], [70.0, 90.0] ]
    #BinID         = ['90-70S', '70-50S', '50-30S', '30-10S', '10S-10N', '10-30N', '30-50N', '50-70N', '70-90N']
    #ClrID         = ['r', 'b', 'g', 'c', 'm', 'y', 'k', 'gray', 'orange', 'brown' ]
    
    #pltFile       = '/data1/Campaign/HIPPO/'+gas.upper()+'_HIPPO_Plts.pdf'

    printTXT       = False

    outtxtpath    = '/data/iortega/pbin/HIPPO_ACE/'


                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    if saveFlg: pdfsav = PdfPages('/data/iortega/pbin/HIPPO_ACE/'+gas.upper()+'_HIPPO_ACEv3p5_Plts.pdf')

    pltFile        = '/data/iortega/pbin/HIPPO_ACE/'+gas.upper()+'_HIPPO_Plts.pdf'
    pltFileACE     = '/data/iortega/pbin/HIPPO_ACE/'+gas.upper()+'_ACEv3.5_Plts.pdf'

    if HIPPOFlg:

        statDataCl = OrderedDict()

        Data     = hc.HIPPOClass(dataDir, file,   outFname= pltFile, saveFlg= saveFlg)
        Data.ReadOutputHIPPO(gasname=gas)
        if pltRFFlg:   Data.pltMapRF()
        if pltPrFlg:   Data.pltPrf()
        if pltPr2Flg:  Data.pltPrf2(BinLat=BinLat, BinID=BinID, ClrID=ClrID)

        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()
        
       
    if ACEFlg:
        DataACE   = hc.ACEClass(dataDirACE, iyear=iyearACE,imnth=imonthACE, iday=1, fyear=fyearACE,fmnth=fmonthACE,fday=31, outFname= pltFileACE, saveFlg= saveFlg)
        ####DataACE.ReadNCDFACE(gasname=gas)      #INITAL NETCDF FROM THE PUBLI RELEASE (v2.2)
        ####DataACE.ReadASCIIACE(gasname=gas)     #ASCII FILES FROM v3.5 (DO NOT CONTAIN IMPORTANT FLAGS)
        DataACE.ReadNCDFACE2(gasname=gas)         #NETCDF FROM v3.5 (CONTAINS IMPORTANT FLAGS)
        if pltPrACEFlg:   DataACE.pltPrfACE(BinLat=BinLat, BinID=BinID, ClrID=ClrID)

        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()



    if HIPPOFlg & ACEFlg:
        #----------------------------------------------------------------------------
        #PROFILES OF BOTH HIPPO AND ACE
        #---------------------------------------------------------------------------

        fig,   ax = plt.subplots(1,2, figsize=(11, 7))
        fig2, ax2 = plt.subplots(1,2, figsize=(11, 7))
        fig3, ax3 = plt.subplots(1,2, figsize=(11, 7))

        PrfFinal   = {}
        SaFinal    = {}
        AltFinal   = {}
        BinFinal   = {}

        for i, bi in enumerate(BinID):

            #---HIPPO
            gas_mean_i = Data.PrfLat[bi+'_Prf_mean'][0]
            gas_var_i  = Data.PrfLat[bi+'_Prf_var'][0]  #np.nanvar(PrfLat[bi+'_Prf'], axis=0)
            gas_std_i  = Data.PrfLat[bi+'_Prf_std'][0]  #np.nanstd(PrfLat[bi+'_Prf'], axis=0)
            gas_med_i  = Data.PrfLat[bi+'_Prf_med'][0]  #np.nanmedian(PrfLat[bi+'_Prf'], axis=0)
            layers_mid = Data.layers_mid


            Fraction = gas_std_i/gas_mean_i

            #---ACE
            gas_mean_ACE_i = DataACE.PrfLat[bi+'_Prf_mean'][0]
            gas_var_ACE_i  = DataACE.PrfLat[bi+'_Prf_var'][0]  #np.nanvar(PrfLat[bi+'_Prf'], axis=0)
            gas_std_ACE_i  = DataACE.PrfLat[bi+'_Prf_std'][0]  #np.nanstd(PrfLat[bi+'_Prf'], axis=0)
            gas_med_ACE_i  = DataACE.PrfLat[bi+'_Prf_med'][0]  #np.nanmedian(PrfLat[bi+'_Prf'], axis=0)
            layers_ACE_mid = DataACE.layers_mid

            Fraction_ACE = gas_std_ACE_i/gas_mean_ACE_i

            #----------------------------------------------------------------------------
            #SMOOTH PROFILES: FIRST PROFILES ARE CONCATENATED, THEN INTERPOLATED TO A SIMILAR GRID AND THEN SMOOTH USING savitzky_golay
            #----------------------------------------------------------------------------

            gas_mean_ACE_i   = gas_mean_ACE_i + gas_mean_ACE_i*0.15
            gas_med_ACE_i    = gas_med_ACE_i + gas_med_ACE_i*0.15

            gas_mean_HA_i   = np.concatenate((gas_mean_i, gas_mean_ACE_i) )
            gas_med_HA_i    = np.concatenate((gas_med_i, gas_med_ACE_i) )
            gas_std_HA_i    = np.concatenate((gas_std_i, gas_std_ACE_i) )
            fraction_HA_i   = np.concatenate((Fraction, Fraction_ACE) ) 
            layers_HA_i     = np.concatenate((layers_mid, layers_ACE_mid) ) 
   
            indsNAN         = np.where(np.isnan(gas_mean_HA_i) == True)
            gas_mean_HA_i   = np.delete(gas_mean_HA_i,indsNAN )
            layers_HA_i     = np.delete(layers_HA_i,  indsNAN )
            gas_std_HA_i    = np.delete(gas_std_HA_i,  indsNAN )
            gas_med_HA_i    = np.delete(gas_med_HA_i,  indsNAN )

            if bi == '90-50S' or bi == '50-90N': layers_HA_i2    = np.arange(0.5, 24, 1.0)
            else: layers_HA_i2    = np.arange(0.5, 28, 1.0)

            #----
            tck_mean = interpolate.splrep(layers_HA_i, gas_mean_HA_i, w= 1/gas_std_HA_i, k=1, s=0)
            tck_med = interpolate.splrep(layers_HA_i, gas_med_HA_i, w= 1/gas_std_HA_i, k=1, s=0)
            tck_std  = interpolate.splrep(layers_HA_i, gas_std_HA_i, k=1, s=0) 
            
            gas_mean_HA_i = interpolate.splev(layers_HA_i2, tck_mean, der=0)
            gas_std_HA_i = interpolate.splev(layers_HA_i2, tck_std, der=0)
            gas_med_HA_i = interpolate.splev(layers_HA_i2, tck_med, der=0)
            
            fraction_HA_i = gas_std_HA_i/gas_mean_HA_i

            #----------------------------------------------------------------------------
            #SMOOTH:savitzky_golay
            #----------------------------------------------------------------------------
            window_size = 9
            poly        = 3
            gas_mean_HA_i  = mf.savitzky_golay(gas_mean_HA_i, window_size, poly) 
            gas_med_HA_i   = mf.savitzky_golay(gas_med_HA_i, window_size, poly) 
            gas_std_HA_i   = mf.savitzky_golay(gas_std_HA_i, window_size, poly)
            fraction_HA_i  = mf.savitzky_golay(fraction_HA_i, window_size, poly)

            #----------------------------------------------------------------------------
            #GO UP TO 120km
            #----------------------------------------------------------------------------
            layers_UP        = np.arange(38.5, 119.5, 1.0)
            gas_mean_UP      = np.zeros(np.shape(layers_UP))

            gas_mean_UP[0]   = 0.138
            gas_mean_UP[1]   = 0.05

            inds             = np.where(layers_HA_i2 >= 18.0)

            gas_mean_HA_up   = np.concatenate((gas_mean_HA_i[inds], gas_mean_UP) )
            layers_HA_UP     = np.concatenate((layers_HA_i2[inds], layers_UP) )

            if i == 0 or i == len(BinID)-1: layers_MID    = np.arange(24.5, 38.0, 1.0)
            else: layers_MID    = np.arange(28.5, 38.0, 1.0)
            
            par = np.poly1d(np.polyfit(layers_HA_UP,gas_mean_HA_up, 11 ))
            fit_y = par(layers_MID)

            gas_mean_all   = np.concatenate((gas_mean_HA_i, fit_y, gas_mean_UP) )
            layers_all     = np.concatenate((layers_HA_i2, layers_MID, layers_UP) )

            gas_mean_all  = mf.savitzky_golay(gas_mean_all, 7, 3)

            gas_mean_all[gas_mean_all< 0. ] = 0.0

            #------
            #Sa
            #------

            SaYutting       = np.flipud(Data.SaTAB)
            zYutting        = np.flipud(Data.zTAB)
            midpointYutting = np.flipud(Data.midpointTAB)

            indsZ   = np.where(midpointYutting > 28.5)[0]
            SaYutting = SaYutting[indsZ]
            midpointYutting = midpointYutting[indsZ]

            if bi == '90-50S':factor = 2.9
            elif bi == '50-20S': factor = 2.2
            elif bi == '20S-20N': factor = 1.0
            elif bi == '20-50N': factor = 1.9
            elif bi == '50-90N': factor = 2.0

            SaYutting = SaYutting*factor

            Fraction_UP    = np.concatenate( (fraction_HA_i, SaYutting))
            layers_HA_UP_f = np.concatenate( (layers_HA_i2, midpointYutting) )


            tck_mean = interpolate.splrep(layers_HA_UP_f, Fraction_UP,  k=3, s=0)
            Fraction_UP = interpolate.splev(layers_all, tck_mean, der=0)


            #Fraction_UP     = interpolate.interp1d(layers_HA_UP_f, Fraction_UP, bounds_error=False, kind='cubic')(layers_UP)
            Fraction_UP     = mf.savitzky_golay(Fraction_UP, 9, 3)

            #------------------------------------------
            #SAVE FINAL PROFILES 
            #------------------------------------------

            #PrfFinal.setdefault('Prf_'+bi,[]).append(np.asarray(gas_mean_all))
            #SaFinal.setdefault('Sa_'+bi,[]).append(np.asarray(Fraction_UP))
            #AltFinal.setdefault('Altitude_'+bi,[]).append(np.asarray(layers_all))
            #BinFinal.setdefault(bi,[]).append(bi)

            PrfFinal.setdefault('Prf_'+bi, np.asarray(gas_mean_all))
            SaFinal.setdefault('Sa_'+bi,np.asarray(Fraction_UP))
            AltFinal.setdefault('Altitude_'+bi,np.asarray(layers_all))
            BinFinal.setdefault(bi,bi)


            #------------------------------------------
            #PRINT OCS PROFILE FOR MLO/TAB/FL0
            #------------------------------------------
            if bi == '50-90N':
                zTAB    = Data.zTAB
                zTAB    = np.flipud(zTAB)
                PrfTAB     = interpolate.interp1d(layers_all, gas_mean_all, fill_value='extrapolate', bounds_error=False, kind='linear')(zTAB)
                PrfTAB = np.flipud(PrfTAB)/1e12
                PrfTAB[PrfTAB==0.0] = 1.472e-14

                if printTXT:

                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_TAB.refprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(PrfTAB, 5):
                            strformat = ','.join('{:>12.3E}' for rrr in row) + ', \n'
                            fopen.write(strformat.format(*row))

                    zTAB_Mid = np.flipud(Data.midpointTAB)

                    SaTAB     = interpolate.interp1d(layers_all, Fraction_UP, fill_value='extrapolate', bounds_error=False, kind='linear')(zTAB_Mid)
                    SaTAB = np.flipud(SaTAB)
                    #SaTAB[SaTAB==0.0] = 1.472e-14

                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_TAB.Saprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(SaTAB, 10):
                            strformat = ' '.join('{:>8.4f}' for rrr in row) +  '\n'
                            fopen.write(strformat.format(*row))

            elif bi == '20S-20N':

                zMLO    = Data.zMLO
                zMLO    = np.flipud(zMLO)
                PrfMLO     = interpolate.interp1d(layers_all, gas_mean_all, fill_value='extrapolate', bounds_error=False, kind='linear')(zMLO)
                PrfMLO = np.flipud(PrfMLO)/1e12
                PrfMLO[PrfMLO==0.0] = 1.472e-14

                if printTXT:
                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_MLO.refprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(PrfMLO, 5):
                            strformat = ','.join('{:>12.3E}' for rrr in row) + ', \n'
                            fopen.write(strformat.format(*row))

                    zMLO_Mid = np.flipud(Data.midpointMLO)

                    SaMLO     = interpolate.interp1d(layers_all, Fraction_UP, fill_value='extrapolate', bounds_error=False, kind='linear')(zMLO_Mid)
                    SaMLO = np.flipud(SaMLO)
                    #SaMLO[SaMLO==0.0] = 1.472e-14

                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_MLO.Saprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(SaMLO, 10):
                            strformat = ' '.join('{:<8.4f}' for rrr in row) +  '\n'
                            fopen.write(strformat.format(*row))

            elif bi == '20-50N':

                zFL0    = Data.zFL0
                zFL0    = np.flipud(zFL0)
                PrfFL0     = interpolate.interp1d(layers_all, gas_mean_all, fill_value='extrapolate', bounds_error=False, kind='linear')(zFL0)
                PrfFL0 = np.flipud(PrfFL0)/1e12
                PrfFL0[PrfFL0==0.0] = 1.472e-14

                if printTXT:
                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_FL0.refprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(PrfFL0, 5):
                            strformat = ','.join('{:>12.3E}' for rrr in row) + ', \n'
                            fopen.write(strformat.format(*row))

                    zFL0_Mid = np.flipud(Data.midpointFL0)

                    SaFL0     = interpolate.interp1d(layers_all, Fraction_UP, fill_value='extrapolate', bounds_error=False, kind='linear')(zFL0_Mid)
                    SaFL0 = np.flipud(SaFL0)
                    #SaMLO[SaMLO==0.0] = 1.472e-14

                    outtxtFile = outtxtpath + 'OCS_Prf_'+bi+'_FL0.Saprfs'

                    with open(outtxtFile,'w') as fopen:
                        for row in segmnt(SaFL0, 10):
                            strformat = ' '.join('{:<8.4f}' for rrr in row) +  '\n'
                            fopen.write(strformat.format(*row))



            #---HIPPO
            ax[0].plot(gas_mean_i, layers_mid,  color=ClrID[i], markersize=4)
            ax[0].errorbar(gas_mean_i, layers_mid,fmt='o',xerr=gas_std_i, color=ClrID[i], ecolor=ClrID[i], label='HIPPO-'+BinID[i])
            ax[0].scatter(gas_med_i, layers_mid, edgecolors=ClrID[i],marker='x', facecolors=ClrID[i], s=35)

            ax[1].plot(Fraction, layers_mid, 'o', linestyle='-', markersize=6, color=ClrID[i], label='HIPPO-'+BinID[i])

            #---ACE
            ax[0].plot(gas_mean_ACE_i, layers_ACE_mid, linestyle='dotted',  color=ClrID[i],  markersize=4)
            ax[0].errorbar(gas_mean_ACE_i, layers_ACE_mid,fmt='s',xerr=gas_std_ACE_i, color=ClrID[i], ecolor=ClrID[i], label='ACE-'+BinID[i] )
            ax[0].scatter(gas_med_ACE_i, layers_ACE_mid, edgecolors=ClrID[i],marker='x', facecolors=ClrID[i], s=35)

            ax[1].plot(Fraction_ACE, layers_ACE_mid, 's', linestyle='dotted', markersize=6, color=ClrID[i], label='ACE-'+BinID[i])

            #---BOTH
            ax2[0].plot(gas_mean_HA_i, layers_HA_i2,  color=ClrID[i], markersize=0)
            ax2[0].errorbar(gas_mean_HA_i, layers_HA_i2,fmt='o',xerr=gas_std_HA_i, color=ClrID[i], ecolor=ClrID[i], label=BinID[i])
            ax2[0].scatter(gas_med_HA_i, layers_HA_i2, edgecolors=ClrID[i],marker='x', facecolors=ClrID[i], s=35)

            ax2[1].plot(fraction_HA_i, layers_HA_i2, 'o', linestyle='-', markersize=6, color=ClrID[i], label=BinID[i])

            #ax2[0].plot(fit_y, layers_MID, linewidth=6,  color=ClrID[i])
            #ax2[0].errorbar(fit_y, layers_MID,fmt='o', color=ClrID[i], ecolor=ClrID[i])
            #ax2[0].scatter(fit_y, layers_MID, edgecolors=ClrID[i],marker='o', facecolors=ClrID[i], s=35)

            
            #ax3[0].scatter(gas_mean_all, layers_all, color=ClrID[i], edgecolors='k',marker='o', facecolors=ClrID[i], s=45)
            #ax3[0].plot(gas_mean_all, layers_all, linewidth=1,  color=ClrID[i], label=BinID[i])
            ax3[0].plot(gas_mean_all, layers_all, 'o', linestyle='-', markersize=6, color=ClrID[i], label=BinID[i])

            ax3[1].plot(Fraction_UP, layers_all, 'o', linestyle='-', markersize=6, color=ClrID[i], label=BinID[i])



        #----------------------------------------------------------------------------
        #PROFILES FROM THULE AND MLO (GOTTEN FROM HIPPO CLASS)
        #---------------------------------------------------------------------------
        aprprfTAB = Data.aprprfTAB
        SaTAB     = Data.SaTAB
        zTAB      = Data.zTAB
        midpointTAB = Data.midpointTAB

        #print aprprfTAB*1e12

        aprprfMLO = Data.aprprfMLO
        SaMLO     = Data.SaMLO
        zMLO     = Data.zMLO
        midpointMLO = Data.midpointMLO

        # ax[0].plot(aprprfTAB*1e12, zTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Toon')
        # ax[1].plot(SaTAB, midpointTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Yutting et al, (2016)')

        # ax[0].plot(aprprfMLO*1e12, zMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Toon')
        # ax[1].plot(SaMLO, midpointMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Yutting et al, (2016)')

        # #----
        # ax2[0].plot(aprprfTAB*1e12, zTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Toon')
        # ax2[1].plot(SaTAB, midpointTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Yutting et al, (2016)')

        # ax2[0].plot(aprprfMLO*1e12, zMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Toon')
        # ax2[1].plot(SaMLO, midpointMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Yutting et al, (2016)')

        # #----
        # ax3[0].plot(aprprfTAB*1e12, zTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Toon')
        # ax3[1].plot(SaTAB, midpointTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Yutting et al, (2016)')

        # ax3[0].plot(aprprfMLO*1e12, zMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Toon')
        # ax3[1].plot(SaMLO, midpointMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Yutting et al, (2016)')


        #----------------------------------------------------------------------------
        ax[0].tick_params(labelsize=14)  
        ax[0].grid(True)
        ax[0].set_xlabel('VMR [ppt]', fontsize = 14)
        ax[0].set_ylabel('Altitude [km]', fontsize = 14)
        ax[0].set_ylim(0,35)
        ax[0].set_xlim(xmin=-10, xmax=600)
        ax[0].legend(prop={'size':10}, loc=3)

        ax[1].tick_params(labelsize=14)
        ax[1].grid(True)
        ax[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
        ax[1].set_ylim(0, 35)
        ax[1].legend(prop={'size':10}, loc='lower right')
        ax[1].set_xlim(xmin=0, xmax=1.0)

        ax2[0].tick_params(labelsize=14)  
        ax2[0].grid(True)
        ax2[0].set_xlabel('VMR [ppt]', fontsize = 14)
        ax2[0].set_ylabel('Altitude [km]', fontsize = 14)
        ax2[0].set_ylim(0,50)
        ax2[0].set_xlim(xmin=-10, xmax=600)
        ax2[0].legend(prop={'size':10}, loc=3)

        ax2[1].tick_params(labelsize=14)
        ax2[1].grid(True)
        ax2[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
        ax2[1].set_ylim(0, 50)
        ax2[1].legend(prop={'size':10}, loc='lower right')
        ax2[1].set_xlim(xmin=0, xmax=1.0)

        ax3[0].tick_params(labelsize=14)  
        ax3[0].grid(True)
        ax3[0].set_xlabel('VMR [ppt]', fontsize = 14)
        ax3[0].set_ylabel('Altitude [km]', fontsize = 14)
        ax3[0].set_ylim(0,50)
        ax3[0].set_xlim(xmin=-10, xmax=600)
        ax3[0].legend(prop={'size':10}, loc=3)

        ax3[1].tick_params(labelsize=14)
        ax3[1].grid(True)
        ax3[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
        ax3[1].set_ylim(0, 50)
        ax3[1].legend(prop={'size':10}, loc='lower right')
        ax3[1].set_xlim(xmin=0, xmax=1.0)
        
        fig.subplots_adjust(left = 0.125, bottom=0.15, top=0.95, right = 0.95)
        fig2.subplots_adjust(left = 0.125, bottom=0.15, top=0.95, right = 0.95)
        fig3.subplots_adjust(left = 0.125, bottom=0.15, top=0.95, right = 0.95)
            
        if saveFlg: 
           pdfsav.savefig(fig,dpi=200)
           pdfsav.savefig(fig2,dpi=200)
           pdfsav.savefig(fig3,dpi=200)
           
        
        else:           
            plt.show(block=False)

    #------------------------------------------
    #PRINT OCS PROFILE BINNED BY LATITUDE (TO DISTRIBUTE)
    #------------------------------------------



        #PrfFinal   = np.asarray(PrfFinal)
        #SaFinal    = np.asarray(SaFinal)  
        #AltFinal   = np.asarray(AltFinal)
      
        #DataAll = [  [AltFinal[i][0],PrfFinal[i][0],SaFinal[i][0]] for i in np.arange(0, Nbins, 1)]
      #  DataAll = zip(*(AltFinal[:][0],PrfFinal[:][0],SaFinal[:][0]))#for i in np.arange(0, Nbins, 1)]
        
        DataAll = zip(*(AltFinal['Altitude_'+BinID[0]], PrfFinal['Prf_'+BinID[0]]/1e12,SaFinal['Sa_'+BinID[0]], PrfFinal['Prf_'+BinID[1]]/1e12,SaFinal['Sa_'+BinID[1]], PrfFinal['Prf_'+BinID[2]]/1e12,SaFinal['Sa_'+BinID[2]], PrfFinal['Prf_'+BinID[3]]/1e12,SaFinal['Sa_'+BinID[3]], PrfFinal['Prf_'+BinID[4]]/1e12,SaFinal['Sa_'+BinID[4]]))

        #DataAll = [zip(*(AltFinal['Altitude_'+bi], PrfFinal['Prf_'+bi],SaFinal['Sa_'+bi])) for bi in BinID]
        
        
        #DataAll = [zip(*(AltFinal,PrfFinal[i],SaFinal[i])) for bi in BinID]

        if printTXT:
            outtxtFile = outtxtpath + 'OCS_HIPPO_ACE.prfs'

            with open(outtxtFile,'w') as fopen:
                fopen.write('#Hannigan, J.W., Ortega, I\n')
                fopen.write('#National Center for Atmospheric Research\n')
                fopen.write('#CONTACT_INFO:   Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
                fopen.write('#CONTACT_INFO:   Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
                fopen.write('#DESCRIPTION:    OCS vertical profiles [mixing ratios] and Sa [unitless] binned by latitudes derived from in-situ HIPPO observation (troposphere) and ACE-FTS (stratosphere)\n')
                fopen.write('#COMMENTS:       Latitude bins: 90-50S, 50-20S, 20S-20N, 20-50N, 50-90N\n')
                #fopen.write('#For details see ***.docx\n')
                #fopen.write('{0:>20s}{1:>30s}{2:>30s}{3:>30s}{4:>30s}\n'.format('Altitude [km]', 'Latitude - '+ BinFinal[BinID[0]],'Latitude - '+BinFinal[BinID[1]], 'Latitude - '+BinFinal[BinID[2]], 'Latitude - '+BinFinal[BinID[3]], 'Latitude - '+BinFinal[BinID[4]]))
                fopen.write('{0:<20s}{1:<20s}{2:<20s}{3:<20s}{4:<20s}{5:<20s}{6:<20s}{7:<20s}{8:<20s}{9:<20s}{10:<20s}\n'.format('Altitude [km]','Profile-'+BinFinal[BinID[0]],'Sa-'+BinFinal[BinID[0]],'Profile-'+BinFinal[BinID[1]],'Sa-'+BinFinal[BinID[1]],'Profile-'+BinFinal[BinID[2]],'Sa-'+BinFinal[BinID[2]],'Profile-'+BinFinal[BinID[3]],'Sa-'+BinFinal[BinID[3]],'Profile-'+BinFinal[BinID[4]],'Sa-'+BinFinal[BinID[4]]))

                for row in DataAll  :
                    strformat = ''.join('{:<20.4E}' for i in row) + '\n'
                    fopen.write(strformat.format(*row))


    if DownFlg: donwloadACE(dataDirACE, iyearACE, fyearACE)


    if saveFlg:
        pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program 

if __name__ == "__main__":
    main()
