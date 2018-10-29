#! /usr/bin/python2.7
#------------------------------------------------------
#SCRIPT TO PLOT IN-SITU Vs FTIR
#------------------------------------------------------

import sys
import datetime as dt
import matplotlib
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
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
import numpy as np

from cycler import cycler

from itertools import izip
import myfunctions as mf
import dataModelOutClass as dm
import dataOutClass as dc
from collections                     import OrderedDict
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

from scipy.io import netcdf
import os
import glob

from scipy import interpolate

import matplotlib.dates as md
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl

import pytz
import csv

from numpy import *

def month_string_to_number(string):
    m = {
        'jan': 1,
        'feb': 2,
        'mar': 3,
        'apr':4,
         'may':5,
         'jun':6,
         'jul':7,
         'aug':8,
         'sep':9,
         'oct':10,
         'nov':11,
         'dec':12
        }
    s = string.strip()[:3].lower()

    try:
        out = m[s]
        return out
    except:
        raise ValueError('Not a month')

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


def main():

    #-----------------------------------------------------------------------------------------
    #                 Initializations for FTIR
    #-----------------------------------------------------------------------------------------

    loc               = 'mlo'
    gasName           = ['h2o', 'hdo']
    ver               = ['Current_v6_50', 'Current_v8_50'] 
    ctlF              = ['sfit4_v6.ctl', 'sfit4_v6.ctl'] 
    maxRMS            = [0.5, 0.5]                      
    minDOF            = [1.0, 1.0]

    saveFlg           = True 

    insFlg            = False   # Read In-situ?

    if insFlg: pltFile  =  '/data/iortega/results/'+loc.lower()+'/fig/dD_'+loc.upper()+'_v8_2016_wInsitu.pdf'
    else: pltFile       =  '/data/iortega/results/'+loc.lower()+'/fig/dD_'+loc.upper()+'_v8_2016_vPartial_50_2.pdf'

    aoi               = 3.4   #Altitude of interest for single grid value
    maxalt            = 130.0  # Max Altitude 

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
    tcMMFlg            = True                  # Flag to filter based on min and max total column amount
    cnvrgFlg           = True                   # Flag to filter profiles that did not converge
    rmsFlg             = True                   # Flag to filter based on max RMS
    chiFlg             = False                  # Flag to filter based on max CHI_2_Y
    mnthFlg            = False                  # Flag to filter based on 
    
    mnths              = [6,7,8]                # Months to filter on (these are the months to include data)
    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E6                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppmv'                 # Name of scale factor for labeling plots

    #pCols = [ [3.0, 5.0], [5.0, 16.0], [3.0, 16.0] ]
    pCols = [ [3.3, 4.1], [4.8, 5.55], [6.3, 7.2], [7.9, 8.9],
              [3.3, 4.9], [5.5, 7.2], [7.9, 9.8], [10.7, 12.8],
              [3.3, 5.55], [6.3, 8.9], [9.7, 12.72], [13.7, 17.2] ]

    if saveFlg: pdfsav = PdfPages(pltFile)

    if insFlg:

        aoi               = 3.4

        #-------------------------------------------------------
        #READ IN-SITU NDATA - H2O (ppm)
        #-------------------------------------------------------
        File = '/data1/ancillary_data/mlo/mlo_vapor_ppmv.dat'    

        utc = pytz.utc
        Hawaii = pytz.timezone('US/Hawaii')

        with open(File ,'r') as fopen: lines = fopen.readlines()

        Data = [row.strip().split() for row in lines]

        Date_in   = []
        insitu   = []

        for i, row in enumerate(Data):

            ho   =  int(row[1][0:2].strip())
            mi   =  int(row[1][3:5].strip())
            ss   =  int(row[1][6:8].strip())

            #if (ho <= 7) or (ho <= 10) or (ho >= 14):
            
            insitu.append(float(row[2]))
            yyyy = int(row[0][7:11].strip())
            smm  = str(row[0][3:6].strip())
            mm   = month_string_to_number(smm)
            dd   = int(row[0][0:2].strip())

            date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
            date_lt  = Hawaii.localize(date_lt)

            Date_in.append(date_lt.astimezone(utc))

        insitu    = np.asarray(insitu)
        Date_in   = np.asarray(Date_in)

        #-------------------------------------------------------
        #READ IN-SITU NDATA - dD
        #-------------------------------------------------------
        File = '/data1/ancillary_data/mlo/dD_LGR_MLO.txt'

        with open(File ,'r') as fopen: lines = fopen.readlines()

        Data = [row.strip().split() for row in lines]

        Date_in2   = []
        insitu2    = []

        for i, row in enumerate(Data):

            ho   =  int(row[1][0:2].strip())
            mi   =  int(row[1][3:5].strip())
            ss   =  int(row[1][6:8].strip())

            insitu2.append(float(row[2]))
            yyyy = int(row[0][7:11].strip())
            smm  = str(row[0][3:6].strip())
            mm   = month_string_to_number(smm)
            dd   = int(row[0][0:2].strip())

            date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
            date_lt  = Hawaii.localize(date_lt)

            Date_in2.append(date_lt.astimezone(utc))

        dDinsitu    = np.asarray(insitu2)
        #Date_in2   = np.asarray(Date_in2)

        #-------------------------------------------------------
        #READ IN-SITU NDATA - HDO
        #-------------------------------------------------------
        File = '/data1/ancillary_data/mlo/mlo_hdo.dat'

        with open(File ,'r') as fopen: lines = fopen.readlines()

        Data = [row.strip().split() for row in lines]

        Date_in3   = []
        insitu3    = []

        for i, row in enumerate(Data):

            ho   =  int(row[1][0:2].strip())
            mi   =  int(row[1][3:5].strip())
            ss   =  int(row[1][6:8].strip())

            insitu3.append(float(row[2]))
            yyyy = int(row[0][7:11].strip())
            smm  = str(row[0][3:6].strip())
            mm   = month_string_to_number(smm)
            dd   = int(row[0][0:2].strip())

            date_lt  = dt.datetime(yyyy, mm, dd, ho, mi, ss)
            date_lt  = Hawaii.localize(date_lt)

            Date_in3.append(date_lt.astimezone(utc))

        HDOinsitu    = np.asarray(insitu3)*2.0   #A CORRECTION OF 2 IS NEEDED

        dDinsitu_r         = (np.divide(HDOinsitu, insitu)/3.1152e-4 - 1.0) *1000.

        #----------------------
        # Date range to process FTIR (depending on insitu dates)
        #----------------------
        # iyear              = Date_in[0].year     
        # imnth              = Date_in[0].month    
        # iday               = Date_in[0].day     
        # fyear              = Date_in[-1].year    
        # fmnth              = Date_in[-1].month   
        # fday               = Date_in[-1].day   

        iyear              = 2016
        imnth              = 5
        iday               = 1
        fyear              = 2016
        fmnth              = 5
        fday               = 31  

    else:    
        iyear              = 2016
        imnth              = 4
        iday               = 1
        fyear              = 2016
        fmnth              = 6
        fday               = 30

    #-------------------------------------------------
    #FTS
    #-------------------------------------------------
    retDir = [ '/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+v+'/' for g,v in izip(gasName,ver)] 
    ctlFile  = ['/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+'x.'+g.lower()+'/'+ ctlF[i] for i, g in enumerate(gasName)]

    #---------------------------
    # Check file and directories
    #---------------------------
    for d in retDir:  ckDir(d,exit=True)
    for c in ctlFile: ckFile(c,exit=True)
    ckDir(os.path.dirname(os.path.realpath(pltFile)),exit=True)

    #-------------------------------------
    # Create instance of output data class   
    #-------------------------------------
    statDataCl = OrderedDict()
    for i,gas in enumerate(gasName):
        statDataCl[gas+'_'+ver[i]] = dc.ReadOutputData(retDir[i],'',ctlFile[i],iyear,imnth,iday,fyear,fmnth,fday)


    #--------------
    # Read profiles
    #--------------
    for g, gasVer in enumerate(statDataCl):

        if gasName[g].lower() == 'h2o':
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'HDO'],retapFlg=1)
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'HDO'],retapFlg=0)
        else: 
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'H2O'],retapFlg=1)
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas, 'H2O'],retapFlg=0)

        if statDataCl[gasVer].empty: 
            print 'No retreivals found for {}. Exiting.......'.format(gasVer)
            sys.exit()

    rPrfVMR      = OrderedDict()
    rPrfMol      = OrderedDict()
    dates        = OrderedDict()
    alt          = OrderedDict()
    Airmass      = OrderedDict()
    waterVMR     = OrderedDict()
    waterMol     = OrderedDict()
    totClmn      = OrderedDict()
    totClmn_e    = OrderedDict()
    TCdryAir     = OrderedDict()
    TCdry        = OrderedDict()
    rPrfDry      = OrderedDict()
    rms          = OrderedDict()
    dofsAvg      = OrderedDict()
    dofsAvg_cs   = OrderedDict()
    LayVMR       = OrderedDict()
    DryLayVMR    = OrderedDict()
    upperHgt     = OrderedDict()
    lowerHgt     = OrderedDict()
    LayThk       = OrderedDict()
    LayDOF       = OrderedDict()
    avkSCF       = OrderedDict()
    avkVMR       = OrderedDict()
    avkSCFav     = OrderedDict()
    avkVMRav     = OrderedDict()
    aPrfVMR      = OrderedDict()
    aPrfMol      = OrderedDict()
    tot_err      = OrderedDict()      
    rand_err     = OrderedDict()
    sys_err      = OrderedDict()
    tot_rnd      = OrderedDict()
    tot_sys      = OrderedDict()
    totClmn_e    = OrderedDict()
    PresPrf      = OrderedDict()
    avkTC        = OrderedDict()

    sys_errvmr   =  OrderedDict()
    rand_errvmr   = OrderedDict()
    tot_errvmr   =  OrderedDict()

    rand_cmpnts  = OrderedDict()
    sys_cmpnts   = OrderedDict()

    rand_cmpnts_vmr  = OrderedDict()
    sys_cmpnts_vmr   = OrderedDict()

    sza           = OrderedDict()
    saa           = OrderedDict()


    for j, gasVer in enumerate(statDataCl):
        statDataCl[gasVer].readsummary()

        rPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        rPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))
        dates[gasVer]    = statDataCl[gasVer].rprfs['date']
        alt[gasVer]      = np.asarray(statDataCl[gasVer].rprfs['Z'][0,:])
        Airmass[gasVer]  = np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])
        waterVMR[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'])
        waterMol[gasVer] = np.asarray(statDataCl[gasVer].aprfs['H2O'] * Airmass[gasVer])
        totClmn[gasVer]  = np.sum(rPrfMol[gasVer],axis=1)
        TCdryAir[gasVer] = np.sum(Airmass[gasVer],axis=1) - np.sum(waterMol[gasVer],axis=1)
        TCdry[gasVer]    = (totClmn[gasVer] / TCdryAir[gasVer]) * sclfct
        aPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        aPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))

        PresPrf[gasVer]  = np.asarray(statDataCl[gasVer].aprfs['PRESSURE'])

        if gasName[j].lower() == 'h2o': rPrfVMR_2  = np.asarray(statDataCl[gasVer].rprfs['HDO']) * sclfct

        #----------------------------------------
        # This is the mixing ratio for DRY AIR!!!
        #----------------------------------------
        rPrfDry[gasVer] = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) / (1.0 - waterVMR[gasVer]) * sclfct       
        
        #----------------------------------
        # Read Summary data (For filtering)
        #----------------------------------
        
        rms[gasVer]     = np.asarray(statDataCl[gasVer].summary[statDataCl[gasVer].PrimaryGas+'_FITRMS'])       

        statDataCl[gasVer].readPbp()

        sza[gasVer]  = np.asarray(statDataCl[gasVer].pbp['sza'])    
        saa[gasVer]  = np.asarray(statDataCl[gasVer].pbp['saa'])    
        
        #--------------------
        # Call to filter data
        #--------------------
        if fltrFlg: statDataCl[gasVer].fltrData(statDataCl[gasVer].PrimaryGas,mxrms=maxRMS[j],rmsFlg=rmsFlg, minDOF=minDOF[j],  dofFlg=dofFlg, 
        tcFlg=tcNegFlg, pcFlg=pcNegFlg , cnvrgFlg=True)
        else:       statDataCl[gasVer].inds = np.array([]) 

        #--------------------------------------------
        # Read Error data to get AVK and profile DOFs
        #-------------------------------------------------
        # Determine if AVK has been created via sfit4 core
        # code or via python error analysis
        #-------------------------------------------------     
        if errorFlg:   # Read AVK from error output
            #statDataCl[gasVer].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=True,avkFlg=True,KbFlg=False)
            statDataCl[gasVer].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False)

            npnts           = np.shape(statDataCl[gasVer].error['Total_Random_Error'])[0]
            nlvls           = np.shape(alt[gasVer])[0]
            
            #-------------------------------------------------
            #Error profiles constructed with the diagonal elements of the covariances matrices
            #-------------------------------------------------
            rand_err[gasVer]   = np.zeros((npnts,nlvls))
            sys_err[gasVer]    = np.zeros((npnts,nlvls))

            for i in range(npnts):
                rand_err[gasVer][i,:] = np.diag(statDataCl[gasVer].error['Total_Random_Error'][i][:,:])
                sys_err[gasVer][i,:]  = np.diag(statDataCl[gasVer].error['Total_Systematic_Error'][i][:,:])

            tot_err[gasVer]     = np.sqrt(rand_err[gasVer] + sys_err[gasVer])            
            sys_err[gasVer]     = np.sqrt(sys_err[gasVer])
            rand_err[gasVer]    = np.sqrt(rand_err[gasVer])

            sys_errvmr[gasVer]  = sys_err[gasVer]/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])*sclfct
            rand_errvmr[gasVer] = rand_err[gasVer]/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS'])*sclfct
            tot_errvmr[gasVer]  = np.sqrt(rand_errvmr[gasVer]**2 + sys_errvmr[gasVer]**2) 

            rand_cmpnts[gasVer] = statDataCl[gasVer].randErrDiag
            sys_cmpnts[gasVer]  = statDataCl[gasVer].sysErrDiag

            rand_cmpnts_vmr[gasVer] = statDataCl[gasVer].randErrDiag 
            sys_cmpnts_vmr[gasVer]  = statDataCl[gasVer].sysErrDiag

            for k in sys_cmpnts_vmr[gasVer]:
                sys_cmpnts_vmr[gasVer][k] = (np.sqrt(sys_cmpnts_vmr[gasVer][k])/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))*sclfct
              
            for k in rand_cmpnts_vmr[gasVer]:
                rand_cmpnts_vmr[gasVer][k] = (np.sqrt(rand_cmpnts_vmr[gasVer][k])/ np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))*sclfct

            #-------------------------------------------------
            #Error in the summary output. Get Total Errors 
            #-------------------------------------------------
            tot_rnd[gasVer]        = np.array(statDataCl[gasVer].error['Total random uncertainty'])
            tot_sys[gasVer]        = np.array(statDataCl[gasVer].error['Total systematic uncertainty'])
            totClmn_e[gasVer]       = np.sqrt(tot_rnd[gasVer]**2 + tot_sys[gasVer]**2)
            
            #err_summary[gasVer]    = statDataCl[gasVer].error
                
            #---------------------
            # Get averaging kernel
            #---------------------   
            avkSCF[gasVer]    = np.asarray(statDataCl[gasVer].error['AVK_scale_factor'])
            avkVMR[gasVer]    = np.asarray(statDataCl[gasVer].error['AVK_vmr']) 
            dofs              = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
            avkSCFav[gasVer]  = np.mean(avkSCF[gasVer],axis=0)    
            avkVMRav[gasVer]  = np.mean(avkVMR[gasVer],axis=0)                   
                
            dofsAvg[gasVer]    = np.diag(avkSCFav[gasVer])
            dofsAvg_cs[gasVer] = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]

            
            
        else:        # Read AVK from sfit4 output (only contains scaled AVK)
            avkSCFi = []
            for d in statDataCl[gasVer].dirLst:
                lines  = dc.tryopen( d + statDataCl[gasVer].ctl['file.out.ak_matrix'][0])
                if not lines: continue
                avkSCFi.append(np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[2:] ] ))
                
            if not statDataCl[gasVer].readPrfFlgApr[statDataCl[gasVer].PrimaryGas]: statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=0)   # Apriori Profiles
                      
            avkSCF[gasVer]  = np.asarray(avkSCFi)
            nobs            = np.shape(avkSCF[gasVer])[0]
            n_layer         = np.shape(avkSCF[gasVer])[1]
            avkVMR[gasVer]  = np.zeros((nobs,n_layer,n_layer))
    
            for obs in range(0,nobs):
                Iapriori        = np.zeros((n_layer,n_layer))
                IaprioriInv     = np.zeros((n_layer,n_layer))
                np.fill_diagonal(Iapriori,statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas.upper()][obs])
                np.fill_diagonal(IaprioriInv, 1.0 / (statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas.upper()][obs]))
                avkVMR[gasVer][obs,:,:] = np.dot(np.dot(Iapriori,np.squeeze(avkSCF[gasVer][obs,:,:])),IaprioriInv)       
                
            
            dofs               = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
            avkSCFav[gasVer]     = np.mean(avkSCF[gasVer],axis=0)
            avkVMRav[gasVer]     = np.mean(avkVMR[gasVer],axis=0)
            
            dofsAvg[gasVer]    = np.diag(avkSCFav[gasVer])
            dofsAvg_cs[gasVer] = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]

        #----------------------------------------
        # Calculate total column averaging kernel
        #----------------------------------------
        nobs            = len(dates[gasVer])
        nlyrs           = len(alt[gasVer])
        avkTC_i           = np.zeros((nobs,nlyrs))
        for i in range(0,nobs):
            AirMtemp  = np.squeeze(Airmass[gasVer][i,:])
            akTemp    = np.squeeze(avkSCF[gasVer][i,:,:])
            AirMinv   = np.diag(1.0/AirMtemp)
            avkTC_i[i,:] = np.dot(np.dot(AirMtemp,akTemp),AirMinv)

        avkTC[gasVer]  = np.asarray(avkTC_i)
            
        #--------------------------------------
        # Remove retrieval data based on filter
        #--------------------------------------
        nfltr             = len(statDataCl[gasVer].inds)
        rms[gasVer]       = np.delete(rms[gasVer],statDataCl[gasVer].inds)
        ntot              = len(rms[gasVer])
        dates[gasVer]     = np.delete(dates[gasVer],statDataCl[gasVer].inds)
        totClmn[gasVer]   = np.delete(totClmn[gasVer],statDataCl[gasVer].inds)

        sza[gasVer]       = np.delete(sza[gasVer],statDataCl[gasVer].inds)
        saa[gasVer]       = np.delete(saa[gasVer],statDataCl[gasVer].inds)

        rPrfVMR[gasVer]   = np.delete(rPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfMol[gasVer]   = np.delete(rPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfDry[gasVer]   = np.delete(rPrfDry[gasVer],statDataCl[gasVer].inds,axis=0)   
        Airmass[gasVer]   = np.delete(Airmass[gasVer],statDataCl[gasVer].inds,axis=0)
        TCdry[gasVer]     = np.delete(TCdry[gasVer],statDataCl[gasVer].inds)

        aPrfVMR[gasVer]   = np.delete(aPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        aPrfMol[gasVer]   = np.delete(aPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)

        PresPrf[gasVer]   = np.delete(PresPrf[gasVer],statDataCl[gasVer].inds,axis=0)   

        if gasName[j].lower() == 'h2o':  rPrfVMR_2 = np.delete(rPrfVMR_2,statDataCl[gasVer].inds,axis=0)

        avkSCF[gasVer]     = np.delete(avkSCF[gasVer],statDataCl[gasVer].inds,axis=0)
        avkVMR[gasVer]     = np.delete(avkVMR[gasVer],statDataCl[gasVer].inds,axis=0) 

        avkTC[gasVer]     = np.delete(avkTC[gasVer],statDataCl[gasVer].inds,axis=0) 

        if errorFlg:
            rand_err[gasVer] = np.delete(rand_err[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_err[gasVer]  = np.delete(sys_err[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_err[gasVer]  = np.delete(tot_err[gasVer],statDataCl[gasVer].inds,axis=0)

            rand_errvmr[gasVer] = np.delete(rand_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_errvmr[gasVer]  = np.delete(sys_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_errvmr[gasVer]  = np.delete(tot_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)

            tot_rnd[gasVer] = np.delete(tot_rnd[gasVer],statDataCl[gasVer].inds)
            tot_sys[gasVer] = np.delete(tot_sys[gasVer],statDataCl[gasVer].inds)
            totClmn_e[gasVer] = np.delete(totClmn_e[gasVer],statDataCl[gasVer].inds)

            
            for k in sys_cmpnts[gasVer]:
                sys_cmpnts[gasVer][k] = np.delete(sys_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                sys_cmpnts_vmr[gasVer][k]  = np.delete(sys_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)
                
            for k in rand_cmpnts[gasVer]:
                rand_cmpnts[gasVer][k] = np.delete(rand_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                rand_cmpnts_vmr[gasVer][k] = np.delete(rand_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)



    #-----------------------------
    #Plot Random error components in VMR
    #-----------------------------

    for j, gasVer in enumerate(statDataCl):
    
        fig, ax1 = plt.subplots(1, 2, figsize=(8,7), sharey=True, sharex=False)

        if statDataCl[gasVer].PrimaryGas.lower() == 'hdo': cf = 3.107e-4
        else: cf = 1.0

        for p in sorted(rand_cmpnts_vmr[gasVer].iterkeys()):

            if npnts > 1:
                errPlt = np.mean(rand_cmpnts_vmr[gasVer][p],axis=0)
                retPrf = np.mean(rPrfVMR[gasVer],axis=0)
            else:
                errPlt = rand_cmpnts_vmr[gasVer][p][0]
                retPrf = rPrfVMR[gasVer][0]

            if p == 'measurement_Random_Error': legend = 'Measurement'
            elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
            elif p == 'temperature_Random_Error': legend = 'Temperature'
            elif p == 'sza_Random_Error': legend = 'SZA'
            elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
            elif p == 'curvature_Random_Error': legend = 'Curvature' #continue
            elif p == 'max_opd_Random_Error': legend = 'Max OPD' #continue     
            else: legend = p
        
            ax1[0].plot(errPlt*cf,alt[gasVer],linewidth=2.0, label=legend)

        ax1[0].plot(np.nanmean(rand_errvmr[gasVer], axis=0)*cf,alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
        ax1[0].set_ylabel('Altitude [km]', fontsize=14)
        ax1[0].set_xlabel('VMR [ppm]', fontsize=14)             
        ax1[0].grid(True,which='both')
        ax1[0].legend(prop={'size':12},loc='upper right')          
        #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
        ax1[0].set_title('Random Components')
        if loc.lower() == 'fl0': ax1[0].set_ylim(1, 15.5)
        if loc.lower() == 'mlo': ax1[0].set_ylim(3, 20)
        ax1[0].tick_params(which='both',labelsize=14)
        #ax1[0].set_xlim(0, 10)

        #-----------------------------
        #systematic error components
        #-----------------------------
        for p in sorted(sys_cmpnts_vmr[gasVer].iterkeys()):

            if npnts > 1:
                errPlt = np.mean(sys_cmpnts_vmr[gasVer][p],axis=0)
                retPrf = np.mean(rPrfVMR[gasVer],axis=0)
            else:
                errPlt = sys_cmpnts_vmr[gasVer][p][0]
                retPrf = rPrfVMR[gasVer][0]
            
            if p == 'temperature_Systematic_Error': legend = 'Temperature'
            elif p == 'smoothing_Systematic_Error':  legend = 'Smoothing'
            #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'            
            elif p == 'linepair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$\gamma$'
            elif p == 'lineint_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$S$'
            elif p == 'linetair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$n$'
            elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn' #continue
            elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
            elif p == 'omega_Systematic_Error': legend = 'Omega' #continue
            else: legend = p

            ax1[1].plot(errPlt*cf, alt[gasVer],linewidth=2.0, label=legend)

        ax1[1].plot(np.nanmean(sys_errvmr[gasVer], axis=0)*cf, alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
        ax1[1].set_xlabel('VMR [ppm]', fontsize=14)            
        ax1[1].grid(True,which='both')
        ax1[1].legend(prop={'size':12}, loc='upper right')
        ax1[1].set_title('Systematic Components')
        if loc.lower() == 'fl0': ax1[1].set_ylim(1, 15.5)
        if loc.lower() == 'mlo': ax1[1].set_ylim(3, 20)
        ax1[1].tick_params(which='both',labelsize=14)
        #ax1[1].set_xlim(0, 10)
        #ax1[1].set_xscale('log')

        plt.suptitle('Uncertainty of {} in VMR'.format(statDataCl[gasVer].PrimaryGas), fontsize=16)

        fig.subplots_adjust(bottom=0.1, top=0.90, left = 0.1, right = 0.95)
        
        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
        else:           
            plt.show(block=False) 
            
    #-----------------------------
    #Plot Random error components in percentage
    #-----------------------------
    for j, gasVer in enumerate(statDataCl):
    
        fig, ax1 = plt.subplots(1, 2, figsize=(8,7), sharey=True, sharex=False)

        #for k, p in enumerate(list_comp_ran):
        for p in sorted(rand_cmpnts_vmr[gasVer].iterkeys()):

            if npnts > 1:
                errPlt = np.mean(rand_cmpnts_vmr[gasVer][p],axis=0)
                retPrf = np.mean(rPrfVMR[gasVer],axis=0)
            else:
                errPlt = rand_cmpnts_vmr[gasVer][p][0]
                retPrf = rPrfVMR[gasVer][0]

            errPlt = errPlt / retPrf *100.

            if p == 'measurement_Random_Error': legend = 'Measurement'
            elif p == 'interfering_species_Random_Error': legend = 'Interfering species'
            elif p == 'temperature_Random_Error': legend = 'Temperature'
            elif p == 'sza_Random_Error': legend = 'SZA'
            elif p == 'retrieval_parameters_Random_Error': legend = 'Retrieval parameters'
            elif p == 'curvature_Random_Error': legend = 'Curvature' #continue
            elif p == 'max_opd_Random_Error': legend = 'Max OPD' #continue     
            else: legend = p
        
            ax1[0].plot(errPlt,alt[gasVer],linewidth=2.0, label=legend)

        ax1[0].plot(np.nanmean(rand_errvmr[gasVer], axis=0)/ retPrf *100.,alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
        ax1[0].set_ylabel('Altitude [km]', fontsize=14)
        #ax1[0].set_xlabel('VMR [ppm]', fontsize=14)
        ax1[0].set_xlabel('% wrt mean', fontsize=14)             
        ax1[0].grid(True,which='both')
        ax1[0].legend(prop={'size':12},loc='upper right')          
        #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
        ax1[0].set_title('Random Components')
        if loc.lower() == 'fl0': ax1[0].set_ylim(1, 15.5)
        if loc.lower() == 'mlo': ax1[0].set_ylim(3, 20)
        ax1[0].tick_params(which='both',labelsize=14)
        ax1[0].set_xlim(0, 10)

        #-----------------------------
        #ERRORS:   Plot systematic error components
        #-----------------------------

        for p in sorted(sys_cmpnts_vmr[gasVer].iterkeys()):

            if npnts > 1:
                errPlt = np.mean(sys_cmpnts_vmr[gasVer][p],axis=0)
                retPrf = np.mean(rPrfVMR[gasVer],axis=0)
            else:
                errPlt = sys_cmpnts_vmr[gasVer][p][0]
                retPrf = rPrfVMR[gasVer][0]

            errPlt = errPlt / retPrf *100.
            
            if p == 'temperature_Systematic_Error': legend = 'Temperature'
            elif p == 'smoothing_Systematic_Error':  legend = 'Smoothing'
            #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'            
            elif p == 'linepair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$\gamma$'
            elif p == 'lineint_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$S$'
            elif p == 'linetair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$n$'
            elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn' #continue
            elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
            elif p == 'omega_Systematic_Error': legend = 'Omega' #continue
            else: legend = p

            ax1[1].plot(errPlt, alt[gasVer],linewidth=2.0, label=legend)

        ax1[1].plot(np.nanmean(sys_errvmr[gasVer], axis=0)/ retPrf *100., alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
        #ax1[1].plot(np.nanmean(tot_errvmr[gasVer], axis=0)/ retPrf *100., alt[gasVer],linestyle='--', color='r',linewidth=2.0, label='Total')
        #ax1[1].set_xlabel('VMR [ppm]', fontsize=14)
        ax1[1].set_xlabel('% wrt mean', fontsize=14)             
        ax1[1].grid(True,which='both')
        ax1[1].legend(prop={'size':12}, loc='upper right')
        ax1[1].set_title('Systematic Components')
        if loc.lower() == 'fl0': ax1[1].set_ylim(1, 15.5)
        if loc.lower() == 'mlo': ax1[1].set_ylim(3, 20)
        ax1[1].tick_params(which='both',labelsize=14)
        ax1[1].set_xlim(0, 20)
            #ax1[1].set_xscale('log')

        fig.subplots_adjust(bottom=0.1, top=0.90, left = 0.1, right = 0.95)

        plt.suptitle('Uncertainty of {} in %'.format(statDataCl[gasVer].PrimaryGas), fontsize=16)
        
        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
        else:           
            plt.show(block=False)


    #---------------------------------
    # Finding exactly the same dates
    #---------------------------------
    doy_h2o  = mf.toYearFraction(dates[gasName[0]+'_'+ver[0]])
    doy_hdo  = mf.toYearFraction(dates[gasName[1]+'_'+ver[1]])
    
    intrsctVals = np.intersect1d(doy_h2o, doy_hdo, assume_unique=True)

    if len(intrsctVals) > 1:
        print 'Number of identical dates = {}'.format(len(intrsctVals))
    else: 
        print 'Error: Dates of measurements are different'
        exit()

    inds1       = np.nonzero( np.in1d( doy_hdo, intrsctVals, assume_unique=False ) )[0]
    inds2       = np.nonzero( np.in1d( doy_h2o, intrsctVals, assume_unique=False ) )[0]

    doy_hdo     = doy_hdo[inds1]
    doy_h2o     = doy_h2o[inds2]

    #---------------------------------
    # Defining variables (FTIR) - H2O
    #---------------------------------

    h2over        = gasName[0]+'_'+ver[0]
    
    indsalt       = np.where(alt[h2over] <= maxalt)[0]

    alt_h2o         = alt[h2over][indsalt]

    indsaoi        = mf.nearestind(aoi, alt_h2o)
    aoi            = alt_h2o[indsaoi]

    TC_h2o          = totClmn[h2over][inds2]
    Dates_h2o       = dates[h2over][inds2]
   
    Prf_h2o         = rPrfVMR[h2over][:, indsalt];                      Prf_h2o = Prf_h2o[inds2, :]

    avkSCF_h2o      = avkSCF[h2over][inds2, indsalt[0]:, indsalt[0]:]
    avkSCav_h2o     = np.mean(avkSCF_h2o, axis=0)
    avkVMR_h2o      = avkVMR[h2over][inds2, indsalt[0]:, indsalt[0]:]
    avkVMRav_h2o     = np.mean(avkVMR_h2o, axis=0)
    aPrf_h2o        = aPrfVMR[h2over][:, indsalt];                       aPrf_h2o    = aPrf_h2o[inds2, :]         
    Airmass_h2o     = Airmass[h2over][:, indsalt];                       Airmass_h2o = Airmass_h2o[inds2, :]                      
    rPrfMol_h2o     = rPrfMol[h2over][:, indsalt];                       rPrfMol_h2o = rPrfMol_h2o[inds2, :]
    aPrfMol_h2o     = aPrfMol[h2over][:, indsalt];                       aPrfMol_h2o = aPrfMol_h2o[inds2, :]


    TCpc_h2o        = np.sum(rPrfMol_h2o,axis=1)
    VMR_ns_h2o      = Prf_h2o[:, indsaoi]
    VMRa_ns_h2o     = aPrf_h2o[:, indsaoi] 

    #Prf_hdoW        = rPrfVMR_2[:, indsalt]*3.107e-4;                   Prf_hdoW = Prf_hdoW[inds2, :]
    #VMR_ns_hdoW     = rPrfVMR_2[inds2, -1]
    #Prf_dDW         = (np.divide(Prf_hdoW, Prf_h2o)/3.1152e-4 - 1) *1000.
    #dDftsW          = Prf_dDW[-1]

    PresPrf_h2o     = PresPrf[h2over][:, indsalt];                      PresPrf_h2o = PresPrf_h2o[inds2, :]

    avkTC_h2o       = avkTC[h2over][inds2, :] 


    if errorFlg: 
        ranErr_h2o     = rand_errvmr[h2over][:,indsalt];                     ranErr_h2o = ranErr_h2o[inds2, :]
        sysErr_h2o     = sys_errvmr[h2over][:, indsalt];                     sysErr_h2o = sysErr_h2o[inds2, :]
        totErr_h2o     = tot_errvmr[h2over][:, indsalt];                     totErr_h2o = totErr_h2o[inds2, :]

        ranErrMol_h2o     = rand_err[h2over][:,indsalt];                     ranErrMol_h2o = ranErrMol_h2o[inds2, :]
        sysErrMol_h2o     = sys_err[h2over][:, indsalt];                     sysErrMol_h2o = sysErrMol_h2o[inds2, :]
        totErrMol_h2o     = tot_err[h2over][:, indsalt];                     totErrMol_h2o = totErrMol_h2o[inds2, :]

        TC_h2o_e        = totClmn_e[h2over][inds2]
        TC_h2o_R_e      = tot_rnd[h2over][inds2]
        #tot_err_h2o     = tot_err[h2over][:, indsalt]; tot_err_h2o = tot_err_h2o[inds2, :]
        #VMR_e_ns_h2o    = tot_errvmr[h2over][inds2, indsaoi]
        VMR_e_ns_h2o    = totErr_h2o[:, indsaoi] #tot_errvmr[h2over][inds2, indsaoi]

    rms_h2o          = rms[h2over][inds2]
    sza_h2o          = sza[h2over][inds2]
    saa_h2o          = saa[h2over][inds2]


    #---------------------------------
    # Defining variables (FTIR) - HDO
    #---------------------------------
    hdover        = gasName[1]+'_'+ver[1]
    hours         = [d.hour for d in dates[hdover]]
    
    TC_hdo          = totClmn[hdover][inds1]*3.107e-4
    Dates_hdo       = dates[hdover][inds1]
    Prf_hdo         = rPrfVMR[hdover][:, indsalt]*3.107e-4 ;             Prf_hdo = Prf_hdo[inds1, :]             
    avkSCF_hdo      = avkSCF[hdover][inds1, indsalt[0]:, indsalt[0]:]
    avkSCav_hdo     = np.mean(avkSCF_hdo, axis=0)
    avkVMR_hdo      = avkVMR[hdover][inds1, indsalt[0]:, indsalt[0]:]
    avkVMRav_hdo    = np.mean(avkVMR_hdo, axis=0)
    aPrf_hdo        = aPrfVMR[hdover][:, indsalt]*3.107e-4 ;            aPrf_hdo = aPrf_hdo[inds1, :]
    Airmass_hdo     = Airmass[hdover][:, indsalt];                      Airmass_hdo = Airmass_hdo[inds1, :]
    rPrfMol_hdo     = rPrfMol[hdover][:, indsalt]*3.107e-4;             rPrfMol_hdo = rPrfMol_hdo[inds1, :]
    aPrfMol_hdo     = aPrfMol[hdover][:, indsalt]*3.107e-4;             aPrfMol_hdo = aPrfMol_hdo[inds1, :]

    TCpc_hdo        = np.sum(rPrfMol_hdo, axis=1)
    VMR_ns_hdo      = Prf_hdo[:, indsaoi]   
    VMRa_ns_hdo     = aPrf_hdo[:, indsaoi]

    avkTC_hdo       = avkTC[hdover][inds1, :] 

    if errorFlg: 
        ranErr_hdo     = rand_errvmr[hdover][:, indsalt]*3.107e-4;          ranErr_hdo = ranErr_hdo[inds1, :]
        sysErr_hdo     = sys_errvmr[hdover][:, indsalt]*3.107e-4;           sysErr_hdo = sysErr_hdo[inds1, :]
        totErr_hdo     = tot_errvmr[hdover][:, indsalt]*3.107e-4;           totErr_hdo = totErr_hdo [inds1, :]

        ranErrMol_hdo     = rand_err[hdover][:,indsalt]*3.107e-4;                     ranErrMol_hdo = ranErrMol_hdo[inds1, :]
        sysErrMol_hdo     = sys_err[hdover][:, indsalt]*3.107e-4;                     sysErrMol_hdo = sysErrMol_hdo[inds1, :]
        totErrMol_hdo     = tot_err[hdover][:, indsalt]*3.107e-4;                     totErrMol_hdo = totErrMol_hdo[inds1, :]

        TC_hdo_e         = totClmn_e[hdover][inds1]*3.107e-4
        TC_hdo_R_e       = tot_rnd[hdover][inds1]*3.107e-4
        #tot_err_hdo     = tot_err[hdover][:, indsalt]*3.107e-4; tot_err_hdo = tot_err_hdo[inds1, :]
        #VMR_e_ns_hdo    = tot_errvmr[hdover][inds1, indsaoi]*3.107e-4

        VMR_e_ns_hdo    = totErr_hdo[:, indsaoi]

    rms_hdo          = rms[hdover][inds1]
    sza_hdo          = sza[hdover][inds1]
    saa_hdo          = saa[hdover][inds1]
    
    #---------------------------------
    # Smoothing --> WE CAN'T DO THIS BECAUSE HDO & H2O ARE DIFFERENT GASES
    #---------------------------------
    # Prf_h2o_s       = []
    # Prf_hdo_s       = []
    # rPrfMol_hdo_s   = []

    # ###TC_hdo_s  = []

    # for itime in range(len(Prf_h2o)):
    #     Prf_h2o_s.append(aPrf_hdo[itime, :]/3.107e-4 + np.dot(avkVMR_hdo[itime, :, :], (Prf_h2o[itime, :] -  aPrf_hdo[itime, :]/3.107e-4)) )
    #     Prf_hdo_s.append(aPrf_h2o[itime, :]*3.107e-4 + np.dot(avkVMR_h2o[itime, :, :], (Prf_hdo[itime, :] -  aPrf_h2o[itime, :]*3.107e-4)) )
    #     ##Prf_hdo_s.append(aPrf_h2o[itime, :]*3.107e-4 + np.dot(avkVMRav_h2o[:, :], (Prf_hdo[itime, :] -  aPrf_h2o[itime, :]*3.107e-4)) )

    #     rPrfMol_hdo_s.append(aPrfMol_h2o[itime, :]*3.107e-4 + np.dot(avkSCF_h2o[itime, :, :], (rPrfMol_hdo[itime, :] -  aPrfMol_h2o[itime, :]*3.107e-4)) )

    #     TC_ap  = np.sum(aPrfMol_h2o[itime,:])

    #     ## print avkTC_h2o[itime, :].shape
    #     ## print TC_ap.shape
    #     ## print TC_hdo[itime].shape
    #     ## exit()

    #     ## print TC_ap*3.107e-4
    #     ## print avkTC_h2o[itime, :]
    #     ## print TC_hdo[itime]
    #     ## print np.dot(avkTC_h2o[itime, :], (TC_hdo[itime] -  TC_ap*3.107e-4))[0]
    #     ## print np.sum(avkTC_h2o[itime, :]*(TC_hdo[itime] -  TC_ap*3.107e-4))
    #     ## print TC_hdo[itime] -  TC_ap*3.107e-4

    #     ## exit()

    #     ## TC_hdo_s.append(TC_ap*3.107e-4+ np.sum(abs(np.dot(avkTC_h2o[itime, :], (TC_hdo[itime] -  TC_ap*3.107e-4)))) )
    #     ## TC_hdo_s.append(TC_ap*3.107e-4+ np.sum(np.dot(avkTC_h2o[itime, :], (TC_hdo[itime] -  TC_ap*3.107e-4))) )

    #     ## Prf_h2o_s.append(aPrf_hdo[itime, :]/3.107e-4 + np.dot(avkVMRav_hdo, (Prf_h2o[itime, :] -  aPrf_hdo[itime, :]/3.107e-4)) )
    #     ## Prf_hdo_s.append(aPrf_h2o[itime, :]*3.107e-4 + np.dot(avkVMRav_h2o, (Prf_hdo[itime, :] -  aPrf_h2o[itime, :]*3.107e-4)) )

    # Prf_h2o_s         = np.asarray(Prf_h2o_s)
    # Prf_hdo_s         = np.asarray(Prf_hdo_s)
    # rPrfMol_hdo_s     = np.asarray(rPrfMol_hdo_s)

    # VMR_ns_h2o_s      = Prf_h2o_s[:,indsaoi]
    # VMR_ns_hdo_s      = Prf_hdo_s[:,indsaoi]

    # TC_hdo_s          = np.sum(rPrfMol_hdo_s, axis=1)      
    # TC_hdo_s          = np.asarray(TC_hdo_s)

    #---------------------------------
    # Defining dD (FTIR)
    #---------------------------------
    Prf_dD      = (np.divide(Prf_hdo, Prf_h2o)/3.1152e-4 - 1.0) *1000.
    #TC_dD       = np.sum(Prf_dD, axis=1)
    TC_dD       = (np.divide(TC_hdo, TC_h2o)/3.1152e-4 - 1.0) *1000.
    dDfts       = (np.divide(VMR_ns_hdo, VMR_ns_h2o)/3.1152e-4 - 1.0) *1000.
    dDfts2      = Prf_dD[:, indsaoi]

    #Prf_dD_s      = (np.divide(Prf_hdo, Prf_h2o_s)/3.1152e-4 - 1.0) *1000.
    #Prf_dD_s      = (np.divide(Prf_hdo_s, Prf_h2o)/3.1152e-4 - 1.0) *1000.
    #TC_dD_s       = np.sum(Prf_dD_s, axis=1)
    #TC_dD_s       = (np.divide(TC_hdo_s, TC_h2o)/3.1152e-4 - 1.0) *1000.
    #dDfts_s       = (np.divide(VMR_ns_hdo_s, VMR_ns_h2o)/3.1152e-4 - 1.0) *1000.

    #---------------------------------
    #Error propagation of dD
    #---------------------------------
    dDfts_RanErr =  np.abs(Prf_dD)*(np.sqrt( (ranErr_hdo/Prf_hdo)**2 + (ranErr_h2o/Prf_h2o)**2 ))
    dDfts_SysErr =  np.abs(Prf_dD)*(np.sqrt( (sysErr_hdo/Prf_hdo)**2 + (sysErr_h2o/Prf_h2o)**2 ))
    dDfts_TotErr =  np.abs(Prf_dD)*(np.sqrt( (totErr_hdo/Prf_hdo)**2 + (totErr_h2o/Prf_h2o)**2 ))

    if errorFlg: 
        dDfts_e           = dDfts_TotErr[:,indsaoi]
        TC_dD_e           = np.abs(TC_dD)*(np.sqrt( (TC_hdo_e/TC_hdo)**2 + (TC_h2o_e/TC_h2o)**2 ))
        TC_dD_R_e         = np.abs(TC_dD)*(np.sqrt( (TC_hdo_R_e/TC_hdo)**2 + (TC_h2o_R_e/TC_h2o)**2 ))

    #---------------------------------
    #Second Filter: positive dD profiles
    #---------------------------------
    indsFilter = []

    rprf_pos  = np.asarray(Prf_dD) >= 100
    indsT     = np.where( np.sum(rprf_pos, axis=1) > 100 )[0]
    print ('Total number observations found with positive partial column = {}'.format(len(indsT)))
    print ('Percent of observations found with positive partial column = {0:.2f}'.format(float(len(indsT))/float(len(Dates_h2o)) * 100.))
    indsFilter = np.union1d(indsT, indsFilter)

    #---------------------------------
    #Filter H2O
    #---------------------------------
    TC_h2o           = np.delete(TC_h2o, indsFilter, axis=0)
    TC_h2o_e         = np.delete(TC_h2o_e, indsFilter, axis=0)
    TC_h2o_R_e       = np.delete(TC_h2o_R_e, indsFilter, axis=0)
    Dates_h2o        = np.delete(Dates_h2o, indsFilter, axis=0)
    Prf_h2o          = np.delete(Prf_h2o, indsFilter, axis=0)
    aPrf_h2o         = np.delete(aPrf_h2o, indsFilter, axis=0)
    ranErr_h2o       = np.delete(ranErr_h2o, indsFilter, axis=0)
    sysErr_h2o       = np.delete(sysErr_h2o, indsFilter, axis=0)
    totErr_h2o       = np.delete(totErr_h2o, indsFilter, axis=0)

    rPrfMol_h2o      = np.delete(rPrfMol_h2o, indsFilter, axis=0)

    ranErrMol_h2o     = np.delete(ranErrMol_h2o, indsFilter, axis=0)
    sysErrMol_h2o     = np.delete(sysErrMol_h2o, indsFilter, axis=0)
    totErrMol_h2o     = np.delete(totErrMol_h2o, indsFilter, axis=0)

    TCpc_h2o         = np.delete(TCpc_h2o, indsFilter, axis=0)
    VMR_ns_h2o       = np.delete(VMR_ns_h2o, indsFilter, axis=0)
    VMRa_ns_h2o      = np.delete(VMRa_ns_h2o, indsFilter, axis=0)
    #Prf_hdoW         = np.delete(Prf_hdoW, indsFilter, axis=0)
    #VMR_ns_hdoW      = np.delete(VMR_ns_hdoW, indsFilter, axis=0)
    #Prf_dDW          = np.delete(Prf_dDW, indsFilter, axis=0)
    #TCe_h2o          = np.delete(TCe_h2o, indsFilter, axis=0)
    #tot_err_h2o      = np.delete(tot_err_h2o, indsFilter, axis=0)
    VMR_e_ns_h2o     = np.delete(VMR_e_ns_h2o, indsFilter, axis=0)
    avkVMR_h2o       = np.delete(avkVMR_h2o, indsFilter, axis=0)
    PresPrf_h2o      = np.delete(PresPrf_h2o, indsFilter, axis=0)

    Airmass_h2o      = np.delete(Airmass_h2o, indsFilter, axis=0)

    #Prf_h2o_s        = np.delete(Prf_h2o_s, indsFilter, axis=0)
    #VMR_ns_h2o_s     = np.delete(VMR_ns_h2o_s, indsFilter, axis=0)


    #---------------------------------
    #Filter HDO
    #---------------------------------
    TC_hdo           = np.delete(TC_hdo, indsFilter, axis=0)
    TC_hdo_e         = np.delete(TC_hdo_e, indsFilter, axis=0)
    TC_hdo_R_e       = np.delete(TC_hdo_R_e, indsFilter, axis=0)
    Dates_hdo        = np.delete(Dates_hdo, indsFilter, axis=0)
    Prf_hdo          = np.delete(Prf_hdo, indsFilter, axis=0)
    aPrf_hdo         = np.delete(aPrf_hdo, indsFilter, axis=0)
    ranErr_hdo       = np.delete(ranErr_hdo, indsFilter, axis=0)
    sysErr_hdo       = np.delete(sysErr_hdo, indsFilter, axis=0)
    totErr_hdo       = np.delete(totErr_hdo, indsFilter, axis=0)

    rPrfMol_hdo      = np.delete(rPrfMol_hdo, indsFilter, axis=0)

    ranErrMol_hdo     = np.delete(ranErrMol_hdo, indsFilter, axis=0)
    sysErrMol_hdo     = np.delete(sysErrMol_hdo, indsFilter, axis=0)
    totErrMol_hdo     = np.delete(totErrMol_hdo, indsFilter, axis=0)


    TCpc_hdo         = np.delete(TCpc_hdo, indsFilter, axis=0)
    VMR_ns_hdo       = np.delete(VMR_ns_hdo, indsFilter, axis=0)
    VMRa_ns_hdo      = np.delete(VMRa_ns_hdo, indsFilter, axis=0)
    #tot_err_hdo      = np.delete(tot_err_hdo, indsFilter, axis=0)
    VMR_e_ns_hdo     = np.delete(VMR_e_ns_hdo, indsFilter, axis=0)
    avkVMR_hdo       = np.delete(avkVMR_hdo, indsFilter, axis=0)

    #Prf_hdo_s          = np.delete(Prf_hdo_s, indsFilter, axis=0)
    #VMR_ns_hdo_s     = np.delete(VMR_ns_hdo_s, indsFilter, axis=0)

    #---------------------------------
    #Filter dD
    #---------------------------------
    Prf_dD           = np.delete(Prf_dD, indsFilter, axis=0)
    dDfts_RanErr     = np.delete(dDfts_RanErr, indsFilter, axis=0)
    dDfts_SysErr     = np.delete(dDfts_SysErr, indsFilter, axis=0)
    dDfts_TotErr     = np.delete(dDfts_TotErr, indsFilter, axis=0)

    dDfts_e          = np.delete(dDfts_e, indsFilter, axis=0)
    dDfts            = np.delete(dDfts, indsFilter, axis=0)
 
    doy_h2o          = mf.toYearFraction(Dates_h2o)
    doy_hdo          = mf.toYearFraction(Dates_hdo)

    #Prf_dD_s         = np.delete(Prf_dD_s, indsFilter, axis=0)
    #dDfts_s          = np.delete(dDfts_s, indsFilter, axis=0)

    TC_dD            = np.delete(TC_dD, indsFilter)
    #TC_dD_s          = np.delete(TC_dD_s, indsFilter)
    TC_dD_e          = np.delete(TC_dD_e, indsFilter)
    TC_dD_R_e        = np.delete(TC_dD_R_e, indsFilter)


    #-------------------------------------------------
    # Calculate weighted VMR H2O, HDO, and dD
    #-------------------------------------------------
    vmrP_hdo     = {}
    vmrP_hdo_e   = {}
    vmrP_h2o     = {}
    vmrP_h2o_e   = {}
    #vmrP_h2o_s   = {}
    #vmrP_hdo_s   = {}
    dDP          = {}
    dDP_e        = {}

    dDP_2        = {}
    dDP_e_2      = {}
    #dDP_s        = {}

    PresP        = {}
    altw         = {}

    tcP_h2o      = {}
    tcP_h2o_e    = {}
    tcP_h2o_R_e = {}
    
    tcP_hdo      = {}
    tcP_hdo_e    = {}
    tcP_hdo_R_e = {}

    #tcP_hdo_s    = {}

    dDP_tcp      = {}
    dDP_tcp_e    = {}
    dDP_tcp_R_e  = {}
    
    #dDP_tcp_s    = {}
    
    for p, pcol in enumerate(pCols):
        
        #ind1             = mf.nearestind(pcol[0], alt_h2o)
        #ind2             = mf.nearestind(pcol[1], alt_h2o)

        inds = np.where( (alt_h2o >= pcol[0]) & (alt_h2o <=pcol[1])  )[0]

        vmrP_h2o.setdefault(str(p) ,[]).append(np.average(Prf_h2o[:,inds],axis=1, weights=Airmass_h2o[:,inds])) 
        vmrP_h2o_e.setdefault(str(p) ,[]).append(np.average(totErr_h2o[:,inds],axis=1, weights=Airmass_h2o[:,inds])) 

        vmrP_hdo.setdefault(str(p) ,[]).append(np.average(Prf_hdo[:,inds],axis=1,weights=Airmass_h2o[:,inds]))
        vmrP_hdo_e.setdefault(str(p) ,[]).append(np.average(totErr_hdo[:,inds],axis=1,weights=Airmass_h2o[:,inds]))

        dDP_2.setdefault(str(p) ,[]).append(np.average(Prf_dD[:,inds],axis=1,weights=Airmass_h2o[:,inds]))
        dDP_e_2.setdefault(str(p) ,[]).append(np.average(dDfts_TotErr[:,inds],axis=1,weights=Airmass_h2o[:,inds]))

        #vmrP_h2o_s.setdefault(str(p) ,[]).append(np.average(Prf_h2o_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        #vmrP_hdo_s.setdefault(str(p) ,[]).append(np.average(Prf_hdo_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        
        #dDP_s.setdefault(str(p) ,[]).append(np.average(Prf_dD_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1])) 

        PresP.setdefault(str(p) ,[]).append(np.average(PresPrf_h2o[:,inds],axis=1, weights=Airmass_h2o[:,inds]))

        altw.setdefault(str(p) ,[]).append(np.average(alt_h2o[inds], weights=np.mean(Airmass_h2o[:, inds], axis=0)))

        tcP_h2o.setdefault(str(p) ,[]).append(np.sum(rPrfMol_h2o[:,inds],axis=1)) 
        tcP_h2o_e.setdefault(str(p) ,[]).append(np.sum(totErrMol_h2o[:,inds],axis=1))
        tcP_h2o_R_e.setdefault(str(p) ,[]).append(np.sum(ranErrMol_h2o[:,inds],axis=1)) 

        tcP_hdo.setdefault(str(p) ,[]).append(np.sum(rPrfMol_hdo[:,inds],axis=1))
        tcP_hdo_e.setdefault(str(p) ,[]).append(np.sum(totErrMol_hdo[:,inds],axis=1)) 
        tcP_hdo_R_e.setdefault(str(p) ,[]).append(np.sum(ranErrMol_hdo[:,inds],axis=1)) 
        #tcP_hdo_s.setdefault(str(p) ,[]).append(np.sum(rPrfMol_hdo_s[:,ind2:ind1],axis=1))



        # vmrP_h2o.setdefault(str(p) ,[]).append(np.average(Prf_h2o[:,ind2:ind1],axis=1, weights=Airmass_h2o[:,ind2:ind1])) 
        # vmrP_h2o_e.setdefault(str(p) ,[]).append(np.average(totErr_h2o[:,ind2:ind1],axis=1, weights=Airmass_h2o[:,ind2:ind1])) 

        # vmrP_hdo.setdefault(str(p) ,[]).append(np.average(Prf_hdo[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        # vmrP_hdo_e.setdefault(str(p) ,[]).append(np.average(totErr_hdo[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))

        # dDP_2.setdefault(str(p) ,[]).append(np.average(Prf_dD[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        # dDP_e_2.setdefault(str(p) ,[]).append(np.average(dDfts_TotErr[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))

        # #vmrP_h2o_s.setdefault(str(p) ,[]).append(np.average(Prf_h2o_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        # #vmrP_hdo_s.setdefault(str(p) ,[]).append(np.average(Prf_hdo_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1]))
        
        # #dDP_s.setdefault(str(p) ,[]).append(np.average(Prf_dD_s[:,ind2:ind1],axis=1,weights=Airmass_h2o[:,ind2:ind1])) 

        # PresP.setdefault(str(p) ,[]).append(np.average(PresPrf_h2o[:,ind2:ind1],axis=1, weights=Airmass_h2o[:,ind2:ind1]))

        # altw.setdefault(str(p) ,[]).append(np.average(alt_h2o[ind2:ind1], weights=np.mean(Airmass_h2o[:, ind2:ind1], axis=0)))

        # tcP_h2o.setdefault(str(p) ,[]).append(np.sum(rPrfMol_h2o[:,ind2:ind1],axis=1)) 
        # tcP_h2o_e.setdefault(str(p) ,[]).append(np.sum(totErrMol_h2o[:,ind2:ind1],axis=1)) 

        # tcP_hdo.setdefault(str(p) ,[]).append(np.sum(rPrfMol_hdo[:,ind2:ind1],axis=1))
        # tcP_hdo_e.setdefault(str(p) ,[]).append(np.sum(totErrMol_hdo[:,ind2:ind1],axis=1)) 
        # #tcP_hdo_s.setdefault(str(p) ,[]).append(np.sum(rPrfMol_hdo_s[:,ind2:ind1],axis=1)) 



    for p, pcol in enumerate(pCols):
        vmrP_h2o[str(p)] = np.asarray(vmrP_h2o[str(p)])
        vmrP_hdo[str(p)] = np.asarray(vmrP_hdo[str(p)])
       
        vmrP_h2o_e[str(p)] = np.asarray(vmrP_h2o_e[str(p)])
        vmrP_hdo_e[str(p)] = np.asarray(vmrP_hdo_e[str(p)])

        dDP_2[str(p)]      = np.asarray(dDP_2[str(p)])
        dDP_e_2[str(p)]      = np.asarray(dDP_e_2[str(p)])

        dDP.setdefault(str(p) ,[]).append((np.divide(vmrP_hdo[str(p)][0], vmrP_h2o[str(p)][0])/3.1152e-4 - 1.0) *1000.)
        dDP_e.setdefault(str(p) ,[]).append(np.abs(dDP[str(p)][0])*(np.sqrt( (vmrP_hdo_e[str(p)][0]/vmrP_hdo[str(p)][0])**2 + (vmrP_h2o_e[str(p)][0]/vmrP_h2o[str(p)][0])**2 ))    )

        dDP[str(p)]        = np.asarray(dDP[str(p)])
        dDP_e[str(p)]      = np.asarray(dDP_e[str(p)])

        #vmrP_h2o_s[str(p)] = np.asarray(vmrP_h2o_s[str(p)])
        #vmrP_hdo_s[str(p)] = np.asarray(vmrP_hdo_s[str(p)])
        

        #dDP_s.setdefault(str(p) ,[]).append((np.divide(vmrP_hdo_s[str(p)][0], vmrP_h2o[str(p)][0])/3.1152e-4 - 1.0) *1000.)
        ###dDP_s.setdefault(str(p) ,[]).append((np.divide(vmrP_hdo[str(p)][0], vmrP_h2o_s[str(p)][0])/3.1152e-4 - 1.0) *1000.)

        #dDP_s[str(p)]      = np.asarray(dDP_s[str(p)])


        PresP[str(p)]      = np.asarray(PresP[str(p)]) 

        altw[str(p)]      = np.asarray(altw[str(p)])


        tcP_h2o[str(p)]       = np.asarray(tcP_h2o[str(p)])
        tcP_hdo[str(p)]       = np.asarray(tcP_hdo[str(p)])
        tcP_h2o_e[str(p)]     = np.asarray(tcP_h2o_e[str(p)])
        tcP_hdo_e[str(p)]     = np.asarray(tcP_hdo_e[str(p)])
        tcP_h2o_R_e[str(p)]   = np.asarray(tcP_h2o_R_e[str(p)])

        #tcP_hdo_s[str(p)] = np.asarray(tcP_hdo_s[str(p)])


        dDP_tcp.setdefault(str(p) ,[]).append((np.divide(tcP_hdo[str(p)][0], tcP_h2o[str(p)][0])/3.1152e-4 - 1.0) *1000.)
        dDP_tcp_e.setdefault(str(p) ,[]).append(np.abs(dDP_tcp[str(p)][0])*(np.sqrt( (tcP_hdo_e[str(p)][0]/tcP_hdo[str(p)][0])**2 + (tcP_h2o_e[str(p)][0]/tcP_h2o[str(p)][0])**2 ))    ) 
        dDP_tcp_R_e.setdefault(str(p) ,[]).append(np.abs(dDP_tcp[str(p)][0])*(np.sqrt( (tcP_hdo_R_e[str(p)][0]/tcP_hdo[str(p)][0])**2 + (tcP_h2o_R_e[str(p)][0]/tcP_h2o[str(p)][0])**2 ))    )   

        #dDP_tcp_s.setdefault(str(p) ,[]).append((np.divide(tcP_hdo_s[str(p)][0], tcP_h2o[str(p)][0])/3.1152e-4 - 1.0) *1000.)

        dDP_tcp[str(p)]      = np.asarray(dDP_tcp[str(p)])
        dDP_tcp_e[str(p)]    = np.asarray(dDP_tcp_e[str(p)])
        dDP_tcp_R_e[str(p)]    = np.asarray(dDP_tcp_R_e[str(p)])
        #dDP_tcp_s[str(p)]    = np.asarray(dDP_tcp_s[str(p)])


        #print dDP_e[str(p)].shape
        #print dDP_s[str(p)].shape

        #exit()

 

    #---------------------------------
    #Error plots of dD - FTIR
    #---------------------------------
    fig, ax1 = plt.subplots(1, 2, figsize=(8,7), sharey=True, sharex=False)

    ax1[0].plot(np.mean(dDfts_RanErr, axis=0), alt_h2o, linewidth=2.0, label='Random Uncertainty')
    ax1[0].plot(np.mean(dDfts_SysErr, axis=0), alt_h2o, linewidth=2.0, label='Systematic Uncertainty')
    #ax1[0].plot(np.mean(dDfts_TotErr, axis=0), alt_h2o, linewidth=2.0, label='Total Uncertainty')
    ax1[0].plot(np.mean(dDfts_TotErr, axis=0), alt_h2o,linestyle='--', color='k',linewidth=2.0, label='Total Uncertainty')

    ax1[1].plot(np.mean(dDfts_RanErr, axis=0)/np.mean(Prf_dD, axis=0) * 100., alt_h2o, linewidth=2.0, label='Random Uncertainty')
    ax1[1].plot(np.mean(dDfts_SysErr, axis=0)/np.mean(Prf_dD, axis=0) * 100., alt_h2o, linewidth=2.0, label='Systematic Uncertainty')
    #ax1[1].plot(np.mean(dDfts_TotErr, axis=0)/np.mean(Prf_dD, axis=0) * 100., alt_h2o, linewidth=2.0, label='Total Uncertainty')
    ax1[1].plot(np.mean(dDfts_TotErr, axis=0)/np.mean(Prf_dD, axis=0) * 100., alt_h2o,linestyle='--', color='k',linewidth=2.0, label='Total Uncertainty')

    
    ax1[0].set_ylabel('Altitude [km]', fontsize=14)
    ax1[0].set_xlabel('Uncertainty in dD', fontsize=14)
    #ax1.set_xlabel('Error [%]', fontsize=14)             
    ax1[0].grid(True,which='both')
    ax1[0].legend(prop={'size':12},loc='upper right')          
    #ax1[0].set_title('Uncertainty [ppm]')
    ax1[0].set_ylim(3, 20)
    ax1[0].tick_params(which='both',labelsize=14)

    ax1[1].grid(True,which='both')
    ax1[1].set_xlabel('% wrt mean]', fontsize=14)        
    #ax1[1].set_title('Uncertainty [%]')
    ax1[1].set_ylim(3.0, 20)
    ax1[1].tick_params(which='both',labelsize=14)
    
    fig.subplots_adjust(bottom=0.1, top=0.90, left = 0.1, right = 0.95)

    plt.suptitle('Uncertainty of dD', fontsize=16)
        
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()   


    #---------------------------------
    # Interpolation of in-situ to FTS
    #---------------------------------
    if insFlg:
        doyin_i            = mf.toYearFraction(Date_in)

        hours              =  [float(d.hour) for d in Dates_h2o]
        hours              = np.asarray(hours)

        insitu_interp      = interpolate.interp1d(doyin_i, insitu, bounds_error=False, fill_value=(insitu[-1], insitu[0]))(doy_hdo)
        HDOinsitu_interp   = interpolate.interp1d(doyin_i, HDOinsitu, bounds_error=False, fill_value=(HDOinsitu[-1], HDOinsitu[0]))(doy_hdo)
        dDinsitu_r_interp  = interpolate.interp1d(doyin_i, dDinsitu_r, bounds_error=False, fill_value=(dDinsitu_r[-1], dDinsitu_r[0]))(doy_hdo)
        dDinsitu_interp    = interpolate.interp1d(doyin_i, dDinsitu,bounds_error=False, fill_value=(dDinsitu[-1], dDinsitu[0]))(doy_hdo)

        #insitu_interp      = interpolate.interp1d(doyin_i, insitu, bounds_error=False)(doy_hdo)
        #HDOinsitu_interp   = interpolate.interp1d(doyin_i, HDOinsitu, bounds_error=False)(doy_hdo)
        #dDinsitu_r_interp  = interpolate.interp1d(doyin_i, dDinsitu_r, bounds_error=False)(doy_hdo)
        #dDinsitu_interp    = interpolate.interp1d(doyin_i, dDinsitu,bounds_error=False)(doy_hdo)

        insitu_interp_s    = []
        HDOinsitu_interp_s = []
        #dDinsitu_interp_s  = []    

        ##Smoothing using FTS AK
        for itime in range(len(insitu_interp)):
            insitu_interp_s.append(aPrf_h2o[itime, -1] + np.dot(avkVMR_h2o[itime, -1, -1], (insitu_interp[itime] -  aPrf_h2o[itime, -1])) )
            HDOinsitu_interp_s.append(aPrf_hdo[itime, -1] + np.dot(avkVMR_hdo[itime, -1, -1], (HDOinsitu_interp[itime] -  aPrf_hdo[itime, -1])) )

            avkVMR_dD = np.dot(avkVMR_h2o[itime,:,:], avkVMR_hdo[itime,:,:])

            #dDinsitu_interp_s.append(aPrf_h2o[itime, -1] + np.dot(avkVMR_dD[itime, -1, -1], (dDinsitu_interp[itime] -  aPrf_h2o[itime, -1])) )

        insitu_interp_s    = np.asarray(insitu_interp_s)
        HDOinsitu_interp_s = np.asarray(HDOinsitu_interp_s)
        #dDinsitu_interp_s = np.asarray(dDinsitu_interp_s)

        dDinsitu_interp_s  = (np.divide(HDOinsitu_interp_s, insitu_interp_s)/3.1152e-4 - 1.0) *1000.

    #---------------------------------
    #plots
    #---------------------------------
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    dayLc        = DayLocator()
    yearsLc      = YearLocator()
    monthLc      = MonthLocator()
    mondays      = WeekdayLocator(MONDAY)
    DateFmt      = DateFormatter('%b %d')

    #--------------------------------
    # FIGURE: H2O Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_h2o])
    fig,(ax1,ax2)  = plt.subplots(1,2,sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    #ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    #ax2.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )

    ax1.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    ax2.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    
    for i in range(len(month)):
        ax1.plot(Prf_h2o[i,:]/1e3,alt_h2o,linewidth=0.75)
        ax2.plot(Prf_h2o[i,:]/1e3,alt_h2o,linewidth=0.75)

        #ax1.plot(aPrf_h2o[i,:],alt_h2o,linewidth=0.75, alpha=0.5)
        #ax2.plot(aPrf_h2o[i,:],alt_h2o,linewidth=0.75, alpha=0.5)

        #ax1.plot(rPrfMol_h2o[i,:],alt_h2o,linewidth=0.75)
        #ax2.plot(rPrfMol_h2o[i,:],alt_h2o,linewidth=0.75)

        #ax1.plot(aPrfMol_h2o[i,:],alt_h2o,linewidth=0.75, alpha=0.5)
        #ax2.plot(aPrfMol_h2o[i,:],alt_h2o,linewidth=0.75, alpha=0.5)
    
    #ax1.plot(np.mean(aPrf_h2o/1e3, axis=0),alt_h2o,'k--',linewidth=4,label='A priori')
    ###ax1.fill_betweenx(alt_h2o,np.mean(aPrf_h2o/1e3, axis=0)-np.std(aPrf_h2o/1e3, axis=0), np.mean(aPrf_h2o/1e3, axis=0)+np.std(aPrf_h2o/1e3, axis=0),alpha=0.5,color='0.75')

    #ax2.plot(np.mean(aPrf_h2o/1e3, axis=0),alt_h2o,'k--',linewidth=4,label='A priori')
    ###ax2.fill_betweenx(alt_h2o,np.mean(aPrf_h2o/1e3, axis=0)-np.std(aPrf_h2o/1e3, axis=0), np.mean(aPrf_h2o/1e3, axis=0)+np.std(aPrf_h2o/1e3, axis=0),alpha=0.5,color='0.75')

    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('VMR [ppm, x10$^3$]', fontsize=14)
    ax2.set_xlabel('Log VMR [ppm, x10$^3$]', fontsize=14)
    
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax2.set_xscale('log')
    
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    ax1.legend(prop={'size':12})
    ax2.legend(prop={'size':12})

    ax2.set_ylim(3.0, 20)
    
    ax1.tick_params(axis='x',which='both',labelsize=14)
    ax2.tick_params(axis='x',which='both',labelsize=14)  
    plt.suptitle('FTIR Profiles: {}'.format('H$_2$O'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()

    #---------------------------------
    # Plot : Averaging Kernel: H2O
    #---------------------------------
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    fig       = plt.figure(figsize=(9,7))
    gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
    ax        = plt.subplot(gs[0])
    axb       = plt.subplot(gs[1])
    axc       = plt.subplot(gs[2])
    cm        = plt.get_cmap(clmap)
    cNorm     = colors.Normalize(vmin=np.min(alt_h2o), vmax=np.max(alt_h2o))
    scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
    scalarMap.set_array(alt_h2o)
    
    #---------------------------------
    #ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_h2o])
    ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in alt_h2o] ) )
    
    for i in range(len(alt_h2o)):
        ax.plot(avkSCav_h2o[i,:], alt_h2o)
        #ax.plot(avkVMRav_h2o [i,:], alt_h2o)
        
    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_xlabel('Averaging Kernels', fontsize=14)
    ax.grid(True)

    #cbaxes = fig.add_axes([0.52, 0.3, 0.03, 0.5]) 
    cbaxes = fig.add_axes([0.4, 0.5, 0.02, 0.35]) 
    cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
    cbar.set_label('Altitude [km]', fontsize=14)
    ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(1, 15)
    if loc.lower() == 'fl0': ax.set_ylim(1, 15)
    if loc.lower() == 'mlo': ax.set_ylim(3, 20)

    #---------------------------------
    axb.plot(np.sum(avkSCav_h2o,axis=0), alt_h2o,color='k')
    #axb.plot(np.sum(avkVMRav_h2o,axis=0), alt_h2o,color='k')
    axb.grid(True)
    axb.set_xlabel('AK Area', fontsize=14)
    axb.tick_params(axis='x',which='both',labelsize=14)
    major_ticks = np.arange(0, 3, 1)
    axb.set_xticks(major_ticks) 
    #axb.set_ylim(1, 15)
    if loc.lower() == 'fl0': axb.set_ylim(1, 15)
    if loc.lower() == 'mlo': axb.set_ylim(3, 20)
    #axb.tick_params(labelsize=14) 

    #---------------------------------
    dofs_cs = np.cumsum(np.diag(avkSCav_h2o)[::-1])[::-1]
    #dofs_cs = np.cumsum(np.diag(avkVMRav_h2o)[::-1])[::-1]
    axc.plot(dofs_cs,alt_h2o,color='k',label='Cumulative Sum of DOFS (starting at surface)')
    xval = range(0,int(np.ceil(max(dofs_cs)))+2)
    axc.set_xlabel('Cumulative\nSum of DOFS', fontsize=14)  
    axc.tick_params(labelsize=14)
    #axc.set_title('DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol), fontsize=9)    
    axc.set_ylim((0,120))
    axc.grid(True,which='both',  alpha=0.5)

    major_ticks = np.arange(0, 3, 1)
    axc.set_xticks(major_ticks)

    if loc.lower() == 'fl0': axc.set_ylim(1, 15)
    if loc.lower() == 'mlo': axc.set_ylim(3, 20)
    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else:       plt.show(block=False)


    #--------------------------------
    # FIGURE: HDO Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_h2o])
    fig,(ax1,ax2)  = plt.subplots(1,2,sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    #ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    #ax2.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )

    ax1.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    ax2.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    
    for i in range(len(month)):
        ax1.plot(Prf_hdo[i,:], alt_h2o,linewidth=0.75)
        ax2.plot(Prf_hdo[i,:], alt_h2o,linewidth=0.75)

        # ax1.plot(aPrf_hdo[i,:],alt_h2o,linewidth=0.75, alpha=0.5)
        # ax2.plot(aPrf_hdo[i,:],alt_h2o,linewidth=0.75, alpha=0.5)

        # ax1.plot(Prf_hdo_s[i,:],alt_h2o,linewidth=0.75, color='r')
        # ax2.plot(Prf_hdo_s[i,:],alt_h2o,linewidth=0.75, color='r')

        

        #ax1.plot(rPrfMol_hdo[i,:],alt_h2o,linewidth=0.75)
        #ax2.plot(rPrfMol_hdo[i,:],alt_h2o,linewidth=0.75)

        #ax1.plot(aPrfMol_hdo[i,:],alt_h2o,linewidth=0.75, alpha=0.5)
        #ax2.plot(aPrfMol_hdo[i,:],alt_h2o,linewidth=0.75, alpha=0.5)

        #ax1.plot(rPrfMol_hdo_s[i,:],alt_h2o,linewidth=0.75, color='r')
        #ax2.plot(rPrfMol_hdo_s[i,:],alt_h2o,linewidth=0.75, color='r')
    
    #ax1.plot(np.mean(aPrf_hdo, axis=0),alt_h2o,'k--',linewidth=4,label='A priori')
    ##ax1.fill_betweenx(alt_h2o,np.mean(aPrf_hdo, axis=0)-np.std(aPrf_hdo, axis=0), np.mean(aPrf_hdo, axis=0)+np.std(aPrf_hdo, axis=0),alpha=0.5,color='0.75')

    #ax2.plot(np.mean(aPrf_hdo, axis=0),alt_h2o,'k--',linewidth=4,label='A priori')
    ##ax2.fill_betweenx(alt_h2o,np.mean(aPrf_hdo, axis=0)-np.std(aPrf_hdo, axis=0), np.mean(aPrf_hdo, axis=0)+np.std(aPrf_hdo, axis=0),alpha=0.5,color='0.75')
    
    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('VMR [ppm]', fontsize=14)
    ax2.set_xlabel('Log VMR [ppm]', fontsize=14)
    
    ax1.grid(True,which='both')
    ax2.grid(True,which='both')
    ax2.set_xscale('log')
    
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    ax1.legend(prop={'size':12})
    ax2.legend(prop={'size':12})

    ax2.set_ylim(3.0, 20)
    
    ax1.tick_params(axis='x',which='both',labelsize=14)
    ax2.tick_params(axis='x',which='both',labelsize=14)  
    plt.suptitle('FTIR Profiles: {}'.format('HDO'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # Plot : Averaging Kernel: HDO
    #---------------------------------
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)
    fig       = plt.figure(figsize=(9,7))
    gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
    ax        = plt.subplot(gs[0])
    axb       = plt.subplot(gs[1])
    axc       = plt.subplot(gs[2])
    cm        = plt.get_cmap(clmap)
    cNorm     = colors.Normalize(vmin=np.min(alt_h2o), vmax=np.max(alt_h2o))
    scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
    scalarMap.set_array(alt_h2o)
    
    #---------------------------------
    #ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_h2o])
    ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in alt_h2o] ) )
    
    for i in range(len(alt_h2o)):
        ax.plot(avkSCav_hdo[i,:], alt_h2o)
        #ax.plot(avkVMRav_hdo[i,:], alt_h2o)
        

    ax.set_ylabel('Altitude [km]', fontsize=14)
    ax.set_xlabel('Averaging Kernels', fontsize=14)
    ax.grid(True)

    #cbaxes = fig.add_axes([0.52, 0.3, 0.03, 0.5]) 
    cbaxes = fig.add_axes([0.4, 0.5, 0.02, 0.35]) 
    cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
    cbar.set_label('Altitude [km]', fontsize=14)
    ax.set_title('HDO Averaging Kernels Scale Factor', fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(1, 15)
    if loc.lower() == 'fl0': ax.set_ylim(1, 15)
    if loc.lower() == 'mlo': ax.set_ylim(3, 20)

    #---------------------------------
    axb.plot(np.sum(avkSCav_hdo,axis=0), alt_h2o,color='k')
    #axb.plot(np.sum(avkVMRav_hdo,axis=0), alt_h2o,color='k')
    axb.grid(True)
    axb.set_xlabel('AK Area', fontsize=14)
    axb.tick_params(axis='x',which='both',labelsize=14)
    major_ticks = np.arange(0, 3, 1)
    axb.set_xticks(major_ticks) 
    #axb.set_ylim(1, 15)
    if loc.lower() == 'fl0': axb.set_ylim(1, 15)
    if loc.lower() == 'mlo': axb.set_ylim(3, 20)
    #axb.tick_params(labelsize=14) 

    #---------------------------------
    dofs_cs = np.cumsum(np.diag(avkSCav_hdo)[::-1])[::-1]
    #dofs_cs = np.cumsum(np.diag(avkVMRav_hdo)[::-1])[::-1]
    axc.plot(dofs_cs,alt_h2o,color='k',label='Cumulative Sum of DOFS (starting at surface)')
    xval = range(0,int(np.ceil(max(dofs_cs)))+2)
    axc.set_xlabel('Cumulative\nSum of DOFS', fontsize=14)  
    axc.tick_params(labelsize=14)
    #axc.set_title('DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol), fontsize=9)    
    axc.set_ylim((0,120))
    axc.grid(True,which='both',  alpha=0.5)

    major_ticks = np.arange(0, 3, 1)
    axc.set_xticks(major_ticks)

    if loc.lower() == 'fl0': axc.set_ylim(1, 15)
    if loc.lower() == 'mlo': axc.set_ylim(3, 20)

    if saveFlg: pdfsav.savefig(fig,dpi=200)
    else:       plt.show(block=False)

    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()


    #--------------------------------
    # FIGURE: dD Profiles as a function of Month
    #--------------------------------
    month = np.array([d.month for d in Dates_h2o])
    fig,(ax1)  = plt.subplots(1, sharey=True)
    cm             = plt.get_cmap(clmap)
    cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    scalarMap.set_array(month)
    
    #ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    ax1.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    
   
    for i in range(len(month)):
        ax1.plot(Prf_dD[i,:],alt_h2o,linewidth=0.75)
    
    #ax1.plot(aprPrf[gas],alt_h2o,'k--',linewidth=4,label='A priori')
    #ax1.fill_betweenx(alt_h2o,np.mean(aPrf_hdo/1e3, axis=0)-np.std(aPrf_hdo/1e3, axis=0), np.mean(aPrf_hdo/1e3, axis=0)+np.std(aPrf_hdo/1e3, axis=0),alpha=0.5,color='0.75')
    
    #ax2.plot(aprPrf[gas],alt,'k--',linewidth=4,label='A priori')
    #ax1.fill_betweenx(alt_h2o,np.mean(aPrf_hdo/1e3, axis=0)-np.std(aPrf_hdo/1e3, axis=0), np.mean(aPrf_hdo/1e3, axis=0)+np.std(aPrf_hdo/1e3, axis=0),alpha=0.5,color='0.75')
    
    
    ax1.set_ylabel('Altitude [km]', fontsize=14)
    ax1.set_xlabel('dD [per mil]', fontsize=14)
    ax1.grid(True,which='both')
    ax1.set_xlim(-1000, 50)
    ax1.set_ylim(3, 20)
   
    cbar = fig.colorbar(scalarMap,orientation='vertical')
    cbar.set_label('Month', fontsize=14)
    
    #ax1.legend(prop={'size':12})
   
    ax1.tick_params(axis='x',which='both',labelsize=14)
   
    plt.suptitle('FTIR Profiles: {}'.format('dD'), fontsize=16)
    
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
        

    # #--------------------------------
    # # FIGURE: dD Profiles as a function of Month
    # #--------------------------------
    # month = np.array([d.month for d in Dates_h2o])
    # fig,(ax1)  = plt.subplots(1, sharey=True)
    # cm             = plt.get_cmap(clmap)
    # cNorm          = colors.Normalize( vmin=np.min(month), vmax=np.max(month) )
    # scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
    
    # scalarMap.set_array(month)
    
    # #ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    # ax1.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    
    # for i in range(len(month)):
    #     ax1.plot(Prf_dD_s[i,:],alt_h2o,linewidth=0.75)
    
    # #ax1.plot(aprPrf[gas],alt_h2o,'k--',linewidth=4,label='A priori')
    # #ax1.fill_betweenx(alt_h2o,np.mean(aPrf_hdo/1e3, axis=0)-np.std(aPrf_hdo/1e3, axis=0), np.mean(aPrf_hdo/1e3, axis=0)+np.std(aPrf_hdo/1e3, axis=0),alpha=0.5,color='0.75')
    
    # #ax2.plot(aprPrf[gas],alt,'k--',linewidth=4,label='A priori')
    # #ax1.fill_betweenx(alt_h2o,np.mean(aPrf_hdo/1e3, axis=0)-np.std(aPrf_hdo/1e3, axis=0), np.mean(aPrf_hdo/1e3, axis=0)+np.std(aPrf_hdo/1e3, axis=0),alpha=0.5,color='0.75')
    
    # ax1.set_ylabel('Altitude [km]', fontsize=14)
    # ax1.set_xlabel('dD [per mil]', fontsize=14)
    # ax1.grid(True,which='both')
    # ax1.set_xlim(-1000, 0)
    # ax1.set_ylim(0, 20)
    
    # cbar = fig.colorbar(scalarMap,orientation='vertical')
    # cbar.set_label('Month', fontsize=14)
    
    # ax1.legend(prop={'size':12})
    
    # ax1.tick_params(axis='x',which='both',labelsize=14)
    
    # plt.suptitle('FTIR Profiles: {} - smoothed'.format('dD'), fontsize=16)
    
    # if saveFlg:     
    #     pdfsav.savefig(fig,dpi=200)
    # else:           
    #     plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit() 

    # #---------------------------------
    # #correlation of dD vd dD_2 AT SEVERAL LAYERS
    # #---------------------------------
    npanels      = int(math.ceil(len(pCols)/3.0))


    # if len(pCols) <= 4: 
    #     outer_grid = gridspec.GridSpec(1, len(pCols),  wspace=0.25, hspace=0.35)
    #     fig = plt.figure(figsize=(12, 5))
    #     fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.92)
    # else: 
    #     outer_grid = gridspec.GridSpec(3, npanels, wspace=0.2, hspace=0.175)
    #     fig = plt.figure(figsize=(18, 12))
    #     fig.subplots_adjust(left=0.1, bottom=0.075, right=0.95, top=0.95)

    # for p, pcol in enumerate(pCols):

    #     ax = plt.Subplot(fig, outer_grid[p])
        
    #     ax.scatter(dDP[str(p)], dDP_2[str(p)], facecolors='white', s=40, color='red')
    #     #ax.scatter(dDP[str(p)], dDP_tcp[str(p)], facecolors='white', s=40, color='green')
              
    #     if p ==0: ax.set_ylabel('dD - 2',fontsize=14)
    #     ax.set_xlabel('dD - original',fontsize=14)
    #     ax.tick_params(labelsize=14)
    #     ax.grid(True)
    #     ax.annotate('{} - {} km'.format(pcol[0], pcol[1]), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=16, ha='left')
    #     #ax.set_xscale('log')
    #     #ax.set_xlim(0, 16)
    #     #ax.set_ylim(0, 16)
    #     #if p ==0: ax.legend(prop={'size':10})

    #     fig.add_subplot(ax)

    #     plt.suptitle('Correlation of weighted dD', fontsize=16)

    #     slope, intercept, r_value, p_value, std_err = stats.linregress(dDP[str(p)], dDP_2[str(p)])

    #     odr, odrErr  = mf.orthoregress(dDP[str(p)], dDP_2[str(p)], xerr=dDP_e[str(p)], yerr=dDP_e_2[str(p)], InError=True)
    #     slope2       = float(odr[0])
    #     intercept2   = float(odr[1])

    #     ax.text(0.025,0.78,"Slope: {0:.2f}".format(slope2),transform=ax.transAxes,  fontsize=14, color='red')
    #     ax.text(0.025,0.71,"Intercept: {:.3f}".format(intercept2),transform=ax.transAxes,  fontsize=14, color='red')
    #     ax.text(0.025,0.64,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='red')

    #     #bias = mf.bias(dDP[str(p)], dDP_s[str(p)])
    #     #ax.text(0.025,0.57,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='r')

    #     o2o = np.arange(np.min(dDP[str(p)]), 0)

    #     ax.plot(o2o, o2o, color= 'k', linestyle='--')
    #     ax.plot(o2o*0.8, o2o, color= 'k', linestyle='--')
    #     ax.plot(o2o, o2o*0.8, color= 'k', linestyle='--')


    # if saveFlg:     
    #     pdfsav.savefig(fig,dpi=200)
    # else:           
    #     plt.show(block=False)

    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit() 
        

    #---------------------------------
    #correlation of dD vd dD smoothed at single altitude
    #---------------------------------

    # fig, ax = plt.subplots(figsize=(8,6))

    # #ax.scatter(Prf_h2o_s[:, -1], Prf_h2o[:, -1], facecolors='white', s=40, color='b')
    # ax.scatter(dDfts, dDfts_s, facecolors='white', s=40, color='r')               
    # ax.set_ylabel('dD - smoothed',fontsize=14)
    # ax.set_xlabel('dD - original',fontsize=14)
    # ax.tick_params(labelsize=14)
    # ax.grid(True)
    # #ax.set_xlim(-850, 400)
    # #ax.set_ylim(-850, 400)
    # ax.set_ylim(top=0)
    # ax.set_xlim(right=0)

    # plt.suptitle('Correlation of weighted mean dD (original vs smoothed)\nAltitude = {0:.1f}'.format(aoi), fontsize=16)

    # slope, intercept, r_value, p_value, std_err = stats.linregress(dDfts, dDfts_s)

    # odr, odrErr  = mf.orthoregress(dDfts, dDfts_s, xerr=dDfts_e, yerr=dDfts_e, InError=True)
    # slope2       = float(odr[0])
    # intercept2   = float(odr[1])

    # ax.text(0.025,0.78,"Slope: {0:.2f}".format(slope2),transform=ax.transAxes,  fontsize=14, color='red')
    # ax.text(0.025,0.71,"Intercept: {:.3f}".format(intercept2),transform=ax.transAxes,  fontsize=14, color='red')
    # ax.text(0.025,0.64,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='red')

    # bias = mf.bias(dDfts, dDfts_s)
    # ax.text(0.025,0.57,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='r')

    # o2o = np.arange(np.min(dDfts), 0)

    # ax.plot(o2o, o2o, color= 'k', linestyle='--')
    # ax.plot(o2o*0.8, o2o, color= 'k', linestyle='--')
    # ax.plot(o2o, o2o*0.8, color= 'k', linestyle='--')

    # if saveFlg:     
    #     pdfsav.savefig(fig,dpi=200)
    # else:           
    #     plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit() 

    # #---------------------------------
    # # FIGURE - TIME SERIES OF H2O
    # #---------------------------------
    # fig, ax = plt.subplots(figsize=(12,6))

    # #insitu    = np.asarray(insitu)
    # #Date_in   = np.asarray(Date_in)
    
    # #ax.scatter(Date_in, insitu/1e3, facecolors='white', s=40, color='b',alpha=0.75, label='in-situ (measured)')

    # if insFlg:
    #     ax.scatter(Dates_hdo, insitu_interp/1e3, facecolors='white', s=40, color='b', label='in-situ (measured)')
    #     ax.scatter(Dates_hdo, insitu_interp_s/1e3, facecolors='white', s=40, color='green', label='in-situ (smoothed by FTIR AK)')
    
    # if errorFlg: ax.errorbar(Dates_h2o,VMR_ns_h2o/1e3, yerr=VMR_e_ns_h2o/1e3, fmt='o', markersize=0, color='r', ecolor='r')
    # ax.scatter(Dates_h2o, VMR_ns_h2o/1e3, facecolors='white', s=40, color='r', label='Retrieved FTIR - near surface')
    # #ax.scatter(Dates_hdo[inds], VMRa_ns_h2o_i[inds]/1e3, facecolors='white', s=40, color='gray', label='a priori (ERA) - near surface')
    
    # ax.grid(True)
    # ax.set_ylabel('VMR [ppm, x10$^3$]',fontsize=14)
    # ax.tick_params(labelsize=14)
    # ax.xaxis.set_minor_locator(dayLc)
    # ax.xaxis.set_major_formatter(DateFmt)
    # ax.legend(prop={'size':14})

    # plt.suptitle('Time series of H2O', fontsize=16)

    # fig.subplots_adjust(left=0.1, right=0.95)

    # if saveFlg:     
    #     pdfsav.savefig(fig,dpi=200)
    # else:           
    #     plt.show(block=False)
            

    # #---------------------------------
    # # FIGURE - TIME SERIES OF HDO
    # #---------------------------------
    # fig, ax = plt.subplots(figsize=(12,6))

    # if insFlg:
    #     ax.scatter(Dates_hdo, HDOinsitu_interp,   facecolors='white', s=40, color='b', label='in-situ (Calculated)')
    #     ax.scatter(Dates_hdo, HDOinsitu_interp_s, facecolors='white', s=40, color='green', label='in-situ (smoothed by FTIR AK)')
    
    # if errorFlg: ax.errorbar(Dates_hdo,VMR_ns_hdo, yerr=VMR_e_ns_hdo, fmt='o', markersize=0, color='r', ecolor='r')
    # ax.scatter(Dates_hdo, VMR_ns_hdo, facecolors='white', s=40, color='r', label='Retrieved FTIR')
    # #ax.scatter(Dates_hdo, VMRa_ns_hdo, facecolors='white', s=40, color='gray', label='a priori (ERA)')

    
    # ax.grid(True)
    # ax.set_ylabel('VMR [ppm]',fontsize=14)
    # ax.tick_params(labelsize=14)
    # ax.xaxis.set_minor_locator(dayLc)
    # ax.xaxis.set_major_formatter(DateFmt)
    # ax.legend(prop={'size':14})

    # plt.suptitle('Time series of HDO', fontsize=16)

    # fig.subplots_adjust(left=0.1, right=0.95)

    # if saveFlg:     
    #     pdfsav.savefig(fig,dpi=200)
    # else:           
    #     plt.show(block=False)
        

    # #---------------------------------
    # # FIGURE - TIME SERIES OF dD
    # #---------------------------------

    fig, ax = plt.subplots(figsize=(12,6))

   
    if errorFlg: 
        ax.errorbar(Dates_h2o, TC_dD, yerr=TC_dD_e,  fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')
        #ax.errorbar(Dates_hdo,dDfts_s, yerr=dDfts_e,  fmt='o', ecolor='r', mec='r', mfc='w', mew=1, label='FTIR-smoothed')
    else:
        ax.scatter(Dates_h2o, TC_dD, facecolors='white', s=40, color='g', label='FTIR')
        #ax.scatter(Dates_hdo, dDfts_s, facecolors='white', s=40, color='blue', label='FTIR-smoothed')
   
    ax.grid(True)
    ax.set_ylabel('dD [per mil]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.set_ylim(-750, 0)
    ax.xaxis.set_minor_locator(dayLc)
    ax.xaxis.set_major_formatter(DateFmt)
    #ax.legend(prop={'size':14})

    #plt.suptitle('Time series of dD', fontsize=16)
    plt.suptitle('Time series of dD (total column'.format(aoi), fontsize=16)

    fig.subplots_adjust(left=0.1, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)


    #---------------------------------
    #dD vs Water vapor (several partial columns) - Random Error
    #---------------------------------
    if len(pCols) <= 4: 
        outer_grid = gridspec.GridSpec(1, len(pCols),  wspace=0.25, hspace=0.35)
        fig = plt.figure(figsize=(12, 5))
        fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.92)
    else: 
        outer_grid = gridspec.GridSpec(3, npanels, wspace=0.2, hspace=0.175)
        fig = plt.figure(figsize=(18, 12))
        fig.subplots_adjust(left=0.1, bottom=0.075, right=0.95, top=0.95)


    for p, pcol in enumerate(pCols):

        ax = plt.Subplot(fig, outer_grid[p])

        ax.errorbar(tcP_h2o[str(p)][0,:]*2.989e-22, dDP_tcp[str(p)][0,:], yerr=dDP_tcp_R_e[str(p)][0, :], xerr=tcP_h2o_R_e[str(p)][0,:]*2.989e-22,  fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')
      
              
        if p == 0: ax.set_ylabel('dD [per mil]',fontsize=14)
        #ax.set_xlabel('H$_2$O VMR [ppm, x10$^3$]',fontsize=14)
        ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
        ax.tick_params(labelsize=14)
        ax.grid(True)
        ax.annotate('{} - {} km'.format(pcol[0], pcol[1]), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=16, ha='left')
        ax.set_ylim(top=0)
        #ax.set_xscale('log')
        #ax.set_xlim(0, 16)
        ax.set_ylim(-1000, 50)
        #if p ==0: ax.legend(prop={'size':10})

        fig.add_subplot(ax)

        plt.suptitle('[dD, H2O] - Partial Columns (Random errors)', fontsize=16)
        
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)


        
    #---------------------------------
    #dD vs Water vapor (several partial columns) - Total Error
    #---------------------------------
    #outer_grid = gridspec.GridSpec(2, 2, wspace=0.2, hspace=0.175)

    if len(pCols) <= 4: 
        outer_grid = gridspec.GridSpec(1, len(pCols),  wspace=0.25, hspace=0.35)
        fig = plt.figure(figsize=(12, 5))
        fig.subplots_adjust(left=0.1, bottom=0.15, right=0.95, top=0.92)
    else: 
        outer_grid = gridspec.GridSpec(3, npanels, wspace=0.2, hspace=0.175)
        fig = plt.figure(figsize=(18, 12))
        fig.subplots_adjust(left=0.1, bottom=0.075, right=0.95, top=0.95)


    for p, pcol in enumerate(pCols):

        ax = plt.Subplot(fig, outer_grid[p])
        
        #ax.errorbar(vmrP_h2o[str(p)][0,:]/1e3, dDP[str(p)][0,:], yerr=dDP_e[str(p)][0, :], xerr=vmrP_h2o_e[str(p)][0,:]/1e3, fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR - 1')
        #ax.errorbar(vmrP_h2o[str(p)][0,:]/1e3, dDP_2[str(p)][0,:], yerr=dDP_e_2[str(p)][0, :], xerr=vmrP_h2o_e[str(p)][0,:]/1e3, fmt='o', ecolor='r', mec='r', mfc='w', mew=1, label='FTIR - 2')
        ###ax.scatter(vmrP_h2o[str(p)]/1e3, dDP[str(p)], facecolors='white', s=40, color='green', label='Original')

        #ax.errorbar(vmrP_h2o[str(p)][0,:]/1e3, dDP_s[str(p)][0,:], yerr=dDP_e[str(p)][0,:], xerr=vmrP_h2o_e[str(p)][0,:]/1e3, fmt='o', ecolor='r', mec='r', mfc='w', mew=1, label='FTIR-smoothed')
        ###ax.scatter(vmrP_h2o[str(p)]/1e3, dDP_s[str(p)], facecolors='white', s=40, color='r',  label='smoothed')
        ###ax.scatter(insitu_interp_s/1e3, VMRa_ns_h2o_i/1e3, facecolors='white', s=40, color='gray', label='A priori (ERA) Vs in-situ (smooth by FTIR AK)')

        ax.errorbar(tcP_h2o[str(p)][0,:]*2.989e-22, dDP_tcp[str(p)][0,:], yerr=dDP_tcp_e[str(p)][0, :], xerr=tcP_h2o_e[str(p)][0,:]*2.989e-22,  fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')
      
              
        if p == 0: ax.set_ylabel('dD [per mil]',fontsize=14)
        #ax.set_xlabel('H$_2$O VMR [ppm, x10$^3$]',fontsize=14)
        ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
        ax.tick_params(labelsize=14)
        ax.grid(True)
        ax.annotate('{} - {} km'.format(pcol[0], pcol[1]), xy=(0.025, 0.9), xycoords='axes fraction', fontsize=16, ha='left')
        ax.set_ylim(top=0)
        #ax.set_xscale('log')
        #ax.set_xlim(0, 16)
        ax.set_ylim(-1000, 50)
        #if p ==0: ax.legend(prop={'size':10})

        fig.add_subplot(ax)

        plt.suptitle('[dD, H2O] - Partial Columns (Total error)', fontsize=16)
        
    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

     #---------------------------------
    # FIGURE - dD vs Water vapor (Total Columns) - Random Error
    #---------------------------------

    fig, ax = plt.subplots(figsize=(8,6))
    
    if errorFlg: 
        ax.errorbar(TC_h2o*2.989e-22,TC_dD, yerr=TC_dD_R_e, xerr= TC_h2o_R_e*2.989e-22,  fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')
        #ax.errorbar(TC_h2o,TC_dD_s, yerr=TC_dD_e, xerr= TC_h2o_e,  fmt='o', ecolor='r', mec='r', mfc='w', mew=1, label='FTIR-smoothed')
    else:
        ax.scatter(TC_h2o*2.989e-22,TC_dD, facecolors='white', s=40, color='r', label='FTIR')
    
    ax.grid(True)
    ax.set_ylabel('dD [total column]',fontsize=14)
    ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(-750, 0)
    #ax.xaxis.set_minor_locator(dayLc)
    #ax.xaxis.set_major_formatter(DateFmt)
    #ax.legend(prop={'size':14})
    ax.set_ylim(-700, -100)

    plt.suptitle('[dD, H2O] - Total Column (Random error)', fontsize=16)

    fig.subplots_adjust(left=0.15, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #---------------------------------
    # FIGURE - dD vs Water vapor (Total Columns) - Total Error
    #--------------------------------- 

    fig, ax = plt.subplots(figsize=(8,6))
    
    #ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
    ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in month] ) )
    
   
    if errorFlg: 
            ax.errorbar(TC_h2o*2.989e-22,TC_dD, yerr=TC_dD_e, xerr= TC_h2o_e*2.989e-22,  fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')

    else:
        ax.scatter(TC_h2o*2.989e-22,TC_dD, facecolors='white', s=40, color='r', label='FTIR')
    
    ax.grid(True)
    ax.set_ylabel('dD [total column]',fontsize=14)
    ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(-750, 0)
    #ax.xaxis.set_minor_locator(dayLc)
    #ax.xaxis.set_major_formatter(DateFmt)
    #ax.legend(prop={'size':14})
    ax.set_ylim(-700, -100)


    plt.suptitle('[dD, H2O] - Total Column (Total error)', fontsize=16)

    fig.subplots_adjust(left=0.15, right=0.95)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()

    #---------------------------------
    # FIGURE - dD vs Water vapor (Total Columns) - Total Error
    #--------------------------------- 


    moi = [ [3, 4, 5], [6, 7, 8], [9, 10, 11], [1, 2, 12]  ]
    clr = ['r', 'g', 'orange', 'b']
    lbl = ['MAM', 'JJA', 'SON', 'DJF']

    fig, ax = plt.subplots(figsize=(8,6))
    
    if errorFlg: 

        for j, mi in enumerate(moi):

            inds = [i for i, m in enumerate(month) if m in mi]
          
            ax.scatter(TC_h2o[inds]*2.989e-22,TC_dD[inds], s=40, color=clr[j], label=lbl[j])
    else:
        ax.scatter(TC_h2o*2.989e-22,TC_dD, facecolors='white', s=40, color='r', label='FTIR')
    
    ax.grid(True)
    ax.set_ylabel('dD [total column]',fontsize=14)
    ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.legend(prop={'size':14})
    #ax.set_ylim(-750, 0)
    #ax.xaxis.set_minor_locator(dayLc)
    #ax.xaxis.set_major_formatter(DateFmt)
    ax.set_ylim(-700, -100)

    #cbar = fig.colorbar(scalarMap,orientation='vertical')
    #cbar.set_label('Month', fontsize=14)

    fig.subplots_adjust(left=0.15, right=0.95)

    plt.suptitle('[dD, H2O] - Total Column', fontsize=16)


    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)

    # #---------------------------------
    # # FIGURE - dD vs Water vapor (Total Columns) - Total Error
    # #--------------------------------- 

    tcks = range(np.min(month),np.max(month)+2)
    norm = colors.BoundaryNorm(tcks,cm.N)   

    fig, ax = plt.subplots(figsize=(8,6))
    
    if errorFlg: 
          
        sc = ax.scatter(TC_h2o*2.989e-22,TC_dD, s=40, c=month, cmap=cm, norm=norm)
    else:
        ax.scatter(TC_h2o*2.989e-22,TC_dD, facecolors='white', s=40, color='r', label='FTIR')
    
    ax.grid(True)
    ax.set_ylabel('dD [total column]',fontsize=14)
    ax.set_xlabel('H$_2$O Column [mm]',fontsize=14)
    ax.tick_params(labelsize=14)
    #ax.set_ylim(-750, 0)
    #ax.xaxis.set_minor_locator(dayLc)
    #ax.xaxis.set_major_formatter(DateFmt)
    ax.set_ylim(-700, -100)

    fig.subplots_adjust(left=0.15, right=0.82)

    cax  = fig.add_axes([0.86, 0.1, 0.03, 0.8])
    cbar = fig.colorbar(sc, cax=cax,format='%3i')
    cbar.set_label('Month')

    plt.suptitle('[dD, H2O] - Total Column', fontsize=16)


    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
        
        

    #---------------------------------
    #dD vs Water vapor near surface
    #---------------------------------
    # indLow = np.where(dDfts_s > -90.)[0]
    # print float(len(indLow))/float(len(dDfts_s))
    # print sza_h2o[indLow]
    # print saa_h2o[indLow]
    # print rms_h2o[indLow]
    # print rms_hdo[indLow]
    # print Dates_hdo[indLow]

    fig, ax = plt.subplots(figsize=(8,6))
    
    #ax.errorbar(vmrP_h2o[str(p)]/1e3, dDP[str(p)], yerr=VMR_e_ns_h2o/1e3, fmt='o',xerr=insitu_interp/1e3*0.05, markersize=0, color='green', ecolor='green')
    if insFlg: ax.plot(insitu/1e3, dDinsitu, 'k.',markersize=4, color='gray', alpha=0.5,  label='in-situ - All', zorder=1)
    if insFlg: ax.scatter(insitu_interp/1e3, dDinsitu_interp, facecolors='white', s=40, color='r', label='in-situ - interpolated to FTIR time', zorder=2)
    if insFlg: ax.scatter(insitu_interp_s/1e3, dDinsitu_interp, facecolors='white', s=40, color='blue', label='in-situ smoothed - interpolated to FTIR time', zorder=2)
    #ax.scatter(VMR_ns_h2o/1e3, dDfts, facecolors='white', s=40, color='green', label='FTIR', zorder=3)
    #ax.scatter(VMR_ns_h2o/1e3, dDfts_s, facecolors='white', s=40, color='k', label='FTIR-smoothed', zorder=3)

    ax.errorbar(VMR_ns_h2o/1e3, dDfts, yerr=dDfts_e, xerr=VMR_e_ns_h2o/1e3, fmt='o', ecolor='g', mec='g', mfc='w', mew=1, label='FTIR')
    #ax.errorbar(VMR_ns_h2o/1e3, dDfts_s, yerr=dDfts_e, xerr=VMR_e_ns_h2o/1e3,fmt='o', ecolor='r', mec='r', mfc='w', mew=1, label='FTIR-smoothed')
    
    ax.set_ylabel('dD [per mil]',fontsize=14)
    ax.set_xlabel('H$_2$O VMR [ppm, x10$^3$]',fontsize=14)
    ax.tick_params(labelsize=14)
    ax.grid(True)
    #ax.set_ylim(top=0)
    ax.set_ylim(-1000, 50)
    if insFlg: ax.legend(prop={'size':12}, loc=3)
    #ax.set_xscale('log')
    #ax.set_ylim(-700, 0)
    #ax.set_xlim(0, 20)


    plt.suptitle('[dD, H2O] pair at altitude = {0:.1f} km'.format(aoi), fontsize=16)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()

    #---------------------------------
    # FIGURE - CORRELATION (H2O)
    #---------------------------------
    if insFlg:
    
        fig, ax = plt.subplots(figsize=(8,6))

        ax.errorbar(insitu_interp/1e3, VMR_ns_h2o/1e3, yerr=VMR_e_ns_h2o/1e3, xerr=insitu_interp/1e3*0.05,  fmt='o', markersize=0, color='b', ecolor='b')
        ax.scatter(insitu_interp/1e3, VMR_ns_h2o/1e3, facecolors='white', s=40, color='b', label='FTIR Vs in-situ')
        
        ax.errorbar(insitu_interp_s/1e3, VMR_ns_h2o/1e3, yerr=VMR_e_ns_h2o/1e3, fmt='o',xerr=insitu_interp/1e3*0.05, markersize=0, color='green', ecolor='green')
        ax.scatter(insitu_interp_s/1e3, VMR_ns_h2o/1e3, facecolors='white', s=40, color='green', label='FTIR Vs in-situ (smoothed by FTIR AK)')
        #ax.scatter(insitu_interp_s/1e3, VMRa_ns_h2o_i/1e3, facecolors='white', s=40, color='gray', label='A priori (ERA) Vs in-situ (smooth by FTIR AK)')
              
        ax.set_ylabel('VMR [ppm, x10$^3$]',fontsize=14)
        ax.set_xlabel('VMR - in-situ [ppm, x10$^3$]',fontsize=14)
        ax.tick_params(labelsize=14)
        ax.grid(True)
        ax.set_xlim(0, 16)
        ax.set_ylim(0, 16)
        ax.legend(prop={'size':12})

        plt.suptitle('Orthogonal Linear Regression of H2O\n(near-surface FTIR Vs in-situ)', fontsize=16)

        slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp/1e3, VMR_ns_h2o/1e3)
       
        odr, odrErr  = mf.orthoregress(insitu_interp/1e3, VMR_ns_h2o/1e3, xerr=insitu_interp*0.05/1e3, yerr=VMR_e_ns_h2o/1e3, InError=True)
        slope       = float(odr[0])
        intercept   = float(odr[1])

        slope_e     = float(odrErr[0])
        intercept_e = float(odrErr[1])
       
        ax.text(0.025,0.92,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='blue')
        ax.text(0.025,0.85,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='blue')
        ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='blue')

        bias = mf.bias(insitu_interp/1e3, VMR_ns_h2o/1e3)
        ax.text(0.025,0.71,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='blue')

        slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp_s/1e3, VMR_ns_h2o/1e3)
      
        odr, odrErr  = mf.orthoregress(insitu_interp_s/1e3, VMR_ns_h2o/1e3, xerr=insitu_interp_s/1e3*0.05, yerr=VMR_e_ns_h2o/1e3, InError=True)
        slope       = float(odr[0])
        intercept   = float(odr[1])

        slope_e     = float(odrErr[0])
        intercept_e = float(odrErr[1])

        ax.text(0.025,0.63,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='green')
        ax.text(0.025,0.56,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='green')
        ax.text(0.025,0.49,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')
        
        bias = mf.bias(insitu_interp_s/1e3, VMR_ns_h2o/1e3)
        ax.text(0.025,0.42,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='green')

        rmse = mf.rmse(insitu_interp_s/1e3, VMR_ns_h2o/1e3)
     
        print 'bias (H2O) = {}'.format(bias) 
        print 'rmse (H2O) = {}'.format(rmse) 

        #slope, intercept, r_value, p_value, std_err = stats.linregress(insitu_interp_s, VMRa_ns_h2o_i)

        #ax.text(0.75,0.7,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='gray')
        #ax.text(0.75,0.63,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='gray')
        #ax.text(0.75,0.56,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='gray')

        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
        else:           
            plt.show(block=False)
            

         #---------------------------------
        # FIGURE - CORRELATION (HDO)
        #---------------------------------

        fig, ax = plt.subplots(figsize=(8,6))


        ax.errorbar(HDOinsitu_interp, VMR_ns_hdo, yerr=VMR_e_ns_hdo, xerr=HDOinsitu_interp*0.05,  fmt='o', markersize=0, color='b', ecolor='b')
        ax.scatter(HDOinsitu_interp, VMR_ns_hdo, facecolors='white', s=40, color='b', label='FTIR Vs in-situ')

        ax.errorbar(HDOinsitu_interp_s, VMR_ns_hdo, yerr=VMR_e_ns_hdo, xerr=HDOinsitu_interp*0.05,  fmt='o', markersize=0, color='green', ecolor='green')
        ax.scatter(HDOinsitu_interp_s, VMR_ns_hdo, facecolors='white', s=40, color='green', label='FTIR Vs in-situ (smoothed by FTIR AK)')

        #ax.errorbar(HDOinsitu_interp_s, VMR_ns_hdo_s, yerr=VMR_e_ns_hdo, xerr=HDOinsitu_interp*0.05,  fmt='o', markersize=0, color='r', ecolor='r')
        #ax.scatter(HDOinsitu_interp_s, VMR_ns_hdo_s, facecolors='white', s=40, color='r', label='FTIR Vs in-situ (smoothed dD and insitu)')
              
        ax.set_ylabel('VMR [ppm]',fontsize=14)
        ax.set_xlabel('VMR - in-situ [ppm]',fontsize=14)
        ax.tick_params(labelsize=14)
        ax.grid(True)
        ax.set_xlim(0, 4)
        ax.set_ylim(0, 4)
        ax.legend(prop={'size':12})

        plt.suptitle('Orthogonal Linear Regression of HDO\n(near-surface FTIR Vs in-situ)', fontsize=16)

        slope, intercept, r_value, p_value, std_err = stats.linregress(HDOinsitu_interp, VMR_ns_hdo)

        odr, odrErr  = mf.orthoregress(HDOinsitu_interp, VMR_ns_hdo, xerr=HDOinsitu_interp*0.05, yerr=VMR_e_ns_hdo, InError=True)
        slope       = float(odr[0])
        intercept   = float(odr[1])

        slope_e     = float(odrErr[0])
        intercept_e = float(odrErr[1])
       
        ax.text(0.025,0.92,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='blue')
        ax.text(0.025,0.85,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='blue')
        ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='blue')

        bias = mf.bias(HDOinsitu_interp, VMR_ns_hdo)
        ax.text(0.025,0.71,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='blue')

        slope, intercept, r_value, p_value, std_err = stats.linregress(HDOinsitu_interp_s, VMR_ns_hdo)

        odr, odrErr  = mf.orthoregress(HDOinsitu_interp_s, VMR_ns_hdo, xerr=HDOinsitu_interp_s*0.05, yerr=VMR_e_ns_hdo, InError=True)
        slope       = float(odr[0])
        intercept   = float(odr[1])

        slope_e     = float(odrErr[0])
        intercept_e = float(odrErr[1])

        ax.text(0.025,0.63,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='green')
        ax.text(0.025,0.56,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='green')
        ax.text(0.025,0.49,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')
        
        bias = mf.bias(HDOinsitu_interp_s, VMR_ns_hdo)
        ax.text(0.025,0.42,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='green')

        rmse = mf.rmse(HDOinsitu_interp_s, VMR_ns_hdo)
     
        print 'bias (HDO) = {}'.format(bias) 
        print 'rmse (HDO) = {}'.format(rmse) 

        #----
        # slope, intercept, r_value, p_value, std_err = stats.linregress(HDOinsitu_interp_s, VMR_ns_hdo_s)

        # odr, odrErr  = mf.orthoregress(HDOinsitu_interp_s, VMR_ns_hdo_s, xerr=HDOinsitu_interp_s*0.05, yerr=VMR_e_ns_hdo, InError=True)
        # slope       = float(odr[0])
        # intercept   = float(odr[1])

        # slope_e     = float(odrErr[0])
        # intercept_e = float(odrErr[1])

        # ax.text(0.025,0.34,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='r')
        # ax.text(0.025,0.27,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='r')
        # ax.text(0.025,0.21,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='r')
        
        # bias = mf.bias(HDOinsitu_interp_s, VMR_ns_hdo_s)
        # ax.text(0.025,0.14,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='r')

        # rmse = mf.rmse(HDOinsitu_interp_s, VMR_ns_hdo_s)
     
        # print 'bias (HDO) = {}'.format(bias) 
        # print 'rmse (HDO) = {}'.format(rmse) 

       
        if saveFlg:     
            pdfsav.savefig(fig,dpi=200)
        else:           
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        
        
        #---------------------------------
        # FIGURE - CORRELATION (dD)
        #---------------------------------
        # fig, ax = plt.subplots(figsize=(8,6))

        # #ax.scatter(dDinsitu_r_interp, dDfts, facecolors='white', s=40, color='b', label='in-situ (Calculated)')


        # ax.errorbar(dDinsitu_interp, dDfts, yerr=dDfts_e, xerr=dDinsitu_interp*0.05,  fmt='o', markersize=0, color='blue', ecolor='blue')
        # ax.scatter(dDinsitu_interp,   dDfts, facecolors='white', s=40, color='blue', label='FTIR Vs in-situ')
        
        # ax.errorbar(dDinsitu_interp_s, dDfts, yerr=dDfts_e, xerr=dDinsitu_interp_s*0.05,  fmt='o', markersize=0, color='green', ecolor='green')
        # ax.scatter(dDinsitu_interp_s, dDfts, facecolors='white', s=40, color='green',label='FTIR Vs in-situ (smoothed by FTIR AK)')

        # #ax.errorbar(dDinsitu_interp_s, dDfts_s, yerr=dDfts_e, xerr=dDinsitu_interp_s*0.05,  fmt='o', markersize=0, color='r', ecolor='r')
        # #ax.scatter(dDinsitu_interp_s, dDfts_s, facecolors='white', s=40, color='r',label='FTIR Vs in-situ (both smoothed)')
              
        # ax.set_ylabel('dD - FTIR [per mil]',fontsize=14)
        # ax.set_xlabel('dD - in-situ [per mil]',fontsize=14)
        # ax.tick_params(labelsize=14)
        # ax.grid(True)
        # #ax.set_xlim(-850, 400)
        # #ax.set_ylim(-850, 400)
        # ax.legend(prop={'size':14})

        # plt.suptitle('Orthogonal Linear Regression of dD\n(near-surface FTIR Vs in-situ)', fontsize=16)

        # slope, intercept, r_value, p_value, std_err = stats.linregress(dDinsitu_interp, dDfts)

        # odr, odrErr  = mf.orthoregress(dDinsitu_interp, dDfts, xerr=dDinsitu_interp*0.05, yerr=dDfts_e, InError=True)
        # slope       = float(odr[0])
        # intercept   = float(odr[1])

        # slope_e     = float(odrErr[0])
        # intercept_e = float(odrErr[1])
       
        # bias = mf.bias(dDinsitu_interp, dDfts)
        # rmse = mf.rmse(dDinsitu_interp, dDfts)
       
        # ax.text(0.025,0.92,"Slope: {0:.2f}".format(slope),transform=ax.transAxes,  fontsize=14, color='blue')
        # ax.text(0.025,0.85,"Intercept: {:.2f}".format(intercept/1e3),transform=ax.transAxes,  fontsize=14, color='blue')
        # ax.text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='blue')

        # ax.text(0.025,0.71,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='blue')


        # print 'bias (dD)= {}'.format(bias) 
        # print 'rmse (dD)= {}'.format(rmse)
        # weights      = np.ones_like(dDinsitu_interp)

        # # res    = mf.fit_driftfourier(dDinsitu_interp, dDfts, weights, 2)
        # # f_drift, f_fourier, f_driftfourier = res[3:6]

        # # print "Fitted slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(dDinsitu_interp)*100.0)
        # # print "Fitted intercept at xmin: {:.3E}".format(res[0])
        # # print "STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(dDinsitu_interp)*100.0)

        # slope, intercept, r_value, p_value, std_err = stats.linregress(dDinsitu_interp_s, dDfts)

        # odr, odrErr  = mf.orthoregress(dDinsitu_interp_s, dDfts, xerr=dDinsitu_interp_s*0.05, yerr=dDfts_e, InError=True)
        # slope       = float(odr[0])
        # intercept   = float(odr[1])

        # slope_e     = float(odrErr[0])
        # intercept_e = float(odrErr[1])
       
        # bias = mf.bias(dDinsitu_interp_s, dDfts)
        # rmse = mf.rmse(dDinsitu_interp_s, dDfts)
        

        # ax.text(0.025,0.63,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='green')
        # ax.text(0.025,0.56,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='green')
        # ax.text(0.025,0.49,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='green')
        # ax.text(0.025,0.42,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='green')

        #--
        # slope, intercept, r_value, p_value, std_err = stats.linregress(dDinsitu_interp_s, dDfts_s)

        # odr, odrErr  = mf.orthoregress(dDinsitu_interp_s, dDfts_s, xerr=dDinsitu_interp_s*0.05, yerr=dDfts_e, InError=True)
        # slope       = float(odr[0])
        # intercept   = float(odr[1])

        # slope_e     = float(odrErr[0])
        # intercept_e = float(odrErr[1])
       
        # bias = mf.bias(dDinsitu_interp_s, dDfts)
        # rmse = mf.rmse(dDinsitu_interp_s, dDfts)
        

        # ax.text(0.025,0.34,"Slope: {0:.2f} +/- {1:.2f}".format(slope, slope_e),transform=ax.transAxes,  fontsize=14, color='red')
        # ax.text(0.025,0.27,"Intercept: {0:.2f} +/- {1:.2f}".format(intercept, intercept_e),transform=ax.transAxes,  fontsize=14, color='red')
        # ax.text(0.025,0.21,"r-value: {0:.2f}".format(r_value),transform=ax.transAxes,  fontsize=14, color='red')
        # ax.text(0.025,0.14,"bias: {0:.2f}".format(bias),transform=ax.transAxes,  fontsize=14, color='red')

        #if saveFlg:     
        #    pdfsav.savefig(fig,dpi=200)
        #else:           
        #    plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
       

   #----------------------------
   #SAVE DATA in ASCII FILE
   #----------------------------

    # if not insFlg:

    #     Data = [Prf_h2o, ranErr_h2o, sysErr_h2o, totErr_h2o, Prf_dD, dDfts_RanErr, dDfts_SysErr, dDfts_TotErr, PresPrf_h2o]
    #     DT   = Dates_h2o

    #     Data  = np.asarray(Data)

    #     #-----------------------------------------
    #     # Open output file for year and write data
    #     #-----------------------------------------
    #     with open('/data/iortega/results/'+loc.lower()+'/WP_dD_'+loc.upper()+'2.ascii', 'w') as fopen:
    #         writer = csv.writer(fopen, delimiter='\t', lineterminator='\n')

    #         fopen.write('#Hannigan, J.W., Ortega, I\n')
    #         fopen.write('#National Center for Atmospheric Research\n')
    #         fopen.write('#Ground Based HR-FTIR Spectrometer at MLO\n')
    #         fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
    #         fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
    #         fopen.write('#DATA_INFO: Tropospheric Water Vapor and dD retrieved from HR-FTIR\n')
    #         fopen.write('#Number of Profiles: {}\n'.format(len(DT)))
    #         alt2print = ["%.2f"%v for v in alt_h2o]
    #         fopen.write('#Altitude [km]: ')
    #         writer.writerows([alt2print])
    #         fopen.write('#-------------------------------------------------\n')
    #         fopen.write('#   LABEL                  Units            Column\n')
    #         fopen.write('#-------------------------------------------------\n')
    #         fopen.write('#   H2O Profile           [ppm]                1  \n')
    #         fopen.write('#   H2O Random Error      [ppm]                2  \n')
    #         fopen.write('#   H2O Systematic Error  [ppm]                3  \n')
    #         fopen.write('#   H2O Total Error       [ppm]                4  \n')
    #         fopen.write('#   dD profile            [per mil]            5  \n')
    #         fopen.write('#   dD Random Error       [per mil]            6  \n')
    #         fopen.write('#   dD Systematic Error   [per mil]            7  \n')
    #         fopen.write('#   dD Total Error        [per mil]            8  \n')
    #         fopen.write('#   Pressure Profile      [mbar]               9  \n')
    #         fopen.write('#-------------------------------------------------\n')

    #         for i, d in enumerate(DT):

    #             fopen.write('# Date: {}\n'.format(d))

    #             Data4dt =  Data[:,i,:]
    #             Data4dt2 = [ ["%.3f"%v for v in row] for row in Data4dt]

    #             writer.writerows(zip(*(k for k in Data4dt2)))


    #    #-----------------------------------------
    #    # Open output file for year and write data
    #    #-----------------------------------------
    #     with open('/data/iortega/results/'+loc.lower()+'/WP_dD_Weighted_'+loc.upper()+'2.ascii', 'w') as fopen:
    #         writer = csv.writer(fopen, delimiter='\t', lineterminator='\n')

    #         fopen.write('#Hannigan, J.W., Ortega, I\n')
    #         fopen.write('#National Center for Atmospheric Research\n')
    #         fopen.write('#Ground Based HR-FTIR Spectrometer at MLO\n')
    #         fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
    #         fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
    #         fopen.write('#DATA_INFO: Tropospheric Water Vapor and dD retrieved from HR-FTIR\n')
    #         fopen.write('#Number of Profiles: {}\n'.format(len(DT)))
    #         alt2print = ["%.2f"%float(v) for v in [4.8, 7.9, 6.6]]
    #         fopen.write('#Altitude Weighted [km]: ')
    #         writer.writerows([alt2print])
    #         fopen.write('#---------------------------------------------------------------\n')
    #         fopen.write('#   LABEL                                Units            Column\n')
    #         fopen.write('#---------------------------------------------------------------\n')
    #         fopen.write('#   H2O Weighted                         [ppm]               1  \n')
    #         fopen.write('#   H2O Error Weighted                   [ppm]               2  \n')
    #         fopen.write('#   Note: fourth row is Total Column     [molec/cm2]         2  \n')
    #         fopen.write('#   dD Weighted                          [per mil]           3  \n')
    #         fopen.write('#   dD Error Weighted                    [per mil]           4  \n')
    #         fopen.write('#   Pressure Weighted                    [mbar]              5  \n')
    #         fopen.write('#---------------------------------------------------------------\n')

    #         for i, d in enumerate(DT):

    #             fopen.write('# Date: {}\n'.format(d))

    #             for p, pcol in enumerate(pCols):
    #                 row = [vmrP_h2o[str(p)][0,i], vmrP_h2o_e[str(p)][0,i], dDP_s[str(p)][0,i], dDP_e[str(p)][0,i], PresP[str(p)][0,i]]
    #                 row2 = [ "%.3f"%k for k in row] 
    #                 writer.writerows([row2])

    #             row = [TC_h2o[i], TC_h2o_e[i], TC_dD[i], TC_dD_e[i], -9999]
    #             row2 = [ "%.3e"%k for k in row] 
    #             writer.writerows([row2])

    if saveFlg:     
        pdfsav.close()
    else:           
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()    


if __name__ == "__main__":
    main()



