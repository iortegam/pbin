#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltFTSWRF.py
#
# Purpose:
#       1) Read and plot netCDF files provided by Gabi Pfister during FRAPPE (WRF-CHEM)
#       2) Integrate Model outputs and FTS
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
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
import numpy.ma as ma
import sys
import glob

from scipy import interpolate

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

from netCDF4 import Dataset
import dataModelOutClass as dm
import dataOutClass as dc
from collections                     import OrderedDict
from scipy import linspace, polyval, polyfit, sqrt, stats, randn


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
    #                     Global Inputs independently of FTS and/or Models
    #-----------------------------------------------------------------------------------------
    loc               = 'fl0'
    pCols             = [1.6, 8.0]       #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR
    
    pltAnalysis       = True
    saveFlg           = False 
    
    pltFile           =  '/data/iortega/results/'+loc.lower()+'/fig/FTS-WRF_'+loc.upper()+'_v2.pdf'

    #-----------------------------------------------------------------------------------------
    #                 Initializations for WRF-Chem during FRAPPE (Jul - Aug 2014)
    #-----------------------------------------------------------------------------------------
    #dataDirWRF        = '/data1/ancillary_data/'+loc.lower()+'/FRAPPE/'
    dataDirWRF        = '/ya4/Campaign/'+loc.upper()+'/FRAPPE/'
    #fileWRF1          = 'BoulderFL.nc'
    fileWRF1          = 'Boulder_FL_WRFV3.6.1_BTX_ECMWF_mozmozaic.nc'  
    fileWRF2          = 'BoulderFL_loc.nc'

    gasNameWRF        = ['c2h6', 'co', 'hcho', 'o3', 'no2', 'ch3oh', 'pan', 'c3h8', 'c2h2' ]              

    ReadWRF           = True
    pltWRF            = False
    calcWind          = False
    saveWRFFlg        = False
    
    pltWRFFile        = '/data/iortega/results/'+loc.lower()+'/fig/WRF-Chem_FRAPPE_v1.pdf'

    #-----------------------------------------------------------------------------------------
    #                 Initializations for CMAQ (EMISSIONS)
    #-----------------------------------------------------------------------------------------
    dataDirCMAQ       = '/ya4/Campaign/'+loc.upper()+'/FRAPPE/'
    fileCMAQ          = 'CMAQ_CB6_20140709_D02.nc'  

    ReadCMAQ          = False
    pltCMAQ           = False
    saveCMAQFlg       = False
    
    pltCMAQFile       = '/data/iortega/results/'+loc.lower()+'/fig/CMAQ_FRAPPE_v1.pdf'


    #-----------------------------------------------------------------------------------------
    #                 Initializations for FTIR
    #-----------------------------------------------------------------------------------------
    #gasName           = ['co', 'c2h2',  'c2h6', 'h2co', 'hcooh', 'o3']#, 'ch4']                   # Name of gas
    gasName           = ['co', 'c2h2', 'c2h6', 'h2co', 'o3', 'ch4', 'nh3']

    ver               = ['Current_v2', 'Current_v2', 'Current_v2', 'Current_WP_v6', 'Current_WP', 'Current_WP', 'Current_v2']          # Name of retrieval version to process
    ctlF              = ['sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v6.ctl', 'sfit4.ctl', 'sfit4_3.ctl', 'sfit4_v2.ctl']            # Name of ctl file
    maxRMS             = [2.0, 1.0, 1.0, 1.0, 4.0, 1.0, 3.0, 1.0, 2.0]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0, 0.5, 0.5, 0.5, 2.0, 0.5, 0.0, 0.0, 0.0]                      # Min DOFs for filtering

    #------
    # Flags - 1
    #------
    ReadFTS            = True   

    #------
    # Flags - 2
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
    mnths              = [6,7,8]                # Months to filter on (these are the months to include data)
    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2014
    imnth              = 7
    iday               = 15
    fyear              = 2014
    fmnth              = 8
    fday               = 20

    #-------------------------------
    # Decide how to create layers
    # above and below the tropopause
    #-------------------------------
    redFct  = 5.0
    hgtTrgU = 10.0       # Upper Layer thickness [km]
    hgtTrgL = False       # Lower Layer thickness [km];Set to False if you want surface value


                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    if loc == 'fl0':
        sLat              = 40.035#    40.4             #--LATITUDE OF BOULDER
        sLon              = -105.245          #--LONGITUDE OF BOULDER
    else:
        print 'Need a location (currently only fl0)'
        sys.exit()

        #lonB, latB = -105.245, 40.035      #LOCATION OF NCAR FOOTHILLS LAB in Boulder

    #-----------------------------------------------------------------------------------------
    #WRF-CHEM
    #-----------------------------------------------------------------------------------------
    if ReadWRF:
        DataWRF = dm.WRFClass(dataDirWRF, fileWRF1, fileWRF2, outFname= pltWRFFile, saveFlg= saveWRFFlg)
        DataWRF.ReadOutputWRF(gasNameWRF, pCols, sLat, sLon, calcWind = calcWind)
        if pltWRF:
            DataWRF.PltWRF()

    #-----------------------------------------------------------------------------------------
    #CMAQ
    #-----------------------------------------------------------------------------------------
    if ReadCMAQ:
        DataCMAQ = dm.CMAQClass(dataDirCMAQ, fileCMAQ, outFname= pltCMAQFile, saveFlg= saveCMAQFlg)
        DataCMAQ.ReadOutputCMAQ(sLat, sLon)
        #if pltCMAQ:
        #    DataCMAQ.PltCMAQ()

    #-------------------------------------------------
    #FTS
    #-------------------------------------------------
    if ReadFTS:
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
        for gasVer in statDataCl:
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=1)
            statDataCl[gasVer].readprfs([statDataCl[gasVer].PrimaryGas],retapFlg=0)
        
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
        avkSCFav       = OrderedDict()
        avkVMRav      = OrderedDict()

        aPrfVMR      = OrderedDict()
        aPrfMol      = OrderedDict()

        for j, gasVer in enumerate(statDataCl):
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
            aPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].aprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].aprfs['AIRMASS']))

            #----------------------------------------
            # This is the mixing ratio for DRY AIR!!!
            #----------------------------------------
            rPrfDry[gasVer] = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) / (1.0 - waterVMR[gasVer]) * sclfct       
            
            #----------------------------------
            # Read Summary data (For filtering)
            #----------------------------------
            statDataCl[gasVer].readsummary()
            rms[gasVer]     = np.asarray(statDataCl[gasVer].summary[statDataCl[gasVer].PrimaryGas+'_FITRMS'])       
        
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
                statDataCl[gasVer].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
                tot_rnd = np.array(statDataCl[gasVer].error['Total random uncertainty'])
                tot_sys = np.array(statDataCl[gasVer].error['Total systematic uncertainty'])
                totClmn_e[gasVer] = np.sqrt(tot_rnd**2 + tot_sys**2)
                    
                #---------------------
                # Get averaging kernel
                #---------------------   
                avkSCF[gasVer]  = np.delete(np.asarray(statDataCl[gasVer].error['AVK_scale_factor']),statDataCl[gasVer].inds,axis=0)
                avkVMR[gasVer]  = np.delete(np.asarray(statDataCl[gasVer].error['AVK_vmr']),statDataCl[gasVer].inds,axis=0)   
                dofs            = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
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
                    
                avkSCF[gasVer]     = np.delete(avkSCF[gasVer],statDataCl[gasVer].inds,axis=0)
                avkVMR[gasVer]     = np.delete(avkVMR[gasVer],statDataCl[gasVer].inds,axis=0)
                dofs               = np.diagonal(avkSCF[gasVer],axis1=1,axis2=2)
                avkSCFav[gasVer]     = np.mean(avkSCF[gasVer],axis=0)
                avkVMRav[gasVer]     = np.mean(avkVMR[gasVer],axis=0)
                
                dofsAvg[gasVer]    = np.diag(avkSCFav[gasVer])
                dofsAvg_cs[gasVer] = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]   
                
            #--------------------------------------
            # Remove retrieval data based on filter
            #--------------------------------------
            nfltr             = len(statDataCl[gasVer].inds)
            rms[gasVer]       = np.delete(rms[gasVer],statDataCl[gasVer].inds)
            ntot              = len(rms[gasVer])
            dates[gasVer]     = np.delete(dates[gasVer],statDataCl[gasVer].inds)
            totClmn[gasVer]   = np.delete(totClmn[gasVer],statDataCl[gasVer].inds)
            totClmn_e[gasVer] = np.delete(totClmn_e[gasVer],statDataCl[gasVer].inds)
            rPrfVMR[gasVer]   = np.delete(rPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
            rPrfMol[gasVer]   = np.delete(rPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)   
            rPrfDry[gasVer]   = np.delete(rPrfDry[gasVer],statDataCl[gasVer].inds,axis=0)   
            Airmass[gasVer]   = np.delete(Airmass[gasVer],statDataCl[gasVer].inds,axis=0)
            TCdry[gasVer]     = np.delete(TCdry[gasVer],statDataCl[gasVer].inds)

            aPrfVMR[gasVer]   = np.delete(aPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
            aPrfMol[gasVer]   = np.delete(aPrfMol[gasVer],statDataCl[gasVer].inds,axis=0) 
        
            #-------------------------------------------------
            # Calculate layer information for each observation
            #-------------------------------------------------
            LayVMR[gasVer]    = np.zeros(len(dates[gasVer]))
            DryLayVMR[gasVer] = np.zeros(len(dates[gasVer]))
            upperHgt[gasVer]  = np.zeros(len(dates[gasVer]))
            lowerHgt[gasVer]  = np.zeros(len(dates[gasVer]))
            LayThk[gasVer]    = np.zeros(len(dates[gasVer]))
            LayDOF[gasVer]    = np.zeros(len(dates[gasVer]))
            
            for i,dateDay in enumerate(dates[gasVer]):        
                indU             = dc.nearestind(hgtTrgU,alt[gasVer])
                if hgtTrgL: indL = dc.nearestind(hgtTrgL,alt[gasVer])
                else:       indL = -1
                
                upperHgt[gasVer][i]  = alt[gasVer][indL]
                lowerHgt[gasVer][i]  = alt[gasVer][indU]
                
                if indL < 0:
                    LayDOF[gasVer][i]    = np.sum(dofs[i,indU+1:])
                    LayVMR[gasVer][i]    = np.average(rPrfVMR[gasVer][i,indU+1:],weights=Airmass[gasVer][i,indU+1:])
                    DryLayVMR[gasVer][i] = np.average(rPrfDry[gasVer][i,indU+1:],weights=Airmass[gasVer][i,indU+1:])
                else:
                    LayDOF[gasVer][i]    = np.sum(dofs[i,indU+1:indL+1])
                    LayVMR[gasVer][i]    = np.average(rPrfVMR[gasVer][i,indU+1:indL+1],weights=Airmass[gasVer][i,indU+1:indL+1])
                    DryLayVMR[gasVer][i] = np.average(rPrfDry[gasVer][i,indU+1:indL+1],weights=Airmass[gasVer][i,indU+1:indL+1])
             
                LayThk[gasVer][i]    = alt[gasVer][indU] - alt[gasVer][indL]
        
        #----------
        # Plot Data
        #----------
    vmrPgas    = {}
    vmrPWRFgas = {}
    Dategas    = {}

    GasStrall = []


    if pltAnalysis:

        if saveFlg: pdfsav = PdfPages(pltFile)
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()
        mondays      = WeekdayLocator(MONDAY)
        DateFmt      = DateFormatter('%b %d')
       

        for k, gasVer in enumerate(statDataCl):

            gas2plt = gasName[k].lower()
            gasSTR  = dm.GasName(gas2plt)
            GasStrall.append(gasSTR)

            if   gas2plt == 'co':    scFactor = 1.0
            elif gas2plt == 'c2h6':  scFactor = 1.0
            else: scFactor = 1.0

            if gas2plt == 'h2co': gas2plt = 'hcho'

            maxalt = 16.0

            #---------------------------------
            # Defining variables (FTIR)
            #---------------------------------
            indsalt      = np.where(alt[gasVer] <= maxalt)[0]
            alt_i        = alt[gasVer][indsalt]
            TC_i         = totClmn[gasVer]
            TCe_i        = totClmn_e[gasVer]
            Dates_i      = dates[gasVer]
            Prf_i        = rPrfVMR[gasVer][:, indsalt]
            avkSCFavs_i  = avkSCFav[gasVer][indsalt[0]:, indsalt[0]:]
            aPrf_i        = aPrfVMR[gasVer][:, indsalt]
            Airmass_i     = Airmass[gasVer][:, indsalt]
            rPrfMol_i     = rPrfMol[gasVer][:,indsalt]
            TCpc_i        = np.sum(rPrfMol_i,axis=1)

            indh1 = mf.nearestind(pCols[0], alt_i)
            indh2 = mf.nearestind(pCols[1], alt_i)               
            vmrP  = np.average(Prf_i[:,indh2:indh1],axis=1,weights=Airmass_i[:,indh2:indh1]) 

            vmrPgas.setdefault(gas2plt,[]).append(vmrP)
            Dategas.setdefault(gas2plt,[]).append(Dates_i)

            #---------------------------------
            # Defining variables (WRF-CHEM) - Data is hourly
            #---------------------------------

            altWRF      = DataWRF.WRF['midpoints'][DataWRF.WRF['indsLoc'][0],DataWRF.WRF['indsLoc'][1], :, 0]
            indsalt     = np.where(altWRF <= maxalt)[0]
            altWRF      = altWRF[indsalt]

            DatesWRF    = np.asarray(DataWRF.WRF['dates'])
            #PrfWRFi     = np.asarray(DataWRF.WRF['GasPrf_'+gas2plt][DataWRF.WRF['indsLoc'][0],DataWRF.WRF['indsLoc'][1], :,:])
            AirmassWRFi = np.asarray(DataWRF.WRF['airmass'][DataWRF.WRF['indsLoc'][0],DataWRF.WRF['indsLoc'][1], :,:])        
            #print PrfWRFi.shape
            
            try:
                TCWRFi      = np.asarray(DataWRF.WRF['GasTC_'+gas2plt][DataWRF.WRF['indsLoc'][0],DataWRF.WRF['indsLoc'][1], :])
                PrfWRFj     = np.asarray(DataWRF.WRF['GasPrf_'+gas2plt][:,:, :,:])

                PrfWRFi     = np.nanmean(PrfWRFj, axis=tuple(range(0, 2)))
                PrfWRFStdi  = np.nanstd(PrfWRFj, axis=tuple(range(0, 2)))
                PrfWRF      = PrfWRFi[indsalt, :]
                PrfStdWRF   = PrfWRFStdi[indsalt, :]
                AirmassWRF  = AirmassWRFi[indsalt, :]

            except KeyError, e:
                npnts       = np.shape(DatesWRF)[0]
                nlvls       = np.shape(altWRF)[0]
            
                PrfWRF      = np.zeros((nlvls, npnts))
                PrfStdWRF   = np.zeros((nlvls, npnts))
                AirmassWRF  = np.zeros((nlvls, npnts))
                PrfWRF[:,:] = float('nan')
                PrfWRF[:,:] = float('nan')
                AirmassWRF[:,:] = float('nan')

            #---------------------------------
            # Smoothing WRF-CHEM using FTIR AK and apriori
            #---------------------------------
            PrfWRF_interp     = interpolate.interp1d(altWRF, PrfWRF, axis=0, fill_value='extrapolate', bounds_error=False)(alt_i)
            PrfWRFStd_interp  = interpolate.interp1d(altWRF, PrfStdWRF, axis=0, fill_value='extrapolate', bounds_error=False)(alt_i)
            
            AirmassWRF_interp = interpolate.interp1d(altWRF, AirmassWRF, axis=0, fill_value='extrapolate', bounds_error=False)(alt_i)

            aPrfm             = np.mean(aPrf_i, axis=0)

            PrfWRFs           = np.zeros( (len(DatesWRF), len(alt_i)) )
            PrfStdWRFs        = np.zeros( (len(DatesWRF), len(alt_i)) )
           
            for itime in range(len(DatesWRF)):
                #The equations below are equivalent
                PrfWRFs[itime, :]  = aPrfm + np.dot(avkSCFavs_i, (PrfWRF_interp[:, itime] -  aPrfm))         # x_s = Apriori(FTS) + AK(Prf(WRF-CHEM) - Apriori(FTS))
                #PrfWRFs[itime, :]  = np.dot(avkSCFavs, PrfWRF_interp[itime, :]) + np.dot(  np.identity(avkSCFavs.shape[0]) - avkSCFavs,    aPrfm )   # x_s = AK*Prf(WRF-CHEM) + (I-AK)*Apriori
                PrfStdWRFs[itime, :]  = aPrfm + np.dot(avkSCFavs_i, (PrfWRFStd_interp[:, itime] -  aPrfm))

            #BELOW I USE THE MEAN AIRMASS FROM THE FTIR, THE AIRMASS FROM THE MODEL IS NOT ADEQUATE
            AirmassFTSMean = np.mean(Airmass_i, axis=0)
            AirmassFTSstd = np.std(Airmass_i, axis=0)

            TCWRF        = np.sum( (PrfWRFs/1e9) *AirmassFTSMean, axis=1)
            TCStdWRF     = np.sum( (PrfStdWRFs/1e9) *AirmassFTSMean, axis=1)
            vmrPWRF      = np.average(PrfWRFs[:,indh2:indh1], axis=1, weights=AirmassFTSMean[indh2:indh1]) 
            vmrPStdWRF   = np.average(PrfStdWRFs[:,indh2:indh1], axis=1, weights=AirmassFTSMean[indh2:indh1]) 

               
 
            #---------------------------------
            # Plot : Averaging Kernel Smoothing Function (row of avk)
            #---------------------------------
            fig       = plt.figure(figsize=(9,9))
            gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
            ax        = plt.subplot(gs[0])
            axb       = plt.subplot(gs[1])
            cm        = plt.get_cmap(clmap)
            cNorm     = colors.Normalize(vmin=np.min(alt_i), vmax=np.max(alt_i))
            scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
            scalarMap.set_array(alt_i)
            ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt_i])
            
            for i in range(len(alt_i)):
                ax.plot(avkSCFavs_i[i,:],alt_i)
                
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('Averaging Kernels', fontsize=14)
            ax.grid(True)
            cbar = fig.colorbar(scalarMap, orientation='vertical')
            cbar.set_label('Altitude [km]', fontsize=14)
            ax.set_title(gasSTR + ' Averaging Kernels Scale Factor', fontsize=14)
            ax.tick_params(labelsize=14)
            
            axb.plot(np.sum(avkSCFavs_i,axis=0), alt_i,color='k')
            axb.grid(True)
            axb.set_xlabel('Averaging Kernel Area', fontsize=14)
            axb.tick_params(axis='x',which='both',labelsize=8)
            #axb.tick_params(labelsize=14) 

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #---------------------------------
            # Plot : Prf Mean
            #---------------------------------

            PrfmeanWRF = np.mean(PrfWRFs, axis=0)
            prfSTDWRF  = np.std(PrfWRFs, axis=0)

            prfMean = np.mean(Prf_i, axis=0)
            prfSTD  = np.std(Prf_i, axis=0)

            PrfmeanWRFi = np.mean(PrfWRF, axis=1)

            fig,ax  = plt.subplots(figsize=(7,9))

            ax.plot(prfMean,alt_i, linewidth = 2.0, color='k', label='FTIR')
            ax.scatter(prfMean,alt_i, facecolors='white', s=60, color='k')
            ax.fill_betweenx(alt_i,prfMean-prfSTD,prfMean+prfSTD, alpha=0.5, color='k') 

            ax.plot(PrfmeanWRFi*scFactor,altWRF, linewidth = 2.0, color='r', label='WRF-Chem x '+' '+str(scFactor))
            ax.scatter(PrfmeanWRFi*scFactor,altWRF, facecolors='white', s=30, color='r') 

            ax.plot(PrfmeanWRF*scFactor,alt_i, linewidth = 2.0, color='green', label='WRF-Chem smoothed x '+' '+str(scFactor))
            ax.scatter(PrfmeanWRF*scFactor,alt_i, facecolors='white', s=60, color='green')
            ax.fill_betweenx(alt_i,PrfmeanWRF*scFactor-prfSTDWRF,PrfmeanWRF*scFactor+prfSTDWRF,alpha=0.5,color='green') 

    
            ax.set_title('Mean Profile of '+gasSTR, fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            ax.grid(True,which='both')
            ax.tick_params(labelsize=14)
            ax.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #-------------------------------------------------
            # Time series of weighted VMR
            #-------------------------------------------------

            fig, (ax, ax2) = plt.subplots(2, figsize=(10,7), sharex=True)

            totErr_frac = TCe_i / TC_i

            ax.plot(Dates_i, vmrP, color='k', label='FTIR')
            ax.errorbar(Dates_i,vmrP, yerr=vmrP*totErr_frac*2.0, fmt='o', markersize=0, color='k', ecolor='k')
            ax.scatter(Dates_i, vmrP, facecolors='white', s=60, color='k')
        
            ax.plot(DatesWRF, vmrPWRF*scFactor, color='r', label='WRF-Chem x '+' '+str(scFactor))
            ax.errorbar(DatesWRF,vmrPWRF, yerr=vmrPStdWRF*2.0, fmt='o', markersize=0, color='r', ecolor='r')
            ax.scatter(DatesWRF, vmrPWRF*scFactor, facecolors='white', s=60, color='r')

            ax.grid(True)
            ax.set_ylabel('VMR [ppb$_v$]',fontsize=14)
            ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                          multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)

            ax.legend(prop={'size':12})

            ax2.plot(Dates_i, TCpc_i, color='k', label='FTIR')
            ax2.errorbar(Dates_i,TCpc_i, yerr=TCpc_i*totErr_frac*2.0, fmt='o', markersize=0, color='k', ecolor='k')
            ax2.scatter(Dates_i, TCpc_i, facecolors='white', s=60, color='k')
        
            ax2.plot(DatesWRF, TCWRF, color='r', label='WRF-Chem x '+' '+str(scFactor))
            ax2.errorbar(DatesWRF,TCWRF, yerr=TCStdWRF*2.0, fmt='o', markersize=0, color='r', ecolor='r')
            ax2.scatter(DatesWRF, TCWRF, facecolors='white', s=60, color='r')

            ax2.grid(True)
            ax2.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax2.set_xlabel('Date',fontsize=14)
            ax2.set_title(gasSTR + ' Partial Column between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                          multialignment='center',fontsize=14)
            ax2.tick_params(labelsize=14)
            ax2.xaxis.set_major_locator(mondays)
            ax2.xaxis.set_minor_locator(dayLc)
            ax2.xaxis.set_major_formatter(DateFmt)
            ax2.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)

            #-------------------------------------------------
            # Time series of weighted VMR - WRF-CHem interpolated to FTIR
            #-------------------------------------------------

            doy    = mf.toYearFraction(Dates_i)
            doyWRF = mf.toYearFraction(DatesWRF) 

            vmrPWRF_interp       = interpolate.interp1d(doyWRF, vmrPWRF, fill_value='extrapolate', bounds_error=False)(doy)
            TCWRF_interp         = interpolate.interp1d(doyWRF, TCWRF, fill_value='extrapolate', bounds_error=False)(doy)

            vmrPWRFgas.setdefault(gas2plt,[]).append(vmrPWRF_interp)

            vmrPStdWRF_interp    = interpolate.interp1d(doyWRF, vmrPStdWRF, fill_value='extrapolate', bounds_error=False)(doy)
            TCStdWRF_interp      = interpolate.interp1d(doyWRF, TCStdWRF, fill_value='extrapolate', bounds_error=False)(doy)

            #fig, (ax, ax2) = plt.subplots(2, figsize=(10,7), sharex=True)
            fig, ax = plt.subplots(figsize=(12,6), sharex=True)

            #ax.plot(Dates_i, vmrP, color='k')
            ax.errorbar(Dates_i,vmrP, yerr=vmrP*totErr_frac, fmt='o', markersize=0, color='k', ecolor='k')
            ax.scatter(Dates_i, vmrP, facecolors='white', s=100, color='k', label='FTIR')
        
            #ax.plot(Dates_i, vmrPWRF_interp*scFactor, color='r')
            ax.errorbar(Dates_i,vmrPWRF_interp, yerr=vmrPStdWRF_interp, fmt='o', markersize=0, color='r', ecolor='r')
            ax.scatter(Dates_i, vmrPWRF_interp*scFactor, facecolors='white', s=100, color='r', label='WRF-Chem')

            ax.grid(True)
            ax.set_ylabel('VMR [ppb$_v$]',fontsize=14)
            ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                          multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.set_xlabel('Date',fontsize=14)

            ax.legend(prop={'size':12})

            #ax2.plot(Dates_i, TCpc_i, color='k')
            # ax2.errorbar(Dates_i,TCpc_i, yerr=TCpc_i*totErr_frac*2.0, fmt='o', markersize=0, color='k', ecolor='k')
            # ax2.scatter(Dates_i, TCpc_i, facecolors='white', s=100, color='k', label='FTIR')
        
            # #ax2.plot(Dates_i, TCWRF_interp, color='r', label='WRF-Chem x '+' '+str(scFactor))
            # ax2.errorbar(Dates_i,TCWRF_interp, yerr=TCStdWRF_interp*2.0, fmt='o', markersize=0, color='r', ecolor='r')
            # ax2.scatter(Dates_i, TCWRF_interp, facecolors='white', s=100, color='r', label='WRF-Chem')

            # ax2.grid(True)
            # ax2.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            # ax2.set_xlabel('Date',fontsize=14)
            # ax2.set_title(gasSTR + ' Partial Column between '+str(pCols[0])+' to '+str(pCols[1])+' km',
            #               multialignment='center',fontsize=14)
            # ax2.tick_params(labelsize=14)
            # ax2.xaxis.set_major_locator(mondays)
            # ax2.xaxis.set_minor_locator(dayLc)
            # ax2.xaxis.set_major_formatter(DateFmt)
            # ax2.legend(prop={'size':12})

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #-------------------------------------------------
            # Correlation Plots
            #-------------------------------------------------

            # slope, intercept, r_value, p_value, std_err = stats.linregress(vmrP,vmrPWRF_interp)

            # print '\n', gasSTR
            # print 'Slope = ', str(slope)
            # print 'Intercept = ', str(intercept)
            # print 'r-value = ', str(r_value)

            # z = np.polyfit(vmrP,vmrPWRF_interp, 1,  full=False, cov=True)

            # fig, (ax, ax2) = plt.subplots(2, figsize=(10, 12), sharex=False)
            # Fit = slope*np.asarray(vmrP) + (intercept)

            # ax.scatter(vmrP, vmrPWRF_interp, s =25, c='red', alpha=0.75, edgecolors='k')
            # ax.errorbar(vmrP,vmrPWRF_interp, yerr=vmrPStdWRF_interp*2.0,  xerr=vmrP*totErr_frac*2.0, fmt='o', markersize=0, color='k', ecolor='k')
            # ax.plot(vmrP, Fit, color='red',  linewidth=2.0, linestyle='--')
            # ax.grid(True)
            # ax.set_ylabel('VMR [ppb$_v$] - WRF-Chem', fontsize=14)
            # ax.set_xlabel('VMR [ppb$_v$] - FTIR', fontsize=14)
            # ax.set_xlim(( np.min(vmrP)- np.min(vmrP)*0.2, np.max(vmrP) +np.max(vmrP)*0.2  ))
            # ax.set_ylim(( np.min(vmrP)- np.min(vmrP)*0.2, np.max(vmrP) +np.max(vmrP)*0.2  ))
            # #ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
            # #              multialignment='center',fontsize=14)
            # ax.tick_params(labelsize=14)

            # slope, intercept, r_value, p_value, std_err = stats.linregress(TCpc_i,TCWRF_interp)

            # print 'Slope = ', str(slope)
            # print 'Intercept = ', str(intercept)
            # print 'r-value = ', str(r_value)

            
            # Fit = slope*np.asarray(TCpc_i) + (intercept)

            # ax2.scatter(TCpc_i, TCWRF_interp, s =25, c='red', alpha=0.75, edgecolors='k')
            # ax2.errorbar(TCpc_i,TCWRF_interp, yerr=TCStdWRF_interp,  xerr=TCWRF_interp*totErr_frac*2.0, fmt='o', markersize=0, color='k', ecolor='k')
            # ax2.plot(TCpc_i, Fit, color='red',  linewidth=2.0, linestyle='--')
            # ax2.grid(True)
            # ax2.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$] - WRF-Chem', fontsize=14)
            # ax2.set_xlabel('Partial Column [molecules$\cdot$cm$^{-2}$] - FTIR', fontsize=14)
            # #ax2.set_title(gasSTR + ' Partial Column between '+str(pCols[0])+' to '+str(pCols[1])+' km',
            # #              multialignment='center',fontsize=14)
            # ax2.set_xlim(( np.min(TCpc_i)- np.min(TCpc_i)*0.2, np.max(TCpc_i) +np.max(TCpc_i)*0.2  ))
            # ax2.set_ylim(( np.min(TCpc_i)- np.min(TCpc_i)*0.2, np.max(TCpc_i) +np.max(TCpc_i)*0.2  ))
            # ax2.tick_params(labelsize=14)

            # if saveFlg:     pdfsav.savefig(fig,dpi=200)
            # else:           plt.show(block=False)

        
        #-------------------------------------------------
        #
        #-------------------------------------------------

        stfile = '/ya4/Campaign/'+loc.upper()+'/FRAPPE//flexpart2/FLX_SCF_12h.txt'

        ckFile(stfile)
        #stfile = file(stfile, 'r')
        cols, indexToName = mf.getColumns(stfile, headerrow=1, delim=',', header=True)

        SCsd = np.asarray(cols['Date'])  
        SCdate = [dt.datetime(int(d[0:4]), int(d[4:6]), int(d[6:8]), int(d[8:10])) for d in SCsd]
        SCdate = np.asarray(SCdate)
        SCNE   = np.asarray(cols['SC_NE'])
        SCNW   = np.asarray(cols['SC_NW'])
        SCSE   = np.asarray(cols['SC_SE'])
        SCSW   = np.asarray(cols['SC_SW'])


        #-------------------------------------------------
        for k, gasVer in enumerate(statDataCl):

            gas2plt = gasName[k].lower()

            if gas2plt != 'co':
                gasSTR  = dm.GasName(gas2plt)

                if gas2plt == 'h2co': gas2plt = 'hcho'

                #-----------------------------------------
                #DEFINE VARIABLES
                #-----------------------------------------
                yWRF   = np.asarray(vmrPWRFgas[gas2plt][0])
                yFTS   = np.asarray(vmrPgas[gas2plt][0])
                yDate  = np.asarray(Dategas[gas2plt][0])
                ydoy   = mf.toYearFraction(yDate)     

                yTC_i   = np.asarray(totClmn[gasVer])
                yTCe_i  = np.asarray(totClmn_e[gasVer])
                
                xWRF   = np.asarray(vmrPWRFgas['co'][0])
                xFTS   = np.asarray(vmrPgas['co'][0])
                xDate  = np.asarray(Dategas['co'][0])
                xdoy   = mf.toYearFraction(xDate)

                xTC_i   = np.asarray(totClmn['co_'+ver[0]])   
                xTCe_i  = np.asarray(totClmn_e['co_'+ver[0]])

                xFTS_interp     = interpolate.interp1d(xdoy, xFTS, axis=0, fill_value='extrapolate', bounds_error=False)(ydoy)
                xWRF_interp     = interpolate.interp1d(xdoy, xWRF, axis=0, fill_value='extrapolate', bounds_error=False)(ydoy)

                xFTS_interp_tc     = interpolate.interp1d(xdoy, xTC_i, axis=0, fill_value='extrapolate', bounds_error=False)(ydoy)
                xFTS_interp_tce    = interpolate.interp1d(xdoy, xTCe_i, axis=0, fill_value='extrapolate', bounds_error=False)(ydoy)
                
                RWRF = yWRF/xWRF_interp
                RFTS = yFTS/xFTS_interp

                RFTSe = RFTS*np.sqrt((xFTS_interp_tce/xFTS_interp_tc)**2 + (yTCe_i/yTC_i)**2)

                #-----------------------------------------
                #CALCULATING ANOMALIES
                #-----------------------------------------
                daysall = [dt.date(da.year, da.month, da.day) for da in yDate]

                dailyVals        = mf.dailyAvg(yFTS, yDate, dateAxis=1, meanAxis=0)
                DT_d             = dailyVals['dates']
                yFTS_d              = dailyVals['dailyAvg']

                dailyVals        = mf.dailyAvg(xFTS_interp, yDate, dateAxis=1, meanAxis=0)
                xFTS_d              = dailyVals['dailyAvg']

                dailyVals        = mf.dailyAvg(yWRF, yDate, dateAxis=1, meanAxis=0)
                yWRF_d              = dailyVals['dailyAvg']

                dailyVals        = mf.dailyAvg(xWRF_interp, yDate, dateAxis=1, meanAxis=0)
                xWRF_d              = dailyVals['dailyAvg']
                
                deltaVMRy      = []
                deltaVMRx      = []
                deltaVMRy_WRF  = []
                deltaVMRx_WRF  = []

                for i, item in enumerate(DT_d):
                    diff = np.asarray(daysall) - item
                    inds = np.where( diff == dt.timedelta(0) )[0]

                    deltaVMRy.append(yFTS[inds] - yFTS_d[i])
                    deltaVMRx.append(xFTS_interp[inds] - xFTS_d[i])

                    deltaVMRy_WRF.append(yWRF[inds] - yWRF_d[i])
                    deltaVMRx_WRF.append(xWRF_interp[inds] - xWRF_d[i])

                deltaVMRy = np.asarray(deltaVMRy)
                dVMRy = [y for x in deltaVMRy for y in x]
                dVMRy = np.asarray(dVMRy)

                deltaVMRx = np.asarray(deltaVMRx)
                dVMRx = [y for x in deltaVMRx for y in x]
                dVMRx = np.asarray(dVMRx)

                deltaVMRy_WRF = np.asarray(deltaVMRy_WRF)
                dVMRy_WRF = [y for x in deltaVMRy_WRF for y in x]
                dVMRy_WRF = np.asarray(dVMRy_WRF)

                deltaVMRx_WRF = np.asarray(deltaVMRx_WRF)
                dVMRx_WRF = [y for x in deltaVMRx_WRF for y in x]
                dVMRx_WRF = np.asarray(dVMRx_WRF)

                #-----------------------------------------
                #FIGURE OF THE RATIO 
                #-----------------------------------------
                ydoysc   = mf.toYearFraction(SCdate)
                SCNE_interp     = interpolate.interp1d(ydoysc, SCNE, axis=0, bounds_error=False)(ydoy)
                SCNW_interp     = interpolate.interp1d(ydoysc, SCNW, axis=0, bounds_error=False)(ydoy)
                SCSE_interp     = interpolate.interp1d(ydoysc, SCSE, axis=0, bounds_error=False)(ydoy)
                SCSW_interp     = interpolate.interp1d(ydoysc, SCSW, axis=0, bounds_error=False)(ydoy)
                SCwEST =  SCSW_interp + SCNW_interp

                slope, intercept, r_value, p_value, std_err = stats.linregress(RFTS, SCNE_interp)
                print '\nr-value (SC) = ', str(r_value)

                fig, ax= plt.subplots(figsize=(12,6), sharex=True)
                #fig , (ax, ax1) = plt.subplots(2, figsize=(10,8), gridspec_kw = {'height_ratios':[3, 1]}, sharex=True)

                ax2 = ax.twinx()
                #ax2.plot(yDate, SCNE_interp, color='gray', linestyle='--', label='FLEXPART', markersize=16, linewidth=3)
                #ax2.scatter(yDate, SCNE_interp, facecolors='white', s=80, color='gray')

                #ax2.fill_between(yDate, 0, SCNE_interp, color='gray', alpha=0.25)
                #ax2.fill_between(yDate, 0, SCNW_interp, color='yellow', alpha=0.25)

                #ax2.bar(yDate-dt.timedelta(days=0.3), SCNE_interp, width=0.1,  align='center', color = 'r', label='North-East', alpha=0.4 )
                #ax2.bar(yDate, SCSE_interp,width=0.1, align='center', color = 'blue',  label='South-East', alpha=0.4)
                #ax2.bar(yDate+dt.timedelta(days=0.3), SCwEST,width=0.1, align='center', color = 'green', label='West', alpha=0.4)

                ax2.bar(yDate, SCNE_interp, width=0.2,  align='center', color = 'r', label='North-East', alpha=0.3 )
                ax2.bar(yDate, SCSE_interp,width=0.2, align='center', color = 'blue',  label='South-East', alpha=0.3,  bottom =sumzip(SCNE_interp))
                ax2.bar(yDate, SCwEST,width=0.2, align='center', color = 'green', label='West', alpha=0.3, bottom =sumzip(SCNE_interp,SCSE_interp))

                ax2.legend(loc=4, ncol=3)
                #ax1.bar(ind+0.27,RateCAM, width = 0.27, align='center', color = 'b', yerr=RateCAMe, ecolor = 'k', label = 'CAM-Chem')
        
                #ax1.set_xticks(ind)
                ax2.set_ylabel('Sector Fraction', fontsize = 18)
                #ax1.set_xticks(ind+0.27/2.0)
            
                #ax1.legend(prop={'size':12}, loc=4)
                ax2.xaxis.set_tick_params(which='major',labelsize=18)
                ax2.yaxis.set_tick_params(which='major',labelsize=18)
                ax2.set_ylim(0,1.05)
                ax2.tick_params(labelsize=18)                

                #for tl in ax2.get_yticklabels():
                #    tl.set_color('gray')
                
                ax.errorbar(yDate,RFTS, yerr=RFTSe, fmt='o', markersize=0, color='red', ecolor='red', elinewidth=2)
                ax.scatter(yDate, RFTS, facecolors='red', s=100, color='k', label='HR-FTIR')
            
                #ax.plot(yDate, RWRF, color='r', label='WRF-Chem', markersize=16, linewidth=3)
                ax.scatter(yDate, RWRF, facecolors='blue', s=100, color='k', label='WRF-Chem')

                ax.grid(True)
                #ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                #              multialignment='center',fontsize=14)
                ax.tick_params(labelsize=18)
                ax.xaxis.set_major_locator(mondays)
                ax.xaxis.set_minor_locator(dayLc)
                ax.xaxis.set_major_formatter(DateFmt)
                ax.set_ylabel(gasSTR + ' to CO ratio - Weighted VMR', fontsize=18)
                ax.set_xlabel('Date', fontsize=18)

                #ax2.stem(yDate, SCNE_interp, linefmt='b-')

                if not ( gas2plt == 'nh3' or gas2plt == 'ch4'): ax.legend(prop={'size':16}, bbox_to_anchor=(0.0, 1.01, 1., 0.102), loc='center', ncol=2, borderaxespad=0.0)


                #ax2 = ax.twinx()
                #ax2.plot(yDate, SCNE_interp, color='gray', linestyle='--', label='FLEXPART', markersize=16, linewidth=3)
                #ax2.scatter(yDate, SCNE_interp, facecolors='white', s=80, color='gray')

                #ax2.fill_between(yDate, 0, SCNE_interp, color='gray', alpha=0.5)
                
                # ax2.plot(yDate, SCNW_interp, color='g', linestyle='--', label='FLEXPART', markersize=16, linewidth=3)
                # ax2.scatter(yDate, SCNW_interp, facecolors='white', s=80, color='g')

                # ax2.plot(yDate, SCSE_interp, color='gray', linestyle='--', label='FLEXPART', markersize=16, linewidth=3)
                # ax2.scatter(yDate, SCSE_interp, facecolors='white', s=80, color='gray')

                # ax1.bar(yDate, SCNE_interp,  align='center', color = 'r', ecolor = 'k', label='North-East', alpha=0.4)
                # ax1.bar(yDate, SCSE_interp, align='center', color = 'blue', ecolor = 'k', label='South-East', alpha=0.4,)
                # ax1.bar(yDate, SCwEST, align='center', color = 'green', ecolor = 'k', label='West', alpha=0.4,)
                # #ax1.bar(ind+0.27,RateCAM, width = 0.27, align='center', color = 'b', yerr=RateCAMe, ecolor = 'k', label = 'CAM-Chem')
                # ax1.yaxis.grid(True)
                # #ax1.set_xticks(ind)
                # ax1.set_ylabel('Fraction', fontsize = 18)
                # #ax1.set_xticks(ind+0.27/2.0)
                # ax1.axhline(0, color='black', lw=1)
                # #ax1.legend(prop={'size':12}, loc=4)
                # ax1.xaxis.set_tick_params(which='major',labelsize=18)
                # ax1.yaxis.set_tick_params(which='major',labelsize=18)
                # ax1.set_xlabel('Date', fontsize=18)

                if saveFlg:     pdfsav.savefig(fig,dpi=200)
                else:           plt.show(block=False)

                #---
                fig , (ax, ax1) = plt.subplots(2, figsize=(10,8), gridspec_kw = {'height_ratios':[3, 1]}, sharex=True)
                #print RFTS
                #print yDate

                #ax.plot(yDate, RFTS, color='k', label='FTIR',  markersize=16, linewidth=3)
                ax.errorbar(yDate,RFTS, yerr=RFTSe, fmt='o', markersize=0, color='red', ecolor='red',  elinewidth=2)
                ax.scatter(yDate, RFTS, facecolors='red', s=100, color='k', label='HR-FTIR')
            
                #ax.plot(yDate, RWRF, color='r', label='WRF-Chem', markersize=16, linewidth=3)
                ax.scatter(yDate, RWRF, facecolors='blue', s=100, color='k', label='WRF-Chem')

                ax.grid(True)
                #ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                #              multialignment='center',fontsize=14)
                ax.tick_params(labelsize=18)
                ax.xaxis.set_major_locator(mondays)
                ax.xaxis.set_minor_locator(dayLc)
                ax.xaxis.set_major_formatter(DateFmt)
                ax.set_ylabel(gasSTR + ' to CO ratio - Weighted VMR', fontsize=18)
                #ax.set_xlabel('Date', fontsize=18)

                #ax2.stem(yDate, SCNE_interp, linefmt='b-')

                #if not ( gas2plt == 'nh3' or gas2plt == 'ch4'): ax.legend(prop={'size':16}, bbox_to_anchor=(0.0, 1.01, 1., 0.102), loc='center', ncol=2, borderaxespad=0.0)
                if not ( gas2plt == 'nh3' or gas2plt == 'ch4'): ax.legend(prop={'size':16})


                ax1.bar(yDate, SCNE_interp, width=0.2,  align='center', color = 'r', label='North-East', alpha=1 )
                ax1.bar(yDate, SCSE_interp,width=0.2, align='center', color = 'blue',  label='South-East', alpha=1,  bottom =sumzip(SCNE_interp))
                ax1.bar(yDate, SCwEST,width=0.2, align='center', color = 'green', label='West', alpha=1, bottom =sumzip(SCNE_interp,SCSE_interp))

                ax1.legend(ncol=3, bbox_to_anchor=(0.0, 1.04, 1., 0.102), loc='center', borderaxespad=0.0)
                #ax1.bar(ind+0.27,RateCAM, width = 0.27, align='center', color = 'b', yerr=RateCAMe, ecolor = 'k', label = 'CAM-Chem')
        
                #ax1.set_xticks(ind)
                ax1.set_ylabel('Sector Fraction', fontsize = 18)
                #ax1.set_xticks(ind+0.27/2.0)
            
                #ax1.legend(prop={'size':12}, loc=4)
                ax1.xaxis.set_tick_params(which='major',labelsize=18)
                ax1.yaxis.set_tick_params(which='major',labelsize=18)
                ax1.set_ylim(0,1.05)
                ax1.tick_params(labelsize=18) 
                ax1.set_xlabel('Date', fontsize=18)
            

                if saveFlg:     pdfsav.savefig(fig,dpi=200)
                else:           plt.show(block=False)

                #-----------------------------------------
                #FIGURE OFANOMALY 
                #-----------------------------------------

                fig2, ax = plt.subplots(figsize=(9,7), sharex=True)

                #ax.plot(dVMRy, dVMRx, color='k', label='FTIR',  markersize=16, linewidth=3)
                #ax.errorbar(yDate,RFTS, yerr=RFTSe*2.0, fmt='o', markersize=0, color='k', ecolor='k')
                ax.scatter(dVMRy, dVMRx, facecolors='white', s=80, color='k')
            
                #ax.plot(dVMRy_WRF, dVMRx_WRF, color='r', label='WRF-Chem x '+' '+str(scFactor), markersize=16, linewidth=3)
                #ax.errorbar(yDate,RFTS, yerr=vmrPStdWRF*2.0, fmt='o', markersize=0, color='r', ecolor='r')
                ax.scatter(dVMRy_WRF, dVMRx_WRF, facecolors='white', s=80, color='r')

                ax.grid(True)
                #ax.set_title(gasSTR + ' VMR Weighted between '+str(pCols[0])+' to '+str(pCols[1])+' km',
                #              multialignment='center',fontsize=14)
                ax.tick_params(labelsize=14)
                ax.set_ylabel('Delta '+gasSTR + '  Weighted VMR [ppb$_v$]', fontsize=18)
                ax.set_xlabel('Delta CO  Weighted VMR [ppb$_v$]', fontsize=18)
                ax.legend(prop={'size':16})

                


        #-------------------------
        # Close output figure file
        #-------------------------
 
        if saveFlg: pdfsav.close()
        else:
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()

    
    
    #if saveFlg: pdfsav.close()
    #else:
    #    user_input = raw_input('Press any key to exit >>> ')
    #    sys.exit()

 
if __name__ == "__main__":
    main()