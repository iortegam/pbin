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
                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

def main():
    #-----------------------------------------------------------------------------------------
    #                             Initializations for FTIR
    #-----------------------------------------------------------------------------------------

    loc               = 'fl0'

    #gasName           = ['co', 'c2h6', 'h2co',   'c2h2', 'ch4', 'hcooh', 'nh3',  'o3']             
    #ver               = ['Current_v2', 'Current_v2', 'Current_WP_v6',    'Current_v2', 'Current_WP', 'Current_v1', 'Current_v2',  'Current_WP']          # Name of retrieval version to process
    #ctlF              = ['sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v6.ctl', 'sfit4_v2.ctl', 'sfit4_3.ctl', 'sfit4_v1.ctl', 'sfit4_v2.ctl', 'sfit4.ctl' ]            # Name of ctl file

    #gasName           = ['co', 'c2h6', 'nh3']             
    #ver               = ['Current_v2', 'Current_v2', 'Current_v2' ]          # Name of retrieval version to process
    #ctlF              = ['sfit4_v2.ctl', 'sfit4_v2.ctl',  'sfit4_v2.ctl']            # Name of ctl file

    gasName            = ['co', 'c2h2',   'c2h6',  'ch4', 'nh3']#,   'h2co', 'hcooh',  'nh3', 'o3']             
    ver                = ['Current_v3', 'Current_v2', 'Current_v2',   'Current_WP', 'Current_v2',   'Current_WP_v6','Current_v1',   'Current_WP']          # Name of retrieval version to process
    ctlF               = ['sfit4_v3.ctl',  'sfit4_v2.ctl',   'sfit4_v2.ctl',  'sfit4_3.ctl', 'sfit4_v2.ctl', 'sfit4_v6.ctl', 'sfit4_v1.ctl',  'sfit4.ctl'] 

    saveFlg           = False 
    pltFile           =  '/data/iortega/results/'+loc.lower()+'/fig/'+loc.upper()+'_Results_b.pdf'
 
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
    maxRMS             = [2.0, 0.65, 1.0, 0.3, 1.0, 1.0, 1.0, 4.0]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0, 0.5, 0.5, 0.9, 0.5,  0.5, 0.5, 0.0]                      # Min DOFs for filtering
    #maxRMS             = [2.0, 1.0, 1.0]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    #minDOF             = [1.0, 0.5, 0.5]                      # Min DOFs for filtering

    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    pCols             = [ [1.6, 8.0] ]       #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2016
    fmnth              = 6
    fday               = 30

    #-------------------------------
    # Decide how to create layers
    # above and below the tropopause
    #-------------------------------
    redFct  = 5.0
    hgtTrgU = 10.0       # Upper Layer thickness [km]
    hgtTrgL = False       # Lower Layer thickness [km];Set to False if you want surface value

    weatherFlag     = False
    weatherDir      = '/data1/ancillary_data/fl0/eol/'
    weatherFileTag  = 'v2'

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

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

    vmrP         = OrderedDict()
    TCp           = OrderedDict()

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
        nfltr           = len(statDataCl[gasVer].inds)
        rms[gasVer]     = np.delete(rms[gasVer],statDataCl[gasVer].inds)
        ntot            = len(rms[gasVer])
        dates[gasVer]   = np.delete(dates[gasVer],statDataCl[gasVer].inds)
        totClmn[gasVer] = np.delete(totClmn[gasVer],statDataCl[gasVer].inds)
        rPrfVMR[gasVer] = np.delete(rPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfMol[gasVer] = np.delete(rPrfMol[gasVer],statDataCl[gasVer].inds,axis=0)   
        rPrfDry[gasVer] = np.delete(rPrfDry[gasVer],statDataCl[gasVer].inds,axis=0)   
        Airmass[gasVer] = np.delete(Airmass[gasVer],statDataCl[gasVer].inds,axis=0)
        TCdry[gasVer]   = np.delete(TCdry[gasVer],statDataCl[gasVer].inds)

        aPrfVMR[gasVer] = np.delete(aPrfVMR[gasVer],statDataCl[gasVer].inds,axis=0)   
        aPrfMol[gasVer] = np.delete(aPrfMol[gasVer],statDataCl[gasVer].inds,axis=0) 
    
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


        #-------------------------------------------------
        # Calculate partial columns and weighted VMR
        #-------------------------------------------------

        for pcol in pCols:
            ind1             = mf.nearestind(pcol[0], alt[gasVer])
            ind2             = mf.nearestind(pcol[1], alt[gasVer])
            vmrP[gasVer]     = np.average(rPrfVMR[gasVer][:,ind2:ind1],axis=1,weights=Airmass[gasVer][:,ind2:ind1])           
            TCp[gasVer]      = np.sum(rPrfMol[gasVer][:,ind2:ind1],axis=1)

    
    #-------------------------------------------------
    #-------------------- PLOTS ---------------------
    #-------------------------------------------------
    if iyear == fyear: yrsFlg = False
    else: yrsFlg = True

    pl = mp.myPltClass(gasName, ver, iyear, imnth, iday, fyear, fmnth, fday, outFname=pltFile, saveFlg= saveFlg, yrsFlg=yrsFlg)

    if weatherFlag: wdir, wspeed, temp, rh, dt_weather = mf.weatherout(loc, weatherDir, weatherFileTag, iyear, fyear )
    else: wdir, wspeed, temp, rh, dt_weather = 0., 0., 0., 0., 0.

    #-------------------------------------------------
    # PLOTs: Time Series
    #-------------------------------------------------
    pl.pltts(dates, totClmn, vmrP, TCp, pCols)

    #pl.pltts_2(dates, totClmn, vmrP, TCp, pCols)      #FOR IGAC POSTER
 
    #-------------------------------------------------
    # PLOTs: BOX PLOTS
    #-------------------------------------------------
    #pl.pltbox(dates, vmrP, idtitle = 'Weighted VMR [ppbv]')

    #-------------------------------------------------
    # PLOTs: Time Series of Ratios, Correlation Analysis, and Q-Q plot and Wind if True
    #-------------------------------------------------
    ##pl.pltcorr(dates, totClmn, gasx = 'co', idtitle = 'Total Column [molecules cm$^{-2}$]', weatherFlag=weatherFlag,
    ##           wdir=wdir, wspeed=wspeed, temp=temp,  dt_weather=dt_weather, qqFlag=True)

    pl.pltcorr(dates, vmrP, gasx = 'co', idtitle = 'Weighted VMR [ppbv]', weatherFlag=weatherFlag,
               wdir=wdir, wspeed=wspeed, temp=temp,  dt_weather=dt_weather, qqFlag=False)

    #pl.pltcorr(dates, vmrP, gasx = 'c2h2', idtitle = 'Weighted VMR [ppbv]', weatherFlag=weatherFlag,
    #           wdir=wdir, wspeed=wspeed, temp=temp,  dt_weather=dt_weather, qqFlag=True)

    #pl.pltcorr(dates, vmrP, gasx = 'ch4', idtitle = 'Weighted VMR [ppbv]', weatherFlag=weatherFlag,
    #           wdir=wdir, wspeed=wspeed, temp=temp,  dt_weather=dt_weather, qqFlag=True)

    #pl.pltcorr(dates, vmrP, gasx = 'nh3', idtitle = 'Weighted VMR [ppbv]', weatherFlag=weatherFlag,
    #           wdir=wdir, wspeed=wspeed, temp=temp,  dt_weather=dt_weather, qqFlag=True)


    if saveFlg: pl.closeFig()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()


 
if __name__ == "__main__":
    main()