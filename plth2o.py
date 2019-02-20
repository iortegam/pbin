#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         plth2o.py
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
import dataOutClass as dc
from collections                     import OrderedDict
import PltClass as mp

sys.path.append('/data/iortega/pbin/sondeFTIR')

import sondeclass as rs

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

    loc        = 'fl0'

    gasName    = ['h2o']             
    ver        = ['Current_ERA']
    #ver        = ['Current']          
    ctlF       = ['sfit4.ctl'] 

    saveFlg    = True 
    pltFile    =  '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+'_Results_H2O.pdf'
    FigDir     = '/data/iortega/Manuscripts/Fig/'


    #------
    # Flags
    #------
    errorFlg   = True                  # Flag to process error data
    fltrFlg    = True                   # Flag to filter the data
    byYrFlg    = False                  # Flag to create plots for each individual year in date range
    szaFlg     = True                   # Flag to filter based on min and max SZA
    dofFlg     = True                   # Flag to filter based on min DOFs
    pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
    tcNegFlg   = True                  # Flagsag to filter profiles with negative total columns
    tcMMFlg    = False                  # Flag to filter based on min and max total column amount
    cnvrgFlg   = False                   # Flag to filter profiles that did not converge
    rmsFlg     = True                   # Flag to filter based on max RMS
    chiFlg     = False                   # Flag to filter based on max CHI_2_Y
    mnthFlg    = False                  # Flag to filter based on 

    mnths      = [6,7,8]                # Months to filter on (these are the months to include data)
    maxRMS     = [3.0]                         # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.0]                   # Min DOFs for filtering
    minSZA     = 0.0                   # Min SZA for filtering
    maxSZA     = 90.0                   # Max SZA for filtering
    maxCHI     = 2.0                    # Max CHI_y_2 value
    maxTC      = 5.0E24                # Max Total column amount for filtering
    minTC      = 0.0                 # Min Total column amount for filtering
    sclfct     = 1.0E6                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName = 'ppmv'                 # Name of scale factor for labeling plots

    #----------------------
    # Date range to process
    #----------------------
    if   loc.lower() == 'fl0':   iyear     = 2010
    elif loc.lower() == 'mlo':   iyear     = 1995
    elif loc.lower() == 'tab':   iyear     = 1999
    else:
        print 'Error: check site name!'
        sys.exit()

    imnth      = 1
    iday       = 1
    fyear      = 2016
    fmnth      = 12
    fday       = 31

    if loc.lower() == 'fl0': pCols      = [ [1.5, 3.0], [3.0, 5.0], [5.0, 7.5], [7.5, 10.0] ]       #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR
    if loc.lower() == 'mlo': pCols      = [ [3.0, 5.0], [5.0, 7.0], [7.0, 9.0], [9.0, 11.0]]
    if loc.lower() == 'tab': pCols      = [ [0.0, 3.0], [3.0, 5.0], [5.0, 7.0], [7.0, 9.0] ] 

    #pCols = [ [1.5, 3.0], [3.0, 5.0], [5.0, 7.5], [7.5, 10.0], [10.0,13.0], [13.0,17.0], [17.0,21.0] ]


#----------------------------------------------------------------------------------------------------------------------

                                    #----------------------------#
                                    #     --- START SONDE ---    #
                                    #----------------------------#
    
    if   loc.lower() != 'tab':
        sondeDir     = '/data1/ancillary_data/sondes/'+loc.lower()+'/'

        sn = rs.ReadOutputsSonde(sondeDir, loc, iyear=2010,imnth=imnth,iday=iday,fyear=fyear,fmnth=fmnth,fday=fday)

        sn.readfilessonde()
        sn.sondePrf()

        sondevmrP      = {}
        sondevmrsdP    = {}
        sondesumP      = {}
        sondesumsdP    = {}
        sondedates     = {}

        for pn, pcol in enumerate(pCols):


            for d, da in enumerate(sn.sonde['date']):

                sondeH2OPrf_a       = np.asarray(sn.sonde['H2Omr_ppmv_a'][d])
                sondeH2OPrfsd_a     = np.asarray(sn.sonde['H2Osd_ppmv_a'][d])
                sondeH2OPrfMol_a    = np.asarray(sn.sonde['H2OMol_a'][d])
                sondeH2OPrfMolsd_a  = np.asarray(sn.sonde['H2OMolsd_a'][d])

                sondealt_a          = np.asarray(sn.sonde['Alt_km_a'][d])
                sondeairmass_a      = np.asarray(sn.sonde['airmass_a'][d])

                ind1 = mf.nearestind(pcol[0], sondealt_a)
                ind2 = mf.nearestind(pcol[1], sondealt_a) 

                sondeH2OPrf_a_i    = np.asarray(sondeH2OPrf_a, dtype=np.float32)
                sondeH2OPrfsd_a_i  = np.asarray(sondeH2OPrfsd_a, dtype=np.float32)
                sondeairmass_a_i   = np.asarray(sondeairmass_a, dtype=np.float32)
                sondeH2OPrfMol_a_i = np.asarray(sondeH2OPrfMol_a, dtype=np.float32)
                sondeH2OPrfMolsd_a_i = np.asarray(sondeH2OPrfMolsd_a, dtype=np.float32)

                sondealt_a  = map(lambda x: float(x),sondealt_a)
                sondealt_a  = np.asarray(sondealt_a)
               
                if ind1 == ind2:
                    continue

                sondevmrP_i = np.average(sondeH2OPrf_a_i[ind1:ind2],  weights=sondeairmass_a_i[ind1:ind2]) 
                sondevmrP.setdefault(pn,[]).append(sondevmrP_i)

                sondevmrsdP_i = np.average(sondeH2OPrfsd_a_i[ind1:ind2], weights=sondeairmass_a_i[ind1:ind2])
                sondevmrsdP.setdefault(pn,[]).append(sondevmrsdP_i)

                sondesumP_i = np.sum(sondeH2OPrf_a_i[ind1:ind2]*sondeairmass_a_i[ind1:ind2]*1e-6)
                sondesumP.setdefault(pn,[]).append(sondesumP_i)

                sondesumsdP_i = np.sum(sondeH2OPrfsd_a_i[ind1:ind2]*sondeairmass_a_i[ind1:ind2]*1e-6)
                sondesumsdP.setdefault(pn,[]).append(sondesumsdP_i)

                sondedates.setdefault(pn,[]).append(da)

   
    #----------------------------------------------------------------------------------------------------------------------


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

    rPrfVMR          = OrderedDict()
    rPrfMol          = OrderedDict()
    dates            = OrderedDict()
    alt              = OrderedDict()
    Airmass          = OrderedDict()
    waterVMR         = OrderedDict()
    waterMol         = OrderedDict()
    totClmn          = OrderedDict()
    TCdryAir         = OrderedDict()
    TCdry            = OrderedDict()
    rPrfDry          = OrderedDict()
    rms              = OrderedDict()
    dofsAvg          = OrderedDict()
    dofsAvg_cs       = OrderedDict()
    LayVMR           = OrderedDict()
    DryLayVMR        = OrderedDict()
    upperHgt         = OrderedDict()
    lowerHgt         = OrderedDict()
    LayThk           = OrderedDict()
    LayDOF           = OrderedDict()

    vmrP             = {}    #OrderedDict()
    TCp              = {}    #OrderedDict()

    avkSCF           = OrderedDict()
    avkVMR           = OrderedDict()
    avkSCFav         = OrderedDict()
    avkVMRav         = OrderedDict()

    aPrfVMR          = OrderedDict()
    aPrfMol          = OrderedDict()
           
    rand_err         = OrderedDict()
    sys_err          = OrderedDict()
    tot_err          = OrderedDict()

    rand_errvmr      = OrderedDict()
    sys_errvmr       = OrderedDict()
    tot_errvmr       = OrderedDict()

    rand_cmpnts      = OrderedDict() 
    sys_cmpnts       = OrderedDict()

    rand_cmpnts_vmr  = OrderedDict()
    sys_cmpnts_vmr   = OrderedDict()

    tot_rnd          = OrderedDict()
    tot_sys          = OrderedDict()
    tot_std          = OrderedDict()

    err_summary      = OrderedDict()


    for j, gasVer in enumerate(statDataCl):
        rPrfVMR[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]) * sclfct
        rPrfMol[gasVer]  = np.asarray(statDataCl[gasVer].rprfs[statDataCl[gasVer].PrimaryGas]  * np.asarray(statDataCl[gasVer].rprfs['AIRMASS']))
        dates[gasVer]    = np.asarray(statDataCl[gasVer].rprfs['date'])
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

        #----------------------------------
        # Read readPbp (observed, fitted, and difference spectra)
        #----------------------------------
        mw    = [str(int(x)) for x in statDataCl[gasVer].ctl['band']]     
        numMW = len(mw)
        statDataCl[gasVer].readPbp()

        #----------------------------------
        # Read Spectra for each gas
        #----------------------------------
        statDataCl[gasVer].readSpectra(statDataCl[gasVer].gasList)

        #-------------------- 
        # Call to filter data
        #--------------------
        if fltrFlg: statDataCl[gasVer].fltrData(statDataCl[gasVer].PrimaryGas,mxrms=maxRMS[j],rmsFlg=rmsFlg, minDOF=minDOF[j],  dofFlg=dofFlg, 
        tcFlg=tcNegFlg, pcFlg=pcNegFlg , cnvrgFlg=cnvrgFlg)
        else:       statDataCl[gasVer].inds = np.array([]) 

        #--------------------------------------------
        # Read Error data to get AVK and profile DOFs
        #-------------------------------------------------
        # Determine if AVK has been created via sfit4 core
        # code or via python error analysis
        #-------------------------------------------------     
        if errorFlg:   # Read AVK from error output
            #statDataCl[gasVer].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
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
            tot_errvmr[gasVer]  = np.sqrt(rand_errvmr[gasVer] + sys_errvmr[gasVer]) 

            #-------------------------------------------------
            #Error profiles of components
            #-------------------------------------------------
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
            tot_std[gasVer]        = np.sqrt(tot_rnd[gasVer]**2 + tot_sys[gasVer]**2)
            
            err_summary[gasVer]    = statDataCl[gasVer].error

            #-------------------------------------------------
            #Error profiles of the layers of interes constructed with the elements of the covariances
            #-------------------------------------------------
            # npntslay       = np.shape(pCols)[0]
            # rand_err_lay   = np.zeros((npnts,npntslay))
            # sys_err_lay    = np.zeros((npnts,npntslay))

            # for i in range(npnts):
            #     for pi, p in enumerate(pCols):
            #         ind1             = mf.nearestind(p[0], alt[gasVer])
            #         ind2             = mf.nearestind(p[1], alt[gasVer])

            #         rand_err_lay[i,pi] = statDataCl[gasVer].error['Total_Random_Error'][i][ind2+1:ind1+1, ind2+1:ind1+1].sum()
            #         sys_err_lay[i,pi]  = statDataCl[gasVer].error['Total_Systematic_Error'][i][ind2+1:ind1+1, ind2+1:ind1+1].sum()

            # tot_err_lay  = np.sqrt(rand_err_lay + sys_err_lay)            
            # rand_err_lay = np.sqrt(rand_err_lay)
            # sys_err_lay  = np.sqrt(sys_err_lay)

            # alt_lay   = np.zeros((npntslay))

            # for pi, p in enumerate(pCols):
            #     alt_lay[pi] = (p[0] + p[1])*0.5

            #-------------------------------------------------
            #Error profiles of components
            #-------------------------------------------------
            # rand_cmpnts_diag = statDataCl[gasVer].randErrDiag
            # sys_cmpnts_diag  = statDataCl[gasVer].sysErrDiag

            # sys_err_comp_lay  = {}

            # #rand_err_i = statDataCl[gasVer].sysErr['linepair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error']

            # fig, ax1       = plt.subplots()

            # for k in sys_cmpnts_diag:
            #     errPlt = np.nanmean(np.sqrt(sys_cmpnts_diag[k]),axis=0)
                
            #     ax1.plot(errPlt,alt[gasVer],linewidth=0.75,linestyle='--', label=k+'_diag')
                
            #     sys_err_comp_i   = statDataCl[gasVer].sysErr[k]

            #     sys_err_comp_lay[k] = np.zeros((npnts,npntslay))

            #     for i in range(npnts):

            #         for pi, p in enumerate(pCols):
            #             ind1             = mf.nearestind(p[0], alt[gasVer])
            #             ind2             = mf.nearestind(p[1], alt[gasVer])

            #             sys_err_comp_lay[k][i,pi] = np.sqrt(sys_err_comp_i[i][ind2+1:ind1+1, ind2+1:ind1+1].sum())

            #     errPlt2 = np.nanmean(sys_err_comp_lay[k],axis=0)
            #     ax1.plot(errPlt2,alt_lay,linewidth=1.75, label=k+'_lay')


            # ax1.plot(np.nanmean(rand_err_lay, axis=0), alt_lay,linestyle='--',linewidth=2.0, label='Ran Err Layer')
            # ax1.plot(np.nanmean(sys_err_lay, axis=0), alt_lay,linestyle='--', linewidth=2.0, label='Sys Err Layer')
            # ax1.plot(np.nanmean(tot_err_lay, axis=0), alt_lay,linestyle='--',linewidth=2.0, label='Tot Err Layer')

            # ax1.plot(np.nanmean(rand_err_diag, axis=0), alt[gasVer],linewidth=2.0, label='Ran Err diag')
            # ax1.plot(np.nanmean(sys_err_diag, axis=0), alt[gasVer], linewidth=2.0, label='Sys Err diag')
            # ax1.plot(np.nanmean(tot_err_diag, axis=0), alt[gasVer],linewidth=2.0, label='Tot Err diag')
            
            # ax1.set_ylabel('Altitude [km]', fontsize=14)
            # ax1.set_xlabel('VMR [ppm]', fontsize=14)             
            # ax1.grid(True,which='both')
            # ax1.legend(prop={'size':12},loc='upper right')          
            # #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
            # ax1.set_title('Random Components')
            # ax1.set_ylim(1, 15)
            # ax1.tick_params(which='both',labelsize=14)
            # plt.show(block=False)
            # user_input = raw_input('Press any key to exit >>> ')
            # sys.exit()


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

        if errorFlg:
            rand_err[gasVer] = np.delete(rand_err[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_err[gasVer]  = np.delete(sys_err[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_err[gasVer]  = np.delete(tot_err[gasVer],statDataCl[gasVer].inds,axis=0)

            rand_errvmr[gasVer] = np.delete(rand_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)
            sys_errvmr[gasVer]  = np.delete(sys_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)  
            tot_errvmr[gasVer]  = np.delete(tot_errvmr[gasVer],statDataCl[gasVer].inds,axis=0)

            tot_rnd[gasVer] = np.delete(tot_rnd[gasVer],statDataCl[gasVer].inds)
            tot_sys[gasVer] = np.delete(tot_sys,statDataCl[gasVer].inds)
            tot_std[gasVer] = np.delete(tot_std,statDataCl[gasVer].inds)

            
            for k in sys_cmpnts[gasVer]:
                sys_cmpnts[gasVer][k]      = np.delete(sys_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                sys_cmpnts_vmr[gasVer][k]  = np.delete(sys_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)
                
            for k in rand_cmpnts[gasVer]:
                rand_cmpnts[gasVer][k]     = np.delete(rand_cmpnts[gasVer][k],statDataCl[gasVer].inds,axis=0)
                rand_cmpnts_vmr[gasVer][k] = np.delete(rand_cmpnts_vmr[gasVer][k],statDataCl[gasVer].inds,axis=0)

            for k in err_summary[gasVer]:
                err_summary[gasVer][k]     = np.delete(err_summary[gasVer][k],statDataCl[gasVer].inds, axis=0)
          

        #--------------------------------------
        # Remove retrieval data based on filter for spectra fit
        #--------------------------------------
        dataSpec  = OrderedDict()
        gasSpec   = OrderedDict()
        gasAbs    = OrderedDict()
        gasAbsSNR = OrderedDict()
        mwList    = OrderedDict()

        for x in mw:  # Loop through micro-windows
            dataSpec['Obs_'+x]        = np.delete(statDataCl[gasVer].pbp['Obs_'+x],statDataCl[gasVer].inds,axis=0)
            dataSpec['Fitted_'+x]     = np.delete(statDataCl[gasVer].pbp['Fitted_'+x],statDataCl[gasVer].inds,axis=0)
            dataSpec['Difference_'+x] = np.delete(statDataCl[gasVer].pbp['Difference_'+x], statDataCl[gasVer].inds, axis=0)
            dataSpec['WaveN_'+x]      = statDataCl[gasVer].spc['MW_'+x]
            dataSpec['All_'+x]        = np.delete(statDataCl[gasVer].spc['All_'+x],statDataCl[gasVer].inds,axis=0)
            if statDataCl[gasVer].solarFlg: dataSpec['Sol_'+x]        = np.delete(statDataCl[gasVer].spc['Solar_'+x],statDataCl[gasVer].inds,axis=0)

        #----------------------------------------
        # Loop through gases and micro-windows
        # Not all gasses are in each micro-window
        # Determine which gasses belong to which
        # micro-window
        #----------------------------------------
                # Use ordered dictionary so that micro-window plots print in order
        for x in mw: mwList[x] = statDataCl[gasVer].ctl['band.'+x+'.gasb']
        
        gasSpec = {gas.upper()+'_'+m:np.delete(statDataCl[gasVer].spc[gas.upper()+'_'+m],statDataCl[gasVer].inds,axis=0) for m in mwList for gas in mwList[m]}

        #---------------------
        # Calculate Statistics
        #---------------------
        for x in mw:  # Loop through micro-windows

            #-------------------------------------------------------------
            # Calculate the total Observed absorption in micro-window
            # This must be done first because below code modifies dataSpec
            #-------------------------------------------------------------
            gasAbs["Total_"+x] = simps(1.0 - dataSpec['Obs_'+x],x=dataSpec['WaveN_'+x],axis=1)               
                        
            if len(dates[gasVer]) > 1:
                dataSpec['Obs_'+x]        = np.mean(dataSpec['Obs_'+x],axis=0)
                dataSpec['Fitted_'+x]     = np.mean(dataSpec['Fitted_'+x],axis=0)
                dataSpec['Difference_'+x] = np.mean(dataSpec['Difference_'+x],axis=0)
                dataSpec['All_'+x]        = np.mean(dataSpec['All_'+x],axis=0)
                if statDataCl[gasVer].solarFlg:dataSpec['Sol_'+x]        = np.mean(dataSpec['Sol_'+x],axis=0)    
                dataSpec['DifSTD_'+x]     = np.std(dataSpec['Difference_'+x],axis=0)
                
            else:
                dataSpec['Obs_'+x]        = dataSpec['Obs_'+x][0]
                dataSpec['Fitted_'+x]     = dataSpec['Fitted_'+x][0]
                dataSpec['Difference_'+x] = dataSpec['Difference_'+x][0]
                dataSpec['All_'+x]        = dataSpec['All_'+x][0]
                if statDataCl[gasVer].solarFlg:dataSpec['Sol_'+x]        = dataSpec['Sol_'+x][0]
                                
            if len(dates[gasVer]) > 1:

                if statDataCl[gasVer].PrimaryGas+"_"+x in gasSpec:  
                #---------------------------------------------------
                # Calculate the integrate absorption for primary gas
                #---------------------------------------------------               
                    gasAbs[statDataCl[gasVer].PrimaryGas+"_"+x] = simps(1.0 - gasSpec[statDataCl[gasVer].PrimaryGas+"_"+x],x=dataSpec['WaveN_'+x],axis=1)       
                
                #-----------------------------------
                # Calculate the peak absorption of 
                # primary gas for each micro-window
                #-----------------------------------
                    gasAbs[statDataCl[gasVer].PrimaryGas+"_trans_"+x] = 1.0 - np.min(gasSpec[statDataCl[gasVer].PrimaryGas+"_"+x],axis=1)
                
                #---------------------------------------------
                # Determine product of SNR and Peak Absorption
                #---------------------------------------------
                    tempSNR  = np.delete(statDataCl[gasVer].summary["SNR_"+x],statDataCl[gasVer].inds)
                    gasAbsSNR[statDataCl[gasVer].PrimaryGas+"_"+x] = gasAbs[statDataCl[gasVer].PrimaryGas+"_"+x] * tempSNR

        if len(dates[gasVer]) > 1:
            gasSpec = {gas.upper()+'_'+x:np.mean(gasSpec[gas.upper()+'_'+x],axis=0) for x in mwList for gas in mwList[x]}   
        else:
            for x in gasSpec: gasSpec[x] = gasSpec[x][0]

        
        #-------------------------------------------------
        # Calculate partial columns and weighted VMR
        #-------------------------------------------------
        for p, pcol in enumerate(pCols):
            
            ind1             = mf.nearestind(pcol[0], alt[gasVer])
            ind2             = mf.nearestind(pcol[1], alt[gasVer])
            #vmrP[gasVer][p]  = np.average(rPrfVMR[gasVer][:,ind2:ind1],axis=1,weights=Airmass[gasVer][:,ind2:ind1])           
            #TCp[gasVer][p]   = np.sum(rPrfMol[gasVer][:,ind2:ind1],axis=1)
            vmrP.setdefault(gasVer+str(p) ,[]).append(np.average(rPrfVMR[gasVer][:,ind2:ind1],axis=1,weights=Airmass[gasVer][:,ind2:ind1]))         
            TCp.setdefault(gasVer+str(p), []).append(np.sum(rPrfMol[gasVer][:,ind2:ind1],axis=1))

        for p, pcol in enumerate(pCols):
            vmrP[gasVer+str(p)] = np.asarray(vmrP[gasVer+str(p)])
            TCp[gasVer+str(p)] = np.asarray(TCp[gasVer+str(p)])
    


    #-------------------------------------------------
    #                     PLOTS  
    #-------------------------------------------------
    if saveFlg: pdfsav = PdfPages(pltFile)
   
    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)              
    yearsLc      = YearLocator()
    monthsAll    = MonthLocator()
    #months       = MonthLocator(bymonth=1,bymonthday=1)
    monthsFmt     = MonthLocator()
    DateFmt      = DateFormatter('%Y')


    for j, gasVer in enumerate(statDataCl):

        #--------------------------------
        # VERTICAL PROFILES (CONTOUR PLOTS)
        #--------------------------------
        indsH = np.where(alt[gasVer] <= 11.0)[0]
        #levels = np.arange(np.round(np.min(rPrfVMR[gasVer][:,indsH]),decimals=0),np.round(np.max(rPrfVMR[gasVer][:,indsH]),decimals=0))
        levels = np.arange(0,np.round(np.max(rPrfVMR[gasVer][:,indsH]/1e3),decimals=1), 0.5)
       
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12,6), sharey=False)
       
        gs   = gridspec.GridSpec(1, 2, width_ratios=[3,1])
        ax  = plt.subplot(gs[0])
        ax2  = plt.subplot(gs[1])

        cax          = ax.contourf(dates[gasVer],alt[gasVer][indsH],np.transpose(rPrfVMR[gasVer][:,indsH])/1e3,levels, cmap=mplcm.jet)

        divider      = make_axes_locatable(ax)
        cb           = divider.append_axes("right",size="3.5%",pad=0.05)
      
        cbar         = plt.colorbar(cax,cax=cb)
        #cbar.ax.tick_params(labelsize=12)
        cbar.set_label('VMR [x10$^3$ ppm]', fontsize=12)#, fontsize=14)
        
        ax.set_xlabel('Year', fontsize=12)#, fontsize=14)
        ax.set_ylabel('Altitude [km]', fontsize=12)#, fontsize=14)
        ax.set_ylim(np.min(alt[gasVer][indsH]), 10)

        #ax.yaxis.set_tick_params(which='major',labelsize=8)
        #ax.xaxis.set_tick_params(which='major',labelsize=8)
        ax.xaxis.set_major_locator(yearsLc)
        ax.xaxis.set_major_formatter(DateFmt)
        ax.xaxis.set_minor_locator(monthsFmt)
        ax.set_title('Time Series of Retrieved H$_2$O Vertical Profiles',multialignment='center')
        ax.xaxis.set_tick_params(which='major',labelsize=14)
        ax.xaxis.set_tick_params(which='minor',labelbottom='off')
        ax.tick_params(labelsize=14)

        #----
        months   = np.asarray([d.month for d in dates[gasVer]])  

        for m in list(set(months)):
            if m == 1: lm = 'Jan'
            if m == 2: lm = 'Feb'
            if m == 3: lm = 'Mar'
            if m == 4: lm = 'Apr'
            if m == 5: lm = 'May'
            if m == 6: lm = 'Jun'
            if m == 7: lm = 'Jul'
            if m == 8: lm = 'Aug'
            if m == 9: lm = 'Sep'
            if m == 10: lm = 'Oct'
            if m == 11: lm = 'Nov'
            if m == 12: lm = 'Dec'
            inds     = np.where(months == m)[0]
            mnthMean = np.mean(rPrfVMR[gasVer][inds,:], axis=0)

            ax2.plot(mnthMean[indsH]/1e3,alt[gasVer][indsH],label=lm)

        #ax2.tick_params(labelsize=14)
        ax2.set_xlabel('VMR [x10$^3$ ppm]', fontsize=12)#, fontsize=14)
        #ax2.axes.get_yaxis().set_visible(False)
        ax2.set_yticklabels([])
        ax2.grid(True)
        ax2.legend(prop={'size':10})
        ax2.set_ylim(np.min(alt[gasVer][indsH]), 10)
        ax2.set_title('Monthly Mean',multialignment='center')
        #ax2.yaxis.set_major_formatter(ax2.NullFormatter())
        #ax2.yticks([])
        
        #----
        fig.subplots_adjust(left=0.075, bottom=0.09,top=0.95, right=0.975)
        
        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(FigDir+'TSPrF_'+loc.upper()+'.pdf', bbox_inches='tight')

        else: 
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        fig1,ax1 = plt.subplots()
        ax1.plot(dates[gasVer],totClmn[gasVer]*2.989e-22,'k.',markersize=6)
        ax1.grid(True)
        ax1.set_ylabel('Retrieved Total Column [mm]',multialignment='center')
        ax1.set_xlabel('Year')
        ax1.set_title('Time Series of Retrieved H$_2$O Total Column [mm]',multialignment='center')
        
        
        ax1.xaxis.set_major_locator(yearsLc)
        ax1.xaxis.set_minor_locator(monthsFmt)
        #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
        ax1.xaxis.set_major_formatter(DateFmt) 
        #ax1.xaxis.set_tick_params(which='major', pad=15)  
        ax1.xaxis.set_tick_params(which='major',labelsize=8)
        ax1.xaxis.set_tick_params(which='minor',labelbottom='off')
        
        
        if saveFlg: pdfsav.savefig(fig1,dpi=200)
        else:       plt.show(block=False)  

        #-------------------------------------------------
        #Fit: Mean   
        #-------------------------------------------------
        fig = plt.figure(figsize=(16,9))
            
        for p, x in enumerate(mwList):
            gs1 = gridspec.GridSpec(3, 2)

            if p == 0: gs1.update(left=0.05, right=0.48, top=0.95, bottom=0.55,  wspace=0.05,  hspace=0.08)
            if p == 1: gs1.update(left=0.05, right=0.48, top = 0.45, bottom=0.1, wspace=0.05, hspace=0.08)
            if p == 2: gs1.update(left=0.55, right=0.98, top=0.95, bottom=0.55,  wspace=0.05,  hspace=0.08)
            if p == 3: gs1.update(left=0.55, right=0.98, top=0.45, bottom=0.1, wspace=0.05,   hspace=0.08)

            #ax1 = fig.add_subplot(gs1[0])
            #ax2 = fig.add_subplot(gs1[1], sharex=ax1)

            ax1 = plt.subplot(gs1[0:1, :])
            ax2 = plt.subplot(gs1[1:3, :])

           
            ax1.plot(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*100.0,  color='k')
            ax1.axhline(y=0,color='r')
            if len(dates[gasVer]) > 1: ax1.fill_between(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*100.0-dataSpec['DifSTD_'+x]*100.0,dataSpec['Difference_'+x]*100.0+dataSpec['DifSTD_'+x]*100.0,alpha=0.5,color='0.75')
            ax1.grid(True, alpha=0.5)
            if (p == 0) or (p == 1): ax1.set_ylabel('% Difference', fontsize=14,)
            ax1.set_title('Micro-Window '+ x)
            ax1.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))
            
            plt.tick_params(which='minor',length=4,color='b')
            
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Obs_'+x],label='Observed')
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Fitted_'+x],label='Fitted')
            if statDataCl[gasVer].solarFlg: ax2.plot(dataSpec['WaveN_'+x],dataSpec['Sol_'+x],label='Solar')
            sclfct = 0.0
            for g in mwList[x]: 
                sclfct += 0.03
                
                ax2.plot(dataSpec['WaveN_'+x],gasSpec[g.upper()+'_'+x]+(gasSpec[g.upper()+'_'+x]*sclfct),label=g.upper())
            
            ax2.grid(True)
            if (p == 1) or (p == 3): ax2.set_xlabel('Wavenumber [cm$^{-1}$]', fontsize=14)
            if (p == 0) or (p == 1): ax2.set_ylabel('Transmission', fontsize=14)
            #ax2.set_ylim(bottom=0.0)
            if (p == 0) or (p == 1):ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))

            ax2.legend(prop={'size':12},loc='upper center', bbox_to_anchor=(0.5, 1.065),
                      fancybox=True, ncol=len(mwList[x])+3)  
            
            ax1.tick_params(axis='x',which='both',labelsize=14, bottom='off', labelbottom='off')
            ax1.tick_params(axis='y',which='both',labelsize=14)
            ax2.tick_params(axis='x',which='both',labelsize=14)
            ax2.tick_params(axis='y',which='both',labelsize=14)

            ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax1.xaxis.set_minor_locator(AutoMinorLocator())

            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())

            ax1.set_ylim(-2,2)
            ax2.set_ylim(0,1.3)

            #major_ticks = np.arange(-2, 2, 1.0)
            #ax1.set_xticks(major_ticks)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(FigDir+'FitMean1_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)


        #-------------------------------------------------
        #Fit: Mean   
        #-------------------------------------------------

        clr = ['maroon', 'gray', 'orange', 'cyan', 'yellow', 'black']
        #clr = ('green', 'red', 'maroon', 'blue', 'gray', 'yellow', 'cyan', 'orange', 'black')
        fig = plt.figure(figsize=(16,9))

            
        for p, x in enumerate(mwList):
          
            gs1 = gridspec.GridSpec(1, 2)

            if p == 0: gs1.update(left=0.06, right=0.49, top=0.95, bottom=0.55)
            if p == 1: gs1.update(left=0.06, right=0.49, top=0.45, bottom=0.1)
            if p == 2: gs1.update(left=0.53, right=0.97, top=0.95, bottom=0.55)
            if p == 3: gs1.update(left=0.53, right=0.97, top=0.45, bottom=0.1)

            ax2 = plt.subplot(gs1[0:1, :])

            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*10.0,  color='k', linewidth=1.5)
            ax2.axhline(y=0, color='k', linestyle='--')
     
            if len(dates[gasVer]) > 1: ax2.fill_between(dataSpec['WaveN_'+x],dataSpec['Difference_'+x]*10.0-dataSpec['DifSTD_'+x]*10.0,dataSpec['Difference_'+x]*10.0+dataSpec['DifSTD_'+x]*10.0,alpha=0.5,color='0.75')

            if statDataCl[gasVer].solarFlg: ax2.plot(dataSpec['WaveN_'+x],dataSpec['Sol_'+x],label='Solar', color='red')
            sclfct = 0.0
            for i,g in enumerate(mwList[x]): 
                sclfct += 0.05

                if g.upper() == 'CH4': 
                    lab = 'CH$_4$'
                    clr = 'maroon'
                elif g.upper() == 'CO2': 
                    lab = 'CO$_2$'
                    clr = 'cyan'
                elif g.upper() == 'N2O': 
                    lab = 'N$_2$O'
                    clr = 'orange'
                elif g.upper() == 'H2O': 
                    lab = 'H$_2$O'
                    clr = 'gray'
                elif g.upper() == 'HCL': 
                    lab = 'HCl'
                    clr = 'yellow'
                elif g.upper() == 'HDO': 
                    lab = 'HDO'
                    clr = 'purple'
                
                else: lab = g.upper() 
                
                ax2.plot(dataSpec['WaveN_'+x],gasSpec[g.upper()+'_'+x]+(gasSpec[g.upper()+'_'+x]*sclfct),label=lab, color=clr)
            
            plt.tick_params(which='minor',length=4,color='b')
            
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Obs_'+x],label='Observed', color='blue')
            ax2.plot(dataSpec['WaveN_'+x],dataSpec['Fitted_'+x],label='Fitted', color='green')

            ax2.legend(prop={'size':11.75},loc='upper center', bbox_to_anchor=(0.5, 1.065),
                      fancybox=True, ncol=len(mwList[x])+3, handletextpad=0.1)  
            
            
            #ax2.grid(True, alpha=0.1)
            if (p == 1) or (p == 3): ax2.set_xlabel('Wavenumber [cm$^{-1}$]', fontsize=18)
            if (p == 0) or (p == 1): ax2.set_ylabel('Transmission', fontsize=18)
            if (p == 0) or (p == 1):ax2.set_xlim((np.min(dataSpec['WaveN_'+x]),np.max(dataSpec['WaveN_'+x])))

            #ax2.legend(prop={'size':10.5},loc='upper center', bbox_to_anchor=(0.5, 1.065),
            #          fancybox=True, ncol=len(mwList[x])+3)  
            
            ax2.tick_params(axis='x',which='both',labelsize=18)
            ax2.tick_params(axis='y',which='both',labelsize=18)

            ax2.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())
            #ax2.set_xticks(ax2.get_xticks()[::2])

            ax2.xaxis.set_major_locator(plt.MaxNLocator(5))
            ax.axis('equal')

            #xmarks=[i for i in range(dataSpec['WaveN_'+x].min(), dataSpec['WaveN_'+x].max(), 5)]

            #plt.xticks(xmarks)

            ax2.annotate('x10', xy=(0.65, 0.12),size=16, xycoords='axes fraction', xytext=(0.7, 0.25), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.05))

            ax2.set_ylim(-0.22,1.3)



        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(FigDir+'FitMean2_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()


        #-----------------------------
        #ERRORS:  
        #-----------------------------
        if errorFlg:
            print '\nMean Total Column = {}'.format(np.mean(totClmn[gasVer]))
            #-----------------------------
            #ERRORS: Calculate mean errors
            #-----------------------------
            #fig, ax1 = plt.subplots()
            for k in sorted(err_summary[gasVer].iterkeys()):
                try:
                    cmpnts_error_i = np.mean(err_summary[gasVer][k])
                    print '{:<45s} = {:15.5E} '.format(k, cmpnts_error_i)
                    print '{:<45s} = {:10.4f} %'.format(k, cmpnts_error_i/np.mean(totClmn[gasVer]) * 100.)
                    #ax1.plot(dates[gasVer],cmpnts_error[gasVer][k], label=k)
                    #if k == 'Total systematic uncertainty lineint_h2o':
             #           ax1.scatter(dates[gasVer], err_summary[gasVer][k], color='r', label=k)
                    #elif k == 'Total systematic uncertainty linepair_h2o':
              #          ax1.scatter(dates[gasVer], err_summary[gasVer][k], color='blue', label=k)
                except Exception as errmsg:
                    print errmsg
            print '\n'
            #ax1.legend(prop={'size':8},loc='upper right') 
            #plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

            #-----------------------------
            #Plot Random error components
            #-----------------------------
            
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
                elif p == 'curvature_Random_Error': continue
                elif p == 'max_opd_Random_Error': continue     
                else: legend = p
            
                ax1[0].plot(errPlt,alt[gasVer],linewidth=2.0, label=legend)

            ax1[0].plot(np.nanmean(rand_errvmr[gasVer], axis=0)/ retPrf *100.,alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
            ax1[0].set_ylabel('Altitude [km]', fontsize=18)
            #ax1[0].set_xlabel('VMR [ppm]', fontsize=14)
            ax1[0].set_xlabel('uncertainty [%]', fontsize=18)             
            #ax1[0].grid(True,which='both')
            ax1[0].legend(prop={'size':12},loc='upper right')          
            #ax2[k,0].tick_params(axis='both',which='both',labelsize=8) 
            ax1[0].set_title('Random Components', fontsize=18)
            if loc.lower() == 'fl0': ax1[0].set_ylim(1, 15.5)
            if loc.lower() == 'mlo': ax1[0].set_ylim(3, 15.5)
            ax1[0].tick_params(which='both',labelsize=18)
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
                elif p == 'smoothing_Systematic_Error': continue # legend = 'Smoothing'
                #elif p == 'smoothing_Systematic_Error': legend = 'Smoothing'            
                elif p == 'linepair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$\gamma$'
                elif p == 'lineint_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$S$'
                elif p == 'linetair_'+statDataCl[gasVer].PrimaryGas.lower()+'_Systematic_Error': legend = '$n$'
                elif p == 'phase_fcn_Systematic_Error': legend = 'Phase fcn'
                elif p == 'apod_fcn_Systematic_Error':  legend = 'Apod fcn'
                else: legend = p

                ax1[1].plot(errPlt, alt[gasVer],linewidth=2.0, label=legend)

            ax1[1].plot(np.nanmean(sys_errvmr[gasVer], axis=0)/ retPrf *100., alt[gasVer],linestyle='--', color='k',linewidth=2.0, label='Total')
            #ax1[1].set_xlabel('VMR [ppm]', fontsize=14)
            ax1[1].set_xlabel('uncertainty [%]', fontsize=18)             
            #ax1[1].grid(True,which='both')
            ax1[1].legend(prop={'size':12}, loc='upper right')
            ax1[1].set_title('Systematic Components', fontsize=18)
            if loc.lower() == 'fl0': ax1[1].set_ylim(1, 15.5)
            if loc.lower() == 'mlo': ax1[1].set_ylim(3, 15.5)
            ax1[1].tick_params(which='both',labelsize=18)
            ax1[1].set_xlim(0, 10)
                #ax1[1].set_xscale('log')

            fig.subplots_adjust(bottom=0.1, top=0.95, left = 0.1, right = 0.95)
           
            if saveFlg: 
                pdfsav.savefig(fig,dpi=200)
                plt.savefig(FigDir+'UncertaintyPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
            else:       
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()


        #---------------------------------
        # Plot : Averaging Kernel Smoothing Function (row of avk)
        #---------------------------------
        indAlt = np.where(alt[gasVer] < 20.)[0]


        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        fig       = plt.figure(figsize=(9,7))
        gs        = gridspec.GridSpec(1,3,width_ratios=[3,1,1])
        ax        = plt.subplot(gs[0])
        axb       = plt.subplot(gs[1])
        axc       = plt.subplot(gs[2])
        cm        = plt.get_cmap(clmap)
        cNorm     = colors.Normalize(vmin=np.min(alt[gasVer][indAlt]), vmax=np.max(alt[gasVer][indAlt]))
        scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        scalarMap.set_array(alt[gasVer][indAlt])

        #---------------------------------
        ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt[gasVer][indAlt]])
        
        for i in range(len(alt[gasVer][indAlt])):
            ax.plot(avkSCFav[gasVer][indAlt[i],indAlt],alt[gasVer][indAlt])
            
        ax.set_ylabel('Altitude [km]', fontsize=16)
        ax.set_xlabel('AK', fontsize=16)
        #ax.grid(True, alpha=0.5)
        #ax.set_title('(a)', loc='left', fontsize=14)
        ax.text(0.025, 0.95,'(a)', fontsize=16,transform=ax.transAxes)   

        cbaxes = fig.add_axes([0.4, 0.55, 0.02, 0.4]) 
        cbar = fig.colorbar(scalarMap, orientation='vertical', cax = cbaxes)
        cbar.set_label('Altitude [km]', fontsize=16)
        #ax.set_title('H2O Averaging Kernels Scale Factor', fontsize=14)
        ax.tick_params(labelsize=16)
        ax.axvline(x=0,color='k', linestyle='--')
        #ax.set_ylim(1, 15)
        if loc.lower() == 'fl0': ax.set_ylim(1, 15)
        if loc.lower() == 'mlo': ax.set_ylim(3, 15)
        if loc.lower() == 'tab': ax.set_ylim(0, 15)

        #----------------------------------------
        # Calculate total column averaging kernel
        #----------------------------------------
        nobs            = len(dates[gasVer])
        nlyrs           = len(alt[gasVer])
        avkTC           = np.zeros((nobs,nlyrs))

        for i in range(0,nobs):
            AirMtemp  = np.squeeze(Airmass[gasVer][i,:])
            akTemp    = np.squeeze(avkSCF[gasVer][i,:,:])
            AirMinv   = np.diag(1.0/AirMtemp)
            avkTC[i,:] = np.dot(np.dot(AirMtemp,akTemp),AirMinv)

        avkTCAv        = np.mean(avkTC, axis=0)

        #---------------------------------        
        #axb.plot(np.sum(avkSCFav[gasVer],axis=0), alt[gasVer],color='k')
        axb.plot(avkTCAv, alt[gasVer],color='k')

        #axb.grid(True,  alpha=0.5)
        axb.set_xlabel('Total Column AK', fontsize=16)
        axb.tick_params(labelsize=16)
        #axb.set_title('(b)', loc='left', fontsize=14)
        axb.text(0.725, 0.95,'(b)', fontsize=16,transform=axb.transAxes)  

        major_ticks = np.arange(0, 3, 1)
        axb.set_xticks(major_ticks) 
        #axb.set_ylim(1, 15)

        if loc.lower() == 'fl0': axb.set_ylim(1, 15)
        if loc.lower() == 'mlo': axb.set_ylim(3, 15)
        if loc.lower() == 'tab': axb.set_ylim(0, 15) 

        #---------------------------------
        dofs_cs = np.cumsum(np.diag(avkSCFav[gasVer])[::-1])[::-1]
        axc.plot(dofs_cs,alt[gasVer],color='k',label='Cumulative Sum of DOFS (starting at surface)')
        #axc.plot(np.sum(avkSCFav[gasVer],axis=0), alt[gasVer],color='k')
        xval = range(0,int(np.ceil(max(dofs_cs)))+2)

        #ind1         = mf.nearestind(Pcol[0], alt[gasVer])
        #ind2         = mf.nearestind(Pcol[1], alt[gasVer])

        #axc.fill_between(xval,alt[gasVer][ind1],alt[gasVer][ind2],alpha=0.5,color='0.75')  
        #axc.axhline(alt[gasVer][ind2],color='k',linestyle='--')
        #dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
        #axc.text(0.15,(alt[idhdf][ind1]+alt[idhdf][ind2])/2.0, 
        #         'DOFs for layer {0:.2f}-{1:.2f}[km] = {2:.3f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol),
        #         fontsize=9)
        #axc.set_title('DOFs Profile - ' +str(pltID[i]))
        #axc.set_ylabel('Altitude [km]')
        
        axc.set_xlabel('Cumulative\nSum of DOFS', fontsize=16)  
        axc.tick_params(labelsize=16)
        #axc.set_title('(c)', loc='left', fontsize=14)
        axc.text(0.025, 0.95,'(c)', fontsize=16,transform=axc.transAxes)  
        #axc.set_title('DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol), fontsize=9)    
        
        axc.set_xticks(major_ticks)

        axc.set_ylim((0,120))
        #axc.grid(True,which='both',  alpha=0.5)
        if loc.lower() == 'fl0': axc.set_ylim(1, 15)
        if loc.lower() == 'mlo': axc.set_ylim(3, 15)
        if loc.lower() == 'tab': axc.set_ylim(0, 15)

        plt.tight_layout()

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(FigDir+'avkSCFavPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
        else:       plt.show(block=False)

        #---------------------------------
        # Plot : VMR Averaging Kernel Smoothing Function (row of avk)
        #---------------------------------
        # clmap        = 'jet'
        # cm           = plt.get_cmap(clmap)
        # fig       = plt.figure(figsize=(7,7))
        # gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
        # ax        = plt.subplot(gs[0])
        # axb       = plt.subplot(gs[1])
        # cm        = plt.get_cmap(clmap)
        # cNorm     = colors.Normalize(vmin=np.min(alt[gasVer]), vmax=np.max(alt[gasVer]))
        # scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
        # scalarMap.set_array(alt[gasVer])
        # ax.set_color_cycle([scalarMap.to_rgba(x) for x in alt[gasVer]])
        
        # for i in range(len(alt[gasVer])):
        #     ax.plot(avkVMRav[gasVer][i,:],alt[gasVer])
            
        # ax.set_ylabel('Altitude [km]', fontsize=14)
        # ax.set_xlabel('Averaging Kernels', fontsize=14)
        # ax.grid(True)
        # ax.set_title('H2O Averaging Kernels [VMR]', fontsize=14)
        # ax.tick_params(labelsize=14)
        
        # axb.plot(np.sum(avkVMRav[gasVer],axis=0), alt[gasVer],color='k')
        # axb.grid(True)
        # axb.set_xlabel('Averaging Kernel Area', fontsize=14)
        # axb.tick_params(axis='x',which='both',labelsize=8)
        # #axb.tick_params(labelsize=14) 

        # cbar = fig.colorbar(scalarMap, orientation='vertical')
        # cbar.set_label('Altitude [km]', fontsize=14)

        # if saveFlg: 
        #     pdfsav.savefig(fig,dpi=200)
        #     plt.savefig(FigDir+'avkVMRavPrf_'+loc.upper()+'.pdf', bbox_inches='tight')
        # else:       plt.show(block=False)

        #---------------------------------
        # Plot :Partial Columns
        #---------------------------------

        #fig = plt.figure(figsize=(16,12))
        fig = plt.figure(figsize=(12,8))

        outer_grid = gridspec.GridSpec(2, 2, wspace=0.17, hspace=0.25)

        for p, pcol in enumerate(pCols):

            ax = plt.Subplot(fig, outer_grid[p])

            #----------------------------------------------------------------------
            #analysis of time series of DAILY averages
            #----------------------------------------------------------------------
            #dailyVals        = mf.dailyAvg(vmrP[gasVer+str(p)][0], dates[gasVer],dateAxis=1, meanAxis=0)
            mnthlyVals       = mf.mnthlyAvg(vmrP[gasVer+str(p)][0],dates[gasVer],dateAxis=1, meanAxis=0)
            #dateYearFrac     = mf.toYearFraction(dailyVals['dates'])
            dateYearFrac     = mf.toYearFraction(mnthlyVals['dates'])
            weights          = np.ones_like(dateYearFrac)

            #res              = mf.fit_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            res          = mf.fit_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            intercept, slope, pfourier = res[0:3]
            f_drift, f_fourier, f_driftfourier = res[3:6]

            #DT_d             = dailyVals['dates']
            DT_d             = mnthlyVals['dates']
            #vmrP_d           = dailyVals['dailyAvg']
            vmrP_d           = mnthlyVals['mnthlyAvg']
            
            #TCsd_d           = dailyVals['std']
            TCsd_d           = mnthlyVals['std']

            int_fourier_d    = intercept
            sl_fourier_d     = slope
            p_fourier_d      = pfourier
            FAT_d            = f_drift(dateYearFrac)
            FATAV_d          = f_driftfourier(dateYearFrac)

            #perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, dailyVals['dailyAvg'], weights, 2)
            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFrac, mnthlyVals['mnthlyAvg'], weights, 2)
            int_boot_d       = intercept_boot
            sl_boot_d        = slope_boot
            p_boot_d         = pfourier_boot


            #Rate_d           = np.divide(res[1], np.mean(dailyVals['dailyAvg']))  *100.0 
            #Rate_e_d         = np.divide(np.std(slope_boot), np.mean(dailyVals['dailyAvg']))*100.0

            Rate_d           = np.divide(res[1], np.mean(mnthlyVals['mnthlyAvg']))  *100.0 
            Rate_e_d         = np.divide(np.std(slope_boot), np.mean(mnthlyVals['mnthlyAvg']))*100.0

            #----------------------------------------------------------------------
            #analysis of time series of Sonde
            #----------------------------------------------------------------------

            datesCFH          = np.asarray(sondedates[p])
            vmrPCFH           = np.asarray(sondevmrP[p])

            #if p == 5:
            #    indmax = np.where(vmrPCFH > 150.)[0]
            #    datesCFH       =  np.delete(datesCFH, indmax, axis=0)
            #    vmrPCFH        =  np.delete(vmrPCFH, indmax, axis=0)

            mnthlyValsCFH       = mf.mnthlyAvg(vmrPCFH,datesCFH, dateAxis=1, meanAxis=0)
           
            dateYearFracCFH   = mf.toYearFraction(mnthlyValsCFH['dates'])
            weights           = np.ones_like(dateYearFracCFH)
            resCFH            = mf.fit_driftfourier(dateYearFracCFH, mnthlyValsCFH['mnthlyAvg'], weights, 2)

            interceptCFH, slopeCFH, pfourierCFH          = resCFH[0:3]
            f_driftCFH, f_fourierCFH, f_driftfourierCFH  = resCFH[3:6]

            FAT_CFH            = f_drift(dateYearFracCFH)
            FATAV_CFH          = f_driftfourier(dateYearFracCFH)

            perc, intercept_boot, slope_boot, pfourier_boot = mf.cf_driftfourier(dateYearFracCFH, mnthlyValsCFH['mnthlyAvg'], weights, 2)
            int_boot_CFH       = intercept_boot
            sl_boot_CFH        = slope_boot
            p_boot_CFH         = pfourier_boot

            Rate_CFH           = np.divide(resCFH[1], np.mean(mnthlyValsCFH['mnthlyAvg']))  *100.0 
            Rate_e_CFH         = np.divide(np.std(slope_boot), np.mean(mnthlyValsCFH['mnthlyAvg']))*100.0

           #  #----------------------------------------------------------------------
        
            ax.scatter(DT_d,vmrP_d,facecolors='red', edgecolors='black', s=35, label='Daily averages')
            ax.plot(DT_d,FAT_d, color='blue', linewidth=2.5)
            ax.plot(DT_d,FATAV_d, color='green', linewidth=2.5)

            ax.scatter(mnthlyValsCFH['dates'],mnthlyValsCFH['mnthlyAvg'],facecolors='blue', edgecolors='black', s=35, label='CFH')

            ax.xaxis.set_major_locator(yearsLc)
            ax.xaxis.set_minor_locator(monthsFmt)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.xaxis.set_tick_params(which='major',labelsize=18)
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')
           
            ax.grid(True)
            if (p == 2) or (p == 3): ax.set_xlabel('Year', fontsize=14)
            if (p == 0) or (p == 2): ax.set_ylabel('VMR [ppm]', fontsize=14)

            ax.set_ylim(ymin=0)
            ax.tick_params(axis='x',which='both',labelsize=14)
            ax.tick_params(axis='y',which='both',labelsize=14)

            #ax.text(0.02,0.94,"Layer: {0:.1f} - {1:.1f} km".format(pcol[0], pcol[1]),transform=ax.transAxes)
            ax.text(0.02, 0.94,"HR-FTIR fitted trend: {0:.2f}$\pm${1:.2f} %/year".format(Rate_d, Rate_e_d), fontsize=12,transform=ax.transAxes)            
            #ax2.text(0.02,0.88,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(mnthlyVals['mnthlyAvg'])*100.0),transform=ax2.transAxes)

            ax.text(0.02, 0.88,"FPH fitted trend: {0:.2f}$\pm${1:.2f} %/year".format(Rate_CFH, Rate_e_CFH), fontsize=12,transform=ax.transAxes)
            
            #ax.text(0.02,0.82,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(resCFH[1],resCFH[1]/np.mean(mnthlyValsCFH['mnthlyAvg'])*100.0),transform=ax.transAxes)
           
           # ax2.text(0.02,0.89,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax2.transAxes)
            #ax2.text(0.02,0.84,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(dailyVals['dailyAvg'])*100.0),transform=ax2.transAxes) 

            ax.set_title("Layer: {0:.1f} - {1:.1f} km".format(pcol[0], pcol[1]))

            fig.add_subplot(ax)
            fig.subplots_adjust(left=0.1, bottom=0.075, right=0.975, top=0.95)
            

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            plt.savefig(FigDir+'TSPC_'+loc.upper()+'.pdf', bbox_inches='tight')

        else:       plt.show(block=False)


    if saveFlg: pdfsav.close()
    else:
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()


 
if __name__ == "__main__":
    main()