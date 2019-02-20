#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         printCAM.py
#
# Purpose:
#         Read and Create a .dat file for input gas in CAM-Chem using Smoothing AK (see ianputs below)
#
# Notes:
#   
#
# Version History:
#       Created, July, 2018  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#

from itertools import izip
import ClassFTS as cfts
import numpy as np

from scipy import interpolate
import datetime as dt

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
from cycler import cycler


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
    loc                = 'fl0'
    gasName            = ['co',           'c2h2',         'c2h6',         'ch4',         'nh3',          'hcn',         'h2co',         'hcooh']             
    ver                = ['Current_v3',   'Current_v2',   'Current_v2',   'Current_WP',  'Current_v2',   'Current_v1', 'Current_v8',   'Current_v3']          # Name of retrieval version to process
    ctlF               = ['sfit4_v3.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_3.ctl', 'sfit4_v2.ctl', 'sfit4.ctl',  'sfit4_v8.ctl', 'sfit4_v3.ctl'] 

    #gasName            = ['h2co',       ]
    #ver                = ['Current_v8'  ]        # Name of retrieval version to process
    #ctlF               = ['sfit4_v8.ctl'] 


    #----------------------
    # First Level Retrieval Directory
    #----------------------
    retDir             = '/data1/ebaumer/'+loc.lower()

    dataPath           = '/data/iortega/results/'+loc.lower()+'/data/'

    saveFlg           = True
    pltDir            =  '/data/iortega/pbin/tropGases/fig/'
    pltFile           =  pltDir + 'printCAM_v2.pdf'
 
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
    tcMMFlg            = False                  # Flag to filter based on min and max total column amount
    cnvrgFlg           = True                   # Flag to filter profiles that did not converge
    rmsFlg             = True                   # Flag to filter based on max RMS
    chiFlg             = False                  # Flag to filter based on max CHI_2_Y
    mnthFlg            = False                  # Flag to filter based on 
    mnths              = [6, 7, 8]

    maxRMS             = [1.5, 0.6, 1.0, 0.3, 1.0, 1.5, 0.5, 1.5]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0, 0.5, 0.5, 0.9, 0.5, 0.5,  0.5, 0.5]                      # Min DOFs for filtering

    #maxRMS             = [0.5]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    #minDOF             = [0.5]  

    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    pColsFlg           = True
    pCols              = [1.6, 8.0]              #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR

    #-----------------------------------------------------------------------------------------
    #                 Initializations for CAM-CHEM
    #-----------------------------------------------------------------------------------------
    #dataDirCAM        = '/data1/ancillary_data/'+loc.lower()+'/'
    #fileCAM           = 'CAM_chem_fmerra_FSDSSOA_2deg_2000_2014_extra_Boulder.nc'       

    #dataDirCAM        = '/net/modeling1/data16a/buchholz/CAM_chem_output/CAM_chem_fmerra2_FCSD_1deg_Boulder/'  # 
    dataDirCAM        = '/net/modeling1/data16a/buchholz/CAM_chem_output/CAM_chem_fmerra2_FCSD_1deg_Boulder_finn/'  # 

    #fileCAM           = 'CAM_chem_fmerra2_FCSD_1deg_Boulder_2009_2017.nc'
    fileCAM           = 'CAM_chem_fmerra2_FCSD_1deg_FINN_Boulder_2009_2017.nc'

    sLat              = 40.4             #--LATITUDE OF BOULDER
    sLon              = -105.24          #--LONGITUDE OF BOULDER

    interpFlg         = False

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2017
    fmnth              = 12
    fday               = 31


                                    #----------------------------#
                                    #                            #
                                    #        --- START ---       #
                                    #                            #
                                    #----------------------------#

    if saveFlg: pdfsav = PdfPages(pltFile)

    #-------------------------------------------------
    #    -- Read CAM-CHEM --
    #-------------------------------------------------
    DataCAM = dm.CAMClass(dataDirCAM, fileCAM,  outFname= '', saveFlg= False)
    DataCAM.ReadOutputCAM(gasName, pCols, sLat, sLon, interpFlg=interpFlg)
    #DataCAM.PltCAM()

    # DatesCAM      = np.asarray(DataCAM.CAM['dates'])
    # PrfCAM        = np.asarray(DataCAM.CAM['GasPrf_'+gasName[0]][:,:,0,0])*1e9
    
    # mnthlyVals            = mf.mnthlyAvg(PrfCAM, DatesCAM, dateAxis=0, meanAxis=0)
    # vmrPCAM_smt_mnth      = mnthlyVals['mnthlyAvg']
    # vmrPCAM_smt_STDmnth   = mnthlyVals['std']
    # DatesCAM_mnth         = mnthlyVals['dates']

    # print mnthlyVals['dates']

    #exit()



    #user_input = raw_input('Press any key to exit >>> ')
    #sys.exit()           

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------
   
    fts = cfts.FTSClass( gasName, retDir,  ctlF, ver, iyear,imnth, iday, fyear, fmnth, fday)

    fts.ReadFTS(fltrFlg=fltrFlg, sclfct=sclfct, sclname=sclfctName, mnthFltr= mnths, mnthFltFlg=mnthFlg,
          errFlg=errorFlg, minSZA=minSZA, maxSZA=maxSZA, maxRMS=maxRMS, minTC=minTC, maxTC= maxTC,
          minDOF=minDOF, maxCHI=maxCHI, dofFlg=dofFlg, rmsFlg=rmsFlg, tcFlg=tcNegFlg,
          pcFlg=pcNegFlg, szaFlg=szaFlg,cnvrgFlg=cnvrgFlg, chiFlg=chiFlg,tcMMflg=tcMMFlg, pColsFlg=pColsFlg, pCols=pCols)

     
    #---------------------------------
    # 
    #---------------------------------
    GasVer = []

    for (g,v) in izip(gasName, ver):
        GasVer.append(g+'_'+v)


    #--------------

    for (g,v) in izip(gasName, ver):

        if DataCAM.CAM['GasFlg_'+g.lower()]:

            gv = g+'_'+v

            #---------------------------------
            # Defining variables (CAM-CHEM) 
            #---------------------------------
            altCAM        = np.asarray(DataCAM.CAM['midpoints'][0, :, DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])
            DatesCAM      = np.asarray(DataCAM.CAM['dates'])
            PrfCAM        = np.asarray(DataCAM.CAM['GasPrf_'+g][:,:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])*1e9
            AirmassCAM    = np.asarray(DataCAM.CAM['AIRMASS_'+g][:,:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])
            pressCAM      = np.asarray(DataCAM.CAM['pressure'][:,:,DataCAM.CAM['indsLoc'][0],DataCAM.CAM['indsLoc'][1]])

        
            #---------------------------------
            # Defining variables (FTIR)
            #---------------------------------
            altFTS       = fts.alt[gv]
            DatesFTS     = fts.dates[gv]
            PrfFTS       = fts.rPrfVMR[gv]
            AirmassFTS   = fts.Airmass[gv]
            aPrfFTS      = fts.aPrfVMR[gv]
            akVMRFTS     = fts.avkVMR[gv]
            vmrPFTS      = fts.vmrP[gv]

            aPrfFTSMean  = np.mean(aPrfFTS, axis=0)
            akVMRFTSMean = np.mean(akVMRFTS, axis=0)

            #----------------------------
            #  Monthly - FTIR
            #----------------------------
            mnthlyVals       = mf.mnthlyAvg(PrfFTS, DatesFTS, dateAxis=0, meanAxis=0)
            PrfFTS_mnth      = mnthlyVals['mnthlyAvg']
            DatesFTS_mnth    = mnthlyVals['dates']

            mnthlyVals       = mf.mnthlyAvg(AirmassFTS, DatesFTS, dateAxis=0, meanAxis=0)
            AirmassFTS_mnth  = mnthlyVals['mnthlyAvg']

            mnthlyVals       = mf.mnthlyAvg(PrfFTS, DatesFTS, dateAxis=0, meanAxis=0)
            PrfFTS_mnth      = mnthlyVals['mnthlyAvg']

            mnthlyVals       = mf.mnthlyAvg(akVMRFTS, DatesFTS, dateAxis=0, meanAxis=0)
            akVMRFTS_mnth    = mnthlyVals['mnthlyAvg']

            mnthlyVals       = mf.mnthlyAvg(vmrPFTS, DatesFTS, dateAxis=0, meanAxis=0)
            vmrPFTS_mnth     = mnthlyVals['mnthlyAvg']
            vmrPFTS_STDmnth  = mnthlyVals['std']


            #----------------------------
            # Altitude Interpolation
            #----------------------------
            PrfCAM_int     = interpolate.interp1d(altCAM, PrfCAM, axis=1, fill_value='extrapolate', bounds_error=False)(altFTS)
            AirmassCAM_int = interpolate.interp1d(altCAM, AirmassCAM, axis=1, fill_value='extrapolate', bounds_error=False)(altFTS)
            pressCAM_int   = interpolate.interp1d(altCAM, pressCAM, axis=1, fill_value='extrapolate', bounds_error=False)(altFTS)

            #----------------------------
            # Smoothing 
            #----------------------------
            PrfCAM_smt = np.zeros( (len(DatesCAM), len(altFTS)) )
            
            for itime in range(len(DatesCAM)):            
                PrfCAM_smt[itime, :]  = aPrfFTSMean + np.dot(akVMRFTSMean, (PrfCAM_int[itime, :] -  aPrfFTSMean))

            AirmassFTSMean = np.mean(AirmassFTS, axis=0)
            AirmassFTSstd = np.std(AirmassFTS, axis=0)

            #----------------------------
            # Weighted VMR 
            #----------------------------
            inds         = np.where( (altFTS >= pCols[0]) & (altFTS <=pCols[1])  )[0]
            vmrPCAM_smt  = np.average(PrfCAM_smt[:,inds],axis=1,weights=pressCAM_int[:,inds])
            vmrPCAM_int  = np.average(PrfCAM_int[:,inds],axis=1,weights=pressCAM_int[:, inds]) 
            vmrPCAM      = np.average(PrfCAM[:,inds],axis=1,weights=pressCAM[:,inds])                 
               
            #---------------------------------
            # Plot : Prf Mean
            #---------------------------------

            PrfmeanCAM     = np.mean(PrfCAM, axis=0)
            prfSTDCAM      = np.std(PrfCAM, axis=0)

            PrfmeanCAM_int = np.mean(PrfCAM_int, axis=0)
            prfSTDCAM_int  = np.std(PrfCAM_int, axis=0)

            PrfmeanCAM_smt = np.mean(PrfCAM_smt, axis=0)
            prfSTDCAM_smt  = np.std(PrfCAM_smt, axis=0)

            prfMeanFTS     = np.mean(PrfFTS_mnth, axis=0)
            prfSTDFTS      = np.std(PrfFTS_mnth, axis=0)

            #----------------------------
            #  Monthly - CAM-Chem
            #---------------------------

            mnthlyVals            = mf.mnthlyAvg(vmrPCAM_smt, DatesCAM, dateAxis=0, meanAxis=0)
            vmrPCAM_smt_mnth      = mnthlyVals['mnthlyAvg']
            vmrPCAM_smt_STDmnth   = mnthlyVals['std']
            DatesCAM_mnth         = mnthlyVals['dates']

            mnthlyVals            = mf.mnthlyAvg(vmrPCAM, DatesCAM, dateAxis=0, meanAxis=0)
            vmrPCAM_mnth          = mnthlyVals['mnthlyAvg']
            vmrPCAM_STDmnth       = mnthlyVals['std']
        

            #----------------------------
            # Print in Dat 
            #---------------------------- 
            with open(dataPath + g+'_CAM-CHEM_v2_finn.dat','w') as fopen:

                YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(d.year, d.month, d.day)    for d in DatesCAM_mnth])
         
                fopen.write('Index, YYYY-MM-DD, wVMRsmth, wVMRsmth_std, wVMR, wVMR_std \n')
                strFormat = '{0:d}, {1:>10s}, {2:.3E}, {3:.3E}, {4:.3E}, {5:.3E}\n'


                for i, sngTime in enumerate(YYYYMMDD):
                  
                    fopen.write(strFormat.format((i+1), YYYYMMDD[i], vmrPCAM_smt_mnth[i], vmrPCAM_smt_STDmnth[i], vmrPCAM_mnth[i], vmrPCAM_STDmnth[i] ) )

            #----------------------------
            #
            #----------------------------
            xmin   = dt.date(iyear, imnth, iday)
            xmax   = dt.date(fyear, fmnth, fday)

            clmap        = 'jet'
            cm           = plt.get_cmap(clmap)
            dayLc        = DayLocator()
            yearsLc      = YearLocator()
            monthLc      = MonthLocator()
            mondays      = WeekdayLocator(MONDAY)
            #DateFmt      = DateFormatter('%b %d')
            DateFmt      = DateFormatter('%Y')

            #---------------------------------
            # Plot : Averaging Kernel Smoothing Function (row of avk)
            #---------------------------------
            fig       = plt.figure(figsize=(9,9))
            gs        = gridspec.GridSpec(1,2,width_ratios=[3,1])
            ax        = plt.subplot(gs[0])
            axb       = plt.subplot(gs[1])
            cm        = plt.get_cmap(clmap)
            cNorm     = colors.Normalize(vmin=np.min(altFTS), vmax=np.max(altFTS))
            scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=clmap)
            scalarMap.set_array(altFTS)
            
            #ax.set_color_cycle([scalarMap.to_rgba(x) for x in altFTS])

            ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in altFTS] ) )
            
            for i in range(len(altFTS)):
                ax.plot(akVMRFTSMean[i,:],altFTS)
                
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('Averaging Kernels', fontsize=14)
            ax.grid(True)
            cbar = fig.colorbar(scalarMap, orientation='vertical')
            cbar.set_label('Altitude [km]', fontsize=14)
            ax.set_title(g.upper() + ' Averaging Kernels Scale Factor', fontsize=14)
            ax.tick_params(labelsize=14)
            
            axb.plot(np.sum(akVMRFTSMean,axis=0), altFTS,color='k')
            axb.grid(True)
            axb.set_xlabel('Averaging Kernel Area', fontsize=14)
            axb.tick_params(axis='x',which='both',labelsize=8)
            #axb.tick_params(labelsize=14) 

            if saveFlg:     pdfsav.savefig(fig,dpi=200)
            else:           plt.show(block=False)

            #---------------------------------
            # Plot : Prf
            #---------------------------------
            fig,ax  = plt.subplots(figsize=(7,9))

            ax.plot(prfMeanFTS,altFTS, linewidth = 2.0, color='k', label='FTIR')
            ax.scatter(prfMeanFTS,altFTS, facecolors='white', s=60, color='k')
            ax.fill_betweenx(altFTS,prfMeanFTS-prfSTDFTS,prfMeanFTS+prfSTDFTS, alpha=0.5, color='k') 

            ax.plot(PrfmeanCAM,altCAM, linewidth = 2.0, color='red', label='CAM-Chem')
            ax.scatter(PrfmeanCAM,altCAM, facecolors='white', s=60, color='red')
            ax.fill_betweenx(altCAM,PrfmeanCAM-prfSTDCAM,PrfmeanCAM+prfSTDCAM,alpha=0.5,color='red') 

            ax.plot(PrfmeanCAM_int,altFTS, linewidth = 2.0, color='green', label='CAM-Chem - FTS grid')
            ax.scatter(PrfmeanCAM_int,altFTS, facecolors='white', s=60, color='green')
            ax.fill_betweenx(altFTS,PrfmeanCAM_int-prfSTDCAM_int,PrfmeanCAM_int+prfSTDCAM_int,alpha=0.5,color='green') 

            ax.plot(PrfmeanCAM_smt,altFTS, linewidth = 2.0, color='blue', label='CAM-Chem - Smoothed')
            #ax.scatter(PrfmeanCAM_smt,altFTS, facecolors='white', s=60, color='blue')
            #ax.fill_betweenx(altFTS,PrfmeanCAM_smt-prfSTDCAM_smt,PrfmeanCAM_smt+prfSTDCAM_smt,alpha=0.5,color='blue') 

            ax.set_title('Mean Profile of '+g.upper(), fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            ax.grid(True,which='both')
            ax.tick_params(labelsize=14)
            ax.legend(prop={'size':12})
            ax.set_ylim(0, 20)
            ax.set_xlim(0, 200)

            if saveFlg: 
                pdfsav.savefig(fig,dpi=200)
                #plt.savefig(pltDir+'Time_Series_all.pdf', bbox_inches='tight')
            else:
                plt.show(block=False) 

            #---------------------------------
            # Plot : Monthly averages of total columns
            #---------------------------------
            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)
            ax.plot(DatesFTS_mnth, vmrPFTS_mnth, color='k', label='FTIR')
            ax.scatter(DatesFTS_mnth, vmrPFTS_mnth, facecolors='white', s=60, color='k')
            ax.fill_between(DatesFTS_mnth, vmrPFTS_mnth - vmrPFTS_STDmnth, vmrPFTS_mnth + vmrPFTS_STDmnth, alpha=0.25, zorder=1)

            ax.plot(DatesCAM_mnth, vmrPCAM_mnth, color='red',  linewidth=3.0, label='CAM-Chem')
            #ax.scatter(DatesCAM, vmrPCAM, facecolors='white', s=60, color='red')

            ax.plot(DatesCAM_mnth, vmrPCAM_smt_mnth, color='blue',  linewidth=3.0, label='CAM-Chem smoothed')
            #ax.scatter(DatesCAM, vmrPCAM_smt, facecolors='white', s=60, color='blue')
            
            ax.grid(True)
            #ax.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax.set_ylabel('Weighted VMR [ppb]',fontsize=14)
            ax.set_title(g.upper() + ' Partial Column [Monthly average], '+ str(altFTS[inds[-1]])+'[km] - '+str(altFTS[inds[0]])+'[km]',multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})
            #ax.set_xlim(xmin, xmax)

            if saveFlg: 
                pdfsav.savefig(fig,dpi=200)
                #plt.savefig(pltDir+'Time_Series_all.pdf', bbox_inches='tight')
            else:
                plt.show(block=False) 

    
    if saveFlg: pdfsav.close()
    else:
        user_input = raw_input('Press any key to exit >>> ')
        exit()       

      
 
if __name__ == "__main__":
    main()