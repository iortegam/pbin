#!/usr/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        print2file.py
#
# Purpose:
#        The purpose of this program is to write total column and associated data
#        to icartt format
#----------------------------------------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import os
import dataOutClass        as dc
import numpy               as np
import datetime            as dt

import matplotlib.pyplot as plt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...





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

    #----------------
    # Initializations
    #----------------
    #loc        = 'fl0'           
    #gasName    = 'ch4'            #'nh3'           #'hcooh'         #'hcn'           #'h2co'            #'ch4'           #'c2h2'         #'c2h6'           #'co'                 
    #ver        = 'Current_WP'    #'Current_v2'    #'Current_v2'    #'Current_WP'    #'Current_WP_v5'   #'Current_WP'    #'Current_v2'   #'Current_v2'     #'Current_v3'             
    #ctlF       = 'sfit4_3.ctl'     #'sfit4_v2.ctl'  #'sfit4_v2.ctl'  #'sfit4.ctl'     #'sfit4_v5.ctl'    #'sfit4_3.ctl'   #'sfit4_v2.ctl' # 'sfit4_v2.ctl'  #'sfit4_v3.ctl'
    
    loc        = 'mlo'
    gasName    = 'hcooh'
    ver        = 'Current_v3'
    ctlF       = 'sfit4_v3.ctl'

    #------
    # Flags
    #------
    fltrFlg    = True                  # Flag to filter the data   
    errorFlg   = True
    saveFlg    = True                  # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
    
    dofFlg     = True                   # Flag to filter based on min DOFs
    pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
    tcNegFlg   = True                   # Flag to filter profiles with negative total columns
    cnvrgFlg   = True                   # Flag to filter profiles that did not converge
    rmsFlg     = True                   # Flag to filter based on max RMS
    
    maxRMS     = 1.0                         # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = 0.1                    # Min DOFs for filtering
    sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName = 'ppbv'                 # Name of scale factor for labeling plots    
  
    #----------------------
    # Date range to process
    #----------------------
    iyear      = 2017
    imnth      = 1
    iday       = 1
    fyear      = 2017
    fmnth      = 12
    fday       = 31


    #------------
    # Directories
    #------------
    retDir  = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'    # Retrieval data directory
    
    #------
    # Files
    #------
    ctlFile     = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF
    outFileDir  = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/' 
    #outFileDir  = '/data/Archive/'+loc.upper()+'/'+gasName.upper()+'/' 
    icarttHdr   = '/data1/ebaumer/icartt_Hdr.txt'

    pltFile     =  outFileDir + loc.upper()+ '-FTS-'+gasName.upper()+'.pdf'

    if saveFlg: pdfsav = PdfPages(pltFile)

    #---------------------------
    # Check file and directories
    #---------------------------
    ckDir(retDir,exit=True)
    ckFile(ctlFile,exit=True)
    ckDir(outFileDir,exit=True)

    #-------------------------------------
    # Create instance of output data class   
    #-------------------------------------
    statDataCl = dc.ReadOutputData(retDir,'',ctlFile,iyear,imnth,iday,fyear,fmnth,fday)
    
    #----------------------------------
    # Read Summary data (For filtering)
    #----------------------------------
    statDataCl.readsummary()
    statDataCl.readPbp()
    totClmn  = np.asarray(statDataCl.summary[statDataCl.PrimaryGas.upper()+'_RetColmn'])   
    dates    = np.asarray(statDataCl.summary['date'])
    dofs     = np.asarray(statDataCl.summary[statDataCl.PrimaryGas.upper()+'_DOFS_TRG']) 
    sza      = np.asarray(statDataCl.pbp['sza'])
    rms      = np.asarray(statDataCl.summary[statDataCl.PrimaryGas.upper()+'_FITRMS'])

    #-------------------------------------
    # Naming output file  
    #-------------------------------------
    idate    = dates[0]
    fdate    = dates[-1]

    idateStr = "{0:04d}{1:02d}{2:02d}T{3:02d}{4:02d}{5:02d}Z".format(idate.year,idate.month,idate.day,idate.hour,idate.minute,idate.second)
    fdateStr = "{0:04d}{1:02d}{2:02d}T{3:02d}{4:02d}{5:02d}Z".format(fdate.year,fdate.month,fdate.day,fdate.hour,fdate.minute,fdate.second)

    outFileName = outFileDir + loc.upper()+ '-FTS-'+gasName.upper()+'_'+idateStr+'_'+fdateStr+'.ascii'

    #----------------
    # Read Error Data
    #----------------
    if errorFlg:
        statDataCl.readError(totFlg=True,sysFlg=False,randFlg=False,vmrFlg=False,avkFlg=False,KbFlg=False)
        tot_rnd  = np.array(statDataCl.error['Total random uncertainty'])
        tot_sys  = np.array(statDataCl.error['Total systematic uncertainty'])
        tot_std  = np.sqrt(tot_rnd**2 + tot_sys**2) 


    
    #----------------------------------
    # Read Profile Data (for filtering)
    #----------------------------------
    statDataCl.readprfs([statDataCl.PrimaryGas],retapFlg=1)         # Retrieved Profiles
    
    #--------------------
    # Call to filter data
    #--------------------
    if fltrFlg: statDataCl.fltrData(statDataCl.PrimaryGas,mxrms=maxRMS,rmsFlg=rmsFlg,tcFlg=tcNegFlg,pcFlg=pcNegFlg,cnvrgFlg=cnvrgFlg, minDOF=minDOF,  dofFlg=dofFlg)
    else:       statDataCl.inds = np.array([])    

    
    if statDataCl.empty: 
        print 'No retreivals found. Exiting.......'
        sys.exit()    
    
    #---------------------------------
    # Remove data based on filter inds
    #---------------------------------
    totClmn = np.delete(totClmn,statDataCl.inds)
    dates   = np.delete(dates,statDataCl.inds)
    dofs    = np.delete(dofs,statDataCl.inds)
    sza     = np.delete(sza,statDataCl.inds)
    rms     = np.delete(rms,statDataCl.inds)
    
    if errorFlg:
        tot_rnd = np.delete(tot_rnd,statDataCl.inds)
        tot_sys = np.delete(tot_sys,statDataCl.inds)
        tot_std = np.delete(tot_std,statDataCl.inds)

    YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(dates[i].year, dates[i].month, dates[i].day)    for i,dum in enumerate(dates)])
    HHMMSS   = np.asarray(['{0:02d}:{1:02d}:{2:02d}'.format(dates[i].hour,dates[i].minute,dates[i].second) for i,dum in enumerate(dates)])
    
    #-----------------------------------------
    # Construct icartt file name based on first 
    # day of data
    #------------------------------------------
    #idateStr    = '{0:04d}{1:02d}{2:02d}'.format(dates[0].year,dates[0].month,dates[0].day)
    
    
    #----------------------------------------
    # Convert times to second since UTC start 
    # of day of measurements (dates[0])
    #----------------------------------------
    times = np.array([(x - dt.datetime(dates[0].year,dates[0].month,dates[0].day,0,0,0)).total_seconds() for x in dates])

    #-----------------------------------
    # Find current date of file creation
    #-----------------------------------
    dateNow = dt.datetime.now()

    if loc.lower() == 'tab':
        locStr = 'Thule, Greenland'
    elif loc.lower() == 'mlo':
        locStr = 'Mauna Loa, HI USA'
    elif loc.lower() == 'fl0':
        locStr = 'Boulder, Colorado USA'

    #--------------
    with open(outFileName,'w') as fopen:
        fopen.write('#Hannigan, J.W.; Ortega, I\n')
        fopen.write('#National Center for Atmospheric Research\n')
        fopen.write('#Ground Based HR-FTIR Spectrometer at {}\n'.format(locStr))
        fopen.write('#CONTACT_INFO: Hannigan, Jim, jamesw@ucar.edu, 303-497-1853, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
        fopen.write('#CONTACT_INFO: Ortega, Ivan, iortega@ucar.edu, 303-497-1861, NCAR, 3090 Center Green Drive, Boulder, CO 80301\n')
        fopen.write('#Gas name, Version, ctlFile, maxRMS, minDOF:{0}, {1}, {2}, {3}, {4}\n'.format( gasName.upper(), ver, ctlF, maxRMS, minDOF  ))
        if errorFlg: 
            fopen.write('Index, YYYY-MM-DD, HH:MM:SS, TC [molecules cm^-2], Random_Err [molecules cm^-2], Systematic_Err [molecules cm^-2], Total_Err [molecules cm^-2], DOF [a.u], SZA [deg], RMS [%]\n')
            #fopen.write('Index, YYYY-MM-DD, HH:MM:SS, TC [mm], Random_Err [mm], Systematic_Err [mm], Total_Err [mm], DOF [a.u], SZA [deg], RMS [%]\n')
            strFormat = '{0:d}, {1:>10s}, {2:>10s}, {3:.3E}, {4:.3E}, {5:.3E}, {6:.3E}, {7:.3f}, {8:.3f}, {9:.3f}\n'
        else:  
            fopen.write('Index, YYYY-MM-DD, HH:MM:SS, TC [molecules cm^-2], DOF [a.u], SZA [deg], RMS [%]\n')
            strFormat = '{0:d}, {1:>10s}, {2:>10s}, {3:.3E}, {4:.3f}, {5:.3f}, {6:.3f}\n'


        for i,sngTime in enumerate(times):
            if errorFlg: fopen.write(strFormat.format((i+1),YYYYMMDD[i], HHMMSS[i] ,totClmn[i],tot_rnd[i],tot_sys[i],tot_std[i],dofs[i],sza[i],rms[i]))
            #if errorFlg: fopen.write(strFormat.format((i+1),YYYYMMDD[i], HHMMSS[i] ,totClmn[i]*2.989e-22,tot_rnd[i]*2.989e-22,tot_sys[i]*2.989e-22,tot_std[i]*2.989e-22,dofs[i],sza[i],rms[i]))
            else: fopen.write(strFormat.format((i+1),YYYYMMDD[i], HHMMSS[i] ,totClmn[i],dofs[i],sza[i],rms[i]))
    

    #-----------------
    # plot time series
    #-----------------
    years   = YearLocator()
    months  = MonthLocator()
    DateFmt = DateFormatter('%Y-%m')
    #DateFmt      = DateFormatter('%b %d')
    #DateFmt      = DateFormatter('%H:%M')
    
    fig,(ax1,ax2,ax3,ax4) = plt.subplots(4, figsize=(10,10), sharex=True)
    
    #ax1.scatter(dates,totClmn*2.989e-22, facecolors='red', s=30, label='data', color='k')
    if errorFlg: ax1.errorbar(dates,totClmn, yerr=tot_std, fmt='o', markersize=5, color='r', ecolor='r')
    else: ax1.errorbar(dates,totClmn, yerr=totClmn*0., fmt='o', markersize=0, color='r', ecolor='r')
    
    ax1.grid(True)
    ax1.set_ylabel('Retrieved Total Column\n[molecules cm$^{-2}$]', multialignment='center', fontsize=14)
    #ax1.set_ylabel('Total Column [mm]', multialignment='center', fontsize=14)
    #ax1.set_title(gasName.upper() + ' For ' + loc.upper(), fontsize=14 )
    ax1.tick_params(labelsize=14)
    
    #ax2.plot(dates,dofs,'r-o')
    ax2.scatter(dates,dofs, facecolors='red', s=30, label='data', color='k')
    ax2.grid(True)
    ax2.set_ylabel('DOF', fontsize=14)
    ax2.tick_params(labelsize=14)
    
    #ax3.plot(dates,sza,'r-o')
    ax3.scatter(dates,sza, facecolors='red', s=30, label='data', color='k')
    ax3.grid(True)
    ax3.set_ylabel('SZA [Degrees]', multialignment='center', fontsize=14)
    ax3.tick_params(labelsize=14)
    
    #ax4.plot(dates,rms,'r-o')
    ax4.scatter(dates,rms, facecolors='red', s=30, label='data', color='k')
    ax4.grid(True)
    ax4.set_ylabel('RMS', fontsize=14)
    ax4.tick_params(labelsize=14)    
    
    #ax4.xaxis.set_major_locator(years)
    #ax4.xaxis.set_minor_locator(months)
    ax4.xaxis.set_major_formatter(DateFmt)
    ax4.set_xlabel('Date', fontsize=14)
    
    fig.autofmt_xdate()
    plt.suptitle('Time series of {} at {}'.format(gasName.upper(), locStr), fontsize=16)
    fig.subplots_adjust(bottom=0.1, top=0.95, left = 0.13, right = 0.97)

    if saveFlg:     
        pdfsav.savefig(fig,dpi=200)
    else:           
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()

    if saveFlg:     
        pdfsav.close()
            
                                                                                    
if __name__ == "__main__":
    main()