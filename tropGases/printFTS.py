#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         print FTS.py
#
# Purpose:
#         Read and Create a .dat file for each input gas (see ianputs below)
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

    #loc                = 'tab'
    #gasName            = ['co',            'c2h6',         'hcn']             
    #ver                = ['Current_B3',    'Current_v2',   'Current_WP' ]          # Name of retrieval version to process
    #ctlF               = ['sfit4_v3.ctl','sfit4_v2.ctl', 'sfit4.ctl'] 

    #----------------------
    # First Level Retrieval Directory
    #----------------------
    retDir             = '/data1/ebaumer/'+loc.lower()

    #dataPath           = '/data/iortega/pbin/tropGases/data/'
    dataPath           = '/data/iortega/results/'+loc.lower()+'/data/'
 
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

    mnths              = [6, 7, 8]                # Months to filter on (these are the months to include data)
    maxRMS             = [1.5, 0.6, 1.0, 0.3, 1.0, 1.5, 0.5, 1.5]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF             = [1.0, 0.5, 0.5, 0.9, 0.5, 0.5,  0.5, 0.5]                      # Min DOFs for filtering

    #maxRMS             = [2.5, 2.2, 1.5]                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    #minDOF             = [1.0, 0.5, 0.9]  

    minSZA             = 0.0                    # Min SZA for filtering
    maxSZA             = 90.0                   # Max SZA for filtering
    maxCHI             = 2.0                    # Max CHI_y_2 value
    maxTC              = 5.0E24                 # Max Total column amount for filtering
    minTC              = 0.0                    # Min Total column amount for filtering
    sclfct             = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
    sclfctName         = 'ppbv'                 # Name of scale factor for labeling plots

    pColsFlg           = True
    pCols              = [1.6, 8.0]         #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR

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
    

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------
   
    fts = cfts.FTSClass( gasName, retDir,  ctlF, ver, iyear,imnth, iday, fyear, fmnth, fday)

    fts.ReadFTS(fltrFlg=fltrFlg, sclfct=sclfct, sclname=sclfctName, mnthFltr= mnths, mnthFltFlg=mnthFlg,
          errFlg=errorFlg, minSZA=minSZA, maxSZA=maxSZA, maxRMS=maxRMS, minTC=minTC, maxTC= maxTC,
          minDOF=minDOF, maxCHI=maxCHI, dofFlg=dofFlg, rmsFlg=rmsFlg, tcFlg=tcNegFlg,
          pcFlg=pcNegFlg, szaFlg=szaFlg,cnvrgFlg=cnvrgFlg, chiFlg=chiFlg,tcMMflg=tcMMFlg, pColsFlg=pColsFlg, pCols=pCols)


    
     
    #---------------------------------
    #  -- Define Global Parameters --
    #---------------------------------
    GasVer = []

    for (g,v) in izip(gasName, ver):
        GasVer.append(g+'_'+v)

    for gv in GasVer:

        print '***************************************'
        print 'Error Analysis for {}'.format(gv)
        print '***************************************'

        print 'Sys Error VMR: {:.3f}'.format(np.mean(fts.vmrP_Syserr[gv]))
        print 'Ran Error VMR: {:.3f}'.format(np.mean(fts.vmrP_Ranerr[gv]))
        print 'Tot Error VMR: {:.3f}'.format(np.mean(fts.vmrP_Toterr[gv]))

        print 'Sys Error %:   {:.1f}'.format(np.mean(fts.vmrP_Syserr[gv])/np.mean(fts.vmrP[gv]) * 100.  )
        print 'Ran Error %:   {:.1f}'.format(np.mean(fts.vmrP_Ranerr[gv])/np.mean(fts.vmrP[gv]) * 100.  )
        print 'Tot Error %:   {:.1f}'.format(np.mean(fts.vmrP_Toterr[gv])/np.mean(fts.vmrP[gv]) * 100.  )

        print 'Mean DOF:      {:.3f}'.format(np.mean(fts.dofs[gv]))

    #--------------

    for gv in GasVer:
        
        YYYYMMDD = np.asarray(['{0:4d}-{1:02d}-{2:02d}'.format(d.year, d.month, d.day)    for d in fts.dates[gv]])
        HHMMSS   = np.asarray(['{0:02d}:{1:02d}:{2:02d}'.format(d.hour,d.minute,d.second) for d in fts.dates[gv]])

        with open(dataPath + gv+'.dat','w') as fopen:
        
            if errorFlg: 
                fopen.write('Index, YYYY-MM-DD, HH:MM:SS, wVMR, TC, TC_rand_err, TC_sys_err, TC_tot_err, SZA, RMS\n')
                #fopen.write('Index, YYYY-MM-DD, HH:MM:SS, TC [mm], Random_Err [mm], Systematic_Err [mm], Total_Err [mm], DOF [a.u], SZA [deg], RMS [%]\n')
                strFormat = '{0:d}, {1:>10s}, {2:>10s}, {3:.3E}, {4:.3E}, {5:.3E}, {6:.3E}, {7:.3f}, {8:.3f}, {9:.3f}\n'
            else:  
                fopen.write('Index, YYYY-MM-DD, HH:MM:SS, wVMR, TC,  SZA, RMS\n')
                strFormat = '{0:d}, {1:>10s}, {2:>10s}, {3:.3E}, {4:.3E}, {5:.3f}, {6:.3f}\n'


            for i, sngTime in enumerate(YYYYMMDD):
                if errorFlg: 
                
                    fopen.write(strFormat.format((i+1),YYYYMMDD[i], HHMMSS[i] , fts.vmrP[gv][i], fts.totClmn[gv][i], fts.tot_rnd[gv][i],fts.tot_sys[gv][i],fts.tot_std[gv][i],fts.sza[gv][i],fts.rms[gv][i]))
                #if errorFlg: fopen.write(strFormat.format((i+1),YYYYMMDD[i], HHMMSS[i] ,totClmn[i]*2.989e-22,tot_rnd[i]*2.989e-22,tot_sys[i]*2.989e-22,tot_std[i]*2.989e-22,dofs[i],sza[i],rms[i]))
                else: fopen.write(strFormat.format((i+1), YYYYMMDD[i], HHMMSS[i], fts.vmrP[gv][i], fts.totClmn[gv][i], fts.sza[gv][i],fts.rms[gv][i] ) )

   


 
if __name__ == "__main__":
    main()