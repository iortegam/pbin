#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        TikOpt.py
#
# Purpose:
#       - Optimize the regularization strength (alpha) used in Tikhonov (SFIT4)
#
# Notes:
#       - See important inputs in main()
#       - This script will run under the folder containing the t15asc.4 and the sfit4.ctl file
#       - Error analysis is optional (See Initializations below)
#       - Layer1Mods, sfitClasses, and dataOutClass classes are needed to plot results
#       - The name of the Sa matrix used in the ctl file needs to be 'Tik.inp'
#       - OPTIONAL = TRANSFORMATION IN CASE IN CASE NON-CONSTANT RETRIEVAL GRID 
#   
# Version History:
#       Created, June, 2017  Ivan Ortega (iortega@ucar.edu)
#
#----------------------------------------------------------------------------------------

                                #-------------------------#
                                #   Import SFIT4 modules  #
                                #-------------------------#
from Layer1Mods import  errAnalysis 
import sfitClasses as sc              
import dataOutClass as dc

                                #-------------------------#
                                # Import standard modules #
                                #-------------------------#

import sys,os
import numpy as np
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

                                #-------------------------#
                                #        Functions        #
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

def getColumns(FileName, headerrow=0, delim=' ', header=True):
    """
    Get columns of data from inFile. The order of the rows is respected
    
    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that 
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in 
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the 
    heading names (i.e. the keys in the cols dict)
    """

    f = open(FileName, 'r')
    inFile = f.readlines()


    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum < headerrow: continue

        if lineNum == headerrow:
            if delim == ' ': headings = line.split()
            else: headings = line.split(delim) 
              
            i = 0
            for heading in headings:
                heading = heading.strip()
         
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            if delim == ' ': cells = line.split()
            else: cells = line.split(delim)

            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1
                
    return cols, indexToName


def tikCov_Trans(gridbound,alpha,pathOut, transflag=True,normalizeflag=True):
    """Create a Tikhonov matrix with strength 1 for a non-equidistant grid

    Input argument: a ndarray of gridboundaries (1D) or a string for the file containing the grid boundaries (increasing height)
    Optional key arguments::
    transflag=True, if true, the matrix is transformed to work on a decreasing grid coordinate system
    normalizeflag=True, if true, the matrix is normalized with mean of the layer thickness (if false, unit is $km^{-2}$);;"""
    ##logger=getlogger(logger,'create_tik')

    if isinstance(gridbound,np.ndarray): 
        
        gridsize=gridbound.shape[0]-1
        boundaries=np.sort(gridbound)
    
    else: 
        
        try: fid=open(gridbound,'r')
        except: 
            return array([])
       
        gridsize=int(fid.readline())
        boundaries=ftir.tools.read_matrixfile(fid)
    
    if boundaries.shape!=(gridsize+1,): 
    	return array([])
    
    l1      = np.zeros((gridsize-1,gridsize))
    mids    = np.array([(boundaries[i+1]+boundaries[i])/2 for i in range(gridsize)])
    
    for i in range(l1.shape[0]): 
        l1[i,i]=-1
        l1[i,i+1]=1

    thickness=np.array([mids[i+1]-mids[i] for i in range(gridsize-1)])
    
    if normalizeflag: thickness/=thickness.mean()
    
    out  = alpha*l1.T.dot(np.diag(thickness**-2)).dot(l1)

 
    if transflag: out=np.eye(gridsize)[::-1,:].dot(out).dot(np.eye(gridsize)[::-1,:])

    np.savetxt(pathOut+'Tik.input',out)

    return out

def tikCov(nlayer, alpha, pathOut):
    '''First Derivative operator'''
    l1=np.zeros((nlayer-1,nlayer),dtype=float)

    for i in range(nlayer-1):
	    l1[i,i]=1.0
	    l1[i,i+1]=-1.0

    out = alpha*np.dot(l1.T,l1)

    np.savetxt(pathOut+'Tik.input',out)

    return out
                                #-------------------------#
                                #         Main            #
                                #-------------------------#

def main():
    	#-----------------------------------------------------------------------------------------
    	#                             Initialization
    	#-----------------------------------------------------------------------------------------
    	loc       = 'tab'                                                                              # LOCATION
        gas       = 'pan'                                                                            # GAS                               
    	alpha     = [0.1, 1, 10, 100, 1000, 10000, 100000, 1000000]                                   # ALPHA VALUES TO TEST
        #alpha     = [10000]#, 1000, 1e4]                                                         # ALPHA VALUES TO TEST
    	TikOut    = '/data1/ebaumer/'+loc.lower()+'/'+gas.lower()+'/x.'+gas.lower()+'/'                # PATH TO SAVE THE TIK MATRIX
        binDir    = '/data/ebaumer/Code/sfit-core-code/src/'                                           # PATH FOR THE SFIT4 SOURCE CODE

        errFlg    = False                                                                               # ERROR ANALYSIS?
        saveFlg   = True                                                                               # SAVE PDF FILE?
        TransFlg  = True                                                                              # OPTIONAL = TRANSFORMATION IN CASE IN CASE NON-CONSTANT RETRIEVAL GRID (e.g., http://www.atmos-meas-tech.net/4/1943/2011/)
        
        
        if saveFlg: pltFile = TikOut + 'TikOpt.pdf'

       
        #-----------------------------------------------------------------------------------------
        #                             START
        #-----------------------------------------------------------------------------------------
        ckDir(TikOut, logFlg=False, exit=True)


        if TransFlg:
        	
            fileLayers = '/data/Campaign/'+loc.upper()+'/local/station.layers'
            ckFile(fileLayers, exit=True)

            cols, indexToName = getColumns(fileLayers, headerrow=2, delim=' ', header=True)
            midpnt   = np.asarray(cols['midpnt'][0:-1]).astype(np.float)
            thick    = np.asarray(cols['thick'][0:-1]).astype(np.float)
            level    = np.asarray(cols['level']).astype(np.float)

        else:

        	if   loc == 'fl0': nlayer  = 44                                                              
        	elif loc == 'mlo': nlayer  = 41
        	elif loc == 'tab': nlayer  = 47
        	else: 
        	    print "nlayer is not known for your location\nModify TikOpt to accomodate nlayer"
        	    exit()


        #-----------------------------------------------------------------------------------------
        #                             Define variable to save
        #-----------------------------------------------------------------------------------------
        totCol   = []
        dof      = []
        rms      = []
        chi2y    = []
        tot_rnd  = []
        tot_sys  = []
        tot_std  = []
        tot_smt  = []
        tot_msr  = []

        #Find Current Directory
        cwd = os.getcwd()

        #-----------------------------------------------------------------------------------------
        #                             Run Sfit
        #-----------------------------------------------------------------------------------------
        for ai in alpha:

            if TransFlg:
            	tikCov_Trans(np.sort(level), ai, TikOut, transflag=True,normalizeflag=False)
            else: 
            	tikCov(nlayer, ai,  TikOut)

            print '************************************************************************'
            print 'Running sfit4 - Tikhonov approach, alpha = {0:.3f}'.format(ai)
            print '************************************************************************'
           
            rtn = sc.subProcRun( [binDir + 'sfit4'] )

            wrkDir    = os.getcwd()
            if not(wrkDir.endswith('/')): wrkDir = wrkDir + '/'

            #-------------------
            # Run error analysis
            #-------------------
            if errFlg:
                print '**********************************************************************************'
                print 'Running error analysis - Tikhonov approach, alpha = {0:.3f}\n'.format(ai)
                print '**********************************************************************************'
                
                ckFile('sfit4.ctl', exit=True)
                ckFile('sb.ctl', exit = True)

                ctlFile = sc.CtlInputFile('sfit4.ctl')
                ctlFile.getInputs()

                sbCtlFile = sc.CtlInputFile('sb.ctl')
                sbCtlFile.getInputs()

                rtn = errAnalysis(ctlFile,sbCtlFile,wrkDir)

            gas = dc.PlotData(wrkDir,wrkDir+'sfit4.ctl')

            gas.readsummary()

            totCol.append(np.asarray(gas.summary[gas.PrimaryGas.upper()+'_RetColmn']))
            rms.append(np.asarray(gas.summary[gas.PrimaryGas.upper()+'_FITRMS'])) 
            dof.append(np.asarray(gas.summary[gas.PrimaryGas.upper()+'_DOFS_TRG']))
            chi2y.append(np.asarray(gas.summary[gas.PrimaryGas.upper()+'_CHI_2_Y']))

            if errFlg:
                gas.readError(totFlg=True,sysFlg=False,randFlg=False,vmrFlg=False,avkFlg=False,KbFlg=False)

                tot_rnd.append(np.array(gas.error['Total random uncertainty']))
                tot_sys.append(np.array(gas.error['Total systematic uncertainty']))
                tot_std.append(np.sqrt((np.array(gas.error['Total random uncertainty']))**2 + (np.array(gas.error['Total systematic uncertainty']))**2))
                tot_smt.append(np.array(gas.error['Smoothing error (Ss)']))          #--> Eric's routine
                #tot_smt.append(np.array(gas.error['Smoothing error (Ss, using sa)'])) #--> Bavo's routine
                tot_msr.append(np.array(gas.error['Measurement error (Sm)']))

        totCol   = np.asarray(totCol)
        dof      = np.asarray(dof)
        rms      = np.asarray(rms)
        chi2y    = np.asarray(chi2y)
        
        if errFlg:
            tot_rnd  = np.asarray(tot_rnd)
            tot_sys  = np.asarray(tot_sys)
            tot_std  = np.asarray(tot_std)
            tot_smt  = np.asarray(tot_smt)
            tot_msr  = np.asarray(tot_msr)

        #-------------------
        # Plots
        #-------------------
        print '**********************'
        print '        Plots        '
        print '**********************'

        if saveFlg: pdfsav = PdfPages(pltFile)
        
        if errFlg: fig, ( (ax, ax2), (ax3, ax4) ) = plt.subplots(2, 2, sharex=True, figsize=(10, 8))
        else: fig, (ax, ax2, ax3) = plt.subplots(3, 1, sharex=True)

        #-------------------Total Column------------------- 
        ax.plot(alpha,totCol,label='totCol')
        ax.scatter(alpha,totCol)
        ax.grid(True,which='both')    
        ax.tick_params(which='both',labelsize=11)
        ax.set_title('Total Column', multialignment='center')
        ax.set_ylabel('Total Column [molec/cm$^2$]')
        plt.xscale('log')

        #-------------------rms------------------- 
        ax2.plot(alpha,rms,label='rms')
        ax2.scatter(alpha,rms)
        ax2.grid(True,which='both')    
        ax2.tick_params(which='both',labelsize=11)
        ax2.set_title('RMS', multialignment='center')
        ax2.set_ylabel('RMS [%]')
        ax2.yaxis.set_major_formatter(FormatStrFormatter('%.4f'))

        #-------------------dof------------------- 
        ax3.plot(alpha,dof,label='dof')
        ax3.scatter(alpha,dof)
        ax3.grid(True,which='both')    
        ax3.tick_params(which='both',labelsize=11)
        ax3.set_title('DOF', multialignment='center')
        ax3.set_ylabel('DOF')
        ax3.set_xlabel('Regularization strength, alpha')

        #-------------------Error------------------- 
        if errFlg:
            ax4.plot(alpha,tot_smt, color='blue', label='smoothing')
            ax4.scatter(alpha,tot_smt, color='blue')
            ax4.plot(alpha,tot_msr, color='green', label='measurement')
            ax4.scatter(alpha,tot_msr, color='green')
            ax4.plot(alpha,np.sqrt(tot_msr**2 + tot_smt**2) , color='gray', label='Total (measurement + smoothing)')
            ax4.scatter(alpha,np.sqrt(tot_msr**2 +tot_smt**2), color='gray')
        
            ax4.grid(True,which='both')    
            ax4.tick_params(which='both',labelsize=11)
            ax4.set_title('Error', multialignment='center')
            ax4.set_ylabel('Error [%]')
            ax4.set_xlabel('Regularization strength, alpha')

            ax4.legend(prop={'size':10}, loc = 1)

        fig.suptitle(cwd)
        fig.tight_layout()

        fig.subplots_adjust(left=0.1, right=0.95, top=0.9)

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            pdfsav.close()
        else:       
            plt.show(block=False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()


if __name__ == "__main__":
    main()