#!/usr/bin/python
##! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        pltvs.py
#
# Purpose:
#
# Plot time series of columns/VMRs/error and linear correlation analyis for comparison
# Usefule when comparing different settings in the retrieval of trace gase
# Input:
#       input_pltvs.py
#----------------------------------------------------------------------------------------


#---------------
# Import modules
#---------------
import sys
import os
import getopt
import dataOutplts as dc
import time
import datetime as dt

import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 

import myfunctions as mf
import numpy as np
from scipy import stats



def usage():
    ''' Prints to screen standard program usage'''
    print 'pltSet.py -i <inputfile> -?'
    print '  -i <file> : Run pltSet.py with specified input file'
    print '  -?        : Show all flags'

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

def main(argv):

    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:?')

    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit()

    #-----------------------------
    # Parse command line arguments
    #-----------------------------
    for opt, arg in opts:
        # Check input file flag and path
        if opt == '-i':

            pltInputs = {}

            ckFile(arg,exit=True)

            try:
                execfile(arg, pltInputs)
            except IOError as errmsg:
                print errmsg + ' : ' + arg
                sys.exit()

            if '__builtins__' in pltInputs:
                del pltInputs['__builtins__']               


        # Show all command line flags
        elif opt == '-?':
            usage()
            sys.exit()

        else:
            print 'Unhandled option: ' + opt
            sys.exit()

    ngases = len(pltInputs['gasName'])

    #--------------------------------
    # Define the Retrieval arguments
    #--------------------------------
    TC       =  []  #TOTAL COLUMN OF MAIN TRACE GAS
    TCApr    =  []  #TOTAL COLUMN OF APRIORI OF TRACE GAS
    DT       =  []  #DATE and TIME
    DOF      =  []  #DEGREES OF FREEDOM
    RMS      =  []  #RMS
    TCRE     =  []  #TOTAL COLUMN RANDOM ERROR
    TCSE     =  []  #TOTAL COLUMN SYSTEMATIC ERROR
    TCTE     =  []  #TOTAL COLUMN TOTAL ERROR
    ID       =  []  #ID version or NAME
    TCall    =  []  #TOTAL COLUMN OF other gases retrieved
    TCallApr =  []  #TOTAL COLUMN OF APRIORI other gases retrieved
    GasList  =  []
    chi2y    =  []



    for k in range(ngases):

        retDir = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ver'][k]+'/'
        ctlFile  = '/data1/ebaumer/'+pltInputs['loc'].lower()+'/'+pltInputs['gasName'][k].lower()+'/'+'x.'+pltInputs['gasName'][k].lower()+'/'+pltInputs['ctlF'][k]
        maxRMS = pltInputs['maxRMS'][k]
        minDOF = pltInputs['minDOF'][k]

    #---------------------------------
    # Check for the existance of files 
    # directories from input file
    #---------------------------------

        ckDir(retDir, exit=True)
        ckFile(ctlFile, exit=True)

        if pltInputs['saveFlg']:  ckDir(os.path.dirname(os.path.realpath(pltInputs['pltFile'])),exit=True)

    #-------------------------
    # Create Instance of Class
    #-------------------------
        gas = dc.PlotData(retDir,ctlFile,iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
        fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'],outFname=pltInputs['pltFile'])

    #----------------------
    # Call to plot profiles
    #----------------------
        ds = gas.pltTotClmn(fltr=pltInputs['fltrFlg'],sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
                   partialCols=pltInputs['pCols'],errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
                   maxSZA=pltInputs['maxSZA'],maxRMS=maxRMS,minDOF=minDOF,maxCHI=pltInputs['maxCHI'],
                   dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
                   pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],chiFlg=pltInputs['chiFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],tcMMflg=pltInputs['tcMMFlg'],
                   allGas = True)


        TC.append(ds['totClmn'])
        TCApr.append(ds['totClmnApr'])
        DT.append(ds['dates'])
        DOF.append(ds['dofs'])
        RMS.append(ds['rms'])
        TCRE.append(ds['tot_rnd'])
        TCSE.append(ds['tot_sys'])
        TCTE.append(ds['tot_std'])
        ID.append(pltInputs['ID'][k])
        TCall.append(ds['TCret'])
        TCallApr.append(ds['TCapr'])
        GasList.append(ds['GasList'])
        chi2y.append(ds['chi2y'])



    clmap        = 'jet'
    cm           = plt.get_cmap(clmap)    
    yearsLc      = YearLocator()
    monthsAll    = MonthLocator()
    months       = MonthLocator()
   
    #DateFmt      = DateFormatter('%m\n%Y')
    DateFmt      = DateFormatter('%Y')

    colors = ['r', 'b', 'g', 'y', 'c', 'k']
    gas4ratio = pltInputs['g4ratio'].upper()

    if pltInputs['errorFlg']: fig, ax = plt.subplots(7,figsize=(8,15), sharex=True)
    else: fig, ax = plt.subplots(5,figsize=(8,10), sharex=True)
    
    
    for k in range(ngases):
        ax[0].scatter(DT[k], TC[k],facecolors='white', s=35, color=colors[k],label=ID[k] )
        ax[0].errorbar(DT[k], TC[k],yerr=TCTE[k],fmt='k.',markersize=0, ecolor=colors[k])  
        ax[0].legend(prop={'size':10})

        ax[1].scatter(DT[k], RMS[k],facecolors='white', s=35, color=colors[k])
        ax[2].scatter(DT[k], DOF[k],facecolors='white', s=35, color=colors[k])
        ax[3].scatter(DT[k], chi2y[k],facecolors='white', s=35, color=colors[k])
        ax[4].scatter(DT[k], TCall[k][gas4ratio]/TCallApr[k][gas4ratio],facecolors='white', s=35, color=colors[k])
        
        if pltInputs['errorFlg']:
            ax[5].scatter(DT[k], TCRE[k],facecolors='white', s=35, color=colors[k])
            ax[6].scatter(DT[k], TCSE[k],facecolors='white', s=35, color=colors[k])
            #ax[5].scatter(DT[k], TCTE[k],facecolors='white', s=35, color=colors[k])
        
    ax[0].set_ylim(np.min(TC[0])-0.1*np.min(TC[0]), np.max(TC[0])+0.1*np.max(TC[0]))
    ax[0].set_ylabel('Total Column\n[molecules$\cdot$cm$^{-2}$]')

    #ax[1].set_ylim(np.min(RMS[2])-0.5*np.min(RMS[2]), np.max(RMS[2])+0.3*np.max(RMS[2]))
    ax[1].set_ylim(0, np.max(RMS[0])+0.3*np.max(RMS[0]))
    ax[1].set_ylabel('RMS')

    ax[2].set_ylim(np.min(DOF[0])-0.5*np.min(DOF[0]), np.max(DOF[0])+0.3*np.max(DOF[0]))
    #ax[2].set_ylim(0, 2.5)
    ax[2].set_ylabel('DOF')

    ax[3].set_ylim(np.min(chi2y[0])-0.5*np.min(chi2y[0]), np.max(chi2y[0])+0.5*np.max(chi2y[0]))
    ax[3].set_ylabel('$\chi_y^{2}$')

    #ax[6].set_ylim(np.min(TCTE[0])-0.5*np.min(TCTE[0]), np.max(TCTE[0])+0.3*np.max(TCTE[0]))
    ax[4].set_ylabel(gas4ratio+'(ret)/'+gas4ratio+'(apr)')

    if pltInputs['errorFlg']:
        ax[5].set_ylim(np.min(TCRE[0])-0.5*np.min(TCRE[0]), np.max(TCRE[0])+0.5*np.max(TCRE[0]))
        ax[5].set_ylabel('Random error\n[molecules$\cdot$cm$^{-2}$]')

        ax[6].set_ylim(np.min(TCSE[0])-0.5*np.min(TCSE[0]), np.max(TCSE[0])+0.5*np.max(TCSE[0]))
        ax[6].set_ylabel('Systematic error\n[molecules$\cdot$cm$^{-2}$]')

        #ax[5].set_ylim(np.min(TCTE[0])-0.5*np.min(TCTE[0]), np.max(TCTE[0])+0.5*np.max(TCTE[0]))
        #ax[5].set_ylabel('Total error\n[molecules$\cdot$cm$^{-2}$]')

    if pltInputs['errorFlg']: nf = 7
    else: nf = 5
    for p in range(0,nf):
        ax[p].grid(True)
            
        if ds['yrsFlg']:
            #plt.xticks(rotation=45)
            ax[p].xaxis.set_major_locator(yearsLc)
            ax[p].xaxis.set_minor_locator(months)
            ax[p].xaxis.set_major_formatter(DateFmt) 
            ax[p].xaxis.set_tick_params(which='major',labelsize=12)
            ax[p].xaxis.set_tick_params(which='minor',labelbottom='off')
            ax[p].set_xlim((dt.date(pltInputs['iyear'],pltInputs['imnth'], 1), dt.date(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'])))
            if k == ngases-1: ax[p].set_xlabel('Year', fontsize = 16)
        
        elif ds['mntsFlg']:
            ax[p].xaxis.set_major_locator(monthsAll)
            ax[p].xaxis.set_major_formatter(DateFmt)
            #ax[p].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax[p].set_xlim((dt.datetime(pltInputs['iyear'],pltInputs['imnth'], pltInputs['iday'], 14), dt.datetime(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'], 23)))
            ax[p].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax[p].set_xlabel('Month/Year')

        elif ds['mdays']:
            ax[p].xaxis.set_major_locator(monthsAll)
            ax[p].xaxis.set_major_formatter(DateFmt)
            #ax[p].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax[p].set_xlim((dt.datetime(pltInputs['iyear'],pltInputs['imnth'], pltInputs['iday'], 14), dt.datetime(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'], 23)))
            ax[p].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax[p].set_xlabel('Month/Year')

        else:
            #ax[p].xaxis.set_major_locator(monthsAll)
            ax[p].xaxis.set_major_formatter(DateFormatter('%H:%M'))
            #ax[p].set_xlim((dt.date(ds['years'][0],1,1), dt.date(ds['years'][0],12,31)))
            ax[p].set_xlim((dt.datetime(pltInputs['iyear'],pltInputs['imnth'], pltInputs['iday'], 14), dt.datetime(pltInputs['fyear'],pltInputs['fmnth'],pltInputs['fday'], 23)))
            ax[p].xaxis.set_minor_locator(AutoMinorLocator())
            if k == ngases-1: ax[p].set_xlabel('Month/Year')
    
    fig.autofmt_xdate()
    plt.tight_layout()
    fig.subplots_adjust(bottom=0.07, top=0.97)
    #fig.subplots_adjust(left=0.07, bottom=0.1, top=0.95, right= 0.95)        
    
    #Linear Correlation analysis
    for k in range(ngases):
        if pltInputs['ID'][k] == pltInputs['ver_vs']: 
            indT = k
        else: indT = 0
    #if not indT: indT = 0
     
    dtfc_fix     = dc.toYearFraction(DT[indT])
    tc_fic       = TC[indT]
    rms_fic      = RMS[indT]
    DOF_fic      = DOF[indT]
    TCRE_fic     = TCRE[indT]
    TCSE_fic     = TCSE[indT]
    ratio_fic    = TCall[indT][gas4ratio]/TCallApr[indT][gas4ratio]

    #one2one = np.arange(np.min(TC[indT]) - 1.0*np.min(TC[indT]), np.max(TC[indT])+1.0*np.max(TC[indT]), np.min(TC[indT])*0.05)
    one2one = np.arange(0, 1e17, 1e15)

    if pltInputs['errorFlg']: fig2, ax = plt.subplots(2,3, figsize=(9,8))
    else: fig2, ax = plt.subplots(2,2, figsize=(9,8))

    for k in range(ngases):

        if k != indT:

            dtfc = dc.toYearFraction(DT[k])
            tci  = np.interp(dtfc_fix, dtfc, TC[k])
            rmsi = np.interp(dtfc_fix, dtfc, RMS[k])
            dofi = np.interp(dtfc_fix, dtfc, DOF[k])
            ratioi = np.interp(dtfc_fix, dtfc, TCall[k][gas4ratio]/TCallApr[k][gas4ratio])

            if pltInputs['errorFlg']: 
                tcei = np.interp(dtfc_fix, dtfc, TCTE[k])
                tcrei = np.interp(dtfc_fix, dtfc, TCRE[k])
                tcsei = np.interp(dtfc_fix, dtfc, TCSE[k])
            
           
            ax[0,0].scatter(tc_fic,  tci, facecolors='white', s=35, color=colors[k], label = ID[k]) 
            #-----linear correlation analysis
            #A = np.vstack([tc_fic, np.ones(len(tc_fic))]).T
            #slope, intercept = np.linalg.lstsq(A, tci)[0]
            #fit = slope*np.array(tc_fic) + intercept

            # Fit
            slope, intercept, r_value, p_value, std_err = stats.linregress(tc_fic, tci)


            #ax[0,0].plot(tc_fic, fit, 'k', linewidth=2, label='Fitted line')
            if ngases ==2: ax[0,0].text(0.5,0.15, "slope: {0:4.2f}".format(slope),transform=ax[0,0].transAxes, fontsize = 11)
            if ngases ==2: ax[0,0].text(0.5,0.1,  "Intercept: {0:.2E}".format(intercept),transform=ax[0,0].transAxes, fontsize = 11)

            print 'Linerar Correlation: '+ pltInputs['ID'][k] + 'vs' + pltInputs['ver_vs']
            print 'Total Column:'
            print "slope: {0:4.2f}".format(slope)
            print "Intercept: {0:.2E}".format(intercept)
            print "r-value: {0:.2E}".format(r_value)

            #---------------------
            #Bias
            #---------------------
            Bias = mf.bias(tc_fic, tci)
            print 'Bias: {0:.2E} (molec/cm$^2$)'.format(Bias)
            print 'Bias [%]: {0:.2f} '.format(Bias/np.mean(tc_fic))
            if ngases ==2: ax[0,0].text(0.5,0.05,  "Bias: {0:.2E} ".format(Bias),transform=ax[0,0].transAxes, fontsize = 11)


            #---------------
            ax[0,1].scatter(rms_fic,  rmsi, facecolors='white', s=35, color=colors[k])
            #-----linear correlation analysis
            #A = np.vstack([rms_fic, np.ones(len(rms_fic))]).T
            #slope, intercept = np.linalg.lstsq(A, rmsi)[0]
            #fit = slope*np.array(rms_fic) + intercept

            slope, intercept, r_value, p_value, std_err = stats.linregress(rms_fic, rmsi)


            #ax[0,1].plot(rms_fic, fit, 'k', linewidth=2, label='Fitted line')
            if ngases ==2: ax[0,1].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[0,1].transAxes, fontsize = 11)
            if ngases ==2: ax[0,1].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[0,1].transAxes, fontsize = 11)
            print 'RMS:'
            print "slope: {0:4.2f}".format(slope)
            print "Intercept: {0:.2E}".format(intercept)
            print "r-value: {0:.2E}".format(r_value)
            #---------------

            if pltInputs['errorFlg']: nr = (0,2)
            else: nr  = (1,0)

            ax[nr].scatter(DOF_fic,  dofi, facecolors='white', s=35, color=colors[k])
            #-----linear correlation analysis
            #A = np.vstack([DOF_fic, np.ones(len(DOF_fic))]).T
            #slope, intercept = np.linalg.lstsq(A, dofi)[0]
            #fit = slope*np.array(DOF_fic) + intercept

            slope, intercept, r_value, p_value, std_err = stats.linregress(DOF_fic, dofi)


            #ax[0,2].plot(DOF_fic, fit, 'k', linewidth=2, label='Fitted line')
            if ngases ==2: ax[nr].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[nr].transAxes, fontsize = 11)
            if ngases ==2: ax[nr].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[nr].transAxes, fontsize = 11)
            print 'DOF:'
            print "slope: {0:4.2f}".format(slope)
            print "Intercept: {0:.2E}".format(intercept)
            print "r-value: {0:.2E}".format(r_value)


        

            if pltInputs['errorFlg']: nr = (1,0)
            else: nr  = (1,1)

            #---------------
            ax[nr].scatter(ratio_fic,  ratioi, facecolors='white', s=35, color=colors[k])
            #-----linear correlation analysis
            #A = np.vstack([ratio_fic, np.ones(len(ratio_fic))]).T
            #slope, intercept = np.linalg.lstsq(A, ratioi)[0]  
            #fit = slope*np.array(ratio_fic) + intercept

            slope, intercept, r_value, p_value, std_err = stats.linregress(ratio_fic, ratioi)
            #ax[1,2].plot(ratio_fic, fit, 'k', linewidth=2, label='Fitted line')
            if ngases ==2: ax[nr].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[nr].transAxes, fontsize = 11)
            if ngases ==2: ax[nr].text(0.4,0.04,  "Intercept: {0:4.3f}".format(intercept),transform=ax[nr].transAxes, fontsize = 11) 


            if pltInputs['errorFlg']:

                ax[1,1].scatter(TCRE_fic,  tcrei, facecolors='white', s=35, color=colors[k])
                #-----linear correlation analysis
                #A = np.vstack([TCRE_fic, np.ones(len(TCRE_fic))]).T
                #slope, intercept = np.linalg.lstsq(A, tcrei)[0]  
                #fit = slope*np.array(TCRE_fic) + intercept

                slope, intercept, r_value, p_value, std_err = stats.linregress(TCRE_fic, tcrei)
                #ax[1,0].plot(TCRE_fic, fit, 'k', linewidth=2, label='Fitted line')
                if ngases ==2: ax[1,1].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[1,1].transAxes, fontsize = 11)
                if ngases ==2: ax[1,1].text(0.4,0.04,  "Intercept: {0:.2E}".format(intercept),transform=ax[1,1].transAxes, fontsize = 11) 
                
                #---------------
                ax[1,2].scatter(TCSE_fic,  tcsei, facecolors='white', s=35, color=colors[k])
                #-----linear correlation analysis
                #A = np.vstack([TCSE_fic, np.ones(len(TCSE_fic))]).T
                #slope, intercept = np.linalg.lstsq(A, tcsei)[0]  
                #fit = slope*np.array(TCSE_fic) + intercept

                slope, intercept, r_value, p_value, std_err = stats.linregress(TCSE_fic, tcsei)
                #ax[1,1].plot(TCSE_fic, fit, 'k', linewidth=2, label='Fitted line')
                if ngases ==2: ax[1,2].text(0.4,0.1, "slope: {0:4.2f}".format(slope),transform=ax[1,2].transAxes, fontsize = 11)
                if ngases ==2: ax[1,2].text(0.4,0.04,  "Intercept: {0:.2E}".format(intercept),transform=ax[1,2].transAxes, fontsize = 11) 
            
            #---------------

    
    ax[0,0].set_xlabel('Total Column [molec$\cdot$cm$^{-2}$] - %s' %str(pltInputs['ver_vs']))
    ax[0,0].set_ylabel('Total Column [molec$\cdot$cm$^{-2}$]')    
    ax[0,0].legend(prop={'size':12})
    ax[0,0].set_ylim(np.min(TC[indT])-0.05*np.min(TC[indT]), np.max(TC[indT])+0.05*np.max(TC[indT]))
    ax[0,0].set_xlim(np.min(TC[indT])-0.05*np.min(TC[indT]), np.max(TC[indT])+0.05*np.max(TC[indT]))
    ax[0,0].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    ax[0,0].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    ax[0,0].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    ax[0,0].grid(True)

    ax[0,1].set_xlabel('RMS - %s' %str(pltInputs['ver_vs']))
    ax[0,1].set_ylabel('RMS')    
    ax[0,1].set_ylim(np.min(RMS[indT])-0.5*np.min(RMS[indT]), np.max(RMS[indT])+0.15*np.max(RMS[indT]))
    ax[0,1].set_xlim(np.min(RMS[indT])-0.5*np.min(RMS[indT]), np.max(RMS[indT])+0.15*np.max(RMS[indT]))
    ax[0,1].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    ax[0,1].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    ax[0,1].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    ax[0,1].grid(True)

    if pltInputs['errorFlg']: nr = (0,2)
    else: nr = (1,0)
    ax[nr].set_xlabel('DOF - %s' %str(pltInputs['ver_vs']))
    ax[nr].set_ylabel('DOF')    
    ax[nr].set_ylim(np.min(DOF[indT])-0.5*np.min(DOF[indT]), np.max(DOF[indT])+0.15*np.max(DOF[indT]))
    ax[nr].set_xlim(np.min(DOF[indT])-0.5*np.min(DOF[indT]), np.max(DOF[indT])+0.15*np.max(DOF[indT]))
    ax[nr].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    ax[nr].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    ax[nr].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    ax[nr].grid(True)

    #ax[1,2].set_xlabel('H$_2$O (ret)/H$_2$O (apr) - %s' %str(pltInputs['ver_vs']))
    #ax[1,2].set_ylabel('H$_2$O (ret)/H$_2$O (apr)')
    if pltInputs['errorFlg']: nr = (1,0)
    else: nr = (1,1)
    ax[nr].set_xlabel(gas4ratio+' (ret)/'+gas4ratio+' (apr) - %s' %str(pltInputs['ver_vs']))
    ax[nr].set_ylabel(gas4ratio+' (ret)/'+gas4ratio+' (apr)')    
    ax[nr].set_ylim(np.min(ratio_fic)-0.5*np.min(ratio_fic), np.max(ratio_fic)+0.15*np.max(ratio_fic))
    ax[nr].set_xlim(np.min(ratio_fic)-0.5*np.min(ratio_fic), np.max(ratio_fic)+0.15*np.max(ratio_fic))
    ax[nr].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
    ax[nr].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
    ax[nr].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
    ax[nr].grid(True)

    if pltInputs['errorFlg']:

        ax[1,1].set_xlabel('Random Error - %s' %str(pltInputs['ver_vs']))
        ax[1,1].set_ylabel('Random Error')    
        ax[1,1].set_ylim(np.min(TCRE[indT])-0.5*np.min(TCRE[indT]), np.max(TCRE[indT])+0.15*np.max(TCRE[indT]))
        ax[1,1].set_xlim(np.min(TCRE[indT])-0.5*np.min(TCRE[indT]), np.max(TCRE[indT])+0.15*np.max(TCRE[indT]))
        ax[1,1].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
        ax[1,1].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
        ax[1,1].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
        ax[1,1].grid(True)

        ax[1,2].set_xlabel('Systematic Error - %s' %str(pltInputs['ver_vs']))
        ax[1,2].set_ylabel('Systematic Error')    
        ax[1,2].set_ylim(np.min(TCSE[indT])-0.5*np.min(TCSE[indT]), np.max(TCSE[indT])+0.15*np.max(TCSE[indT]))
        ax[1,2].set_xlim(np.min(TCSE[indT])-0.5*np.min(TCSE[indT]), np.max(TCSE[indT])+0.15*np.max(TCSE[indT]))
        ax[1,2].plot(one2one, one2one, ls ='--', color='gray', linewidth=2)
        ax[1,2].plot(one2one, one2one*0.8, ls ='--', color='gray', linewidth=2)
        ax[1,2].plot(one2one*0.8, one2one, ls ='--', color='gray', linewidth=2) 
        ax[1,2].grid(True)

    
    plt.tight_layout()
    fig2.subplots_adjust(left=0.075, bottom=0.1, top=0.95, right= 0.95)

    if gas.pdfsav: 
        gas.pdfsav.savefig(fig,dpi=200)
        gas.pdfsav.savefig(fig2,dpi=200)
       
    else:           
        plt.show(block=False)  


    #-------------------------------------------------------------------------
    #END PLOTS
    #-------------------------------------------------------------------------
    if pltInputs['saveFlg']: gas.closeFig()

    print('\nFinished Plots.......\n')
       
      #--------------------------------
      # Pause so user can look at plots
      #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program        


if __name__ == "__main__":
    main(sys.argv[1:])
