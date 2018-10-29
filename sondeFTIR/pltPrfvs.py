#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltPrfvs.py
#
# Purpose:
#
# Plot comparison of Profiles using different versions
#
# Input:
#       input_pltPrfvs.py
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
import classSondeFTS as rs
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
from scipy import linspace, polyval, polyfit, sqrt, stats, randn

import pylab as P

from scipy import interpolate

from math import acos, asin, atan2, cos, hypot
from math import degrees, pi as PI, radians, sin, tan


def usage():
    ''' Prints to screen standard program usage'''
    print 'pltPrfvs.py -i <inputfile> -?'
    print '  -i <file> : Run pltPrfvs.py with specified input file'
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

def getidver(version):

    labstr = version.strip().split('/')
    if labstr[1] == 'Current_ERA': idver = 'ERA-d'
    if labstr[1] == 'Current_ERA_v66': idver = 'ERA-6' 
    if labstr[1] == 'Current_WACCM':  idver = 'WACCM'
    if labstr[1] == 'Current_NCEP': idver = 'NCEP-d'

    return idver

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


    if pltInputs['saveFlg']: pdfsav = PdfPages(pltInputs['pltFile'])
    

    #---------------------------------
    #
    #     -- Read FTS --
    #
    #---------------------------------
    
    fts = rs.FTSClass( pltInputs['gasName'], pltInputs['retDir'],  pltInputs['ctlF'], pltInputs['loc'], pltInputs['ver'], iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
    fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'], incr=1)

    fts.ReadFTS(fltrFlg=pltInputs['fltrFlg'],allGas=False,sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
           errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=pltInputs['maxRMS'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
           minDOF=pltInputs['minDOF'],maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
           pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],chiFlg=pltInputs['chiFlg'],tcMMflg=pltInputs['tcMMFlg'])

    clr     = mf.clrplt()
    vers    = pltInputs['ver']
    loc     = pltInputs['loc']
    idvers  = pltInputs['ID']
    pCols   = pltInputs['pCols']

    #---------------------------------
    #Defining the reference profile
    #---------------------------------
    for i, d in enumerate(idvers):
        if d == pltInputs['ref']: 
            indT = i
        else: indT = 0
     
    dt_ref       = fts.dates[vers[indT]]   
    doy_ref      = dc.toYearFraction(dt_ref)
    rPrf_ref     = fts.rPrfVMR[vers[indT]]
    if pltInputs['errorFlg']: rPrf_ref_e   = fts.tot_errvmr[vers[indT]]
    alt          = fts.alt[vers[indT]]
    rms_ref      = fts.rms[vers[indT]]
    dof_ref      = fts.dofsAvg[vers[indT]]

    bias         = {}
    bias_e       = {}

    prec         = {}
    prec_e       = {}

    bias_perc    = {}
    bias_perc_e  = {}

    prec_perc    = {}
    prec_perc_e  = {}

    slope        = {}
    slope_e      = {}
    
    intercept    = {}
    intercept_e  = {}
    
    rvalue       = {}

    rms          = {}
    dof          = {}

    for k in range(len(vers)):

        if k != indT:



            dt_i     = fts.dates[vers[k]]   
            doy_i    = dc.toYearFraction(dt_i)
            rPrf_i   = fts.rPrfVMR[vers[k]]
            if pltInputs['errorFlg']: rPrf_i_e   = fts.tot_errvmr[vers[k]]
            rms_i      = fts.rms[vers[k]]
            dof_i      = fts.dofsAvg[vers[k]] 

            intrsctVals = np.intersect1d(doy_ref, doy_i, assume_unique=False)
            
            inds1       = np.nonzero( np.in1d( doy_i, intrsctVals, assume_unique=False ) )[0]
            inds2       = np.nonzero( np.in1d( doy_ref, intrsctVals, assume_unique=False ) )[0]

            print '\nTotal Number of coincident dates between {} and {} = {}'.format(pltInputs['ref'], idvers[k], str(len(intrsctVals)))

            odr, odrErr  = mf.orthoregress(rms_ref[inds2], rms_i[inds1], xerr=rms_ref[inds2]*0., yerr=rms_i[inds1]*0., InError=True)
            print 'RMS - slope +/- intercept: {0:.3f} +/- {1:.3f}'.format(float(odr[0]), float(odr[1]))
            #rms.setdefault(k, []).append(prec_n)
            

            for p, pcol in enumerate(pCols[0:-1]):

                pcolstr = str(pcol[0])

                print 'Column: '+pcolstr

                indsH = np.where( (alt >= pcol[0]) & (alt <pcol[1]))[0]

                rPrf_i2       = rPrf_i[inds1, :]
                rPrf_i2       = rPrf_i2[:, indsH]

                rPrf_ref2     = rPrf_ref[inds2, :]
                rPrf_ref2     = rPrf_ref2[:,indsH]

                if pltInputs['errorFlg']: 

                    rPrf_i2_e     = rPrf_i_e[inds1, :]
                    rPrf_i2_e     = rPrf_i2_e[:, indsH]

                    rPrf_ref2_e   = rPrf_ref_e[inds2, :]
                    rPrf_ref2_e   = rPrf_ref2_e[:,indsH]

                #--------------
                #Bias and precision in Retrieval
                #--------------
                bias_n = np.sum(rPrf_i2 - rPrf_ref2, axis=1 )/float(len(indsH))
                bias_n = bias_n[~np.isnan(bias_n)]

                mu = np.nanmean(bias_n)
                me = np.nanmedian(bias_n)
                sigma = np.std(bias_n)
                stdE = sigma/np.sqrt(len(bias_n))   
                prec_n = sigma/np.sqrt(len(bias_n)) * 2.

                print 'Mean Bias (vmr)   = {0:.5f} +/- {1:.5f}'.format(me, stdE)
                print 'Mean Bias (%)     = {0:.5f} +/- {1:.5f}'.format(me/np.nanmean(rPrf_ref2) * 100., stdE/np.nanmean(rPrf_ref2) * 100.)
                bias.setdefault(k, []).append(me)
                bias_e.setdefault(k, []).append(stdE)

                bias_perc.setdefault(k, []).append(me/np.nanmean(rPrf_ref2) * 100.)
                bias_perc_e.setdefault(k, []).append(stdE/np.nanmean(rPrf_ref2) * 100.)

                prec.setdefault(k, []).append(prec_n)
                prec_e.setdefault(k, []).append(stdE * 0.71)

                prec_perc.setdefault(k, []).append(prec_n/np.nanmean(rPrf_ref2) * 100.)
                prec_perc_e.setdefault(k, []).append((stdE * 0.71)/np.nanmean(rPrf_ref2) * 100.)

                #--------------
                #Orthogonal Regression in Retrieval
                #--------------
                xx = rPrf_ref2
                xx = xx.flatten()

                yy = rPrf_i2
                yy = yy.flatten()

                if pltInputs['errorFlg']:
                    xx_e = rPrf_ref2_e
                    xx_e = xx_e.flatten()
                
                    yy_e = rPrf_i2_e
                    yy_e = yy_e.flatten()
                else:
                    xx_e = xx*0.
                    yy_e = yy*0.


                odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx_e, yerr=yy_e, InError=True)
                slope.setdefault(k, []).append(float(odr[0]))
                intercept.setdefault(k, []).append(float(odr[1]))

                slope_e.setdefault(k, []).append(float(odrErr[0]))
                intercept_e.setdefault(k, []).append(float(odrErr[1]))

                slope2lr, intercept2lr, r_value2lr, p_value2lr, std_err2lr = stats.linregress(xx, yy)
                rvalue.setdefault(k, []).append(float(r_value2lr))
                #slope.setdefault(k, []).append(float(slope2lr))
                #intercept.setdefault(k, []).append(float(intercept2lr))

                print 'Slope         = {0:.2f} +/- {1:.2f}'.format(float(odr[0]), float(odrErr[0]))
                print 'Intercept     = {0:.3f} +/- {1:.3f}'.format(float(odr[1]), float(odrErr[1]))
                print 'r-value       = {0:.2f}'.format(float(r_value2lr))

    colors2 = 'ykrgbcm'
    colors3 = 'mkrgbcm'

    #---------------------------------
    #Figure: SINGLE selected profile
    #---------------------------------   
    if pltInputs['sdoi']:

        doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in pltInputs['sdoi']]
        doidt = np.asarray(doidt) 

        fig, ax0   = plt.subplots(1, 2, figsize=(8,6), sharey=True, sharex=False)
        fig2, ax02 = plt.subplots(figsize=(4.5,6), sharey=True, sharex=False)

        for vi, v in enumerate(pltInputs['ver']):

            indsH = np.where(fts.alt[v] < 16)[0]

            days = np.asarray([dt.date(dd.year,dd.month,dd.day) for dd in fts.dates[v]])

            for doi_i in doidt:

                deltadoi = days - doi_i

                indsTime   = np.where(deltadoi == dt.timedelta(0))[0]
                
                PrfMean    =  np.nanmean(fts.rPrfVMR[v][indsTime],axis=0)
                PrfStd     =  np.nanstd(fts.rPrfVMR[v][indsTime],axis=0)

                PrfMeanWP  =  np.nanmean(fts.waterVMR[v][indsTime],axis=0) *1e-6
                PrfStdWP   =  np.nanstd(fts.waterVMR[v][indsTime],axis=0)

                PrfRef        = np.nanmean(fts.rPrfVMR[pltInputs['ver'][0]][indsTime],axis=0)

                ax02.plot(PrfMeanWP[indsH],fts.alt[v][indsH], color=colors3[vi], linewidth=2.0, zorder=0)
                ax02.scatter(PrfMeanWP[indsH],fts.alt[v][indsH],facecolors='white', s=45, color=colors3[vi], zorder=1, label=idvers[vi])
                #ax0[0].fill_betweenx(fts.alt[v],PrfMeanWP-PrfStdWP,PrfMeanWP+PrfStdWP,alpha=0.25,color=clr[vi])

                ax0[0].plot(PrfMean[indsH],fts.alt[v][indsH], color=colors3[vi], linewidth=2.0,zorder=0)
                ax0[0].scatter(PrfMean[indsH],fts.alt[v][indsH],facecolors='white', s=45, color=colors3[vi], zorder=1, label=idvers[vi])
                #ax0[1].scatter(PrfMean,fts.alt[v],facecolors='white', s=35, color=clr[vi], zorder=1)
                #ax0[1].fill_betweenx(fts.alt[v],PrfMean-PrfStd,PrfMean+PrfStd,alpha=0.25,color=clr[vi])

        
                if v != pltInputs['ver'][0]: 
                    ax0[1].plot((PrfMean[indsH] - PrfRef[indsH])/PrfRef[indsH] * 100.,fts.alt[v][indsH], color=colors3[vi], linewidth=2.0,zorder=0)
                    ax0[1].scatter((PrfMean[indsH] - PrfRef[indsH])/PrfRef[indsH] * 100.,fts.alt[v][indsH], facecolors='white', s=45, color=colors3[vi], zorder=1, label=idvers[vi])
                #ax0[2].scatter(PrfMean,fts.alt[v],facecolors='white', s=35, color=clr[vi], zorder=1)

        ax02.set_xlabel('VMR [x10$^3$ ppm]', fontsize=16)
        ax02.set_ylabel('Altitude [km]', fontsize=16)
        ax02.set_ylim(1, 15)                   
        ax02.grid(True,which='both', alpha=0.35)
        ax02.tick_params(which='both',labelsize=16)
        ax02.set_title('Water Vapor', fontsize=16)   
        ax02.set_xlim(xmin=0, xmax=22)
        ax02.legend(prop={'size':12})

        if pltInputs['gasName'].lower() == 'hcn': 
            gaslabel = 'HCN'
        elif pltInputs['gasName'].lower() == 'c2h6': 
            gaslabel = 'C$_2$H$_6$'
        elif pltInputs['gasName'].lower() == 'co': 
            gaslabel = 'CO'
        else: 
            gaslabel = pltInputs['gasName'].upper()


        ax0[0].set_xlabel('VMR [ppb]', fontsize=16)                 
        ax0[0].grid(True,which='both', alpha=0.35)
        ax0[0].tick_params(which='both',labelsize=16)
        ax0[0].set_title(gaslabel, fontsize=16)  
        ax0[0].set_ylim(1, 15)
        #ax0[0].set_yticklabels([])
        ax0[0].set_ylabel('Altitude [km]', fontsize=16)
        #ax0[0].yaxis.label.set_visible(False)
        #ax0[1].set_xlim(xmin=0)

        ax0[1].set_xlabel('Relative Difference [%]', fontsize=16)                 
        ax0[1].grid(True,which='both', alpha=0.35)
        ax0[1].tick_params(which='both',labelsize=16)
        ax0[1].set_title(gaslabel + ' Difference', fontsize=16) 
        #ax0[1].set_xlim(xmin=0)
        #ax0[0].set_title(doi_i, fontsize=14)   
        #ax0[0].legend(prop={'size':11})

                #ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(len(indsTime)), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)
                
        fig.subplots_adjust(bottom=0.1,top=0.95, left=0.08, right=0.975)
        fig2.subplots_adjust(bottom=0.1,top=0.95, left=0.15, right=0.975)  

                        
        if pltInputs['saveFlg']: 
            pdfsav.savefig(fig,dpi=200)
            fig.savefig(pltInputs['pltDir']+pltInputs['gasName']+'_singlePrf_'+pltInputs['sdoi'][0]+'_'+pltInputs['loc']+'.pdf', bbox_inches='tight')
            fig2.savefig(pltInputs['pltDir']+pltInputs['gasName']+'_singlePrf_WV_'+pltInputs['sdoi'][0]+'_'+pltInputs['loc']+'.pdf', bbox_inches='tight')
        
        else:       
            plt.show(block=False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()           # Exit program     

    
    #-------------------
    #Slope and Intercept Bars
    #-------------------
    for pi, pcol in enumerate(pCols[0:-1]):

        if pi == 0:

            pcolstr = str(pcol[0])

            fig, (ax, ax2, ax3)  = plt.subplots(3, 1, figsize=(7.5,9), sharex=True)

            ind = np.arange(len(pCols[0:-1]))                        
            width = 1. / (len(vers))

            labels     = []
            bar_groups = []
            
            for k in range(len(vers)):

                if k != indT:        

                    bars = ax.bar(ind+k*width, slope[k], width, yerr=slope_e[k], ecolor='k', color=colors2[k % len(colors2)])
                    bar_groups.append(bars)  
                    labels.append(idvers[k])    

                    ax2.bar(ind+k*width, intercept[k], width, yerr=intercept_e[k], ecolor='k', color=colors2[k % len(colors2)])

                    ax3.bar(ind+k*width, rvalue[k], width, color=colors2[k % len(colors2)])

                         
            ax.set_ylabel('Slope', fontsize=16)
            ax.legend([b[0] for b in bar_groups], labels, fontsize=11, ncol=len(labels), loc='upper center', bbox_to_anchor=(0.5, 1.1))
            ax.tick_params(which='both',labelsize=16)
            ax.axhline(y=1.0, linestyle='--', linewidth=1.0, color='k', alpha=0.5)
            ax.set_ylim(0.85, 1.1)

            ax2.set_ylabel('Intercept [{}]'.format(pltInputs['sclfctName']), fontsize=16)
            ax2.tick_params(which='both',labelsize=16)
            ax2.axhline(y=0.0, linestyle='-', linewidth=1.0, color='k')

            ax3.set_xlim(-width*0.5,len(ind)+width)
            ax3.set_ylabel('r-value', fontsize=16)
            xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
            ax3.set_xticks(ind+width*3.5)
            xtickNames = ax3.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=0, fontsize=12)
            ax3.tick_params(which='both',labelsize=16)
            ax3.set_xlabel('Layer [km]', fontsize=16)
            ax3.set_ylim(0.9, 1.005)
            ax3.axhline(y=1.0, linestyle='--', linewidth=1.0, color='k', alpha=0.5)
            #ax3.text(0.01, 0.95, '(c)', va='center', ha='left', transform=ax3.transAxes,fontsize=14)

            fig.subplots_adjust(bottom=0.075,top=0.975, left=0.15, right=0.95)

            if pltInputs['saveFlg']: 
                pdfsav.savefig(fig,dpi=200)
                plt.savefig(pltInputs['pltDir']+pltInputs['gasName']+'_OLR_'+pltInputs['loc']+'.pdf', bbox_inches='tight')
    
            else:       
                plt.show(block=False)

    #-------------------
    #Bias and Precision Bars
    #-------------------
    
    for pi, pcol in enumerate(pCols[0:-1]):

        if pi == 0:

            pcolstr = str(pcol[0])

            fig, (ax, ax2, ax3, ax4)  = plt.subplots(4, 1, figsize=(7.5,9), sharex=True)

            ind = np.arange(len(pCols[0:-1]))                        
            width = 1. / (len(vers))

            labels     = []
            bar_groups = []
            
            for k in range(len(vers)):

                if k != indT:        

                    bars = ax.bar(ind+k*width, bias[k], width, yerr=bias_e[k], ecolor='k', color=colors2[k % len(colors2)])
                    bar_groups.append(bars)  
                    labels.append(idvers[k])    

                    ax2.bar(ind+k*width, bias_perc[k], width, yerr=bias_perc_e[k], ecolor='k', color=colors2[k % len(colors2)])

                    ax3.bar(ind+k*width, prec[k], width, color=colors2[k % len(colors2)])
                    #ax3.bar(ind+k*width, prec[k], width, yerr=prec_e[k], ecolor='k', color=colors2[k % len(colors2)])

                    ax4.bar(ind+k*width, prec_perc[k], width, color=colors2[k % len(colors2)])
                    #ax4.bar(ind+k*width, prec_perc[k], width, yerr=prec_perc_e[k], ecolor='k', color=colors2[k % len(colors2)])                    

            ax.set_ylabel('Bias [{}]'.format(pltInputs['sclfctName']), fontsize=16)
            #ax.set_ylim(0,1.5)
            #ax.axhline(y=1.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
            ax.legend([b[0] for b in bar_groups], labels, fontsize=11, ncol=len(labels), loc='upper center', bbox_to_anchor=(0.5, 1.1))
            ax.tick_params(which='both',labelsize=16)
            ax.axhline(y=0.0, linestyle='-', linewidth=1.0, color='k')

            ax2.set_ylabel('Bias [%]', fontsize=16)
            ax2.tick_params(which='both',labelsize=16)
            #ax2.set_ylim(-15, 30)
            ax2.axhline(y=0.0, linestyle='-', linewidth=1.0, color='k')

            ax3.set_ylabel('Precision [{}]'.format(pltInputs['sclfctName']), fontsize=16)
            ax3.tick_params(which='both',labelsize=16)
            ax3.axhline(y=0.0, linestyle='--', linewidth=1.0, color='k', alpha=0.5)

            ax4.set_ylabel('Precision [%]', fontsize=16)
            xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
            ax4.set_xticks(ind+width*3.5)
            ax4.set_xlim(-width*0.5,len(ind)+width)
            xtickNames = ax4.set_xticklabels(xTickMarks)
            plt.setp(xtickNames, rotation=0, fontsize=12)
            ax4.tick_params(which='both',labelsize=16)
            ax4.set_xlabel('Layer [km]', fontsize=16)
            #ax3.text(0.01, 0.95, '(c)', va='center', ha='left', transform=ax3.transAxes,fontsize=14)

            fig.subplots_adjust(bottom=0.075,top=0.975, left=0.15, right=0.95)


            if pltInputs['saveFlg']: 
                pdfsav.savefig(fig,dpi=200)
                plt.savefig(pltInputs['pltDir']+pltInputs['gasName']+'_Bias_'+pltInputs['loc']+'.pdf', bbox_inches='tight')
    
            else:       
                plt.show(block=False)



    #---------------------------------
    #Figure: selected profiles (raw grid)
    #---------------------------------   
    if pltInputs['dois']:

        doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in pltInputs['dois']]
        doidt = np.asarray(doidt) 

        fig, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)

        for vi, v in enumerate(pltInputs['ver']):

            ndoi = 0

            days = np.asarray([dt.date(dd.year,dd.month,dd.day) for dd in fts.dates[v]])

            for doi_i in doidt:

                deltadoi = days - doi_i

                indsTime = np.where(deltadoi == dt.timedelta(0))[0]
                PrfMean  =  np.nanmean(fts.rPrfVMR[v][indsTime],axis=0)
                PrfStd   =  np.nanstd(fts.rPrfVMR[v][indsTime],axis=0)

                if ndoi < len(doidt):

                    if ndoi<=4:

                        if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                        ax0[0, ndoi].plot(PrfMean,fts.alt[v], color=clr[vi], linewidth=2.0, label=idvers[vi], zorder=5)
                        ax0[0, ndoi].scatter(PrfMean,fts.alt[v],facecolors='white', s=35, color=clr[vi], zorder=6)
                        ax0[0, ndoi].fill_betweenx(fts.alt[v],PrfMean-PrfStd,PrfMean+PrfStd,alpha=0.25,color=clr[vi])

                        ax0[0, ndoi].set_ylim(1, 15)
                        if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                        if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                           
                        ax0[0, ndoi].grid(True,which='both')
                        ax0[0, ndoi].tick_params(which='both',labelsize=14)
                        ax0[0, ndoi].set_title(doi_i, fontsize=14)   

                        #ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(len(indsTime)), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)

                        ax0[0, ndoi].set_xlim(xmin=0)
                        if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                    if ndoi>=5:

                        ax0[1, (ndoi-5)].plot(PrfMean,fts.alt[v], color=clr[vi],  linewidth=2.0, label=idvers[vi], zorder=5)
                        ax0[1, (ndoi-5)].scatter(PrfMean,fts.alt[v],facecolors='white', s=35, color=clr[vi], zorder=6)
                        ax0[1, (ndoi-5)].fill_betweenx(fts.alt[v],PrfMean-PrfStd,PrfMean+PrfStd,alpha=0.25,color=clr[vi])

                        if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)

                        if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                        if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)
    
                        if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [{}]'.format(pltInputs['sclfctName']), fontsize=14)
                        ax0[1, (ndoi-5)].grid(True,which='both')
                        ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                        ax0[1, (ndoi-5)].set_title(doi_i, fontsize=14)

                        #ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(len(indsTime)), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)

                        ax0[1, (ndoi-5)].set_xlim(xmin=0)

                    ndoi+=1

        fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 

        #idver = getidver(v)
                        
        if pltInputs['saveFlg']: 
            pdfsav.savefig(fig,dpi=200)
        
        else:       
            plt.show(block=False)
    
     
    if pltInputs['saveFlg']: pdfsav.close()

    print('\nFinished Plots.......\n')
       
 #    #--------------------------------
 #    # Pause so user can look at plots
 #    #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program       

if __name__ == "__main__":
    main(sys.argv[1:])