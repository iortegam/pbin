import datetime as dt
import time
import math
import sys
import numpy as np
import os
import csv
import itertools
from collections import OrderedDict
import os
from os import listdir
from os.path import isfile, join
import re
from scipy.integrate import simps
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import matplotlib.dates as md

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import glob
import myfunctions as mf

import dataOutClass as dc
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
import pylab as P

from math import acos, asin, atan2, cos, hypot
from math import degrees, pi as PI, radians, sin, tan
from itertools import izip




#------------------------------------------------------------------------------------------------------------------------------    
class _DateRange(object):
    '''
    This is an extension of the datetime module.
    Adds functionality to create a list of days.
    '''
    def __init__(self,iyear,imnth,iday,fyear,fmnth,fday, incr=1):
        self.i_date   = dt.date(iyear,imnth,iday)                                                     # Initial Day
        self.f_date   = dt.date(fyear,fmnth,fday)                                                     # Final Day
        self.dateList =[self.i_date + dt.timedelta(days=i) for i in range(0, self.numDays(), incr)]   # Incremental day list from initial to final day

    def numDays(self):
        '''Counts the number of days between start date and end date'''
        return (self.f_date + dt.timedelta(days=1) - self.i_date).days

    def inRange(self,crntyear,crntmonth,crntday):
        '''Determines if a specified date is within the date ranged initialized'''
        crntdate = dt.date(crntyear,crntmonth,crntday)
        if self.i_date <= crntdate <= self.f_date:
            return True
        else:
            return False

    def nearestDate(self, year, month, day=1, daysList=False):
        ''' Finds the nearest date from a list of days based on a given year, month, and day'''
        testDate = dt.date(year, month, day)
        if not daysList:
            daysList = self.dateList
        return min( daysList, key=lambda x:abs(x-testDate) )

    def yearList(self):
        ''' Gives a list of unique years within DateRange '''
        years = [ singDate.year for singDate in self.dateList]               # Find years for all date entries
        years = list(set(years))                                             # Determine all unique years
        years.sort()
        return years

    def daysInYear(self,year):
        ''' Returns an ordered list of days from DateRange within a specified year '''
        if isinstance(year,int):
            newyears = [inYear for inYear in self.dateList if inYear.year == year]
            return newyears
        else:
            print 'Error!! Year must be type int for daysInYear'
            return False

class FTSClass():

    def __init__(self, gasName, retDir, ctlF, ver, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,outFname='', fleoutFlg=False, version=''):
    

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

        if not( retDir.endswith('/') ): retDir = retDir + '/'

        self.ver      = ver
        self.gasName  = gasName
        self.retDir   = [retDir+g.lower()+'/'+v+'/' for g,v in izip(gasName,ver)] 
        self.ctlFile  = [retDir+g.lower()+'/'+'x.'+g.lower()+'/'+ c for g,c in izip(gasName,ctlF)] 
        
        #---------------------------
        # Check file and directories
        #---------------------------
        for d in self.retDir:  dc.ckDir(d,exitFlg=True)
        for c in self.ctlFile: dc.ckFile(c,exitFlg=True)

        self.iyear = iyear
        self.imnth = imnth
        self.iday  = iday

        self.fyear = fyear
        self.fmnth = fmnth
        self.fday  = fday


    def ReadFTS(self, fltrFlg=False,minSZA=0.0,maxSZA=90.0,maxRMS=10.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,dofFlg=False,rmsFlg=True,tcFlg=True,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
               pcFlg=True,cnvrgFlg=True, allGas=False,sclfct=1.0,sclname='ppv',pltStats=True,szaFlg=False,errFlg=False,chiFlg=False,tcMMflg=False,mnthFltFlg=False, pColsFlg=False, pCols=[[0,8], [8,16]]):

        #-------------------------------------
        # Create instance of output data class   
        #-------------------------------------
        statDataCl = OrderedDict()
        #for i, v in enumerate(self.ver):
        #    statDataCl[v] = dc.ReadOutputData(self.retDir[i],'',self.ctlFile[i],self.iyear,self.imnth,self.iday,self.fyear,self.fmnth,self.fday)

        for i, (g,v) in enumerate(izip(self.gasName,self.ver)):
            statDataCl[g+'_'+v] = dc.ReadOutputData(self.retDir[i],'',self.ctlFile[i], self.iyear,self.imnth,self.iday,self.fyear,self.fmnth,self.fday)

        #--------------
        # Read profiles
        #--------------
        for GasVer in statDataCl:
            statDataCl[GasVer].readprfs([statDataCl[GasVer].PrimaryGas],retapFlg=1)
            statDataCl[GasVer].readprfs([statDataCl[GasVer].PrimaryGas],retapFlg=0)
        
            if statDataCl[GasVer].empty: 
                print 'No retreivals found for {}. Exiting.......'.format(GasVer)
                sys.exit()

            #statLyrFile  = '/data/Campaign/FL0/local/station.layers'
            #statDataCl[GasVer].readStatLyrs(statLyrFile)  

        self.rPrfVMR          = OrderedDict()
        self.rPrfMol          = OrderedDict()
        self.dates            = OrderedDict()
        self.alt              = OrderedDict()
        self.Airmass          = OrderedDict()
        self.waterVMR         = OrderedDict()
        self.waterMol         = OrderedDict()
        self.totClmn          = OrderedDict()
        self.TCdryAir         = OrderedDict()
        self.TCdry            = OrderedDict()
        self.rPrfDry          = OrderedDict()
        self.rms              = OrderedDict()
        self.dofs             = OrderedDict()
        self.dofsAvg          = OrderedDict()
        self.dofsAvg_cs       = OrderedDict()
        self.LayVMR           = OrderedDict()
        self.DryLayVMR        = OrderedDict()
        self.upperHgt         = OrderedDict()
        self.lowerHgt         = OrderedDict()
        self.LayThk           = OrderedDict()
        self.LayDOF           = OrderedDict()

        self.avkSCF           = OrderedDict()
        self.avkVMR           = OrderedDict()
        self.avkSCFav         = OrderedDict()
        self.avkVMRav         = OrderedDict()

        self.aPrfVMR          = OrderedDict()
        self.aPrfMol          = OrderedDict()

        self.sza              = OrderedDict()
        self.saa              = OrderedDict()

        if errFlg:          
            self.rand_err         = OrderedDict()
            self.sys_err          = OrderedDict()
            self.tot_err          = OrderedDict()

            self.rand_errvmr      = OrderedDict()
            self.sys_errvmr       = OrderedDict()
            self.tot_errvmr       = OrderedDict()

            self.rand_cmpnts      = OrderedDict() 
            self.sys_cmpnts       = OrderedDict()

            self.rand_cmpnts_vmr  = OrderedDict()
            self.sys_cmpnts_vmr   = OrderedDict()

            self.tot_rnd          = OrderedDict()
            self.tot_sys          = OrderedDict()
            self.tot_std          = OrderedDict()

            self.err_summary      = OrderedDict()

        if pColsFlg:
            self.vmrP             = OrderedDict()
            self.TCp              = OrderedDict()

            if errFlg:
                self.vmrP_Toterr     = OrderedDict()
                self.TCp_Toterr      = OrderedDict()

                self.vmrP_Ranerr     = OrderedDict()
                self.TCp_Ranerr      = OrderedDict()

                self.vmrP_Syserr     = OrderedDict()
                self.TCp_Syserr      = OrderedDict()

        for j, GasVer in enumerate(statDataCl):
 
            self.rPrfVMR[GasVer]  = np.asarray(statDataCl[GasVer].rprfs[statDataCl[GasVer].PrimaryGas]) * sclfct
            self.rPrfMol[GasVer]  = np.asarray(statDataCl[GasVer].rprfs[statDataCl[GasVer].PrimaryGas]  * np.asarray(statDataCl[GasVer].rprfs['AIRMASS']))
            self.dates[GasVer]    = np.asarray(statDataCl[GasVer].rprfs['date'])
            self.alt[GasVer]      = np.asarray(statDataCl[GasVer].rprfs['Z'][0,:])
           
            #self.LayThk[GasVer]   = np.asarray(statDataCl[GasVer].thckAlt[:-1]) 
            #print statDataCl[GasVer].thckAlt
            #print statDataCl[GasVer].midAlt
            

            self.Airmass[GasVer]  = np.asarray(statDataCl[GasVer].rprfs['AIRMASS'])
            self.waterVMR[GasVer] = np.asarray(statDataCl[GasVer].aprfs['H2O']) * sclfct
            self.waterMol[GasVer] = np.asarray(statDataCl[GasVer].aprfs['H2O'] * self.Airmass[GasVer] )
            self.totClmn[GasVer]  = np.sum(self.rPrfMol[GasVer],axis=1)
            #TCdryAir[GasVer] = np.sum(Airmass[GasVer],axis=1) - np.sum(waterMol[ver],axis=1)
            #TCdry[GasVer]    = (totClmn[GasVer] / TCdryAir[GasVer]) * sclfct

            self.aPrfVMR[GasVer]  = np.asarray(statDataCl[GasVer].aprfs[statDataCl[GasVer].PrimaryGas]) * sclfct
            self.aPrfMol[GasVer]  = np.asarray(statDataCl[GasVer].aprfs[statDataCl[GasVer].PrimaryGas]  * np.asarray(statDataCl[GasVer].aprfs['AIRMASS']))

            #----------------------------------------
            # This is the mixing ratio for DRY AIR!!!
            #----------------------------------------
            #rPrfDry[ver] = np.asarray(statDataCl[GasVer].rprfs[statDataCl[GasVer].PrimaryGas]) / (1.0 - waterVMR[ver]) * sclfct       
            
            #----------------------------------
            # Read Summary data (For filtering)
            #----------------------------------
            statDataCl[GasVer].readsummary()
            self.rms[GasVer]     = np.asarray(statDataCl[GasVer].summary[statDataCl[GasVer].PrimaryGas+'_FITRMS'])

            #----------------------------------
            # Read readPbp (observed, fitted, and difference spectra)
            #----------------------------------
            mw    = [str(int(x)) for x in statDataCl[GasVer].ctl['band']]     
            numMW = len(mw)
            statDataCl[GasVer].readPbp()

            self.sza[GasVer]   = statDataCl[GasVer].pbp['sza']
            self.saa[GasVer]   = statDataCl[GasVer].pbp['saa']

            #----------------------------------
            # Read Spectra for each gas
            #----------------------------------
            statDataCl[GasVer].readSpectra(statDataCl[GasVer].gasList)

            #-------------------- 
            # Call to filter data
            #--------------------
            if fltrFlg: statDataCl[GasVer].fltrData(statDataCl[GasVer].PrimaryGas,mxrms=maxRMS[j],rmsFlg=rmsFlg, minDOF=minDOF[j],  dofFlg=dofFlg, tcFlg=tcFlg, pcFlg=pcFlg , cnvrgFlg=cnvrgFlg)
            else:       statDataCl[GasVer].inds = np.array([]) 

                #--------------------------------------------
            # Read Error data to get AVK and profile DOFs
            #-------------------------------------------------
            # Determine if AVK has been created via sfit4 core
            # code or via python error analysis
            #-------------------------------------------------     
            if errFlg:   # Read AVK from error output
                #statDataCl[GasVer].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
                statDataCl[GasVer].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False)

                npnts           = np.shape(statDataCl[GasVer].error['Total_Random_Error'])[0]
                nlvls           = np.shape(self.alt[GasVer])[0]
                
                #-------------------------------------------------
                #Error profiles constructed with the diagonal elements of the covariances matrices
                #-------------------------------------------------
                self.rand_err[GasVer]   = np.zeros((npnts,nlvls))
                self.sys_err[GasVer]    = np.zeros((npnts,nlvls))

                for i in range(npnts):
                    self.rand_err[GasVer][i,:] = np.diag(statDataCl[GasVer].error['Total_Random_Error'][i][:,:])
                    self.sys_err[GasVer][i,:]  = np.diag(statDataCl[GasVer].error['Total_Systematic_Error'][i][:,:])

                self.tot_err[GasVer]     = np.sqrt(self.rand_err[GasVer] + self.sys_err[GasVer])            
                self.sys_err[GasVer]     = np.sqrt(self.sys_err[GasVer])
                self.rand_err[GasVer]    = np.sqrt(self.rand_err[GasVer])

                self.sys_errvmr[GasVer]  = self.sys_err[GasVer]/ np.asarray(statDataCl[GasVer].rprfs['AIRMASS']) * sclfct
                self.rand_errvmr[GasVer] = self.rand_err[GasVer]/ np.asarray(statDataCl[GasVer].rprfs['AIRMASS']) * sclfct
                self.tot_errvmr[GasVer]  = np.sqrt(self.rand_errvmr[GasVer]**2 + self.sys_errvmr[GasVer]**2) 

                #-------------------------------------------------
                #Error profiles of components
                #-------------------------------------------------
                self.rand_cmpnts[GasVer] = statDataCl[GasVer].randErrDiag
                self.sys_cmpnts[GasVer]  = statDataCl[GasVer].sysErrDiag

                self.rand_cmpnts_vmr[GasVer] = statDataCl[GasVer].randErrDiag 
                self.sys_cmpnts_vmr[GasVer]  = statDataCl[GasVer].sysErrDiag

                for k in self.sys_cmpnts_vmr[GasVer]:
                    self.sys_cmpnts_vmr[GasVer][k] = (np.sqrt(self.sys_cmpnts_vmr[GasVer][k])/ np.asarray(statDataCl[GasVer].rprfs['AIRMASS']))*sclfct
                  
                for k in self.rand_cmpnts_vmr[GasVer]:
                    self.rand_cmpnts_vmr[GasVer][k] = (np.sqrt(self.rand_cmpnts_vmr[GasVer][k])/ np.asarray(statDataCl[GasVer].rprfs['AIRMASS']))*sclfct

                #-------------------------------------------------
                #Error in the summary output. Get Total Errors 
                #-------------------------------------------------
                self.tot_rnd[GasVer]        = np.array(statDataCl[GasVer].error['Total random uncertainty'])
                self.tot_sys[GasVer]        = np.array(statDataCl[GasVer].error['Total systematic uncertainty'])
                self.tot_std[GasVer]        = np.sqrt(self.tot_rnd[GasVer]**2 + self.tot_sys[GasVer]**2)
                
                self.err_summary[GasVer]    = statDataCl[GasVer].error

                #---------------------
                # Get averaging kernel
                #---------------------   
                self.avkSCF[GasVer]    = np.delete(np.asarray(statDataCl[GasVer].error['AVK_scale_factor']),statDataCl[GasVer].inds,axis=0)
                self.avkVMR[GasVer]    = np.delete(np.asarray(statDataCl[GasVer].error['AVK_vmr']),statDataCl[GasVer].inds,axis=0)   
                diagAK                 = np.diagonal(self.avkSCF[GasVer],axis1=1,axis2=2)
                self.dofs[GasVer]      = np.sum(diagAK, axis=1)        
                self.avkSCFav[GasVer]  = np.mean(self.avkSCF[GasVer],axis=0)    
                self.avkVMRav[GasVer]  = np.mean(self.avkVMR[GasVer],axis=0)                   
                    
                self.dofsAvg[GasVer]    = np.diag(self.avkSCFav[GasVer])
                self.dofsAvg_cs[GasVer] = np.cumsum(np.diag(self.avkSCFav[GasVer])[::-1])[::-1]

         
            
            else:        # Read AVK from sfit4 output (only contains scaled AVK)
                avkSCFi = []
                for d in statDataCl[GasVer].dirLst:
                    lines  = dc.tryopen( d + statDataCl[GasVer].ctl['file.out.ak_matrix'][0])
                    if not lines: continue
                    avkSCFi.append(np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[2:] ] ))
                    
                if not statDataCl[GasVer].readPrfFlgApr[statDataCl[GasVer].PrimaryGas]: statDataCl[GasVer].readprfs([statDataCl[GasVer].PrimaryGas],retapFlg=0)   # Apriori Profiles
                          
                self.avkSCF[GasVer]  = np.asarray(avkSCFi)
                nobs            = np.shape(self.avkSCF[GasVer])[0]
                n_layer         = np.shape(self.avkSCF[GasVer])[1]
                self.avkVMR[GasVer]  = np.zeros((nobs,n_layer,n_layer))
        
                for obs in range(0,nobs):
                    Iapriori        = np.zeros((n_layer,n_layer))
                    IaprioriInv     = np.zeros((n_layer,n_layer))
                    np.fill_diagonal(Iapriori,statDataCl[GasVer].aprfs[statDataCl[GasVer].PrimaryGas.upper()][obs])
                    np.fill_diagonal(IaprioriInv, 1.0 / (statDataCl[GasVer].aprfs[statDataCl[GasVer].PrimaryGas.upper()][obs]))
                    self.avkVMR[GasVer][obs,:,:] = np.dot(np.dot(Iapriori,np.squeeze(self.avkSCF[GasVer][obs,:,:])),IaprioriInv)       
                    
                self.avkSCF[GasVer]     = np.delete(self.avkSCF[GasVer],statDataCl[GasVer].inds,axis=0)
                self.avkVMR[GasVer]     = np.delete(self.avkVMR[GasVer],statDataCl[GasVer].inds,axis=0)
                self.dofs[GasVer]       = np.diagonal(self.avkSCF[GasVer],axis1=1,axis2=2)
                self.avkSCFav[GasVer]   = np.mean(self.avkSCF[GasVer],axis=0)
                self.avkVMRav[GasVer]   = np.mean(self.avkVMR[GasVer],axis=0)
                
                self.dofsAvg[GasVer]    = np.diag(self.avkSCFav[GasVer])
                self.dofsAvg_cs[GasVer] = np.cumsum(np.diag(self.avkSCFav[GasVer])[::-1])[::-1]


            #--------------------------------------
            # Remove retrieval data based on filter
            #--------------------------------------
            nfltr              = len(statDataCl[GasVer].inds)
            self.rms[GasVer]      = np.delete(self.rms[GasVer],statDataCl[GasVer].inds)
            self.sza[GasVer]      = np.delete(self.sza[GasVer],statDataCl[GasVer].inds)
            self.saa[GasVer]      = np.delete(self.saa[GasVer],statDataCl[GasVer].inds)
            ntot               = len(self.rms[GasVer])
            self.dates[GasVer]    = np.delete(self.dates[GasVer],statDataCl[GasVer].inds)
            self.totClmn[GasVer]  = np.delete(self.totClmn[GasVer],statDataCl[GasVer].inds)
            self.rPrfVMR[GasVer]  = np.delete(self.rPrfVMR[GasVer],statDataCl[GasVer].inds,axis=0)   
            self.rPrfMol[GasVer]  = np.delete(self.rPrfMol[GasVer],statDataCl[GasVer].inds,axis=0)
            self.waterVMR[GasVer] = np.delete(self.waterVMR[GasVer],statDataCl[GasVer].inds,axis=0)   
            self.waterMol[GasVer] = np.delete(self.waterMol[GasVer],statDataCl[GasVer].inds,axis=0) 

    

            #self.rPrfConc[GasVer]  = np.delete(self.rPrfMol[GasVer],statDataCl[GasVer].inds,axis=0)  

            #rPrfDry[GasVer] = np.delete(rPrfDry[GasVer],statDataCl[GasVer].inds,axis=0)   
            self.Airmass[GasVer] = np.delete(self.Airmass[GasVer],statDataCl[GasVer].inds,axis=0)
            #TCdry[GasVer]   = np.delete(TCdry[GasVer],statDataCl[GasVer].inds)

            self.aPrfVMR[GasVer] = np.delete(self.aPrfVMR[GasVer],statDataCl[GasVer].inds,axis=0)   
            self.aPrfMol[GasVer] = np.delete(self.aPrfMol[GasVer],statDataCl[GasVer].inds,axis=0)

            if errFlg:
                self.rand_err[GasVer] = np.delete(self.rand_err[GasVer],statDataCl[GasVer].inds,axis=0)
                self.sys_err[GasVer]  = np.delete(self.sys_err[GasVer],statDataCl[GasVer].inds,axis=0)  
                self.tot_err[GasVer]  = np.delete(self.tot_err[GasVer],statDataCl[GasVer].inds,axis=0)

                self.rand_errvmr[GasVer] = np.delete(self.rand_errvmr[GasVer],statDataCl[GasVer].inds,axis=0)
                self.sys_errvmr[GasVer]  = np.delete(self.sys_errvmr[GasVer],statDataCl[GasVer].inds,axis=0)  
                self.tot_errvmr[GasVer]  = np.delete(self.tot_errvmr[GasVer],statDataCl[GasVer].inds,axis=0)

                self.tot_rnd[GasVer] = np.delete(self.tot_rnd[GasVer],statDataCl[GasVer].inds)
                self.tot_sys[GasVer] = np.delete(self.tot_sys[GasVer],statDataCl[GasVer].inds)
                self.tot_std[GasVer] = np.delete(self.tot_std[GasVer],statDataCl[GasVer].inds)

                
                for k in self.sys_cmpnts[GasVer]:
                    self.sys_cmpnts[GasVer][k]      = np.delete(self.sys_cmpnts[GasVer][k],statDataCl[GasVer].inds,axis=0)
                    self.sys_cmpnts_vmr[GasVer][k]  = np.delete(self.sys_cmpnts_vmr[GasVer][k],statDataCl[GasVer].inds,axis=0)
                    
                for k in self.rand_cmpnts[GasVer]:
                    self.rand_cmpnts[GasVer][k]     = np.delete(self.rand_cmpnts[GasVer][k],statDataCl[GasVer].inds,axis=0)
                    self.rand_cmpnts_vmr[GasVer][k] = np.delete(self.rand_cmpnts_vmr[GasVer][k],statDataCl[GasVer].inds,axis=0)

                for k in self.err_summary[GasVer]:
                    self.err_summary[GasVer][k]     = np.delete(self.err_summary[GasVer][k],statDataCl[GasVer].inds, axis=0)

            #-------------------------------------------------
            # Calculate partial columns and weighted VMR
            #-------------------------------------------------
            if pColsFlg:
                #ind1                  = mf.nearestind(pCols[0], self.alt[GasVer])
                #ind2                  = mf.nearestind(pCols[1], self.alt[GasVer])

                inds = np.where( (self.alt[GasVer] >= pCols[0]) & (self.alt[GasVer] <=pCols[1])  )[0]


                self.vmrP[GasVer]     = np.average(self.rPrfVMR[GasVer][:,inds],axis=1,weights=self.Airmass[GasVer][:,inds])           
                self.TCp[GasVer]      = np.sum(self.rPrfMol[GasVer][:,inds],axis=1)

                if errFlg:
                    self.vmrP_Syserr[GasVer]   = np.average(self.sys_errvmr[GasVer][:,inds],axis=1,weights=self.Airmass[GasVer][:,inds])
                    self.vmrP_Ranerr[GasVer]   = np.average(self.rand_errvmr[GasVer][:,inds],axis=1,weights=self.Airmass[GasVer][:,inds])  
                    self.vmrP_Toterr[GasVer]   = np.average(self.tot_errvmr[GasVer][:,inds],axis=1,weights=self.Airmass[GasVer][:,inds])           
                    
                    self.TCp_Syserr[GasVer]      = np.sum(self.sys_err[GasVer][:,inds],axis=1)
                    self.TCp_Ranerr[GasVer]      = np.sum(self.rand_err[GasVer][:,inds],axis=1)
                    self.TCp_Toterr[GasVer]      = np.sum(self.tot_err[GasVer][:,inds],axis=1)



