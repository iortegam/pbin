#----------------------------------------------------------------------------------------
# Name:
#        HIPPO_ACE_outClass.py
#
# Purpose:

# Notes:
#
#
# Version History:
#       Created, Sep, 2016  Ivan Ortega (iortega@ucar.edu)
#
#
#----------------------------------------------------------------------------------------

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
from scipy import interpolate

import matplotlib
#
from numpy import fromstring, vectorize, ndarray, array, genfromtxt

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

import myfunctions as mf
from netCDF4 import Dataset, Group

from mpl_toolkits.basemap import Basemap
import sfitClasses as sc
import glob
from datetime import datetime, timedelta
from geographiclib.geodesic import Geodesic
from mpl_toolkits.basemap import Basemap
import srtm
import dataOutClass as dc




def ckDir(dirName,logFlg=False,exitFlg=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if logFlg: logFlg.error('Directory %s does not exist' % dirName)
        if exitFlg: sys.exit()
        return False
    else:
        return True   

def ckFile(fName,logFlg=False,exitFlg=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if logFlg: logFlg.error('Unable to find file: %s' % fName)
        if exitFlg: sys.exit()
        return False
    else:
        return True

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

class HIPPOClass():
    '''
    This class deals with reading HIPPO outputs.
    '''

    def __init__(self,dataDir, file,  outFname='', saveFlg= False):
        
        # check if '/' is included at end of path
        if not( dataDir.endswith('/') ):
            dataDir = dataDir + '/'           
            
        # Check if directory exits
        ckDir(dataDir, exitFlg=True)
        
        self.dirLst  = dataDir + file

        if saveFlg:
            self.pdfsav = PdfPages(outFname)
            self.saveFlg   = True
        else: self.saveFlg = False

        #---------------
        # ReadOutputData
        #---------------
    def ReadOutputHIPPO(self, gasname=''):
        '''Function to Read out the FLX-CEHM net CDF files.'''

        #infile = file(self.dirLst, 'r')
        cols, indextoname = mf.getColumns(self.dirLst, headerrow=0, delim=' ', header=True)
        #infile.close()

        #---------------
        # ReadOutputData
        #---------------
        instr        = np.asarray(cols[indextoname[0]])
        jd           = np.asarray(cols[indextoname[1]], dtype=float)
        year         = np.asarray(cols[indextoname[2]], dtype=int)
        Hno          = np.asarray(cols[indextoname[3]], dtype=int)
        flt          = np.asarray(cols[indextoname[4]], dtype=int)
        doy          = np.asarray(cols[indextoname[5]], dtype=float)
        utc          = np.asarray(cols[indextoname[6]], dtype=float)
        ut_mid       = np.asarray(cols[indextoname[7]])
        ut_start     = np.asarray(cols[indextoname[8]])
        ut_stop      = np.asarray(cols[indextoname[9]])

        #---------------
        # Aircraft information
        #---------------
        alt          = np.asarray(cols[indextoname[92]], dtype=float)
        lat          = np.asarray(cols[indextoname[93]], dtype=float)
        lon          = np.asarray(cols[indextoname[94]], dtype=float)

        #---------------
        # Gases
        #---------------
        
        if gasname.lower() == 'ocs': ncol = 29
        elif gasname.lower() == 'c2h6': ncol = 51
        elif gasname.lower() == 'c2h2': ncol = 52
        else: 
            print 'Gas {} is not setup in ReadOutputHIPPO. Go to ReadOutputHIPPO.py and modify it accordingly'.format(gasname)
            exit()

        gas = cols[indextoname[ncol]]
        gas = np.asarray(gas)

        #---------------
        # Filtering NAN
        #---------------
        inds    = np.where(gas == 'NA')[0]

        self.gasname   = gasname

        self.gas       = np.delete(gas, inds)
        self.instr     = np.delete(instr, inds)
        self.jd        = np.delete(jd, inds)
        self.year      = np.delete(year, inds)
        self.Hno       = np.delete(Hno, inds)
        self.flt       = np.delete(flt, inds)
        self.doy       = np.delete(doy, inds)
        self.utc       = np.delete(utc, inds)
        self.ut_mid    = np.delete(ut_mid, inds)
        self.ut_start  = np.delete(ut_start, inds)
        self.ut_stop   = np.delete(ut_stop, inds)
        self.alt       = np.delete(alt, inds)
        self.lat       = np.delete(lat, inds)
        self.lon       = np.delete(lon, inds)

        self.date      = [dt.datetime(self.year[i], 1, 1) + dt.timedelta(self.doy[i] - 1) + dt.timedelta(seconds=self.utc[i])  for i, da in enumerate(self.doy)]
        self.date      = np.asarray(self.date)

        inds    = np.where(lat == 'NA')[0]

        self.alt_all      = np.delete(alt, inds)
        self.lat_all      = np.delete(lat, inds)
        self.lon_all      = np.delete(lon, inds)
        self.Hno_all      = np.delete(Hno, inds)


    def pltPrf(self):

        HiPPO_Number = list(set(self.Hno))
        HiPPO_Number = np.asarray(HiPPO_Number, dtype=int)
        HiPPO_Number.sort()

        DateFmt      = DateFormatter('%Y-%m-%d')
        dayLc        = DayLocator()
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)

        layers     = np.arange(0, 18, 1.0)
        layers_mid = np.asarray([(layers[l] + layers[l+1])*0.5  for l in np.arange(0, len(layers)-1, 1)])
        
        for hn in HiPPO_Number:

            inds1     = np.where(self.Hno == hn)[0]
            
            lat       = self.lat[inds1]
            lon       = self.lon[inds1]
            alt       = self.alt[inds1]
            doy       = self.doy[inds1]
            dt        = self.date[inds1]
            RFHN      = self.flt[inds1]
            gas       = self.gas[inds1]

            RFs = list(set(RFHN))
            RFs = np.asarray(RFs, dtype=int)
            RFs.sort()

            for rf in RFs:
                inds2     = np.where(RFHN == rf)[0]

                gas_i     = np.asarray(gas[inds2], dtype=float)
                alt_i     = np.asarray(alt[inds2])/1000.

                gas_mean_i  = []
                gas_std_i   = []

                for l, la in enumerate(layers_mid):
                    inds3     = np.where( (alt_i >= layers[l]) & (alt_i < layers[l+1]))[0]
                    if len(inds3) >=2 : 
                        gas_mean_i.append(np.mean(gas_i[inds3]))
                        gas_std_i.append(np.std(gas_i[inds3]))
                    elif len(inds3) ==1 : 
                        gas_mean_i.append(gas_i[inds3])
                        gas_std_i.append(float('nan'))
                    else: 
                        gas_mean_i.append(float('nan'))
                        gas_std_i.append(float('nan'))

                gas_mean_i = np.asarray(gas_mean_i)
                gas_std_i  = np.asarray(gas_std_i)

                fig, ax = plt.subplots(1)

                ax.plot(gas_i, alt_i, color='k', linestyle='none', marker ='.', markersize=4)
                ax.scatter(gas_mean_i, layers_mid, facecolors='red', edgecolors='black', s=35)
                #ax.fill_betweenx(layers_mid,gas_mean_i-gas_std_i,gas_mean_i+gas_std_i,alpha=0.5,color='0.75')
                ax.errorbar(gas_mean_i, layers_mid,fmt=None,xerr=gas_std_i,ecolor='r')

                ax.tick_params(labelsize=16)
                
                ax.grid(True)
                ax.set_xlabel('Gas', fontsize = 16)
                ax.set_ylabel('Altitude [km]', fontsize = 16)
            
                fig.subplots_adjust(left = 0.125, bottom=0.15, top=0.95, right = 0.95)
                
                if self.saveFlg: 
                   self.pdfsav.savefig(fig,dpi=200)
               
                else:           
                    plt.show(block=False)
                    #user_input = raw_input('Press any key to exit >>> ')
                    #sys.exit()

    def pltPrf2(self, BinLat='', BinID='', ClrID=''):

        HiPPO_Number = list(set(self.Hno))
        HiPPO_Number = np.asarray(HiPPO_Number, dtype=int)
        HiPPO_Number.sort()

        DateFmt      = DateFormatter('%Y-%m-%d')
        dayLc        = DayLocator()
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)

        layers     = np.arange(0, 15, 1.0)
        layers_mid = np.asarray([(layers[l] + layers[l+1])*0.5  for l in np.arange(0, len(layers)-1, 1)])
        self.layers_mid = layers_mid

        PrfLat    =   {}
        self.PrfLat = {}

        #----------------------------------------------------------------------------------------------
        #THE CTL FILES AND APRIORI ARE READ BELOW 
        #----------------------------------------------------------------------------------------------

        #-----------------------------------------------
        #READING THE CTL FILE AND DEFINE PARAMETERS: THULE
        #-----------------------------------------------
        #ctlFile = '/data1/ebaumer/tab/ocs/x.ocs/sfit4_v2.ctl'
        ctlFile = '/data1/ebaumer/tab/'+self.gasname.lower()+'/x.'+self.gasname.lower()+'/sfit4.ctl'

        (PrimaryGas, ctlinfo) = dc.readCtlF(ctlFile)

        stfile = ctlinfo['file.in.stalayers'][0]
        PrimaryGas = ctlinfo['gas.profile.list'][0]     
        SaTAB = ctlinfo['gas.profile.'+str(PrimaryGas)+'.sigma']
        SaTAB =np.asarray(SaTAB)
        scl = ctlinfo['gas.profile.'+str(PrimaryGas)+'.scale']
        

        #-----------------------------------------------
        #READING THE STATION LAYER FILE: THULE
        #-----------------------------------------------
        ckFile(stfile)
        #stfile = file(stfile, 'r')
        cols, indexToName = mf.getColumns(stfile, headerrow=2, delim=' ', header=True)
        midpointTAB = np.asarray(cols['midpnt'][0:-1]).astype(np.float)

        WACCMfile = '/data/Campaign/TAB/waccm/reference.prf.REFC1.3'
        ckFile(WACCMfile)
        prfwaccm = mf.readRefPrf(fname=WACCMfile, parms = ['ALTITUDE','PRESSURE','TEMPERATURE', PrimaryGas.upper()])
        aprprfTAB =  np.asarray(prfwaccm[PrimaryGas.upper()][0]) * scl 
        zTAB = np.asarray(prfwaccm['ALTITUDE'][0])

        self.SaTAB = SaTAB
        self.aprprfTAB = aprprfTAB
        self.zTAB = zTAB
        self.midpointTAB = midpointTAB

        #-----------------------------------------------
        #READING THE CTL FILE AND DEFINE PARAMETERS: MLO
        #-----------------------------------------------
        ctlFile = '/data1/ebaumer/mlo/ocs/x.ocs/sfit4_v2.ctl'
        #ctlFile = '/data1/ebaumer/mlo/'+self.gasname.lower()+'/x.'+self.gasname.lower()+'/sfit4.ctl'
        (PrimaryGas, ctlinfo) = dc.readCtlF(ctlFile)

        stfile = ctlinfo['file.in.stalayers'][0]
        PrimaryGas = ctlinfo['gas.profile.list'][0]     
        SaMLO = ctlinfo['gas.profile.'+str(PrimaryGas)+'.sigma']
        SaMLO =np.asarray(SaMLO)
        scl = ctlinfo['gas.profile.'+str(PrimaryGas)+'.scale']

        #-----------------------------------------------
        #READING THE STATION LAYER FILE: MLO
        #-----------------------------------------------
        ckFile(stfile)
        #stfile = file(stfile, 'r')
        cols, indexToName = mf.getColumns(stfile, headerrow=2, delim=' ', header=True)
        midpointMLO = np.asarray(cols['midpnt'][0:-1]).astype(np.float)

        WACCMfile = '/data/Campaign/MLO/waccm/reference.prf.REFC1.3'
        ckFile(WACCMfile)
        prfwaccm = mf.readRefPrf(fname=WACCMfile, parms = ['ALTITUDE','PRESSURE','TEMPERATURE', PrimaryGas.upper()])
        aprprfMLO =  np.asarray(prfwaccm[PrimaryGas.upper()][0]) * scl 
        zMLO = np.asarray(prfwaccm['ALTITUDE'][0])

        self.SaMLO = SaMLO
        self.aprprfMLO = aprprfMLO
        self.zMLO = zMLO
        self.midpointMLO = midpointMLO

        #-----------------------------------------------
        #READING THE CTL FILE AND DEFINE PARAMETERS: FL0
        #-----------------------------------------------
        ctlFile = '/data1/ebaumer/fl0/ocs/x.ocs/sfit4_v8.ctl'
        #ctlFile = '/data1/ebaumer/mlo/'+self.gasname.lower()+'/x.'+self.gasname.lower()+'/sfit4.ctl'
        (PrimaryGas, ctlinfo) = dc.readCtlF(ctlFile)

        stfile = ctlinfo['file.in.stalayers'][0]
        PrimaryGas = ctlinfo['gas.profile.list'][0]     
        SaFL0 = ctlinfo['gas.profile.'+str(PrimaryGas)+'.sigma']
        SaFL0 =np.asarray(SaFL0)
        scl = ctlinfo['gas.profile.'+str(PrimaryGas)+'.scale']

        #-----------------------------------------------
        #READING THE STATION LAYER FILE: MLO
        #-----------------------------------------------
        ckFile(stfile)
        #stfile = file(stfile, 'r')
        cols, indexToName = mf.getColumns(stfile, headerrow=2, delim=' ', header=True)
        midpointFL0 = np.asarray(cols['midpnt'][0:-1]).astype(np.float)

        WACCMfile = '/data/Campaign/FL0/waccm/reference.prf.REFC1.3'
        ckFile(WACCMfile)
        prfwaccm = mf.readRefPrf(fname=WACCMfile, parms = ['ALTITUDE','PRESSURE','TEMPERATURE', PrimaryGas.upper()])
        aprprfFL0 =  np.asarray(prfwaccm[PrimaryGas.upper()][0]) * scl 
        zFL0 = np.asarray(prfwaccm['ALTITUDE'][0])

        self.SaMLO = SaMLO
        self.aprprfFL0 = aprprfFL0
        self.zFL0 = zFL0
        self.midpointFL0 = midpointFL0

        #----------------------------------------------------------------------------------------------

    
        for i, b in enumerate(BinLat):

            numprfs = 0
            
            for hn in HiPPO_Number:

                inds1     = np.where(self.Hno == hn)[0]
                
                lat       = self.lat[inds1]
                lon       = self.lon[inds1]
                alt       = self.alt[inds1]
                doy       = self.doy[inds1]
                dt        = self.date[inds1]
                RFHN      = self.flt[inds1]
                gas       = self.gas[inds1]

                RFs = list(set(RFHN))
                RFs = np.asarray(RFs, dtype=int)
                RFs.sort()

                for rf in RFs:
                    inds2     = np.where(RFHN == rf)[0]

                    gas_i     = np.asarray(gas[inds2], dtype=float)
                    alt_i     = np.asarray(alt[inds2])/1000.
                    lat_i     = np.asarray(lat[inds2])


                    if ((np.mean(lat_i) >= b[0]) & (np.mean(lat_i) < b[1]) ):

                        gas_mean_i  = []
                        gas_std_i   = []

                        for l, la in enumerate(layers_mid):
                            inds3     = np.where( (alt_i >= layers[l]) & (alt_i < layers[l+1]))[0]
                            if len(inds3) >=2 : 
                                gas_mean_i.append(np.mean(gas_i[inds3]))
                                gas_std_i.append(np.std(gas_i[inds3]))
                            elif len(inds3) ==1 : 
                                gas_mean_i.append(gas_i[inds3])
                                gas_std_i.append(float('nan'))
                            else: 
                                gas_mean_i.append(float('nan'))
                                gas_std_i.append(float('nan'))

                        gas_mean_i = np.asarray(gas_mean_i, dtype=float)
                        gas_std_i  = np.asarray(gas_std_i, dtype=float)

                        PrfLat.setdefault(BinID[i]+'_Prf',[]).append(gas_mean_i)

                    numprfs += 1
                    
            PrfLatbi   = np.asarray(PrfLat[BinID[i]+'_Prf'])
            gas_mean_i = np.nanmean(PrfLatbi, axis=0)
            gas_var_i  = np.nanvar(PrfLatbi, axis=0)
            gas_std_i  = np.nanstd(PrfLatbi, axis=0)
            gas_med_i  = np.nanmedian(PrfLatbi, axis=0)

            print 'Number of HIPPO profiles in {0}: {1} '.format(b, PrfLatbi.shape)
            print 'Total of HIPPO profiles in {0}: {1}'.format(b, numprfs)

            self.PrfLat.setdefault(BinID[i]+'_Prf_mean',[]).append(gas_mean_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_std',[]).append(gas_std_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_med',[]).append(gas_med_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_var',[]).append(gas_var_i)


        fig, ax = plt.subplots(1,2, figsize=(11, 7))
        fig2, ax2 = plt.subplots(1,2, figsize=(11, 7))

        for i, bi in enumerate(BinID):
            
            PrfLat[bi+'_Prf'] = np.asarray(PrfLat[bi+'_Prf'])
                  
            gas_mean_i = self.PrfLat[bi+'_Prf_mean'][0]  #np.nanmean(PrfLat[bi+'_Prf'], axis=0)
            gas_var_i  = self.PrfLat[bi+'_Prf_var'][0]  #np.nanvar(PrfLat[bi+'_Prf'], axis=0)
            gas_std_i  = self.PrfLat[bi+'_Prf_std'][0]  #np.nanstd(PrfLat[bi+'_Prf'], axis=0)
            gas_med_i  = self.PrfLat[bi+'_Prf_med'][0]  #np.nanmedian(PrfLat[bi+'_Prf'], axis=0)


            ax[0].plot(gas_mean_i, layers_mid, color=ClrID[i], markersize=4)            
            #ax.fill_betweenx(layers_mid,gas_mean_i-gas_std_i,gas_mean_i+gas_std_i,alpha=0.5,color='0.75')
            ax[0].errorbar(gas_mean_i, layers_mid,fmt='o',xerr=gas_std_i, color=ClrID[i], label=bi)
            ax[0].scatter(gas_med_i, layers_mid, edgecolors=ClrID[i], marker='x', facecolors=ClrID[i], s=35)

            Fraction = gas_std_i/gas_mean_i

            ax[1].plot(Fraction, layers_mid, 'o', linestyle='-', color=ClrID[i], markersize=6)


            ax2[0].plot(gas_mean_i, layers_mid, color=ClrID[i], markersize=4)            
            ax2[0].errorbar(gas_mean_i, layers_mid,fmt='o',xerr=gas_std_i, color=ClrID[i], label=bi)
            ax2[0].scatter(gas_med_i, layers_mid, edgecolors=ClrID[i], marker='x', facecolors=ClrID[i], s=35)
  
            ax2[1].plot(Fraction, layers_mid, 'o', linestyle='-', color=ClrID[i], markersize=6, label=bi)
        
        #ax2[0].plot(aprprfTAB*1e12, zTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Thule - Toon')
        #ax2[1].plot(SaTAB, midpointTAB, color='k', linestyle='-', marker ='o', markersize=6, label='Yutting et al, (2016)')

        #ax2[0].plot(aprprfMLO*1e12, zMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='MLO - Toon')
        #ax2[1].plot(SaMLO, midpointMLO, color='gray', linestyle='-', marker ='o', markersize=6, label='Yutting et al, (2016)')

        ax[0].tick_params(labelsize=14)
        ax[0].grid(True)
        ax[0].set_xlabel('VMR [ppt]', fontsize = 16)
        ax[0].set_ylabel('Altitude [km]', fontsize = 16)
        #ax[0].legend(prop={'size':14}, loc=3)
        ax[0].legend(prop={'size':14})
        ax[0].set_ylim(0, 16)
        ax[0].set_xlim(xmin=0.0)
        #ax.title('Vertical Profile of '+ str(self.gasname.upper() )+'\nHIPPO', fontsize=14)

        ax[1].tick_params(labelsize=14)
        ax[1].grid(True)
        ax[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
        ax[1].set_ylim(0, 16)

        ax2[0].tick_params(labelsize=14)
        ax2[0].grid(True)
        ax2[0].set_xlabel('VMR [ppt]', fontsize = 16)
        ax2[0].set_ylabel('Altitude [km]', fontsize = 16)
        ax2[0].set_ylim(0, 20)
        #ax2[0].legend(prop={'size':12}, loc=3)
        ax2[0].legend(prop={'size':12})
        #ax2[0].set_yscale("log", nonposx='clip')
        ax2[0].set_xlim(xmin=0.0)

        ax2[1].tick_params(labelsize=14)
        ax2[1].grid(True)
        ax2[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
        ax2[1].set_ylim(0, 20)
        ax2[1].legend(prop={'size':12}, loc=1)
        #ax2[1].set_xlim(xmin=0.0, xmax=0.15)

       
        fig.subplots_adjust(left = 0.075, bottom=0.1, top=0.95, right = 0.95)
        fig2.subplots_adjust(left = 0.075, bottom=0.1, top=0.95, right = 0.95)
            
        if self.saveFlg: 
           self.pdfsav.savefig(fig,dpi=200)
           self.pdfsav.savefig(fig2,dpi=200)
           self.pdfsav.close()
            
        else:           
            plt.show(block=False)

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()   
               

    def pltMapRF(self):

        HiPPO_Number = list(set(self.Hno))
        HiPPO_Number = np.asarray(HiPPO_Number, dtype=int)
        HiPPO_Number.sort()

        DateFmt      = DateFormatter('%Y-%m-%d')
        dayLc        = DayLocator()
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap) 

        for hn in HiPPO_Number:

            fig, ax = plt.subplots(2, gridspec_kw = {'height_ratios':[3, 1]}, figsize=(11,10))

            #-----------------
            #MAP
            #-----------------
            #m  = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-180,urcrnrlon=180,resolution='c', ax=ax[0])
            # m  = Basemap(projection='mill',llcrnrlat=-90,urcrnrlat=90, llcrnrlon=-180,urcrnrlon=180,resolution='l', ax=ax[0])
            # # plot coastlines, draw label meridians and parallels.
            # m.drawcoastlines()
            # m.drawparallels(np.arange(-90,90,30),labels=[1,0,0,0], fontsize=16)
            # m.drawmeridians(np.arange(m.lonmin,m.lonmax+30,60),labels=[0,0,0,1], fontsize=16)
            # # fill continents 'coral' (with zorder=0), color wet areas 'aqua'
            # m.drawmapboundary(fill_color='aqua')
            # m.fillcontinents(color='coral',lake_color='aqua')

            m = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0, lat_0=0, lon_0=180, ax=ax[0])
            m.drawcoastlines(linewidth =0.25, color = 'gray', )
            #eq_map.drawcountries()
            m.fillcontinents(color = 'gray', alpha=0.5)
            m.drawmapboundary()
            m.drawmeridians(np.arange(0, 360, 30), color='gray', linewidth=0.5, alpha = 0.5)
            m.drawparallels(np.arange(-90, 90, 30), color='gray', linewidth=0.5, labels=[1,0,0,0], alpha = 0.5)

            inds1     = np.where(self.Hno_all == hn)[0]
            
            lat_all   = self.lat_all[inds1]
            lon_all   = self.lon_all[inds1]

            x, y = m(lon_all,lat_all)
            #m.plot(x, y, marker ='o', markersize=4, color='r', linestyle='None', linewidth=3)
            #m.scatter(x, y, facecolors='red', edgecolors='none', s=40, zorder=2)

            inds2  = np.where(self.Hno == hn)[0]
            lat    = self.lat[inds2]
            lon    = self.lon[inds2]
            alt    = self.alt[inds2]
            rf     = self.flt[inds2]

            #norm = colors.Normalize( vmin=np.nanmin(rf), vmax=np.nanmax(rf) )
            bounds = []
            bounds.extend([int(ii+1) for ii in range(np.max(rf))])

            norm = colors.BoundaryNorm(bounds, cm.N )

       
            xx, yy = m(lon,lat)
            #m.plot(xx, yy, marker ='o', markersize=4, color='b', linestyle='None', linewidth=3)
            m.scatter(xx, yy, edgecolors='none', s=30, c=rf,cmap=cm, norm=norm, zorder=2)

            ax[0].tick_params(labelsize=18)

            #-----------------
            #Altitude profile
            #-----------------
            alt   = self.alt[inds2]
            doy   = self.doy[inds2]
            dt    = self.date[inds2]
            
            #ax[1].plot(dt,alt, color='k', linestyle='-', marker ='.', markersize=0)
            sc = ax[1].scatter(dt,alt/1000.,facecolors='blue', edgecolors='none', s=30, c=rf,cmap=cm, norm=norm)
            ax[1].tick_params(labelsize=16)
            
            ax[1].grid(True)
            #ax[1].set_ylim(np.min(vmrP[gv])-0.1*np.min(vmrP[gv]), np.max(vmrP[gv])+0.13*np.max(vmrP[gv]))
            ax[1].set_xlabel('Date', fontsize = 16)
            ax[1].set_ylabel('Altitude [km]', fontsize = 16)
            ax[1].xaxis.set_major_formatter(DateFmt)
            ax[1].xaxis.set_minor_locator(dayLc)
            datemin = datetime.date(dt.min())
            datemax = datetime.date(dt.max())
            ax[1].set_xlim(datemin, datemax)
            fig.autofmt_xdate()

            cax  = fig.add_axes([0.88, 0.15, 0.03, 0.8])
                
            cbar = fig.colorbar(sc, cax=cax,  orientation='vertical')

            cbar.set_ticks(np.asarray(bounds)+0.5)
            cbar.set_ticklabels(bounds)
            cbar.set_label('Research Flight', fontsize=16)  
            
            
            fig.suptitle('HIPPO-'+str(hn)+ ' ('+str(datemin) + ' to '+str(datemax)+')', fontsize=20)
            fig.subplots_adjust(left = 0.12, bottom=0.15, top=0.95, right = 0.85)
            
            if self.saveFlg: 
               self.pdfsav.savefig(fig,dpi=200)
               
            else:           
                plt.show(block=False)

        if self.saveFlg: self.pdfsav.close()

class ACEClass(_DateRange):

    def __init__(self,dataDir, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,  outFname='', saveFlg= False):

        self.idate = dt.date( int(iyear), int(imnth), int(iday))
        self.fdate = dt.date( int(fyear), int(fmnth), int(fday))

        if not( dataDir.endswith('/') ):
            dataDir = dataDir + '/' 

        self.dataDir   = dataDir

        #-------------------------------------------------------------------------
        #BELOW IS COMMENTED OU BECUASE IT WAS THE INITIAL ATTEMP TO READ INDIVIDUAL ASCCI/NETCDF FILES
        #THE CURRENT APPROACH IS TO READ SINGLE NETCDF FILES FROM V3.5 
        #WEBSITE: https://databace.scisat.ca/level2/ace_v3.5/
        #username: iortega, pw: grindcore
        #-------------------------------------------------------------------------
        
        # check if '/' is included at end of path
        # self.dataDir   = dataDir
        # self.dirLst    = []  
        # #---------------------------------
        # # Test if date range given. If so 
        # # create a list of directories to 
        # # read data from
        # #---------------------------------
        # if all([iyear,imnth,iday,fyear,fmnth,fday]):

        #     self.dirDate = []
            
        #     # check if '/' is included at end of path
        #     if not( dataDir.endswith('/') ):
        #         dataDir = dataDir + '/'           
                
        #     # Check if directory exits
        #     dc.ckDir(dataDir,exitFlg=True)
            
        #     _DateRange.__init__(self, iyear, imnth, iday, fyear, fmnth, fday, incr=1)

        #     #--------------------------------------------
        #     # Walk through first level of directories
        #     #--------------------------------------------
        #     for drs in os.walk(dataDir).next()[1]: 
                
        #         #-------------------------------------------
        #         # Test directory to make sure it is a number
        #         #-------------------------------------------
        #         try:    int(drs[0:4])
        #         except: continue

        #         if _DateRange.inRange(self, int(drs[0:4]), int(drs[5:7]), 1 ):
                
        #             self.dirDate.append(dt.date(int(drs[0:4]), int(drs[5:7]), 1 ) )

        #             #idfile = '*.nc'
        #             idfile = '*v3.5.asc'
        #             files = glob.glob( dataDir + drs +'/' + idfile)

        #             files.sort()
                    
        #             self.dirLst.append(files)            
            
        #     #----------------------------------------------------------------------------
        #     # This is important in keeping all the gathered data in order
        #     #----------------------------------------------------------------------------
        #     self.dirLst.sort()
            

        if saveFlg:
            self.pdfsav = PdfPages(outFname)
            self.saveFlg   = True
        else: self.saveFlg = False

    def ReadNCDFACE(self, gasname=''):
        '''Function to Read out the net CDF files.'''

        self.lat     = []
        self.lon     = []
        self.gas     = []
        self.alt     = []
    

        for indMain, Dir in enumerate(self.dirLst):

            lat  = []
            lon  = []
            gas  = []
            alt  = []


            for sngDir in Dir:
                

                try:
                    #print sngDir
                    nc_fid = Dataset(sngDir , 'r')  # Dataset is the class behavior to open the file
                    #nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=True)

                    vgrp =  str(nc_fid.groups)

                    #print vgrp
                    #print vgrp[746:752]
                    #print vgrp[767:773]


                    try:
                        lat_i   =  float(vgrp[745:752])
                        lon_i   = float(vgrp[767:773])
                    except ValueError:

                        try:
                            lat_i    = float(vgrp[746:752])
                            lon_i   = float(vgrp[768:773])
                        except ValueError:
                            lat_i = 'not good'
                            lon_i = 'not good'

                    if isinstance(lat_i, float):

                        gas_i   =  np.asarray(nc_fid['/ACE-FTS-v2.2/Data-L2_1km_grid/'+gasname.upper()])
                        alt_i   =  np.asarray(nc_fid['/ACE-FTS-v2.2/Data-L2_1km_grid/z'])

                        ind = np.where(gas_i != -999)[0] #

                        gas_i  = gas_i[ind]*1e12   #vmr to ppt
                        alt_i  = alt_i[ind]

                        lat.append(lat_i)
                        lon.append(lon_i)
                        gas.append(gas_i)
                        alt.append(alt_i)

                except IOError as e:
                    #print "I/O error({0}): {1}".format(e.errno, e.strerror)
                    pass

            self.lat.extend(lat)
            self.lon.extend(lon)
            self.gas.extend(gas)
            self.alt.extend(alt)

        self.lat   = np.asarray(self.lat)
        self.lon   = np.asarray(self.lon)
        self.gas   = np.asarray(self.gas)
        self.alt   = np.asarray(self.alt)

        self.lat.flatten()
        self.lon.flatten()
        self.gas.flatten()
        self.alt.flatten()


    def ReadNCDFACE2(self, gasname=''):
        
        '''Function to Read out the net CDF files.'''
        #-------------------------------------------------------------------------
        #THE CURRENT APPROACH IS TO READ SINGLE NETCDF FILES FROM V3.5 
        #WEBSITE: https://databace.scisat.ca/level2/ace_v3.5/
        #username: iortega, pw: grindcore
        #-------------------------------------------------------------------------
        
        self.gas        = []
        self.alt        = []
        self.date       = []
        self.pres       = []
        self.temp       = []

        if gasname.lower() == 'ocs': 
            nc_fid    = Dataset(self.dataDir+'ACEFTS_L2_v3p5_OCS.nc' , 'r')  # Dataset is the class behavior to open the file
        else:
            print 'Error: Look if the NETCDF of {0} is in the directory'.format(gasname)
        
        c_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=False)

        orbit     = np.asarray(nc_fid.variables['orbit'])
        lat       = np.asarray(nc_fid.variables['latitude'])
        lon       = np.asarray(nc_fid.variables['longitude'])
        qfl       = np.asarray(nc_fid.variables['quality_flag'])
        year      = np.asarray(nc_fid.variables['year'])
        month     = np.asarray(nc_fid.variables['month'])
        day       = np.asarray(nc_fid.variables['day'])
        hour      = np.asarray(nc_fid.variables['day'])

        date      = [dt.date(y, month[i], day[i]) for i, y in enumerate(year)]
        date      = np.asarray(date)

        alt       = np.asarray(nc_fid.variables['altitude'])
        gas       = np.asarray(nc_fid.variables[gasname.upper()])
        Temp      = np.asarray(nc_fid.variables['temperature'])
        Pres      = np.asarray(nc_fid.variables['pressure'])

        indD      = np.where((date >= self.idate) & (date <= self.fdate))[0]

        print 'Number of ACE-FTS within the dates ({0} - {1}) of intereste: {2}'.format(date[indD][0], date[indD][-1], len(indD))
        print 'Total number of ACE-FTS within the dates ({0} - {1}) of intereste: {2}'.format(date[0], date[-1], len(date))

        self.gas_all    = np.asarray(gas[indD])*1e12
        self.alt_all    = np.asarray(alt)
        self.temp_all   = np.asarray(Temp[indD])
        self.pres_all   = np.asarray(Pres[indD])

        self.orbit      = np.asarray(orbit[indD])
        self.lat        = np.asarray(lat[indD])
        self.lon        = np.asarray(lon[indD])
        self.qfl        = np.asarray(qfl[indD])
        self.date       = np.asarray(date[indD])

        self.gas        = np.asarray(gas[indD])*1e12
        self.alt        = np.asarray(alt)
        self.temp       = np.asarray(Temp[indD])
        self.pres       = np.asarray(Pres[indD])


        for ii, qf in enumerate(self.qfl):
            #-------------------------------------------------------------------------
            #UPDATE WITH NAN BAD VALUES IN GAS PROFILES
            #SEE Data Usage and File Format Document: https://databace.scisat.ca/level2/ace_v3.5/
            #username: iortega, pw: grindcore
            #-------------------------------------------------------------------------
            ind = np.where( (qf > 3) )[0]

            self.gas[ii][ind]   = float('nan')
            self.temp[ii][ind]  = float('nan') 
            self.pres[ii][ind]  = float('nan')


        self.gas        = np.asarray(self.gas)
        self.alt        = np.asarray(self.alt)
        self.temp       = np.asarray(self.temp)
        self.pres       = np.asarray(self.pres)

        #-------------------------------------------------
        # Calculating density profile, layer thickness and midpoint
        #-------------------------------------------------
        self.rho_all = mf.AirDensity(self.pres_all, self.temp_all,  Punits='atm', Tunits='K' )
        self.rho     = mf.AirDensity(self.pres, self.temp,  Punits='atm', Tunits='K' )
        #self.rho     = [mf.AirDensity(self.pres[i], self.temp[i],  Punits='atm', Tunits='K' ) for i, p in enumerate(self.gas)]

        #-------------------------------------------------
        # Calculating layer thickness, airmass and profiles in molec/cm2.
        #-------------------------------------------------
        self.dz_all     = np.absolute ((self.alt_all[1:] - self.alt_all[:-1]))*1000.*100.  #in cm
        self.dz_all    = np.append(self.dz_all, self.dz_all[-1])
        self.airmass_all  = self.rho_all*self.dz_all
        self.rPrfMol_all  = np.asarray(self.gas_all/1e12) * self.airmass_all
      
        self.airmass  = self.rho*(1.0*1000.*100.)
        self.rPrfMol  = np.asarray(self.gas/1e12) * self.airmass


    def ReadASCIIACE(self, gasname=''):
        '''----------------------------------------------------------------------------
           Below is a script to read the ASCCII FILES: NOT USED AFTER I FOUND SINGLE NETCDF FILES
        ---------------------------------------------------------------------------'''  
        retrvdAll   = ['z', 'T', 'P', 'dens', gasname.upper(), gasname.upper()+'_err']
       
        self.lat        = []
        self.lon        = []
        self.gas        = []
        self.alt        = []
        self.name       = []
        self.start_time = []
        self.end_time   = []
        self.date       = []

        self.deflt   = {}


        for indMain, Dir in enumerate(self.dirLst):


            lat  = []
            lon  = []
            gas  = []
            alt  = []

            for sngDir in Dir:
                try:
                    with open(sngDir , 'r') as fopen:
                        
                        defltLines = fopen.readlines()


                        self.name.append(defltLines[0].strip().split()[2])
                        self.start_time.append(defltLines[3].strip().split()[2])
                        self.end_time.append( defltLines[4].strip().split()[2])
                        
                        self.lat.append(float(defltLines[6].strip().split()[2]))
                        self.lon.append(float(defltLines[7].strip().split()[2]))

                        strdate = defltLines[5].strip().split()[2]
                        strdate = strdate.strip().split('-')

                        self.date.append(dt.date(int(strdate[0]), int(strdate[1]), int(strdate[2])))

                        #--------------------------------
                        # Get Names of profiles retrieved
                        #--------------------------------
                        defltParm = defltLines[11].strip().split()

                        #gas_i  = [ float(row.strip().split()[defltParm.index( gasname.upper())]) for row in defltLines[12:] ]
                        gas_i  = [ float(row.strip().split()[25]) for row in defltLines[12:]  ] 
                        alt_i  = [ float(row.strip().split()[defltParm.index( 'z')]) for row in defltLines[12:] ]

                        gas_i = np.asarray(gas_i)*1e12
                        alt_i = np.asarray(alt_i)

                        ind = np.where(gas_i >= 0.0)[0]

                        self.gas.append(gas_i[ind])
                        self.alt.append(alt_i[ind])

                        for rtrvdSing in retrvdAll:

                            df = [ float(row.strip().split()[defltParm.index(rtrvdSing)]) for row in defltLines[12: ] ]
                            df = np.asarray(df)
                            df = df[ind]
                            self.deflt.setdefault(rtrvdSing,[]).append(df )       

                except Exception as errmsg:
                    print errmsg
                    continue


        self.lat   = np.asarray(self.lat)
        self.lon   = np.asarray(self.lon)
        self.gas   = np.asarray(self.gas)
        self.alt   = np.asarray(self.alt)
        self.date   = np.asarray(self.date)

        for rtrvdSing in retrvdAll:
            self.deflt[rtrvdSing] =  np.asarray(self.deflt[rtrvdSing])



    def pltPrfACE(self, BinLat='', BinID='', ClrID=''):

        #----------------------------------------------------------------------------
        #CALCULATING MEAN PROFILES, AND CREATING PLOTS
        #---------------------------------------------------------------------------

        layers     = np.arange(14, 30, 1.0)
        layers_mid = np.asarray([(layers[l] + layers[l+1])*0.5  for l in np.arange(0, len(layers)-1, 1)])
        self.layers_mid = layers_mid

        self.PrfLat = {}

        yearsLc      = YearLocator()
        #months       = MonthLocator(bymonth=1,bymonthday=1)
        months       = MonthLocator()
        DateFmt      = DateFormatter('%m\n%Y')   

        fig, ax = plt.subplots(1,2, figsize=(11, 7))   #MEAN PROFILES BINNED BY LATITUDE
        fig2, ax2 = plt.subplots(figsize=(10, 6))      #TIME SERIES OF PARTIAL COLUMNS
        fig3, ax3 = plt.subplots(figsize=(10, 6))      #PARTIAL COLUMN BY MONTH

        for i, bi in enumerate(BinLat):

            indS       = np.where( (self.lat >= bi[0]) &  (self.lat < bi[1]))[0]
            print 'Number of ACE profiles in {0}: {1} '.format(bi, len(indS))

            gas_int  = []

            #alt = self.alt[indS]
            alt = self.alt
            gas = self.gas[indS]


            for ii in range(len(indS)):
                try:
                    #gint  = interpolate.interp1d(alt[ii], gas[ii], bounds_error=False, axis=0)(layers_mid)
                    gint  = interpolate.interp1d(alt, gas[ii], bounds_error=False, axis=0)(layers_mid)
                    gas_int.append(gint)
                except ValueError:
                    pass

            gas_int = np.asarray(gas_int)

            gas_mean_i = np.nanmean(gas_int, axis=0)
            gas_var_i  = np.nanvar(gas_int, axis=0)
            gas_std_i  = np.nanstd(gas_int, axis=0)
            gas_med_i  = np.nanmedian(gas_int, axis=0)

            self.PrfLat.setdefault(BinID[i]+'_Prf_mean',[]).append(gas_mean_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_std',[]).append(gas_std_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_med',[]).append(gas_med_i)
            self.PrfLat.setdefault(BinID[i]+'_Prf_var',[]).append(gas_var_i)
            #----------------------------------------------------------------------------
            #MEAN PROFILES
            #---------------------------------------------------------------------------

            ax[0].plot(gas_mean_i, layers_mid,  color=ClrID[i], markersize=4)
            ax[0].errorbar(gas_mean_i, layers_mid,fmt='o',xerr=gas_std_i, color=ClrID[i], label=BinID[i])
            ax[0].scatter(gas_med_i, layers_mid, edgecolors=ClrID[i],marker='x', facecolors=ClrID[i], s=35)

            Fraction = gas_std_i/gas_mean_i

            ax[1].plot(Fraction, layers_mid, 'o', linestyle='-', markersize=6, color=ClrID[i], label=BinID[i])

            ax[0].tick_params(labelsize=14)  
            ax[0].grid(True)
            ax[0].set_xlabel('VMR [ppt]', fontsize = 14)
            ax[0].set_ylabel('Altitude [km]', fontsize = 14)
            ax[0].set_ylim(0,35)
            ax[0].legend(prop={'size':14}, loc=3)

            ax[1].tick_params(labelsize=14)
            ax[1].grid(True)
            ax[1].set_xlabel('Standard deviation / Mean ', fontsize = 16)
            ax[1].set_ylim(0, 40)
            ax[1].legend(prop={'size':14}, loc=1)
            
            fig.subplots_adjust(left = 0.125, bottom=0.15, top=0.95, right = 0.95)

            #----------------------------------------------------------------------------
            #PARTIAL COLUMNS: ALTITUDE IS HARDCODED BELOW
            #---------------------------------------------------------------------------
            ind1 = mf.nearestind(10, self.alt)
            ind2 = mf.nearestind(25, self.alt)

            sumP = np.sum(self.rPrfMol[indS,ind1:ind2], axis=1)

            ax2.plot(self.date[indS], sumP,'k.', markersize=7, color=ClrID[i],  label=BinID[i])
            ax2.grid(True)
            ax2.set_ylabel('Partial Column [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            ax2.set_title('Altitude Layer '+str(self.alt[ind1])+'[km] - '+str(self.alt[ind2])+'[km]',
                          multialignment='center',fontsize=14)
            ax2.set_xlabel('Date', fontsize=14)
            ax2.xaxis.set_major_locator(yearsLc)
            ax2.xaxis.set_minor_locator(months)
            ax2.xaxis.set_major_formatter(DateFmt)
            ax2.tick_params(labelsize=14)
            ax2.legend(prop={'size':14})

            fig2.subplots_adjust(left = 0.1, bottom=0.15, top=0.9, right = 0.95)

            #----------------------------
            # Plot total columns by month
            #----------------------------
            month    = np.array([d.month for d in self.date[indS]])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))
            
            for ii,m in enumerate(mnthSort):
                indsM        = np.where(month == m)[0]
                mnthMean[ii]  = np.nanmean(sumP[indsM])
                mnthSTD[ii]   = np.nanstd(sumP[indsM])   
                
            
            ax3.plot(mnthSort,mnthMean, linestyle='-', color=ClrID[i],  label=BinID[i])
            ax3.errorbar(mnthSort,mnthMean,yerr=mnthSTD,fmt='k.',markersize=10, color= ClrID[i], ecolor=ClrID[i])     
            ax3.grid(True,which='both')
            ax3.set_ylabel('Partial Column [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            ax3.set_xlabel('Month', fontsize=14)
            ax3.set_title('Retrieved Monthly Mean with Standard Deviation\n Altitude Layer '+str(self.alt[ind1])+'[km] - '+str(self.alt[ind2])+'[km]', fontsize=14)
            ax3.set_xlim((0,13))
            ax3.set_xticks(range(1,13))
            ax3.tick_params(labelsize=14)
            ax3.legend(prop={'size':14})
            fig3.subplots_adjust(left = 0.1, bottom=0.15, top=0.9, right = 0.95)


            
        if self.saveFlg: 
           self.pdfsav.savefig(fig,dpi=200)
           self.pdfsav.savefig(fig2,dpi=200)
           self.pdfsav.savefig(fig3,dpi=200)
        
        else:           
            plt.show(block=False)


        # #----------------------------------------------------------------------------
        # #CONTOUR PLOT: TIME SERIES OF VERTICAL PROFILES
        # #---------------------------------------------------------------------------

        # levels1 = np.arange(0., 600,  50.)        
        # fig,ax1 = plt.subplots(figsize=(10, 6))
        
        # cax1          = ax1.contourf(self.date, self.alt, np.transpose(self.gas), levels1,  cmap=mplcm.jet)
        
        # divider1      = make_axes_locatable(ax1)
        # cb1           = divider1.append_axes("right",size="5%",pad=0.075)

        # cbar1         = plt.colorbar(cax1,cax=cb1)
        # cbar1.set_label('VMR [ppt]', fontsize=12)
        # cbar1.ax.tick_params(labelsize=12)
       
        # ax1.grid(True)
        # ax1.set_xlabel('Date',  fontsize=14)
        # ax1.set_ylabel('Altitude [km]',  fontsize=14)   
        # ax1.set_ylim((0,35))     
        # ax1.tick_params(labelsize=14)
        # ax1.xaxis.set_major_locator(yearsLc)
        # ax1.xaxis.set_minor_locator(months)
        # ax1.xaxis.set_major_formatter(DateFmt)
        # ax1.set_title('Time Series of vertical Profiles', fontsize=14)
        # fig.subplots_adjust(left = 0.1, bottom=0.15, top=0.9, right = 0.92)

        # if self.saveFlg: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.show(block=False)  

        # #----------------------------------------------------------------------------
        # #LATITUDE/ALTITUDE SECTIONS OF PROFILES
        # #---------------------------------------------------------------------------
        # fig,ax1 = plt.subplots(figsize=(10, 6))
       
        # cax1          = ax1.contourf(self.lat[:,0], self.alt, np.transpose(self.gas), levels1,  cmap=mplcm.jet)
        
        # divider1      = make_axes_locatable(ax1)
        # cb1           = divider1.append_axes("right",size="5%",pad=0.075)

        # cbar1         = plt.colorbar(cax1,cax=cb1)
        # cbar1.set_label('VMR [ppt]', fontsize=12)
        # cbar1.ax.tick_params(labelsize=12)
       
        # ax1.grid(True)
        # ax1.set_xlabel('Latitude',  fontsize=14)
        # ax1.set_ylabel('Altitude [km]',  fontsize=14)   
        # ax1.set_ylim((0,35))     
        # ax1.tick_params(labelsize=14)
        # ax1.set_title('Latitude/Altitude sections', fontsize=14)
        # fig.subplots_adjust(left = 0.1, bottom=0.1, top=0.9, right = 0.92)
       
        
        # if self.saveFlg: 
        #     self.pdfsav.savefig(fig,dpi=200)
        # else:
        #     plt.show(block=False)


        # if self.saveFlg: 
        #     self.pdfsav.close()
        # else:
        #     plt.show(block=False)   
                






               

            