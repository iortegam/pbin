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


#------------------------------------------------------------------------------------------------------------------------------                     
class SondeClass(_DateRange):

    def __init__(self, dataDirSonde,loc, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1, fleoutFlg=False):       
        
        self.dirLstsonde           = []   #DIRECTORY WITH THE SONDE FILES

        self.loc = loc.lower()
        
        #---------------------------------
        # Test if date range given. If so 
        # create a list of directories to 
        # read data from
        #---------------------------------
        if all([iyear,imnth,iday,fyear,fmnth,fday]):
            
            # check if '/' is included at end of path
            if not( dataDirSonde.endswith('/') ):
                dataDirSonde = dataDirSonde + '/'  

            self.dataDirSonde          = dataDirSonde         
                
            # Check if directory exits
            dc.ckDir(dataDirSonde,exitFlg=True)
            
            _DateRange.__init__(self, iyear, imnth, iday, fyear, fmnth, fday, incr=1)


            #--------------------------------------------
            # Walk through first level of directories and
            # collect directory names for processing
            #--------------------------------------------

            for dd in self.dateList:
                # Find year month and day strings
                yrstr   = "{0:02d}".format(dd.year)
                mnthstr = "{0:02d}".format(dd.month)
                daystr  = "{0:02d}".format(dd.day)

                if loc.lower() == 'mlo':   
                    filename = 'HIH_H2O_'+yrstr+mnthstr+daystr+'*.txt'
                    if fleoutFlg: filenamefleout =  '*'+yrstr+'_'+mnthstr+'_'+daystr+'*'
                elif loc.lower() == 'fl0': 
                    filename = 'BLD_H2O_'+yrstr+mnthstr+daystr+'*.txt'
                    if fleoutFlg: filenamefleout =  '*'+yrstr+'_'+mnthstr+'_'+daystr+'*'

                
                if fleoutFlg: sondefilefleout = glob.glob( dataDirSonde + filenamefleout )

                sondefile = glob.glob( dataDirSonde + filename )

                if not sondefile: continue

                for k in sondefile:
                    self.dirLstsonde.append(k)

            self.dirLstsonde.sort()
            self.fleoutFlg = fleoutFlg


    def readfilessonde(self):
        ''' Reads in txt files from sondes '''
        self.sonde  = {}
        self.fleout = {}

        #------------------------------------
        # Loop through collected directories
        # self.dirLst has already been sorted
        #------------------------------------
        for indMain,sngDir in enumerate(self.dirLstsonde):

            #-----------------------------------------
            #Reading FPH Sonde
            #-----------------------------------------
            if not os.path.isfile(sngDir): continue

            with open(sngDir ,'r') as fopen: lines = fopen.readlines()
            info = [ row.strip().split() for row in lines if 'Water Vapor Flight Date' in row]
            info2 = [ row.strip().split() for row in lines if 'Number of header lines' in row] 
            info3 = [ row.strip().split() for row in lines if 'Number of variables' in row]

            nheaders = int(info2[-1][-1])
            nvariable = int(info3[-1][-1])

            lines[:] = [ row.strip().split() for row in lines[nheaders:-1] ] 
           
            #lines[:] = [ row.strip().split() for row in lines if len(row.strip().split()) == 13 and not 'to' in row and not 'Level' in row and not 'Traditionally' in row and not 'Number' in row]              

            npoints = len(lines)

            time = [row[7].split()[0] for row in info]
            
            strfile = [x for x in sngDir.split('/')]
            filename = strfile[-1]
            yyyy = filename[8:12]
            mm   = filename[12:14]
            dd   = filename[14:16]
           
            time2 = time[0].split(':')
            
            hh = int(time2[0])
            mi = int(time2[1])
            se = int(time2[2])  

            testQ = [row[7] for row in lines if float(row[7]) > 0.0]
            testQ = np.asarray(testQ)

            if len(testQ) > 2:

                self.sonde.setdefault('Level_Number',[]).append([row[0] for row in lines if float(row[7]) > 0.0])  #and filtering when h2o is negative 
                self.sonde.setdefault('Alt_km',[]).append([row[1] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Press_hPa',[]).append([row[2] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Temp_degC',[]).append([row[3] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('Theta_degK',[]).append([row[4] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('RH_FPH',[]).append([row[5] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('RH_RS',[]).append([row[6] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Omr_ppmv',[]).append([row[7] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Osd_ppmv',[]).append([row[8] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('H2Omr_orig_ppmv',[]).append([row[9] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('TFP_FPH_decC',[]).append([row[10] for row in lines if float(row[7]) > 0.0])
                self.sonde.setdefault('TFP_RS_degC',[]).append([row[11] for row in lines if float(row[7]) > 0.0])
                if nvariable > 13: self.sonde.setdefault('O3mr_ppbv',[]).append([row[12] for row in lines if float(row[7]) > 0.0])
                else: self.sonde.setdefault('O3mr_ppbv',[]).append([])
                self.sonde.setdefault('dt',[]).append(dt.datetime( int(yyyy), int(mm), int(dd), hh, mi, se  ))
                self.sonde.setdefault('date',[]).append(dt.date(int(yyyy),  int(mm), int(dd) ))

                if self.loc.lower() == 'fl0':

                    #-----------------------------------------
                    #Reading fleout Sonde
                    #-----------------------------------------
                    filename =  '*'+yyyy+'_'+mm+'_'+dd+'*'

                    sondefilefleout = glob.glob( self.dataDirSonde + filename )

                    if sondefilefleout:

                        with open(sondefilefleout[0] ,'r') as fname: 
                            lines = fname.readlines()
                        
                            info = [ row.strip().split() for row in lines if 'Header lines' in row]
                            info2 = [ row.strip().split() for row in lines if 'Data columns' in row] 

                            nheaders = int(info[-1][-1])
                            nvariable = int(info2[-1][-1])

                            lines[:] = [ row.strip().split(',') for row in lines[nheaders-2:-1] ] 
                            headings = lines[0]

                            headings = [f.strip() for f in headings]
                           
                            indsAlt = [i for i, s in enumerate(headings) if s == 'Alt']#np.where(headings == 'Alt')[0]
                            indsLat = [i for i, s in enumerate(headings) if s == 'GPS lat']
                            indsLon = [i for i, s in enumerate(headings) if s == 'GPS lon']

                            
                            if len(indsLat) == 1:
                                self.fleout.setdefault('Alt',[]).append([row[int(indsAlt[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])
                                self.fleout.setdefault('GPS lat',[]).append([row[int(indsLat[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])
                                self.fleout.setdefault('GPS lon',[]).append([row[int(indsLon[0])] for row in lines[2:-1] if float(row[int(indsLat[0])]) > 0.0])

                            else:

                                self.fleout.setdefault('Alt',[]).append([np.nan for row in lines[2:-1]])
                                self.fleout.setdefault('GPS lat',[]).append([np.nan for row in lines[2:-1]])
                                self.fleout.setdefault('GPS lon',[]).append([np.nan for row in lines[2:-1]])

                    else: 

                        self.fleout.setdefault('Alt',[]).append([np.nan for row in lines if float(row[7]) > 0.0])
                        self.fleout.setdefault('GPS lat',[]).append([np.nan for row in lines if float(row[7]) > 0.0])
                        self.fleout.setdefault('GPS lon',[]).append([np.nan for row in lines if float(row[7]) > 0.0])


    def sondePrf(self):
        ''' Post analysis of sonde profiles, i.e., ascent/descent profiles, concentration profiles '''


        for a, alt in enumerate(self.sonde['Alt_km']):


            #------------------------------------------------
            #First, calculate density and concentration profiles in molec/cm2
            #------------------------------------------------
            pres         = np.asarray(self.sonde['Press_hPa'][a], dtype=np.float32) 
            Temp         = np.asarray(self.sonde['Temp_degC'][a], dtype=np.float32)
            H2Omr        = np.asarray(self.sonde['H2Omr_ppmv'][a], dtype=np.float32)
            H2Osd        = np.asarray(self.sonde['H2Osd_ppmv'][a], dtype=np.float32)
            altsonde     = np.asarray(alt, dtype=np.float32)
            Temp         = Temp + 273.15
            rho          = (pres*100.00) / ((Temp)*8.314) * 6.022e23 /100. /100. /100.
            self.sonde.setdefault('rho',[]).append(rho)
            
            
            #------------------------------------------------
            #Midpoints
            #------------------------------------------------
            rho_mid      = np.absolute ((rho[1:] + rho[:-1]) / 2.0)
            altsonde_mid = np.absolute ((altsonde[1:] + altsonde[:-1]) / 2.0)
            H2Omr_mid    = np.absolute ((H2Omr[1:] + H2Omr[:-1]) / 2.0)
            H2Osd_mid    = np.absolute ((H2Osd[1:] + H2Osd[:-1]) / 2.0)

            dz           = np.absolute ((altsonde[1:] - altsonde[:-1])) *1000.*100.
          
            airmass      = rho_mid*dz
            self.sonde.setdefault('airmass',[]).append(airmass)
            CPrf_mid     = (H2Omr_mid*1e-6)*rho_mid
            CPrf_sd_mid  =  (H2Osd_mid*1e-6) *rho_mid
            
            H2OMol_mid   = (H2Omr_mid*1e-6)*airmass
            H2OMolsd_mid = (H2Osd_mid*1e-6)*airmass

            self.sonde.setdefault('H2Omr_ppmv_mid',[]).append(H2Omr_mid)
            self.sonde.setdefault('H2Osd_ppmv_mid',[]).append(H2Osd_mid)
            self.sonde.setdefault('H2OMol_mid',[]).append(H2OMol_mid)
            self.sonde.setdefault('H2OMolsd_mid',[]).append(H2OMolsd_mid)
            self.sonde.setdefault('CPrf_mid',[]).append(CPrf_mid)
            self.sonde.setdefault('CPrf_sd_mid',[]).append(CPrf_sd_mid)
            self.sonde.setdefault('Alt_km_mid',[]).append(altsonde_mid)
            
    
            #------------------------------------------------
            #Calculate ascent and descent profiles
            #----------------------------------------------
            inds = np.where(altsonde_mid == np.amax(altsonde_mid))[0]
            
            altsonde_a    = np.asarray(altsonde_mid[0:inds[0]] , dtype=np.float32)
            altsonde_d    = np.asarray(altsonde_mid[inds[0]+1:] , dtype=np.float32)           

            self.sonde.setdefault('H2Omr_ppmv_a',[]).append(H2Omr_mid[0:inds[0]])
            self.sonde.setdefault('H2Osd_ppmv_a',[]).append(H2Osd_mid[0:inds[0]])
            self.sonde.setdefault('Alt_km_a',[]).append(altsonde_a)
            self.sonde.setdefault('CPrf_a',[]).append(CPrf_mid[0:inds[0]])
            self.sonde.setdefault('rho_a',[]).append(rho_mid[0:inds[0]])
            self.sonde.setdefault('H2OMol_a',[]).append(H2OMol_mid[0:inds[0]])
            self.sonde.setdefault('H2OMolsd_a',[]).append(H2OMolsd_mid[0:inds[0]])
            self.sonde.setdefault('airmass_a',[]).append(airmass[0:inds[0]])


            self.sonde.setdefault('H2Omr_ppmv_d',[]).append(H2Omr_mid[inds[0]+1:])
            self.sonde.setdefault('H2Osd_ppmv_d',[]).append(H2Osd_mid[inds[0]+1:])
            self.sonde.setdefault('Alt_km_d',[]).append(altsonde_d)
            self.sonde.setdefault('CPrf_d',[]).append(CPrf_mid[inds[0]+1:])
            self.sonde.setdefault('rho_d',[]).append(rho_mid[inds[0]+1:])
            self.sonde.setdefault('H2OMol_d',[]).append(H2OMol_mid[inds[0]+1:])
            self.sonde.setdefault('H2OMolsd_d',[]).append(H2OMolsd_mid[inds[0]+1:])
            self.sonde.setdefault('airmass_d',[]).append(airmass[inds[0]+1:])
            
            #------------------------------------------------
            #Interpolate fleout lat/lon to FPH altitude
            #----------------------------------------------
            if self.fleoutFlg:

                try: 
                
                    AltFleout    = np.asarray(self.fleout['Alt'][a], dtype=np.float32)
                    latFleout    = np.asarray(self.fleout['GPS lat'][a], dtype=np.float32)  
                    lonFleout    = np.asarray(self.fleout['GPS lon'][a], dtype=np.float32) 

                    latFleout_a       = interpolate.interp1d(AltFleout, latFleout, bounds_error=False)(altsonde_a)
                    latFleout_d       = interpolate.interp1d(AltFleout, latFleout, bounds_error=False)(altsonde_d)

                    self.sonde.setdefault('latFleout_a',[]).append(latFleout_a)
                    self.sonde.setdefault('latFleout_d',[]).append(latFleout_d)

                    lonFleout_a       = interpolate.interp1d(AltFleout, lonFleout, bounds_error=False)(altsonde_a)
                    lonFleout_d       = interpolate.interp1d(AltFleout, lonFleout, bounds_error=False)(altsonde_d)

                    self.sonde.setdefault('lonFleout_a',[]).append(lonFleout_a)
                    self.sonde.setdefault('lonFleout_d',[]).append(lonFleout_d)

                except Exception as errmsg:
                    #print errmsg
                  
                    nNan    = len(altsonde_a)
                    tempNan = np.zeros([nNan])
                    tempNan[:] = np.nan

                    
                    self.sonde.setdefault('latFleout_a',[]).append(tempNan)
                    self.sonde.setdefault('latFleout_d',[]).append(tempNan)

                    self.sonde.setdefault('lonFleout_a',[]).append(tempNan)
                    self.sonde.setdefault('lonFleout_d',[]).append(tempNan)


class FTSClass():

    def __init__(self, gasName, retDir, ctlF, loc, ver, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,outFname='', fleoutFlg=False, version=''):\
    

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

        self.ver      = ver
        self.retDir   = [ retDir+v+'/' for v in ver] 
        self.ctlFile  = [retDir+'x.'+gasName.lower()+'/'+ cf for cf in ctlF]
        
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
               pcFlg=True,cnvrgFlg=True,allGas=False,sclfct=1.0,sclname='ppv',pltStats=True,szaFlg=False,errFlg=False,chiFlg=False,tcMMflg=False,mnthFltFlg=False, pCols=False, smthFlg=False,
               doiplt=False, loc='', latFTS=False, lonFTS=False,
               altCG=False, diffT= False, dfValue=False):

        #-------------------------------------
        # Create instance of output data class   
        #-------------------------------------
        statDataCl = OrderedDict()
        for i, v in enumerate(self.ver):
            statDataCl[v] = dc.ReadOutputData(self.retDir[i],'',self.ctlFile[i],self.iyear,self.imnth,self.iday,self.fyear,self.fmnth,self.fday)

        #--------------
        # Read profiles
        #--------------
        for ver in statDataCl:
            statDataCl[ver].readprfs([statDataCl[ver].PrimaryGas],retapFlg=1)
            statDataCl[ver].readprfs([statDataCl[ver].PrimaryGas],retapFlg=0)
        
            if statDataCl[ver].empty: 
                print 'No retreivals found for {}. Exiting.......'.format(ver)
                sys.exit()


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

        self.sza              = OrderedDict()
        self.saa              = OrderedDict()

        for j, ver in enumerate(statDataCl):
            self.rPrfVMR[ver]  = np.asarray(statDataCl[ver].rprfs[statDataCl[ver].PrimaryGas]) * sclfct
            self.rPrfMol[ver]  = np.asarray(statDataCl[ver].rprfs[statDataCl[ver].PrimaryGas]  * np.asarray(statDataCl[ver].rprfs['AIRMASS']))
            self.dates[ver]    = np.asarray(statDataCl[ver].rprfs['date'])
            self.alt[ver]      = np.asarray(statDataCl[ver].rprfs['Z'][0,:])
            self.Airmass[ver]  = np.asarray(statDataCl[ver].rprfs['AIRMASS'])
            self.waterVMR[ver] = np.asarray(statDataCl[ver].aprfs['H2O']) * sclfct
            self.waterMol[ver] = np.asarray(statDataCl[ver].aprfs['H2O'] * self.Airmass[ver] )
            self.totClmn[ver]  = np.sum(self.rPrfMol[ver],axis=1)
            #TCdryAir[ver] = np.sum(Airmass[ver],axis=1) - np.sum(waterMol[ver],axis=1)
            #TCdry[ver]    = (totClmn[ver] / TCdryAir[ver]) * sclfct

            self.aPrfVMR[ver]  = np.asarray(statDataCl[ver].aprfs[statDataCl[ver].PrimaryGas]) * sclfct
            self.aPrfMol[ver]  = np.asarray(statDataCl[ver].aprfs[statDataCl[ver].PrimaryGas]  * np.asarray(statDataCl[ver].aprfs['AIRMASS']))

            #----------------------------------------
            # This is the mixing ratio for DRY AIR!!!
            #----------------------------------------
            #rPrfDry[ver] = np.asarray(statDataCl[ver].rprfs[statDataCl[ver].PrimaryGas]) / (1.0 - waterVMR[ver]) * sclfct       
            
            #----------------------------------
            # Read Summary data (For filtering)
            #----------------------------------
            statDataCl[ver].readsummary()
            self.rms[ver]     = np.asarray(statDataCl[ver].summary[statDataCl[ver].PrimaryGas+'_FITRMS'])

            #----------------------------------
            # Read readPbp (observed, fitted, and difference spectra)
            #----------------------------------
            mw    = [str(int(x)) for x in statDataCl[ver].ctl['band']]     
            numMW = len(mw)
            statDataCl[ver].readPbp()

            self.sza[ver]   = statDataCl[ver].pbp['sza']
            self.saa[ver]   = statDataCl[ver].pbp['saa']

            #----------------------------------
            # Read Spectra for each gas
            #----------------------------------
            statDataCl[ver].readSpectra(statDataCl[ver].gasList)

            #-------------------- 
            # Call to filter data
            #--------------------
            if fltrFlg: statDataCl[ver].fltrData(statDataCl[ver].PrimaryGas,mxrms=maxRMS,rmsFlg=rmsFlg, minDOF=minDOF,  dofFlg=dofFlg, tcFlg=tcFlg, pcFlg=pcFlg , cnvrgFlg=cnvrgFlg)
            else:       statDataCl[ver].inds = np.array([]) 

                #--------------------------------------------
            # Read Error data to get AVK and profile DOFs
            #-------------------------------------------------
            # Determine if AVK has been created via sfit4 core
            # code or via python error analysis
            #-------------------------------------------------     
            if errFlg:   # Read AVK from error output
                #statDataCl[ver].readError(totFlg=False,sysFlg=False,randFlg=False,vmrFlg=True,avkFlg=True,KbFlg=False)
                statDataCl[ver].readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False)

                npnts           = np.shape(statDataCl[ver].error['Total_Random_Error'])[0]
                nlvls           = np.shape(self.alt[ver])[0]
                
                #-------------------------------------------------
                #Error profiles constructed with the diagonal elements of the covariances matrices
                #-------------------------------------------------
                self.rand_err[ver]   = np.zeros((npnts,nlvls))
                self.sys_err[ver]    = np.zeros((npnts,nlvls))

                for i in range(npnts):
                    self.rand_err[ver][i,:] = np.diag(statDataCl[ver].error['Total_Random_Error'][i][:,:])
                    self.sys_err[ver][i,:]  = np.diag(statDataCl[ver].error['Total_Systematic_Error'][i][:,:])

                self.tot_err[ver]     = np.sqrt(self.rand_err[ver] + self.sys_err[ver])            
                self.sys_err[ver]     = np.sqrt(self.sys_err[ver])
                self.rand_err[ver]    = np.sqrt(self.rand_err[ver])

                self.sys_errvmr[ver]  = self.sys_err[ver]/ np.asarray(statDataCl[ver].rprfs['AIRMASS']) * sclfct
                self.rand_errvmr[ver] = self.rand_err[ver]/ np.asarray(statDataCl[ver].rprfs['AIRMASS']) * sclfct
                self.tot_errvmr[ver]  = np.sqrt(self.rand_errvmr[ver]**2 + self.sys_errvmr[ver]**2) 

                #-------------------------------------------------
                #Error profiles of components
                #-------------------------------------------------
                self.rand_cmpnts[ver] = statDataCl[ver].randErrDiag
                self.sys_cmpnts[ver]  = statDataCl[ver].sysErrDiag

                self.rand_cmpnts_vmr[ver] = statDataCl[ver].randErrDiag 
                self.sys_cmpnts_vmr[ver]  = statDataCl[ver].sysErrDiag

                for k in self.sys_cmpnts_vmr[ver]:
                    self.sys_cmpnts_vmr[ver][k] = (np.sqrt(self.sys_cmpnts_vmr[ver][k])/ np.asarray(statDataCl[ver].rprfs['AIRMASS']))*sclfct
                  
                for k in self.rand_cmpnts_vmr[ver]:
                    self.rand_cmpnts_vmr[ver][k] = (np.sqrt(self.rand_cmpnts_vmr[ver][k])/ np.asarray(statDataCl[ver].rprfs['AIRMASS']))*sclfct

                #-------------------------------------------------
                #Error in the summary output. Get Total Errors 
                #-------------------------------------------------
                self.tot_rnd[ver]        = np.array(statDataCl[ver].error['Total random uncertainty'])
                self.tot_sys[ver]        = np.array(statDataCl[ver].error['Total systematic uncertainty'])
                self.tot_std[ver]        = np.sqrt(self.tot_rnd[ver]**2 + self.tot_sys[ver]**2)
                
                self.err_summary[ver]    = statDataCl[ver].error

                #---------------------
                # Get averaging kernel
                #---------------------   
                self.avkSCF[ver]  = np.delete(np.asarray(statDataCl[ver].error['AVK_scale_factor']),statDataCl[ver].inds,axis=0)
                self.avkVMR[ver]  = np.delete(np.asarray(statDataCl[ver].error['AVK_vmr']),statDataCl[ver].inds,axis=0)   
                self.dofs            = np.diagonal(self.avkSCF[ver],axis1=1,axis2=2)
                self.avkSCFav[ver]  = np.mean(self.avkSCF[ver],axis=0)    
                self.avkVMRav[ver]  = np.mean(self.avkVMR[ver],axis=0)                   
                    
                self.dofsAvg[ver]    = np.diag(self.avkSCFav[ver])
                self.dofsAvg_cs[ver] = np.cumsum(np.diag(self.avkSCFav[ver])[::-1])[::-1]            
            
            else:        # Read AVK from sfit4 output (only contains scaled AVK)
                avkSCFi = []
                for d in statDataCl[ver].dirLst:
                    lines  = dc.tryopen( d + statDataCl[ver].ctl['file.out.ak_matrix'][0])
                    if not lines: continue
                    avkSCFi.append(np.array( [ [ float(x) for x in line.strip().split() ] for line in lines[2:] ] ))
                    
                if not statDataCl[ver].readPrfFlgApr[statDataCl[ver].PrimaryGas]: statDataCl[ver].readprfs([statDataCl[ver].PrimaryGas],retapFlg=0)   # Apriori Profiles
                          
                self.avkSCF[ver]  = np.asarray(avkSCFi)
                nobs            = np.shape(self.avkSCF[ver])[0]
                n_layer         = np.shape(self.avkSCF[ver])[1]
                self.avkVMR[ver]  = np.zeros((nobs,n_layer,n_layer))
        
                for obs in range(0,nobs):
                    Iapriori        = np.zeros((n_layer,n_layer))
                    IaprioriInv     = np.zeros((n_layer,n_layer))
                    np.fill_diagonal(Iapriori,statDataCl[ver].aprfs[statDataCl[ver].PrimaryGas.upper()][obs])
                    np.fill_diagonal(IaprioriInv, 1.0 / (statDataCl[ver].aprfs[statDataCl[ver].PrimaryGas.upper()][obs]))
                    self.avkVMR[ver][obs,:,:] = np.dot(np.dot(Iapriori,np.squeeze(self.avkSCF[ver][obs,:,:])),IaprioriInv)       
                    
                self.avkSCF[ver]     = np.delete(self.avkSCF[ver],statDataCl[ver].inds,axis=0)
                self.avkVMR[ver]     = np.delete(self.avkVMR[ver],statDataCl[ver].inds,axis=0)
                self.dofs               = np.diagonal(self.avkSCF[ver],axis1=1,axis2=2)
                self.avkSCFav[ver]     = np.mean(self.avkSCF[ver],axis=0)
                self.avkVMRav[ver]     = np.mean(self.avkVMR[ver],axis=0)
                
                self.dofsAvg[ver]    = np.diag(self.avkSCFav[ver])
                self.dofsAvg_cs[ver] = np.cumsum(np.diag(self.avkSCFav[ver])[::-1])[::-1]


            #--------------------------------------
            # Remove retrieval data based on filter
            #--------------------------------------
            nfltr              = len(statDataCl[ver].inds)
            self.rms[ver]      = np.delete(self.rms[ver],statDataCl[ver].inds)
            self.sza[ver]      = np.delete(self.sza[ver],statDataCl[ver].inds)
            self.saa[ver]      = np.delete(self.saa[ver],statDataCl[ver].inds)
            ntot               = len(self.rms[ver])
            self.dates[ver]    = np.delete(self.dates[ver],statDataCl[ver].inds)
            self.totClmn[ver]  = np.delete(self.totClmn[ver],statDataCl[ver].inds)
            self.rPrfVMR[ver]  = np.delete(self.rPrfVMR[ver],statDataCl[ver].inds,axis=0)   
            self.rPrfMol[ver]  = np.delete(self.rPrfMol[ver],statDataCl[ver].inds,axis=0)
            self.waterVMR[ver] = np.delete(self.waterVMR[ver],statDataCl[ver].inds,axis=0)   
            self.waterMol[ver] = np.delete(self.waterMol[ver],statDataCl[ver].inds,axis=0)   

            #rPrfDry[ver] = np.delete(rPrfDry[ver],statDataCl[ver].inds,axis=0)   
            self.Airmass[ver] = np.delete(self.Airmass[ver],statDataCl[ver].inds,axis=0)
            #TCdry[ver]   = np.delete(TCdry[ver],statDataCl[ver].inds)

            self.aPrfVMR[ver] = np.delete(self.aPrfVMR[ver],statDataCl[ver].inds,axis=0)   
            self.aPrfMol[ver] = np.delete(self.aPrfMol[ver],statDataCl[ver].inds,axis=0)

            if errFlg:
                self.rand_err[ver] = np.delete(self.rand_err[ver],statDataCl[ver].inds,axis=0)
                self.sys_err[ver]  = np.delete(self.sys_err[ver],statDataCl[ver].inds,axis=0)  
                self.tot_err[ver]  = np.delete(self.tot_err[ver],statDataCl[ver].inds,axis=0)

                self.rand_errvmr[ver] = np.delete(self.rand_errvmr[ver],statDataCl[ver].inds,axis=0)
                self.sys_errvmr[ver]  = np.delete(self.sys_errvmr[ver],statDataCl[ver].inds,axis=0)  
                self.tot_errvmr[ver]  = np.delete(self.tot_errvmr[ver],statDataCl[ver].inds,axis=0)

                self.tot_rnd[ver] = np.delete(self.tot_rnd[ver],statDataCl[ver].inds)
                self.tot_sys[ver] = np.delete(self.tot_sys,statDataCl[ver].inds)
                self.tot_std[ver] = np.delete(self.tot_std,statDataCl[ver].inds)

                
                for k in self.sys_cmpnts[ver]:
                    self.sys_cmpnts[ver][k]      = np.delete(self.sys_cmpnts[ver][k],statDataCl[ver].inds,axis=0)
                    self.sys_cmpnts_vmr[ver][k]  = np.delete(self.sys_cmpnts_vmr[ver][k],statDataCl[ver].inds,axis=0)
                    
                for k in self.rand_cmpnts[ver]:
                    self.rand_cmpnts[ver][k]     = np.delete(self.rand_cmpnts[ver][k],statDataCl[ver].inds,axis=0)
                    self.rand_cmpnts_vmr[ver][k] = np.delete(self.rand_cmpnts_vmr[ver][k],statDataCl[ver].inds,axis=0)

                for k in self.err_summary[ver]:
                    self.err_summary[ver][k]     = np.delete(self.err_summary[ver][k],statDataCl[ver].inds, axis=0)



    


#------------------------------------------------------------------------------------------------------------------------------        
class Pltsonde(SondeClass, dc.ReadOutputData):

    def __init__(self, dataDir, ctlF, dataDirSonde, loc, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1,outFname='', fleoutFlg=False, version=''):\
        
        primGas = ''
        if outFname: self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.ver = version
        self.loc = loc


        #---------------
        # ReadOutputsSonde
        #---------------
        ReadOutputsSonde.__init__(self,dataDirSonde,loc,iyear,imnth,iday,fyear,fmnth,fday,incr, fleoutFlg)    

        #------------
        # ReadOutputData
        #------------
        dc.ReadOutputData.__init__(self,dataDir,primGas,ctlF,iyear,imnth,iday,fyear,fmnth,fday,incr) 

   
    def closeFig(self):
        self.pdfsav.close()

    def plt(self,fltr=False,minSZA=0.0,maxSZA=80.0,maxRMS=10.0,minDOF=1.0,maxCHI=2.0,minTC=1.0E15,maxTC=1.0E16,dofFlg=False,rmsFlg=True,tcFlg=True,mnthFltr=[1,2,3,4,5,6,7,8,9,10,11,12],
               pcFlg=True,cnvrgFlg=True,allGas=False,sclfct=1.0,sclname='ppv',pltStats=True,szaFlg=False,errFlg=False,chiFlg=False,tcMMflg=False,mnthFltFlg=False, pCols=False, smthFlg=False, doiplt=False, loc='', latFTS=False, lonFTS=False,
               altCG=False, diffT= False, dfValue=False):
        
        #-------------------------------------------------------
        #READING ALL THE SONDE FILES
        #-------------------------------------------------------
        self.readfilessonde()

        print 'Total Number of sonde observations = ' +str(len(self.dirLstsonde))

        self.datefts = []
        self.dirDateTime.sort()

        for dd in self.dirDateTime:
            self.datefts.append(dt.date( int(dd.year), int(dd.month), int(dd.day) ) )

        self.datefts = np.asarray(self.datefts)
        self.sonde['date'] = np.asarray(self.sonde['date'])


        doy_fts   = dc.toYearFraction(self.datefts)
        doy_sonde = dc.toYearFraction(self.sonde['date'])

        #-------------------------------------------------------
        #FINDING COINCIDENT DATES
        #-------------------------------------------------------
        intrsctVals = np.intersect1d(doy_fts, doy_sonde, assume_unique=False)
    
        inds1       = np.nonzero( np.in1d( doy_sonde, intrsctVals, assume_unique=False ) )[0]
        inds2       = np.nonzero( np.in1d( doy_fts, intrsctVals, assume_unique=False ) )[0]

        print 'Total Number of coincident dates between FTS and sondes = ' +str(len(intrsctVals))+'\n'

        #-------------------------------------------------------
        #CREATING NEW DIR LIST FOR COINCIDENT DATES
        #-------------------------------------------------------
        self.dirLstsonde = np.array(self.dirLstsonde)
        self.dirLst      = np.array(self.dirLst)

        self.dirLstsonde = self.dirLstsonde[inds1]
        self.dirLst      = self.dirLst[inds2]

        #-------------------------------------------------------
        #READING THE SONDE FILES FOR COINCIDENT DATES
        #-------------------------------------------------------   

        self.readfilessonde()
        self.sondePrf()

        #-------------------------------------------------------
        #Test sonde and calculate center of gravity
        #-------------------------------------------------------
        #sondelat    = np.asarray(self.sonde['latFleout_a'])
        #sondelon    = np.asarray(self.sonde['lonFleout_a'])

        sondealt  = np.asarray(self.sonde['Alt_km_a'])
        sondeprf  = np.asarray(self.sonde['H2Omr_ppmv_a'])

        opp   = []
        for ii, r in enumerate(sondeprf):

            if self.loc.upper() == 'MLO': inds = np.where( sondealt[ii] > 3.5)[0]
            if self.loc.upper() == 'FL0': inds = np.where( sondealt[ii] > 1.5)[0]

            opp.append(np.average(sondealt[ii][inds], weights= r[inds]))

        opp = np.asarray(opp)

        print 'Center of gravity = {0:.2f} +/- {1:.2f}'.format(np.mean(opp), np.std(opp))

        
        # fig, ax1       = plt.subplots()

        # for i in range(sondelat.shape[0]):
    
        #     ax1.plot(sondelon[i][:], sondelat[i][:], linewidth=0.75)
        
        # ax1.set_ylabel('Latitude')
        # ax1.set_xlabel('Longitude')        
        # ax1.grid(True,which='both')
        # ax1.set_title('Location', fontsize=12)
        # ax1.tick_params(labelsize=12)
        # ax1.set_xlim(-106, -102)
        # ax1.set_ylim(38, 42)

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig,dpi=200)
        #     plt.close(fig)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #     sys.exit()


        #-------------------------------------------------------
        #READING THE FTS DATA ONLY FOR COINCIDENT DATES
        #-------------------------------------------------------
        aprPrf       = {}
        rPrf         = {}
        localGasList = [self.PrimaryGas]
        
        #-------------------------------------------------------
        # Get profile, summary for Primary gas.... for filtering
        #-------------------------------------------------------
        if not self.readPrfFlgRet[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=1)   # Retrieved Profiles
        if self.empty: return False
        if not self.readPrfFlgApr[self.PrimaryGas]: self.readprfs([self.PrimaryGas],retapFlg=0)   # Apriori Profiles
        
        aprPrf      = np.asarray(self.aprfs[self.PrimaryGas]) * sclfct
        rPrf        = np.asarray(self.rprfs[self.PrimaryGas]) * sclfct
        rPrfMol     = np.asarray(self.rprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        aprPrffMol  = np.asarray(self.aprfs[self.PrimaryGas]) * np.asarray(self.rprfs['AIRMASS'])
        Airmass     = np.asarray(self.rprfs['AIRMASS'])

        if len(self.dirLst) > 1: 
            dates                   = self.rprfs['date']

        alt = np.asarray(self.rprfs['Z'][0,:])

        
        if not self.readsummaryFlg: self.readsummary()                                      # Summary File info
        rms     = np.asarray(self.summary[self.PrimaryGas+'_FITRMS'])
        dofs    = np.asarray(self.summary[self.PrimaryGas+'_DOFS_TRG'])
        totClmn = np.asarray(self.summary[self.PrimaryGas.upper()+'_RetColmn'])
        
        if not self.readPbpFlg: self.readPbp()                                              # Pbp file info
        sza   = self.pbp['sza']
        saa   = self.pbp['saa']
                 
        if errFlg:                                    # Error info
            
            if not all((self.readErrorFlg['totFlg'],self.readErrorFlg['sysFlg'],self.readErrorFlg['randFlg'])):
                self.readError(totFlg=True,sysFlg=True,randFlg=True,vmrFlg=False,avkFlg=True,KbFlg=False) 
            
            npnts    = np.shape(self.error['Total_Random_Error'])[0]
            nlvls    = np.shape(alt)[0]

            #-------------------------------------------------
            #Error profiles constructed with the diagonal elements of the covariances matrices
            #-------------------------------------------------
            rand_err = np.zeros((npnts,nlvls))
            sys_err  = np.zeros((npnts,nlvls))
            
            for i in range(npnts):
                rand_err[i,:] = np.diag(self.error['Total_Random_Error'][i][:,:])
                sys_err[i,:]  = np.diag(self.error['Total_Systematic_Error'][i][:,:])
                
            tot_err  = np.sqrt(rand_err + sys_err ) #* sclfct         
            rand_err = np.sqrt(rand_err) #* sclfct
            sys_err  = np.sqrt(sys_err)  #* sclfct

            rand_errvmr  = (rand_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            sys_errvmr   = (sys_err/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
            tot_errvmr   = np.sqrt(rand_errvmr**2 + sys_errvmr**2)
            
            #-------------------------------------------------
            #Error profiles of components
            #-------------------------------------------------
            rand_cmpnts      = self.randErrDiag 
            sys_cmpnts       = self.sysErrDiag

            rand_cmpnts_vmr  = self.randErrDiag  
            sys_cmpnts_vmr   = self.sysErrDiag

            for k in sys_cmpnts_vmr:
                 sys_cmpnts_vmr[k] = (np.sqrt(sys_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct
                 
            for k in rand_cmpnts_vmr:
                 rand_cmpnts_vmr[k] = (np.sqrt(rand_cmpnts_vmr[k])/ np.asarray(self.rprfs['AIRMASS'][0:npnts]))*sclfct

            #-------------------------------------------------
            #Error in the summary output. Get Total Errors 
            #-------------------------------------------------
            rand_errTC  = np.asarray(self.error['Total random uncertainty'])
            sys_errTC   = np.asarray(self.error['Total systematic uncertainty'])
            tot_errTC   = np.sqrt(rand_errTC**2 + sys_errTC**2)  
    
        #--------------------
        # Call to filter data
        #--------------------
        if fltr: self.fltrData(self.PrimaryGas,mxrms=maxRMS,minsza=minSZA,mxsza=maxSZA,minDOF=minDOF,maxCHI=maxCHI,minTC=minTC,maxTC=maxTC,mnthFltr=mnthFltr,
                               dofFlg=dofFlg,rmsFlg=rmsFlg,tcFlg=tcFlg,pcFlg=pcFlg,szaFlg=szaFlg,cnvrgFlg=cnvrgFlg,chiFlg=chiFlg,tcMinMaxFlg=tcMMflg,mnthFltFlg=mnthFltFlg)
        else:    self.inds = np.array([]) 
        
        if self.empty: return False
        
        #--------------------------------------------------------
        # Get profile data for other gases if allGas flag is True
        #--------------------------------------------------------
        if allGas:
            nonPgas = (gas for gas in self.gasList if gas != self.PrimaryGas)
            for gas in nonPgas:
                if not self.readPrfFlgRet[gas.upper()]: self.readprfs([gas.upper()],retapFlg=1)   # Retrieved Profiles
                if not self.readPrfFlgApr[gas.upper()]: self.readprfs([gas.upper()],retapFlg=0)   # Apriori Profiles
                aprPrf[gas.upper()] = np.asarray(self.aprfs[gas.upper()]) * sclfct
                rPrf[gas.upper()]   = np.asarray(self.rprfs[gas.upper()]) * sclfct
                localGasList.append(gas)
          
        #----------------------------
        # Remove inds based on filter
        #----------------------------
        nfltr   = len(self.inds)
        rms     = np.delete(rms, self.inds)
        ntot    = len(rms)
        sza     = np.delete(sza, self.inds)
        saa     = np.delete(saa, self.inds)
        if len(self.dirLst) > 1:
            dates   = np.delete(dates, self.inds)
            datefts = []
            for dd in dates:
                datefts.append(dt.date( int(dd.year), int(dd.month), int(dd.day) ) )
            datefts = np.asarray(datefts)
        dofs       = np.delete(dofs,self.inds)
        totClmn    = np.delete(totClmn,self.inds)
        rPrfMol    = np.delete(rPrfMol,self.inds,axis=0)
        aprPrffMol = np.delete(aprPrffMol,self.inds,axis=0)
        Airmass    = np.delete(Airmass,self.inds,axis=0)
        rPrf       = np.delete(rPrf,self.inds,axis=0)
        aprPrf     = np.delete(aprPrf,self.inds,axis=0)

            
        if errFlg:
            rand_err    = np.delete(rand_err,self.inds,axis=0)
            sys_err     = np.delete(sys_err,self.inds,axis=0)  
            tot_err     = np.delete(tot_err,self.inds,axis=0)

            rand_errvmr = np.delete(rand_errvmr,self.inds,axis=0)
            sys_errvmr  = np.delete(sys_errvmr,self.inds,axis=0)  
            tot_errvmr  = np.delete(tot_errvmr,self.inds,axis=0)

            rand_errTC  = np.delete(rand_errTC,self.inds)
            sys_errTC   = np.delete(sys_errTC,self.inds)
            tot_errTC   = np.delete(tot_errTC,self.inds)

            
            for k in sys_cmpnts:
                sys_cmpnts[k]      = np.delete(sys_cmpnts[k],self.inds,axis=0)
                sys_cmpnts_vmr[k]  = np.delete(sys_cmpnts_vmr[k],self.inds,axis=0)
                
            for k in rand_cmpnts:
                rand_cmpnts[k]     = np.delete(rand_cmpnts[k],self.inds,axis=0)
                rand_cmpnts_vmr[k] = np.delete(rand_cmpnts_vmr[k],self.inds,axis=0)

            #---------------------
            # Get averaging kernel
            #---------------------   
            avkSCF      = np.delete(np.asarray(self.error['AVK_scale_factor']),self.inds,axis=0)
            avkVMR      = np.delete(np.asarray(self.error['AVK_vmr']), self.inds,axis=0)   
            avkSCFav    = np.mean(avkSCF,axis=0)    
            avkVMRav    = np.mean(avkVMR,axis=0)

        #--------------------------------------------------------
        # Get location of FTS based on altitude profile weighted from sone, azimuth (bearing) and SZA
        #--------------------------------------------------------
        if altCG:

            sza    = np.asarray(sza)
            saa    = np.asarray(saa)
            szarad = [radians(s) for s in sza]

            distFTS = [altCG*tan(s) for s in   szarad]

            latFTS2  = []
            lonFTS2  = []

            for i, a in enumerate(saa):
                lat2, lon2 = mf.destination(latFTS, lonFTS, distFTS[i], a, radius=6371.008771415)
                latFTS2.append(lat2)
                lonFTS2.append(lon2)

            latFTS2 = np.asarray(latFTS2)
            lonFTS2 = np.asarray(lonFTS2)

        Dist  = [mf.haversine(lonFTS, latFTS, lo, la) for (lo, la) in zip(lonFTS2, latFTS2) ]

        print 'Mean horizontal positions of the center of gravity = {0:.2f} $\pm$ {1:.2f}'.format(np.mean(Dist), np.std(Dist))


        #----------------------------
        # Determine if multiple years
        #----------------------------
        if len(self.dirLst) > 1:
            years = [ singDate.year for singDate in dates]      # Find years for all date entries
            if len(list(set(years))) > 1: yrsFlg = True         # Determine all unique years
            else:                         yrsFlg = False

        #---------------------------------
        # Re - Defining variables (FTIR) using Max/Min altitude
        #---------------------------------
        maxalt      = 25. #km
        indsalt     = np.where(alt<= maxalt)[0]
        alt         = alt[indsalt]
        rPrf        = rPrf[:, indsalt]
        aprPrf      = aprPrf[:, indsalt]
        tot_err     = tot_err[:, indsalt]
        rand_errvmr = rand_errvmr[:, indsalt]
        sys_errvmr  = sys_errvmr[:, indsalt]
        rPrfMol     = rPrfMol[:, indsalt]
        tot_errvmr  = tot_errvmr[:, indsalt]
        aprPrffMol  = aprPrffMol[:, indsalt]
        avkSCFav    = avkSCFav[indsalt[0]:, indsalt[0]:]
        avkVMRav    = avkVMRav[indsalt[0]:, indsalt[0]:]
        avkSCF      = avkSCF[:, indsalt[0]:, indsalt[0]: ]
        avkVMR      = avkVMR[:, indsalt[0]:, indsalt[0]: ]
        Airmass     = Airmass[:, indsalt]

        #---------------------------
        # Define Delta time,  plot Profiles, and calculate columns/VRM(weighted) for both FTS and sonde
        #---------------------------

        prfmean_day            = {}
        prferr_day             = {}
        prfSTD_day             = {}
        aPrf_day               = {}
        Airmassmean            = {}
        NobsFTS                = {}

        latmean                = {}
        lonmean                = {}

        avkSCFday              = {}
        avkSCFVMRday           = {}
        rPrfMolmean_day        = {}
        aPrfMol_day            = {}
        prferrMol_day          = {}

        sondePrf               = {}
        sondePrfsd             = {}
        sondeairmass           = {}
        sondePrfMol            = {}
        sondePrfMolsd          = {}
        sondePrflat            = {}
        sondePrflon            = {}

        Prfsonde_interp        = {}
        Prfsonde_sd_interp     = {}
        PrfsondeMol_interp     = {}
        PrfsondeMol_sd_interp  = {}
        sondeairmass_a_interp  = {}
        Prfsondelat_interp     = {}
        Prfsondelon_interp     = {}
        sondealt               = {}

        PrfDiff                = {}
        PrfDiffApr             = {}
        PrfDiffRel             = {}
        PrfDiffAprRel          = {}
        PrfDist                = {}

        FTSvmrP                = {}
        FTSsumP                = {}
        FTSvmrsdP              = {}
        FTSsumsdP              = {}
        FTSdates               = {}

        FTSAprvmrP             = {}
        FTSAprsumP             = {}

        sondevmrP              = {}
        sondevmrsdP            = {}
        sondesumP              = {}
        sondesumsdP            = {}
        sondedates             = {}

        sondelatP              = {}
        sondelonP              = {}
        distP                  = {}


        #---------------------------------
        #
        #---------------------------------

        for t, df in enumerate(diffT):

            for d, da in enumerate(self.sonde['date']):

                #---------------------------------
                #Index for same days - FTIR
                #---------------------------------
                dtsonde   = self.sonde['dt'][d]
                deltaT    = dtsonde -   dates 
                deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.
                
                #inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= df])

                if t == 0: inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= df])
                else: inds = np.array([x for x, dff in enumerate(deltaTmin) if( (abs(dff) >= diffT[t-1]) and ( abs(dff) < df))])
                ###inds = np.array([x for x, df in enumerate(datefts) if df == da])

                if len(inds) >= 1:
                    
                    #prfmean_day   = np.mean(rPrf[inds],axis=0)
                    ##WEIGTHED MEAN
                    prfmean_day.setdefault(df, []).append(np.average(rPrf[inds],axis=0, weights=tot_errvmr[inds]) ) # ==> prfmean_day   = np.sum(rPrf[inds]/(tot_errvmr[inds]**2), axis=0 )/ np.sum(1/(tot_err[inds]**2), axis=0); http://www.physics.umd.edu/courses/Phys261/F06/ErrorPropagation.pdf
                    #prfmean_day.setdefault(df, []).append(np.sum(rPrf[inds]/(tot_errvmr[inds]**2), axis=0 )/ np.sum(1/(tot_errvmr[inds]**2), axis=0))
                    #prfmean_day.setdefault(df, []).append(np.mean(rPrf[inds],axis=0 ))
                    ## uncertainty in the weighted mean can be evaluated using standard error propagation  
                    ##WEIGTHED UNCERTAINTY (DOUBLE CHECK THIS)
                    #prferr_day.setdefault(df, []).append(np.sqrt (1.0/np.sum(1.0/tot_errvmr[inds]**2, axis=0)))
                    ##UNCERTAINTY USING standard error propagation
                    prferr_day.setdefault(df, []).append(np.sqrt (np.sum(tot_errvmr[inds]**2, axis=0))) 

                    prfSTD_day.setdefault(df, []).append(np.std(rPrf[inds], axis=0))
                    
                    aPrf_day.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0))

                    Airmassmean.setdefault(df, []).append(np.mean(Airmass[inds],axis=0))
                    avkSCFday.setdefault(df, []).append(np.mean(avkSCF[inds],axis=0))
                    avkSCFVMRday.setdefault(df, []).append(np.mean(avkVMR[inds],axis=0))
                    
                    #rPrfMolmean_day   = np.mean(rPrfMol[inds],axis=0)
                    rPrfMolmean_day.setdefault(df, []).append(np.average(rPrfMol[inds],axis=0, weights=tot_err[inds]))
                    prferrMol_day.setdefault(df, []).append(np.sqrt (1.0/np.sum(1.0/tot_err[inds]**2, axis=0)))

                    aPrfMol_day.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0))

                    NobsFTS.setdefault(df, []).append(len(inds))

                    latmean.setdefault(df, []).append(np.mean(latFTS2[inds]))
                    lonmean.setdefault(df, []).append(np.mean(lonFTS2[inds]))

                    sondeH2OPrf_a      = np.asarray(self.sonde['H2Omr_ppmv_a'][d], dtype=np.float32)
                    sondeH2OPrfsd_a    = np.asarray(self.sonde['H2Osd_ppmv_a'][d], dtype=np.float32)
                    sondeH2OPrfMol_a   = np.asarray(self.sonde['H2OMol_a'][d], dtype=np.float32)
                    sondeH2OPrfMolsd_a = np.asarray(self.sonde['H2OMolsd_a'][d], dtype=np.float32)

                    sondealt_a          = np.asarray(self.sonde['Alt_km_a'][d], dtype=np.float32)
                    sondeairmass_a      = np.asarray(self.sonde['airmass_a'][d], dtype=np.float32)

                    if self.fleoutFlg: 
                        sondelat_a            = np.asarray(self.sonde['latFleout_a'][d])
                        sondelon_a            = np.asarray(self.sonde['lonFleout_a'][d])

                else: continue
                
                #ind1 = mf.nearestind(pcol[0], sondealt_a)
                #ind2 = mf.nearestind(pcol[1], sondealt_a) 

                sondePrf.setdefault(df, []).append(sondeH2OPrf_a)
                sondePrfsd.setdefault(df, []).append(sondeH2OPrfsd_a)
                sondeairmass.setdefault(df, []).append(sondeairmass_a)
                sondePrfMol.setdefault(df, []).append(sondeH2OPrfMol_a)
                sondePrfMolsd.setdefault(df, []).append(sondeH2OPrfMolsd_a)

                if self.fleoutFlg: 
                    sondePrflat.setdefault(df, []).append(sondelat_a)
                    sondePrflon.setdefault(df, []).append(sondelon_a)

                sondealt_a  = map(lambda x: float(x),sondealt_a)
                sondealt_a  = np.asarray(sondealt_a)

                sondedates.setdefault(df, []).append(da)
                sondealt.setdefault(df, []).append(sondealt_a)

                day = np.asarray([dt.date(d.year,d.month,d.day) for d in dates[inds]])
                uniqueDay = list(set(day))          # Find a list of unique days   
                FTSdates.setdefault(df, []).append(uniqueDay)


                Prfsonde_interp_i           = interpolate.interp1d(sondealt_a, sondeH2OPrf_a, axis=0, fill_value=sondeH2OPrf_a[0], bounds_error=False, kind='linear')(alt)
                #Prfsonde_interp_i           = interpolate.interp1d(sondealt_a, sondeH2OPrf_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                #Prfsonde_sd_interp_i        = interpolate.interp1d(sondealt_a, sondeH2OPrfsd_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt)
                Prfsonde_sd_interp_i        = interpolate.interp1d(sondealt_a, sondeH2OPrfsd_a, axis=0, fill_value=sondeH2OPrfsd_a[0], bounds_error=False, kind='linear')(alt)
                PrfsondeMol_interp_i        = interpolate.interp1d(sondealt_a, sondeH2OPrfMol_a, axis=0, fill_value=sondeH2OPrfsd_a[0], bounds_error=False, kind='linear')(alt)
                PrfsondeMol_sd_interp_i     = interpolate.interp1d(sondealt_a, sondeH2OPrfMolsd_a, axis=0, fill_value=sondeH2OPrfsd_a[0], bounds_error=False, kind='linear')(alt)
                #sondeairmass_a_interp.setdefault(df, []).append(interpolate.interp1d(sondealt_a, sondeairmass_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='linear')(alt))
                sondeairmass_a_interp.setdefault(df, []).append(interpolate.interp1d(sondealt_a, sondeairmass_a, axis=0, fill_value=sondeairmass_a[0], bounds_error=False, kind='linear')(alt))

                
                if self.fleoutFlg:
                    #sondePrflon_interp_i            = interpolate.interp1d(sondealt_a, sondelon_a, axis=0,  bounds_error=False, kind='linear')(alt)
                    #sondePrflat_interp_i            = interpolate.interp1d(sondealt_a, sondelat_a, axis=0,  bounds_error=False, kind='linear')(alt)

                    sondePrflon_interp_i            = interpolate.interp1d(sondealt_a, sondelon_a, axis=0, fill_value=sondelon_a[0], bounds_error=False, kind='linear')(alt)
                    sondePrflat_interp_i            = interpolate.interp1d(sondealt_a, sondelat_a, axis=0, fill_value=sondelat_a[0], bounds_error=False, kind='linear')(alt)

                    Prfsondelon_interp.setdefault(df, []).append(sondePrflon_interp_i)
                    Prfsondelat_interp.setdefault(df, []).append(sondePrflat_interp_i)


                if smthFlg:
                    #---------------------------------
                    # Smoothing using FTIR AK and apriori
                    #---------------------------------
                    Prfsonde_interp.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0)/1e3 + np.dot(np.mean(avkVMR[inds],axis=0), (Prfsonde_interp_i -  np.mean(aprPrf[inds], axis=0)/1e3)))
                    Prfsonde_sd_interp.setdefault(df, []).append(np.mean(aprPrf[inds], axis=0)/1e3 + np.dot(np.mean(avkVMR[inds],axis=0), (Prfsonde_sd_interp_i -  np.mean(aprPrf[inds], axis=0)/1e3)))

                    PrfsondeMol_interp.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0)/1e3 + np.dot(np.mean(avkSCF[inds],axis=0), (PrfsondeMol_interp_i -  np.mean(aprPrffMol[inds],axis=0)/1e3)))
                    PrfsondeMol_sd_interp.setdefault(df, []).append(np.mean(aprPrffMol[inds],axis=0)/1e3 + np.dot(np.mean(avkSCF[inds],axis=0), (PrfsondeMol_sd_interp_i -  np.mean(aprPrffMol[inds],axis=0)/1e3)))

                else:
                    Prfsonde_interp.setdefault(df, []).append(Prfsonde_interp_i)
                    Prfsonde_sd_interp.setdefault(df, []).append(Prfsonde_sd_interp_i)
                    PrfsondeMol_interp.setdefault(df, []).append(PrfsondeMol_interp_i)
                    PrfsondeMol_sd_interp.setdefault(df, []).append(PrfsondeMol_sd_interp_i)
                    #sondeairmass_a_interp.setdefault(df, []).append(interpolate.interp1d(sondealt_a, sondeairmass_a, axis=0, fill_value='extrapolate', bounds_error=False, kind='nearest')(alt))
                    

            prfmean_day[df]           = np.asarray(prfmean_day[df])
            prferr_day[df]            = np.asarray(prferr_day[df])
            prfSTD_day[df]            = np.asarray(prfSTD_day[df])
            aPrf_day[df]              = np.asarray(aPrf_day[df])
            Airmassmean[df]           = np.asarray(Airmassmean[df])
            avkSCFday[df]             = np.asarray(avkSCFday[df])
            rPrfMolmean_day[df]       = np.asarray(rPrfMolmean_day[df])
            prferrMol_day[df]         = np.asarray(prferrMol_day[df])
            aPrfMol_day[df]           = np.asarray(aPrfMol_day[df])
            FTSdates[df]              = np.asarray(FTSdates[df])
            NobsFTS[df]               = np.asarray(NobsFTS[df])
            latmean[df]               = np.asarray(latmean[df])
            lonmean[df]               = np.asarray(lonmean[df])

            Prfsonde_interp[df]       = np.asarray(Prfsonde_interp[df])
            Prfsonde_sd_interp[df]    = np.asarray(Prfsonde_sd_interp[df])
            PrfsondeMol_interp[df]    = np.asarray(PrfsondeMol_interp[df])
            PrfsondeMol_sd_interp[df] = np.asarray(PrfsondeMol_sd_interp[df])
            sondeairmass_a_interp[df] = np.asarray(sondeairmass_a_interp[df])
            sondedates[df]            = np.asarray(sondedates[df])
            sondealt[df]              = np.asarray(sondealt[df])
            PrfDiff[df]               = np.asarray(prfmean_day[df] - Prfsonde_interp[df])
            PrfDiffApr[df]            = np.asarray(aPrf_day[df] - Prfsonde_interp[df])

            PrfDiffRel[df]            = np.true_divide(PrfDiff[df], Prfsonde_interp[df])*100.
            PrfDiffAprRel[df]         = np.true_divide(PrfDiffApr[df], Prfsonde_interp[df])*100.

            if self.fleoutFlg:
                Prfsondelon_interp[df]    = np.asarray(Prfsondelon_interp[df])
                Prfsondelat_interp[df]    = np.asarray(Prfsondelat_interp[df])
                PrfDist[df]  = [ [mf.haversine(np.nanmean(lonmean[df]), np.nanmean(latmean[df]), lo2, la2 ) for (lo2, la2) in  zip(lo,la)] for (lo, la) in zip(Prfsondelon_interp[df], Prfsondelat_interp[df]) ]
                PrfDist[df]    = np.asarray(PrfDist[df])


            #------------------------------------
            #For ssome reason at ~10k some values does not make sense
            #------------------------------------
            #PrfDist[df][PrfDist[df] > 500.] = np.nan
             
      
            sondeairmass_a_i = Airmassmean[df]

            #------------------------------------
            #
            #------------------------------------
            # for pn, pcol in enumerate(pCols):

            #     inds = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

            #     sondevmrP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(Prfsonde_interp[df][:,inds]), axis=1,  weights=sondeairmass_a_i[:,inds])
            #     sondevmrsdP[str(df)+'_'+str(pn)]  = np.ma.average(np.ma.masked_invalid(Prfsonde_sd_interp[df][:, inds]), axis=1, weights=sondeairmass_a_i[:,inds])
            #     sondesumP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(Prfsonde_interp[df][:,inds])*np.ma.masked_invalid(sondeairmass_a_i[:,inds])*1e-6, axis=1)
            #     sondesumsdP[str(df)+'_'+str(pn)]  = np.nansum(np.ma.masked_invalid(Prfsonde_sd_interp[df][:,inds])*np.ma.masked_invalid(sondeairmass_a_i[:,inds])*1e-6, axis=1)
                
            #     FTSvmrP[str(df)+'_'+str(pn)]      = np.ma.average(np.ma.masked_invalid(prfmean_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSvmrsdP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(prferr_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSAprvmrP[str(df)+'_'+str(pn)]   = np.ma.average(np.ma.masked_invalid(aPrf_day[df][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,inds]))
            #     FTSsumP[str(df)+'_'+str(pn)]      = np.nansum(np.ma.masked_invalid(rPrfMolmean_day[df][:,inds]), axis=1)
            #     FTSsumsdP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(prferrMol_day[df][:,inds]), axis=1)
            #     FTSAprsumP[str(df)+'_'+str(pn)]   = np.nansum(np.ma.masked_invalid(aPrfMol_day[df][:,inds]), axis=1)

            #     sondelatP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[df][:,inds]), axis=1)
            #     sondelonP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[df][:,inds]), axis=1)

            #     #distP[str(df)+'_'+str(pn)]        = [mf.haversine(np.nanmean(lonmean[df]), np.mean(latmean[df]), lo, la ) for (lo, la) in  zip(sondelonP[str(df)+'_'+str(pn)], sondelatP[str(df)+'_'+str(pn)])]
            #     distP[str(df)+'_'+str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[df][:,inds]), axis=1)
        
        
        #--------------------------------
        # Profiles as a function of Month
        #--------------------------------
        month          = np.array([d.month for d in self.sonde['date']])
        fig, ax1       = plt.subplots()
        cm             = plt.get_cmap('jet')
        cNorm          = colors.Normalize( vmin=1, vmax=12 )
        scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )
        
        scalarMap.set_array(month)
        
        ax1.set_color_cycle( [scalarMap.to_rgba(x) for x in month] )
        
        for i in range(len(month)):
            ax1.plot(np.asarray(self.sonde['H2Omr_ppmv_a'][i])/1e3, np.asarray(self.sonde['Alt_km_a'][i]), linewidth=0.75)
        
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=15)        
        ax1.grid(True,which='both')
        ax1.set_ylim(1,12)
        ax1.set_title('CFH profiles (2010 - 2016)', fontsize=16)
        ax1.tick_params(labelsize=14)
        
        cbar = fig.colorbar(scalarMap,orientation='vertical')#, format='%2i')
        cbar.set_label('Month')

        labels = np.arange(1, 13 , 1)
        location = labels 
        cbar.set_ticks(location)
        cbar.set_ticklabels(labels)

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        
        #---------------------------------
        #Figure: Individual Prfs
        #---------------------------------
        
        nprfs   = prfmean_day[dfValue].shape[0]

        if self.pdfsav:

            for d, da in enumerate(self.sonde['date']):

                #---------------------------------
                #Index for same days - FTIR
                #---------------------------------
                dtsonde   = self.sonde['dt'][d]
                deltaT    = dtsonde -   dates 
                deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.
                
                inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= dfValue])
                #inds = np.array([x for x, df in enumerate(datefts) if df == da])

                if len(inds) >= 1:

                    fig,(ax1,ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(15,7), sharey=True, sharex=False)
                                  
                    for i in inds:
                        #ax1.plot(rPrf[i]/1e3,alt, label='{:%H:%M}'.format(dates[i]),  linewidth=2.0)  #
                        ax1.errorbar(rPrf[i]/1e3,alt, xerr=tot_errvmr[i]/1e3, markersize=35, linestyle='-', label='{:%H:%M}'.format(dates[i]), linewidth=2.0)
                        ax1.scatter(rPrf[i]/1e3,alt,facecolors='white', s=35, color='k')

                        ax1.plot(aprPrf[i]/1e3,alt, '-', label='{:%H:%M}'.format(dates[i]) )
                        ax1.scatter(aprPrf[i]/1e3,alt,facecolors='white', s=35, color='k')
                      
                    prfmean_i    = np.average(rPrf[inds],axis=0, weights=tot_errvmr[inds])
                    Aprfmean_i   = np.average(aprPrf[inds], axis=0)
                    prferrmean_i = np.sqrt (1.0/np.sum(1.0/tot_errvmr[inds]**2, axis=0))
                    prferrmean_i    = np.sqrt (np.sum(tot_errvmr[inds]**2, axis=0))
                    prfstd_i     = np.std(rPrf[inds], axis=0)


                    ax1.plot(Aprfmean_i/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                    ax1.scatter(Aprfmean_i/1e3, alt,facecolors='white', s=35, color='gray')
                    
                    ax1.plot(np.asarray(self.sonde['H2Omr_ppmv_a'][d])/1e3, np.asarray(self.sonde['Alt_km_a'][d]),color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                    ax1.scatter(np.asarray(self.sonde['H2Omr_ppmv_a'][d])/1e3, np.asarray(self.sonde['Alt_km_a'][d]),facecolors='white', s=35, color='k')
                    ax1.legend(prop={'size':11})

                    ax2.plot(prfmean_i/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                    ax2.scatter(prfmean_i/1e3,alt,facecolors='white', s=35, color='b')
                    ax2.fill_betweenx(alt,prfmean_i/1e3-prferrmean_i/1e3,prfmean_i/1e3+prferrmean_i/1e3,alpha=0.25,color='blue')

                    ax2.plot(Aprfmean_i/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                    ax2.scatter(Aprfmean_i/1e3, alt,facecolors='white', s=35, color='gray')

                    Prfsonde_interp_i = interpolate.interp1d(np.asarray(self.sonde['Alt_km_a'][d]), np.asarray(self.sonde['H2Omr_ppmv_a'][d]), axis=0, fill_value=np.asarray(self.sonde['H2Omr_ppmv_a'][d][0]), bounds_error=False)(alt)
                    
                    #if smthFlg:
                    ax2.plot(Prfsonde_interp_i/1e3, alt, color='k',  linewidth=2.0, label='sonde\n({:%H:%M})'.format(dtsonde))
                    ax2.scatter(Prfsonde_interp_i/1e3, alt,facecolors='white', s=35, color='k')
                        
                    rd                 = (Prfsonde_interp_i - prfmean_i)
                    rda                 = (Prfsonde_interp_i - Aprfmean_i)

                    ax3.plot(rd/1e3, alt, '-', color='blue', linewidth=2.0, label='Retrieved')
                    ax3.scatter(rd/1e3, alt,facecolors='white', s=35, color='b')
                    
                    ax3.plot(rda/1e3, alt, '-', color='gray', linewidth=2.0, label='a priori')
                    ax3.scatter(rda/1e3, alt,facecolors='white', s=35, color='gray')
                    ax3.legend(prop={'size':11})
                    ax3.axvline(x=0, linestyle='--', linewidth=2, color='k')

                    rd = rd/Prfsonde_interp_i *100.
                    rda = rda/Prfsonde_interp_i *100.

                    ax4.plot(rd, alt, '-', color='blue', linewidth=2.0, label='Retrieved')
                    ax4.scatter(rd, alt,facecolors='white', s=35, color='b')
                    
                    ax4.plot(rda, alt, '-', color='gray', linewidth=2.0, label='a priori')
                    ax4.scatter(rda, alt,facecolors='white', s=35, color='gray')
                    ax4.legend(prop={'size':11})
                    ax4.axvline(x=0, linestyle='--', linewidth=2, color='k')

                    ax2.legend(prop={'size':12})

                    ax1.set_ylabel('Altitude [km]')
                    if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
                    if loc.lower() == 'mlo': ax1.set_ylim(3, 15)
                    #ax1.set_ylim(1, 15)
                    ax1.set_xlim(xmin=0)

                    ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]')
                    ax2.set_xlabel('VMR [x10$^3$ ppm$_v$]')
                    ax2.set_xlim(xmin=0)
                    ax3.set_xlabel('Residuals [x10$^3$ ppm$_v$]')
                    ax4.set_xlabel('Relative difference [%]')

                    ax1.grid(True,which='both')
                    ax2.grid(True,which='both')
                    ax3.grid(True,which='both')
                    ax4.grid(True,which='both')
                    #print mf.bias(prfmean_day/1e3,Prfsonde_interp)

                    ax1.tick_params(which='both',labelsize=12)
                    ax2.tick_params(which='both',labelsize=12)
                    ax3.tick_params(which='both',labelsize=12)
                    ax4.tick_params(which='both',labelsize=12)

                    plt.suptitle(da, fontsize=16)    

                    fig.subplots_adjust(bottom=0.1,top=0.94, left=0.05, right=0.95)    
                    #ax1.text(-0.1,1.1, 'Number of Obs Filtered        = '+str(nfltr),ha='left',va='center',transform=ax1.transAxes,fontsize=8)
                    #ax1.text(-0.1,1.05,'Number of Obs After Filtering = '+str(ntot), ha='left',va='center',transform=ax1.transAxes,fontsize=8)
               
                    if self.pdfsav: 
                        self.pdfsav.savefig(fig,dpi=200)
                        plt.close(fig)
                    else:
                        plt.show(block=False)
                        #user_input = raw_input('Press any key to exit >>> ')
                        #sys.exit()

        #---------------------------------
        #Figure: selected profiles (No smoothed)
        #---------------------------------
        if doiplt:
            doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in doiplt]
            doidt = np.asarray(doidt)
            
            fig0, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
            ndoi = 0 

            for d, da in enumerate(sondedates[dfValue]):

                deltadoi = doidt - da

                if ndoi < len(doidt):
                
                    if deltadoi[ndoi] == dt.timedelta(0):

                        if ndoi < len(doidt):

                            if ndoi<=4:

                                if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                                ax0[0, ndoi].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='HR-FTIR')
                                ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                ax0[0, ndoi].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[0, ndoi].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[0, ndoi].plot(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],color='k', label='CFH', linewidth=2.0)
                                ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=20, color='k')

                                if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                
                                ax0[0, ndoi].grid(True,which='both')
                                ax0[0, ndoi].tick_params(which='both',labelsize=14)
                                ax0[0, ndoi].set_title(da, fontsize=14)   

                                ax0[0, ndoi].set_xlim(xmin=0)

                                #ax0[0, ndoi].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')

                                ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)
                                if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                            if ndoi>=5:

                                ax0[1, (ndoi-5)].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                                ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[1, (ndoi-5)].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                ax0[1, (ndoi-5)].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[1, (ndoi-5)].plot(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=20, color='k')

                                if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)
                                #ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)                                

                                if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                                ax0[1, (ndoi-5)].grid(True,which='both')
                                ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                                ax0[1, (ndoi-5)].set_title(da, fontsize=14)

                                ax0[1, (ndoi-5)].set_xlim(xmin=0)

                                

                                ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)
                                #ax0[1, (ndoi-5)].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')
                                #ax0[1, (ndoi-5)].legend(prop={'size':11})

                            fig0.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 
                            ndoi+=1

            if self.pdfsav: 
                self.pdfsav.savefig(fig0,dpi=200)
                plt.savefig('/data/iortega/Manuscripts/Fig/PrFSelected_'+self.loc.upper()+'.pdf', bbox_inches='tight')
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

        #---------------------------------
        #Figure: selected profiles (same grid, smoothed or not)
        #---------------------------------
        if doiplt:
            doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in doiplt]
            doidt = np.asarray(doidt)
            
            fig0, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
            ndoi = 0 

            for d, da in enumerate(sondedates[dfValue]):

                deltadoi = doidt - da

                if ndoi < len(doidt):
                
                    if deltadoi[ndoi] == dt.timedelta(0):

                        if ndoi < len(doidt):

                            if ndoi<=4:

                                if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                                ax0[0, ndoi].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='HR-FTIR')
                                ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                #ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                
                                ax0[0, ndoi].grid(True,which='both')
                                ax0[0, ndoi].tick_params(which='both',labelsize=14)
                                ax0[0, ndoi].set_title(da, fontsize=14)   

                                ax0[0, ndoi].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[0, ndoi].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[0, ndoi].plot(Prfsonde_interp[dfValue][d]/1e3, alt,color='k', label='CFH', linewidth=2.0)
                                ax0[0, ndoi].scatter(Prfsonde_interp[dfValue][d]/1e3, alt,facecolors='white', s=20, color='k')

                                #ax0[0, ndoi].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[0, ndoi].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')

                                ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)

                                ax0[0, ndoi].set_xlim(xmin=0)

                                if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                            if ndoi>=5:

                                ax0[1, (ndoi-5)].plot(prfmean_day[dfValue][d]/1e3,alt, color='b',  linewidth=2.0, label='Mean (Retrieved)')
                                ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
                                ax0[1, (ndoi-5)].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

                                if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)
                                #ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)

                                
                                if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                                ax0[1, (ndoi-5)].grid(True,which='both')
                                ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                                ax0[1, (ndoi-5)].set_title(da, fontsize=14)

                                ax0[1, (ndoi-5)].plot(aPrf_day[dfValue][d]/1e3, alt, '-', color='gray', label='a priori', linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue][d]/1e3, alt,facecolors='white', s=20, color='gray')
                                
                                ax0[1, (ndoi-5)].plot(Prfsonde_interp[dfValue][d]/1e3, alt,color='k', label='sonde\n({:%H:%M})'.format(dtsonde), linewidth=2.0)
                                ax0[1, (ndoi-5)].scatter(Prfsonde_interp[dfValue][d]/1e3, alt,facecolors='white', s=20, color='k')

                                #ax0[1, (ndoi-5)].errorbar(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d], xerr=sondePrf[dfValue][d]*0.1/1e3, color='k', markersize=0, linestyle='-', label='CFH', linewidth=2.0)
                                #ax0[1, (ndoi-5)].scatter(sondePrf[dfValue][d]/1e3, sondealt[dfValue][d],facecolors='white', s=35, color='k')
                                #ax0[1, (ndoi-5)].legend(prop={'size':11})

                                ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)

                                ax0[1, (ndoi-5)].set_xlim(xmin=0)

                            fig0.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 
                            ndoi+=1

            if self.pdfsav: 
                self.pdfsav.savefig(fig0,dpi=200)
                plt.savefig('/data/iortega/Manuscripts/Fig/PrFSelected_'+self.loc.upper()+'_interp.pdf', bbox_inches='tight')
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

        #------------------------------------
        #HISTOGRAM OF RESIDUALS BEFORE FILTERING
        #------------------------------------
        Qpercent = 95
        fig = plt.figure(figsize=(10,7.5))

        print '\nBefore Filtering:'

        for p, pcol in enumerate(pCols[0:-1]):

            gs1 = gridspec.GridSpec(1, 3)

            if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
            if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
            if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
            
            if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
            if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
            if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

            ax1 = plt.subplot(gs1[0:1, :])

            indsH = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]
            x = PrfDiff[dfValue][:,indsH]
            x = x[~np.isnan(x)]
      
            n, bins, patches = ax1.hist(x/1e3, bins=20, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
            P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
            mu = np.nanmean(x/1e3)
            me = np.nanmedian(x/1e3)
            sigma = np.std(x/1e3)
            #add a line showing the expected distribution
            y = P.normpdf( bins, mu, sigma)
            l = ax1.plot(bins, y, 'k--', linewidth=3)

            ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent), color='k', linestyle='--')
            ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent)*-1., color='k', linestyle='--')

            ax1.axvline(x=mu, color='blue', linestyle='--')
            ax1.axvline(x=me, color='green', linestyle='--')

            ax1.axvline(x=mu + sigma, color='r', linestyle='--')
            ax1.axvline(x=mu - sigma, color='r', linestyle='--')


            ax1.grid(True,which='both')
            ax1.tick_params(which='both',labelsize=12)#, labelbottom='off')
            #ax1.set_ylabel('Probability', fontsize=12)
            ax1.set_title(str(pcol[0])+' - '+str(pcol[1])+' km',
                      horizontalalignment='center', verticalalignment='baseline', fontsize=14)
            #ax1.set_xlim(-1, 1)

            if (p == 4) or (p == 5): 
                ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                ax1.tick_params(labelbottom='on')
            if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

           
            print '\nAltitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
            #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
            print 'Mean Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, sigma)
            print 'Median Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, sigma)
            print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100.)
            print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100.)

            print 'Percentile 95 in Delta T {} min : {}'.format(dfValue, np.percentile(np.abs(x/1e3), Qpercent))


        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/Hist_Apriori'+self.loc.upper()+'.pdf', bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        #------------------------------------
        #
        #------------------------------------
        Prct        = {}

        for df in diffT:

            indsBad = []

            #print 'Percentile 95 in Delta T {} min : {}'.format(df, np.percentile(np.abs(PrfDiffRel[df][:,:]), Qpercent))

            for pn, pcol in enumerate(pCols):

                indsH = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

                print 'Percentile 95 in Delta T {} min and Pcol {} : {}'.format(df, pcol, np.percentile(np.abs(PrfDiff[df][:,indsH]), Qpercent))

                Prct.setdefault(df, []).append(np.percentile(np.abs(PrfDiffRel[df][:,indsH]), Qpercent))

                if self.loc.lower() == 'fl0': iBad = np.where( (np.abs(PrfDiff[df][:,indsH]) >= np.percentile(np.abs(PrfDiff[df][:,indsH]), Qpercent)) | (PrfDist[df][:,indsH] > 300.0) )
                if self.loc.lower() == 'mlo': iBad = np.where( (np.abs(PrfDiff[df][:,indsH]) >= np.percentile(np.abs(PrfDiff[df][:,indsH]), Qpercent)))


                #for d, da in enumerate(FTSdates[df]): 
                #if self.fleoutFlg: iBad = np.where( (np.abs(PrfDiffRel[df][d]) >= 300.0) | (PrfDist[df][d] > 300.0) )[0]
                #if self.fleoutFlg: iBad = np.where( (np.abs(PrfDiffRel[df]) >= np.percentile(np.abs(PrfDiffRel[df]), 95)) | (PrfDist[df] > 300.0) )
                #else:              iBad = np.where( (np.abs(PrfDiffRel[df][d]) >= 300.0) )[0]
                #iBad = np.where( (np.abs(PrfDist[df][d]) >= 300.0)  )[0]

                iBad = np.asarray(iBad)
                #print iBad.shape
                # print iBad.shape[1]
                #print iBad
                # print PrfDiffRel[df].shape
                # print PrfDiffRel[df].shape[0]*PrfDiffRel[df].shape[1]
                # exit()
                #print prfmean_day[df].shape

                #print iBad
                #print iBad[0]
                #print iBad[1]
        
                    #if len(iBad) >= 1:

                    #    indsBad.append(d)

                #print 'Percentile 95 in Delta T = {} minutes: {}'.format(df, np.percentile(np.abs(PrfDiffRel[df]), Qpercent))
                #print float(iBad.shape[1])
                #print float(PrfDiffRel[df][:,indsH].shape[0]) * (PrfDiffRel[df][:,indsH].shape[1])
                print 'Percent of bad points in Delta T = {} minutes: {}'.format(df, float(iBad.shape[1])/float((PrfDiffRel[df][:,indsH].shape[0] * PrfDiffRel[df][:,indsH].shape[1]))*100.)
                #print PrfDiffRel[df]
                prfmean_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(prfmean_day[df], indsBad, axis=0)
                
                prferr_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(prferr_day[df], indsBad, axis=0)
                prfSTD_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(prfSTD_day[df], indsBad, axis=0)
                aPrf_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                        = np.delete(aPrf_day[df], indsBad, axis=0)
                Airmassmean[df][iBad[0], indsH[iBad[1]]] = np.nan#                     = np.delete(Airmassmean[df], indsBad, axis=0)
                avkSCFday[df][iBad[0], indsH[iBad[1]]] = np.nan#                       = np.delete(avkSCFday[df], indsBad, axis=0)
                rPrfMolmean_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                 = np.delete(rPrfMolmean_day[df], indsBad, axis=0)
                prferrMol_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(prferrMol_day[df], indsBad, axis=0)
                aPrfMol_day[df][iBad[0], indsH[iBad[1]]] = np.nan#                     = np.delete(aPrfMol_day[df], indsBad, axis=0)
                #FTSdates[df][indsBad][:] = np.nan#                        = np.delete(FTSdates[df], indsBad)
                #NobsFTS[df][indsBad][:] = np.nan#                         = np.delete(NobsFTS[df], indsBad)

                Prfsonde_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#                 = np.delete(Prfsonde_interp[df], indsBad, axis=0)
                Prfsonde_sd_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsonde_sd_interp[df], indsBad, axis=0)
                PrfsondeMol_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(PrfsondeMol_interp[df], indsBad, axis=0)
                PrfsondeMol_sd_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#           = np.delete(PrfsondeMol_sd_interp[df], indsBad, axis=0)
                sondeairmass_a_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#           = np.delete(sondeairmass_a_interp[df], indsBad, axis=0)
                #sondedates[df][indsBad][:] = np.nan#                      = np.delete(sondedates[df], indsBad)
                
                PrfDiff[df][iBad[0], indsH[iBad[1]]] = np.nan#                         = np.delete(PrfDiff[df], indsBad, axis=0)
                PrfDiffApr[df][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(PrfDiffApr[df], indsBad, axis=0)
                PrfDiffRel[df][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(PrfDiffRel[df], indsBad, axis=0)
                #PrfDiffRel[df]                      = np.delete(PrfDiffRel[df], indsBad, axis=0)
                #print '\n'
                #print PrfDiffRel[df]
                #exit()
                PrfDiffAprRel[df][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(PrfDiffAprRel[df], indsBad, axis=0)

                if self.fleoutFlg:
                    Prfsondelon_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsondelon_interp[df], indsBad, axis=0)
                    Prfsondelat_interp[df][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsondelat_interp[df], indsBad, axis=0)
                    PrfDist[df][iBad[0], indsH[iBad[1]]] = np.nan#                         = np.delete(PrfDist[df], indsBad, axis=0)

            #for pn, pcol in enumerate(pCols):

            #    inds = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

                # sondevmrP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(Prfsonde_interp[df][:,indsH]), axis=1,  weights=Airmassmean[df][:,indsH])
                # sondevmrsdP[str(df)+'_'+str(pn)]  = np.ma.average(np.ma.masked_invalid(Prfsonde_sd_interp[df][:, indsH]), axis=1, weights=Airmassmean[df][:,indsH])
                # sondesumP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(Prfsonde_interp[df][:,indsH])*np.ma.masked_invalid(Airmassmean[df][:,indsH])*1e-6, axis=1)
                # sondesumsdP[str(df)+'_'+str(pn)]  = np.nansum(np.ma.masked_invalid(Prfsonde_sd_interp[df][:,indsH])*np.ma.masked_invalid(Airmassmean[df][:,indsH])*1e-6, axis=1)
                
                # FTSvmrP[str(df)+'_'+str(pn)]      = np.ma.average(np.ma.masked_invalid(prfmean_day[df][:,indsH]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,indsH]))
                # FTSvmrsdP[str(df)+'_'+str(pn)]    = np.ma.average(np.ma.masked_invalid(prferr_day[df][:,indsH]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,indsH]))
                # FTSAprvmrP[str(df)+'_'+str(pn)]   = np.ma.average(np.ma.masked_invalid(aPrf_day[df][:,indsH]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[df][:,indsH]))
                # FTSsumP[str(df)+'_'+str(pn)]      = np.nansum(np.ma.masked_invalid(rPrfMolmean_day[df][:,indsH]), axis=1)
                # FTSsumsdP[str(df)+'_'+str(pn)]    = np.nansum(np.ma.masked_invalid(prferrMol_day[df][:,indsH]), axis=1)
                # FTSAprsumP[str(df)+'_'+str(pn)]   = np.nansum(np.ma.masked_invalid(aPrfMol_day[df][:,indsH]), axis=1)

                if self.fleoutFlg:
                    sondelatP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[df][:,indsH]), axis=1)
                    sondelonP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[df][:,indsH]), axis=1)

                    #distP[str(df)+'_'+str(pn)]        = [mf.haversine(np.nanmean(lonmean[df]), np.mean(latmean[df]), lo, la ) for (lo, la) in  zip(sondelonP[str(df)+'_'+str(pn)], sondelatP[str(df)+'_'+str(pn)])]
                    distP[str(df)+'_'+str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[df][:,indsH]), axis=1)

                    # for pn, pcol in enumerate(pCols):

                    #     sondevmrP[str(df)+'_'+str(pn)]    = np.delete(sondevmrP[str(df)+'_'+str(pn)], indsBad)
                    #     sondevmrsdP[str(df)+'_'+str(pn)]  = np.delete(sondevmrsdP[str(df)+'_'+str(pn)], indsBad)
                    #     sondesumP[str(df)+'_'+str(pn)]    = np.delete(sondesumP[str(df)+'_'+str(pn)], indsBad)
                    #     sondesumsdP[str(df)+'_'+str(pn)]  = np.delete(sondesumsdP[str(df)+'_'+str(pn)], indsBad)
                        
                    #     FTSvmrP[str(df)+'_'+str(pn)]      = np.delete(FTSvmrP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSvmrsdP[str(df)+'_'+str(pn)]    = np.delete(FTSvmrsdP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSAprvmrP[str(df)+'_'+str(pn)]   = np.delete(FTSAprvmrP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSsumP[str(df)+'_'+str(pn)]      = np.delete(FTSsumP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSsumsdP[str(df)+'_'+str(pn)]    = np.delete(FTSsumsdP[str(df)+'_'+str(pn)], indsBad)
                    #     FTSAprsumP[str(df)+'_'+str(pn)]   = np.delete(FTSAprsumP[str(df)+'_'+str(pn)] , indsBad)
                    #     distP[str(df)+'_'+str(pn)]        = np.delete(distP[str(df)+'_'+str(pn)] , indsBad)
            #print PrfDiffRel[df]
            #exit()

        #------------------------------------
        #HISTOGRAM OF RESIDUALS AFTER FILTERING
        #------------------------------------
        fig = plt.figure(figsize=(10,7.5))

        for p, pcol in enumerate(pCols[0:-1]):

            gs1 = gridspec.GridSpec(1, 3)

            if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
            if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
            if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
            
            if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
            if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
            if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

            ax1 = plt.subplot(gs1[0:1, :])

            indsH = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]
            x = PrfDiff[dfValue][:,indsH]
            x = x[~np.isnan(x)]
        
            n, bins, patches = ax1.hist(x/1e3, bins=20, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
            P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
            mu = np.nanmean(x/1e3)
            me = np.nanmedian(x/1e3)
            sigma = np.std(x/1e3)
            #add a line showing the expected distribution
            y = P.normpdf( bins, me, sigma)
            l = ax1.plot(bins, y, 'k--', linewidth=3)

            ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent), color='k', linestyle='--')
            ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent)*-1., color='k', linestyle='--')

            ax1.axvline(x=mu, color='blue', linestyle='--')
            ax1.axvline(x=me, color='green', linestyle='--')

            ax1.axvline(x=mu + sigma, color='r', linestyle='--')
            ax1.axvline(x=mu - sigma, color='r', linestyle='--')


            ax1.grid(True,which='both')
            ax1.tick_params(which='both',labelsize=12)#, labelbottom='off')
            #ax1.set_ylabel('Probability', fontsize=12)
            ax1.set_title(str(pcol[0])+' - '+str(pcol[1])+' km',
                      horizontalalignment='center', verticalalignment='baseline', fontsize=14)
            #ax1.set_xlim(-1, 1)

            if (p == 4) or (p == 5): 
                ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                ax1.tick_params(labelbottom='on')
            if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

            print '\nAltitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
            #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
            print 'Mean Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, sigma)
            print 'Median Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, sigma)
            print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100.)
            print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue][:, indsH]/1e3) * 100.)

            # print 'Percentile 95 in Delta T {} min : {}'.format(dfValue, np.percentile(np.abs(x/1e3), Qpercent))


        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/Hist_Apriori'+self.loc.upper()+'.pdf', bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #------------------------------------
        #Fig: Stats with Coincident Interval
        #------------------------------------
        fig1, ax1  = plt.subplots(2, figsize=(9,8), sharex=True)
        clr = mf.clrplt()

        if self.loc.lower() == 'fl0':
            fileName = 'Variability_Temporal_FL0.ascii'

        
        for pn, pcol in enumerate(pCols):

            stdDF        = []
            Npnts        = []
            Npnts2       = []
            biasDF       = []
            rmseDF       = []

            for df in diffT:

                inds = np.where( (alt > pcol[0]) & (alt <=pcol[1])  )[0]


                x = PrfDiffRel[df][:, inds]
                x = np.asarray(x)
                x = x.flatten()

                biasDF.append(np.nansum( x) / len(x))

                ftsPrf = [i[inds] for i in prfmean_day[df]]
                sonPrf = [i[inds] for i in Prfsonde_interp[df]]

                st = prfSTD_day[df][:, inds]
                st = np.asarray(st)
                st = st.flatten()

                stdDF.append(float(np.nanmean(st)/np.nanmean(ftsPrf)))


                ftsPrf = np.asarray(ftsPrf)
                ftsPrf = ftsPrf.flatten()

                sonPrf = np.asarray(sonPrf)
                sonPrf = sonPrf.flatten()

                #------------------------------
                # Root Mean Square Error (RMSE)
                #------------------------------
                ss_res = np.nansum( x**2)
                rmseDF.append(np.sqrt( ss_res / len(x) ))

                Npnts.append(float(np.sum(NobsFTS[df])))
                Npnts2.append(float(len(NobsFTS[df])))

            Npnts  = np.asarray(Npnts)
            Npnts2 = np.asarray(Npnts2)
            stdDF  = np.asarray(stdDF)


            #with open('/data/iortega/Manuscripts/Fig/'+fileName, 'wb') as fopen:
            #    fopen.write('{0:.4f}\t{1:.4f}\t{2:.4f}\n'.format(Npnts2, Npnts, np.asarray(stdDF)*100.))


            if pn == 0:
                ax1[0].scatter(diffT, Npnts2, color='k',s=60, label='# of Dates')
                ax1[0].plot(diffT, Npnts2,  color='k')
                ax1[0].scatter(diffT, Npnts, color='blue',s=60, label='# of Profiles')
                ax1[0].plot(diffT, Npnts,  color='blue')

                print '# of Profiles = {}'.format(Npnts)
                print '# of Dates = {}'.format(Npnts2)
                print 'diffT = {}'.format(diffT)

            ax1[1].scatter(diffT, np.asarray(stdDF)*100., color=clr[pn], s=60, label= str(pcol[0])+'-'+str(pcol[1])+' km')
            ax1[1].plot(diffT, np.asarray(stdDF)*100.,   color=clr[pn])
    

        ax1[0].grid(True) 
        ax1[0].tick_params(which='both',labelsize=14)
        ax1[0].set_ylabel('Number of observations', fontsize=14)
        ax1[0].set_ylim(bottom=0)
        ax1[0].legend(prop={'size':11}, loc=2)

        ax1[1].grid(True)
        ax1[1].set_ylabel('Variability [%]', fontsize=14)
        ax1[1].tick_params(which='both',labelsize=14)
        ax1[1].set_ylim(bottom=0)
        ax1[1].legend(prop={'size':11}, loc=2)
        ax1[1].set_xlim(left=0, right=np.max(diffT) + 10)
        ax1[1].set_xlabel('$\Delta$t [min]', fontsize=14)


        if self.pdfsav: 
            self.pdfsav.savefig(fig1,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/StatTime_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        # if self.loc.lower() == 'fl0':
        #     fileName = 'Variability_Temporal_FL0.ascii'

        #     with open('/data/iortega/Manuscripts/Fig/'+fileName, 'wb') as fopen:
            
        #         fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('diffT', 'Num_Dates', '', 'Num_Prf', 'Var_perc'))
        #         for i, s in enumerate(diffT):
                    
        #             fopen.write('{0:.4f}\t{1:.4f}\t{2:.4f}\n'.format(diffT[i], Npnts2[i], Npnts[i], bias_e[i]))



        #---------------------------------
        #CALCULATION OF BIAS (MEDIAN OF DIFFERENCES), PRECISION (STDV OF RESIDUALS), CORRELATION AND HISTOGRAM OF BIAS
        #---------------------------------
        fig  = plt.figure(figsize=(10,7.5))
        fig2  = plt.figure(figsize=(10,7.5))
        outer_grid = gridspec.GridSpec(3, 2, wspace=0.15, hspace=0.15)

        for df in diffT:

            if df == dfValue:

                bias         = []
                bias_e       = []

                prec         = []
                prec_e       = []

                bias_perc    = []
                bias_perc_e  = []

                prec_perc    = []
                prec_perc_e  = []


                biasApr      = []
                bias_eApr    = []

                precApr      = []
                prec_eApr    = []

                slope        = []
                slope_e      = []
                
                intercept    = []
                intercept_e  = []
                
                rvalue       = []

                slopeApr     = []
                slope_eApr     = []
                
                interceptApr = []
                intercept_eApr = []
                rvalueApr    = []

                #---------------------------------
                # Figure: HISTOGRAMS FOR DIFFERENT ALTITUDE LAYERS
                #---------------------------------
                print '\nBias and Correlation:'
                for p, pcol in enumerate(pCols[0:-1]):

                    inds = np.where( (alt >= pcol[0]) & (alt <pcol[1])  )[0]

                    #--------------
                    #Bias and precision in Retrieval
                    #--------------
                    bias_n = np.sum(prfmean_day[df][:, inds] - Prfsonde_interp[df][:, inds], axis=1 )/float(len(inds))#Prfsonde_interp[df][:, inds].shape[1] 
                    #biasCalc = prfmean_day[df][:, inds] - Prfsonde_interp[df][:, inds]
                    bias_n = bias_n[~np.isnan(bias_n)]/1e3

                    mu = np.nanmean(bias_n)
                    me = np.nanmedian(bias_n)
                    sigma = np.std(bias_n)
                    stdE = sigma/sqrt(len(bias_n))      #http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf
                    prec_n = sigma/sqrt(len(bias_n)) * 2.

                    bias.append(me)
                    bias_e.append(stdE)

                    bias_perc.append(me/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)
                    bias_perc_e.append(stdE/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)

                    prec.append(prec_n)
                    prec_e.append(stdE * 0.71)

                    prec_perc.append(prec_n/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)
                    prec_perc_e.append((stdE * 0.71)/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)

                    #--------------
                    #Bias and precision in Apriori
                    #--------------
                    bias_nApr = np.sum(aPrf_day[df][:, inds] - Prfsonde_interp[df][:, inds], axis=1 )/float(len(inds))#Prfsonde_interp[df][:, inds].shape[1] 
                    bias_nApr = bias_nApr[~np.isnan(bias_nApr)]/1e3

                    muApr = np.nanmean(bias_nApr)
                    meApr = np.nanmedian(bias_nApr)
                    sigmaApr = np.std(bias_nApr)
                    stdEApr = sigmaApr/sqrt(len(bias_nApr)) * 2.  #http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf
                    #add a line showing the expected distribution



                    biasApr.append(meApr)
                    bias_eApr.append(sigmaApr/np.sqrt(len(bias_nApr)))

                    precApr.append(sigmaApr/np.sqrt(len(bias_nApr)) * 2.)
                    prec_eApr.append(sigmaApr/np.sqrt(len(bias_nApr)) * 0.71)

                    #--------------
                    #Orthogonal Regression in Retrieval
                    #--------------
                    xx = Prfsonde_interp[df][:, inds]
                    xx = xx[~np.isnan(xx)]/1e3

                    xx_e = Prfsonde_sd_interp[df][:, inds]
                    xx_e = xx_e[~np.isnan(xx_e)]/1e3

                    yy = prfmean_day[df][:, inds]
                    yy = yy[~np.isnan(yy)]/1e3

                    yy_e = prferr_day[df][:, inds]
                    yy_e = yy_e[~np.isnan(yy_e)]/1e3

                    odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx_e, yerr=yy_e, InError=True)
                    slope.append(float(odr[0]))
                    intercept.append(float(odr[1]))

                    slope_e.append(float(odrErr[0]))
                    intercept_e.append(float(odrErr[1]))

                    slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(xx, yy)
                    rvalue.append(float(r_value2))
                    

                    #--------------
                    #Orthogonal Regression in Apriori
                    #--------------
                    yyApr = aPrf_day[df][:, inds]
                    yyApr = yyApr[~np.isnan(yyApr)]/1e3

                    odrApr, odrErrApr  = mf.orthoregress(xx, yyApr, xerr=xx_e, yerr=yyApr*0.,InError=True)
                    slopeApr.append(float(odrApr[0]))
                    interceptApr.append(float(odrApr[1]))

                    slope_eApr.append(float(odrErrApr[0]))
                    intercept_eApr.append(float(odrErrApr[1]))

                    slope2Apr, intercept2Apr, r_value2Apr, p_value2Apr, std_err2Apr = stats.linregress(xx, yyApr)
                    
                    rvalueApr.append(float(r_value2Apr))
                    
                    #-------------- 
                    #Plot: Histogram
                    #--------------                  
                    ax1 = plt.Subplot(fig, outer_grid[p])

                    #normed=True means that the total area under the histogram is equal to 1 but the sum of heights is not equal to 1    
                    n, bins, patches = ax1.hist(bias_n, bins=20, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
                    P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)

                    y = P.normpdf(bins, me, sigma)
                    l = ax1.plot(bins, y, 'k--', linewidth=3)

                    ax1.grid(True,which='both')
                    ax1.tick_params(which='both',labelsize=12)#, labelbottom='off')
                    ax1.axvline(x=mu, color='blue', linestyle='--')
                    ax1.axvline(x=me, color='aqua', linestyle='--')
                    ax1.axvline(x=mu + sigma, color='r', linestyle='--')
                    ax1.axvline(x=mu - sigma, color='r', linestyle='--')
                   
                    if (p == 4) or (p == 5): 
                        ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                        ax1.tick_params(labelbottom='on')
                    if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

                    ax1.annotate(str(pcol[0])+' - '+str(pcol[1])+' km',  xycoords='axes fraction', xy=(0.035, 0.875), color='k',  fontsize=14, ha='left')

                    fig.add_subplot(ax1)
                    fig.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95)

                    #--------------
                    #Print in Terminal 
                    #--------------
                    
                    print '\nAltitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
                    #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
                    print 'Mean Bias Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, stdE)
                    print 'Median Bias Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, stdE)
                    print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100., stdE/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)
                    print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100., stdE/np.nanmean(Prfsonde_interp[df][:, inds]/1e3) * 100.)
        
                    #-------------- 
                    #Plot: Histogram
                    #--------------                  
                    ax2 = plt.Subplot(fig2, outer_grid[p])

                    ax2.errorbar(xx, yy, yerr=yy_e, xerr=xx_e, color='b', markersize=60, linestyle='none')
                    ax2.scatter(xx, yy, color='b', s=60, edgecolor='k', label= str(pcol[0])+'-'+str(pcol[1])+' km')

                    o2o = np.arange(0, 50)

                    ax2.plot(o2o, o2o, color= 'k', linestyle='--')
                    ax2.plot(o2o*0.8, o2o, color= 'k', linestyle='--')
                    ax2.plot(o2o, o2o*0.8, color= 'k', linestyle='--')
                    

                    ax2.grid(True,which='both')
                    ax2.tick_params(which='both',labelsize=12)#, labelbottom='off')

                    ax2.set_ylim(0, np.max(xx) + np.max(xx)*0.15)
                    ax2.set_xlim(0, np.max(xx) + np.max(xx)*0.15)

                    if (p == 4) or (p == 5): 
                        ax2.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                        #ax1.tick_params(labelbottom='on')
                    if (p == 0) or (p == 2) or (p == 4): ax2.set_ylabel('VMR [x10$^3$ ppm]', fontsize=14)


                    ax2.annotate(str(pcol[0])+' - '+str(pcol[1])+' km',  xycoords='axes fraction', xy=(0.035, 0.875), color='k',  fontsize=14, ha='left')
                    ax2.annotate('slope = {0:.2f} $\pm$ {1:.2f}'.format(float(odr[0]), float(odrErr[0])) ,  xycoords='axes fraction', xy=(0.4, 0.125), color='k',  fontsize=12, ha='left')
                    ax2.annotate('intercept = {0:.3f} $\pm$ {1:.3f}'.format(float(odr[1]), float(odrErr[1])) ,  xycoords='axes fraction', xy=(0.4, 0.05), color='k',  fontsize=12, ha='left')

                    fig2.add_subplot(ax2)
                    fig2.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95)

                    #--------------
                    #Print in Terminal 
                    #--------------
                    print 'Slope: {0:.2f} +/- {1:.2f}'.format(float(odr[0]), float(odrErr[0]))
                    print 'Intercept = {0:.3f} +/- {1:.3f}'.format(float(odr[1]), float(odrErr[1]))
                    print 'R value = {0:.2f}'.format(float(r_value2))

                    print 'Slope Apriori: {0:.2f} +/- {1:.2f}'.format(float(odrApr[0]), float(odrErrApr[0]))
                    print 'Intercept Apriori = {0:.3f} +/- {1:.3f}'.format(float(odrApr[1]), float(odrErrApr[1]))
                    print 'R value Apriori = {0:.2f}'.format(float(r_value2Apr))
                    


                if self.pdfsav: 
                    self.pdfsav.savefig(fig,dpi=200)
                    plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)
        
        #----------------------------------------
        #
        #----------------------------------------

        if self.ver == 'Current_ERA_v66' : fileName = 'Stats_ERA_v66_'+self.loc.lower()+'.ascii'
        elif self.ver == 'Current_ERA' : fileName = 'Stats_ERA_'+self.loc.lower()+'.ascii'
        elif self.ver == 'Current_NCEP' : fileName = 'Stats_NCEP_'+self.loc.lower()+'.ascii'
        elif self.ver == 'Current_WACCM' : fileName = 'Stats_WACCM_'+self.loc.lower()+'.ascii'
        else: sys.error()

        if self.loc.lower() == 'fl0':

            with open('/data/iortega/Manuscripts/Fig/'+fileName, 'wb') as fopen:
            
                if self.ver == 'Current_ERA_v66' : fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_ERAv66', 'int_ERAv66', 'rP_ERAv66', 'slp_e_ERAv66', 'int_e_ERAv66', 'slpApr_ERAv66', 'intApr_ERAv66', 'rPApr_ERAv66', 'bias_ERAv66', 'precision_ERAv66', 'bias_e_ERAv66', 'precision_e_ERAv66', 'bias_percent_ERAv66', 'bias_percent_e_ERAv66', 'precision_percent_ERAv66', 'precision_percent_e_ERAv66'))
                elif self.ver == 'Current_ERA' :   fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_ERA', 'int_ERA', 'rP_ERA', 'slp_e_ERA', 'int_e_ERA', 'slpApr_ERA', 'intApr', 'rPApr_ERA', 'bias_ERA', 'precision_ERA', 'bias_e_ERA', 'precision_e_ERA', 'bias_percent_ERA', 'bias_percent_e_ERA', 'precision_percent_ERA', 'precision_percent_e_ERA'))
                elif self.ver == 'Current_NCEP' :    fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_NCEP', 'int_NCEP', 'rP_NCEP', 'slp_e_NCEP', 'int_e_NCEP', 'slpApr_NCEP', 'intApr_NCEP', 'rPApr_NCEP', 'bias_NCEP', 'precision_NCEP', 'bias_e_NCEP', 'precision_e_NCEP', 'bias_percent_NCEP', 'bias_percent_e_NCEP', 'precision_percent_NCEP', 'precision_percent_e_NCEP'))
                elif self.ver == 'Current_WACCM' : fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_WACCM', 'int_WACCM', 'rP_WACCM', 'slp_e_WACCM', 'int_e_WACCM', 'slpApr_WACCM', 'intApr_WACCM', 'rPApr_WACCM', 'bias_WACCM', 'precision_WACCM', 'bias_e_WACCM', 'precision_e_WACCM', 'bias_percent_WACCM', 'bias_percent_e_WACCM', 'precision_percent_WACCM', 'precision_percent_e_WACCM'))
                for i, s in enumerate(slope):
                    
                    fopen.write('{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}\t{13:.4f}\t{14:.4f}\t{15:.4f}\n'.format(s, intercept[i], rvalue[i], slope_e[i], intercept_e[i], slopeApr[i], interceptApr[i], rvalueApr[i], bias[i], prec[i], bias_e[i], prec_e[i], bias_perc[i], bias_perc_e[i], prec_perc[i], prec_perc_e[i] ))
        
        elif self.loc.lower() == 'mlo':

            with open('/data/iortega/Manuscripts/Fig/'+fileName, 'wb') as fopen:
            
                if self.ver == 'Current_ERA_v66' : fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_ERAv66_mlo', 'int_ERAv66_mlo', 'rP_ERAv66_mlo', 'slp_e_ERAv66_mlo', 'int_e_ERAv66_mlo', 'slpApr_ERAv66_mlo', 'intApr_ERAv66_mlo', 'rPApr_ERAv66_mlo', 'bias_ERAv66_mlo', 'precision_ERAv66_mlo', 'bias_e_ERAv66_mlo', 'precision_e_ERAv66_mlo', 'bias_percent_ERAv66_mlo', 'bias_percent_e_ERAv66_mlo', 'precision_percent_ERAv66_mlo', 'precision_percent_e_ERAv66_mlo'))
                elif self.ver == 'Current_ERA' :   fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_ERA_mlo', 'int_ERA_mlo', 'rP_ERA_mlo', 'slp_e_ERA_mlo', 'int_e_ERA_mlo', 'slpApr_ERA_mlo', 'intApr_mlo', 'rPApr_ERA_mlo', 'bias_ERA_mlo', 'precision_ERA_mlo', 'bias_e_ERA_mlo', 'precision_e_ERA_mlo', 'bias_percent_ERA_mlo', 'bias_percent_e_ERA_mlo', 'precision_percent_ERA_mlo', 'precision_percent_e_ERA_mlo'))
                elif self.ver == 'Current_NCEP' :    fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_NCEP_mlo', 'int_NCEP_mlo', 'rP_NCEP_mlo', 'slp_e_NCEP_mlo', 'int_e_NCEP_mlo', 'slpApr_NCEP_mlo', 'intApr_NCEP_mlo', 'rPApr_NCEP_mlo', 'bias_NCEP_mlo', 'precision_NCEP_mlo', 'bias_e_NCEP_mlo', 'precision_e_NCEP_mlo', 'bias_percent_NCEP_mlo', 'bias_percent_e_NCEP_mlo', 'precision_percent_NCEP_mlo', 'precision_percent_e_NCEP_mlo'))
                elif self.ver == 'Current_WACCM' : fopen.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('slp_WACCM_mlo', 'int_WACCM_mlo', 'rP_WACCM_mlo', 'slp_e_WACCM_mlo', 'int_e_WACCM_mlo', 'slpApr_WACCM_mlo', 'intApr_WACCM_mlo', 'rPApr_WACCM_mlo', 'bias_WACCM_mlo', 'precision_WACCM_mlo', 'bias_e_WACCM_mlo', 'precision_e_WACCM_mlo', 'bias_percent_WACCM_mlo', 'bias_percent_e_WACCM_mlo', 'precision_percent_WACCM_mlo', 'precision_percent_e_WACCM_mlo'))
                for i, s in enumerate(slope):
                    
                    fopen.write('{0:.4f}\t{1:.4f}\t{2:.4f}\t{3:.4f}\t{4:.4f}\t{5:.4f}\t{6:.4f}\t{7:.4f}\t{8:.4f}\t{9:.4f}\t{10:.4f}\t{11:.4f}\t{12:.4f}\t{13:.4f}\t{14:.4f}\t{15:.4f}\n'.format(s, intercept[i], rvalue[i], slope_e[i], intercept_e[i], slopeApr[i], interceptApr[i], rvalueApr[i], bias[i], prec[i], bias_e[i], prec_e[i], bias_perc[i], bias_perc_e[i], prec_perc[i], prec_perc_e[i] ))
        
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()

        #----------------------------------------
        

        #-------------------
        #Bar Plot with slope, intercept and Rvalue
        #-------------------
        fig, (ax, ax2, ax3)  = plt.subplots(3,1,figsize=(7.5,9), sharex=True)

        ind = np.arange(len(slope))                # the x locations for the groups
        width = 0.35                      # the width of the bars
        rects1 = ax.bar(ind, slope, width, color='green', label='Retrieved', yerr=slope_e, error_kw=dict(elinewidth=1.0, ecolor='k'))
        rects2 = ax.bar(ind+width, slopeApr, width, color='blue', label='A priori')

        ax.set_xlim(-width,len(ind)+width*0.5)
        ax.set_ylabel('Slope', fontsize=14)
        ax.set_ylim(0,1.5)

        ax.axhline(y=1.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)

        #xTickMarks = [str(pcol[0])+'-'+str(pcol[1])+ ' km' for pcol in pCols]
        ax.set_xticks(ind+width)
        #xtickNames = ax.set_xticklabels(xTickMarks)
        #plt.setp(xtickNames, rotation=45, fontsize=12)
        ax.tick_params(which='both',labelsize=14)
        ax.legend(prop={'size':12}, loc=1)

        rects1 = ax2.bar(ind, intercept, width, color='green', yerr=intercept_e, error_kw=dict(elinewidth=1.0, ecolor='k'))
        rects2 = ax2.bar(ind+width, interceptApr, width, color='blue')

        ax2.set_xlim(-width,len(ind)+width*0.5)
        ax2.set_ylabel('Intercept [x10$^3$ ppm$_v$]', fontsize=14)
        #xTickMarks = [str(pcol[0])+'-'+str(pcol[1])+ ' km' for pcol in pCols]
        ax2.set_xticks(ind+width)
        #xtickNames = ax2.set_xticklabels(xTickMarks)
        #plt.setp(xtickNames, rotation=45, fontsize=12)
        ax2.tick_params(which='both',labelsize=14)

        ax2.axhline(y=0.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)

        rects1 = ax3.bar(ind, rvalue, width, color='green')
        rects2 = ax3.bar(ind+width, rvalueApr, width, color='blue')

        ax3.set_xlim(-width,len(ind)+width*0.5)
        ax3.set_ylabel('r-value', fontsize=14)
        xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
        ax3.set_xticks(ind+width)
        xtickNames = ax3.set_xticklabels(xTickMarks)
        plt.setp(xtickNames, rotation=0, fontsize=11)
        ax3.tick_params(which='both',labelsize=14)
        ax3.set_xlabel('Layer [km]', fontsize=14)

        fig.subplots_adjust(bottom=0.075,top=0.975, left=0.15, right=0.95)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/StatsBars_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()
        
                
        #--------------------------------
        # Percentile plot
        #--------------------------------
        #60 MIN
        #Prct_ERA6    = [32.34,  48.72,  150.48, 186.52,  109.21, 119.23, 54.83]
        #Prct_ERA     = [39.61, 102.59,  251.55, 217.72, 91.24, 102.54, 59.65]
        #Prct_NCEP    = [43.57, 130.14, 262.20,  199.32,  306.48, 360.813,  89.95499549]
        #Prct_WACCM   = [44.94, 117.74,  266.45,  336.41,  145.82,  116.28,   56.13]

        #30 MIN
        Prct_ERA6    = [38.26473116,   90.38224722,  218.17328925,  120.93111681,   99.24577035, 116.7481478,    75.58856323]
        Prct_ERA     = [40.45348291,  199.46864549,  262.56195281,  192.43043034,   79.89317067, 97.266399,     59.4740131]
        Prct_NCEP    = [53.00947181,  285.95735394,  285.22738118,  182.02674434,  240.68429733,  271.58422531,   90.11004353]
        Prct_WACCM   = [56.01401609,  266.68757646,  464.99212307,  550.66765635,  201.97047404 ,105.91502458,   76.86148022]
        
        midP = [np.mean(p) for p in  pCols]
    
        clr = mf.clrplt()

        #----------------------------------------
        #
        #----------------------------------------
        ind = np.arange(len(Prct_ERA6))
        fig, ax = plt.subplots(figsize=(9, 6))

        xlabel = [str(p[0])+'-'+str(p[1]) for p in pCols]

        ax.bar(ind-0.2, Prct_ERA6, 0.2, align='center', color = 'r', ecolor = 'k', label='ERA-I-6')
        ax.bar(ind, Prct_ERA, 0.2, align='center', color = 'b', ecolor = 'k', label='ERA-I-d')
        ax.bar(ind+0.2, Prct_NCEP, 0.2, align='center', color = 'green', ecolor = 'k', label='NCEP-d')
        ax.bar(ind+0.2*2, Prct_WACCM, 0.2,  align='center', color = 'gray', ecolor = 'k', label='WACCM')

        ax.legend(prop={'size':11}, loc=2)
        
        ax.yaxis.grid(True, alpha=0.5)
        ax.set_ylabel('Difference [%], Percentile 95%', fontsize=14)
        ax.set_xlabel('Layer [km]', fontsize=14)    
        ax.set_xticks(ind)
        ax.set_xticklabels(xlabel, rotation=45)
        ax.tick_params(labelsize=14)
        #ax.set_xlabel('Site')
        #ax.axvline(0, color='black', lw=1)

        #plt.gca().invert_yaxis()

        fig.subplots_adjust(left=0.1, bottom=0.2)

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/Percentile_Apriori'+self.loc.upper()+'.pdf', bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #
        #--------------------------------
        # horizontal positions of the center of gravity of water vapor (HR-FTIR)
        #--------------------------------
        fig, ax1 = plt.subplots()

        ax1.plot(lonFTS2, latFTS2, 'k.', markersize=8)
        ax1.plot(lonFTS, latFTS, 'k.', markersize=20, color='r')
        ax1.annotate('{}'.format('HR-FTIR'), xy=(lonFTS, latFTS), color='r',  fontsize=14, ha='left')
        
        if self.loc.upper() == 'MLO':
            ax1.plot(-155.049, 19.717, 'k.', markersize=20, color='blue')
            Dist2  = mf.haversine(np.mean(lonFTS2), np.mean(latFTS2), -155.049, 19.717)
            Dist  = mf.haversine(lonFTS, latFTS, -155.049, 19.717)
            ax1.annotate('{}'.format('FPH'), xy=(-155.049, 19.717), color='blue',  fontsize=14, ha='left')
        else:
            ax1.plot(-105.1973, 39.949, 'k.', markersize=20, color='blue')
            Dist2  = mf.haversine(np.mean(lonFTS2), np.mean(latFTS2), -105.1973, 39.949)
            ax1.annotate('{}'.format('FPH'), xy=(-105.1973, 39.949), color='blue', fontsize=14, ha='left')
            Dist  = mf.haversine(lonFTS, latFTS, -105.1973, 39.949)

        print 'Mean Distance from the sonde launch to horizontal HR-FTIR = {0:.2f}'.format(Dist2)
        print 'Mean Distance from the sonde launch to HR-FTIR site= {0:.2f}'.format(Dist)
        
        ax1.set_ylabel('Latitude [$^{\circ}$]', fontsize=14)
        ax1.set_xlabel('Longitude [$^{\circ}$]', fontsize=14)        
        ax1.grid(True,which='both')
        #ax1.set_title('CFH profiles (2010 - 2016)', fontsize=16)
        ax1.tick_params(labelsize=14)

        plt.ticklabel_format(useOffset=False)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/Distance_'+self.loc.upper()+'.pdf', bbox_inches='tight')
            plt.close(fig)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


        #-------------------------------------------------
        # Profile position diference at specific Delta T
        #-------------------------------------------------
        if self.fleoutFlg:
            clmap     = 'jet'
            fig, ax   =plt.subplots(figsize=(7, 7))
            cm        = plt.get_cmap(clmap)
          
            norm = colors.Normalize(vmin=-15, vmax=30)

            #times    =  [15, 30, 60, 120, 240]
            times    =  [30, 60, 90, 120, 150]
            markers  =  ['o', '>', 'D', '*', 's', '<']

            for ti, t in enumerate(times):

                if (ti == 0) or (ti ==3):

                    inds = np.where(np.asarray(diffT) == t)[0]

                    #PrfDiff_i = np.asarray(PrfDiffRel[diffT[inds]])
            
                    #------------------------------
                    # Root Mean Square Error (RMSE)
                    #------------------------------
                    #SS_res = np.nansum( (Prfsonde_interp[diffT[inds]] - prfmean_day[diffT[inds]])**2, axis=0 )
                    #SS_res = np.asarray(SS_res)                
                    #rmse = np.sqrt( SS_res / len(Prfsonde_interp[diffT[inds]]) )

                    #------------------------------
                    # Bias
                    #------------------------------
                    biasCalc = np.nansum(prfmean_day[diffT[inds]] - Prfsonde_interp[diffT[inds]], axis=0 ) /Prfsonde_interp[diffT[inds]].shape[0] 
                    PrfDiff_i = biasCalc/np.nanmean(Prfsonde_interp[diffT[inds]], axis=0) * 100.


                    biasCalc_e = np.nansum(prfSTD_day[diffT[inds]], axis=0 ) /Prfsonde_interp[diffT[inds]].shape[0] 
                    PrfDiff_i_e = biasCalc_e/np.nanmean(Prfsonde_interp[diffT[inds]], axis=0) * 100.

                    Nx = prfmean_day[diffT[inds]].shape[0]
                    Ny = prfmean_day[diffT[inds]].shape[1]

                    Nx = int(np.sum(NobsFTS[diffT[inds]]))

                    PrfDist_i = np.asarray(PrfDist[diffT[inds]])

                    PrfDistMean  = np.nanmean(PrfDist_i, axis=0)
                    PrfDistStd   = np.nanstd(PrfDist_i, axis=0)
                    
                    #PrfDifftMean = np.nanmean(PrfDiff_i, axis=0)
                    #PrfDifftStd  = np.nanmean(PrfDiff_i, axis=0)

                    print '\nMean Relative difference for Delta T = {} minutes'.format(t)

                    for pn, pcol in enumerate(pCols):
                        
                        inds = np.where( (alt >= pcol[0]) & (alt <pcol[1])  )[0]
                        print 'Mean Horizontal DIfference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmean(PrfDistMean[inds]), np.nanmean(PrfDistStd[inds]), Nx)
                        print 'Median Horizontal DIfference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmedian(PrfDistMean[inds]), np.nanmean(PrfDistStd[inds]), Nx)
                        print 'Mean Bias Difference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmean(PrfDiff_i[inds]), np.nanstd(PrfDiff_i[inds]), Nx)
                        print 'Median Bias Difference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}\n'.format(pcol,  np.nanmedian(PrfDiff_i[inds]), np.nanstd(PrfDiff_i[inds]), Nx)

                    #ax.plot(PrfDistMean, alt, color='gray', alpha=0.5, zorder=1)
                    if ti == 0: lab = '0-'+str(t) + ' min (N = {})'.format(Nx)
                    else: lab = str(times[ti-1])+'-'+str(t) + ' min (N = {})'.format(Nx)
                    #sc1 = ax.scatter(PrfDistMean, alt, c=PrfDiff_i, s=50, cmap=cm, norm=norm, zorder=2, marker=markers[ti], label=lab)
                    #sc1 = ax.scatter(PrfDistMean, alt, s=50,  marker=markers[ti], label=lab)

                    ax.scatter(PrfDistMean, alt, color=clr[ti], s=60, label= lab)
                    ax.plot(PrfDistMean, alt,   color=clr[ti])


                    #print np.transpose([PrfDiff_i, alt])

                #for pi in range(Nx):
                #    #ax.plot(PrfDist_i[pi, :], alt, color='gray', alpha=0.5, zorder=1)
                #    sc1 = ax.scatter(PrfDist_i[pi, :], alt, c=PrfDiff_i[pi, :], s=30, cmap=cm, norm=norm, zorder=2)
                
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('$\Delta$x [km]', fontsize=14)
            ax.grid(True)
            ax.legend(prop={'size':11}, loc=2)
            ax.set_ylim(0, 20)

            #cax  = fig.add_axes([0.86, 0.1, 0.03, 0.8])
            
            #cbar = fig.colorbar(sc1, cax=cax,format='%3i')
            #cbar.set_label('Bias [%]', fontsize=13)
            #ax.set_title('Averaging Kernels Scale Factor')
            
            ax.tick_params(which='both',labelsize=14)       

            #fig.subplots_adjust(right=0.82) 

            if self.pdfsav:
                    self.pdfsav.savefig(fig,dpi=200)
                    plt.savefig('/data/iortega/Manuscripts/Fig/PrfDist_'+self.loc.upper()+'.pdf', bbox_inches='tight')
            else:           
                plt.show(block=False) 
                user_input = raw_input('Press any key to exit >>> ')
                sys.exit()  
            

        #-------------------------------------------------
        #Profile position diference at every coincidence time   
        #-------------------------------------------------
        if self.fleoutFlg:
            fig = plt.figure(figsize=(14,9))

            outer_grid = gridspec.GridSpec(3, 3, wspace=0.15, hspace=0.1)
                
            for di, df in enumerate(diffT):
                ax = plt.Subplot(fig, outer_grid[di])

                for pr in PrfDist[df]:
                    ax.plot(pr, alt,'k')

                ax.set_xlim(0, 200)
                ax.grid(True, color='gray', alpha=0.5)
                ax.tick_params(which='both',labelsize=10)
                ax.annotate('{} min'.format((df)), xy=(0.8, 0.9), xycoords='axes fraction', fontsize=12, ha='left')

                ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
          
                ax.xaxis.set_tick_params(which='minor',labelbottom='off')

                fig.add_subplot(ax)

            all_axes = fig.get_axes()
            #show only the outside spines
            for ax in all_axes:
                for sp in ax.spines.values():
                    sp.set_visible(False)
                    plt.setp(ax.get_xticklabels(), visible=False)
                if ax.is_first_row():
                    ax.spines['top'].set_visible(True)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                if ax.is_last_col():
                    ax.spines['right'].set_visible(True)

            if (9 % 2 == 1): #even
                all_axes[-2].spines['bottom'].set_visible(True)
                plt.setp(all_axes[-2].get_xticklabels(), visible=True)
                all_axes[-2].set_zorder(1)

            all_axes[-1].set_xlabel('$\Delta$x [km]')
            all_axes[-2].set_xlabel('$\Delta$x [km]')
            all_axes[-3].set_xlabel('$\Delta$x [km]')
            
            #fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
            fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

            if self.pdfsav:
                self.pdfsav.savefig(fig,dpi=200)
                
            else:
                plt.show(block=False)
                #user_input = raw_input('Press any key to exit >>> ')
                #sys.exit()

        #---------------------------------
        # Figure: Profile with the difference as function of Time
        #---------------------------------

        clr = mf.clrplt()

        #fig1, (ax1, ax2)  = plt.subplots(1, 2,figsize=(10,8), sharey=True)
        fig1, ax1  = plt.subplots(figsize=(6,7))
        
        for ii, df in enumerate(diffT):

            PrfDiffMean       = np.nanmean(PrfDiff[df], axis=0)
            PrfDiffAprMean    = np.nanmean(PrfDiffApr[df], axis=0)
            PrfDifRelfMean    = np.nanmean(PrfDiffRel[df], axis=0)
            PrfDiffAprRelMean = np.nanmean(PrfDiffAprRel[df], axis=0)

            #ax1.plot(PrfDiffMean/1e3, alt, '-', color=clr[ii], label=str(df) + '-'+str(len(NobsFTS[df])) + '-'+str(np.sum(NobsFTS[df])) , linewidth=3)
            ax1.plot(PrfDiffMean/1e3, alt, '-', color=clr[ii], label=str(df)+' min', linewidth=3)
            ax1.scatter(PrfDiffMean/1e3, alt, facecolors='white', color=clr[ii], s=45)

            #ax0[0, ndoi].scatter(prfmean_day[dfValue][d]/1e3,alt,facecolors='white', s=35, color='b')
            #ax0[0, ndoi].fill_betweenx(alt,prfmean_day[dfValue][d]/1e3-prferr_day[dfValue][d]/1e3,prfmean_day[dfValue][d]/1e3+prferr_day[dfValue][d]/1e3,alpha=0.25,color='blue')

            #ax2.plot(PrfDifRelfMean, alt, '-', color=clr[ii], linewidth=3)
            #ax1.plot(PrfDiffAprMean/1e3, alt, '--', color=clr[ii], linewidth=3)
            #ax2.plot(PrfDiffAprRelMean, alt, '--',  color=clr[ii], linewidth=3)
            
        if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
        if loc.lower() == 'mlo': ax1.set_ylim(3, 15)

        ax1.set_xlim(xmin=-2, xmax=1.5)
        ax1.grid(True,which='both')
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=14)   
        ax1.tick_params(which='both',labelsize=14)
        ax1.axvline(x=0, linestyle='--', linewidth=2, color='k')

        ax1.legend(prop={'size':14}, loc=2)

        xval = range(-2,3)
        for pn, pcol in enumerate(pCols):
            ax1.axhline(y=pcol[0], linewidth=2, color='gray', alpha=0.5)  
            ax1.axhline(y=pcol[1], linewidth=2, color='gray', alpha=0.5)
            ax1.fill_between(xval,pcol[0],pcol[1],alpha=0.15, color='0.5')  


        #ax2.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax2.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax2.set_ylim(3, 15)
        # ax2.set_xlim(-50, 100)
        # ax2.grid(True,which='both')
        # ax2.set_xlabel('Relative Difference [%]', fontsize=14)   
        # ax2.tick_params(which='both',labelsize=14)
        # ax2.axvline(x=0, linestyle='--', linewidth=2, color='k')

        if self.pdfsav: 
            self.pdfsav.savefig(fig1,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        #---------------------------------
        # Figure: Profile with the difference of all profiles for fixed time
        #---------------------------------
        # fig1, (ax1, ax2)  = plt.subplots(1, 2,figsize=(10,8), sharey=True)
        
        # PrfDiffMean       = np.nanmean(PrfDiff[dfValue], axis=0)
        # PrfDiffAprMean    = np.nanmean(PrfDiffApr[dfValue], axis=0)
        # PrfDifRelfMean    = np.nanmean(PrfDiffRel[dfValue], axis=0)
        # PrfDiffAprRelMean = np.nanmean(PrfDiffAprRel[dfValue], axis=0)


        # for i, diff in enumerate(PrfDiff[dfValue]):
        #     ax1.plot(diff/1e3, alt, '-', color='gray', linewidth=1, alpha=0.5)
        #     ax2.plot(PrfDiffRel[dfValue][i], alt, '-', color='gray', linewidth=1, alpha=0.5)
            

        # ax1.plot(PrfDiffMean/1e3, alt, '-', color='b', label='Mean (retrieved)', linewidth=3)
        # ax1.plot(PrfDiffAprMean/1e3, alt, '--', color='b',label='Mean (a priori)', linewidth=3)
        # #ax1.plot(tot_err_Mean, alt, '-', color='r', label='a priori', linewidth=3)
        # #ax1.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax1.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax1.set_ylim(3, 15)

        # #ax1.set_xlim(xmin=0)
        # ax1.grid(True,which='both')
        # ax1.set_ylabel('Altitude [km]', fontsize=14)
        # ax1.set_xlabel('VMR [x10$^3$ ppm$_v$]', fontsize=14)   
        # ax1.tick_params(which='both',labelsize=14)
        # ax1.axvline(x=0, linestyle='--', linewidth=2, color='k')

        # ax1.legend(prop={'size':14})

        # ax2.plot(PrfDifRelfMean, alt, '-',  color='b',  linewidth=3)
        # ax2.plot(PrfDiffAprRelMean, alt, '--',  color='b', linewidth=3)
        # #ax2.set_ylim(0, 15)
        # if loc.lower() == 'fl0': ax2.set_ylim(1, 15)
        # if loc.lower() == 'mlo': ax2.set_ylim(3, 15)
        # #ax2.set_xlim(-100, 100)
        # ax2.grid(True,which='both')
        # ax2.set_xlabel('Relative Difference [%]', fontsize=14)   
        # ax2.tick_params(which='both',labelsize=14)
        # ax2.axvline(x=0, linestyle='--', linewidth=2, color='k')

        # if self.pdfsav: 
        #     self.pdfsav.savefig(fig1,dpi=200)
        # else:
        #     plt.show(block=False)
        #     user_input = raw_input('Press any key to exit >>> ')
        #    sys.exit()

        #---------------------------------
        #
        #---------------------------------
        if self.loc == 'mlo': altHist  = pCols[0:]# [ [3, 5], [5, 7], [7, 9], [9, 11], [11, 13] ]
        else: altHist  = pCols[0:-1]# [[0, 3], [3, 5], [5, 7], [7, 9], [9,11] ]
        
        numBins = 20

        if self.loc == 'fl0': dfValue = 30
        elif self.loc == 'mlo': dfValue = 60

        
        #---------------------------------
        #
        #---------------------------------

        # fig = plt.figure(figsize=(10,7.5))

        # for df in diffT:

        #     if df == dfValue:

        #         #---------------------------------
        #         # Figure: bias as function of Distance
        #         #---------------------------------
        #         for p, h in enumerate(altHist):

        #             res = np.true_divide( (FTSvmrP[str(df)+'_'+str(p)] - sondevmrP[str(df)+'_'+str(p)]), sondevmrP[str(df)+'_'+str(p)])*100.
        #             #res = FTSsumP[str(df)+'_'+str(p)] - sondesumP[str(df)+'_'+str(p)]

        #             gs1 = gridspec.GridSpec(1, 3)

        #             if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
        #             if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
        #             if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                    
        #             if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
        #             if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
        #             if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

        #             ax1 = plt.subplot(gs1[0:1, :])

        #             weights = np.ones_like(x)/len(x)
        #             #plt.hist(myarray, weights=weights)
        #             #normed=True means that the total area under the histogram is equal to 1 but the sum of heights is not equal to 1  

        #             ax1.plot(distP[str(df)+'_'+str(p)], res, '-',  linewidth=3)
        #             ax1.scatter(distP[str(df)+'_'+str(p)], res, facecolors='white', s=45)  

                    
        #             ax1.grid(True,which='both')
        #             ax1.tick_params(which='both',labelsize=12, labelbottom='off')
        #             #ax1.set_ylabel('Probability', fontsize=12)
        #             ax1.set_title(str(h[0])+' - '+str(h[1])+' km',
        #                       horizontalalignment='center', verticalalignment='baseline', fontsize=14)
        #             ax1.set_xlim(0, 150)
        #             ax1.set_ylim(-100, 100)

        #         #------
        #             #ax1.set_xlabel('Relative Difference [%]', fontsize=12)
        #             #ax1.text(0.55,0.875,"Bias: {0:.3f}$\pm${1:.3}".format(mu, sigma),transform=ax1.transAxes,  fontsize=13)
                   
        #             if (p == 4) or (p == 5): 
        #                 ax1.set_xlabel('$\Delta$x [km]', fontsize=14)
        #                 ax1.tick_params(labelbottom='on')
        #             if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)

        
        #         if self.pdfsav: 
        #             self.pdfsav.savefig(fig,dpi=200)
        #             plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        #         else:
        #             plt.show(block=False)
                    #user_input = raw_input('Press any key to exit >>> ')
                    #sys.exit()
                    
        #------------------------------------
        # Plot time series of partial columns
        #------------------------------------
        #DateFmt      = DateFormatter('%m/%d/%Y')
        DateFmt      = DateFormatter('%Y')
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()

        #------------------------------------
        #
        #------------------------------------
        fig  = plt.figure(figsize=(10,7.5))

        outer_grid = gridspec.GridSpec(3, 2, wspace=0.15, hspace=0.15)

        for pn, pcol in enumerate(pCols[0:-1]):

            inds = np.where( (alt >= pcol[0]) & (alt <=pcol[1])  )[0]

            sondevmrP[str(pn)]    = np.ma.average(np.ma.masked_invalid(Prfsonde_interp[dfValue][:,inds]), axis=1,  weights=np.ma.masked_invalid(sondeairmass_a_interp[dfValue][:,inds]))
            sondevmrsdP[str(pn)]  = np.ma.average(np.ma.masked_invalid(Prfsonde_sd_interp[dfValue][:, inds]), axis=1, weights=np.ma.masked_invalid(sondeairmass_a_interp[dfValue][:,inds]))
            sondesumP[str(pn)]    = np.nansum(np.ma.masked_invalid(Prfsonde_interp[dfValue][:,inds])*np.ma.masked_invalid(sondeairmass_a_interp[dfValue][:,inds])*1e-6, axis=1)
            sondesumsdP[str(pn)]  = np.nansum(np.ma.masked_invalid(Prfsonde_sd_interp[dfValue][:,inds])*np.ma.masked_invalid(sondeairmass_a_interp[dfValue][:,inds])*1e-6, axis=1)
            
            FTSvmrP[str(pn)]      = np.ma.average(np.ma.masked_invalid(prfmean_day[dfValue][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[dfValue][:,inds]))
            FTSvmrsdP[str(pn)]    = np.ma.average(np.ma.masked_invalid(prferr_day[dfValue][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[dfValue][:,inds]))
            
            FTSAprvmrP[str(pn)]   = np.ma.average(np.ma.masked_invalid(aPrf_day[dfValue][:,inds]), axis=1,  weights=np.ma.masked_invalid(Airmassmean[dfValue][:,inds]))
            FTSsumP[str(pn)]      = np.nansum(np.ma.masked_invalid(rPrfMolmean_day[dfValue][:,inds]), axis=1)
            FTSsumsdP[str(pn)]    = np.nansum(np.ma.masked_invalid(prferrMol_day[dfValue][:,inds]), axis=1)
            FTSAprsumP[str(pn)]   = np.nansum(np.ma.masked_invalid(aPrfMol_day[dfValue][:,inds]), axis=1)



            ax1 = plt.Subplot(fig, outer_grid[pn])

            ax1.scatter(FTSdates[dfValue], FTSvmrP[str(pn)],facecolors='white', s=60, color='red', label='HR-FTIR')
            ax1.plot(FTSdates[dfValue], FTSvmrP[str(pn)], color='red')

            ax1.scatter(sondedates[dfValue], sondevmrP[str(pn)],facecolors='white', s=60, color='gray', label='FPH')
            ax1.plot(sondedates[dfValue], sondevmrP[str(pn)], color='gray')

            ax1.grid(True,which='both')
            ax1.tick_params(which='both',labelsize=12)

            ax1.grid(True)
            
            ax1.xaxis.set_major_formatter(DateFmt)
            ax1.xaxis.set_minor_locator(monthLc)
            
            
            if (p == 4) or (p == 5): 
                ax1.set_xlabel('Dates', fontsize=14)
                ax1.tick_params(labelbottom='on')
            if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Partial Column\n[molecules cm$^{-2}$]',multialignment='center', fontsize=14)

            ax1.annotate(str(pcol[0])+' - '+str(pcol[1])+' km',  xycoords='axes fraction', xy=(0.035, 0.875), color='k',  fontsize=14, ha='left')

            
            fig.add_subplot(ax1)
            fig.subplots_adjust(left=0.1, bottom=0.1, right=0.975, top=0.95)

            #sondelatP[str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[dfValue][:,inds]), axis=1)
            #sondelonP[str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[dfValue][:,inds]), axis=1)

            #distP[str(df)+'_'+str(pn)]        = [mf.haversine(np.nanmean(lonmean[df]), np.mean(latmean[df]), lo, la ) for (lo, la) in  zip(sondelonP[str(df)+'_'+str(pn)], sondelatP[str(df)+'_'+str(pn)])]
            #distP[str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[dfValue][:,inds]), axis=1)


        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            plt.savefig('/data/iortega/Manuscripts/Fig/HistLayers_'+self.loc.upper()+'.pdf', bbox_inches='tight')
        else:
            plt.show(block=False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()
        #fig1, ax1  = plt.subplots(len(pCols),1,figsize=(9,12), sharex=True)
        #fig2, ax2  = plt.subplots(len(pCols),1,figsize=(9,12), sharex=True)

        #fig3, ax3  = plt.subplots(len(pCols),1,figsize=(7,13))
        #fig4, ax4  = plt.subplots(len(pCols),1,figsize=(7,13))

        #print '\n'

        #for df in diffT:

        #    if df == dfValue:
       
        #        for pn, pcol in enumerate(pCols):

                    #-------------------------------------------
                    #
                    #-------------------------------------------
  

                    #ax3[pn].text(0.025,0.92,"slope: {0:.2f}+/-{1:.2}".format(slope, sd_slope),transform=ax3[pn].transAxes,  fontsize=11)
                    #ax3[pn].text(0.025,0.85,"Intercept: {0:.2f}+/-{1:.2}".format(intercept, sd_intercept),transform=ax3[pn].transAxes,  fontsize=11)
                    #ax3[pn].text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax3[pn].transAxes,  fontsize=11)        

            #         ax1[pn].scatter(FTSdates[df], FTSAprvmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')
            #         ax1[pn].plot(FTSdates[df], FTSAprvmrP[str(df)+'_'+str(pn)], color='gray')
            
            #         ax1[pn].errorbar(FTSdates[df],FTSvmrP[str(df)+'_'+str(pn)],yerr=FTSvmrsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='-', color='k', ecolor='k')
            #         ax1[pn].scatter(FTSdates[df], FTSvmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='k', label='Retrieved')
                    
            #         ax1[pn].errorbar(sondedates[df], sondevmrP[str(df)+'_'+str(pn)],yerr=sondevmrsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='r', ecolor='r')
            #         ax1[pn].scatter(sondedates[df], sondevmrP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='r', label='Sonde')
            #         ax1[0].legend(prop={'size':12}, loc=1)

            #         ax2[pn].scatter(FTSdates[df], FTSAprsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')
            #         ax2[pn].plot(FTSdates[df], FTSAprsumP[str(df)+'_'+str(pn)], color='gray')
                    
            #         ax2[pn].errorbar(FTSdates[df],FTSsumP[str(df)+'_'+str(pn)],yerr=FTSsumsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='k', ecolor='k', label='Retrieved')
            #         ax2[pn].scatter(FTSdates[df], FTSsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='k')
            #         ax2[pn].errorbar(sondedates[df], sondesumP[str(df)+'_'+str(pn)],yerr=sondesumsdP[str(df)+'_'+str(pn)],markersize=60, linestyle='-', color='r', ecolor='r', label='Sonde')
            #         ax2[pn].scatter(sondedates[df], sondesumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='r')
            #         ax2[0].legend(prop={'size':12}, loc=1)


            #         ax1[pn].grid(True)
            #         ax1[pn].set_ylabel('Weighted VMR ['+sclname+']', fontsize=14)
            #         ax1[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km',
            #                   multialignment='center',fontsize=14)
            #         ax1[pn].tick_params(which='both',labelsize=14)
            #         ax1[pn].xaxis.set_major_formatter(DateFmt)
            #         ax1[pn].xaxis.set_minor_locator(monthLc)


            #         ax2[pn].grid(True)
            #         ax2[pn].set_ylabel('Partial Column\n[molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax2[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km',
            #                   multialignment='center',fontsize=14)
            #         ax2[pn].tick_params(which='both',labelsize=14)
            #         ax2[pn].xaxis.set_major_formatter(DateFmt)
            #         ax2[pn].xaxis.set_minor_locator(monthLc)
                   
            #         ax1[pn].set_ylim(bottom=0)
            #         ax2[pn].set_ylim(bottom=0)
            #         fig1.autofmt_xdate()
            #         fig2.autofmt_xdate()

            #         fig1.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)
            #         fig2.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)

                    
            #         ax3[pn].errorbar(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)],xerr=sondevmrsdP[str(df)+'_'+str(pn)], yerr=FTSvmrsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='none', color='k', ecolor='k')
            #         ax3[pn].scatter(sondevmrP[str(df)+'_'+str(pn)], FTSvmrP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='k', label='Retrieved')
            #         ax3[pn].scatter(sondevmrP[str(df)+'_'+str(pn)], FTSAprvmrP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='gray', label='a priori')            
            #         ax3[0].legend(prop={'size':12}, loc=1)

            #         # # o2o = np.arange(0, 1e4, 50)

            #         # # ax1.plot(o2o, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # # ax1.plot(o2o*0.8, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # # ax1.plot(o2o, o2o*0.8, linestyle='--',linewidth=2,  color='gray')

                    
                    
            #         ax3[pn].text(0.65,0.2,"slope: {0:.2f}".format(slope),transform=ax3[pn].transAxes, color='gray', fontsize=14)
            #         ax3[pn].text(0.65,0.13,"Intercept: {:.2f}".format(intercept),transform=ax3[pn].transAxes, color='gray', fontsize=14)
            #         ax3[pn].text(0.65,0.06,"r-value: {0:.2f}".format(r_value),transform=ax3[pn].transAxes, color='gray', fontsize=14)

            #         ax4[pn].errorbar(sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)],yerr=FTSsumsdP[str(df)+'_'+str(pn)], xerr=sondesumsdP[str(df)+'_'+str(pn)], markersize=60, linestyle='none', color='k', ecolor='k')
            #         ax4[pn].scatter( sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)], facecolors='white', s=60, color='k', label='Retrieved')
            #         ax4[pn].scatter( sondesumP[str(df)+'_'+str(pn)], FTSAprsumP[str(df)+'_'+str(pn)],facecolors='white', s=60, color='gray', label='a priori')

            #         ax4[pn].legend(prop={'size':12}, loc=1)

            #         # o2o = np.arange(0, 1e23, 1e22)

            #         # ax2.plot(o2o, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # ax2.plot(o2o*0.8, o2o, linestyle='--',linewidth=2,  color='gray')
            #         # ax2.plot(o2o, o2o*0.8, linestyle='--',linewidth=2,  color='gray')


            #         slope, intercept, r_value, p_value, std_err = stats.linregress(sondesumP[str(df)+'_'+str(pn)], FTSsumP[str(df)+'_'+str(pn)])


            #         ax4[pn].text(0.025,0.92,"slope: {0:.2f}".format(slope),transform=ax4[pn].transAxes,  fontsize=14)
            #         ax4[pn].text(0.025,0.85,"Intercept: {:.2E}".format(intercept),transform=ax3[pn].transAxes,  fontsize=14)
            #         ax4[pn].text(0.025,0.78,"r-value: {0:.2f}".format(r_value),transform=ax4[pn].transAxes,  fontsize=14)

            #         slope, intercept, r_value, p_value, std_err = stats.linregress(sondesumP[str(df)+'_'+str(pn)], FTSAprsumP[str(df)+'_'+str(pn)])
                
            #         ax4[pn].text(0.65,0.2,"slope: {0:.2f}".format(slope),transform=ax4[pn].transAxes, color='gray', fontsize=14)
            #         ax4[pn].text(0.65,0.13,"Intercept: {:.2E}".format(intercept),transform=ax4[pn].transAxes, color='gray', fontsize=14)
            #         ax4[pn].text(0.65,0.06,"r-value: {0:.2f}".format(r_value),transform=ax4[pn].transAxes, color='gray', fontsize=14)

            #         ax3[pn].grid(True)
            #         ax3[pn].grid(True)
                
            #         ax3[pn].set_ylabel('Weighted VMR ['+sclname+']', fontsize=14)
            #         ax3[-1].set_xlabel('Weighted VMR - Sonde ['+sclname+']', fontsize=14)
            #         ax3[pn].tick_params(which='both',labelsize=14)
            #         ax3[pn].set_xlim(xmin=0)
            #         ax3[pn].set_ylim(ymin=0)
            #         ax3[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km', multialignment='center',fontsize=14)

            #         ax4[pn].grid(True)
            #         ax4[pn].grid(True)
                
            #         ax4[pn].set_ylabel('partial Column [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax4[-1].set_xlabel('Partial Column - sonde [molecules cm$^{-2}$]',multialignment='center', fontsize=14)
            #         ax4[pn].tick_params(which='both',labelsize=14)
            #         ax4[pn].set_xlim(xmin=0)
            #         ax4[pn].set_ylim(ymin=0)
            #         ax4[pn].set_title('Altitude Layer '+str(pcol[0])+' - '+str(pcol[1])+' km', multialignment='center',fontsize=14)


            #         fig3.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)
            #         fig4.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)

            # if self.pdfsav: 
            #     self.pdfsav.savefig(fig1,dpi=200)
            #     self.pdfsav.savefig(fig2,dpi=200)
            #     self.pdfsav.savefig(fig3,dpi=200)
            #     self.pdfsav.savefig(fig4,dpi=200)

            # else:
            #     plt.show(block=False)





        

       
      
    