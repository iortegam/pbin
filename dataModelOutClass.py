#----------------------------------------------------------------------------------------
# Name:
#        dataModelOutClasss.py
#
# Purpose:

# Notes:
#
#
# Version History:
#       Created, July, 2016  Ivan Ortega (iortega@ucar.edu)
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
#import statsmodels.api as sm
from scipy.integrate import simps

#MODIFIED (IVAN)
import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
#

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
from netCDF4 import Dataset

from mpl_toolkits.basemap import Basemap
import sfitClasses as sc
import glob

from scipy.interpolate import interp2d, interp1d

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

def findCls(dataArray, val):
    ''' Returns the indice and closest value in dataArray to val'''
    return np.argmin(abs(val-dataArray))

def findCls2(dataArray, val):
    return min(dataArray, key=lambda p:dist_sq(val, p))

def dist_sq(a, b): # distance squared (don't need the square root)
  return (a[0] - b[0])**2 + (a[1] - b[1])**2

def comb(a,b):
    c = []
    for i in a:
        for j in b:
            c.append(r_[i,j])
    return c

def closeFig(self):
    self.pdfsav.close()

def GasName(gas):

    if gas.lower() == 'c2h6':     gasSTR = 'C$_2$H$_6$'
    elif gas.lower() == 'o3':     gasSTR = 'O$_3$'
    elif gas.lower() == 'co':     gasSTR = 'CO'
    elif gas.lower() == 'hcho':   gasSTR = 'HCHO'
    elif gas.lower() == 'no2':    gasSTR = 'NO$_2$'
    elif gas.lower() == 'ch3oh':  gasSTR = 'CH$_3$OH'
    elif gas.lower() == 'c3h8':   gasSTR = 'C$_3$H$_8$'
    elif gas.lower() == 'pan':    gasSTR = 'PAN'
    elif gas.lower() == 'ch2o':   gasSTR = 'HCHO'
    elif gas.lower() == 'c2h2':   gasSTR = 'C$_2$H$_2$'
    elif gas.lower() == 'h2co':   gasSTR = 'HCHO'
    elif gas.lower() == 'ch4':    gasSTR = 'CH$_4$'
    elif gas.lower() == 'hcooh':  gasSTR = 'HCOOH'
    elif gas.lower() == 'nh3':    gasSTR = 'NH$_3$'
    else: gasSTR = gas

    return gasSTR

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
class WRFClass():
    '''
    This class deals with reading WRF-CHEM outputs.
    '''  
    def __init__(self,dataDir, file1, file2, outFname='', saveFlg= False):

        if (not dataDir) and (not file1) and (not file2):
            print 'Either Directory or Files of WRF-CHEM need to be specified'
            return False

        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'    

        self.WRFnc_f1   = dataDir + file1  # Your filename
        self.WRFnc_f2   = dataDir + file2  # Your filename

        ckDir(dataDir, exitFlg=True)

        ckFile(self.WRFnc_f1, exitFlg=True)
        ckFile(self.WRFnc_f1, exitFlg=True)

        self.WRF     = {}

        if saveFlg:
            self.WRFpdfsav = PdfPages(outFname)
            self.WRF['saveFlg'] = True
        else: self.WRF['saveFlg'] = False


        #---------------
        # ReadOutputData
        #---------------
    def ReadOutputWRF(self, Gases, pcol, sLat , sLon, calcWind=False):
        '''Function to Read out the WRF-CEHM net CDF files.'''

        self.WRF['calcWind'] = calcWind
        self.WRF['pcol'] = pcol

        nc_fid1 = Dataset(self.WRFnc_f1 , 'r')  # Dataset is the class behavior to open the file
        nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=True)
        
        nc_fid2 = Dataset(self.WRFnc_f2 , 'r')  # Dataset is the class behavior to open the file
        nc_attrs2, nc_dims2, nc_vars2 = mf.ncdump(nc_fid1, verb=True)

        self.WRF['Gases'] = Gases
    
        #-------------------------------------------------
        #Get variables from file 2
        #-------------------------------------------------
        cosalpha      = np.squeeze(np.transpose(nc_fid2.variables['COSALPHA'][:]))
        sinalpha      = np.squeeze(np.transpose(nc_fid2.variables['SINALPHA'][:]))
        xlon          = np.squeeze(np.transpose(nc_fid2.variables['XLONG'][:]))
        xlat          = np.squeeze(np.transpose(nc_fid2.variables['XLAT'][:]))

        LatLon = [zip(x,y) for x,y in zip(*(xlat,xlon))]
        LatLon = np.squeeze(np.asarray(LatLon))


        coord = (sLat, sLon)
        LatLon = np.reshape(LatLon, (9,2))

        LatLonclose  = findCls2(LatLon, coord)
        LatLon       = np.reshape(LatLon, (3,3,2))
        indsLoc      = np.where(LatLonclose == LatLon)[0]

        # #
        # print indsLoc
        # print LatLonclose

        # print sLat
        # print sLon

        # print xlon[indsLoc[0], indsLoc[1]]
        # print xlat[indsLoc[0], indsLoc[1]]
       
        # exit()

        #

        LatLon = np.reshape(LatLon, (9,2))
        coord = (np.max(xlat), np.max(xlon))
        LatLonclose  = findCls2(LatLon, coord)
        LatLon       = np.reshape(LatLon, (3,3,2))
        indsONG      = np.where(LatLonclose == LatLon)[0]

        LatLon = np.reshape(LatLon, (9,2))
        coord = (np.min(xlat), np.min(xlon))
        LatLonclose  = findCls2(LatLon, coord)
        LatLon       = np.reshape(LatLon, (3,3,2))
        indsBKG      = np.where(LatLonclose == LatLon)[0]
        
        print '\nWRF-Chem cell close to FL0:' + str(indsLoc)

        self.WRF['LatLon']     = LatLon
        self.WRF['xlon']       = xlon
        self.WRF['xlat']       = xlat
        self.WRF['indsLoc']    = indsLoc
        self.WRF['indsONG']    = indsONG
        self.WRF['indsBKG']    = indsBKG



        #-------------------------------------------------
        # Get variables from file 1 (Transposes to be consistent with IDL snippets from Gabi)
        #-------------------------------------------------
        time        = nc_fid1.variables['Times'][:]
        a1          = np.transpose(nc_fid1.variables['PB'][:])
        a2          = np.transpose(nc_fid1.variables['P'][:])   
        pressure    = (a1 + a2)/100.0
        a1          = np.transpose(nc_fid1.variables['HGT'][:])
        a2          = np.transpose(nc_fid1.variables['PHB'][:])
        a3          = np.transpose(nc_fid1.variables['PH'][:])
        Ttmp        = np.transpose(nc_fid1.variables['T'][:])
        p00         = nc_fid1.variables['P00'][:]

        ntimes      = len(time)
        nz          = pressure.shape[2] 
        ni          = pressure.shape[0]
        nj          = pressure.shape[1]

        print '\nNumber of times:'+ str(ntimes)
        print 'Number of layers:'+ str(nz)
        print 'Number of boxes:'+ str(ni)
        print 'Number of dimensions:' + str(pressure.shape)

        self.WRF['pressure']  = pressure

      
        #-------------------------------------------------
        # Calculating Date/Time
        #-------------------------------------------------
        dates = []
        for itime in time:
            yyyy =  int(''.join(itime[0:4])) 
            mm =  int(''.join(itime[5:7])) 
            dd =  int(''.join(itime[8:10])) 
            ho =  int(''.join(itime[11:13])) 
            mi =  int(''.join(itime[14:16]))
            se =  int(''.join(itime[17:19]))
            dates.append(dt.datetime(yyyy, mm, dd, ho, mi, se) )
        
        dates = np.asarray(dates)

        self.WRF['dates']  = dates
        
        #-------------------------------------------------
        # Calculating altitude in Km
        #-------------------------------------------------
        altitude    = np.zeros( (ni, nj, nz, ntimes)  )

        for itime in range(ntimes):
            for i in range(nz):
                atemp = (a2[:,:,i, itime]+a3[:,:,i,itime])/9.8/1000.
                atemp2 = (a2[:,:,i+1,itime]+a3[:,:,i+1,itime])/9.8/1000.
                altitude[:,:,i,itime]= (atemp+atemp2)/2.

        self.WRF['altitude']  = altitude

        #-------------------------------------------------
        # Calculating Temperature From Potential Temperature. Pressure, and p0 (Farenheit?)
        #-------------------------------------------------
        temperature = (Ttmp+300.)* ((np.true_divide(pressure,1000.0))**0.286)

        self.WRF['temperature']  = temperature

        #-------------------------------------------------
        # Calculating density profile, layer thickness and midpoint
        #-------------------------------------------------
        rho = (pressure*100.00) / ((temperature)*8.314) * 6.022e23 /100. /100. /100.  ##in molec/cm3
        self.WRF['rho']  = rho

        #-------------------------------------------------
        # Calculating layer thickness and midpoint. These values are 1 layer shorter
        #-------------------------------------------------

        #dz = np.zeros((ni, nj, nz, ntimes))
        #dz[:,:,:-1,:]     = np.absolute ((altitude[:,:,1:,:] - altitude[:,:,:-1,:]))*1000.*100.  #in cm
        dz     = np.absolute ((altitude[:,:,1:,:] - altitude[:,:,:-1,:]))*1000.*100.  #in cm
        self.WRF['dz']  = dz

        #midpoint = np.zeros((ni, nj, nz, ntimes))
        #midpoint[:,:,:-1,:]     = np.absolute ((altitude[:,:,1:,:] + altitude[:,:,:-1,:])/2.0)  #in km
        midpoint     = np.absolute ((altitude[:,:,1:,:] + altitude[:,:,:-1,:])/2.0)  #in km
        self.WRF['midpoints']  = midpoint

        #-------------------------------------------------
        # Calculating density values and air mass (same dimensions as dz, and midpoint)
        #-------------------------------------------------
        rho_interp = np.zeros( (ni, nj, nz-1, ntimes))
        for itime in range(ntimes):
            for ix in range(ni):
                for iy in range(nj):
                    rho_interp[ix, iy, :, itime] = np.interp(midpoint[ix, iy, :, itime], altitude[ix, iy, :, itime], rho[ix, iy, :, itime])

        airmass = rho_interp*dz 
        self.WRF['airmass']  = airmass
        #-------------------------------------------------
        # Getting Gases
        #-------------------------------------------------

        GasPrf = {}
        for g in Gases:
            GasPrf.setdefault(g, []).append(np.transpose(nc_fid1.variables[g][:]))
            GasPrf[g] = np.squeeze(np.asarray(GasPrf[g]))
            GasPrf[g] = np.array(GasPrf[g])*1000.    #ppb

            self.WRF.setdefault('GasPrf_'+g, []).append(np.transpose(nc_fid1.variables[g][:]))
            self.WRF['GasPrf_'+g] = np.squeeze(np.asarray(self.WRF['GasPrf_'+g]))
            self.WRF['GasPrf_'+g] = np.array(self.WRF['GasPrf_'+g])*1000.    #ppb


        #-------------------------------------------------
        # Calculating VMR weighted in total/partial columns
        #-------------------------------------------------
        GasWvmr  = {}  
        GasTC    = {}
        GasPTC   = {}   
        for g in Gases:
            Prf    = np.squeeze(np.asarray(GasPrf[g]))
            for itime in range(ntimes):
                for ix in range(ni):
                    for iy in range(nj):
                        ind1   = mf.nearestind(pcol[0], midpoint[ix, iy, :, itime])
                        ind2   = mf.nearestind(pcol[1], midpoint[ix, iy, :, itime])
                        alt1   = midpoint[ix, iy, ind1, itime]
                        alt2   = midpoint[ix, iy, ind2, itime]
                        Prf_i    = np.asarray(Prf[ix, iy,:, itime])                 
                        Prf_interp  = np.interp(midpoint[ix, iy, :, itime], altitude[ix, iy, :, itime], Prf_i)
                        self.WRF.setdefault('GasWvmr_'+g, []).append(np.average(Prf_interp[ind1:ind2], weights=airmass[ix, iy,ind1:ind2, itime]))
                        self.WRF.setdefault('GasTC_'+g, []).append(np.sum(Prf_interp[:]*airmass[ix, iy,:, itime]))
                        self.WRF.setdefault('GasPTC_'+g, []).append(np.sum(Prf_interp[ind1:ind2]*airmass[ix, iy,ind1:ind2, itime]))

            

            self.WRF['GasWvmr_'+g]   = np.asarray(self.WRF['GasWvmr_'+g])
            self.WRF['GasWvmr_'+g]   = np.reshape(self.WRF['GasWvmr_'+g],(ni, nj, ntimes), order='F') # Re-shape to boxes (tested), the F is FORTRAN indexing
            self.WRF['GasWvmr_'+g]   = np.array(self.WRF['GasWvmr_'+g])

            self.WRF['GasTC_'+g]     = np.asarray(self.WRF['GasTC_'+g])
            self.WRF['GasTC_'+g]     = np.reshape(self.WRF['GasTC_'+g],(ni, nj, ntimes), order='F') # Re-shape to boxes (tested), the F is FORTRAN indexing
            self.WRF['GasTC_'+g]     = np.array(self.WRF['GasTC_'+g])/1e9   #molec/cm2

            self.WRF['GasPTC_'+g]    = np.asarray(self.WRF['GasPTC_'+g])
            self.WRF['GasPTC_'+g]    = np.reshape(self.WRF['GasPTC_'+g],(ni, nj, ntimes), order='F') # Re-shape to boxes (tested), the F is FORTRAN indexing
            self.WRF['GasPTC_'+g]    = np.array(self.WRF['GasPTC_'+g])/1e9   #molec/cm2


        #-------------------------------------------------
        # Calculate Wind from Gabi IDL snippets
        #-------------------------------------------------
        if calcWind:
            print 'Rotate winds!'

            u  = np.transpose(nc_fid1.variables['U'][:])
            v  = np.transpose(nc_fid1.variables['V'][:])

            wd = np.zeros((ni,nj,nz,ntimes)) 
            ws = np.zeros((ni,nj,nz,ntimes))

            for itime in range(ntimes):
                utmp = np.zeros( (ni, nj, nz) )
                vtmp = np.zeros( (ni, nj, nz) )
                
                for iz in range(nz):
                    for i in range(ni):
                        #print i
                        for j in range(nj):
                            #print j
                            u2 = 0.5*(u[i,j,iz,itime]+u[i+1,j,iz,itime])
                            v2 = 0.5*(v[i,j,iz,itime]+v[i,j+1,iz,itime])
                            utmp[i,j, iz] =  u2*cosalpha[i,j] - v2*sinalpha[i,j]
                            vtmp[i,j, iz] =  v2*cosalpha[i,j] + u2*sinalpha[i,j]

                for iz in range(nz):
                    for ix in range(ni):
                        for iy in range(nj):

                            theta=0.
                            a1=utmp[ix,iy, iz]
                            a2=vtmp[ix,iy, iz]
                            if (a1 >= 0) & (a2 < 0): theta=360
                            if (a2 >= 0):  theta=180
                            tmp=np.arctan(a1/a2)*360/(2*np.pi)+theta
                            wd[ix,iy,iz,itime] = tmp
                            ws[ix,iy,iz,itime] = np.sqrt(a1*a1+a2*a2)

            self.WRF['wd'] = wd
            self.WRF['ws'] = wd

        
        dist = mf.haversine(xlon[1,2], xlat[1,2], xlon[1,1],xlat[1,1])

        

    def PltWRF(self):
        #-------------------------------------------------
        #                       Plots
        #-------------------------------------------------
        print '\nStarting WRF-CHEM plots'

        #---------------------------
        # Date locators for plotting
        #---------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()
        #DateFmt      = DateFormatter('%Y/%m/%d')
        mondays      = WeekdayLocator(MONDAY)
        DateFmt      = DateFormatter('%b %d')
        
        maxalt       = 20.0

        indsLoc = self.WRF['indsLoc']
        indsONG = self.WRF['indsONG']
        indsBKG = self.WRF['indsBKG']
        pcol = self.WRF['pcol']

        alt = self.WRF['altitude'][indsLoc[0],indsLoc[1],:,0]
        indsalt = np.where(alt <= maxalt)[0]
        alt = alt[indsalt]

        #-------------------------------------------------
        # Wind
        #-------------------------------------------------

        if self.WRF['calcWind']:

            fig0, (ax1, ax2) = plt.subplots(2, figsize=(10,7), sharex=True)
            ax1.plot(self.WRF['dates'], self.WRF['wd'][indsLoc[0],indsLoc[1],0,:],'k.',markersize=4)
            ax1.grid(True)
            ax1.set_ylabel('Wind direction',fontsize=14)
            ax1.set_title('Near-surface Wind direction\nWRF-Chem (FRAPPE - 2014)',
                      multialignment='center',fontsize=14)
            ax1.tick_params(labelsize=14)
            ax1.xaxis.set_major_locator(mondays)
            ax1.xaxis.set_minor_locator(dayLc)
            ax1.xaxis.set_major_formatter(DateFmt)

            ax2.plot(dates, ws[indsLoc[0],indsLoc[1],0,:],'k.',markersize=4)
            ax2.grid(True)
            ax2.set_ylabel('Wind Speed [m/s]',fontsize=14)
            ax2.set_xlabel('Date',fontsize=14)
            ax2.set_title('Near-surface Wind speed\nWRF-Chem (FRAPPE - 2014)',
                      multialignment='center',fontsize=14)
            ax2.tick_params(labelsize=14)
            ax2.xaxis.set_major_locator(mondays)
            ax2.xaxis.set_minor_locator(dayLc)
            ax2.xaxis.set_major_formatter(DateFmt)

            # ''' Wind Profile '''
            # #x = np.linspace(0,len(dates),len(dates)+1)
            # #y = np.linspace(0,alt[:,0], alt[:,0]+1)

            # U, V = np.meshgrid(wd[indsLoc[0],indsLoc[1],:,:], ws[indsLoc[0],indsLoc[1],:,:])
            # fig00, ax = plt.subplots(1, figsize=(10,7), sharex=True)
            # ax.barbs(dates, alt[:,], U, V)

            if  self.WRF['saveFlg']:
                self.WRFpdfsav.savefig(fig0,dpi=200)
                plt.close()
            else:
                plt.show(block=False)

        fig0, ax = plt.subplots(figsize=(10,7))

        m = Basemap(llcrnrlon=-105.3,llcrnrlat=40.0,urcrnrlon=-105.17,urcrnrlat=40.075,
                   projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        #N=1000
        N=500       
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=15)

        x, y = m(self.WRF['xlon'], self.WRF['xlat'])
        ax.plot(x,y,'.k', markersize=10)
        ax.plot(x[0,0],y[0,0],'.g', markersize=10)
        ax.plot(x[2,2],y[2,2],'.r', markersize=10)
        ax.plot(x[1,1],y[1,1],'.b', markersize=10)

        p = m.contourf(x, y, self.WRF['GasWvmr_c2h6'][:,:,N])

        cbar = fig0.colorbar(p, orientation='vertical', fraction=0.035, pad=0.05)
        cbar.set_label('C$_2$H$_6$ [ppb]', fontsize=14)

        parallels = np.arange(40.,40.1, 0.02)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 0.025)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')

        m.drawmapscale(-105.255, 40.07, -105.0, 40.0, 3, barstyle='fancy')
        ax.set_title('WRF-Chem Domain, example: '+ str(self.WRF['dates'][N]))

        if self.WRF['saveFlg']: self.WRFpdfsav.savefig(fig0,dpi=200)



        for gas in self.WRF['Gases']:

            gas2plt = gas

            gasSTR = GasName(gas2plt)

            #-------------------------------------------------
            # Time series of near-surface and weighted VMR
            #-------------------------------------------------

            fig, (ax, ax2) = plt.subplots(2, figsize=(10,7), sharex=True)
            
            ax.plot(self.WRF['dates'], self.WRF['GasPrf_'+gas2plt][indsLoc[0],indsLoc[1],0,:],'k.', markersize=8, label='NCAR-FL')
            ax.plot(self.WRF['dates'], self.WRF['GasPrf_'+gas2plt][indsONG[0],indsONG[1],0,:],'r.', markersize=8, label='North-East')
            ax.plot(self.WRF['dates'], self.WRF['GasPrf_'+gas2plt][indsBKG[0],indsBKG[1],0,:],'g.', markersize=8, label='South-West')

            ax.grid(True)
            ax.set_ylabel('VMR [ppb$_v$]',fontsize=14)
            ax.set_title(gasSTR + ' near-Surface VMR\nWRF-Chem (FRAPPE - 2014)',
                          multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)

            ax.legend(prop={'size':12})

            ax2.plot(self.WRF['dates'], self.WRF['GasWvmr_'+gas2plt][indsLoc[0],indsLoc[1],:],'k.', markersize=8, label='NCAR-FL')
            ax2.plot(self.WRF['dates'], self.WRF['GasWvmr_'+gas2plt][indsONG[0],indsONG[1],:],'r.', markersize=8, label='North-East')
            ax2.plot(self.WRF['dates'], self.WRF['GasWvmr_'+gas2plt][indsBKG[0],indsBKG[1],:],'g.', markersize=8, label='South-West')


            ax2.grid(True)
            ax2.set_ylabel('VMR [ppb$_v$]',fontsize=14)
            ax2.set_xlabel('Date',fontsize=14)
            ax2.set_title(gasSTR + ' VMR Weighted between '+str(pcol[0])+' to '+str(pcol[1])+' km',
                          multialignment='center',fontsize=14)
            ax2.tick_params(labelsize=14)
            ax2.xaxis.set_major_locator(mondays)
            ax2.xaxis.set_minor_locator(dayLc)
            ax2.xaxis.set_major_formatter(DateFmt)
            ax2.legend(prop={'size':12})
            #

            #-------------------------------------------------
            # Time series of Total/Partial Columns
            #-------------------------------------------------

            fig1, (ax,ax2) = plt.subplots(2, figsize=(10,7), sharex=True)
            ax.plot(self.WRF['dates'], self.WRF['GasTC_'+gas2plt][indsLoc[0],indsLoc[1],:],'k.', markersize=8, label='NCAR-FL')
            ax.plot(self.WRF['dates'], self.WRF['GasTC_'+gas2plt][indsONG[0],indsONG[1],:],'r.', markersize=8, label='North-East')
            ax.plot(self.WRF['dates'], self.WRF['GasTC_'+gas2plt][indsBKG[0],indsBKG[1],:],'g.', markersize=8, label='South-West')
            ax.grid(True)
            ax.set_ylabel('Total Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax.set_title(gasSTR + ' Total Column\nWRF-Chem (FRAPPE - 2014)',multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)
            ax.legend(prop={'size':12})

            ax2.plot(self.WRF['dates'], self.WRF['GasPTC_'+gas2plt][indsLoc[0],indsLoc[1],:],'k.', markersize=8, label='NCAR-FL')
            ax2.plot(self.WRF['dates'], self.WRF['GasPTC_'+gas2plt][indsONG[0],indsONG[1],:],'r.', markersize=8, label='North-East')
            ax2.plot(self.WRF['dates'], self.WRF['GasPTC_'+gas2plt][indsBKG[0],indsBKG[1],:],'g.', markersize=8, label='South-West')
            ax2.grid(True)
            ax2.set_ylabel('Partial Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax2.set_xlabel('Date',fontsize=14)
            ax2.set_title(gasSTR + ' Partial Column between '+str(pcol[0])+' to '+str(pcol[1])+' km',
                multialignment='center',fontsize=14)
            ax2.tick_params(labelsize=14)
            ax2.xaxis.set_major_locator(mondays)
            ax2.xaxis.set_minor_locator(dayLc)
            ax2.xaxis.set_major_formatter(DateFmt)




            #-------------------------------------------------
            # Diurnal profiles
            #-------------------------------------------------

            Prf2plt = self.WRF['GasPrf_'+gas2plt][indsLoc[0],indsLoc[1],:,:]
            Prf2ONG = self.WRF['GasPrf_'+gas2plt][indsONG[0],indsONG[1],:,:]
            Prf2BKG = self.WRF['GasPrf_'+gas2plt][indsBKG[0],indsBKG[1],:,:]
            
            Prf2plt_interp = np.zeros( (len(alt), len(self.WRF['dates'])))
            Prf2ONG_interp = np.zeros( (len(alt), len(self.WRF['dates'])))
            Prf2BKG_interp = np.zeros( (len(alt), len(self.WRF['dates'])))

            
            for itime in range(len(self.WRF['dates'])):
                indsalt = np.where(self.WRF['altitude'][indsLoc[0],indsLoc[1],:,itime] <= maxalt)[0]
                Prf2plt_interp[:,itime] = np.interp(alt, self.WRF['altitude'][indsLoc[0],indsLoc[1],indsalt,itime], Prf2plt[indsalt,itime])
                Prf2ONG_interp[:,itime] = np.interp(alt, self.WRF['altitude'][indsONG[0],indsONG[1],indsalt,itime], Prf2ONG[indsalt,itime])
                Prf2BKG_interp[:,itime] = np.interp(alt, self.WRF['altitude'][indsBKG[0],indsBKG[1],indsalt,itime], Prf2BKG[indsalt,itime])


            levels1 = np.arange(np.round(np.min(Prf2plt_interp), decimals=1),np.round(np.max(Prf2plt_interp),decimals=1))
            if len(levels1) < 3:
                 levels1 = np.arange(np.round(np.min(Prf2plt_interp), decimals=1)-0.1,np.round(np.max(Prf2plt_interp),decimals=1)+0.1,0.1)
          
            fig2,ax = plt.subplots(1, figsize=(10,6))
            cax1          = ax.contourf(self.WRF['dates'],alt,Prf2plt_interp,levels1,cmap=mplcm.jet)
            
            divider1      = make_axes_locatable(ax)
            cb1           = divider1.append_axes("right",size="5%",pad=0.04)
            cbar1         = plt.colorbar(cax1,cax=cb1)
            cbar1.ax.tick_params(labelsize=14)
            cbar1.ax.set_ylabel('VMR [ppb$_v$]', fontsize=14)
            
            ax.set_xlabel('Date', fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)        
            ax.yaxis.set_tick_params(which='major',labelsize=12)
            ax.xaxis.set_tick_params(which='major',labelsize=12)      
            ax.set_title('Diurnal Profile of ' + gasSTR+'\n' + 'WRF-Chem (FRAPPE - 2014)', fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)

            #-------------------------------------------------
            # Average Profiles
            #-------------------------------------------------

            Prfmean    = np.mean(Prf2plt_interp, axis=1)
            prfSTD     = np.std(Prf2plt_interp,axis=1)

            PrfmeanONG = np.mean(Prf2ONG_interp, axis=1)
            PrfmeanBKG = np.mean(Prf2BKG_interp, axis=1)
           


            fig3,ax  = plt.subplots(figsize=(7,9))
            ax.plot(Prfmean,alt, linewidth = 2.0, color='k', label='NCAR-FL')
            ax.fill_betweenx(alt,Prfmean-prfSTD,Prfmean+prfSTD,alpha=0.5,color='0.75')  
            ax.plot(PrfmeanONG,alt, linewidth = 2.0, color='r', label='North-East')
            ax.plot(PrfmeanBKG,alt, linewidth = 2.0, color='g', label='South-West')
           
            ax.set_title('Mean Profile of '+gasSTR+'\n' + 'WRF-Chem (FRAPPE - 2014)', fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            ax.grid(True,which='both')
            ax.tick_params(labelsize=14)
            ax.legend(prop={'size':12})

            #-------------------------------------------------
            # Anomalies
            #-------------------------------------------------

            daysall = [dt.date(da.year, da.month, da.day) for da in self.WRF['dates']]

            dailyVals        = mf.dailyAvg(self.WRF['GasWvmr_'+gas2plt][indsLoc[0],indsLoc[1],:], self.WRF['dates'], dateAxis=1, meanAxis=0)
            DT_d             = dailyVals['dates']
            y_d              = dailyVals['dailyAvg']
            deltaVMR = []

            for i, item in enumerate(DT_d):
                diff = np.asarray(daysall) - item
                inds = np.where( diff == dt.timedelta(0) )[0]

                deltaVMR.append(self.WRF['GasWvmr_'+gas2plt][indsLoc[0],indsLoc[1],inds] - y_d[i])
            deltaVMR = np.asarray(deltaVMR)
            dVMR = [y for x in deltaVMR for y in x]
            dVMR = np.asarray(dVMR)

            fig4, ax = plt.subplots(1, figsize=(10,7), sharex=True)
            
            ax.plot(self.WRF['dates'], dVMR,'k.', markersize=8, label='NCAR-FL')
            
            ax.grid(True)
            ax.set_ylabel('Delta - VMR [ppb$_v$]',fontsize=14)
            ax.set_title('Anomaly of '+gasSTR + ' in VMR\nWRF-Chem (FRAPPE - 2014)',
                          multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(dayLc)
            ax.xaxis.set_major_formatter(DateFmt)

            ax.legend(prop={'size':12})

            if self.WRF['saveFlg']:
                self.WRFpdfsav.savefig(fig,dpi=200)
                self.WRFpdfsav.savefig(fig1,dpi=200)
                self.WRFpdfsav.savefig(fig2,dpi=200)
                self.WRFpdfsav.savefig(fig3,dpi=200)
                plt.close()
                plt.clf()
            else:
                plt.show(block=False)
    
        if self.WRF['saveFlg']: self.WRFpdfsav.close()
        else:
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()


class CAMClass():
    '''
    This class deals with reading WRF-CHEM outputs.
    '''  
    def __init__(self,dataDir, file1, outFname='', saveFlg= False):

        if (not dataDir) and (not file1):
            print 'Either Directory or Files of CAM-CHEM need to be specified'
            return False

        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'    

        self.CAMnc_f1   = dataDir + file1  # Your filename

        ckDir(dataDir, exitFlg=True)
        ckFile(self.CAMnc_f1, exitFlg=True)


        self.CAM     = {}

        if saveFlg:
            self.CAMpdfsav = PdfPages(outFname)
            self.CAM['saveFlg'] = True
        else: self.CAM['saveFlg'] = False


        #---------------
        # ReadOutputData
        #---------------
    def ReadOutputCAM(self, Gases, pcol, sLat, sLon, interpFlg=False):
        '''Function to Read out the CAM-CHEM net CDF files.'''

        #-------------------------------------------------
        # CONSTANTS 
        #-------------------------------------------------
        NAv      = 6.0221415e+23                    #--- Avogadro's number
        g        = 9.81                             #--- m/s - gravity
        H        = (8.314*240.0)/(0.0289751*9.8)      #--- scale height
        MWair    = 28.94                            #--- g/mol
        xp_const = (NAv* 10.0)/(MWair*g)              #--- scaling factor for turning vmr into pcol

        self.CAM['pcol'] = pcol

        nc_fid1 = Dataset(self.CAMnc_f1 , 'r')  # Dataset is the class behavior to open the file
        nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=True)

        Gases = [g.upper() for g in Gases]
        
        self.CAM['Gases'] = Gases

        #-------------------------------------------------
        # First Method to find the closest pixel
        #-------------------------------------------------
        xlon          = np.squeeze(np.transpose(nc_fid1.variables['lon'][:]))
        xlat          = np.squeeze(np.transpose(nc_fid1.variables['lat'][:]))

        xlon = xlon - 360.0

        self.slat = sLat
        self.slon = sLon

        self.CAM['xlon']       = xlon
        self.CAM['xlat']       = xlat




        #-------------------------------------------------
        # Get variables from file 
        #-------------------------------------------------
        time      = nc_fid1.variables['time'][:] 
        p0        = nc_fid1.variables['P0'][:]
        pSurf     = nc_fid1.variables['PS'][:]
        hyai      = nc_fid1.variables['hyai'][:]
        hybi      = nc_fid1.variables['hybi'][:]

        lev       = nc_fid1.variables['lev'][:]
        dat       = nc_fid1.variables['date'][:]   

        ntimes      = len(time)
        nz          = lev.shape

        p0    = p0*0.01
        pSurf = pSurf*0.01

        print '\nNumber of times:'+ str(ntimes)
        print 'Number of dimensions:' + str(lev.shape)

        exit()

        # -------------------------------
        # Hybrid levels to pressure levels
        # -------------------------------
        pi = np.zeros((ntimes, len(hyai), len(xlat), len(xlon)))
        for i in range(len(hyai)):
            pi[:,i,:,:] = hyai[i]*p0 + hybi[i]*pSurf

        

        # -------------------------------
        # calculating dp
        # -------------------------------
        dp    = np.zeros((ntimes, len(lev), len(xlat), len(xlon)))
        pimid = np.zeros((ntimes, len(lev), len(xlat), len(xlon)))
        for i in range(len(lev)):
            dp[:,i,:,:] = pi[:,i+1,:,:] - pi[:,i,:,:]
            pimid[:,i,:,:] = (pi[:,i+1,:,:] + pi[:,i,:,:])*0.5


        #-------------------------------------------------
        # Calculating Date/Time
        #-------------------------------------------------
        dates = []
        for itime in dat:

            number_string = str(itime)
            
            yyyy =  int(''.join(number_string[0:4])) 
            mm =  int(''.join(number_string[4:6])) 
            dd =  int(''.join(number_string[6:8])) 

            dates.append(dt.date(yyyy, mm, dd) )
        
        dates = np.asarray(dates)

        self.CAM['dates']  = dates

        #-------------------------------------------------
        # Getting Gases
        #-------------------------------------------------

        GasPrf   = {}
        GasTC    = {}
        GasFlg   = {}
        for g in Gases:

            try:
                
                if g =='H2CO': 
                    GasPrf.setdefault(g, []).append(nc_fid1.variables['CH2O'][:])
               
                else: 
                    GasPrf.setdefault(g, []).append(nc_fid1.variables[g][:])
                
                GasPrf[g] = np.squeeze(np.asarray(GasPrf[g]))
                GasPrf[g] = np.array(GasPrf[g])

                GasTC[g]  = np.sum(GasPrf[g]*xp_const*dp, axis=1)

                self.CAM.setdefault('GasPrf_'+g.lower(), []).append(GasPrf[g])
               
                self.CAM['GasPrf_'+g.lower()] = np.squeeze(np.asarray(self.CAM['GasPrf_'+g.lower()]))
                self.CAM['GasPrf_'+g.lower()] = np.array(self.CAM['GasPrf_'+g.lower()])

                self.CAM['GasTC_'+g.lower()]  = np.sum(GasPrf[g]*xp_const*dp, axis=1)
                self.CAM['AIRMASS_'+g.lower()] = np.asarray(xp_const*dp)

                self.CAM['GasFlg_'+g.lower()] = True

            except:
                self.CAM['GasFlg_'+g.lower()] = False
                continue



        #-------------------------------------------------
        # Calculating altitude in Km. In absence of Temp we use the following equation 
        # (http://www.srh.noaa.gov/images/epz/wxcalc/pressureAltitude.pdf)
        #-------------------------------------------------
        altitude = ( 1.0 - (pi/p0)**0.190284   ) * 145366.45
        altitude = altitude * 0.3048 /1000.0

        self.CAM['altitude'] = altitude

        midpoints = np.zeros((ntimes, len(lev), len(xlat), len(xlon)))
        dz        = np.zeros((ntimes, len(lev), len(xlat), len(xlon)))
        pressure  = np.zeros((ntimes, len(lev), len(xlat), len(xlon)))

        for i in range(len(lev)):
            midpoints[:,i,:,:] = (altitude[:,i+1,:,:] + altitude[:,i,:,:]) * 0.5
            dz[:,i,:,:]        = (altitude[:,i,:,:] - altitude[:,i+1,:,:])
            pressure[:,i,:,:]  = (pi[:,i+1,:,:] + pi[:,i,:,:]) * 0.5

        self.CAM['midpoints'] = midpoints
        self.CAM['pressure']  = pressure

        
        if interpFlg:
            print 'Interpolating --> in construction, for now exit!'
            exit()

            #-------------------------------------------------
            # Getting Gases
            #-------------------------------------------------

            #dp_interp  = interp2d(xlon[:],xlat[:], dp,kind='linear', bounds_error=True)(sLon,sLat)

            #print dp.shape
            #print dp_interp.shape
            #exit()

            GasPrf   = {}
            
            for g in Gases:

                try:
                    
                    if g =='H2CO': 
                        GasPrf.setdefault(g, []).append(nc_fid1.variables['CH2O'][:])
                   
                    else: 
                        GasPrf.setdefault(g, []).append(nc_fid1.variables[g][:])
                    
                    GasPrf[g] = np.squeeze(np.asarray(GasPrf[g]))

                    GasPrf[g] = np.array(GasPrf[g])

                    Prf2dinterp = interp2d(xlon[:],xlat[:],GasPrf[g],kind='linear', bounds_error=True)(sLon,sLat)


                    self.CAM.setdefault('GasPrf_site_'+g.lower(), []).append(Prf2dinterp)
                   
                    self.CAM['GasPrf_site_'+g.lower()] = np.squeeze(np.asarray(self.CAM['GasPrf_site_'+g.lower()]))
                    self.CAM['GasPrf_site_'+g.lower()] = np.array(self.CAM['GasPrf_site_'+g.lower()])


                    self.CAM['GasTC_site_'+g.lower()] = np.sum(Prf2dinterp*xp_const*dp_interp, axis=1)
                    self.CAM['AIRMASS__site'+g.lower()] = np.asarray(xp_const*dp_interp)

                    self.CAM['GasFlg_site_'+g.lower()] = True

                except:
                    self.CAM['GasFlg_site_'+g.lower()] = False
                    continue


        else:

            latind = findCls(xlat[:],sLat)
            lonind = findCls(xlon[:],sLon)

            indsLoc = [latind, lonind]


            print 'CAM-Chem: Latitude of %s and Longitude of %s close to FL0:'%(str(xlat[latind]), str(xlon[lonind]) ) 
            distance = mf.haversine(sLon, sLat, xlon[lonind],  xlat[latind])
            print 'Distance from FL0 to the closest Pixel %s Km:'%(str(distance) ) 
            

            #-------------------------------------------------
            # Second Method to find the closest pixel
            #-------------------------------------------------
            xlon2, xlat2 = np.meshgrid(xlon, xlat)
            
            xlon2 = np.asarray(xlon2)
            xlat2 = np.asarray(xlat2)

            distance = np.zeros((xlon2.shape))

            for xi, x in enumerate(xlon):
                for yi, y in enumerate(xlat):

                    distance[yi, xi] = mf.haversine(sLon, sLat, x,  y)

            for xi, x in enumerate(xlon):
                for yi, y in enumerate(xlat):

                    if distance[yi, xi] == np.min(distance):
                        latind = yi
                        lonind = xi

                        print 'CAM-Chem: Latitude of %s and Longitude of %s close to FL0:'%(str(xlat[latind]), str(xlon[lonind]) ) 
                        print 'Distance from FL0 to the closest Pixel %s Km:'%(str(distance[yi, xi]) )

            #-------------------------------------------------
            #
            #-------------------------------------------------
            
            self.CAM['indsLoc']    = indsLoc

          
            


        # dpmean  = np.mean(pimid[:,:,indsLoc[0], indsLoc[1]], axis=0)
        # altCAM = midpoints[0,:, indsLoc[0], indsLoc[1]]

        # fig,ax  = plt.subplots(figsize=(7,9))
        # ax.plot(dpmean ,altCAM, linewidth = 2.0, color='k')
        # ax.scatter(dpmean,altCAM, facecolors='white', s=60, color='k')
    
        # ax.set_ylabel('Altitude [km]', fontsize=14)
        # ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
        # ax.grid(True,which='both')
        # ax.tick_params(labelsize=14)

           
        # plt.show(block=False)


        # Prfmairmass_1  = np.mean(self.CAM['AIRMASS_co'][:,:,indsLoc[0], indsLoc[1]], axis=0)
        # Prfmairmass_1  = np.mean(dz[:,:,indsLoc[0], indsLoc[1]], axis=0)
        
        # altCAM = midpoints[0,:, indsLoc[0], indsLoc[1]]

        # fig,ax  = plt.subplots(figsize=(7,9))
        # ax.plot(Prfmairmass_1 ,altCAM, linewidth = 2.0, color='k')
        # ax.scatter(Prfmairmass_1,altCAM, facecolors='white', s=60, color='k')
    
        # ax.set_ylabel('Altitude [km]', fontsize=14)
        # ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
        # ax.grid(True,which='both')
        # ax.tick_params(labelsize=14)

           
        # plt.show(block=False)
        # user_input = raw_input('Press any key to exit >>> ')
        # sys.exit()
        
        
    def PltCAM(self):
        #-------------------------------------------------
        #                       Plots
        #-------------------------------------------------
        print '\nStarting CAM-CHEM plots'

        #---------------------------
        # Date locators for plotting
        #---------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()
        
        mondays      = WeekdayLocator(MONDAY)
        #DateFmt      = DateFormatter('%b %d')
        DateFmt      = DateFormatter('%Y')
        
        #---------------------------
        # Defining variable
        #---------------------------
        maxalt       = 20.0
        indsLoc      = self.CAM['indsLoc']
        pcol         = self.CAM['pcol']
        alt          = self.CAM['midpoints'][0, :, indsLoc[0],indsLoc[1]]
        indsalt      = np.where(alt <= maxalt)[0]
        alt          = alt[indsalt]


        



        #for gas in self.CAM['Gases']:

        #    gasSTR = GasName(gas)

            

       
        for gas in self.CAM['Gases']:

            gas2plt = gas.lower()

            gasSTR = GasName(gas2plt)
            #-------------------------------------------------
            #
            #-------------------------------------------------

            #---------------------------------
            # Defining variables (CAM-CHEM) - Data is monthly mean already
            #---------------------------------

            Dates         = np.asarray(self.CAM['dates'])
            PrfCAM        = np.asarray(self.CAM['GasPrf_'+gas.lower()][:,:,self.CAM['indsLoc'][0],self.CAM['indsLoc'][1]])*1e9
            AirmassCAM    = np.asarray(self.CAM['AIRMASS_'+gas.lower()][:,:,self.CAM['indsLoc'][0],self.CAM['indsLoc'][1]])
            pressCAM      = np.asarray(self.CAM['pressure'][:,:,self.CAM['indsLoc'][0],self.CAM['indsLoc'][1]])
            TotCol        = np.asarray(self.CAM['GasTC_'+gas.lower()][:,self.CAM['indsLoc'][0],self.CAM['indsLoc'][1]])

            mnthlyVals            = mf.mnthlyAvg(PrfCAM, Dates, dateAxis=0, meanAxis=0)
            Prf_mnth              = mnthlyVals['mnthlyAvg']
            Prf_mnth_std          = mnthlyVals['std']
            Dates_mnth            = mnthlyVals['dates']

            mnthlyVals            = mf.mnthlyAvg(TotCol, Dates, dateAxis=0, meanAxis=0)
            TotCol_mnth           = mnthlyVals['mnthlyAvg']

            data = np.mean(self.CAM['GasPrf_'+gas.lower()][:, -1, :, :], axis=0)*1e9


            xmin = self.CAM['dates'][0]
            xmax = self.CAM['dates'][-1]

            #-------------------------------------------------
            #
            #-------------------------------------------------

            fig, ax = plt.subplots(figsize=(10,6))

            #m = Basemap(llcrnrlon=-160,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            #           projection='lcc',lat_1=33,lat_2=45,lon_0=-95, resolution='l')

            m = Basemap(llcrnrlat=10,urcrnrlat=60, llcrnrlon=-130,urcrnrlon=-55, rsphere=(6378137.00,6356752.3142),
                resolution='l',area_thresh=1000.,projection='lcc', lat_1=10,lon_0=-55)

            xfts, yfts = m(self.slon, self.slat)
            
            m.drawcoastlines(color='black')
            m.drawcountries(color='lightgrey')
            #m.fillcontinents(color='gray')
            m.drawstates(color='gray')

            m.drawparallels(np.arange(-80., 81., 5.), labels=[1,0,0,0], alpha=0.5)
            m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1],alpha=0.5)

            ax.plot(xfts, yfts, 'r^', markersize=5)

            xlon2, xlat2 = np.meshgrid(self.CAM['xlon'], self.CAM['xlat'])
            x, y = m(xlon2, xlat2)

            ax.plot(x,y,'.k', markersize=2.5, alpha=0.5)

            cs = m.contourf(x,y, data, cmap=mplcm.jet,extend='both') #,alpha=0.5)
            cs.axis='tight'

            cbar = fig.colorbar(cs, orientation='vertical', fraction=0.035, pad=0.05)
            cbar.set_label(gasSTR +' [ppb]', fontsize=14)

            ax.set_title('Near-surface vmr in CAM-Chem - {}\n'.format(gasSTR) + str(xmin)+' - ' +str(xmax), fontsize=14)
            #ax.set_title('Near-surface vmr in CAM-Chem - {}\n'.format(self.CAM['midpoints'][0,-1], gasSTR) + str(xmin)+' - ' +str(xmax), fontsize=14)

            if self.CAM['saveFlg']:
                self.CAMpdfsav.savefig(fig,dpi=200)

            #-------------------------------------------------
            #
            #-------------------------------------------------

            fig2, ax = plt.subplots(figsize=(10,8))

            #m = Basemap(llcrnrlon=-160,llcrnrlat=22,urcrnrlon=-64,urcrnrlat=49,
            #           projection='lcc',lat_1=33,lat_2=45,lon_0=-95, resolution='l')


            m = Basemap(llcrnrlat=30,urcrnrlat=45, llcrnrlon=-115,urcrnrlon=-90, rsphere=(6378137.00,6356752.3142),
                resolution='l',area_thresh=1000.,projection='lcc', lat_1=30,lon_0=-115)

            xfts, yfts = m(self.slon, self.slat)


            xclose, yclose = m(self.CAM['xlon'][self.CAM['indsLoc'][1]],  self.CAM['xlat'][self.CAM['indsLoc'][0]])

            
            m.drawcoastlines(color='black')
            m.drawcountries(color='lightgrey')
            #m.fillcontinents(color='gray')
            m.drawstates(color='gray')

            #m.drawparallels(np.arange(-80., 81., 5.), labels=[1,0,0,0], alpha=0.5)
            #m.drawmeridians(np.arange(-180., 181., 10.), labels=[0,0,0,1],alpha=0.5)

            ax.plot(xfts, yfts, 'r^', markersize=10)
            ax.plot(xclose, yclose, 'ko', markersize=10)


            xlon2, xlat2 = np.meshgrid(self.CAM['xlon'], self.CAM['xlat'])
            x, y = m(xlon2, xlat2)

            ax.plot(x,y,'.k', markersize=5, alpha=0.5)

            cs = m.contourf(x,y,data, cmap=mplcm.jet,extend='both') #,alpha=0.5)
            cs.axis='tight'

            cbar = fig2.colorbar(cs, orientation='vertical', fraction=0.035, pad=0.05)
            cbar.set_label(gasSTR +' [ppb]', fontsize=14)

            ax.set_title('Near-surface vmr in CAM-Chem - {}\n'.format(gasSTR) + str(xmin)+' - ' +str(xmax), fontsize=14)


            if self.CAM['saveFlg']:
                self.CAMpdfsav.savefig(fig2,dpi=200)


            #-------------------------------------------------
            # Time series of near-surface
            #-------------------------------------------------

            fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)

            ax.plot(Dates_mnth, Prf_mnth[:,-1]*1e9, color='k')
            ax.scatter(Dates_mnth, Prf_mnth[:,-1]*1e9, facecolors='white', s=60, color='k')
            
            ax.grid(True)
            ax.set_ylabel('VMR [ppb$_v$]',fontsize=14)
            ax.set_title(gasSTR + ' near-Surface VMR [Monthly average]\nCAM-Chem',
                          multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            #ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)


            #-------------------------------------------------
            # Time series of Total Columns
            #-------------------------------------------------
            dateYearFrac = mf.toYearFraction(Dates_mnth)
            weights      = np.ones_like(dateYearFrac)
            res          = mf.fit_driftfourier(dateYearFrac, TotCol_mnth, weights, 2)
            f_drift, f_fourier, f_driftfourier = res[3:6]

            fig1, ax = plt.subplots(1, figsize=(10,6), sharex=True)
            ax.plot(Dates_mnth, TotCol_mnth, color='k')
            ax.scatter(Dates_mnth, TotCol_mnth, facecolors='white', s=60, color='k')
            ax.plot(Dates_mnth,f_drift(dateYearFrac),label='Fitted Anual Trend')
            ax.plot(Dates_mnth,f_driftfourier(dateYearFrac),label='Fitted Anual Trend + intra-annual variability')
            ax.grid(True)
            ax.set_ylabel('Total Column\n[molecules$\cdot$cm$^{-2}$]',fontsize=14)
            ax.set_title(gasSTR + ' Total Column [Monthly average]\nCAM-Chem',multialignment='center',fontsize=14)
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Year', fontsize=14)
            #ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)

            ax.text(0.02,0.94,"Fitted trend -- slope: {0:.3E} ({1:.3f}%)".format(res[1],res[1]/np.mean(TotCol_mnth)*100.0),transform=ax.transAxes)
            ax.text(0.02,0.9,"Fitted intercept at xmin: {:.3E}".format(res[0]),transform=ax.transAxes)
            ax.text(0.02,0.86,"STD of residuals: {0:.3E} ({1:.3f}%)".format(res[6],res[6]/np.mean(TotCol_mnth)*100.0),transform=ax.transAxes) 

            #-------------------------------------------------
            # Diurnal profiles
            #-------------------------------------------------

            Prf2plt = Prf_mnth

            Prf2plt_interp = np.zeros( (len(alt), len(Dates_mnth) ) )
            
            for itime in range(len(Dates_mnth)):
                indsalt = np.where(self.CAM['midpoints'][itime,:,indsLoc[0],indsLoc[1]] <= maxalt)[0]
                #Prf2plt_interp[:, itime] = np.interp(alt, self.CAM['midpoints'][itime, indsalt,indsLoc[0],indsLoc[1]], Prf2plt[itime, indsalt], right=Prf2plt[itime, indsalt[-1]])
                Prf2plt_interp[:, itime] = np.transpose(Prf2plt[itime, indsalt])

            #print Prf2plt_interp.shape
            #print np.min(Prf2plt_interp), np.max(Prf2plt_interp)
            #print gasSTR

            levels1 = np.arange(np.round(np.min(Prf2plt_interp), decimals=1)-0.1,np.round(np.max(Prf2plt_interp),decimals=1)+0.1,0.1)
            #levels1 = np.arange(np.round(np.min(Prf2plt_interp)) ,np.round(np.max(Prf2plt_interp)) )
            
          
            fig2,ax = plt.subplots(1, figsize=(10,6))
            cax1          = ax.contourf(Dates_mnth,alt,Prf2plt_interp, levels1,cmap=mplcm.jet)
            
            divider1      = make_axes_locatable(ax)
            cb1           = divider1.append_axes("right",size="5%",pad=0.04)
            cbar1         = plt.colorbar(cax1,cax=cb1)
            cbar1.ax.tick_params(labelsize=14)
            cbar1.ax.set_ylabel('VMR [ppb$_v$]', fontsize=14)
            
            ax.set_xlabel('Year', fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)        
            ax.yaxis.set_tick_params(which='major',labelsize=12)
            ax.xaxis.set_tick_params(which='major',labelsize=12)      
            ax.set_title('Diurnal Profile of ' + gasSTR+'\n' + 'CAM-Chem', fontsize=14)
            ax.tick_params(labelsize=14)
            #ax.xaxis.set_major_locator(mondays)
            ax.xaxis.set_minor_locator(monthLc)
            ax.xaxis.set_major_formatter(DateFmt)

            #-------------------------------------------------
            # Average Profiles
            #-------------------------------------------------

            Prfmean = np.mean(Prf2plt_interp, axis=1)
            prfSTD  = np.std(Prf2plt_interp, axis=1)

            fig3,ax  = plt.subplots(figsize=(7,9))
            ax.plot(Prfmean,alt, linewidth = 2.0, color='k')
            ax.fill_betweenx(alt,Prfmean-prfSTD,Prfmean+prfSTD,alpha=0.5,color='0.75')  
            ax.set_title('Mean Profile of '+gasSTR+'\n' + 'CAM-Chem', fontsize=14)
            ax.set_ylabel('Altitude [km]', fontsize=14)
            ax.set_xlabel('VMR [ppb$_v$]', fontsize=14)    
            ax.grid(True,which='both')
            ax.tick_params(labelsize=14)

            #--------
            # Monthly
            #--------
            #----------------------------
            # Plot total columns by month
            #----------------------------
            month    = np.array([d.month for d in Dates_mnth])
            mnthSort = list(set(month))
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))
            
            for i,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[i] = np.mean(TotCol_mnth[inds])
                mnthSTD[i]  = np.std(TotCol_mnth[inds])   
                
            fig4,ax1  = plt.subplots(figsize=(10,6))
            ax1.plot(mnthSort,mnthMean,'k.',markersize=12)
            ax1.errorbar(mnthSort,mnthMean,yerr=mnthSTD,fmt='k.',markersize=12, ecolor='red')     
            ax1.grid(True,which='both')
            ax1.set_ylabel('Total Column\n[molecules cm$^{-2}$]', fontsize=14, multialignment='center')
            ax1.set_xlabel('Month', fontsize=14)
            ax1.set_title(gasSTR +' Monthly Mean with Standard Deviation\nCAM-Chem', fontsize=14)
            ax1.set_xlim((0,13))
            ax1.set_xticks(range(1,13))
            ax.tick_params(labelsize=14)     


            if self.CAM['saveFlg']:
                self.CAMpdfsav.savefig(fig,dpi=200)
                self.CAMpdfsav.savefig(fig1,dpi=200)
                self.CAMpdfsav.savefig(fig2,dpi=200)
                self.CAMpdfsav.savefig(fig3,dpi=200)
                self.CAMpdfsav.savefig(fig4,dpi=200)
                plt.close()
                plt.clf()
            else:
                plt.show(block=False)
    
        if self.CAM['saveFlg']: self.CAMpdfsav.close()
        else:
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()


class FLXClass(_DateRange):
    '''
    This class deals with reading FLEXPART outputs during FRAPPE.
    '''
    def __init__(self,dataDir, iyear=False,imnth=False,iday=False,fyear=False,fmnth=False,fday=False,incr=1, outFname='', saveFlg= False):
        
        #---------------------------------
        # Test if date range given. If so 
        # create a list of directories to 
        # read data from
        #---------------------------------
        if all([iyear,imnth,iday,fyear,fmnth,fday]): 
            
            # check if '/' is included at end of path
            if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'           
                
            # Check if directory exits
            ckDir(dataDir,exitFlg=True)
            
            self.endtime      = []
            self.starttime    = []
            self.dirFlg       = True
            self.dirLst       = []

            #---------------------
            # Establish date range
            #---------------------
            _DateRange.__init__(self, iyear, imnth, iday, fyear, fmnth, fday, incr=1)

            self.FLX     = {}

            #--------------------------------------------
            # Walk through first level of directories and
            # collect directory names for processing
            #--------------------------------------------
            for drs in os.walk(dataDir).next()[1]:

                if drs[0:8] == 'FLEXPART':

                    if _DateRange.inRange(self, int(drs[24:28]), int(drs[28:30]), int(drs[30:32]) ):

                        self.endtime.append(dt.datetime(int(drs[13:17]), int(drs[17:19]), int(drs[19:21]), int(drs[21:23])))
                        self.starttime.append(dt.datetime(int(drs[24:28]), int(drs[28:30]), int(drs[30:32]), int(drs[32:34])))
                        self.dirLst.append(dataDir+drs+'/output/')     
            
            #----------------------------------------------------------------------------
            # This is important in keeping all the gathered data in order (i.e. profiles,
            # summary, etc.). The sorted directory is carried through and dates are
            # collected in the summary and the profile gather
            #----------------------------------------------------------------------------
            self.dirLst.sort()
            self.starttime.sort()
            self.endtime.sort()

            self.starttime = np.asarray(self.starttime)
            self.endtime   = np.asarray(self.endtime)

            if saveFlg:
                #self.FLXpdfsav = PdfPages(outFname)
                self.FLX['saveFlg']   = True
            else: self.FLX['saveFlg'] = False

        #---------------
        # ReadOutputData
        #---------------
    def ReadOutputFLX(self, pltFLX=False):
        '''Function to Read out the FLX-CEHM net CDF files.'''

        self.FLX['SC'] = {}

        #-----------------------------------
        # Loop through collected directories
        #-----------------------------------
        for i, sngDir in enumerate(self.dirLst):

            if self.FLX['saveFlg']: self.FLXpdfsav = PdfPages(sngDir+'flx.pdf')

            nc_fid1 = Dataset(sngDir+'header_d01.nc' , 'r') 
            nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=False)

            #self.FLX['starttime'].append()

            self.FLX['xlon']         = np.asarray(nc_fid1.variables['XLONG'])
            self.FLX['xlat']        = np.asarray(nc_fid1.variables['XLAT'])
            self.FLX['ztop']         = np.asarray(nc_fid1.variables['ZTOP'])
            #self.species      = np.squeeze(np.asarray(nc_fid1.variables['SPECIES'])  )
            self.FLX['species']      = np.array(nc_fid1.variables['SPECIES'])
            self.FLX['species']      = np.array(self.FLX['species'] , dtype='c')

            self.FLX['ageclass']      = np.array(nc_fid1.variables['AGECLASS'])
            #self.FLX['ageclass']      = np.array(self.FLX['ageclass'] , dtype='c')

            self.FLX['times']        = nc_fid1.variables['Times']
            #self.FLX['times']        = np.array(self.FLX['times'])
            self.FLX['ReleaseName']  = np.squeeze(np.asarray(nc_fid1.variables['ReleaseName']))
            self.FLX['ReleaseName']  = np.array(self.FLX['ReleaseName'], dtype='S')
            
            self.FLX['gridarea']     = np.asarray(nc_fid1.variables['GRIDAREA'])
            self.FLX['topography']     = np.asarray(nc_fid1.variables['TOPOGRAPHY'])

            self.FLX['ReleaseZstart_end']     = np.asarray(nc_fid1.variables['ReleaseZstart_end'])
            self.FLX['ReleaseTstart_end']     = np.asarray(nc_fid1.variables['ReleaseTstart_end'])
            self.FLX['ReleaseXstart_end']     = np.asarray(nc_fid1.variables['ReleaseXstart_end'])

            listnc = glob.glob(sngDir+'flxout*.nc')

            rt2   = []
            
            for nc in listnc:

                tim = dt.datetime(int(nc[101:105]), int(nc[105:107]), int(nc[107:109]), int(nc[110:112]), int(nc[112:114]))
                deltat =  self.starttime[i] - tim

                if deltat <= dt.timedelta(hours=24):

                    nfile = nc.split('/')[-1].split('.')[0]
                    
                    nc_fid1   = Dataset(nc , 'r')  
                    nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid1, verb=False)

                    Times     = nc_fid1.variables['Times']
                    conc      = np.asarray(nc_fid1.variables['CONC'])
                    conc      = conc[0,0,6,0,:,:]

                    #rt        = conc/( (self.FLX['gridarea'])* (self.FLX['ztop'][0]))
                    rt        = conc/np.max(conc)
                    rt2.append(rt)

                ##-----Plot
                #self.PltFLX2(conc, title=nfile, ztitle='Concentration [s m$^3$ kg$^{-1}$]')

            rt2  = np.asarray(rt2)
            
            rts  = np.sum(rt2, axis=0)
            rtsa = np.sum(rt2)

            # dz = np.zeros(list(self.FLX['ztop'].shape)[0])

            # dz[0] = self.FLX['ztop'][0]
            # for x, z in enumerate(self.FLX['ztop'][1:]):
            #     dz[x+1] = z - self.FLX['ztop'][x]

            # for x, adz in enumerate(dz):
            #     rts[x] /= adz


            x = self.FLX['xlon']
            y = self.FLX['xlat']

        
            indNE = np.logical_and(y >=40.05003, x >= -105.0038)
            xNE   = x[indNE]
            yNE   = y[indNE]
            rtsNE = rts[indNE]
            rtsNE_Sum = np.sum(rtsNE)
            FracArea =  np.sum(self.FLX['gridarea'][indNE])/np.sum(self.FLX['gridarea'])

            indNW = np.logical_and(y>=40.05003, x <= -105.0038)
            xNW   = x[indNW]
            yNW   = y[indNW]
            rtsNW = rts[indNW]
            rtsNW_Sum = np.sum(rtsNW)
            FracArea =  np.sum(self.FLX['gridarea'][indNW])/np.sum(self.FLX['gridarea'])

            indSW = np.logical_and(y<=40.05003, x <= -105.0038)
            xSW   = x[indSW]
            ySW   = y[indSW]
            rtsSW = rts[indSW]
            rtsSW_Sum = np.sum(rtsSW)
            FracArea =  np.sum(self.FLX['gridarea'][indSW])/np.sum(self.FLX['gridarea'])

            indSE = np.logical_and(y<=40.05003, x >= -105.0038)
            xSE   = x[indSE]
            ySE   = y[indSE]
            rtsSE = rts[indSE]
            rtsSE_Sum = np.sum(rtsSE)
            FracArea =  np.sum(self.FLX['gridarea'][indSE])/np.sum(self.FLX['gridarea'])


            self.FLX['SC'].setdefault('NE',[]).append(rtsNE_Sum/rtsa)
            self.FLX['SC'].setdefault('NW',[]).append(rtsNW_Sum/rtsa)
            self.FLX['SC'].setdefault('SW',[]).append(rtsSW_Sum/rtsa)
            self.FLX['SC'].setdefault('SE',[]).append(rtsSE_Sum/rtsa)

            #print 'Fraction NE: ', rtsNE_Sum/rtsa
            #print 'Fraction NW: ', rtsNW_Sum/rtsa
            #print 'Fraction SW: ', rtsSW_Sum/rtsa
            #print 'Fraction SE: ', rtsSE_Sum/rtsa

            #----------------------
            #EXAMPLE OF MAP 
            #----------------------
            # fig, ax = plt.subplots()
            # ax.plot(x,y,'k.', markersize=2)
            # ax.plot(xNE,yNE,'r.', markersize=2)
            # ax.plot(xSE,ySE,'g.', markersize=2)
            # plt.show(block=False)
            # user_input = raw_input('Press any key to exit >>> ')
            # sys.exit()
            #----------------------
            #PLOT OF THE SUM
            #----------------------
            #self.PltFLX2(rts, title='24h-sum of Normalized Concentration', ztitle='Normalized Concentration')

            #----------------------
            #PLOT OF THE MEAN
            #----------------------           
            #rtm = np.mean(rt2, axis=0)
            #self.PltFLX2(rtm, title='24h-mean of Normalized Concentration', ztitle='Normalized Concentration')
            

            #if self.FLX['saveFlg']:
            #    self.FLXpdfsav.close()
            #else:
            #    user_input = raw_input('Press any key to exit >>> ')
            #    sys.exit()

            #user_input = raw_input('Press any key to exit >>> ')
            #exit()

        self.FLX['SC']['NE'] = np.asarray(self.FLX['SC']['NE'])
        self.FLX['SC']['NW'] = np.asarray(self.FLX['SC']['NW'])
        self.FLX['SC']['SW'] = np.asarray(self.FLX['SC']['SW'])
        self.FLX['SC']['SE'] = np.asarray(self.FLX['SC']['SE'])

    def PltFLX2(self, data, title='', ztitle=''):

        #----------------------

        fig, ax = plt.subplots(figsize=(7,6))
        
        m = Basemap(llcrnrlon=-110,llcrnrlat=37,urcrnrlon=-102,urcrnrlat=43,
                   projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=10)

        xbao, ybao = m(-105.0038, 40.05003)
        ax.plot(xbao,ybao,'g^', markersize=10)

        x, y = m(self.FLX['xlon'], self.FLX['xlat'])

        data[ data==0. ] = np.nan

        p = m.contourf(x, y, data)
        
        #plt.clabel(p, inline=1, fontsize=10)
        #plt.title('Simplest default with labels')

        cbar = fig.colorbar(p, orientation='vertical', fraction=0.035, pad=0.02)
        cbar.set_label(ztitle, fontsize=12)

        parallels = np.arange(35, 45, 1)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 2)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')

        ax.set_title(title)

        fig.subplots_adjust(left = 0.12, bottom=0.075, top=0.95, right = 0.87)

        if self.FLX['saveFlg']:
            self.FLXpdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


    def PltFLX(self):
        #-------------------------------------------------
        #                       Plots
        #-------------------------------------------------
        print '\nStarting FLEXPART plots'

        #---------------------------
        # Date locators for plotting
        #---------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)
        dayLc        = DayLocator()
        yearsLc      = YearLocator()
        monthLc      = MonthLocator()
        mondays      = WeekdayLocator(MONDAY)
        DateFmt      = DateFormatter('%b %d')


        minlo=np.min(self.FLX['xlon'])
        minla=np.min(self.FLX['xlat'])
        maxlo=np.max(self.FLX['xlon'])
        maxla=np.max(self.FLX['xlat'])

        conc_i    = np.asarray(self.FLX['conc_i'])
        conc      = np.sum(self.FLX['conc_i'], axis=0)


        fig, ax = plt.subplots(figsize=(7,6))

        #m = Basemap(llcrnrlon=-105.3,llcrnrlat=40.0,urcrnrlon=-105.17,urcrnrlat=40.075,
        #           projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        #m = Basemap(llcrnrlon=minlo,llcrnrlat=minla,urcrnrlon=maxlo,urcrnrlat=maxla,
        #           projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)
        
        m = Basemap(llcrnrlon=-110,llcrnrlat=37,urcrnrlon=-102,urcrnrlat=43,
                   projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

   
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=10)

        xbao, ybao = m(-105.0038, 40.05003)
        ax.plot(xbao,ybao,'g^', markersize=10)

        x, y = m(self.FLX['xlon'], self.FLX['xlat'])

        conc[ conc==0. ] = np.nan

        levels = np.arange(np.min(conc[0,0,6,0,:,:]),    np.max(conc[0,0,6,0,:,:]), 0.2)

        p = m.contourf(x, y, conc[0,0,6,0,:,:])
        
        #plt.clabel(p, inline=1, fontsize=10)
        #plt.title('Simplest default with labels')

        cbar = fig.colorbar(p, orientation='vertical', fraction=0.035, pad=0.02)
        cbar.set_label('Concentration [s m$^3$ kg$^{-1}$]', fontsize=12)

        parallels = np.arange(35, 45, 1)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 2)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')

        #m.drawmapscale(-105.255, 40.07, -105.0, 40.0, 3, barstyle='fancy')
        #ax.set_title('Flexpart, example: '+ str(n))

        # load the shapefile, use the name 'states'
        #m.readshapefile('st99_d00', name='states', drawbounds=True)

        #state_names = []
        #for shape_dict in map.states_info:
        #    state_names.append(shape_dict['NAME'])

        #print state_names


        if self.FLX['saveFlg']:
            self.FLXpdfsav.savefig(fig,dpi=200)
        else:
            plt.show(block=False)
        
        
       
##------------------------------------------------------------------------------------------------------------------------------    
class CMAQClass():
    '''
    This class deals with reading WRF-CHEM outputs.
    '''  
    def __init__(self,dataDir, file, outFname='', saveFlg= False):

        if (not dataDir) and (not file1) and (not file2):
            print 'Either Directory or Files of WRF-CHEM need to be specified'
            return False

        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'    

        self.CMAQnc   = dataDir + file

        ckDir(dataDir, exitFlg=True)

        ckFile(self.CMAQnc, exitFlg=True)

        self.CMAQ     = {}

        if saveFlg:
            self.CMAQpdfsav = PdfPages(outFname)
            self.CMAQ['saveFlg']   = True
        else: self.CMAQ['saveFlg'] = False


        #---------------
        # ReadOutputData
        #---------------
    def ReadOutputCMAQ(self, sLat, sLon):
        '''Function to Read out the WRF-CEHM net CDF files.'''

        Gases = ['NH3', 'NH3_FERT', 'CO', 'CH4', 'ETH', 'ETHA', 'ETHY', 'ETOH', 'FORM', 'FORM_PRIMARY', 'MEOH', 'HCL']

        self.CMAQ = {}

        #-------------------------------------------------
        #Reading the NetCDF File 
        #-------------------------------------------------
       
        nc_fid = Dataset(self.CMAQnc , 'r')  # Dataset is the class behavior to open the file
        nc_attrs, nc_dims, nc_vars = mf.ncdump(nc_fid, verb=True)

        TFLAG        = nc_fid.variables['TFLAG'][:]
        TFLAG        = np.asarray(TFLAG)

        TFLAG2     = [str(t) for t in TFLAG[:, 0, 0]]
        YYYY       = [int(t2[0:4]) for t2 in TFLAG2]
        d          = [int(t2[4:7]) for t2 in TFLAG2]

        TFLAG2     = [str(t) for t in TFLAG[:, 0, 1]]
        h          = np.zeros(len(d), dtype=np.int) #[int(t2[0:2]) for t2 in TFLAG2]
        for i, t in enumerate(TFLAG2):
          
            if len(t) == 5: h[i] = int(t[0:1])
            else: h[i] = int(t[0:2])

        self.CMAQ['dates'] = [ dt.datetime(y, 1, 1, h[i]) + dt.timedelta(d[i] - 1) for i, y in enumerate(YYYY)]
        self.CMAQ['dates'] = np.asarray( self.CMAQ['dates'])

        for g in Gases:
            self.CMAQ.setdefault(g, []).append(np.transpose(nc_fid.variables[g][:]))
            self.CMAQ[g] = np.squeeze(np.asarray(self.CMAQ[g]))

        #print self.CMAQ['NH3'].shape
        #print self.CMAQ['dates'].shape

        #-------------------------------------------------
        #Calculating LAT and LON (Hard coded from Gabi) 
        #-------------------------------------------------
        #lat/lon for Sojin domain
        nx        = 138
        ny        = 123

        truelat1 = 33.
        truelat2 = 45.
        std_lon  = -97.
        dx       = 4.0
        xloc     = -972             
        yloc     = -312. 

        lon_data = np.zeros((nx,ny)) 
        lat_data = np.zeros((nx,ny))


        for i in range(nx):
            for j in range(ny):
                lat_i, lon_i = self.ijll_lc_camx(xloc+i*dx+0.5*dx, yloc+j*dx+0.5*dx, truelat1, truelat2, 1, std_lon,dx)
            
                lon_data[i,j]  = lon_i     
                lat_data[i,j]  = lat_i

        #---

        #---
        print self.CMAQ['CO'].shape
        print self.CMAQ['dates'].shape

        lat_data          = np.squeeze(lat_data)
        lon_data          = np.squeeze(lon_data)

        LatLon = [zip(x,y) for x,y in zip(*(lon_data, lat_data))]
        LatLon = np.squeeze(np.asarray(LatLon))

        coord = (sLon, sLat)
        LatLon = np.reshape(LatLon, (ny*nx, 2))

        LatLonclose  = findCls2(LatLon, coord)
        print 'Lon and Lat close to FL0: {}'.format(LatLonclose)

        #indsLoc  =  np.unravel_index(LatLonclose, (nx, ny, 2))
        
        LatLon       = np.reshape(LatLon, (nx, ny, 2))
        indsLoc      = np.where(LatLonclose == LatLon)

        indsLoc = [indsLoc[0][0], indsLoc[1][0]]
      
        #print lat_data[indsLoc[0], indsLoc[1]]
        #print lon_data[indsLoc[0], indsLoc[1]]

        hours = [d.hour for d in self.CMAQ['dates']]
        hours = np.asarray(hours)
        indH  = np.where(hours >= 12)[0]

        print 'Mean Emission of CO: {}'.format(np.mean(self.CMAQ['CO'][indsLoc[0], indsLoc[1], indH]) )
        print 'Mean Emission of C2H6: {}'.format(np.mean(self.CMAQ['ETHA'][indsLoc[0], indsLoc[1], indH]) )


        #------------------
        clmap = 'jet'

        fig, (ax, ax2, ax3) = plt.subplots(3,1, figsize=(10,8), sharex=True)

        # for i in range(nx):
        #     for j in range(ny):
       
        #          ax.plot(self.CMAQ['dates'], self.CMAQ['CO'][i, j, :], 'k', alpha=0.1)
        #          ax2.plot(self.CMAQ['dates'], self.CMAQ['ETHA'][i, j, :], 'k', alpha=0.1)

        ax.scatter(self.CMAQ['dates'][indH], self.CMAQ['CO'][indsLoc[0], indsLoc[1], indH], color='red', s=6)
        ax2.scatter(self.CMAQ['dates'][indH], self.CMAQ['ETH'][indsLoc[0], indsLoc[1], indH], color='red', s=6)
        ax3.scatter(self.CMAQ['dates'][indH], self.CMAQ['NH3'][indsLoc[0], indsLoc[1], indH], color='red', s=6)

        ax.set_ylabel('CO [molecules/s]')
        #ax.set_xlabel('Date')
        ax.grid(True)

        ax2.set_ylabel('C$_2$H$_6$ [molecules/s]')
        #ax2.set_xlabel('Date')
        ax2.grid(True)

        ax3.set_ylabel('NH$_3$ [molecules/s]')
        ax3.set_xlabel('Date')
        ax3.grid(True)

        plt.show(block=False)
        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()


        #--------------C2H6
        fig, ax = plt.subplots(figsize=(8,7))
        #m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
        #           projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
                   projection='lcc', lat_0 = 40.2, lon_0 = -105.0, ax=ax)

   
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=10)

        xbao, ybao = m(-105.0038, 40.05003)
        ax.plot(xbao,ybao,'g^', markersize=10)

        x2, y2 = m(lon_data[indsLoc[0], indsLoc[1]],lat_data[indsLoc[0], indsLoc[1]])
        ax.plot(x2, y2,'k.', markersize=10)

        x, y = m(lon_data, lat_data)
        ctrmean = np.mean(self.CMAQ['ETHA'], axis=2)
        Levels = np.arange(0.1, 13, 0.5)

        p = m.contourf(x, y, ctrmean, Levels, cmap=mplcm.jet)

        cbar = fig.colorbar(p, orientation='vertical', fraction=0.04, pad=0.1)
        cbar.set_label('C$_2$H$_6$ [moles/s]', fontsize=14)

        parallels = np.arange(38, 43, 0.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 0.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')

        #--------------NH3
        fig, ax = plt.subplots(figsize=(8,7))
        #m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
        #           projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
                   projection='lcc', lat_0 = 40.2, lon_0 = -105.0, ax=ax)

   
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=10)

        xbao, ybao = m(-105.0038, 40.05003)
        ax.plot(xbao,ybao,'g^', markersize=10)

        x2, y2 = m(lon_data[indsLoc[0], indsLoc[1]],lat_data[indsLoc[0], indsLoc[1]])
        ax.plot(x2, y2,'k.', markersize=10)

        x, y = m(lon_data, lat_data)
        ctrmean = np.mean(self.CMAQ['NH3'], axis=2)
        Levels = np.arange(0.05, 0.2, 0.005)

        p = m.contourf(x, y, ctrmean, Levels, cmap=mplcm.jet)

        cbar = fig.colorbar(p, orientation='vertical', fraction=0.04, pad=0.1)
        cbar.set_label('NH$_3$ [moles/s]', fontsize=14)

        parallels = np.arange(38, 43, 0.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 0.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')

        #--------------CO
        fig, ax = plt.subplots(figsize=(8,7))
        #m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
        #           projection='lcc', lat_0 = 40.0, lon_0 = -105.0, ax=ax)

        m = Basemap(llcrnrlon=-106,llcrnrlat=39,urcrnrlon=-104,urcrnrlat=41,
                   projection='lcc', lat_0 = 40.2, lon_0 = -105.0, ax=ax)

   
        m.drawcounties()
        x, y = m(-105.245, 40.035)
        ax.plot(x,y,'r^', markersize=10)

        xbao, ybao = m(-105.0038, 40.05003)
        ax.plot(xbao,ybao,'g^', markersize=10)

        x2, y2 = m(lon_data[indsLoc[0], indsLoc[1]],lat_data[indsLoc[0], indsLoc[1]])
        ax.plot(x2, y2,'k.', markersize=10)

        x, y = m(lon_data, lat_data)
        ctrmean = np.mean(self.CMAQ['CO'], axis=2)
        Levels = np.arange(1.0, 7.0, 0.05)

        p = m.contourf(x, y, ctrmean, Levels, cmap=mplcm.jet)

        cbar = fig.colorbar(p, orientation='vertical', fraction=0.04, pad=0.1)
        cbar.set_label('CO [moles/s]', fontsize=14)

        parallels = np.arange(38, 43, 0.5)
        m.drawparallels(parallels,labels=[1,0,0,0],fontsize=13, linewidth=1.0, color='gray')

        meridians = np.arange(180.,360., 0.5)
        m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=13, linewidth=1.0, color='gray')



        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()
        

    def ijll_lc_camx(self, xloc, yloc, truelat1, truelat2, hemi, stdlon,dx):

        #-------------------------------------------------
        #Hard coded from Gabi
        #-------------------------------------------------
        iway        = 0    #; hard-wired for EPA/SMOKE domains!!!!!!!!!!!!!!!
        phic        = 40.
        xlonc       = -97.

        conv        = 57.29578
        a           = 6370.

        if phic < 0: sign = -1. 
        else: sign = 1.

        pole        = 90.
        
        if np.abs(truelat1) > 90.:
            truelat1 = 60.
            truelat2 = 30.
            truelat1 = sign*truelat1
            truelat2 = sign*truelat2

        if (truelat1 == truelat2):
            xn = np.sin(np.abs(truelat2)/conv)
        else:
            xn = np.log10(np.cos(truelat1/conv)) - np.log10(np.cos(truelat2/conv))
            xn = xn/(np.log10(np.tan((45. - np.abs(truelat1)/2.)/conv)) - np.log10(np.tan((45. - np.abs(truelat2)/2.)/conv)))

        psi1 = 90. - np.abs(truelat1)
        psi1 = psi1/conv
       
        if (phic < 0.):
            psi1 = -psi1
            pole = -pole
        
        psi0 = (pole - phic)/conv
        xc = 0.
        yc = -a/xn*np.sin(psi1)*(np.tan(psi0/2.)/np.tan(psi1/2.))**xn

        #-------------------------------------------------
        #     Calculate lat/lon of the point (xloc,yloc)
        #-------------------------------------------------

        xloc = xloc + xc
        yloc = yloc + yc
        if (yloc == 0.):
           if (xloc >= 0.): flp = 90./conv
           if (xloc < 0.):  flp = -90./conv
        else:
           if (phic < 0.): flp = np.arctan2(xloc,yloc) 
           else: flp = np.arctan2(xloc,-yloc)

        flpp = (flp/xn)*conv + xlonc

        if (flpp < -180.): flpp = flpp + 360.
        if (flpp >  180.): flpp = flpp - 360.
        xlon = flpp
        
        r = np.sqrt(xloc*xloc + yloc*yloc)
        if (phic < 0.): r = -r
       
        if (truelat1 == truelat2): cell = r/(a*np.tan(psi1))
        else: cell = (r*xn)/(a*np.sin(psi1))
       
        rxn  = 1.0/xn
        cel1 = np.tan(psi1/2.)*cell**rxn
        cel2 = np.arctan(cel1)
        psx  = 2.*cel2*conv
        ylat = pole - psx


        return ylat, xlon
           

