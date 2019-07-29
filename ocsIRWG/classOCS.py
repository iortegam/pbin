#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         pltHDFGases.py
#
# Purpose:
#         Plot time series of multiple species using HDF GEOMS
#
# Notes:  
#         Initially created for OCS (using trop height) and see trend in trop/strat OCS
#   
#
# Version History:
#       Created, Feb, 2017  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#
from scipy.io import netcdf
import os
import datetime as dt
import numpy as np
import numpy.ma as ma
import sys
import glob
import math

from scipy import interpolate
from scipy.integrate import simps

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
from collections                     import OrderedDict
import PltClass as mp
#import sondeclass as rs
import HDFClassRead as hc
import scipy.stats.mstats as stats
from mpl_toolkits.basemap import Basemap

try:
    from itertools import product
except ImportError:
    # product is new in v 2.6
    def product(*args, **kwds):
        pools = map(tuple, args) * kwds.get('repeat', 1)
        result = [[]]
        for pool in pools:
            result = [x+[y] for x in result for y in pool]
        for prod in result:
            yield tuple(prod)


    
                                    #----------------------------#
                                    #                            #
                                    #        --- Main---         #
                                    #                            #
                                    #----------------------------#

#------------------------------------------------------------------------------------------------------------------------------                     
class ReadData():

    def __init__(self, gasName, dataDir, locID, pltID, locs, errorFlg=False, pColsFlg=False, fltrFlg=False, iyear=False, imonth=False, iday=False, fyear=False, fmonth=False, fday=False):

        #-------------------------------------
        # 
        #-------------------------------------
        self.sclfct        = 1.0E12                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
        self.sclfctName    = 'ppt'                 # Name of scale factor for labeling plots
        self.TCsclfct      = 1.0e15
        self.TCsclfctName  = 'x10$^{15}$'

        #-------------------------------------
        # Create instance of output data class   
        #-------------------------------------
        statDataCl  = OrderedDict()
        
        self.Group  = zip(dataDir,locID, pltID, locs)
        self.Group.sort(key=lambda Group: Group[2])

        self.pltID  = pltID
        self.pltID.sort()
        
        self.locs   = [l for dd, id, pl, l in self.Group]


        for dd, id, pl, l in self.Group:
            #-------------------------------------
            # Some HDF files are in specific folder: change here accordingly
            #-------------------------------------
            if pl == 'Wollongong':      dd = dd + 'ocs_hippov2/'
            elif pl == 'Jungfraujoch' : dd = dd + 'hdf_OCS/'
            elif pl == 'Toronto' :      dd = dd + 'OCS/'
            elif pl == 'Rikubetsu':     dd = dd + 'HDF_Fil4/'
            elif pl == 'Tsukuba' :      dd = dd + 'HDFfiles/'
            elif pl == 'Zugspitze':     dd = dd + 'OCS_Zugspitze/'
            elif pl == 'Kiruna':        dd = dd + 'OCS_2019/'
            elif pl == 'Izana':         dd = dd + 'OCS_2019/'
            elif pl == 'St Petersburg': dd = dd + 'HDF_OCS_SPb_O3_atm16/'
            elif pl == 'Paris':         dd = dd + '2019_Paris/'
            else: dd = dd

            statDataCl[pl] = hc.ReadHDFData(dd, id, gasName)

        #-------------------------------------
        # Variables from HDF files 
        #-------------------------------------
        self.datesJD2K    = OrderedDict()
        self.rPrf         = OrderedDict()   #retrieved Prf in mixing ratio
        self.aPrf         = OrderedDict()   #apriori Prf in mixing ratio
        self.rPrfMol      = OrderedDict()   #retrieved Prf partial Column (molec/cm2)
        self.aPrfMol      = OrderedDict()   #apriori Prf partial Column (molec/cm2)
        self.totClmn      = OrderedDict()   #retrieved total column (molec/cm2)
        self.atotClmn     = OrderedDict()   #apriori total column (molec/cm2)
        self.avkVMR       = OrderedDict()   #Averaging kernel (VMR)
        self.avkTC        = OrderedDict()   #Averaging kernel total column
        self.alt          = OrderedDict()   #Altitude 
        self.sza          = OrderedDict()   #Solar Zenith Angle
        self.TempPrf      = OrderedDict()   #Temperature Profile
        self.PresPrf      = OrderedDict()   #Pressure Profile

        #-------------------------------------
        # Variables calculated 
        #-------------------------------------
        #alt_orig     = OrderedDict()
        self.dates        = OrderedDict()
        self.avkSCF       = OrderedDict()   #Averaging kernel (scale factor)
        self.dofs         = OrderedDict()   #degrees of freedom
        self.AirMPrf      = OrderedDict()   #Airmass
        self.rPrfMol      = OrderedDict()   #retrieved Prf in molec/cm2
        self.aPrfMol      = OrderedDict()   #apriori Prf in molec/cm2

        self.totWvmr      = OrderedDict()    #Weightet VMR A priori
        self.atotWvmr     = OrderedDict()

        self.alttpp       = OrderedDict()
        self.alttpp2      = OrderedDict()

        self.altbl1       = OrderedDict()
        self.altbl2       = OrderedDict()

        self.altft1       = OrderedDict()
        self.altft2       = OrderedDict()

        self.altst1       = OrderedDict()
        self.altst2       = OrderedDict()

        self.Lat          = []
        self.Lon          = []

        if errorFlg:
            self.tot_rnd       = OrderedDict()
            self.tot_sys       = OrderedDict()
            self.tot_std       = OrderedDict()
            self.vmr_rnd_err   = OrderedDict()
            self.vmr_sys_err   = OrderedDict()
            self.vmr_tot_err   = OrderedDict()

        if pColsFlg:
            self.dtp           = OrderedDict()
            self.datesdtp      = OrderedDict()
            
            self.PcolStrat     = OrderedDict()   #partial columns
            self.PcolTrop1      = OrderedDict()
            self.PcolTrop2     = OrderedDict()

            self.PcolStratapr  = OrderedDict()   #partial columns A priori
            self.PcolTropapr1   = OrderedDict()
            self.PcolTropapr2  = OrderedDict()

            self.WvmrStrat     = OrderedDict()   #Weighted VMR
            self.WvmrTrop1      = OrderedDict()
            self.WvmrTrop2     = OrderedDict()

            self.WvmrStratapr  = OrderedDict()    #Weighted VMR A priori
            self.WvmrTropapr1   = OrderedDict()
            self.WvmrTropapr2  = OrderedDict()

            self.rPcol         = OrderedDict() 
            self.aPcol         = OrderedDict()

            self.rPvmr         = OrderedDict()
            self.aPvmr         = OrderedDict()


        for ii, idhdf in enumerate(self.pltID):

            self.datesJD2K[idhdf]    = statDataCl[idhdf].HDF[statDataCl[idhdf].getDatetimeName()]
            self.dates[idhdf]        = hc.jdf_2_datetime(self.datesJD2K[idhdf])
            self.alt[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeName()]
            self.sza[idhdf]          = statDataCl[idhdf].HDF[statDataCl[idhdf].getAngleSolarZenithAstronomicalName()]
            
            conv                     = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()+'VAR_SI_CONVERSION']            
            self.rPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarName()]*float(conv[0][1])*self.sclfct
            self.aPrf[idhdf]         = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAprioriName()]*float(conv[0][1])*self.sclfct

            conv                      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()+'VAR_SI_CONVERSION']
            self.rPrfMol[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarName()]*float(conv[0][1])*(6.02e23/100./100.)
            self.aPrfMol[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnPartialAbsorptionSolarAprioriName()]*float(conv[0][1])*(6.02e23/100./100.)

            conv                     = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()+'VAR_SI_CONVERSION']
            self.totClmn[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarName()]*float(conv[0][1]) * (6.02e23) /100./100. / self.TCsclfct
            self.atotClmn[idhdf]     = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarAprioriName()]*float(conv[0][1]) * (6.02e23) /100./100. / self.TCsclfct
            
            self.PresPrf[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].getPressureIndependentName()]
            self.TempPrf[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].getTemperatureIndependentName()]

            AltBo                    = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeBoundariesName()]
                 
            nobs                     = self.rPrf[idhdf].shape[0]
            n_layer                  = self.rPrf[idhdf].shape[1]

            if statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName() in statDataCl[idhdf].HDF.keys():
                self.avkVMR[idhdf]       = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName()]
                self.avkTC[idhdf]        = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarAvkName()]
            else:
                self.avkVMR[idhdf]  = np.empty([nobs,n_layer,n_layer])
                self.avkTC[idhdf]   = np.empty([nobs,n_layer,n_layer])
                self.avkVMR[idhdf].fill('nan')
                self.avkTC[idhdf].fill('nan')

            #----------------------------------------
            #CALCULATED AIR MASS
            #----------------------------------------
            self.AirMPrf[idhdf]     =  np.divide(self.rPrfMol[idhdf], self.rPrf[idhdf])*self.sclfct

            #----------------------------------------
            #EXTRACT SINGLE ALTITUDE VECTOR
            #----------------------------------------
            if (idhdf == 'Kiruna') or (idhdf == 'Zugspitze') or (idhdf == 'Izana') or (idhdf == 'Paris'):
                self.alt[idhdf]          = self.alt[idhdf][0, :]
            else:
                self.alt[idhdf]          = self.alt[idhdf][0:n_layer]

            #----------------------------------------
            #READ LAT/LON/HEIGHT OF INSTRUMENT
            #----------------------------------------
            Lat_i           = statDataCl[idhdf].HDF[statDataCl[idhdf].getLatitudeInstrumentName()]
            Lon_i           = statDataCl[idhdf].HDF[statDataCl[idhdf].getLongitudeInstrumentName()]
            alt_instru      = statDataCl[idhdf].HDF[statDataCl[idhdf].getAltitudeInstrumentName()]

            self.Lat.append(float(Lat_i[0]))
            self.Lon.append(float(Lon_i[0]))

            print '\n'
            print idhdf
            print 'Latitude          = {0:.2f}'.format(Lat_i[0])
            print 'Longitude         = {0:.2f}'.format(Lon_i[0])
            print 'Altitude of Instr = {0:.2f}'.format(alt_instru[0])        

            #----------------------------------------
            #CALCULATE SCALING FACTOR AK
            #----------------------------------------
            if statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarAvkName() in statDataCl[idhdf].HDF.keys():
                self.avkSCF[idhdf]  = np.zeros((nobs,n_layer,n_layer))

                for obs in range(0,nobs):
                    Iapriori        = np.zeros((n_layer,n_layer))
                    IaprioriInv     = np.zeros((n_layer,n_layer))
                    np.fill_diagonal(Iapriori, self.aPrf[idhdf][obs])
                    np.fill_diagonal(IaprioriInv, 1.0 / (self.aPrf[idhdf][obs]))
                    self.avkSCF[idhdf][obs,:,:] = np.dot(np.dot(IaprioriInv,np.squeeze(self.avkVMR[idhdf][obs,:,:])),Iapriori)

                self.dofs[idhdf]         = np.asarray([np.trace(aki) for aki in self.avkSCF[idhdf]])
            else:
                self.avkSCF[idhdf]  = np.zeros((nobs,n_layer,n_layer))
                self.avkSCF[idhdf].fill('nan')

            #----------------------------------------
            #OBTAIN ERROR VARIABLES
            #---------------------------------------- 
            if errorFlg:                               
                self.tot_rnd[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintyRandomName()]
                self.tot_sys[idhdf]      = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getColumnAbsorptionSolarUncertaintySystematicName()]
                self.tot_std[idhdf]      = np.sqrt(self.tot_rnd[idhdf]**2 + self.tot_sys[idhdf]**2)

                npnts               = np.shape(self.tot_std[idhdf])[0]
                nlvls               = np.shape(self.alt[idhdf])[0]
                
                self.vmr_rnd_err[idhdf]  = np.zeros((npnts,nlvls))
                self.vmr_sys_err[idhdf]  = np.zeros((npnts,nlvls))

                for i in range(npnts):
                    conv    = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()+'VAR_SI_CONVERSION']  
                    cov_rnd = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintyRandomName()]
                    cov_sys = statDataCl[idhdf].HDF[statDataCl[idhdf].PrimaryGas.upper()+'.'+statDataCl[idhdf].getMixingRatioAbsorptionSolarUncertaintySystematicName()]

                    #vmr_rnd_err[idhdf][i,:] = np.diag(cov_rnd[i][:,:])*float(conv[0][1])*sclfct**2
                    #vmr_sys_err[idhdf][i,:] = np.diag(cov_sys[i][:,:])*float(conv[0][1])*sclfct**2

                    self.vmr_rnd_err[idhdf][i,:] = np.diag(cov_rnd[i][:,:])*float(conv[0][1])
                    self.vmr_sys_err[idhdf][i,:] = np.diag(cov_sys[i][:,:])*float(conv[0][1])

                self.vmr_tot_err[idhdf]  = np.sqrt(self.vmr_rnd_err[idhdf]+ self.vmr_sys_err[idhdf]) *self.sclfct
                self.vmr_rnd_err[idhdf]  = np.sqrt(self.vmr_rnd_err[idhdf])*self.sclfct
                self.vmr_sys_err[idhdf]  = np.sqrt(self.vmr_sys_err[idhdf])*self.sclfct

            #----------------------------------------
            # FILTER DATA
            #----------------------------------------
            if fltrFlg: statDataCl[idhdf].fltrData(statDataCl[idhdf].PrimaryGas,iyear=iyear, imonth=imonth, iday=iday, fyear=fyear, fmonth=fmonth, fday=fday, dateFlg=True, tcFlg=True,pcFlg=True)
            else:    statDataCl[idhdf].inds = np.array([]) 
            
            try:
                self.dates[idhdf]    = np.delete(self.dates[idhdf], statDataCl[idhdf].inds)
                self.sza[idhdf]      = np.delete(self.sza[idhdf], statDataCl[idhdf].inds)
                self.totClmn[idhdf]  = np.delete(self.totClmn[idhdf], statDataCl[idhdf].inds)
                self.atotClmn[idhdf] = np.delete(self.atotClmn[idhdf], statDataCl[idhdf].inds)
                self.rPrf[idhdf]     = np.delete(self.rPrf[idhdf], statDataCl[idhdf].inds, axis=0)
                self.rPrfMol[idhdf]  = np.delete(self.rPrfMol[idhdf], statDataCl[idhdf].inds, axis=0)
                self.aPrf[idhdf]     = np.delete(self.aPrf[idhdf], statDataCl[idhdf].inds, axis=0)
                self.aPrfMol[idhdf]  = np.delete(self.aPrfMol[idhdf], statDataCl[idhdf].inds, axis=0)
                self.avkVMR[idhdf]   = np.delete(self.avkVMR[idhdf], statDataCl[idhdf].inds, axis=0)
                self.avkSCF[idhdf]   = np.delete(self.avkSCF[idhdf], statDataCl[idhdf].inds, axis=0)
                self.avkTC[idhdf]    = np.delete(self.avkTC[idhdf], statDataCl[idhdf].inds, axis=0)
                self.AirMPrf[idhdf]  = np.delete(self.AirMPrf[idhdf], statDataCl[idhdf].inds, axis=0)

            except Exception as errmsg:
                print '\nError: ', errmsg


            #if pColsFlg: 
            #    dtp[idhdf]      = np.delete(dtp[idhdf], statDataCl[idhdf].inds)
            #    datesdtp[idhdf] = np.delete(datesdtp[idhdf], statDataCl[idhdf].inds)

            if errorFlg:
                try:
                    self.vmr_rnd_err[idhdf]  = np.delete(self.vmr_rnd_err[idhdf],statDataCl[idhdf].inds,axis=0)
                    self.vmr_sys_err[idhdf]  = np.delete(self.vmr_sys_err[idhdf],statDataCl[idhdf].inds,axis=0)  
                    self.vmr_tot_err[idhdf]  = np.delete(self.vmr_tot_err[idhdf],statDataCl[idhdf].inds,axis=0)

                    self.tot_rnd[idhdf]      = np.delete(self.tot_rnd[idhdf],statDataCl[idhdf].inds)
                    self.tot_sys[idhdf]      = np.delete(self.tot_sys[idhdf],statDataCl[idhdf].inds)
                    self.tot_std[idhdf]      = np.delete(self.tot_std[idhdf],statDataCl[idhdf].inds)
                except Exception as errmsg:
                    print '\nError: ', errmsg


            if pColsFlg:

                #---------------------------------------------------
                #STATISTICS OF TROPOPAUSE HEIGHT BASED ON DAILY AVERAGES
                #---------------------------------------------------
                # AvgTpp       = mf.dailyAvg(dtp[idhdf], dates[idhdf], dateAxis=1, meanAxis=0)
                # AvgTpp       = AvgTpp['dailyAvg']

                # maxTpp       = np.max(AvgTpp)
                # minTpp       = np.min(AvgTpp)
                # meanTpp      = np.mean(AvgTpp)
                # stdTpp       = np.std(AvgTpp)

                # print '\nMean TPH: {0:.2f} +/- {1:.2f}'.format(meanTpp, stdTpp)

                #----------------------------------------------------
                #
                #----------------------------------------------------
                if float(Lat_i[0]) >=70.: 
                    meanTpp = 8.8
                    stdTpp  = 1.2
                
                elif (float(Lat_i[0]) >= 60.0) & (float(Lat_i[0]) < 70.0):
                    meanTpp = 9.8
                    stdTpp  = 1.3
                
                elif (float(Lat_i[0]) >= 50.0) & (float(Lat_i[0]) < 60.0):
                    meanTpp = 10.9
                    stdTpp  = 1.2

                elif (float(Lat_i[0]) >= 40.0) & (float(Lat_i[0]) < 50.0):
                    meanTpp = 11.6
                    stdTpp  = 1.6

                elif (float(Lat_i[0]) >= 30.0) & (float(Lat_i[0]) < 40.0):
                    meanTpp = 12.9 #12.58
                    stdTpp  = 2.4  #2.72

                elif (float(Lat_i[0]) >= 20.0) & (float(Lat_i[0]) < 30.0):
                    meanTpp = 15.0
                    stdTpp  = 1.3

                elif (float(Lat_i[0]) >= -25.0) & (float(Lat_i[0]) < 20.0):
                    meanTpp = 16.5
                    stdTpp  = 0.4

                elif (float(Lat_i[0]) >= -40.0) & (float(Lat_i[0]) < -25.0):
                    meanTpp = 12.3
                    stdTpp  = 2.2

                elif (float(Lat_i[0]) >= -50.0) & (float(Lat_i[0]) < -40.0):
                    meanTpp = 11.1
                    stdTpp  = 1.3

                elif float(Lat_i[0]) < -50:
                    meanTpp = 8.8
                    stdTpp  = 1.7


                partialCols  = [ [0.0, 4.0], [4.0, (meanTpp - stdTpp*2.)], [(meanTpp+stdTpp*2.), 40.] ]

                self.alttpp[idhdf]       = np.zeros((len(self.rPrfMol[idhdf][:,0])))
                self.alttpp2[idhdf]      = np.zeros((len(self.rPrfMol[idhdf][:,0])))

                ind1         = mf.nearestind(partialCols[1][1], self.alt[idhdf])
                ind2         = mf.nearestind(partialCols[2][0], self.alt[idhdf])

                inds = np.where( (self.alt[idhdf] >= partialCols[1][1]) & (self.alt[idhdf] <=partialCols[2][0])  )[0]

                #alttpp[idhdf][:]      = alt[idhdf][ind1]
                #alttpp2[idhdf][:]     = alt[idhdf][ind2]

                self.alttpp[idhdf][:]      = self.alt[idhdf][inds[-1]]
                self.alttpp2[idhdf][:]     = self.alt[idhdf][inds[0]]

            

                for ii, pc in enumerate(partialCols):

                    #ind1         = mf.nearestind(pc[0], alt[idhdf])
                    #ind2         = mf.nearestind(pc[1], alt[idhdf])

                    inds = np.where( (self.alt[idhdf] >= pc[0]) & (self.alt[idhdf] <= pc[1])  )[0]


                    #---------------------------------------------------
                    #THESE SITES REPORT INCREASING ALTITUDE
                    #---------------------------------------------------
                    if (idhdf == 'Kiruna') or (idhdf == 'Izana') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):       
                        
                        #rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,ind1:ind2], axis=1)
                        #aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,ind1:ind2], axis=1)

                        self.rPcol[idhdf+str(pc)]  = np.sum(self.rPrfMol[idhdf][:,inds], axis=1)
                        self.aPcol[idhdf+str(pc)]  = np.sum(self.aPrfMol[idhdf][:,inds], axis=1)

                        try:
                            #rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,ind1:ind2], weights=AirMPrf[idhdf][:,ind1:ind2],axis=1)
                            #aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,ind1:ind2], weights=AirMPrf[idhdf][:,ind1:ind2],axis=1)

                            self.rPvmr[idhdf+str(pc)]  = np.average(self.rPrf[idhdf][:,inds], weights=self.AirMPrf[idhdf][:,inds],axis=1)
                            self.aPvmr[idhdf+str(pc)]  = np.average(self.aPrf[idhdf][:,inds], weights=self.AirMPrf[idhdf][:,inds],axis=1)
                        
                        except Exception as errmsg:
                            self.rPvmr[idhdf+str(pc)]    = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                            self.rPvmr[idhdf+str(pc)][:] = float('nan')

                            self.aPvmr[idhdf+str(pc)]    = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                            self.aPvmr[idhdf+str(pc)][:] = float('nan')

                    else: 
                        #rPcol[idhdf+str(pc)]  = np.sum(rPrfMol[idhdf][:,ind2:ind1], axis=1)
                        #aPcol[idhdf+str(pc)]  = np.sum(aPrfMol[idhdf][:,ind2:ind1], axis=1)

                        self.rPcol[idhdf+str(pc)]  = np.sum(self.rPrfMol[idhdf][:,inds], axis=1)
                        self.aPcol[idhdf+str(pc)]  = np.sum(self.aPrfMol[idhdf][:,inds], axis=1)

                        try:
                            #rPvmr[idhdf+str(pc)]  = np.average(rPrf[idhdf][:,ind2:ind1], weights=AirMPrf[idhdf][:,ind2:ind1],axis=1)
                            #aPvmr[idhdf+str(pc)]  = np.average(aPrf[idhdf][:,ind2:ind1], weights=AirMPrf[idhdf][:,ind2:ind1],axis=1)

                            self.rPvmr[idhdf+str(pc)]  = np.average(self.rPrf[idhdf][:,inds], weights=self.AirMPrf[idhdf][:,inds],axis=1)
                            self.aPvmr[idhdf+str(pc)]  = np.average(self.aPrf[idhdf][:,inds], weights=self.AirMPrf[idhdf][:,inds],axis=1)

                        except Exception as errmsg:
                            self.rPvmr[idhdf+str(pc)]    = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                            self.rPvmr[idhdf+str(pc)][:] = float('nan')

                            self.aPvmr[idhdf+str(pc)]    = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                            self.aPvmr[idhdf+str(pc)][:] = float('nan')

                    if ii == 0:
                        self.PcolTrop1[idhdf]     = np.asarray(self.rPcol[idhdf+str(pc)])/self.TCsclfct
                        self.PcolTropapr1[idhdf]  = np.asarray(self.aPcol[idhdf+str(pc)])/self.TCsclfct

                        self.WvmrTrop1[idhdf]     = np.asarray(self.rPvmr[idhdf+str(pc)])
                        self.WvmrTropapr1[idhdf]  = np.asarray(self.aPvmr[idhdf+str(pc)])

                        self.altbl1[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altbl1[idhdf][:]    = np.asarray(alt[idhdf][ind1])
                        self.altbl1[idhdf][:]    = np.asarray(self.alt[idhdf][inds[-1]])

                        self.altbl2[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altbl2[idhdf][:]    = np.asarray(alt[idhdf][ind2])
                        self.altbl2[idhdf][:]    = np.asarray(self.alt[idhdf][inds[0]])


                    elif ii == 1:
                        self.PcolTrop2[idhdf]     = np.asarray(self.rPcol[idhdf+str(pc)])/self.TCsclfct
                        self.PcolTropapr2[idhdf]  = np.asarray(self.aPcol[idhdf+str(pc)])/self.TCsclfct

                        self.WvmrTrop2[idhdf]     = np.asarray(self.rPvmr[idhdf+str(pc)])
                        self.WvmrTropapr2[idhdf]  = np.asarray(self.aPvmr[idhdf+str(pc)])

                        self.altft1[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altft1[idhdf][:]    = np.asarray(alt[idhdf][ind1])
                        self.altft1[idhdf][:]    = np.asarray(self.alt[idhdf][inds[-1]])

                        self.altft2[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altft2[idhdf][:]    = np.asarray(alt[idhdf][ind2])
                        self.altft2[idhdf][:]    = np.asarray(self.alt[idhdf][inds[0]])

                    elif ii == 2:
                        self.PcolStrat[idhdf]    = np.asarray(self.rPcol[idhdf+str(pc)])/self.TCsclfct
                        self.PcolStratapr[idhdf] = np.asarray(self.aPcol[idhdf+str(pc)])/self.TCsclfct

                        self.WvmrStrat[idhdf]    = np.asarray(self.rPvmr[idhdf+str(pc)])
                        self.WvmrStratapr[idhdf] = np.asarray(self.aPvmr[idhdf+str(pc)])

                        self.altst1[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altst1[idhdf][:]    = np.asarray(alt[idhdf][ind1])
                        self.altst1[idhdf][:]    = np.asarray(self.alt[idhdf][inds[-1]])

                        self.altst2[idhdf]       = np.zeros(len(self.rPrfMol[idhdf][:,0]))
                        #altst2[idhdf][:]    = np.asarray(alt[idhdf][ind2])
                        self.altst2[idhdf][:]    = np.asarray(self.alt[idhdf][inds[0]])


            self.totWvmr[idhdf]  = np.average(self.rPrf[idhdf], axis=1, weights=self.AirMPrf[idhdf])
            #self.atotWvmr[idhdf] = np.average(self.aPrf[idhdf], axis=1, weights=self.AirMPrf[idhdf])

        #----------------------------
        #CONCATENATE st denis and maido
        #---------------------------- 
        self.dates2         = {}
        
        self.PcolTrop12     = {}
        self.PcolTrop22     = {}
        self.PcolStrat2     = {}
        
        self.PcolTropapr12  = {}
        self.PcolTropapr22  = {}
        self.PcolStratapr2  = {}

        self.WvmrTrop12     = {}
        self.WvmrTrop22     = {}
        self.WvmrStrat2     = {}
         
        self.WvmrTropapr12  = {}
        self.WvmrTropapr22  = {}
        self.WvmrStratapr2  = {}

        self.totClmn2       = {}

        self.dtp2           = {}
        self.alttpp12       = {}
        self.alttpp22       = {}

        self.pltID2         = []
        self.Lat2           = []
        self.Lon2           = []
        self.locID2         = []

        for i, idhdf in enumerate(self.pltID):

            if idhdf == 'Maido': continue
            
            if idhdf  == 'St Denis':
                
                self.dates2[idhdf]          = np.concatenate( (self.dates[idhdf], self.dates['Maido']))

                self.totClmn2[idhdf]       = np.concatenate( (self.totClmn[idhdf], self.totClmn['Maido']))

                self.PcolTrop12[idhdf]       = np.concatenate( (self.PcolTrop1[idhdf], self.PcolTrop1['Maido']))
                self.PcolTrop22[idhdf]       = np.concatenate( (self.PcolTrop2[idhdf], self.PcolTrop2['Maido']))
                self.PcolStrat2[idhdf]      = np.concatenate( (self.PcolStrat[idhdf], self.PcolStrat['Maido']))
                
                self.PcolTropapr12[idhdf]    = np.concatenate( (self.PcolTropapr1[idhdf], self.PcolTropapr1['Maido']))
                self.PcolTropapr22[idhdf]    = np.concatenate( (self.PcolTropapr2[idhdf], self.PcolTropapr2['Maido']))
                self.PcolStratapr2[idhdf]   = np.concatenate( (self.PcolStratapr[idhdf], self.PcolStratapr['Maido']))

                self.WvmrTrop12[idhdf]       = np.concatenate( (self.WvmrTrop1[idhdf], self.WvmrTrop1['Maido']))
                self.WvmrTrop22[idhdf]       = np.concatenate( (self.WvmrTrop2[idhdf], self.WvmrTrop2['Maido']))
                self.WvmrStrat2[idhdf]      = np.concatenate( (self.WvmrStrat[idhdf], self.WvmrStrat['Maido']))
                
                self.WvmrTropapr12[idhdf]    = np.concatenate( (self.WvmrTropapr1[idhdf], self.WvmrTropapr1['Maido']))
                self.WvmrTropapr22[idhdf]    = np.concatenate( (self.WvmrTropapr2[idhdf], self.WvmrTropapr2['Maido']))
                self.WvmrStratapr2[idhdf]   = np.concatenate( (self.WvmrStratapr[idhdf], self.WvmrStratapr['Maido']))

                #dtp2[idhdf]              = np.concatenate( (dtp[idhdf], dtp['Maido']))
                #alttpp12[idhdf]           = np.concatenate( (alttpp[idhdf], alttpp['Maido']))
                #alttpp22[idhdf]           = np.concatenate( (alttpp2[idhdf], alttpp2['Maido']))

                self.pltID2.append(self.pltID[i])
                self.Lat2.append(self.Lat[i])
                self.Lon2.append(self.Lon[i])
                self.locID2.append('STD-MAI')

            else:

                self.dates2[idhdf]         = self.dates[idhdf]
                self.totClmn2[idhdf]       = self.totClmn[idhdf]
                self.PcolTrop12[idhdf]     = self.PcolTrop1[idhdf] 
                self.PcolTrop22[idhdf]     = self.PcolTrop2[idhdf] 
                self.PcolStrat2[idhdf]     = self.PcolStrat[idhdf]
                self.PcolTropapr12[idhdf]   = self.PcolTropapr1[idhdf]
                self.PcolTropapr22[idhdf]   = self.PcolTropapr2[idhdf]
                self.PcolStratapr2[idhdf]  = self.PcolStratapr[idhdf]

                self.WvmrTrop12[idhdf]      = self.WvmrTrop1[idhdf] 
                self.WvmrTrop22[idhdf]      = self.WvmrTrop2[idhdf] 
                self.WvmrStrat2[idhdf]     = self.WvmrStrat[idhdf]
                self.WvmrTropapr12[idhdf]   = self.WvmrTropapr1[idhdf]
                self.WvmrTropapr22[idhdf]   = self.WvmrTropapr2[idhdf]
                self.WvmrStratapr2[idhdf]  = self.WvmrStratapr[idhdf]

                #dtp2[idhdf]            = dtp[idhdf]
                #alttpp12[idhdf]        = alttpp[idhdf]
                #alttpp22[idhdf]        = alttpp2[idhdf]

                self.pltID2.append(self.pltID[i])
                self.Lat2.append(self.Lat[i])
                self.Lon2.append(self.Lon[i])
                self.locID2.append(self.locs[i])

#------------------------------------------------------------------------------------------------------------------------------ 
#------------------------------------------------------------------------------------------------------------------------------                     
class pltOCS(ReadData):

    def __init__(self, gasName, dataDir, locID, pltID, locs, iyear=False, imonth=False, iday=False, fyear=False, fmonth=False, fday=False, saveFlg=False, outFname='', errorFlg=False, pColsFlg=False, fltrFlg=False):

        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg: self.pdfsav = PdfPages(outFname)
        else:       self.pdfsav = False

        self.clr = mf.clrplt()
        
        ReadData.__init__(self,gasName, dataDir, locID, pltID, locs, errorFlg=errorFlg, pColsFlg=pColsFlg, fltrFlg=fltrFlg, iyear=iyear, imonth=imonth, iday=iday, fyear=fyear, fmonth=fmonth, fday=fday,)

        #---------------------------------------------------
        # Determine if multiple years
        #---------------------------------------------------
        self.years  = []
        self.Nyears = []

        for i, idhdf in enumerate(self.pltID):
            self.years.append([ singDate.year for singDate in self.dates[idhdf]] )
        
        self.years = np.asarray(self.years)

        self.Nyears = [len(list(set(y))) for y in self.years]
        indY = np.max(self.Nyears)
        
        if indY > 1: self.yrsFlg = True         
        else:        self.yrsFlg = False

        self.clmap        = 'jet'
        self.cm           = plt.get_cmap(self.clmap)    
        self.yearsLc      = YearLocator()
        self.monthsAll    = MonthLocator()
        #months       = MonthLocator(bymonth=1,bymonthday=1)
        self.months       = MonthLocator()

        if self.yrsFlg: self.DateFmt      = DateFormatter('%Y')
        else:      self.DateFmt      = DateFormatter('%m\n%Y')


        #---------------------------------------------------
        # Order data based on +Lat to -Lat
        #---------------------------------------------------
        self.pltID  = [y for (x, y) in sorted(zip(self.Lat,self.pltID), reverse=True)]
        self.locsID = [y for (x, y) in sorted(zip(self.Lat,self.locs), reverse=True)]
        self.Lon    = [y for (x, y) in sorted(zip(self.Lat,self.Lon), reverse=True)]
        
        self.Lat    = sorted(self.Lat, reverse=True)
        self.pltID  = np.asarray(self.pltID)
        self.locsID = np.asarray(self.locsID)

        self.pltID2 = [y for (x,y) in sorted(zip(self.Lat2,self.pltID2), reverse=True)]
        self.Lon2   = [y for (x,y) in sorted(zip(self.Lat2,self.Lon2), reverse=True)]
        self.locID2 = [y for (x,y) in sorted(zip(self.Lat2,self.locID2), reverse=True)]
        self.Lat2   = sorted(self.Lat2, reverse=True)
        self.pltID2 = np.asarray(self.pltID2)
        self.locID2 = np.asarray(self.locID2)

        #---------------------------------------------------
        # Defining variable for plots
        #---------------------------------------------------
        self.npanels   = len(self.pltID)
        self.npanels2  = int(math.ceil(self.npanels/4.0))
        self.npanels3  = int(math.ceil(self.npanels/3.0))

        
        self.xmin      = dt.date(iyear, imonth, iday)
        self.xmax      = dt.date(fyear, fmonth, fday)

    def closeFig(self):
        self.pdfsav.close()


    def pltMap(self):

        #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        #                                                           PLOTS
        #----------------------------------------------------------------------------------------------------------------------------------------------------------------
        print '\n Plotting map.......\n'
            
        #---------------------------------------------------
        # Map with Sites
        #---------------------------------------------------
        fig = plt.figure(figsize=(11,7))

        #eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0, lat_0=0, lon_0=-130)
        eq_map = Basemap(projection='robin', resolution = 'l', area_thresh = 1000.0, lat_0=0, lon_0=0)
        eq_map.drawcoastlines(linewidth =0.25, color = 'gray', )
        #eq_map.drawcountries()
        eq_map.fillcontinents(color = 'gray', alpha=0.5)
        eq_map.drawmapboundary()
        eq_map.drawmeridians(np.arange(0, 360, 30), color='gray', linewidth=0.5, alpha = 0.5)
        eq_map.drawparallels(np.arange(-90, 90, 30), color='gray', linewidth=0.5, labels=[1,0,0,0], alpha = 0.5)
         
        for lo, la, idhdf in zip(self.Lon, self.Lat, self.locsID):

            if la >= 50.:     
                self.clr = 'lightcyan'
                
            elif (la >= 20.) & (la < 50.):
                self.clr = 'lightgreen'
               
            elif (la >= -20.)  & (la < 20.):
                self.clr = 'mistyrose'
                
            elif (la >= -50.)  & (la < -20.):
                self.clr = 'cornsilk'
                
            elif (la < -50.):
                self.clr = 'lightgrey'
            
            else:
                self.clr = 'lightgrey'

            x,y = eq_map(lo, la)
            eq_map.plot(x, y, 'yo', markersize=0)
            eq_map.scatter(x, y, s=80, facecolor='yellow', edgecolor='k', zorder=7)

            if idhdf.upper() == 'BRE': plt.text(x-1000000,y+300000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'EUR': plt.text(x-1800000,y-100000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'KIR': plt.text(x+50000,y+200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'MAI': plt.text(x-30000,y+300000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'PAR': plt.text(x-1800000,y, idhdf.upper(), fontsize=12, color='r', weight='bold')
            #elif idhdf.upper() == 'AHS': plt.text(x+250000,y-200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'JFJ': plt.text(x-1000000,y-800000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'ZGP': plt.text(x+100000,y+200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'STP': plt.text(x+30000,y+200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'NYA': plt.text(x+80000,y+100000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'SPU': plt.text(x+30000,y+300000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'RKB': plt.text(x+100000,y+200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'LDR': plt.text(x-1800000,y-100000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'WLG': plt.text(x+100000,y+200000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            elif idhdf.upper() == 'AHS': plt.text(x-1800000,y-100000, idhdf.upper(), fontsize=12, color='r', weight='bold')
            else: plt.text(x+10000,y-800000, idhdf.upper(), fontsize=12, color='r', weight='bold')
                  
        #plt.title("IRWG-NDACC Sites participating in the OCS global study")

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            #fig.savefig('/data1/projects/ocs/figures/fig/'+'Map_OCS.pdf', bbox_inches='tight')
            fig.savefig('/data1/projects/ocs/figures/fig/'+'Map.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)
        
        #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

    def pltPrf(self):

        print '\n Plotting Profiles.......\n'

        #---------------------------------------------------
        # Mean vertical profiles
        #---------------------------------------------------
        fig = plt.figure(figsize=(12,12))

        outer_grid = gridspec.GridSpec(self.npanels2, 4, wspace=0.11, hspace=0.085)

        for i, idhdf in enumerate(self.pltID):

            ax = plt.Subplot(fig, outer_grid[i])

            Lat = float(self.Lat[i])

            if Lat >= 50.:     
                ax.set_facecolor('lightcyan')
                
            elif (Lat >= 20.) & (Lat < 50.):
                #ax.set_axis_bgcolor('lightgreen')
                ax.set_facecolor('lightgreen')
               
            elif (Lat >= -20.)  & (Lat < 20.):
                ax.set_facecolor('mistyrose')
                
            elif (Lat >= -50.)  & (Lat < -20.):
                ax.set_facecolor('cornsilk')
                
            elif (Lat < -50.):
                ax.set_facecolor('lightgrey')
            
            else:
                ax.set_facecolor('lightgrey')


            prfMean    = np.nanmean(self.rPrf[idhdf],axis=0)
            prfSTD     = np.nanstd(self.rPrf[idhdf],axis=0)
            avkVMRAv   = np.nanmean(self.avkVMR[idhdf], axis=0)
            avkSCFAv   = np.nanmean(self.avkSCF[idhdf], axis=0)
            ax.plot(prfMean,self.alt[idhdf],label='Retrieved')
            ax.fill_betweenx(self.alt[idhdf],prfMean-prfSTD,prfMean+prfSTD,alpha=0.25)

            #dtpMean    = np.nanmean(self.dtp[idhdf])
            #dtpStd     = np.nanstd(self.dtp[idhdf])

            #ax.axhline(y=dtpMean, linewidth=1.5, color='gray', alpha=0.5)  
            #ax.axhline(y=dtpMean - dtpStd*2., linewidth=1.5, color='gray', alpha=0.5)  
            #ax.axhline(y=dtpMean + dtpStd*2, linewidth=1.5, color='gray', alpha=0.5)

            #altftmean  =   np.nanmean(self.alttpp[idhdf])
            #altstmean  =   np.nanmean(self.alttpp2[idhdf])

            #ax.axhline(y=altftmean, linewidth=1.5, color='green', alpha=0.5) 
            #ax.axhline(y=altstmean, linewidth=1.5, color='green', alpha=0.5)


            ax.plot(self.aPrf[idhdf][0],self.alt[idhdf],linestyle='--', label='-Apriori', color='r')

            #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
            ax.annotate(self.locsID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat)), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

            ax.grid(True,which='both', alpha=0.5)    
            #ax.legend(prop={'size':10})    
            #ax.set_ylabel('Altitude [km]')
            ax.tick_params(which='both',labelsize=10)
            ax.set_ylim(0,50)
            ax.set_xlim(0,600) 
            #ax.set_title('Mean profile of'+ gasName.upper())
          
            fig.add_subplot(ax)

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            # for sp in ax.spines.values():
            #     sp.set_visible(False)
            #     plt.setp(ax.get_xticklabels(), visible=False)
            # if ax.is_first_row():
            #     ax.spines['top'].set_visible(True)
            # if ax.is_last_row():
            #     ax.spines['bottom'].set_visible(True)
            #     plt.setp(ax.get_xticklabels(), visible=True)
            # if ax.is_first_col():
            #     ax.spines['left'].set_visible(True)
            # if ax.is_last_col():
            #     ax.spines['right'].set_visible(True)

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)

                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)

        if (self.npanels % 2 == 1): #even

            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)


        all_axes[-3].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-3].get_xticklabels(), visible=True)
        all_axes[-3].set_zorder(1)

        all_axes[-4].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-4].get_xticklabels(), visible=True)
        all_axes[-4].set_zorder(1)

        all_axes[-1].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-2].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-3].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-4].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)


        fig.text(0.01, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
        #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
        fig.subplots_adjust(left=0.06, bottom=0.05, right=0.975, top=0.975)

        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            fig.savefig('/data1/projects/ocs/figures/fig/'+'Profiles.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)

        
    def pltPrfErr(self):

        print '\n Plotting error profiles.......\n'
        #---------------------------------------------------
        # Mean Error vertical profiles
        #---------------------------------------------------
        fig = plt.figure(figsize=(12,12))

        outer_grid = gridspec.GridSpec(self.npanels2, 4, wspace=0.11, hspace=0.085)

        for i, idhdf in enumerate(self.pltID):

            ax = plt.Subplot(fig, outer_grid[i])

            Lat = float(self.Lat[i])

            if Lat >= 50.:     
                ax.set_facecolor('lightcyan')
                
            elif (Lat >= 20.) & (Lat < 50.):
                #ax.set_axis_bgcolor('lightgreen')
                ax.set_facecolor('lightgreen')
               
            elif (Lat >= -20.)  & (Lat < 20.):
                ax.set_facecolor('mistyrose')
                
            elif (Lat >= -50.)  & (Lat < -20.):
                ax.set_facecolor('cornsilk')
                
            elif (Lat < -50.):
                ax.set_facecolor('lightgrey')
            
            else:
                ax.set_facecolor('lightgrey')       

            vmr_rnd_err_Mean  = np.nanmean(self.vmr_rnd_err[idhdf],axis=0)   
            vmr_sys_err_Mean  = np.nanmean(self.vmr_sys_err[idhdf],axis=0)
            vmr_tot_err_Mean  = np.nanmean(self.vmr_tot_err[idhdf],axis=0)

            vmr_tot_err_std   = np.nanstd(self.vmr_tot_err[idhdf],axis=0)

           
            ax.plot(vmr_rnd_err_Mean,self.alt[idhdf],label='Random Error', color='r')
            ax.plot(vmr_sys_err_Mean,self.alt[idhdf],label='Systematic Error', color='b')
            ax.plot(vmr_tot_err_Mean,self.alt[idhdf],label='Total Error', color='k')
                #ax.fill_betweenx(alt[idhdf],vmr_tot_err_Mean-vmr_tot_err_std,vmr_tot_err_Mean+vmr_tot_err_std,alpha=0.25, color='k')
            

            #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
            ax.annotate(self.locsID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat)), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')

            ax.grid(True,which='both')    
            #ax.legend(prop={'size':10})    
            #ax.set_ylabel('Altitude [km]')
            ax.tick_params(which='both',labelsize=10)
            ax.set_ylim(0,50)
            ax.set_xlim(0,60) 
            #ax.set_title('Mean profile of'+ gasName.upper())
          
            fig.add_subplot(ax)

        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:
            # for sp in ax.spines.values():
            #     sp.set_visible(False)
            #     plt.setp(ax.get_xticklabels(), visible=False)
            # if ax.is_first_row():
            #     ax.spines['top'].set_visible(True)
            # if ax.is_last_row():
            #     ax.spines['bottom'].set_visible(True)
            #     plt.setp(ax.get_xticklabels(), visible=True)
            # if ax.is_first_col():
            #     ax.spines['left'].set_visible(True)
            # if ax.is_last_col():
            #     ax.spines['right'].set_visible(True)

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
                if ax.is_last_row():
                    ax.spines['bottom'].set_visible(True)
                    plt.setp(ax.get_xticklabels(), visible=True)

                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)

        if (npanels % 2 == 1): #even

            all_axes[-2].spines['bottom'].set_visible(True)
            plt.setp(all_axes[-2].get_xticklabels(), visible=True)
            all_axes[-2].set_zorder(1)


        all_axes[-3].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-3].get_xticklabels(), visible=True)
        all_axes[-3].set_zorder(1)

        all_axes[-4].spines['bottom'].set_visible(True)
        plt.setp(all_axes[-4].get_xticklabels(), visible=True)
        all_axes[-4].set_zorder(1)

        all_axes[-1].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-2].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-3].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)
        all_axes[-4].set_xlabel('VMR ['+self.sclfctName+']', fontsize=14)


        fig.text(0.01, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle('{} Vertical Profiles'.format(gasName.upper()), fontsize=16  )
        #plt.tight_layout(h_pad=0.25) #w_pad=1.75 pad=1.75,
        fig.subplots_adjust(left=0.06, bottom=0.05, right=0.975, top=0.975)

        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            fig.savefig('/data1/projects/ocs/figures/fig/'+'Profiles_Error.pdf', bbox_inches='tight')

        else: 
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()


    def pltAK(self, pltWvmr=False):

        print '\n Plotting AKs.......\n'

        #---------------------------------------------------
        # Averaging Kernel Smoothing Function (row of avk)
        #--------------------------------------------------
        fig = plt.figure(figsize=(12,12))
        outer_grid = gridspec.GridSpec(self.npanels2, 4, wspace=0.12, hspace=0.085)

        for i, idhdf in enumerate(self.pltID):

            indsAlt   = np.where(self.alt[idhdf] < 50.)[0]

            if pltWvmr: pltAKAv   = np.nanmean(self.avkVMR[idhdf], axis=0)
            else: pltAKAv   = np.nanmean(self.avkSCF[idhdf], axis=0)

            #if pltWvmr: pltAKAv   = avkVMR[idhdf][indsAlt]
            #else: pltAKAv   = avkSCF[idhdf][indsAlt]
            
            avkSCFAv   = np.nanmean(self.avkSCF[idhdf], axis=0)

            #---------
            #Total Column AK
            #---------
            avkTCAv    = np.nanmean(self.avkTC[idhdf], axis=0)
            
            #dtpMean    = np.nanmean(self.dtp[idhdf])
            #dtpStd     = np.nanstd(self.dtp[idhdf])
            #Pcol       = [dtpMean, 120.]

            Lat = float(self.Lat[i])
            
            gs        = gridspec.GridSpecFromSubplotSpec(1,3,subplot_spec=outer_grid[i], width_ratios=[3,1,1])
            axa       = plt.subplot(gs[0])
            axb       = plt.subplot(gs[1])
            axc       = plt.subplot(gs[2])
            cm        = plt.get_cmap(self.clmap)
            cNorm     = colors.Normalize(vmin=np.min(self.alt[idhdf][indsAlt]), vmax=np.max(self.alt[idhdf][indsAlt]))
            
            scalarMap = mplcm.ScalarMappable(norm=cNorm,cmap=self.clmap)
            scalarMap.set_array(self.alt[idhdf][indsAlt])
            axa.set_color_cycle([scalarMap.to_rgba(x) for x in self.alt[idhdf][indsAlt]])

            if Lat >= 50.:     
                axa.set_facecolor('lightcyan')
                axb.set_facecolor('lightcyan')
                axc.set_facecolor('lightcyan')
                
            elif (Lat >= 20.) & (Lat < 50.):
                #ax.set_axis_bgcolor('lightgreen')
                axa.set_facecolor('lightgreen')
                axb.set_facecolor('lightgreen')
                axc.set_facecolor('lightgreen')
               
            elif (Lat >= -20.)  & (Lat < 20.):
                axa.set_facecolor('mistyrose')
                axb.set_facecolor('mistyrose')
                axc.set_facecolor('mistyrose')
                
            elif (Lat >= -50.)  & (Lat < -20.):
                axa.set_facecolor('cornsilk')
                axb.set_facecolor('cornsilk')
                axc.set_facecolor('cornsilk')
                
            elif (Lat < -50.):
                axa.set_facecolor('lightgrey')
                axb.set_facecolor('lightgrey')
                axc.set_facecolor('lightgrey')
            
            else:
                axa.set_facecolor('lightgrey')
                axb.set_facecolor('lightgrey')
                axc.set_facecolor('lightgrey')

            for j in range(len(self.alt[idhdf][indsAlt])):
                axa.plot(pltAKAv[indsAlt[j],indsAlt],self.alt[idhdf][indsAlt])
                #ax.plot(avkSCFav[gasVer][indAlt[i],indAlt],alt[gasVer][indAlt])

            #axa.axhline(y=dtpMean, linewidth=1.5, color='r')  
            #axa.axhline(y=dtpMean - dtpStd, linewidth=1.5, color='r', alpha=0.5)  
            #axa.axhline(y=dtpMean + dtpStd, linewidth=1.5, color='r', alpha=0.5)  
            axa.set_ylim(0, 40)
            axa.set_xlim(-0.15, 0.25)
            axa.xaxis.set_ticks(np.arange(-0.15, 0.25, 0.1))
            axa.tick_params(axis = 'both', which = 'major', labelsize = 8)
            axa.tick_params(axis = 'both', which = 'minor', labelsize = 0)

            if pltWvmr:
                if i == self.npanels-1: axa.set_xlabel('VMR AK')
                if i == self.npanels-2: axa.set_xlabel('VMR AK')
                if i == self.npanels-3: axa.set_xlabel('VMR AK')
                if i == self.npanels-4: axa.set_xlabel('VMR AK')
            else:
                if i == self.npanels-1: axa.set_xlabel('AK')
                if i == self.npanels-2: axa.set_xlabel('AK')
                if i == self.npanels-3: axa.set_xlabel('AK')
                if i == self.npanels-4: axa.set_xlabel('AK')

            axa.annotate(self.locsID[i].upper(), xy=(0.025, 0.875), xycoords='axes fraction', fontsize=10, ha='left')

            #---------------------------------------------------
            # Total Column AK
            #---------------------------------------------------
            axb.plot(avkTCAv,self.alt[idhdf],color='k')
            #axb.axhline(y=dtpMean, linewidth=1.5, color='r')  
            #axb.axhline(y=dtpMean - dtpStd, linewidth=1.5, color='r', alpha=0.5)  
            #axb.axhline(y=dtpMean + dtpStd, linewidth=1.5, color='r', alpha=0.5)
            axb.set_ylim(0, 40)
            axb.set_xlim(0, 2)
            axb.xaxis.set_ticks(np.arange(0,2,1))
            axb.tick_params(axis = 'both', which = 'major', labelsize = 8)
            axb.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            if i == self.npanels-1: axb.set_xlabel('TC AK')
            if i == self.npanels-2: axb.set_xlabel('TC AK')
            if i == self.npanels-3: axb.set_xlabel('TC AK')
            if i == self.npanels-4: axb.set_xlabel('TC AK')

            #---------------------------------------------------
            # Cum Sum of DOF
            #---------------------------------------------------
            if self.alt[idhdf][indsAlt][0] > self.alt[idhdf][indsAlt][-1]:
                dofs_cs = np.cumsum(np.diag(avkSCFAv)[::-1])[::-1]
            else:
                dofs_cs = np.cumsum(np.diag(avkSCFAv))

            axc.plot(dofs_cs,self.alt[idhdf],color='k',label='Cumulative Sum of DOFS (starting at surface)')

            #axc.axhline(y=dtpMean, linewidth=1.5, color='r')  
            #axc.axhline(y=dtpMean - dtpStd, linewidth=1.5, color='r', alpha=0.5)  
            #axc.axhline(y=dtpMean + dtpStd, linewidth=1.5, color='r', alpha=0.5) 
            
            try:
                xval = range(0,int(np.ceil(max(dofs_cs)))+3)

            except Exception as errmsg:
                print '\nError: ', errmsg
                xval = range(0, 2)

            #ind1         = mf.nearestind(Pcol[0], self.alt[idhdf])
            #ind2         = mf.nearestind(Pcol[1], self.alt[idhdf])

            #dofsPcol = dofs_cs[ind2] - dofs_cs[ind1]
        
            axc.set_ylim(0, 40)
            axc.set_xlim(0, 3)
            axc.xaxis.set_ticks(np.arange(0,3,1))
            axc.tick_params(axis='x',which='both',labelsize=8)
            axc.tick_params(axis = 'both', which = 'major', labelsize = 8)
            axc.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            if i == self.npanels-1: axc.set_xlabel('Cum\nsum of DOFS')
            if i == self.npanels-2: axc.set_xlabel('Cum\nsum of DOFS')
            if i == self.npanels-3: axc.set_xlabel('Cum\nsum of DOFS')
            if i == self.npanels-4: axc.set_xlabel('Cum\nsum of DOFS')

            print idhdf
            #print 'DOFs for layer {0:.1f}-{1:.1f}[km] = {2:.2f}'.format(alt[idhdf][ind1],alt[idhdf][ind2],dofsPcol)
            print 'DOFs all = {0:.2f}'.format(np.trace(avkSCFAv))
          
            fig.add_subplot(axa)
            fig.add_subplot(axb)
            fig.add_subplot(axc)  

        all_axes = fig.get_axes()
        #show only the outside spines
        for i, ax in enumerate(all_axes):

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)

        for i in range(-1, -13, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True)
            all_axes[i].set_zorder(1)

        for i in range(0, 64, 12):
            all_axes[i].spines['left'].set_visible(True)
            plt.setp(all_axes[i].get_yticklabels(), visible=True)
            all_axes[i].set_zorder(1)

        cbaxes = fig.add_axes([0.55, 0.1, 0.4, 0.015]) 
        cbar   = fig.colorbar(scalarMap, orientation='horizontal', cax = cbaxes)
        cbar.set_label('Altitude [km]', fontsize=14)


        fig.text(0.01, 0.5, 'Altitude [km]', fontsize=16, va='center', rotation='vertical')
        fig.subplots_adjust(left=0.06, bottom=0.05, right=0.975, top=0.975)

        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            if pltWvmr: fig.savefig('/data1/projects/ocs/figures/fig/AK_VMR.pdf', bbox_inches='tight')
            else: fig.savefig('/data1/projects/ocs/figures/fig/AK_SC.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)

    def pltMnthlyVMR(self):

        #---------------------------------------------------
        # Plot total Weighted VMR by month
        #---------------------------------------------------
        print '\nPlot: Monthly mean Partial Columns:\n' 
        fig, ax  = plt.subplots(3, 5, figsize=(17, 10), sharey=False, sharex=True)
        
        for i, idhdf in enumerate(self.pltID2):

            month    = np.array([d.month for d in self.dates2[idhdf]])
            mnthSort = list(set(month))
            
            mnthMean = np.zeros(len(mnthSort))
            mnthSTD  = np.zeros(len(mnthSort))

            mnthMean2 = np.zeros(len(mnthSort))
            mnthSTD2  = np.zeros(len(mnthSort))

            mnthMean3 = np.zeros(len(mnthSort))
            mnthSTD3  = np.zeros(len(mnthSort))
            
            for j,m in enumerate(mnthSort):
                inds        = np.where(month == m)[0]
                mnthMean[j] = np.mean(self.WvmrTrop12[idhdf][inds])
                mnthSTD[j]  = np.std(self.WvmrTrop12[idhdf][inds])

                mnthMean2[j] = np.mean(self.WvmrTrop22[idhdf][inds])
                mnthSTD2[j]  = np.std(self.WvmrTrop22[idhdf][inds])

                mnthMean3[j] = np.mean(self.WvmrStrat2[idhdf][inds])
                mnthSTD3[j]  = np.std(self.WvmrStrat2[idhdf][inds])

            Lat2 = float(self.Lat2[i])

            if Lat2 >= 50.:     
                ax[2,0].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                ax[1,0].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                ax[0,0].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                ax[0,0].grid(True,which='both')
                ax[1,0].grid(True,which='both')
                ax[2,0].grid(True,which='both')

                ax[2,0].set_ylabel('Low Tropospheric wVMR ['+self.sclfctName+']',multialignment='center')
                ax[1,0].set_ylabel('Free Tropospheric wVMR ['+self.sclfctName+']',multialignment='center')
                ax[0,0].set_ylabel('Stratospheric wVMR ['+self.sclfctName+']',multialignment='center')

                ax[0,0].set_xlim((0,13))
                ax[0,0].set_xticks(range(1,13))
                ax[1,0].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                ax[2,0].set_xlabel('Month')
                ax[0, 0].set_title('50$^{\circ}$N - 90$^{\circ}$N', fontsize=14)

                ax[2, 0].set_ylim(350, 580)
                ax[1, 0].set_ylim(350, 580)
                ax[0, 0].set_ylim(100, 450)
                
            elif (Lat2 >= 20.) & (Lat2 < 50.):
                ax[2,1].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                ax[1,1].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                ax[0,1].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                ax[0,1].grid(True,which='both')
                ax[1,1].grid(True,which='both')
                ax[2,1].grid(True,which='both')

                ax[1,1].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                ax[2,1].set_xlabel('Month')

                ax[0,1].set_title('20$^{\circ}$N - 50$^{\circ}$N', fontsize=14)

                ax[2, 1].set_ylim(350, 580)
                ax[1, 1].set_ylim(350, 580)
                ax[0, 1].set_ylim(100, 450)
               
            elif (Lat2 >= -20.)  & (Lat2 < 20.):
                ax[2,2].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                ax[1,2].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                ax[0,2].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                ax[0,2].grid(True,which='both')
                ax[1,2].grid(True,which='both')
                ax[2,2].grid(True,which='both')

                ax[1,2].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                ax[2,2].set_xlabel('Month')
                ax[0,2].set_title('20$^{\circ}$S - 20$^{\circ}$N', fontsize=14)

                ax[2, 2].set_ylim(350, 580)
                ax[1, 2].set_ylim(350, 580)
                ax[0, 2].set_ylim(100, 450)
                
            elif (Lat2 >= -50.)  & (Lat2 < -20.):
                ax[2,3].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                ax[1,3].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                ax[0,3].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                ax[0,3].grid(True,which='both')
                ax[1,3].grid(True,which='both')
                ax[2,3].grid(True,which='both')

                ax[1,3].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                ax[2,3].set_xlabel('Month')

                ax[0,3].set_title('50$^{\circ}$S - 20$^{\circ}$S', fontsize=14)

                ax[2, 3].set_ylim(350, 580)
                ax[1, 3].set_ylim(350, 580)
                ax[0, 3].set_ylim(100, 450)

                
            elif (Lat2 < -50.):
                ax[2,4].errorbar(mnthSort,mnthMean, fmt='o', yerr=mnthSTD ,markersize=6)
                ax[1,4].errorbar(mnthSort,mnthMean2, fmt='o', yerr=mnthSTD2 ,markersize=6, label=idhdf)
                ax[0,4].errorbar(mnthSort,mnthMean3, fmt='o', yerr=mnthSTD3 ,markersize=6)

                ax[0,4].grid(True,which='both')
                ax[1,4].grid(True,which='both')
                ax[2,4].grid(True,which='both')

                ax[1,4].legend(ncol=2, prop={'size':8}, loc = 'upper center', bbox_to_anchor=[0.5, 1.2])

                ax[2,4].set_xlabel('Month')
                ax[0,4].set_title('90$^{\circ}$S - 50$^{\circ}$S', fontsize=14)

                ax[2, 4].set_ylim(350, 580)
                ax[1, 4].set_ylim(350, 580)
                ax[0, 4].set_ylim(100, 450)

        fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.98, wspace=0.15, hspace=0.15)
        
        if self.pdfsav: 
            self.pdfsav.savefig(fig,dpi=200)
            fig.savefig('/data1/projects/ocs/figures/fig/mnthly_wVMR2.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)
            #user_input = raw_input('Press any key to exit >>> ')
            #sys.exit()

        

def AnalTSAnom(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1, qboFlg=False, xmin=False, xmax=False):
    
    #--------------------
    #Slope and Amplitude of time series
    #--------------------
    slope       = []   #slope
    slope_e     = []   #slope error

    slope_p1    = []   #slope 1995 - 2002
    slope_p1_e  = []   #slope error 

    slope_p2    = []   #slope 2002 - 2008
    slope_p2_e  = []   #slope error

    slope_p3    = []   #slope  2008 - 2016
    slope_p3_e  = []   #slope error

    amp         = []   #Amplitude 

    avgD        = []
    stdD        = []

    avgD_p1     = []
    stdD_p1     = []

    avgD_p2     = []
    stdD_p2     = []

    avgD_p3     = []
    stdD_p3     = []

    #xmin      = dt.date(1993, 1, 1)
    #xmax      = dt.date(2018, 12, 31)


    if pltFig:    
        
        fig  = plt.figure(figsize=(18,13))
        fig2 = plt.figure(figsize=(18,13))  
        
        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.05, hspace=0.075)

    fi = 0

    for i, idhdf in enumerate(pltID):

        ax   = plt.Subplot(fig, outer_grid[i])
        ax2  = plt.Subplot(fig2, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_facecolor('lightcyan')
            ax2.set_facecolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_facecolor('lightgreen')
            ax2.set_facecolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_facecolor('mistyrose')
            ax2.set_facecolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_facecolor('cornsilk')
            ax2.set_facecolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_facecolor('lightgrey')
            ax2.set_facecolor('lightgrey')
        
        else:
            ax.set_facecolor('lightgrey')
            ax2.set_facecolor('lightgrey')
        
        fig.add_subplot(ax)

        fig2.add_subplot(ax2)

        #----------------------------
        if AvgType == 'Daily':
            Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      = Avg['dailyAvg']

            
        elif AvgType == 'Monthly':
            Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      =  Avg['mnthlyAvg']

        elif AvgType == 'none':
            Dates        = xDates[idhdf]
            dateYearFrac = mf.toYearFraction(Dates)
            AvgData      = yData[idhdf]

        else:
            print 'Error: Define average type: Daily, Monthly, or none'
            exit()

        indsnan = [ni for ni, n in enumerate(AvgData) if np.isnan(n)]

        if len(indsnan) == len(AvgData):

            slope.append(float('nan'))
            slope_e.append(float('nan'))

            slope_p1.append(float('nan'))
            slope_p1_e.append(float('nan'))

            slope_p2.append(float('nan'))
            slope_p2_e.append(float('nan'))

            slope_p3.append(float('nan'))
            slope_p3_e.append(float('nan'))

            amp.append(float('nan'))

            avgD.append(float('nan'))
            stdD.append(float('nan'))

            avgD_p1.append(float('nan'))
            stdD_p1.append(float('nan'))

            avgD_p2.append(float('nan'))
            stdD_p2.append(float('nan'))

            avgD_p3.append(float('nan'))
            stdD_p3.append(float('nan'))


        else:
            #continue

    
            #---------------------------------------------------
            #To make a continuous fit
            #---------------------------------------------------
            numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
            dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
            dates2  = np.asarray(dates2)
            dateYearFrac2 = mf.toYearFraction(dates2)
            #---------------------------------------------------

            yyyy = [sngdate.year for sngdate in Dates]
            yyyy = np.asarray(yyyy)

            months = [sngdate.month for sngdate in Dates]
            months = np.asarray(months)


            #--------------------
            # Apply savitzky golay filter (Comment out if not wanted)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            #--------------------
            # Anomalies
            #--------------------
            if AvgType == 'Monthly':

                weights      = np.ones_like(dateYearFrac)

                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                AnnualCy        = f_fourier(dateYearFrac)
                Anomaly         = AvgData - AnnualCy


                # Anomaly  = np.zeros(len(Dates))

                # month    = np.array([d.month for d in Dates])
                # mnthSort = list(set(month))
                
                # mnthMean = np.zeros(len(mnthSort))
                # mnthSTD  = np.zeros(len(mnthSort))
                
                # for j,m in enumerate(mnthSort):
                #     inds        = np.where(month == m)[0]
                #     mnthMean[j] = np.mean(AvgData[inds])
                #     mnthSTD[j]  = np.std(AvgData[inds])

                #     Anomaly[inds]   = AvgData[inds]  - mnthMean[j]

            if AvgType == 'Daily':
                weights      = np.ones_like(dateYearFrac)

                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                AnnualCy        = f_fourier(dateYearFrac)
                Anomaly         = AvgData - AnnualCy

            #--------------------
            # Removing Additional Modulation (QBO, solar cycle, etc) --> TO BE IMPROVED
            #--------------------
            if qboFlg:
                weights      = np.ones_like(dateYearFrac)

                #res          = mf.fit_driftfourier_qbo(dateYearFrac, Anomaly, weights, 2, half_period=2.3)
                res          = mf.fit_qbo(dateYearFrac, Anomaly+300., weights, 2, half_period=2.33)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                Mdl         = f_fourier(dateYearFrac)
                Anomaly2    = Anomaly - Mdl + 300.

            else:

                Anomaly2    = Anomaly

            #--------------------
            # Plot Initial Anomalies (removing seasonal cycle)
            #--------------------
            if pltFig:
                
                ax.plot(Dates, Anomaly,'k.',markersize=0)
                ax.scatter(Dates, Anomaly, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)
                #print len(Dates)
                #print len(f_fourier(dateYeardateYearFracFrac))
                if qboFlg: ax.plot(dates2,f_fourier(dateYearFrac2)-300.,label='Fitted Anual Trend',linewidth=2.0)

                ax.axhline(y=0, linestyle='--', linewidth=1.5, color='gray', alpha=0.5)


                #ax2 = plt.Subplot(fig2, outer_grid[fi])
                fi += 1
                ax2.plot(Dates, Anomaly2,'k.',markersize=0)
                ax2.scatter(Dates, Anomaly2, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                ax2.axhline(y=0, linestyle='--', linewidth=1.5, color='gray', alpha=0.5)

                #--------------------
                #
                #--------------------

                if qboFlg:
                    qboDir        = '/data1/projects/ocs/'                 # Name of station location
                    qbofile       = 'QBO.dat' 

                    cols, indexToName = mf.getColumns(qboDir + qbofile, headerrow=8, delim=' ', header=True)

                    YYMM  = cols[indexToName[1]]
                    qbo50 = cols[indexToName[3]]
                    qbo30 = cols[indexToName[5]]

                    qboDate = []

                    for q in YYMM:
                        YY = q[0:2]
                        MM = q[2:4]

                        if int(YY) <= 20 :
                            YYYY = '20'+YY
                        else:
                            YYYY = '19'+YY

                        qboDate.append(dt.date(int(YYYY), int(MM), 15))

                    qboDate    = np.asarray(qboDate)
                    qbo50   = np.asarray(qbo50, dtype=float)
                    qbo30   = np.asarray(qbo30, dtype=float)

                    axx = ax.twinx()
                    axx.plot(qboDate, qbo30, color= 'gray')
                #fig.add_subplot(axx)

            #--------------------

            if fits:
                
                #---------------------------------------------------
                #LINEAR FIT AND DRIFT FOR ALL YEARS
                #---------------------------------------------------
                #if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Paris'):
                yoi    = [1996, 2018]

                if idhdf == 'AHTS':
                    

                    indx1  = np.where( (yyyy >= yoi[0]) &  (yyyy <= yoi[1]) & (months >= 1) & (months <= 12) )[0]

                else:
                    indx1  = np.where( (yyyy >= yoi[0]) &  (yyyy <= yoi[1]))[0]

                weights      = np.ones_like(dateYearFrac)

                res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A, df_drift = res[3:9]

                res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yoi[0], yoi[-1])

                if (idhdf == 'Zugspitze') or (idhdf == 'Jungfraujoch') or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Thule'):

                    #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                    

                    slope.append(res[1]/np.mean(AvgData[indx1])*100.0)
                    slope_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                    #print res[1]/np.mean(AvgData)*100.0

                    Amp   = np.sum(res[2]**2)
                    Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                    amp.append(Amp)

                    avgD.append(np.mean(AvgData[indx1]))
                    stdD.append(np.std(AvgData[indx1]))

                    ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color='k')

                else:

                    slope.append(float('nan'))
                    slope_e.append(float('nan'))

                    avgD.append(float('nan'))
                    stdD.append(float('nan'))

                

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, Anomaly, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, Anomaly, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                #---------------------------------------------------
                # POLY FIT FOR STATIONS WITH LONG TIME SERIES (TO KNOW INFLECTION POINTS)
                #---------------------------------------------------
                if (idhdf != 'Eureka') & (idhdf != 'St Petersburg') & (idhdf != 'Boulder') & (idhdf != 'Maido') &  (idhdf != 'St Denis') & (idhdf != 'Bremen') & (idhdf != 'Paris') & (idhdf != 'Altzomoni'):

                    res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly2, weights, 2, half_period=period)
                    f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]

            #---------------------------------------------------
            # Finding local min & Max (decrease or increase rate of change)
            #---------------------------------------------------
                    if fits: 
                        #roc  = df_drift(dateYearFrac2)/np.mean(AvgData)*100.0
                        #roc  = np.asarray(roc)
                        #zero_crossings = np.where(np.diff(np.sign(roc)))[0]

                        a = diff(sign(diff(f_driftfourier(dateYearFrac2)))).nonzero()[0] + 1 # local min+max
                        b = (diff(sign(diff(f_driftfourier(dateYearFrac2)))) > 0).nonzero()[0] + 1 # local min
                        c = (diff(sign(diff(f_driftfourier(dateYearFrac2)))) < 0).nonzero()[0] + 1 # local max
                        
                        #ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                        #ax.plot(dates2[a], f_driftfourier(dateYearFrac2)[a], "o", color='green')
                        ##ax.plot(dates2[c], f_driftfourier(dateYearFrac2)[c],  "o", color='red')

                
                        print idhdf
                        print 'Number of crossings = {}'.format(len(dates2[a]))
                        print 'Dates of crossings (minimums)  = {}'.format(dates2[b])
                        print 'Dates of crossings (maximums)  = {}'.format(dates2[c])

                # if pltFig:
                #     #ax.plot(dates2,f_fourier(dateYearFrac2),label='Fitted Anual Trend', linewidth=2.0)
                #     ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #     ax2.plot(dates2,f_drift(dateYearFrac2),label='Fitted Anual Trend',linewidth=2.0)

                #     if idhdf.lower() == 'jungfraujoch':

                #         fig, ax = plt.subplots(2, figsize=(10,6), sharex=True)
                #         ax[0].plot(dates2, f_driftfourier(dateYearFrac2), linewidth=2.0)
                #         ax[0].plot(dates2[b], f_driftfourier(dateYearFrac2)[b], "o", color='green')
                #         ax[0].plot(dates2[c], f_driftfourier(dateYearFrac2)[c],  "o", color='red')

                #         ax[1].plot(dates2, df_drift(dateYearFrac2), linewidth=2.0)

                #         plt.show(block= False)

                #         user_input = raw_input('Press any key to exit >>> ')
                #         sys.exit()   
                
            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    
                # ax = plt.Subplot(fig, outer_grid[i])
                # ax.plot(Dates, AvgData,'k.',markersize=0)
                # ax.scatter(Dates,AvgData, s=10, facecolor='lightgray', edgecolor='k',alpha=0.85)

                # if yData2:
                #     ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                # if yData3:
                #     ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                # if yData4:
                #     dtpStd      = np.nanstd(yData4[idhdf])
                #     dtpMean     = np.nanmean(yData4[idhdf])

                #     ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                #     ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  \

                

                


                #---------------------------------------------------


                # if (idhdf == 'Paris') or (idhdf == 'Altzomoni'):


                #     slope_p1.append(float('nan'))
                #     slope_p1_e.append(float('nan'))

                #     slope_p2.append(float('nan'))
                #     slope_p2_e.append(float('nan'))

                #     slope_p3.append(float('nan'))
                #     slope_p3_e.append(float('nan'))

                #     avgD_p1.append(float('nan'))
                #     stdD_p1.append(float('nan'))

                #     avgD_p2.append(float('nan'))
                #     stdD_p2.append(float('nan'))

                #     avgD_p3.append(float('nan'))
                #     stdD_p3.append(float('nan'))



                if (idhdf == 'Eureka')  or (idhdf == 'St Petersburg') or (idhdf == 'Bremen') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Tsukuba') or (idhdf == 'Paramaribo') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                    yoi  = [[2009, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color='red')


                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund') or (idhdf == 'Izana'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    clr_yoi = ['green', 'red']

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
        
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                      
                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])
                  
                        ##ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Rikubetsu') :

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))

                    clr_yoi = ['blue', 'green']
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])



                elif ( (idhdf == 'Jungfraujoch')  or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze') or (idhdf == 'Mauna Loa') or (idhdf == 'Kiruna') ):

                    yoi  = [[1996, 2002], [2002, 2008], [2009, 2018]]

                    clr_yoi = ['blue', 'green', 'red']
                
                    for ii, y in enumerate(yoi):

                        
                        if idhdf == 'AHTS':
                        
                            indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]) & (months >= 1) & (months <= 12) )[0]

                            #ax.plot(Dates[indx1], Anomaly[indx1],'k.',markersize=0)
                            #ax.scatter(Dates[indx1], Anomaly[indx1], s=40, facecolor='red', edgecolor='k',alpha=0.85)

                        else:
    
                            indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        #res          = mf.fit_driftfourier_poly(dateYearFrac, Anomaly, weights, 2, half_period=period)
                        #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], Anomaly2[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                        ax.plot(Dates[indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend',linewidth=2.0, color=clr_yoi[ii])

                        
                    ##ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)


        if pltFig:
            if (xmin and xmax): ax.set_xlim(xmin, xmax)

            ax.grid(True, color='gray', alpha=0.25)
            ax.tick_params(which='both',labelsize=10)
            #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
            ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
            #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
            #start, end = ax1[i].get_xlim()
            #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
            #ax1[i].set_ylim(bottom=0)

            yearsLc1      = YearLocator(2)
            yearsLc2      = YearLocator(1)
            months        = MonthLocator()
            DateFmt       = DateFormatter('%Y')

            #plt.xticks(rotation=45)
            ax.xaxis.set_major_locator(yearsLc1)
            ax.xaxis.set_minor_locator(yearsLc2)
            #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax.xaxis.set_major_formatter(DateFmt)
            #ax.set_xlabel('Year')
            #ax1.xaxis.set_tick_params(which='major', pad=15)  
            ax.xaxis.set_tick_params(which='minor',labelbottom='off')
            ax.tick_params(which='both')

            if (ymin and ymax):
                ax.set_ylim(ymin, ymax)
          
            #ax.set_xticks([])
            #ax.set_yticks([])
            #fig.add_subplot(ax)


            if (xmin and xmax): ax2.set_xlim(xmin, xmax)
            #ax2.set_ylim(-7, 7)
            #ax2.axhline(y=0, linestyle='--', color='k')
            ax2.grid(True, color='gray', alpha=0.25)
            ax2.tick_params(which='both',labelsize=10)
            ax2.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')

            ax2.xaxis.set_major_locator(yearsLc1)
            ax2.xaxis.set_minor_locator(yearsLc2)
            #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
            #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
            ax2.xaxis.set_major_formatter(DateFmt)
            #ax.set_xlabel('Year')
            #ax1.xaxis.set_tick_params(which='major', pad=15)  
            ax2.xaxis.set_tick_params(which='minor',labelbottom='off')
            ax2.tick_params(which='both')

        

            if (ymin and ymax):
                ax2.set_ylim(ymin, ymax)
    
            #fig2.add_subplot(ax2)


            #---------------------------------------------------

    if pltFig:
    
        # all_axes = fig.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')
        
        # fig.autofmt_xdate()
        # fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        # fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        # fig.suptitle(ytypeStr, fontsize=16  )
        
        #-----
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
            
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)

        
        for i in range(-1, -5, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        ytypeStr = ytypeStr.replace(' ', '')

        #fig.text(0.02, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

       

        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
            else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)
        

    #---------------------------------------------------
    #
    #---------------------------------------------------
    if pltFig:
    
        # all_axes = fig2.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')
        
        # fig2.autofmt_xdate()
        # fig2.text(0.03, 0.5, 'Rate of change [%/y]', fontsize=16, va='center', rotation='vertical')
        # fig2.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)
        # fig2.suptitle('Rate of change of VMR Anomalies', fontsize=16  )

        #-----
        all_axes = fig2.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
            
                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)

        for i in range(-1, -4, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        ytypeStr = ytypeStr.replace(' ', '')

        #fig2.autofmt_xdate()

        #fig2.text(0.02, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #fig2.text(0.0075, 0.5, 'Rate of change [%/y]', fontsize=16, va='center', rotation='vertical')
        fig2.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig2.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

        
        if saveFlg: pdfsav.savefig(fig2,dpi=200)
        else: 
            plt.show(block=False)

        

        #if fits: fig2.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'ROC_wFit.pdf', bbox_inches='tight')
        #else: fig2.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'ROC.pdf', bbox_inches='tight')


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)   

        
def AnalTS(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1, xmin=False, xmax=False):
    
    #--------------------
    #Slope and Amplitude of time series
    #--------------------
    slope       = []   #slope
    slope_e     = []   #slope error

    slope_p1    = []   #slope 1995 - 2002
    slope_p1_e  = []   #slope error 

    slope_p2    = []   #slope 2002 - 2008
    slope_p2_e  = []   #slope error

    slope_p3    = []   #slope  2008 - 2016
    slope_p3_e  = []   #slope error

    amp         = []   #Amplitude 

    avgD        = []
    stdD        = []

    avgD_p1     = []
    stdD_p1     = []

    avgD_p2     = []
    stdD_p2     = []

    avgD_p3     = []
    stdD_p3     = []

    #xmin      = dt.date(1993, 1, 1)
    #xmax      = dt.date(2018, 12, 31)

    if pltFig:    
        fig = plt.figure(figsize=(18,13))  

        outer_grid = gridspec.GridSpec(npanels, 3, wspace=0.05, hspace=0.075)
    #--------------
    # with open('/data1/projects/ocs/figures/fig/'+ytypeStr+'.txt','w') as fopen:
    #     fopen.write('#Site\tLatitude\tMean TPH\tstd TPH\tMax TPH\tMin TPH\tAmplitude\n')

    
    #     strformat = '{:<15}\t'+'\t'.join('{:>5.1f}' for i in range(6)) + '\n'
        
    fi = 0
    for i, idhdf in enumerate(pltID):

        ax   = plt.Subplot(fig, outer_grid[i])

        Lat[i] = float(Lat[i])

        if Lat[i] >= 50.:     
            ax.set_facecolor('lightcyan')
            
        elif (Lat[i] >= 20.) & (Lat[i] < 50.):
            ax.set_facecolor('lightgreen')
           
        elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
            ax.set_facecolor('mistyrose')
            
        elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
            ax.set_facecolor('cornsilk')
            
        elif (Lat[i] < -50.):
            ax.set_facecolor('lightgrey')
        
        else:
            ax.set_facecolor('lightgrey')
        
        fig.add_subplot(ax)

        
        #----------------------------
        if AvgType == 'Daily':
            Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      = Avg['dailyAvg']

            
        elif AvgType == 'Monthly':
            Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
            Dates        = Avg['dates']
            dateYearFrac = mf.toYearFraction(Avg['dates'])
            AvgData      =  Avg['mnthlyAvg']

        elif AvgType == 'none':
            Dates        = xDates[idhdf]
            dateYearFrac = mf.toYearFraction(Dates)
            AvgData      = yData[idhdf]

        else:
            print 'Error: Define average type: Daily, Monthly, or none'
            exit()

        indsnan = [ni for ni, n in enumerate(AvgData) if np.isnan(n)]



        if len(indsnan) == len(AvgData):

    
        #if (idhdf == 'Mauna Loa') or (idhdf == 'Altzomoni') or (idhdf == 'Jungfraujoch') :

            slope.append(float('nan'))
            slope_e.append(float('nan'))

            slope_p1.append(float('nan'))
            slope_p1_e.append(float('nan'))

            slope_p2.append(float('nan'))
            slope_p2_e.append(float('nan'))

            slope_p3.append(float('nan'))
            slope_p3_e.append(float('nan'))

            amp.append(float('nan'))

            avgD.append(float('nan'))
            stdD.append(float('nan'))

            avgD_p1.append(float('nan'))
            stdD_p1.append(float('nan'))

            avgD_p2.append(float('nan'))
            stdD_p2.append(float('nan'))

            avgD_p3.append(float('nan'))
            stdD_p3.append(float('nan'))

        else:
            #continue



            #--------------------
            #Apply savitzky golay filter (Comment out if not wated)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            if fits:

                weights      = np.ones_like(dateYearFrac)

                #---------------------------------------------------
                # To make a continuous fit
                #---------------------------------------------------
                numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
                dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
                dates2  = np.asarray(dates2)
                dateYearFrac2 = mf.toYearFraction(dates2)
                #---------------------------------------------------

                yyyy = [sngdate.year for sngdate in Dates]
                yyyy = np.asarray(yyyy)
                
                #---------------------------------------------------
                # Calculate the linear regression of all years and estimate Fourier seris
                #---------------------------------------------------
                
                res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                slope.append(res[1]/np.mean(AvgData)*100.0)
                slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                Amp   = np.sum(res[2]**2)
                Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                amp.append(Amp)

                avgD.append(np.mean(AvgData))
                stdD.append(np.std(AvgData))


                #fopen.write(strformat.format(idhdf, float(Lat[i]), np.mean(AvgData), np.std(AvgData), np.max(AvgData), np.min(AvgData), float(Amp)))


                # if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Bremen') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                    #res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                    #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]


            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    

                #ax = plt.Subplot(fig, outer_grid[fi])
                fi +=1
                ax.plot(Dates, AvgData,'k.',markersize=0)
                ax.scatter(Dates,AvgData, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                if yData2:
                    ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                if yData3:
                    ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                if yData4:
                    dtpStd      = np.nanstd(yData4[idhdf])
                    dtpMean     = np.nanmean(yData4[idhdf])

                    ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                    ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  

                # Lat[i] = float(Lat[i])

                # if Lat[i] >= 50.:     
                #     ax.set_facecolor('lightcyan')
                    
                # elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                #     ax.set_facecolor('lightgreen')
                   
                # elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                #     ax.set_facecolor('mistyrose')
                    
                # elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                #     ax.set_facecolor('cornsilk')
                    
                # elif (Lat[i] < -50.):
                #     ax.set_facecolor('lightgrey')
                
                # else:
                #     ax.set_facecolor('lightgrey')
     
            if fits:

                if pltFig:
                    #ax.plot(dates2,f_fourier(dateYearFrac2),label='Intra-annual variability', linewidth=1.0)
                    ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #---------------------------------------------------

                if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Bremen') or  (idhdf == 'Paris') or (idhdf == 'Altzomoni') or (idhdf == 'Tsukuba')   or (idhdf == 'Paramaribo'):

                    yoi  = [[2010, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))



                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund') or (idhdf == 'Izana'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                        roc    =  np.mean(roc)

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])


                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                elif (idhdf == 'Rikubetsu'):

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                  
                        #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Jungfraujoch')  or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze') or (idhdf == 'Kiruna'):

                    yoi  = [[1995, 2002], [2002, 2008], [2009, 2018]]
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                    #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

        if pltFig:
                if (xmin and xmax): ax.set_xlim(xmin, xmax)

                ax.grid(True, color='gray', alpha=0.25)
                ax.tick_params(which='both',labelsize=10)
                #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
                ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
                #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
                #start, end = ax1[i].get_xlim()
                #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
                #ax1[i].set_ylim(bottom=0)

                yearsLc1      = YearLocator(2)
                yearsLc2      = YearLocator(1)
                months        = MonthLocator()
                DateFmt       = DateFormatter('%Y')


                #plt.xticks(rotation=45)
                ax.xaxis.set_major_locator(yearsLc1)
                ax.xaxis.set_minor_locator(yearsLc2)
                #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)
                #ax.set_xlabel('Year')
                #ax1.xaxis.set_tick_params(which='major', pad=15)  
                ax.xaxis.set_tick_params(which='minor',labelbottom='off')
                ax.tick_params(which='both')

                if (ymin and ymax):
                    ax.set_ylim(ymin, ymax)
              
                #ax.set_xticks([])
                #ax.set_yticks([])
                #fig.add_subplot(ax)

    if pltFig:
    
        # all_axes = fig.get_axes()
        # #show only the outside spines
        # for ax in all_axes:
        #     for sp in ax.spines.values():
        #         sp.set_visible(False)
        #         plt.setp(ax.get_xticklabels(), visible=False)
        #     if ax.is_first_row():
        #         ax.spines['top'].set_visible(True)
        #     if ax.is_last_row():
        #         ax.spines['bottom'].set_visible(True)
        #         plt.setp(ax.get_xticklabels(), visible=True)
        #     if ax.is_first_col():
        #         ax.spines['left'].set_visible(True)
        #     if ax.is_last_col():
        #         ax.spines['right'].set_visible(True)

        # if (npanels % 2 == 1): #even
        #     all_axes[-2].spines['bottom'].set_visible(True)
        #     plt.setp(all_axes[-2].get_xticklabels(), visible=True)
        #     all_axes[-2].set_zorder(1)

        
        # all_axes[-1].set_xlabel('Year')
        # all_axes[-2].set_xlabel('Year')
        # all_axes[-3].set_xlabel('Year')

        
        # fig.autofmt_xdate()

        # fig.text(0.03, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        # plt.suptitle(ytypeStr, fontsize=16  )

        # fig.subplots_adjust(left=0.08, bottom=0.05, right=0.975, top=0.95)

        #-----
        all_axes = fig.get_axes()
        #show only the outside spines
        for ax in all_axes:

            for sp in ax.spines.values():
                sp.set_visible(False)
                plt.setp(ax.get_xticklabels(), visible=False)
                plt.setp(ax.get_yticklabels(), visible=False)
                ax.spines['top'].set_visible(True)
                ax.spines['left'].set_visible(True)
                ax.spines['right'].set_visible(True)
                ax.spines['bottom'].set_visible(True)
                #ax.set_tick_params(which='major',labelsize=16)
                # if ax.is_last_row():
                #     ax.spines['bottom'].set_visible(True)
                #     plt.setp(ax.get_xticklabels(), visible=True)

                if ax.is_first_col():
                    ax.spines['left'].set_visible(True)
                    plt.setp(ax.get_yticklabels(), visible=True)
                    ax.tick_params(labelsize = 14)


        for i in range(-1, -4, -1):

            all_axes[i].spines['bottom'].set_visible(True)
            plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
            all_axes[i].set_zorder(1)

            all_axes[i].set_xlabel('Year', fontsize=14)
            all_axes[i].tick_params(labelsize = 14)

        #fig.autofmt_xdate()

        fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
        #plt.suptitle(ytypeStr, fontsize=16  )

        fig.subplots_adjust(left=0.05, bottom=0.075, right=0.975, top=0.975)

        ytypeStr = ytypeStr.replace(' ', '')
        
        if saveFlg: 
            pdfsav.savefig(fig,dpi=200)
            if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
            else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')
        else: 
            plt.show(block=False)


    return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
            avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)


def AnalTS2(npanels, xDates, yData, pltID, Lat, ID, fits=True, AvgType='Daily', smthFlg=True, pltFig=False, saveFlg=False, pdfsav=' ', ytypeStr=' ', unitsStr=' ', ymin=False, ymax=False, yData2=False, yData3=False, yData4=False, period=1):
        
        #--------------------
        #Slope and Amplitude of time series
        #--------------------
        slope       = []   #slope
        slope_e     = []   #slope error

        slope_p1    = []   #slope 1995 - 2002
        slope_p1_e  = []   #slope error 

        slope_p2    = []   #slope 2002 - 2008
        slope_p2_e  = []   #slope error

        slope_p3    = []   #slope  2008 - 2016
        slope_p3_e  = []   #slope error

        amp         = []   #Amplitude 

        avgD        = []
        stdD        = []

        avgD_p1     = []
        stdD_p1     = []

        avgD_p2     = []
        stdD_p2     = []

        avgD_p3     = []
        stdD_p3     = []

        xmin      = dt.date(1985, 1, 1)
        xmax      = dt.date(2019, 12, 31)

        if pltFig:    
            #fig = plt.figure(figsize=(18,13))  
            fig = plt.figure(figsize=(18,13))

            outer_grid = gridspec.GridSpec(npanels, 4, wspace=0.075, hspace=0.075)
            #outer_grid = gridspec.GridSpec(npanels2, 4, wspace=0.11, hspace=0.085) 

        for i, idhdf in enumerate(pltID):
            
            #----------------------------
            if AvgType == 'Daily':
                Avg          = mf.dailyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
                Dates        = Avg['dates']
                dateYearFrac = mf.toYearFraction(Avg['dates'])
                AvgData      = Avg['dailyAvg']

                
            elif AvgType == 'Monthly':
                Avg          = mf.mnthlyAvg(yData[idhdf],xDates[idhdf], dateAxis=1, meanAxis=0)
                Dates        = Avg['dates']
                dateYearFrac = mf.toYearFraction(Avg['dates'])
                AvgData      =  Avg['mnthlyAvg']

            elif AvgType == 'none':
                Dates        = xDates[idhdf]
                dateYearFrac = mf.toYearFraction(Dates)
                AvgData      = yData[idhdf]

            else:
                print 'Error: Define average type: Daily, Monthly, or none'
                exit()

            #--------------------
            #Apply savitzky golay filter (Comment out if not wated)
            #--------------------
            if smthFlg: AvgData = mf.savitzky_golay(AvgData, 7, 3)

            if fits:

                weights      = np.ones_like(dateYearFrac)

                #---------------------------------------------------
                # To make a continuous fit
                #---------------------------------------------------
                numdays = (Dates.max() + dt.timedelta(days=1) - Dates.min()).days
                dates2  = [Dates.min() + dt.timedelta(days=x) for x in range(0, numdays)]
                dates2  = np.asarray(dates2)
                dateYearFrac2 = mf.toYearFraction(dates2)
                #---------------------------------------------------

                yyyy = [sngdate.year for sngdate in Dates]
                yyyy = np.asarray(yyyy)
                
                #---------------------------------------------------
                # Calculate the linear regression of all years and estimate Fourier seris
                #---------------------------------------------------
                
                #res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                #f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 3, half_period=period)
                f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=period)
                perc, intercept_b, slope_b, pfourier_b = res_b

                #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                slope.append(res[1]/np.mean(AvgData)*100.0)
                slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                Amp   = np.sum(res[2]**2)
                Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                amp.append(Amp)

                avgD.append(np.mean(AvgData))
                stdD.append(np.std(AvgData))


                # if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'StD-Maido') or (idhdf == 'Bremen') or (idhdf == 'Paris') or (idhdf == 'Altzomoni'):

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                # else:

                #     res          = mf.fit_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     f_drift, f_fourier, f_driftfourier,  res_std, A= res[3:8]

                #     res_b        = mf.cf_driftfourier(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                #     perc, intercept_b, slope_b, pfourier_b = res_b

                #     #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[0], yyyy[-1])
                #     #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData)*100.0, np.std(slope_b)/np.mean(AvgData)*100.0, yyyy[0], yyyy[-1])

                #     slope.append(res[1]/np.mean(AvgData)*100.0)
                #     slope_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)

                #     Amp   = np.sum(res[2]**2)
                #     Amp   = np.sqrt(Amp)*2.0##/np.mean(f_driftfourier(dateYearFrac)) * 100.0
                #     amp.append(Amp)

                #     avgD.append(np.mean(AvgData))
                #     stdD.append(np.std(AvgData))

                    #res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=1.0)
                    #f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]


            #---------------------------------------------------
            #start plot
            #---------------------------------------------------
            if pltFig:    
                ax = plt.Subplot(fig, outer_grid[i])
                ax.plot(Dates, AvgData,'k.',markersize=0)
                ax.scatter(Dates,AvgData, s=40, facecolor='lightgray', edgecolor='k',alpha=0.85)

                if yData2:
                    ax.plot(xDates[idhdf], yData2[idhdf], color='green')

                if yData3:
                    ax.plot(xDates[idhdf], yData3[idhdf], color='green')

                if yData4:
                    dtpStd      = np.nanstd(yData4[idhdf])
                    dtpMean     = np.nanmean(yData4[idhdf])

                    ax.axhline(y=dtpMean + (2.0*dtpStd), color='red', alpha=0.5)
                    ax.axhline(y=dtpMean - (2.0*dtpStd), color='red', alpha=0.5)  

                Lat[i] = float(Lat[i])

                if Lat[i] >= 50.:     
                    ax.set_facecolor('lightcyan')
                    
                elif (Lat[i] >= 20.) & (Lat[i] < 50.):
                    ax.set_facecolor('lightgreen')
                   
                elif (Lat[i] >= -20.)  & (Lat[i] < 20.):
                    ax.set_facecolor('mistyrose')
                    
                elif (Lat[i] >= -50.)  & (Lat[i] < -20.):
                    ax.set_facecolor('cornsilk')
                    
                elif (Lat[i] < -50.):
                    ax.set_facecolor('lightgrey')
                
                else:
                    ax.set_facecolor('lightgrey')
     
            if fits:

                if pltFig:
                    #ax.plot(dates2,f_fourier(dateYearFrac2),label='Intra-annual variability', linewidth=1.0)
                    ax.plot(dates2,f_driftfourier(dateYearFrac2),label='Fitted Anual Trend + intra-annual variability',linewidth=2.0)
                #---------------------------------------------------

                if (idhdf == 'Eureka') or (idhdf == 'St Petersburg') or (idhdf == 'Boulder') or (idhdf == 'Maido') or (idhdf == 'St Denis') or (idhdf == 'Bremen') or  (idhdf == 'Paris') or (idhdf == 'Altzomoni') or (idhdf == 'Tsukuba')   or (idhdf == 'Paramaribo'):

                    yoi  = [[2010, 2018]]

                    for y in yoi:
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        slope_p1.append(float('nan'))
                        slope_p1_e.append(float('nan'))

                        slope_p2.append(float('nan'))
                        slope_p2_e.append(float('nan'))

                        slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                        slope_p3_e.append(np.std(slope_b)/np.mean(AvgData)*100.0)


                        avgD_p1.append(float('nan'))
                        stdD_p1.append(float('nan'))

                        avgD_p2.append(float('nan'))
                        stdD_p2.append(float('nan'))

                        avgD_p3.append(np.mean(AvgData[indx1]))
                        stdD_p3.append(np.std(AvgData[indx1]))



                elif (idhdf == 'Thule') or (idhdf == 'Lauder') or (idhdf == 'Toronto')  or (idhdf == 'Ny Alesund') or (idhdf == 'Izana'):

                    yoi  = [[2002, 2008], [2009, 2018]]

                    slope_p1.append(float('nan'))
                    slope_p1_e.append(float('nan'))

                    avgD_p1.append(float('nan'))
                    stdD_p1.append(float('nan'))

                    for ii, y in enumerate(yoi):
                        
                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                        roc    =  np.mean(roc)

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])


                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                elif (idhdf == 'Rikubetsu'):

                    yoi  = [[1996, 2002], [2002, 2008]]

                    slope_p3.append(float('nan'))
                    slope_p3_e.append(float('nan'))

                    avgD_p3.append(float('nan'))
                    stdD_p3.append(float('nan'))
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                  
                        #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

                elif (idhdf == 'Jungfraujoch')  or (idhdf == 'Mauna Loa') or (idhdf == 'Wollongong') or (idhdf == 'AHTS') or (idhdf == 'Zugspitze') or (idhdf == 'Kiruna'):

                    yoi  = [[1995, 2002], [2002, 2008], [2009, 2018]]
                
                    for ii, y in enumerate(yoi):

                        indx1  = np.where( (yyyy >= y[0]) &  (yyyy <= y[1]))[0]

                        res          = mf.fit_driftfourier_poly(dateYearFrac, AvgData, weights, 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A,  df_drift = res[3:9]
                        roc    = df_drift(dateYearFrac[indx1[0]:indx1[-1]])
                       
                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Poly)".format(pltID[i], np.mean(roc), np.std(roc), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Poly)".format(pltID[i], np.mean(roc)/np.mean(AvgData[indx1])*100.0, np.std(roc)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        res    = mf.fit_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        f_drift, f_fourier, f_driftfourier,  res_std, A = res[3:8]
                        df_drift     = res[1]
                        roc          = df_drift

                        res_b        = mf.cf_driftfourier(dateYearFrac[indx1], AvgData[indx1], weights[indx1], 2, half_period=period)
                        perc, intercept_b, slope_b, pfourier_b = res_b

                        #print "Rate of Change ({}) = {:.3f} +/- {:.3f} molec/cm2 for years: {} - {} (Linear)".format(pltID[i], res[1], np.std(slope_b), yyyy[indx1[0]], yyyy[indx1[-1]])
                        #print "Rate of Change ({}) = {:.2f} +/- {:.3f}% for years: {} - {} (Linear)".format(pltID[i], res[1]/np.mean(AvgData[indx1])*100.0, np.std(slope_b)/np.mean(AvgData[indx1])*100.0, yyyy[indx1[0]], yyyy[indx1[-1]])

                        if ii == 0:
                            slope_p1.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p1_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p1.append(np.mean(AvgData[indx1]))
                            stdD_p1.append(np.std(AvgData[indx1]))

                        if ii == 1:
                            slope_p2.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p2_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p2.append(np.mean(AvgData[indx1]))
                            stdD_p2.append(np.std(AvgData[indx1]))

                        if ii == 2:
                            slope_p3.append(res[1]/np.mean(AvgData[indx1])*100.0)
                            slope_p3_e.append(np.std(slope_b)/np.mean(AvgData[indx1])*100.0)

                            avgD_p3.append(np.mean(AvgData[indx1]))
                            stdD_p3.append(np.std(AvgData[indx1]))

                    #ax.plot(dailyVals['dates'][indx1],f_drift(dateYearFrac[indx1]),label='Fitted Anual Trend', linewidth=2.0)

            if pltFig:
                ax.set_xlim(xmin, xmax)

                ax.grid(True, color='gray', alpha=0.25)
                ax.tick_params(which='both',labelsize=10)
                #ax.annotate(pltID[i] + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.025, 0.8), xycoords='axes fraction', fontsize=16, ha='left')
                ax.annotate(ID[i].upper() + ' ({0:.2f}$^\circ$)'.format(float(Lat[i])), xy=(0.015, 0.85), xycoords='axes fraction', fontsize=16, ha='left')
                #if i == 0: ax.set_title('{} Total Columns'.format(gasName.upper()),multialignment='center')
                #start, end = ax1[i].get_xlim()
                #ax1[i].xticks.set_ticks(np.arange(min(totClmn[idhdf]), max(totClmn[idhdf])))
                #ax1[i].set_ylim(bottom=0)

                yearsLc1      = YearLocator(2)
                yearsLc2      = YearLocator(1)
                months        = MonthLocator()
                DateFmt       = DateFormatter('%Y')


                #plt.xticks(rotation=45)
                ax.xaxis.set_major_locator(yearsLc1)
                ax.xaxis.set_minor_locator(yearsLc2)
                #ax.tick_params(axis = 'both', which = 'minor', labelsize = 0)
                #ax1.xaxis.set_minor_formatter(DateFormatter('%m'))
                ax.xaxis.set_major_formatter(DateFmt)
                #ax.set_xlabel('Year')
                #ax1.xaxis.set_tick_params(which='major', pad=15)  
                ax.xaxis.set_tick_params(which='minor',labelbottom='off')
                ax.tick_params(which='both')

                if (ymin and ymax):
                    ax.set_ylim(ymin, ymax)
              
                #ax.set_xticks([])
                #ax.set_yticks([])
                fig.add_subplot(ax)

        if pltFig:
        
            #-----
            all_axes = fig.get_axes()
            #show only the outside spines
            for i, ax in enumerate(all_axes):

                for sp in ax.spines.values():
                    sp.set_visible(False)
                    plt.setp(ax.get_xticklabels(), visible=False)
                    plt.setp(ax.get_yticklabels(), visible=False)
                    ax.spines['top'].set_visible(True)
                    ax.spines['left'].set_visible(True)
                    ax.spines['right'].set_visible(True)
                    ax.spines['bottom'].set_visible(True)
                
                    if ax.is_first_col():
                        ax.spines['left'].set_visible(True)
                        plt.setp(ax.get_yticklabels(), visible=True)
                        ax.tick_params(labelsize = 14)

            #fig.autofmt_xdate()

            for i in range(-1, -5, -1):

                all_axes[i].spines['bottom'].set_visible(True)
                plt.setp(all_axes[i].get_xticklabels(), visible=True, rotation=45)
                all_axes[i].set_zorder(1)

                all_axes[i].set_xlabel('Year', fontsize=14)
                all_axes[i].tick_params(labelsize = 14)

            #fig.autofmt_xdate()

            fig.text(0.0075, 0.5, ytypeStr +' ['+unitsStr+']', fontsize=16, va='center', rotation='vertical')
            #plt.suptitle(ytypeStr, fontsize=16  )

            fig.subplots_adjust(left=0.05, bottom=0.075, right=0.99, top=0.99)

            
            if saveFlg: pdfsav.savefig(fig,dpi=200)
            else: 
                plt.show(block=False)

            ytypeStr = ytypeStr.replace(' ', '')


            if fits: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'_wFit.pdf', bbox_inches='tight')
            else: fig.savefig('/data1/projects/ocs/figures/fig/'+ytypeStr+'.pdf', bbox_inches='tight')


        return (slope, slope_e, slope_p1, slope_p1_e, slope_p2, slope_p2_e, slope_p3, slope_p3_e, amp,
                avgD, stdD, avgD_p1, stdD_p1 , avgD_p2, stdD_p2 , avgD_p3, stdD_p3)

def hbarplt3(b1, b2, b3, pltID, b1_label='', b2_label='', b3_label='', subtitle='', saveFlg=False, pdfsav=' '):

    pltID        = np.asarray(pltID)
    pltID        = [p.upper() for p in pltID]

    #-------------------------------------------------
    # Calculate partial columns and weighted VMR
    #-------------------------------------------------

    with open('/data1/projects/ocs/figures/fig/'+'ROC_lowTrop.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b1[0][pi], b1[1][pi], b1[2][pi], b1[3][pi], b1[4][pi], b1[5][pi], b1[6][pi], b1[7][pi]     ) )

    with open('/data1/projects/ocs/figures/fig/'+'ROC_freeTrop.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b2[0][pi], b2[1][pi], b2[2][pi], b2[3][pi], b2[4][pi], b2[5][pi], b2[6][pi], b2[7][pi]     ) )

    with open('/data1/projects/ocs/figures/fig/'+'ROC_Strat.dat','w') as fopen:

        fopen.write('# Site, roc1, roc1_e, roc2, roc2_e, roc3, roc3_e, roc4, roc4_e\n')

        strFormat = '{0:<10s} & {1:>.2f} $\pm$ {2:>.2f} & {3:.2f} $\pm$ {4:.2f} & {5:.2f} $\pm$ {6:.2f} & {7:.2f} $\pm$ {8:.2f}\n'

        for pi, p in enumerate(pltID):

            fopen.write(strFormat.format(p.upper(), b3[0][pi], b3[1][pi], b3[2][pi], b3[3][pi], b3[4][pi], b3[5][pi], b3[6][pi], b3[7][pi]     ) )


    
    #---------------------------------------------------
    # Bar plot: Three different periods ==> Retrieval
    #---------------------------------------------------
    

    ind = np.arange(len(b1[0]))
    
    fig, (ax, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True, figsize=(15, 8.5))
    #fig, (ax2, ax3, ax4) = plt.subplots(1, 3, sharey=True, figsize=(15, 8.5))
    
    ax.barh(ind-0.27, b1[0], 0.27, xerr=b1[1], align='center', color = 'r', ecolor = 'k', label = b1_label)
    ax.barh(ind, b2[0], 0.27, xerr=b2[1], align='center', color = 'b', ecolor = 'k', label = b2_label)
    ax.barh(ind+0.27, b3[0], 0.27, xerr=b3[1], align='center', color = 'g', ecolor = 'k', label = b3_label)  #yerr=slope_TC*0
    ax.xaxis.grid(True)
    ax.yaxis.grid(True)
    ax.set_xlabel('Rate of change [%/y]', fontsize=14)
    ax.set_yticks(ind)
    ax.set_yticklabels(np.transpose(pltID), rotation=0)
    ax.set_title('1996 - 2016', multialignment='center', fontsize=14)
    #ax.set_xlabel('Site')
    ax.set_xlim(-2.0, 2.0)
    ax.axvline(0, color='black', lw=1)
    ax.legend(prop={'size':8}, loc = 1)
    ax.tick_params(labelsize = 14)
    #ax.invert_yaxis()

    ax2.barh(ind-0.27, b1[2], 0.27, xerr=b1[3], align='center', color = 'r', ecolor = 'k', label = b1_label)
    ax2.barh(ind, b2[2], 0.27, xerr=b2[3], align='center', color = 'b', ecolor = 'k', label = b2_label)
    ax2.barh(ind+0.27, b3[2], 0.27, xerr=b3[3], align='center', color = 'g', ecolor = 'k', label = b3_label)  
    ax2.xaxis.grid(True)
    ax2.yaxis.grid(True)
    ax2.set_xlabel('Rate of change [%/y]', fontsize=14)
    ax2.set_title('1996 - 2002', multialignment='center', fontsize=14)
    ax2.set_yticks(ind)
    ax2.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    #ax2.legend(prop={'size':8}, loc = 1)
    ax2.set_xlim(-2.0, 2.0)
    ax2.axvline(0, color='black', lw=1)
    ax2.tick_params(labelsize = 14)
    #ax2.invert_yaxis()

    #ax2.invert_yaxis()
    #ax2.yticks([])

    ax3.barh(ind-0.27, b1[4], 0.27, xerr=b1[5], align='center', color = 'r', ecolor = 'k')
    ax3.barh(ind, b2[4], 0.27, xerr=b2[5], align='center', color = 'b', ecolor = 'k')
    ax3.barh(ind+0.27, b3[4], 0.27, xerr=b3[5], align='center', color = 'g', ecolor = 'k')  
    ax3.xaxis.grid(True)
    ax3.yaxis.grid(True)
    ax3.set_xlabel('Rate of change [%/y]', fontsize=14)
    ax3.set_title('2002 - 2008', multialignment='center', fontsize=14)
    ax3.set_yticks(ind)
    #ax3.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    ax3.set_xlim(-2.0, 2.0)
    ax3.axvline(0, color='black', lw=1)
    ax3.tick_params(labelsize = 14)

    #ax3.invert_yaxis()

    ax4.barh(ind-0.27, b1[6], 0.27, xerr=b1[7], align='center', color = 'r', ecolor = 'k')
    ax4.barh(ind, b2[6], 0.27, xerr=b2[7], align='center', color = 'b', ecolor = 'k')
    ax4.barh(ind+0.27, b3[6], 0.27, xerr=b3[7], align='center', color = 'g', ecolor = 'k')  
    ax4.xaxis.grid(True)
    ax4.yaxis.grid(True)
    ax4.set_xlabel('Rate of change [%/y]', fontsize=14)
    ax4.set_title('2009 - 2018', multialignment='center', fontsize=14)
    ax4.set_yticks(ind)
    #ax4.set_yticklabels(np.transpose(pltID), rotation=0)
    #ax.set_xlabel('Site')
    ax4.set_xlim(-2.0, 2.0)
    ax4.axvline(0, color='black', lw=1)
    ax4.tick_params(labelsize = 14)

    #ax4.invert_yaxis()

    plt.gca().invert_yaxis()
    #fig.tight_layout()
    fig.subplots_adjust(left=0.1, right=0.97, top=0.95, bottom=0.1)


    subtitle = subtitle.replace(' ', '')

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        fig.savefig('/data1/projects/ocs/figures/fig/'+subtitle+'.pdf', bbox_inches='tight')
    else: 
        plt.show(block=False)


def latPlt(resvmrLC, resvmrLC2, resvmrSC, Lat2, saveFlg=False, pdfsav=' '):

    print '\nPlot: Hemispheric Differences:\n'
    Lat2 = np.asarray(Lat2)

    fig, (ax, ax2, ax3, ax4) = plt.subplots(4, 1, sharey=True, sharex=True, figsize=(7, 10))

    #ax.errorbar(Lat, resvmrTC[9], yerr=resvmrTC[10], fmt='o', color='red', ecolor='red', label ='Total')
    ax.errorbar(Lat2, resvmrLC[9], yerr=resvmrLC[10], fmt='o', color='red', ecolor='red', label ='Low Tropospheric')
    ax.errorbar(Lat2, resvmrLC2[9], yerr=resvmrLC2[10], fmt='o', color='blue', ecolor='blue', label= 'Free Tropospheric')
    ax.errorbar(Lat2, resvmrSC[9], yerr=resvmrSC[10], fmt='o', color='green', ecolor='green', label = 'Stratospheric')
    ax.grid(True)
    ax.set_title('Years submitted', multialignment='center', fontsize=14)
    ax.tick_params(labelsize = 14)
    #ax.legend(prop={'size':8}, loc=4)


    ax2.errorbar(Lat2, resvmrLC[11], yerr=resvmrLC[12], fmt='o', color='red', ecolor='red', label ='Low Tropospheric')
    ax2.errorbar(Lat2, resvmrLC2[11], yerr=resvmrLC2[12], fmt='o', color='blue', ecolor='blue', label= 'Free Tropospheric')
    ax2.errorbar(Lat2, resvmrSC[11], yerr=resvmrSC[12], fmt='o', color='green', ecolor='green', label = 'Stratospheric')
    ax2.grid(True)
    ax2.set_title('1998 - 2002', multialignment='center', fontsize=14)
    ax2.tick_params(labelsize = 14)
    
    ax2.legend(ncol=1, prop={'size':10}, loc = 'upper right', bbox_to_anchor=[1.05, 1.2])


    ax3.errorbar(Lat2, resvmrLC[13], yerr=resvmrLC[14], fmt='o', color='red', ecolor='red')
    ax3.errorbar(Lat2, resvmrLC2[13], yerr=resvmrLC2[14], fmt='o', color='blue', ecolor='blue')
    ax3.errorbar(Lat2, resvmrSC[13], yerr=resvmrSC[14], fmt='o', color='green', ecolor='green')
    ax3.grid(True)
    ax3.set_title('2002 - 2008', multialignment='center', fontsize=14)
    ax3.tick_params(labelsize = 14)

    ax4.errorbar(Lat2, resvmrLC[15], yerr=resvmrLC[16], fmt='o', color='red', ecolor='red')
    ax4.errorbar(Lat2, resvmrLC2[15], yerr=resvmrLC2[16], fmt='o', color='blue', ecolor='blue')
    ax4.errorbar(Lat2, resvmrSC[15], yerr=resvmrSC[16], fmt='o', color='green', ecolor='green')
    ax4.grid(True)
    ax4.set_title('2008 - 2016', multialignment='center', fontsize=14)
    ax4.set_xlabel('Latitude', fontsize=14)
    ax4.set_xlim(-90, 90)
    ax2.tick_params(labelsize = 14)

    fig.text(0.015, 0.5, 'Weighted VMR [ppt]', fontsize=16, va='center', rotation='vertical')
    
    fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)

    if saveFlg: 
        pdfsav.savefig(fig,dpi=200)
        fig.savefig('/data1/projects/ocs/figures/fig/Hemispheric_wVMR.pdf', bbox_inches='tight')

    else:  plt.show(block=False)

    


