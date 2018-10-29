#----------------------------------------------------------------------------------------
# Name:
#        MLSHDFClass.py
#
# Purpose:
#       This is a collection of classes and functions used for processing and ploting HDF files
#
#
# Notes:
#
#
# Version History:
#       Created, Sep, 2017  Ivan Ortega (iortega@ucar.edu)
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
import matplotlib.animation as animation
import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
import hdfBaseRetDat
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.cm as mplcm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.gridspec as gridspec
import glob
from cycler import cycler

from pyhdf.SD import SD, SDC
from pyhdf.SD import *
from pyhdf.HDF import *
import coda


from mpl_toolkits.basemap import Basemap

import h5py
import myfunctions as mf

from netCDF4 import Dataset, Group
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA

from pyhdf.HDF import *
from pyhdf.V   import *
from pyhdf.VS  import *



def describevg(refnum, v, vs):
 
    # Describe the vgroup with the given refnum.

    # Open vgroup in read mode.
    vg = v.attach(refnum)
    print "----------------"
    print "name:", vg._name, "class:",vg._class, "tag,ref:", vg._tag, vg._refnum

    # Show the number of members of each main object type.
    print "members: ", vg._nmembers,
    print "datasets:", vg.nrefs(HC.DFTAG_NDG),
    print "vdatas:  ", vg.nrefs(HC.DFTAG_VH),
    print "vgroups: ", vg.nrefs(HC.DFTAG_VG)

    # Read the contents of the vgroup.
    members = vg.tagrefs()

    # Display info about each member.
    index = -1
    for tag, ref in members:
        index += 1
        print "member index", index
        # Vdata tag
        if tag == HC.DFTAG_VH:
            vd = vs.attach(ref)
            nrecs, intmode, fields, size, name = vd.inquire()
            print "  vdata:",name, "tag,ref:",tag, ref
            print "    fields:",fields
            print "    nrecs:",nrecs
            vd.detach()

        # SDS tag
        elif tag == HC.DFTAG_NDG:
            sds = sd.select(sd.reftoindex(ref))
            name, rank, dims, type, nattrs = sds.info()
            print "  dataset:",name, "tag,ref:", tag, ref
            print "    dims:",dims
            print "    type:",type
            sds.endaccess()

        # VS tag
        elif tag == HC.DFTAG_VG:
            vg0 = v.attach(ref)
            print "  vgroup:", vg0._name, "tag,ref:", tag, ref
            vg0.detach()

        # Unhandled tag
        else:
            print "unhandled tag,ref",tag,ref

    # Close vgroup
    vg.detach()

def HDFread(refnum, v, vs, section='', group=''):
    """
    Extract the data for non-scientific data in V mode of hdf file
    """

    vg = v.attach(refnum)
    print "----------------"
    print "name:", vg._name, "class:",vg._class, "tag,ref:", vg._tag, vg._refnum

    if vg._name in section:

        # Show the number of members of each main object type.
        print "members: ", vg._nmembers,
        print "datasets:", vg.nrefs(HC.DFTAG_NDG),
        print "vdatas:  ", vg.nrefs(HC.DFTAG_VH),
        print "vgroups: ", vg.nrefs(HC.DFTAG_VG)

        # Read the contents of the vgroup.
        members = vg.tagrefs()

        # Display info about each member.
        index = -1
        for tag, ref in members:
            index += 1
            print "member index", index
            # Vdata tag
            if tag == HC.DFTAG_VH:
                vd = vs.attach(ref)
                nrecs, intmode, fields, size, name = vd.inquire()

                if vd._name in group:
                
                    print "  vdata:",vd._name, "tag,ref:",tag, ref
                    print "  fields:",fields
                    print "  nrecs:",nrecs

                    data = np.asarray(vd[:])

                vd.detach()

            # SDS tag
            elif tag == HC.DFTAG_NDG:
                sds = sd.select(sd.reftoindex(ref))
                name, rank, dims, type, nattrs = sds.info()
                print "  dataset:",name, "tag,ref:", tag, ref
                print "    dims:",dims
                print "    type:",type
                sds.endaccess()

            # VS tag
            elif tag == HC.DFTAG_VG:
                vg0 = v.attach(ref)
                print "  vgroup:", vg0._name, "tag,ref:", tag, ref
                vg0.detach()

            # Unhandled tag
            else:
                print "unhandled tag,ref",tag,ref

        return nrecs, intmode, fields, size, name, data

    # Close vgroup
    vg.detach()

    




    

                                                #----------------#
                                                # Define classes #
                                                #----------------#

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
#------------------------------------------------------------------------------------------------------------------------------                     
class ReadSAGE3():


    def __init__(self, dataDir):

        #---------------------------------
        #
        #---------------------------------
        if not( dataDir.endswith('/') ):
                dataDir = dataDir + '/'

        self.dataDir    = dataDir

        self.HDF          = {}
        self.ReadSAGE3Flg = False

        #-------------------------
        # Create Directory w files
        #-------------------------
        DirFiles = glob.glob(dataDir + 'g3b.ssp*')
        DirFiles.sort()

        self.dirLst     = DirFiles

        self.showinfodFlg  = False
        self.readFlg       = False

    
    def showHDFinfo(self):

        if len(self.dirLst) > 1: drs = self.dirLst[0]
        else: drs = self.dirLst

        hdf = HDF(drs, HC.READ)
        sd  = SD(drs)

        vs  = hdf.vstart()
        v   = hdf.vgstart()

        # Scan all vgroups in the file.
        ref = -1
        while 1:
            try:
                ref = v.getid(ref)
                describevg(ref, v, vs)
           
            except HDF4Error,msg:    # no more vgroup
                    break

        self.showinfodFlg  = True

        #vs.detach()

    #--------------------------------------
    # Read all HDF files from a list
    #--------------------------------------
    def hdfget(self, id =False, section=False, group=False):
  
        #for i, drs in enumerate(self.dirLst[0:500]):
        for i, drs in enumerate(self.dirLst[:]):

            #print drs

            hdf = HDF(drs, HC.READ)
            sd  = SD(drs)

            vs  = hdf.vstart()
            v   = hdf.vgstart()

            # Scan all vgroups in the file.
            ref = -1
            while 1:
                try:
                    ref = v.getid(ref)

                    vg = v.attach(ref)
                
                    if vg._name == section:

                        # Read the contents of the vgroup.
                        members = vg.tagrefs()

                        # Display info about each member.
                        index = -1
                        for tag, ref2 in members:
                            index += 1
                            #print "member index", index
                            # Vdata tag
                            if tag == HC.DFTAG_VH:
                                vd = vs.attach(ref2)
                                nrecs, intmode, fields, size, name = vd.inquire()

                                if vd._name == group:

                                    self.HDF.setdefault(id+'_fields',[]).append(fields)
                                    self.HDF.setdefault(id+'_nrecs',[]).append(nrecs)
                                    self.HDF.setdefault(id+'_data',[]).append(np.asarray(vd[:]))

                                vd.detach()

                            # SDS tag
                            elif tag == HC.DFTAG_NDG:
                                sds = sd.select(sd.reftoindex(ref2))
                                name, rank, dims, type, nattrs = sds.info()

                                sds.endaccess()

                            elif tag == HC.DFTAG_VG:
                                vg0 = v.attach(ref2)
                                vg0.detach()

                            else:
                                print "unhandled tag,ref",tag,ref2

                except HDF4Error,msg:    # no more vgroup
                    break
                #describevg(ref, v, vs)
                #HDFread(ref, v, vs, section='Section 5.0 - Altitude-based Data', group='Tropopause Data')

    def ReadVariables(self):

        if not self.showinfodFlg: self.showHDFinfo()

        print "----------------" 
        print '\nTotal Number of files in SAGE folder: {}\n'.format(len(self.dirLst ))

        #------------------------------------------------------------
        # ID variables defined by the user
        #------------------------------------------------------------
        id1          = 'altitude'
        id2          = 'o3'
        id3          = 'track'
        id4          = 'start'
        id5          = 'end'

        #------------------------------------------------------------
        # Below are sections and groups (These are case sensitive and must be present in the HDF files)
        #------------------------------------------------------------
        id1InHDF     = ['Section 5.0 - Altitude-based Data',                 'Tropopause Data' ]
        id2InHDF     = ['Section 5.2 - Composite Ozone profiles',            'Ozone profiles' ]
        id3InHDF     = ['Section 4.3 - Ground Track Data Over This Event',   'Section 4.3 - Ground Track Data Over This Event' ]
        id4InHDF     = ['Section 4.1 - Science Data Start Information',      'Section 4.1 - Science Data Start Information' ]
        id5InHDF     = ['Section 4.2 - Science Data End Information',        'Section 4.2 - Science Data End Information' ]

        #------------------------------------------------------------
        # Call hdfget
        #------------------------------------------------------------
        self.hdfget(id=id1, section=id1InHDF[0], group=id1InHDF[1])
        self.hdfget(id=id2, section=id2InHDF[0], group=id2InHDF[1])        
        self.hdfget(id=id3, section=id3InHDF[0], group=id3InHDF[1])        
        self.hdfget(id=id4, section=id4InHDF[0], group=id4InHDF[1])        
        self.hdfget(id=id5, section=id5InHDF[0], group=id5InHDF[1])

        print "----------------" 
        print '\nTrack fields: {}'.format(self.HDF['track_fields'][0])
        print '\no3 fields: {}'.format(self.HDF['o3_fields'][0])
        print '\naltitude fields: {}'.format(self.HDF['altitude_fields'][0])
        print '\nstart fields: {}'.format(self.HDF['start_fields'][0])
        print '\nend fields: {}'.format(self.HDF['start_fields'][0])
        print "----------------"
    

        #------------------------------------------------------------
        #
        #------------------------------------------------------------        
        startDate    = np.asarray(self.HDF['start_data'], dtype=np.int32)[:, 0, 0]
        endDate      = np.asarray(self.HDF['end_data'], dtype=np.int32)[:, 0, 0]

        startTime    = np.asarray(self.HDF['start_data'], dtype=np.int32)[:, 0, 1]
        endTime      = np.asarray(self.HDF['end_data'], dtype=np.int32)[:, 0, 1]

        Dates         = [np.mean([startDate[i], endDate[i]]) for i in range(len(startDate))]
        Times         = [np.mean([startTime[i], endTime[i]]) for i in range(len(startDate))]

        Dates         = np.asarray(Dates, dtype='i')
        Times         = np.asarray(Times, dtype='i')

        #Time2        = ["{:06d}".format(t) for t in Time  ]
        #Dates2       = ["{:08d}".format(t) for t in Dates  ]

        #Altitude = np.asarray(self.HDF['altitude_data'])[:, :, 1]
        #Dates        = np.asarray(self.HDF['track_data'], dtype=np.int32)[:, 1:, 0]
        #Time         = np.asarray(self.HDF['track_data'], dtype=np.int32)[:, 1:, 1]
        

        #Dates        = np.asarray(np.mean(np.asarray(Dates) , dtype='i', axis=1))
        #Time         = np.asarray(np.nanmean(np.asarray(Time) , dtype='i', axis=1))

        #------------------------------------------------------------
        #
        #------------------------------------------------------------  
        Lat          = np.asarray(self.HDF['track_data'])[:, 1:, 2]
        Lon          = np.asarray(self.HDF['track_data'])[:, 1:, 3]
        Lat          = np.asarray(np.mean(Lat, axis=1))
        Lon          = np.asarray(np.mean(Lon, axis=1))

        #------------------------------------------------------------
        #
        #------------------------------------------------------------  

        Times2        = ["{:06d}".format(t) for t in Times  ]
        Dates2       = ["{:08d}".format(t) for t in Dates  ]

        #------------------------------------------------------------
        #
        #------------------------------------------------------------  

        dateTime    = []

        index   = -1
        indxBad = []

        indsT   = []

        for i, d in enumerate(Dates2):
            index += 1
            
            try:
                dateTime.append(dt.datetime(int(str(d)[0:4]), int(str(d)[4:6]), int(str(d)[6:8]), int(str(Times2[i])[0:2]), int(str(Times2[i])[2:4])) )
            except:
                indsT.append(index)
                dateTime.append(dt.datetime(int(str(Dates2[index-1])[0:4]), int(str(Dates2[index-1])[4:6]), int(str(Dates2[index-1])[6:8]), int(str(Times2[index-1])[0:2]), 0) )

        indsT   = np.asarray(indsT)
        print ('Total number observations found with non sense dates/times = {}'.format(len(indsT)))
        indxBad = np.union1d(indsT, indxBad)

        #------------------------------------------------------------
        #
        #------------------------------------------------------------  
        O3_prf_conc  = np.asarray(self.HDF['o3_data'])[:, :, 0] # --> Concentration Profile
        O3_prf_cold  = np.asarray(self.HDF['o3_data'])[:, :, 2] # --> Column Density Profile
        O3_QAFlg     = np.asarray(self.HDF['o3_data'])[:, :, 4] # --> Quality
        altitude     = np.asarray(self.HDF['altitude_data'])[:, :, 1]

        #----------------------------------------
        #
        #----------------------------------------
        O3_Pcol  = []

        for i in range(len(dateTime)):

            inds = np.where( (altitude[i, :] >= 15.0) &  (altitude[i, :] <= 30.0) )
            O3_Pcol.append(np.sum(O3_prf_cold[i, inds]) )

        print 'Number of Profiles: {}'.format(len(dateTime))
        print 'Number of Grid Points: {}'.format(altitude.shape[1])

        #----------------------------------------
        #
        #----------------------------------------
        #badPrf = np.asarray(O3_QAFlg) != 0
        #indsT  = np.where( np.sum(badPrf, axis=1) != 0 )[0]
        #print ('Total number observations found with bad Flag = {}'.format(len(indsT)))
        #indxBad = np.union1d(indsT, indxBad)

        #rprf_neg = np.asarray(O3_prf_conc) <= 0
        #indsT = np.where( np.sum(rprf_neg,axis=1) > 0 )[0]
        #print ('Total number observations found with negative partial column = {}'.format(len(indsT)))
        
        #indxBad = np.union1d(indsT, indxBad)

        #----------------------------------------
        # Assigning NAN to non-sense values in profiles
        #----------------------------------------
        rprf_pos = np.asarray(O3_prf_conc) >= 1e15
        O3_prf_conc[rprf_pos] = np.nan

        rprf_neg = np.asarray(O3_prf_conc) < 0
        O3_prf_conc[rprf_neg] = np.nan
        
        #indsT = np.where( np.sum(rprf_neg,axis=1) > 1e13 )[0]
        #print ('Total number observations found with negative partial column = {}'.format(len(indsT)))
        
        #indxBad = np.union1d(indsT, indxBad)
        #---------------------------------------------
        # Find total column amount < minTC and > maxTC
        #---------------------------------------------
        
        # indsT1 = np.where(np.asarray(O3_Pcol) < 0)[0]
        # indsT2 = np.where(np.asarray(O3_Pcol) > 1e30)[0]
        # indsT  = np.union1d(indsT1,indsT2)
        # print "Total number of observations found with partial column < minTotalColumn = {}".format(len(indsT1))
        # print "Total number of observations found with partial column > maxTotalColumn = {}".format(len(indsT2))

        # indxBad = np.union1d(indsT, indxBad)  

        #----------------------------------------
        #
        #----------------------------------------
        self.dateTime         = np.delete(dateTime,indxBad)
        self.Lat              = np.delete(Lat,indxBad)
        self.Lon              = np.delete(Lon,indxBad)
        self.O3_prf_conc      = np.delete(O3_prf_conc,indxBad, axis=0)
        self.O3_prf_cold      = np.delete(O3_prf_cold,indxBad, axis=0)
        self.O3_Pcol          = np.delete(O3_Pcol,indxBad)
        self.altitude         = np.delete(altitude,indxBad, axis=0)
        dates                 = [dt.date(d.year, d.month, d.day) for d in self.dateTime]
        self.date             = np.asarray(dates)

        #---------------------------------------------
        # 
        #---------------------------------------------
        print 'First Date: {}'.format(dates[0])
        print 'Last Date:  {}'.format(dates[-1])

        self.readFlg = True



class SAGE3Cls(ReadSAGE3):

    def __init__(self,dataDir, saveFlg=False, outFname=''):

        print "----------------" 
        print " SAGE III/ISS"
        print "----------------" 
        #------------------------------------------------------------
        # If outFname is specified, plots will be saved to this file,
        # otherwise plots will be displayed to screen
        #------------------------------------------------------------
        if saveFlg:  self.pdfsav = PdfPages(outFname)
        else:        self.pdfsav = False

        self.dataDir = dataDir

        #---------------
        # ReadOutputData
        #---------------
        ReadSAGE3.__init__(self, dataDir)   
       
        
    def closeFig(self):
        self.pdfsav.close()



    def pltSAGE3(self):

        #------------------------------------------------------------
        # Show General Information Stored in one HDF file
        #------------------------------------------------------------
        if not self.showinfodFlg: self.showHDFinfo()

        #------------------------------------------------------------
        # Get HDF Specific Variables 
        #------------------------------------------------------------
        if not self.readFlg: self.ReadVariables()

        #----------------------------------------
        #
        #----------------------------------------
        clmap        = 'jet'
        cm           = plt.get_cmap(clmap)  

        #fig, ax = plt.subplots()

        fig,ax  = plt.subplots(figsize=(8,5))


        m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
                    llcrnrlon=-180,urcrnrlon=180,resolution='c')

        m.drawcoastlines()
        m.drawcountries()

        m.drawparallels([-90,-60,-30,0,30,60,90],labels=[1,0,0,0],fontsize=10)
        m.drawmeridians([-180,-120,-60,0,60,120,180],labels=[0,0,0,1],fontsize=10)

        m.drawmapboundary(fill_color='white')
        #m.fillcontinents(color='gray')

        #x, y = m(*np.meshgrid(Lon, Lat))
        x, y = np.meshgrid(self.Lon, self.Lat)
        x, y = m(x, y)

        #plt.scatter(xx,yy,c=zz,cmap=cm.Reds)
        #ax.scatter(x,y,c=O3_Pcol, cmap=mplcm.jet)
        #ax.scatter(x,y, 'k.')
        #z = z.reshape(xx.shape)
        #-- contour levels
        #clevs = np.arange(np.min(O3_Pcol),np.max(O3_Pcol), )

        #cs = m.contourf(x,y,O3_Pcol, cmap=mplcm.jet,extend='both') #,alpha=0.5)
        #cs.axis='tight'

        ax.plot(x,y,'.r', markersize=6, alpha=0.15)

        #plt.title('Track Fields')
        plt.title('Track - SAGE III/ISS - {} - {} '.format(self.date[0], self.date[-1]), fontsize=14)


        if self.pdfsav: self.pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)  

        #----------------------------------------
        #
        #----------------------------------------
        fig, ax  = plt.subplots(figsize=(6,6))

        cm             = plt.get_cmap('jet')
        cNorm          = colors.Normalize( vmin=np.nanmin(self.Lat), vmax=np.nanmax(self.Lon) )
        scalarMap      = mplcm.ScalarMappable( norm=cNorm, cmap=cm )        
        scalarMap.set_array(self.Lat)

        ax.set_prop_cycle( cycler('color', [scalarMap.to_rgba(x) for x in self.Lat] ) )

        for i in range(len(self.dateTime)):
        #for i in inds:
            ax.plot(self.O3_prf_conc[i,:]/1e13, self.altitude[i,:],  linewidth = 1.0)
            #ax.plot(O3_prf_conc[i,:]/1e13,altitude[i,:], 'k.')
        
        #ax.fill_betweenx(alt,Prfmean-prfSTD,Prfmean+prfSTD,alpha=0.5,color='0.75')  
        ax.set_title('O$_3$ - SAGE III/ISS - {} - {}'.format(self.date[0], self.date[-1]), fontsize=14)
        ax.set_ylabel('Altitude [km]', fontsize=14)
        ax.set_xlabel('Concentration [x10$^{13}$ molec/cm$^{3}$]', fontsize=14) 
        #ax.set_xlim(xmin=1e10, xmax=1e13)    
        ax.grid(True,which='both')
        ax.tick_params(labelsize=14)

        cbar = fig.colorbar(scalarMap,orientation='vertical')
        cbar.set_label('Latitude', fontsize=14)

        fig.subplots_adjust(left=0.12, bottom=0.1,top=0.94, right=0.95)

        if self.pdfsav: self.pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)  

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()

        #print altitude[0, :].shape, O3_prf_conc.shape, Lat.shape

        #----------------------------------------
        #
        #----------------------------------------
       

        fig,ax1  = plt.subplots(figsize=(8,5))

        levels        = np.arange(0, 0.6, 0.05)

        cax          = ax1.contourf(self.Lat,self.altitude[0, :],np.transpose(self.O3_prf_conc)/1e13, levels, cmap=mplcm.jet)
        #cax          = ax1.contourf(altitude[0, :], Lat, O3_prf_conc, cmap=mplcm.jet)

        divider      = make_axes_locatable(ax1)
        cb           = divider.append_axes("right",size="3.5%",pad=0.05)
        
        cbar         = plt.colorbar(cax,cax=cb)
        #cbar.ax.tick_params(labelsize=12)
        #cbar.set_label('VMR [x10$^3$ ppm]', fontsize=12)#, fontsize=14)
        cbar.set_label('Concentration [x10$^{13}$ molec/cm$^{3}$]', fontsize=14) 

        #ax1.set_title('Ozone Partial Column')
        ax1.set_title('O$_3$ - SAGE III/ISS - {} - {} '.format(self.date[0], self.date[-1]), fontsize=14)
        ax1.set_ylabel('Altitude [km]', fontsize=14)
        ax1.set_xlabel('Latitude', fontsize=14)    
        ax1.set_xlim((-55,60))
        ax1.grid(True,which='both')
        ax1.tick_params(labelsize=14)

        fig.subplots_adjust(left=0.1, bottom=0.1,top=0.94, right=0.92)
         

        if self.pdfsav: self.pdfsav.savefig(fig,dpi=200)
        else:           plt.show(block=False)  


    

