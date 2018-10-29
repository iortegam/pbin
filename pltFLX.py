#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltFTSWRF.py
#
# Purpose:
#       1) Read and plot netCDF files provided by Gabi Pfister during FRAPPE (WRF-CHEM)
#       2) Integrate Model outputs and FTS
#
# Notes:
#   
#
# Version History:
#       Created, Sep, 2016  Ivan Ortega (iortega@ucar.edu)

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

from scipy import interpolate

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

from netCDF4 import Dataset
import dataModelOutClass as dm
import dataOutClass as dc
from collections                     import OrderedDict
from scipy import linspace, polyval, polyfit, sqrt, stats, randn


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
    #                     Global Inputs independently of FTS and/or Models
    #-----------------------------------------------------------------------------------------
    loc               = 'fl0'
    pCols             = [1.6, 8.0]       #--ALTITUDE TO CALCULATE PARTIAL COLUMNS AND WEIGHTED VMR
    
    pltAnalysis       = False
    saveFlg           = False 
    
    pltFile           =  '/data/iortega/results/'+loc.lower()+'/fig/FTS-FLEXPART_'+loc.upper()+'.pdf'

    #-----------------------------------------------------------------------------------------
    #                 Initializations for WRF-Chem during FRAPPE (Jul - Aug 2014)
    #-----------------------------------------------------------------------------------------
    #dataDirFLX        = '/data1/ancillary_data/'+loc.lower()+'/FRAPPE/flexpart2/'
    dataDirFLX        = '/ya4/Campaign/'+loc.upper()+'/FRAPPE/flexpart2/'
    
    ReadFLX           = True
    pltFLX            = False
    saveFLXFlg        = False
    
    pltFLXFile        = '/data/iortega/results/'+loc.lower()+'/fig/FLEXPART__FRAPPE.pdf'

    #----------------------
    # Date range to process
    # #----------------------
    # iyear              = 2014
    # imnth              = 7
    # iday               = 22
    # fyear              = 2014
    # fmnth              = 7
    # fday               = 22

    iyear              = 2014
    imnth              = 7
    iday               = 14
    fyear              = 2014
    fmnth              = 8
    fday               = 20

                                    #----------------------------#
                                    #        --- START ---       #
                                    #----------------------------#

    if loc == 'fl0':
        sLat              = 40.035#    40.4             #--LATITUDE OF BOULDER
        sLon              = -105.245          #--LONGITUDE OF BOULDER
    else:
        print 'Need a location (currently only fl0)'
        sys.exit()


    #-----------------------------------------------------------------------------------------
    #FLEXPART
    #-----------------------------------------------------------------------------------------

    if ReadFLX:
        DataFLX = dm.FLXClass(dataDirFLX, iyear,imnth,iday,fyear,fmnth,fday,  outFname= pltFLXFile, saveFlg= saveFLXFlg)
        DataFLX.ReadOutputFLX(pltFLX)
        #if pltFLX:
        #    DataFLX.PltFLX()

        with open(dataDirFLX + 'FLX_SCF_12h.txt','w') as fopen:
            hdr = '#Results of the Fractional Contribution calculated with pltFLX.py\n'
            fopen.write(hdr)
            hdr = 'Date,  SC_NE,  SC_NW,  SC_SE,  SC_SW\n'
            fopen.write(hdr)


            strformat = '{0:<12}, {1:.4f}, {2:.4f}, {3:<.4f}, {4:<.4f}\n'
            for i,indDay in enumerate(DataFLX.starttime):
                daystr = "{0:04d}{1:02d}{2:02d}{3:02d}".format(indDay.year,indDay.month,indDay.day, indDay.hour)
                temp   = [daystr, DataFLX.FLX['SC']['NE'][i], DataFLX.FLX['SC']['NW'][i], DataFLX.FLX['SC']['SE'][i], DataFLX.FLX['SC']['SW'][i]]
                fopen.write(strformat.format(*temp))

        fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)

        ax.plot(DataFLX.starttime, DataFLX.FLX['SC']['NE'], color='k')
        ax.scatter(DataFLX.starttime, DataFLX.FLX['SC']['NE'], facecolors='white', s=60, color='k')
        
        ax.grid(True)
        ax.set_ylabel('SC',fontsize=14)
        ax.tick_params(labelsize=14)
        ax.set_xlabel('Date', fontsize=14)
    


    if saveFlg:
        DataFLX.pdfsav.close()
    else:
        plt.show(block=False)
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()


 
if __name__ == "__main__":
    main()