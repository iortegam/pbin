#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#         printCAM.py
#
# Purpose:
#         Read and Create a .dat file for input gas in CAM-Chem using Smoothing AK (see ianputs below)
#
# Notes:
#   
#
# Version History:
#       Created, July, 2018  Ivan Ortega (iortega@ucar.edu)

#----------------------------------------------------------------------------------------

                                    #-------------------------#
                                    # Import Standard modules #
                                    #-------------------------#

from itertools import izip
import numpy as np

from scipy import interpolate
import datetime as dt

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
from cycler import cycler


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
    #                 Initializations for CAM-CHEM
    #-----------------------------------------------------------------------------------------
    #gasName          = ['h2co']#, 'co', 'c2h6']
    #gasName             = ['co',           'c2h2',         'c2h6',         'ch4',          'h2co',         'hcooh']   
    gasName            = ['co',           'c2h2',         'c2h6',         'ch4',         'nh3',          'hcn',         'h2co',         'hcooh']            

    #dataDirCAM          = '/data1/ancillary_data/fl0/'
    #fileCAM             = 'CAM_chem_fmerra_FSDSSOA_2deg_2000_2014_extra_Boulder.nc'     

    dataDirCAM        = '/net/modeling1/data16a/buchholz/CAM_chem_output/CAM_chem_fmerra2_FCSD_1deg_Boulder/'
    fileCAM           = 'CAM_chem_fmerra2_FCSD_1deg_Boulder_2009_2017.nc'

    saveFlg             = True  

    sLat                = 40.4             #--LATITUDE OF BOULDER
    sLon                = -105.24          #--LONGITUDE OF BOULDER

    pCols              = [1.6, 8.0]

    #----------------------
    # Date range to process
    #----------------------
    iyear              = 2010
    imnth              = 1
    iday               = 1
    fyear              = 2017
    fmnth              = 12
    fday               = 31

    pltFile           =  '/data/iortega/pbin/tropGases/fig/CAM-Chem.pdf'


                                    #----------------------------#
                                    #                            #
                                    #        --- START ---       #
                                    #                            #
                                    #----------------------------#

    #-------------------------------------------------
    #    -- Read CAM-CHEM --
    #-------------------------------------------------
    DataCAM = dm.CAMClass(dataDirCAM, fileCAM,  outFname= pltFile, saveFlg= saveFlg)
    DataCAM.ReadOutputCAM(gasName, pCols, sLat, sLon)  
    DataCAM.PltCAM()  

    #---------------------------------
    #     -- Read FTS --
    #---------------------------------


      
 
if __name__ == "__main__":
    main()