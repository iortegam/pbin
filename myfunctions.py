import datetime as dt
from datetime import datetime, timedelta
import time
import math
import sys
import numpy as np
from numpy import fromstring, vectorize, ndarray, array, genfromtxt
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
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
#from scipy.interpolate import interp1d

#MODIFIED (IVAN)
#import matplotlib
# Force matplotlib to not use any Xwindows backend.
#matplotlib.use('Agg')
#



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
import matplotlib.colorbar as colorbar
import glob
from math import factorial

from math import radians, cos, sin, asin

from scipy.odr import Model, Data, RealData, ODR
from scipy.stats import linregress
import numpy as np

from math import acos, asin, atan2, hypot
from math import degrees, pi as PI, tan
from collections import deque,Counter
from itertools import islice

from scipy import interpolate
from pylab import setp





                                            #------------------#
                                            # Define functions #
                                            #------------------#

def nearestDate(daysList, year, month, day=1):
    ''' Finds the nearest date from a list of days based on a given year, month, and day'''
    testDate = dt.date(year, month, day)
    return min( daysList, key=lambda x:abs(x-testDate) )

def nearestind(val,array):
    ''' Returns the index in array for closest value to val'''
    return (np.abs(array-val)).argmin()

def daysList(timeList):
    ''' Finds a unique set of days within a listof dates '''
    dateList = [dt.date(x.year,x.month,x.day) for x in timeList]
    dateList = np.unique(dateList)
    
    return np.sort(dateList)

def sortDict(DataDict,keyval,excld=[]):
    ''' Sort all values of dictionary based on values of one key'''
    base = DataDict[keyval]
    for k in DataDict:
        if k not in excld: DataDict[k] = [y for (x,y) in sorted(zip(base,DataDict[k]))]
    return DataDict

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
    
def tryopen(fname):
    ''' Try to open a file and read contents'''
    try:
        with open(fname, 'r' ) as fopen:
            return fopen.readlines()
    except IOError as errmsg:
        print errmsg
        return False

def convrtD(rhs):
    '''This identifies numbers specified in scientific notation with d 
       for double precision and converts them to floats'''
    if 'd' in rhs.lower():                               # Test and handle strings containing D (ie 1.0D0)
        mant,trsh,expo = rhs.lower().partition('d')
        try:
            rhs = float(mant)*10**int(expo)
        except ValueError:
            pass
    else:                                                # Handle strings => number
        try:
            rhs = float(rhs)
        except ValueError:
            pass        

    return rhs

def dbFilterUL(fltrParam,lwrVal=False,uprVal=False):
    ''' Filter the output data based on upper an lower value and
        return the corresponding indicies'''
    inds = []
    #-------------------------------------------
    # At least one bound must be set. Check this
    #-------------------------------------------
    #if (not lwrVal) or (not uprVal):
    #    print "Must specify at least one bound in dbFilterUL"
    #    return

    #----------------------------------------
    # Case where lower value is not specified
    #----------------------------------------
    if not lwrVal:
        for ind,val in enumerate(fltrParam):
            if val <= uprVal: inds.append(ind)

    #----------------------------------------
    # Case where upper value is not specified
    #----------------------------------------
    elif not uprVal:
        for ind,val in enumerate(fltrParam):
            if val <= lwrVal: inds.append(ind)

    #---------------------------------------------
    # Case where lower and upper val are specified
    #---------------------------------------------
    else:
        for ind,val in enumerate(fltrParam):
            if ( val >= lwrVal and val <= uprVal ): inds.append(ind)

    return inds


def dbFilterTF(fltrParam,fltrVal):
    ''' Filter the output data based on True False value
        and return the corresponding indicies'''
    inds = []

    #------------------------------
    # Find cases that match fltrVal
    #------------------------------
    for ind,val in enumerate(fltrParam):
        if val == fltrVal: inds.append(ind)

    return inds


def getintrsct(vals1,vals2):
    ''' Return the indicies of vals1 and vals2 for the intersection '''

    intrsctVals = np.intersect1d(vals1,vals2)
    inds1       = np.nonzero( np.in1d( vals1, intrsctVals ) )[0]
    inds2       = np.nonzero( np.in1d( vals2, intrsctVals ) )[0]

    return (inds1,inds2)


def bias(xData,yData):
    ''' Cacluates the mean value of the residuals. 
    Where xData is observation and yData is the model'''

    #--------------------------------
    # Make sure arrays are numpy type
    #--------------------------------
    yData = np.asarray(yData)
    xData = np.asarray(xData)

    #-------------------------------
    # Make sure arrays are same size
    #-------------------------------
    if ( yData.shape[0] != xData.shape[0] ): 
        print 'Data sets must have same size in corrcoef: xData = {}, yData = {}'.format(xData.shape[0],yData.shape[0])
        sys.exit() 

    biasCalc = np.nansum( yData - xData ) / len(yData) 

    return biasCalc    

def rmse(xData, yData):
    '''Calcuates the RMSE and Persons correlation coefficient, 
       xData = observed data, yData = modeled data'''

    #--------------------------------
    # Make sure arrays are numpy type
    #--------------------------------
    yData = np.asarray(yData)
    xData = np.asarray(xData)

    #-------------------------------
    # Make sure arrays are same size
    #-------------------------------
    if ( yData.shape[0] != xData.shape[0] ): 
        print 'Data sets must have same size in corrcoef: xData = {}, yData = {}'.format(xData.shape[0],yData.shape[0])
        sys.exit() 

    #---------------------------------
    # Find the residual sum of squares 
    #---------------------------------
    SS_res = np.nansum( (yData - xData)**2 )

    #------------------------------
    # Root Mean Square Error (RMSE)
    #------------------------------
    rmse = np.sqrt( SS_res / len(yData) )

    return rmse

def yearList(self):
    ''' Gives a list of unique years within DateRange '''
    years = [ singDate.year for singDate in self.dateList]               # Find years for all date entries
    years = list(set(years))                                             # Determine all unique years
    years.sort()
    return years

def daysInYear(dateList, year):
    ''' Returns an ordered list of days from DateRange within a specified year '''
    if isinstance(year,int):
        newyears = [inYear for inYear in dateList if inYear.year == year]
        return newyears
    else:
        print 'Error!! Year must be type int for daysInYear'
        return False


def toYearFraction(dates):
    ''' Convert datetime to year and fraction of year'''

    #def sinceEpoch(date): # returns seconds since epoch
        #return time.mktime(date.timetuple())
    #s = sinceEpoch
    ep_fnc = lambda x: time.mktime(x.timetuple())
    
    retrnDates = np.zeros(len(dates))
    
    for i,sngDate in enumerate(dates):
        year = sngDate.year
        startOfThisYear = dt.datetime(year=year, month=1, day=1)
        startOfNextYear = dt.datetime(year=year+1, month=1, day=1)
    
        yearElapsed = ep_fnc(sngDate) - ep_fnc(startOfThisYear)
        yearDuration = ep_fnc(startOfNextYear) - ep_fnc(startOfThisYear)
        fraction = yearElapsed/yearDuration
        retrnDates[i] = sngDate.year + fraction

    return retrnDates

def t2dt(atime):
    """
    Convert atime (a float) to dt.datetime
    This is the inverse of dt2t.
    assert dt2t(t2dt(atime)) == atime
    """
    year = int(atime)
    remainder = atime - year
    boy = dt.datetime(year, 1, 1)
    eoy = dt.datetime(year + 1, 1, 1)
    seconds = remainder * (eoy - boy).total_seconds()
    return boy + dt.timedelta(seconds=seconds)

def perdelta(start, end, delta):
    curr = start
    while curr < end:
        yield curr
        curr += delta

def findCls(dataArray, val):
    ''' Returns the indice and closest value in dataArray to val'''
    return np.argmin(abs(val-dataArray))

    
def findCls2(dataArray, val):
    return min(dataArray, key=lambda p:dist_sq(val, p))


def ndaysAvg(data, dates, dateAxis=1, meanAxis=0, quad=0, deltaDays = dt.timedelta(days=7)):
    ''' Creates daily averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []
    dataMin = []
    days    = []

    #-----------------------------------------
    # Convert date time input to strictly date
    #-----------------------------------------
    dates_daily    = np.asarray([dt.date(d.year,d.month,d.day) for d in dates])

    dates_interval = [d for d in perdelta(dates_daily[0], dates_daily[-1], deltaDays)]
    Ndates_interval = len(dates_interval)




    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)


    #-------------------------
    # Loop through unique days
    #-------------------------
    for indDay in range(Ndates_interval-1):

        inds = np.where( (dates_daily >= dates_interval[indDay]) &  (dates_daily < dates_interval[indDay+1]) )[0]

        #print dates_daily[inds]
        #exit()

        if len(inds) >1:    
                            
       
            if data.ndim == 1:                                                              # For single dimension array   
                if quad == 1: dataAvg.append( np.sqrt(np.sum(data[inds]**2)) / len(inds) )  # Find daily average (sqrt(sum(x^2))/N)
                else:         dataAvg.append( np.mean(data[inds]) )                         # Find daily average (straight mean)
                std.append(np.std(data[inds]))                                              # Find standard deviation
                dataMin.append(np.min(data[inds]))
                days.append(dates_daily[inds[int(len(inds)/2)]])
                
            else:                                                                           # For multi-Dimension array
                s           = [slice(None)] * data.ndim
                s[dateAxis] = inds
                tempMat     = data[s]
                if quad == 1: dataAvg.append( np.sqrt(np.sum(tempMat,axis=meanAxis) / len(inds)))
                else:         dataAvg.append( np.mean(tempMat,axis=meanAxis) )
                std.append(np.std(tempMat,axis=meanAxis))                                  # Find standard deviation
                dataMin.append(np.min(tempMat,axis=meanAxis))
                days.append(dates_daily[inds[int(len(inds)/2)]])
            
            cnt.append(len(inds))                                                         # Number of observations used in daily average

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)
    dataMin = np.asarray(dataMin)
    days = np.asarray(days)


    outdata['Avg'] = dataAvg
    outdata['dates']    = days
    outdata['cnt']      = cnt
    outdata['std']      = std
    outdata['Min'] = dataMin

    return outdata


def dailyAvg(data, dates, dateAxis=1, meanAxis=0, quad=0):
    ''' Creates daily averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []
    dataMin = []

    #-----------------------------------------
    # Convert date time input to strictly date
    #-----------------------------------------
    dates_daily = np.asarray([dt.date(d.year,d.month,d.day) for d in dates])

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #--------------------------------------
    # Create a list of days to loop through
    #--------------------------------------
    uniqueDays = list(set(dates_daily))          # Find a list of unique days
    uniqueDays.sort()
    uniqueDays = np.asarray(uniqueDays)    

    #-------------------------
    # Loop through unique days
    #-------------------------
    for indDay in uniqueDays:

        inds = np.where(dates_daily == indDay)[0]                                       # Index for values to use as daily average
       
        if data.ndim == 1:                                                              # For single dimension array   
            if quad == 1: dataAvg.append( np.sqrt(np.sum(data[inds]**2)) / len(inds) )  # Find daily average (sqrt(sum(x^2))/N)
            else:         dataAvg.append( np.mean(data[inds]) )                         # Find daily average (straight mean)
            std.append(np.std(data[inds]))                                              # Find standard deviation
            dataMin.append(np.min(data[inds]))
            
        else:                                                                           # For multi-Dimension array
            s           = [slice(None)] * data.ndim
            s[dateAxis] = inds
            tempMat     = data[s]
            if quad == 1: dataAvg.append( np.sqrt(np.sum(tempMat,axis=meanAxis) / len(inds)))
            else:         dataAvg.append( np.mean(tempMat,axis=meanAxis) )
            std.append(np.std(tempMat,axis=meanAxis))                                  # Find standard deviation
            dataMin.append(np.min(tempMat,axis=meanAxis))
        
        cnt.append(len(inds))                                                         # Number of observations used in daily average

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)
    dataMin = np.asarray(dataMin)


    outdata['dailyAvg'] = dataAvg
    outdata['dates']    = uniqueDays
    outdata['cnt']      = cnt
    outdata['std']      = std
    outdata['dailyMin'] = dataMin

    return outdata


def mnthlyAvg(data, dates, dateAxis=1, meanAxis=0, quad=0):
    ''' Creates monthly averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []

    #-----------------------------------------
    # Convert date time input to strictly date
    #-----------------------------------------
    dates_mnthly = np.asarray([dt.date(d.year,d.month,1) for d in dates])

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #----------------------------------------
    # Create a list of months to loop through
    #----------------------------------------
    uniqueMnths = list(set(dates_mnthly))          # Find a list of unique months
    uniqueMnths.sort()
    uniqueMnths = np.asarray(uniqueMnths)    

    #---------------------------
    # Loop through unique months
    #---------------------------
    for indDay in uniqueMnths:

        inds = np.where(dates_mnthly == indDay)[0]                                      # Index for values to use as monhtly average
        
        if data.ndim == 1:                                                              # For single dimension array   
            if quad == 1: dataAvg.append( np.sqrt(np.sum(data[inds]**2)) / len(inds) )  # Find monthly average (sqrt(sum(x^2))/N)
            else:         dataAvg.append( np.mean(data[inds]) )                         # Find monthly average (straight mean)
            std.append(np.std(data[inds]))                                              # Find standard deviation
            
        else:                                                                           # For multi-Dimension array
            s           = [slice(None)] * data.ndim
            s[dateAxis] = inds
            tempMat     = data[s]
            if quad == 1: dataAvg.append( np.sqrt(np.sum(tempMat,axis=meanAxis) / len(inds)))
            else:         dataAvg.append( np.mean(tempMat,axis=meanAxis) )
            std.append(np.std(tempMat,axis=meanAxis))                                  # Find standard deviation    

        cnt.append(len(inds))                                                          # Number of observations used in monhtly average

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)

    outdata['mnthlyAvg'] = dataAvg
    outdata['dates']     = uniqueMnths
    outdata['cnt']       = cnt
    outdata['std']       = std

    return outdata


def allMnthAvg(data, dates):
    ''' Creates monthly averages of specified quantity'''

    #-----------------------------
    # Initialize output dictionary
    #-----------------------------
    outdata = {}
    dataAvg = []
    cnt     = []
    std     = []
    mnths   = np.array(range(1,13))

    #-----------------------------
    # Ensure data is a numpy array
    #-----------------------------
    data = np.asarray(data)

    #------------------------------
    # Loop through months (Jan-Dec)
    #------------------------------
    monthlist = np.array([d.month for d in dates])   # Create a month list

    for indMnth in mnths:
        inds = np.where(monthlist == indMnth)[0]     # Index for values to use as monhtly average
        dataAvg.append(np.mean(data[inds]))          # Find monhtly average
        cnt.append(len(inds))                        # Number of observations used in monhtly average
        std.append(np.std(data[inds]))               # Find standard deviation

    cnt     = np.asarray(cnt)
    dataAvg = np.asarray(dataAvg)
    std     = np.asarray(std)

    outdata['mnthlyAvg'] = dataAvg
    outdata['month']     = mnths
    outdata['cnt']       = cnt
    outdata['std']       = std

    return outdata

def clrplt():
    ''' Colors for plots in order for for loops'''
    
    clr = ('green', 'red', 'maroon', 'blue', 'gray', 'yellow', 'cyan', 'orange', 'black')
    return clr


def readParm(fname):
    ''' Function to read parm.output file '''

    #--------------
    # Read all file
    #--------------
    with open(fname, 'r' ) as fopen:
        lines = fopen.readlines()

    params = lines[3].strip().split()
    data   = np.array([[float(x) for x in row.split()[1:]] for row in lines[4:]])

    #-----------------------------------------------------------
    # Create a dictionary where keys are header of parm.out file
    # A list of numpy arrays is created for repeated keys
    #-----------------------------------------------------------
    parm = {}
    for k in set(params):
        inds = [i for i, val in enumerate(params) if val == k]
        parm.setdefault(k,[]).append(data[:,inds])

    #--------------------------------------
    # Un-nest numpy arrays in parm dictionary
    #--------------------------------------
    for k in parm: parm[k] = parm[k][0]   # Check this

    return parm


def wet_to_dry_vmr(f_wet,f_H2O):
    ''' Convert wet VMR to dry VMR '''
    
    #----------------------------------------
    # Convert wet vmr to dry:
    #   --- f_gas(dry)= f_gas(wet)/ (1-f_H2O)
    #----------------------------------------
    f_dry    = f_wet / (1.0 - f_H2O)    
    
    return f_dry


def lyrVMR(lAlt,hAlt,Z,prf,airmass,sclfct):
    ''' Find layer mixing ratio and convert to ppx (based on scale factor sclfct)'''
    
    Zlwr     = (np.abs(Z-lAlt)).argmin()
    Zupr     = (np.abs(Z-hAlt)).argmin()
    part_VMR = np.asarray( np.average(rprf1[Zlwr:Zupr,:], axis=0, weights=Airmass1[Zlwr:Zupr,:] ) ) * sclfct 
    
    return part_VMR

def readCtlF(ctlF):
    ''' Read .ctl file'''
    
    ctl = {}
    lines = tryopen(ctlF)
    
    if lines:
        gas_flg = True
        
        for line in lines:
            line = line.strip()  

            #---------------------------------------------
            # Lines which start with comments and comments
            # embedded within lines
            #---------------------------------------------                
            if line.startswith('#'): continue                # For lines which start with comments
            if '#' in line: line = line.partition('#')[0]    # For comments embedded in lines

            #-------------------
            # Handle empty lines
            #-------------------
            if len(line) == 0: continue

            #--------------------------
            # Populate input dictionary
            #--------------------------
            if '=' in line:
                lhs,_,rhs = line.partition('=')
                lhs       = lhs.strip()
                rhs       = rhs.strip().split()


                ctl[lhs] = [convrtD(val) for val in rhs]
                #for rhs_part in rhs:                          
                    #if 'd' in rhs_part.lower():                               # Test and handle strings containing D (ie 1.0D0)
                        #mant,trsh,expo = rhs_part.lower().partition('d')
                        #try:
                            #rhs_part = float(mant)*10**int(expo)
                        #except ValueError:
                            #pass
                    #else:                                                # Handle strings => number
                        #try:
                            #rhs_part = [float(ind) for ind in rhs_part]
                        #except ValueError:
                            #pass
                    #ctl[lhs] = rhs

            else:                                                        # Handle multiple line assignments
                rhs = line.strip().split()             

                try:
                    rhs = [float(ind) for ind in rhs]
                except ValueError:
                    pass

                ctl[lhs] += rhs          

            #----------------------    
            # Primary Retrieval Gas
            #----------------------
            # Search for primary retrieval gas in gas.column.list or gas.profile.list
            # The first gas listed is considered the primary retrieval gas
            if gas_flg:
                match = re.match(r'\s*gas\.\w+\.list\s*=\s*(\w+)', line)
                if match:
                    PrimaryGas = match.group(1)
                    gas_flg = False

            #----------------   
            # Spectral Region
            #----------------                
            # Find number of bands
            #match = re.match(r'\s*band\s*=\s*(.+)',line)
            #if match:
                #bands = match.group(1).split()
                #bnd_flg = True

            # Find the upper and lower bands for set including all regions    
            #if bnd_flg:
                #for band in bands:
                    #match = re.match(r'\s*band.' + band + '.nu_\w+\s*=\s*(.*)',line)
                    #if match:
                        #nu.append(float(match.group(1)))

    #nu.sort()              # Sort wavenumbers
    #nuUpper = nu[-1]  # Pull off upper wavenumber
    #nuLower = nu[0]   # Pull off lower wavenumber  
    
    if not gas_flg: return (PrimaryGas,ctl)
    else:           return ctl
    

#---------------------------------
# A simple code for finding linear 
# trend with prediciton intervals
#---------------------------------



#-------------------------------------
# Code to do boot strap trend analysis
#-------------------------------------
def fourier_basis(x, degree, half_period):
    """Returns a 2-d array of fourier basis."""
    A = np.ones((x.size, 2 * degree + 1))
    
    for d in range(1, degree + 1):
        A[:, 2*d-1] = np.cos(d * np.pi * x / half_period)
        A[:, 2*d]   = np.sin(d * np.pi * x / half_period)

    return A

def basis_qbo(x, degree, half_period, x2):
    """Returns a 2-d array of fourier basis."""
    A = np.ones((x.size, 2 * degree + 1))
    
    for d in range(1, degree + 1):
        A[:, 2*d-1] = np.cos(d * np.pi * x / half_period)
        A[:, 2*d]   = np.sin(d * np.pi * x / half_period)

    #-------------------------------------
    # QBO
    #-------------------------------------

    Dir        = '/data1/projects/ocs/'                 # Name of station location
    file       = 'QBO.dat' 

    fileName  = Dir + file

    cols, indexToName = getColumns(fileName, headerrow=8, delim=' ', header=True)

    site  = cols[indexToName[0]]
    YYMM  = cols[indexToName[1]]
    qbo50 = cols[indexToName[3]]
    qbo30 = cols[indexToName[5]]

    qboDate = []

    for i in YYMM:
        YY = i[0:2]
        MM = i[2:4]

        if int(YY) <= 20 :
            YYYY = '20'+YY
        else:
            YYYY = '19'+YY

        qboDate.append(dt.date(int(YYYY), int(MM), 15))

    qboDate     = np.asarray(qboDate)
    qbo50       = np.asarray(qbo50, dtype=float) + 500.
    qbo30       = np.asarray(qbo30, dtype=float) + 500.

    qboDateFrac = toYearFraction(qboDate)

    xminqbo     = qboDateFrac.min()
    xnormqbo    = qboDateFrac - xminqbo

    qbo30Int = np.asarray(interpolate.interp1d(xnormqbo, qbo30, kind='linear')(x))

    return qbo30Int

    #-------------------------------------
    # Return
    #-------------------------------------
    # B = np.ones((x.size, 2 * degree + 1))

    # for i in range(0, 2 * degree + 1):
    #     B[:, i] = A[:, i]*qbo30Int
    
    # return B

def fit_driftfourier_poly(x, data, weights, degree, half_period=0.5):
    """
    Fit y = f(x - x.min()) to data where f is given by
    fourier series + drift.
    
    Parameters
    ----------
    x : 1-d array
        x-coordinates
    data : 1-d array
        data values
    weights : 1-d array
        weights (>=0)
    degree : int
        degree of fourier series
    half_period : float
        half period
    
    Returns
    -------
    intercept : float
        intercept at x.min()
    slope : float
        slope (drift) for the normalized data
        (x - x.min())
    pfourier : 1-d array
        Fourier series parameters for the
        normalized data
    f_drift : callable
        Can be used to calculate the drift
        given any (non-normalized) x
    f_fourier : callable
        Can be used to calculate fourier series
    f_driftfourier : callable
        Can be used to calculate drift + fourier
    residual_std : float
        estimated standard deviation of residuals
    A : 2-d array
        matrix of "coefficients"
    
    """
    xmin = x.min()
    xnorm = x - xmin
    
    # coefficient matrix
    A = np.ones((x.size, 2 * degree + 6))
    A[:, 1] = xnorm
    A[:, 2] = xnorm**2
    A[:, 3] = xnorm**3
    A[:, 4] = xnorm**4
    A[:, 5] = xnorm**5
    A[:, 6:] = fourier_basis(xnorm, degree, half_period)[:, 1:]
    
    # linear weighted least squares
    results = np.linalg.lstsq(A * weights[:, np.newaxis],
                              data * weights)

    params    = results[0]
    
    intercept = params[0]
    slope     = params[1]
    poly      = params[2]
    poly2     = params[3]
    poly3     = params[4]
    poly4     = params[5]
    pfourier  = params[6:]

    
    #f_drift   = lambda t: slope * (t - xmin) + poly * (t - xmin)**2 + poly2 * (t - xmin)**3 + intercept
    #f_drift   = lambda t: (slope *t - slope*xmin) +  (poly*( t**2 - 2.*t*xmin + xmin**2)) +  (poly2*( t**3 - 3.*t**2*xmin + 3.*t*xmin**2 - xmin**3))  + intercept #+(poly3*( t**4 - 4.*t**3*xmin + 6.*t**2xmin**2 - 4*t*xmin**3 + xmin**4)) +(poly4*( t**5 - 5.*t**4*xmin + 10.*t**3*xmin**2 - 10*t**2*xmin**3 + 5*t*xmin**4 + xmin**5)) 

    f_drift   = lambda t: slope * (t - xmin) + poly * (t - xmin)**2 + poly2 * (t - xmin)**3 + poly3 * (t - xmin)**4 +  poly4 * (t - xmin)**5 + intercept 
    
    ####df_drift  = lambda t: slope + (poly*(2.*t.max() - 2.*t.min())) + (poly2*(3.*t.max()**2 - 6.*t.max()*t.min() + 3.*t.min()**2))
    #df_drift  = lambda t: slope + (poly*(2.*t - 2.*xmin)) + (poly2*(3.*t**2 - 6.*t*xmin + 3.*xmin**2))
    df_drift  = lambda t: slope + (poly*(2.*t - 2.*xmin)) + (poly2*(3.*t**2 - 6.*t*xmin + 3.*xmin**2)) + (poly3*(4.*t**3 - 12.*t**2*xmin + 12.*t*xmin**2 - 4*xmin**3)) + (poly4*(5.*t**4 - 20.*t**3*xmin + 30.*t**2*xmin**2 - 20*t*xmin**3 +5*xmin**4))

    f_fourier = lambda t: np.sum(fourier_basis(t - xmin, degree,
                                               half_period)[:, 1:]
                                 * pfourier[np.newaxis, :],
                                 axis=1) + intercept

    f_driftfourier = lambda t: f_drift(t) + f_fourier(t) - intercept


    
    residual_std = np.sqrt(results[1][0] / (x.size - 2 * degree + 2)) 
    
    return (intercept, slope, pfourier,
            f_drift, f_fourier, f_driftfourier,
            residual_std, A, df_drift)


def fit_driftfourier(x, data, weights, degree, half_period=0.5):
    """
    Fit y = f(x - x.min()) to data where f is given by
    fourier series + drift.
    
    Parameters
    ----------
    x : 1-d array
        x-coordinates
    data : 1-d array
        data values
    weights : 1-d array
        weights (>=0)
    degree : int
        degree of fourier series
    half_period : float
        half period
    
    Returns
    -------
    intercept : float
        intercept at x.min()
    slope : float
        slope (drift) for the normalized data
        (x - x.min())
    pfourier : 1-d array
        Fourier series parameters for the
        normalized data
    f_drift : callable
        Can be used to calculate the drift
        given any (non-normalized) x
    f_fourier : callable
        Can be used to calculate fourier series
    f_driftfourier : callable
        Can be used to calculate drift + fourier
    residual_std : float
        estimated standard deviation of residuals
    A : 2-d array
        matrix of "coefficients"
    
    """
    xmin = x.min()
    xnorm = x - xmin
    
    # coefficient matrix
    A = np.ones((x.size, 2 * degree + 2))
    A[:, 1] = xnorm
    A[:, 2:] = fourier_basis(xnorm, degree, half_period)[:, 1:]
    
    # linear weighted least squares
    results = np.linalg.lstsq(A * weights[:, np.newaxis],
                              data * weights)

    params    = results[0]
    intercept = params[0]
    slope     = params[1]
    pfourier  = params[2:]
    
    f_drift   = lambda t: slope * (t - xmin) + intercept

    df_drift  = lambda t: slope*t/t
    
    f_fourier = lambda t: np.sum(fourier_basis(t - xmin, degree,
                                               half_period)[:, 1:]
                                 * pfourier[np.newaxis, :],
                                 axis=1) + intercept
    f_driftfourier = lambda t: f_drift(t) + f_fourier(t) - intercept
    
    residual_std = np.sqrt(results[1][0] / (x.size - 2 * degree + 2)) 
    
    return (intercept, slope, pfourier,
            f_drift, f_fourier, f_driftfourier,
            residual_std, A, df_drift)

def fit_qbo(x, data, weights, degree, half_period=0.5):
    """
    Fit y = f(x - x.min()) to data where f is given by
    fourier series + drift.
    
    Parameters
    ----------
    x : 1-d array
        x-coordinates
    data : 1-d array
        data values
    weights : 1-d array
        weights (>=0)
    degree : int
        degree of fourier series
    half_period : float
        half period
    
    Returns
    -------
    intercept : float
        intercept at x.min()
    slope : float
        slope (drift) for the normalized data
        (x - x.min())
    pfourier : 1-d array
        Fourier series parameters for the
        normalized data
    f_drift : callable
        Can be used to calculate the drift
        given any (non-normalized) x
    f_fourier : callable
        Can be used to calculate fourier series
    f_driftfourier : callable
        Can be used to calculate drift + fourier
    residual_std : float
        estimated standard deviation of residuals
    A : 2-d array
        matrix of "coefficients"
    
    """
    xmin  = x.min()
    xnorm = x - xmin
    
    # coefficient matrix
    #A = np.ones((x.size, 2 * degree + 2))
    A       = np.ones((x.size, 1 + 2))
    A[:, 1] = xnorm
    A[:, 2] = basis_qbo(xnorm, degree, half_period, x)
    
    # linear weighted least squares
    results = np.linalg.lstsq(A * weights[:, np.newaxis],
                              data * weights)


    params    = results[0]
    intercept = params[0]
    slope     = params[1]
    pfourier  = params[2]
    
    f_drift   = lambda t: slope * (t - xmin) + intercept

    df_drift  = lambda t: slope*t/t
    
    f_fourier = lambda t: basis_qbo(t - xmin, degree,half_period, t)* pfourier + intercept

    #f_fourier = lambda t: np.sum(fourier_basis(t - xmin, degree,
    #                                           half_period)[:, 1:]
    #                             * pfourier[np.newaxis, :],
    #                             axis=1) + intercept
    
    f_driftfourier = lambda t: f_drift(t) + f_fourier(t) - intercept
    
    residual_std = np.sqrt(results[1][0] / (x.size - 2 * degree + 2)) 
    
    return (intercept, slope, pfourier,
            f_drift, f_fourier, f_driftfourier,
            residual_std, A, df_drift)



def cf_driftfourier(x, data, weights, degree,
                    half_period=0.5, nboot=5000,
                    percentiles=(2.5, 50., 97.5)):
    """
    Calculate confidence intervals for the fitted
    parameters from fourier series + drift modelling,
    using bootstrap resampling.
    
    Parameters
    ----------
    nboot : int
        number of bootstrap replicates
    percentiles : sequence of floats
        percentiles of parameter estimate
        distributions to return 
    
    Returns
    -------
    perc : dict
        percentiles for of each parameter
        distribution
    intercept : 1-d array
        intercept estimates from bootstraped
        datasets.
    slope : 1-d array
        slope estimates
    pfourier : 2-d array
        fourier parameters estimates
    
    See Also
    --------
    :func:`fit_driftfourier`
    """
    
    # 1st fit without bootstraping
    results = fit_driftfourier(x, data, weights,
                               degree, half_period)
    f_driftfourier = results[5]
    A = results[7]
    model = f_driftfourier(x)
    residuals = data - model
    
    # generate bootstrap resamples of residuals
    # and new datasets from these resamples
    boot_dataset = np.empty((x.size, nboot))
    for i in range(nboot):
        resample_i = np.floor(np.random.rand(x.size) * x.size).astype(int)
        resample_residuals = residuals[resample_i]
        boot_dataset[:, i] = model + resample_residuals
    
    # fit all bootstrap datasets
    results_boot = np.linalg.lstsq(A * weights[:, np.newaxis],
                                   boot_dataset * weights[:, np.newaxis])
    
    params_boot = results_boot[0]
    
    # compute percentiles
    perc_boot = np.column_stack(np.percentile(params_boot,
                                              percentiles, axis=1))
    
    perc = {'intercept' : perc_boot[0],
            'slope' : perc_boot[1],
            'pfourier' : perc_boot[2:]}
    
    intercept = params_boot[0]
    slope = params_boot[1]
    pfourier = params_boot[2:]
    
    return perc, intercept, slope, pfourier


def weatherout(loc, datadir, dataFileTag, iyear, fyear):
    ''' Read met data (yearly) data from eol at FL0  (initial files are created with the read_FL0_EOL_data.py'''

    tempList       = []
    rhList         = []
    pressList      = []
    tempDPList     = []
    wdirList       = []
    wspdList       = []
    wmaxList       = []
    datet_weather  = []


    nyears = int(fyear) - int(iyear) + 1

    for i in range(nyears):

        #------------------------
        # Open file and read data
        #------------------------
        #fileName = loc + '_met_data_'+str(iyear + i)+'_'+dataFileTag+'.txt'
        fileName = glob.glob(datadir + loc + '_met_data_'+str(iyear + i)+'_'+dataFileTag+'.txt')

        if not fileName:
            continue
        
        fileName = fileName[0]
        print 'reading met data: %s' %fileName

        f = open(fileName, 'r')
        header1 = f.readline()
       
        temp_i          = []
        rh_i            = []
        press_i         = []
        tempdp_i        = []
        wdir_i          = []
        wspd_i          = []
        wmax_i          = []
        datet_weather_i = []
        
        for line in f:
            line = line.strip()
            columns = line.split()
            
            year = int(columns[0])
            month = int(columns[1])
            day = int(columns[2])
            hour = int(columns[3])
            minute = int(columns[4])
            datet_weather_i.append(dt.datetime(year,month,day,hour,minute))
            temp_i.append(float(columns[5]))
            rh_i.append(float(columns[6]))
            press_i.append(float(columns[7]))
            tempdp_i.append(float(columns[8]))
            wdir_i.append(float(columns[9]))
            wspd_i.append(float(columns[10]))
            wmax_i.append(float(columns[11]))

        tempList.extend(temp_i)
        rhList.extend(rh_i)
        pressList.extend(press_i)
        tempDPList.extend(tempdp_i)
        wdirList.extend(wdir_i)
        wspdList.extend(wspd_i)
        wmaxList.extend(wmax_i)
        datet_weather.extend(datet_weather_i)

    return wdirList, wspdList, tempList, rhList, datet_weather
             
   

def getgasname(rawname):

    if rawname.upper() == 'C2H6': gasname = 'C$_2$H$_6$'
    elif rawname.upper() == 'CH4': gasname = 'CH$_4$'
    elif rawname.upper() == 'C2H2': gasname = 'C$_2$H$_2$'
    elif rawname.upper() == 'NH3': gasname = 'NH$_3$'
    elif rawname.upper() == 'O3': gasname = 'O$_3$'
    elif rawname.upper() == 'H2CO': gasname = 'CH$_2$O'
    elif rawname.upper() == 'HCL': gasname = 'HCl'
    elif rawname.upper() == 'HNO3': gasname = 'HNO$_3$'
    elif rawname.upper() == 'N2O': gasname = 'N$_2$O'
    elif rawname.upper() == 'CLONO2': gasname = 'ClONO$_2$'
    elif rawname.upper() == 'H2O': gasname = 'H$_2$O'
    elif rawname.upper() == 'NO2': gasname = 'NO$_2$'
    elif rawname.upper() == 'CCL4': gasname = 'CCl$_4$'
    elif rawname.upper() == 'CO': gasname = 'CO'
    elif rawname.upper() == 'HCHO': gasname = 'H2CO'
    elif rawname.upper() == 'H2CO': gasname = 'H2CO'
    elif rawname.upper() == 'HCOOH': gasname = 'HCOOH'
    elif rawname.upper() == 'HCN': gasname = 'HCN'
    elif rawname.upper() == 'COF2': gasname = 'COF$_2$'
    else: gasname = rawname.upper()

    return gasname

def getColumns(FileName, headerrow=0, delim=' ', header=True):
    """
    Get columns of data from inFile. The order of the rows is respected
    
    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that 
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in 
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the 
    heading names (i.e. the keys in the cols dict)
    """

    f = open(FileName, 'r')
    inFile = f.readlines()


    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum < headerrow: continue

        if lineNum == headerrow:
            if delim == ' ': headings = line.split()
            else: headings = line.split(delim) 
              
            i = 0
            for heading in headings:
                heading = heading.strip()
         
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            if delim == ' ': cells = line.split()
            else: cells = line.split(delim)

            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1
                
    return cols, indexToName

def read_ICARTT(path):
    """
    Get ICARTT information and data.
    
    Note: some icartt file may deiffer, hence error may occur
    """

    print 'Reading File:', str(path)

    PI_LINE         = 2
    ORG_LINE        = 3
    PLAT_LINE       = 4
    MISSION_LINE    = 5
    VOL_LINE        = 6
    DATE_LINE       = 7
    TIME_INT_LINE   = 8
    UNIT_LINE       = 9
    DATE_VAR_LINE   = 10
    SCALE_LINE      = 11
    MISSING_LINE    = 12

    f = open(path, 'rU')
    missing = []
    units = []
    l = f.readline()
    if ',' in l:
        delim = ','
    else:
        delim = None
    split = lambda s: list(map(str.strip, s.split(delim)))

    if split(l)[-1] != '1001':
        raise TypeError("File is the wrong format.  Expected 1001; got %s" % (split(l)[-1],))
    
    n, fmt = split(l)
    n_user_comments = 0
    n_special_comments = 0
    n_header_lines = int(n)

    try:
        for li in range(n_header_lines-1):
            li += 2
            l = f.readline()
            LAST_VAR_DESC_LINE = 12+len(missing)
            SPECIAL_COMMENT_COUNT_LINE = LAST_VAR_DESC_LINE + 1
            LAST_SPECIAL_COMMENT_LINE = SPECIAL_COMMENT_COUNT_LINE + n_special_comments
            USER_COMMENT_COUNT_LINE = 12+len(missing)+2+n_special_comments
            if li == PI_LINE:
                PI_NAME = l.strip()
            elif li == ORG_LINE:
                ORGANIZATION_NAME = l.strip()
            elif li == PLAT_LINE:
                SOURCE_DESCRIPTION = l.strip()
            elif li == MISSION_LINE:
                MISSION_NAME = l.strip()
            elif li == VOL_LINE:
                VOLUME_INFO = l.strip()
            elif li == DATE_LINE:
                l = l.replace(',', '').split()

                SDATE = "".join(l[:3])
                WDATE = "".join(l[3:])
                SDATE = SDATE
                WDATE = WDATE
                SDATE = datetime.strptime(SDATE, '%Y%m%d')
                WDATE = datetime.strptime(WDATE, '%Y%m%d')
            elif li == TIME_INT_LINE:
                TIME_INTERVAL = l.strip()
            elif li == UNIT_LINE:
                units.append(l.replace('\n', '').replace('\r', '').strip())
                INDEPENDENT_VARIABLE = units[-1]
            elif li == SCALE_LINE:
                scales = [eval(i) for i in split(l)]
                if set([float(s) for s in scales]) != set([1.]):
                    raise ValueError("Unsupported: scaling is unsupported.  data is scaled by %s" % (str(scales),))
            elif li == MISSING_LINE:
                missing = [eval(i) for i in split(l)]
            elif li > MISSING_LINE and li <= LAST_VAR_DESC_LINE:
                nameunit = l.replace('\n','').split(',')
                name = nameunit[0].strip()
                if len(nameunit) > 1:
                    units.append(nameunit[1].strip())
                elif re.compile('(.*)\((.*)\)').match(nameunit[0]):
                    desc_groups = re.compile('(.*)\((.*)\).*').match(nameunit[0]).groups()
                    name = desc_groups[0].strip()
                    units.append(desc_groups[1].strip())
                elif '_' in name:
                    units.append(name.split('_')[1].strip())
                else:
                    warn('Could not find unit in string: "%s"' % l)
                    units.append(name.strip())
            elif li == SPECIAL_COMMENT_COUNT_LINE:
                n_special_comments = int(l.replace('\n', ''))
            elif li > SPECIAL_COMMENT_COUNT_LINE and li <= LAST_SPECIAL_COMMENT_LINE:
                pass
            elif li == USER_COMMENT_COUNT_LINE:
                n_user_comments = int(l.replace('\n',''))
            elif li > USER_COMMENT_COUNT_LINE and li < n_header_lines:
                colon_pos = l.find(':')
                k = l[:colon_pos].strip()
                v = l[colon_pos+1:].strip()
                #setattr(self,k,v)
            elif li == n_header_lines:
                variables = l.replace(',','').split()
                TFLAG = variables[0]
    except Exception as e:
        raise SyntaxError("Error parsing icartt file %s: %s" % (path, repr(e)))


    infile = file(path, 'r')
    cols, indextoname = getColumns(path, headerrow=n_header_lines-1, delim=',', header=True)
    infile.close()

    start_utc = np.asarray(cols[indextoname[0]], dtype=np.float32)
    start_utc.astype(int)

    dateutc = SDATE + vectorize(lambda s: timedelta(seconds = int(s), microseconds = (s - int(s)) * 1.E6 ))(start_utc).view(type = ndarray)
    
    keys = {'PI_NAME':PI_NAME, 'ORGANIZATION_NAME':ORGANIZATION_NAME, 'SOURCE_DESCRIPTION':SOURCE_DESCRIPTION, 'MISSION_NAME':MISSION_NAME,
            'VOLUME_INFO':VOLUME_INFO, 'TIME_INTERVAL':TIME_INTERVAL, 'MISSING':missing}
    
    return keys, cols, dateutc, indextoname



def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def saveFlg(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    #if verb:
    #    print "NetCDF Global Attributes:"
    #    for nc_attr in nc_attrs:
    #        print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            saveFlg(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                saveFlg(var)
    return nc_attrs, nc_dims, nc_vars


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2.)**2 + cos(lat1) * cos(lat2) * sin(dlon/2.)**2
    c = 2.0 * asin(sqrt(a)) 
    r = 6371.0 # Radius of earth in kilometers. Use 3956 for miles
    return c * r

def destination(lat, lon, distance, bearing, radius=6371.008771415):
    '''Return the destination LatLon from this point after having
       travelled the given distance on the given initial bearing.
       @param {number} distance - Distance travelled.
       @param {degrees} bearing - Initial bearing in degrees from North.
       @param {number} [radius=R_M] - Mean radius of earth (default.

    Notes: see https://github.com/mrJean1/PyGeodesy/blob/master/geodesy/sphericalTrigonometry.py
    '''
    d        = float(distance)/float(radius) 
    t        = radians(bearing)
    latRad   = radians(lat)
    lonRad   = radians(lon)

    ca1, sa1 = cos(latRad), sin(latRad)
    cd, sd   = cos(d), sin(d) 

    a2 = asin(sa1 * cd + ca1 * sd * cos(t))
    b2 = atan2(sin(t) * sd * ca1, cd - sa1 * sin(a2)) + lonRad

    lat2  = degrees(a2)
    lon2  = degrees(b2)

    return lat2, lon2

def readRefPrf(fname='', parms=''):

        ''' Reads in reference profile, an input file for sfit4 (raytrace) '''
        refPrf = {}

        
        try:
            with open(fname,'r') as fopen: lines = fopen.readlines()
                            
            #----------------------------------------
            # Get Altitude, Pressure, and Temperature
            # from reference.prf file
            #----------------------------------------
            nlyrs  = int(lines[0].strip().split()[1])
            nlines = int(np.ceil(nlyrs/5.0))
            
            for ind,line in enumerate(lines):
                if any(p in line for p in parms):

                    val = [x for x in parms if x in line][0]
                   
                    refPrf.setdefault(val,[]).append([float(x[:-1]) for row in lines[ind+1:ind+nlines+1] for x in row.strip().split()])

        except Exception as errmsg:
            print errmsg
        
        #------------------------
        # Convert to numpy arrays
        # and sort based on date
        #------------------------
        for k in refPrf:
            refPrf[k] = np.asarray(refPrf[k])

        return refPrf

def AirDensity(Press='', Temp='', Punits=False, Tunits=False):

    
    if Punits == 'hPa':      Press  = Press*100.0
    elif Punits == 'atm':    Press  = Press*101325.0
    elif Punits == 'Pa':     Press  = Press
    else:
        'Units of Pressure are not given'
        exit()
        
    if Tunits == 'C': Temp  = Temp + 273.15
    elif Tunits == 'K': Temp  = Temp 
    else:
        'Units of Temperature are not given'
        exit()

    #-------------------------------------------------
    # Calculating density profile, layer thickness and midpoint
    #-------------------------------------------------
    rho =  np.true_divide(Press, (Temp*8.314)  ) * 6.022e23 /100. /100. /100.  ##in molec/cm3
    return rho


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
        r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
   3     The Savitzky-Golay filter removes high frequency noise from data.
   4     It has the advantage of preserving the original shape and
   5     features of the signal better than other types of filtering
   6     approaches, such as moving averages techniques.
   7     Parameters
   8     ----------
   9     y : array_like, shape (N,)
  10         the values of the time history of the signal.
  11     window_size : int
  12         the length of the window. Must be an odd integer number.
  13     order : int
  14         the order of the polynomial used in the filtering.
  15         Must be less then `window_size` - 1.
  16     deriv: int
  17         the order of the derivative to compute (default = 0 means only smoothing)
  18     Returns
  19     -------
  20     ys : ndarray, shape (N)
  21         the smoothed signal (or it's n-th derivative).
  22     Notes
  23     -----
  24     The Savitzky-Golay is a type of low-pass filter, particularly
  25     suited for smoothing noisy data. The main idea behind this
  26     approach is to make for each point a least-square fit with a
  27     polynomial of high order over a odd-sized window centered at
  28     the point.
  29     Examples
  30     --------
  31     t = np.linspace(-4, 4, 500)
  32     y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
  33     ysg = savitzky_golay(y, window_size=31, order=4)
  34     import matplotlib.pyplot as plt
  35     plt.plot(t, y, label='Noisy signal')
  36     plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
  37     plt.plot(t, ysg, 'r', label='Filtered signal')
  38     plt.legend()
  39     plt.show()
  40     References
  41     ----------
  42     .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
  43        Data by Simplified Least Squares Procedures. Analytical
  44        Chemistry, 1964, 36 (8), pp 1627-1639.
  45     .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
  46        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
  47        Cambridge University Press ISBN-13: 9780521880688
  48     """
    
        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError, msg:
            raise ValueError("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2
        # precompute coefficients
        b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
        lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve( m[::-1], y, mode='valid')





def orthoregress(x, y,  xerr=False, yerr=False, InError=False):
    """Perform an Orthogonal Distance Regression on the given data,
    using the same interface as the standard scipy.stats.linregress function.
    Arguments:
    x:    x data
    y:    y data
    sx:   The data, with weightings as actual standard deviations and/or covariances
    sy:   The data, with weightings as actual standard deviations and/or covariances
    Returns:
    [m, c, nan, nan, nan]
    Uses standard ordinary least squares to estimate the starting parameters
    then uses the scipy.odr interface to the ODRPACK Fortran code to do the
    orthogonal distance calculations.

    Sources: https://docs.scipy.org/doc/scipy-0.9.0/reference/odr.html
             https://docs.scipy.org/doc/scipy/reference/odr.html
    """

    if ( x.shape[0] != y.shape[0] ):
        print 'ODR error: Data sets must have same size in: xData = {}, yData = {}'.format(x.shape[0],y.shape[0])
        sys.exit()

    if InError:
        if ( x.shape[0] != xerr.shape[0] ):
            print 'ODR error: Data sets must have same size in: xData = {}, xError = {}'.format(x.shape[0],xerr.shape[0])
            sys.exit()

        if ( y.shape[0] != yerr.shape[0] ):
            print 'ODR error: Data sets must have same size in: yData = {}, yError = {}'.format(y.shape[0],yerr.shape[0])
            sys.exit()

        if ( xerr.shape[0] != yerr.shape[0] ):
            print 'ODR error: Data sets must have same size in: xError = {}, yError = {}'.format(xerr.shape[0],yerr.shape[0])
            sys.exit()

   
    linreg  = linregress(x, y)
    mod     = Model(f)
    
    if InError: 
        dat = RealData(x, y, sx=xerr, sy=yerr)
    else:   
        dat = Data(x, y)
   
    od = ODR(dat, mod, beta0=linreg[0:2])
    out = od.run()

    #out.pprint()  Pretty-print important results

    if InError:
        return out.beta, out.sd_beta
    else:
        return out.beta

    #return list(out.beta) + [np.nan, np.nan, np.nan]


def f(p, x):
    """Basic linear regression 'model' for use with ODR"""
    return (p[0] * x) + p[1]


def running_mean(x, N):
     cumsum = np.cumsum(np.insert(x, 0, 0)) 
     return (cumsum[N:] - cumsum[:-N]) / N


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    #return y
    return y[(window_len/2-1):-(window_len/2)]

def getseason(dates):
    
    #seas2get = {'Summer':(dt.date(2014,6,21), dt.date(2014,9,22)),
    #           'Autumn':(dt.date(2014,9,23), dt.date(2014,12,20)),
    #           'Spring':(dt.date(2014,3,21), dt.date(2014,6,20))}

    seas2get = {'Summer':(dt.date(2014,5,1), dt.date(2014,11,1))}
               #'Autumn':(dt.date(2014,9,23), dt.date(2014,12,20)),
               #'Spring':(dt.date(2014,3,21), dt.date(2014,6,20))}
    
    dd = np.asarray([dt.date(2014,d.month,d.day) for d in dates])
    
    season = []

    for indDay in dd: 

        for seas,(season_start, season_end) in seas2get.items():
            if indDay>=season_start and indDay<= season_end:
                season.append(seas)
            else:
                season.append('Winter')

    return season

def setBoxColors(bp):
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    #setp(bp['fliers'][0], color='blue')
    #setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')

    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    #setp(bp['fliers'][2], color='red')
    #setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')


# Total Column to mm
def TC2mm(TC):

    mm = TC*2.989e-22

    return mm



                                                #----------------#
                                                # Define classes #
                                                #----------------#
    


#------------------------------------------------------------------------------------------------------------------------------    
class DateRange(object):
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
