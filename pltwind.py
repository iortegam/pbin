import sys
import os
import datetime                                            as dt
import dataOutClass                                        as dc
import numpy                                               as np
from scipy.interpolate import InterpolatedUnivariateSpline as intrpUniSpl
#import matplotlib.pyplot                                   as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import glob
from matplotlib import pyplot as plt
import glob 

from windrose import WindroseAxes
from windrose import WindAxes


def ckDir(dirName,exitFlg=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if exitFlg: sys.exit()
        return False
    else:
        return True

def ckFile(fName,exitFlg=False):
    '''Check if a file exists'''
    if not os.path.isfile(fName):
        print 'File %s does not exist' % (fName)
        if exitFlg: sys.exit()
        return False
    else:
        return True

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

    return np.asarray(wdirList), np.asarray(wspdList), np.asarray(tempList), np.asarray(rhList), np.asarray(datet_weather)

                            #----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#


def main():

    loc = 'fl0'
    weatherDir     = '/data1/ancillary_data/fl0/eol/'
    weatherFileTag = 'v2'
    iyear = 2010
    fyear = 2010


    wdir, wspeed, temp, rh, dtw = weatherout(loc, weatherDir, weatherFileTag, iyear, fyear )
   
    hours = []
    for k in dtw:
        hours.append(k.hour)

    hours = np.asarray(hours)
    
    inds = np.where( (wspeed <= 10.0) & (hours > 8) & (hours < 18) )[0]
    

    fig = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    rect = [0.1, 0.1, 0.8, 0.8]
    ax = WindroseAxes(fig, rect, axisbg='w')
    fig.add_axes(ax)

    ax.bar(wdir[inds], wspeed[inds], normed=True, opening=0.9, edgecolor='white')
    ax.set_legend()

    #fig2 = plt.figure(figsize=(8, 8), dpi=80, facecolor='w', edgecolor='w')
    ax2 = WindAxes.from_ax()
    bins = np.arange(0, 10 , 0.5)
    bins = bins[1:]
    ax2, params = ax2.pdf(wspeed[inds], bins=bins)

    ax3 = WindAxes.from_ax()
    bins = np.arange(0, 360, 15)
    bins = bins[1:]
    ax3, params = ax3.pdf(wdir[inds], bins=bins)

    # fig2,  ax2   = plt.subplots(figsize=(8,6))
    # ax2.scatter(wdir, wspeed, facecolors='red', edgecolors='black', s=35)
    # ax2.grid(True)        

    plt.show(block=False)

    pdfsav = PdfPages('/data/iortega/results/fl0/windrose.pdf')
    pdfsav.savefig(fig,dpi=200)
    pdfsav.close()
    user_input = raw_input('Press any key to exit >>> ')
    sys.exit()    


if __name__ == "__main__":
    main()