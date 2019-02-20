#! /usr/bin/python2.7
#----------------------------------------------------------
#READ THE  DAYLY TROPOSPHERIC HEIGHT USING YEARLY FILES (ERIC CREATED)
#----------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages 
import numpy as np
from scipy.interpolate import interp2d, interp1d
from scipy import stats
import myfunctions as mf

from sklearn import linear_model, datasets

def main(argv):
	#----------------------------------------------------------
	#INPUTS
	#----------------------------------------------------------
	
	loc        = 'par'

	iyear      = 2011
	fyear      = 2016

	savePDF    = False
	pltFile    = '/data1/ancillary_data/NCEPdata/'+loc.upper()+'_TropHeight.pdf'

	if savePDF: pdfsav = PdfPages(pltFile)

	#----------------------------------------------------------
	#                           MAIN
	#----------------------------------------------------------

	TH_ncep = []
	TP_ncep = []
	dt_ncep = []
	doy_ncep = []

	nyears = int(fyear) - int(iyear) + 1

	if loc.lower() ==  'ahs':
		NCEPhgtDir  = '/data1/projects/ocs/ldr/'
	else:
		NCEPhgtDir  = '/data1/projects/ocs/'+loc.lower()+'/'


	if loc.lower() == 'tab':
		loc2 = 'Thule'
	
	elif loc.lower() == 'mlo':
		loc2 = 'Mauna Loa'

	elif loc.lower() == 'bld':
		loc2 = 'Boulder'

	elif loc.lower() == 'eur':
		loc2 = 'Eureka'

	elif loc.lower() == 'stp':
		loc2 = 'St Petersburg'

	elif loc.lower() == 'jfj':
		loc2 = 'Jungfraujoch'

	elif loc.lower() == 'rkb':
		loc2 = 'Rikubetsu'

	elif loc.lower() == 'tor':
		loc2 = 'Toronto'

	elif loc.lower() == 'tsk':
		loc2 = 'Tsukuba'

	elif loc.lower() == 'wlo':
		loc2 = 'Wollongong'

	elif loc.lower() == 'ldr':
		loc2 = 'Lauder'

	elif loc.lower() == 'ahs':
		loc2 = 'AHTS'

	elif loc.lower() == 'mai':
		loc2 = 'Maido'

	elif loc.lower() == 'std':
		loc2 = 'Stdenis'

	elif loc.lower() == 'zgp':
	    loc2 = 'Zugspitze'

	elif loc.lower() == 'kir':
	    loc2 = 'Kiruna'

	elif loc.lower() == 'iza':
	    loc2 = 'Izana'

	elif loc.lower() == 'par':
	    loc2 = 'Paris'

	elif loc.lower() == 'nya':
	    loc2 = 'Ny Alesund'

	elif loc.lower() == 'bre':
	    loc2 = 'Bremen'

	elif loc.lower() == 'pmb':
	    loc2 = 'Paramaribo'

	elif loc.lower() == 'alt':
	    loc2 = 'Altzomoni'

	else:
	    'Site does not exist!'
	    exit()


	for i in range(nyears):

	    f2 = open( NCEPhgtDir +  'NCEP_trpp/TropHght_'+loc.lower()+'_'+str( int(iyear) +i)+ '.dat', 'r')
	    

	    lines = f2.readlines()
	    f2.close()

	    c1 = []
	    c2 = []
	    c3 = []
	    dt_temp = []

	    count = 0
	    for line in lines:
	        if count >=1:
	            p = line.split()
	            yyyy = int(p[0][0:4])
	            mm   = int(p[0][4:6])
	            dd   = int(p[0][6:8])        
	            c1.append(str(p[0]))
	            c2.append(float(p[1]))
	            c3.append(float(p[1]))
	            dt_temp.append(dt.date(yyyy, mm, dd))
	        count+= 1

	    TH_ncep.extend(np.asarray(c2))
	    TP_ncep.extend(np.asarray(c3))
	    dt_ncep.extend(np.asarray(dt_temp))
	    doy_ncep.extend(mf.toYearFraction(np.asarray(dt_temp)))

	TH_ncep = np.asarray(TH_ncep)
	TP_ncep = np.asarray(TP_ncep)
	dt_ncep = np.asarray(dt_ncep)
	doy_ncep = np.asarray(doy_ncep)

	#----------------------------------------------
	#
	#----------------------------------------------
	
	MeasData  = NCEPhgtDir+'/Dates_'+loc2+'.ascii'

	f2 = open(MeasData, 'r')
	lines = f2.readlines()
	f2.close()

	dateMeas   = []
	hh         = []

	for i, line in enumerate(lines):
		if i >= 5:
			p = line.split()
			yyyy = int(p[0][0:4])
			mm   = int(p[0][4:6])
			dd   = int(p[0][6:8])

			hh.append(p[1])

			dateMeas.append(dt.date(yyyy, mm, dd))
		

	dateMeas    = np.asarray(dateMeas)
	doyMeas     = mf.toYearFraction(dateMeas)
	hh          = np.asarray(hh)

	TP_ncep_interp = interp1d(doy_ncep,TP_ncep, kind='linear')(doyMeas)


	if loc.lower() ==  'ahs':
		outData  = '/data1/projects/ocs/ldr/Dates_'+loc2+'_dNCEP.ascii'
	else: 
		outData  = '/data1/projects/ocs/'+loc.lower()+'/Dates_'+loc2+'_dNCEP.ascii'




	#outData  = '/data1/projects/ocs/'+loc.lower()+'/Dates_Mauna Loa_dNCEP.ascii'

	#----------------------------------------
	#CREATE A FILE WITH DATE AND TIME (TO SEND AND GET BACK THE TROPOPAUSE HEIGHT)
	#----------------------------------------

	j = 0 

	with open(outData, 'wb') as fopen:
		for i, lin in enumerate(lines[0:5]):
			fopen.write('{}'.format(lin))

		for i, dd in enumerate(dateMeas):
			YYYYMMDD = '{0:4d}{1:02d}{2:02d}'.format(dd.year, dd.month, dd.day)
			fopen.write('{0:13}{1:13}{2:.3f}\n'.format(YYYYMMDD, hh[i], TP_ncep_interp[i]))
				

if __name__ == "__main__":
    main(sys.argv[1:])
