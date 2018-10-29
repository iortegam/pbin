#----------------------------------------------------------
#READ AND PLOT THE  YEARLY HOUSEKEEPING FILES
#----------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages 
import numpy as np
from scipy import stats
import myfunctions as mf

def main(argv):
	#----------------------------------------------------------
	#INPUTS
	#----------------------------------------------------------
	loc        = 'mlo'

	Dir  = '/data/Campaign/'+loc.upper()+'/House_Log_Files/'

	iyear      = 2018

	ipolardays  = dt.datetime(iyear, 3, 1, 0, 0, 0)
	fpolardays  = dt.datetime(iyear, 10, 15, 0, 0, 0)

	savePDF    = False
	pltFile    = '/data/iortega/results/'+loc.lower()+'/fig/'+'HK.pdf'

	if savePDF: pdfsav = PdfPages(pltFile)

	#----------------------------------------------------------
	#                           MAIN
	#----------------------------------------------------------
	dt_hk     = []     	#Date
	T_hk      = []		#temperature
	wdir_hk   = []		#wind direction
	wspeed_hk = []		#wind speed
	P_hk      = []		#pressure
	RH_hk     = []		#rel humidity
	ES_hk     = []		#External solar sensor
	Npnts     = []
	NpntsTot  = []

	f2 = open( Dir + loc.upper()+'_HouseData_'+str(int(iyear))+ '.dat', 'r')
	lines = f2.readlines()
	f2.close()
	count = 0

	for line in lines:
		#if line.startswith('#'): continue
		if not line.startswith("#"):

			pp = line.split()
			yyyy = int(pp[0][0:4])
			mm   = int(pp[0][4:6])
			dd   = int(pp[0][6:8])
			hh = int(pp[1][0:2])
			mi   = int(pp[1][3:5])
			se   = int(pp[1][6:8])

			if (dt.datetime(yyyy, mm, dd, hh, mi, se) >= ipolardays) & (dt.datetime(yyyy, mm, dd, hh, mi, se) <= fpolardays):
				if float(pp[26]) > -5:

					T_hk.append(float(pp[11]))
					wdir_hk.append(float(pp[12]))
					wspeed_hk.append(float(pp[13]))
					P_hk.append(float(pp[14]))
					RH_hk.append(float(pp[15]))
					ES_hk.append(float(pp[26]))
					dt_hk.append(dt.datetime(yyyy, mm, dd, hh, mi, se))
					count+= 1

	T_hk      = np.asarray(T_hk)
	wdir_hk   = np.asarray(wdir_hk)
	wspeed_hk = np.asarray(wspeed_hk)
	P_hk      = np.asarray(P_hk)
	RH_hk     = np.asarray(RH_hk)
	ES_hk     = np.asarray(ES_hk)
	dt_hk     = np.asarray(dt_hk)
	Npnts     = len(dt_hk)
	NpntsTot  = len(lines)

	inds = np.where(ES_hk >= 5.0)[0]

	print 'Total points = ' + str(NpntsTot)
	print 'Total points (polardays) = ' + str(Npnts)
	print 'Number of points above =' + str(len(inds))
	print 'Fraction = '+ str(np.divide(float(len(inds)), float(Npnts)))

	Dateslist = mf.daysList(dt_hk)
	print 'Number of days = '+ str(len(Dateslist))

    #-----------------------------------------------
    #Monthly averages
    #-----------------------------------------------
	#outdata = mf.mnthlyAvg(TH_ncep, dt_ncep, dateAxis=1, meanAxis=0, quad=0)
	#TH_ncep_mth = np.asarray(outdata['mnthlyAvg'])
	#TH_ncep_std_mth = np.asarray(outdata['std'])
	#dt_ncep_mth = np.asarray(outdata['dates'] )

    #-----------------------------------------------
    #plots
    #-----------------------------------------------

	clmap = 'jet'
	cm           = plt.get_cmap(clmap)
	yearsLc      = YearLocator()
	daysLc       = DayLocator()
	months       = MonthLocator()
	DateFmt      = DateFormatter('%m')

	fig, ax = plt.subplots(1, figsize=(10,6), sharex=True)
	ax.plot(dt_hk, ES_hk, color='k', linestyle='None', marker ='', markersize=4)
	ax.scatter(dt_hk, ES_hk, facecolors='b', s=35, color='b', edgecolor='k')
	ax.xaxis.set_major_locator(months)
	ax.xaxis.set_minor_locator(daysLc)
	ax.xaxis.set_major_formatter(DateFmt) 
	ax.xaxis.set_tick_params(which='major',labelsize=12)
	ax.xaxis.set_tick_params(which='minor',labelbottom='off')
	#ax[0].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
	ax.set_ylabel('External Solar Sensor [V]', fontsize=16)
	ax.set_xlabel('Month - ' + str(iyear), fontsize=16)
	ax.tick_params(axis='both', which='major', labelsize=14)
	ax.grid(True)

	fig.autofmt_xdate()
	fig.subplots_adjust(left = 0.1, bottom=0.12, top=0.96, right = 0.95)


	if savePDF: 
	    pdfsav.savefig(fig,dpi=200)
	    pdfsav.close() 
	else:           
		plt.show(block=	False)
		user_input = raw_input('Press any key to exit >>> ')
		sys.exit()   

if __name__ == "__main__":
    main(sys.argv[1:])
