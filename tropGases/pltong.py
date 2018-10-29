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
import numpy as np
import datetime as dt

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


pltDir            =  '/data/iortega/pbin/tropGases/ong/'
file              = 'production-aggregated-data.csv'

cols, indexToName = mf.getColumns(pltDir + file, headerrow=0, delim=',', header=False) #"Category","Oil Production","Gas Production","Water Production","Number of Wells"


SCsd    = np.asarray(cols[0])[1:]
yyyystr = [int(s[1:5]) for s in SCsd]
mnthstr = [int(s[6:8]) for s in SCsd]


OilProd   = np.asarray(cols[1])[1:]; OilProd = np.asarray(OilProd,dtype=np.float32)
GasProd   = np.asarray(cols[2])[1:]; GasProd = np.asarray(GasProd,dtype=np.float32)
WatProd   = np.asarray(cols[3])[1:]; WatProd = np.asarray(WatProd,dtype=np.float32)
WellNum   = np.asarray(cols[4])[1:]; WellNum = np.asarray(WellNum,dtype=np.float32)

date = [dt.date(int(d[1:5]), int(d[6:8]), 15) for d in SCsd]

#ind = np.where(OilProd == np.max(OilProd))[0]
#print date[ind[0]]

#ind1 = np.where(OilProd < 1e6)[0]
#print date[ind1[-1]]
#print OilProd[ind] / OilProd[ind1[-1]]
#print GasProd[ind] / GasProd[ind1[-1]]
#exit()


#----------------------------
# PLOT:
#----------------------------

yearsLc   = YearLocator()
months    = MonthLocator()
DateFmt   = DateFormatter('%Y')

fig, ax = plt.subplots(figsize=(12,6))

ax1 = ax.twinx()
ax2 = ax.twinx()

# Offset the right spine of par2.  The ticks and label have already been
# placed on the right by twinx above.
ax2.spines["right"].set_position(("axes", 1.1))
make_patch_spines_invisible(ax2)
ax2.spines["right"].set_visible(True)

p,  = ax.plot(date, WellNum/1e3, label = 'Number of wells', color='k', linewidth=4.0)
p1, = ax1.plot(date, OilProd/1e6, label = 'Oil Production', color='green', linewidth=4.0)
p2, = ax2.plot(date, GasProd/1e6, label = 'Gas Production', color='red', linewidth=4.0)


ax.tick_params(axis='y',labelsize=16)
ax.tick_params(axis='x', labelsize=16)
#ax.grid(True,which='both', alpha=0.35)  
ax.set_ylabel('Number of wells (x10$^3$)', fontsize=16, color = 'k')
ax.xaxis.set_major_locator(yearsLc)
ax.xaxis.set_major_formatter(DateFmt)
ax.set_xlabel('Year', fontsize=16)
ax.set_xlim(date[0],date[-1])

ax1.tick_params(axis='y',labelsize=16, color='green' , labelcolor="green")
ax1.tick_params(axis='x', labelsize=16)
#ax1.grid(True,which='both', alpha=0.35)  
ax1.set_ylabel('Oild Production (million of barrels)', fontsize=16, color = 'green')
ax1.xaxis.set_major_locator(yearsLc)
ax1.xaxis.set_major_formatter(DateFmt)


ax2.tick_params(labelsize=16, color='red', labelcolor="red")
ax2.tick_params(axis='x', labelsize=16)
ax2.set_ylabel('Gas Production (MCM in Millions)', fontsize=16, color='red')

lines = [p, p1, p2]

#ax.legend(prop={'size':11}, loc='upper center', bbox_to_anchor=(0.92, 1.075), fancybox=True, ncol=1)
ax.legend(lines, [l.get_label() for l in lines], prop={'size':16})

fig.subplots_adjust(left=0.075, bottom=0.1, right=0.85)


fig.autofmt_xdate()

plt.savefig('/data/iortega/pbin/tropGases/fig/ONGProd.pdf', bbox_inches='tight')


#if i == 0: 
#    ax.set_title('Detrended Weigthed VMR & Enhancements\nDotted line represents natural variability', fontsize=18)
#    ax.legend(prop={'size':11}, loc='upper center', bbox_to_anchor=(0.92, 1.075), fancybox=True, ncol=1) 
plt.show(block=False)
user_input = raw_input('Press any key to exit >>> ')
exit() 