#----------------------------------------------------------------------------------------
# Name:
#        input_pltsonde.py --> 

# Purpose:
#        This is the input file for pltsonde.py   --> Plot sonde vs FTS 
#----------------------------------------------------------------------------------------
loc        = 'fl0'
gasName    = 'co'               
    
if loc.lower() == 'fl0':

    if gasName    == 'c2h6':
        ver        = ['sonde/Current_v2_sonde', 'sonde/Current_v2_Ret', 'sonde/Current_v2_ERA_v66', 'sonde/Current_v2_ERA', 'sonde/Current_v2_NCEP', 'sonde/Current_v2_WACCM']
        ctlF       = ['sfit4_v2.ctl','sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl']
        ID         = ['FPH', 'Ret', 'ERA-6', 'ERA-d', 'NCEP', 'WACCM']
        ref        = 'FPH'

    elif gasName    == 'hcn':
        ver        = ['sonde/Current_WP_sonde', 'sonde/Current_WP_Ret', 'sonde/Current_WP_ERA_v66', 'sonde/Current_WP_ERA', 'sonde/Current_WP_NCEP', 'sonde/Current_WP_WACCM']
        ctlF       = ['sfit4.ctl','sfit4.ctl', 'sfit4.ctl', 'sfit4.ctl', 'sfit4.ctl', 'sfit4.ctl']
        ID         = ['FPH', 'Ret', 'ERA-6', 'ERA-d', 'NCEP', 'WACCM']
        ref        = 'FPH'

    elif gasName    == 'co':
        ver        = ['sonde/Current_v3_sonde', 'sonde/Current_v3_Ret', 'sonde/Current_v3_ERA_v66', 'sonde/Current_v3_ERA', 'sonde/Current_v3_NCEP', 'sonde/Current_v3_WACCM']
        ctlF       = ['sfit4_v3.ctl', 'sfit4_v3.ctl', 'sfit4_v3.ctl', 'sfit4_v3.ctl', 'sfit4_v3.ctl', 'sfit4_v3.ctl']
        ID         = ['FPH', 'Ret', 'ERA-6', 'ERA-d',  'NCEP', 'WACCM']
        ref        = 'FPH'
    else:
        print 'Gas Not Found!!'
        exit()

    dois        = ['20100914', '20101105', '20110728', '20120430', '20130117', '20130816', '20140722', '20150629', '20150805', '20160706']
    sdoi        = ['20140722']

else:
    print 'You need to input a location!'
    exit()

#------
# Flags
#------
saveFlg    = True	               # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
errorFlg   = True                  # Flag to process error data
fltrFlg    = True                  # Flag to filter the data
byYrFlg    = False                  # Flag to create plots for each individual year in date range
szaFlg     = True                   # Flag to filter based on min and max SZA
dofFlg     = True                   # Flag to filter based on min DOFs
pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
tcNegFlg   = True                  # Flagsag to filter profiles with negative total columns
tcMMFlg    = False                  # Flag to filter based on min and max total column amount
cnvrgFlg   = False                   # Flag to filter profiles that did not converge
rmsFlg     = True                   # Flag to filter based on max RMS
chiFlg     = False     	            # Flag to filter based on max CHI_2_Y
mnthFlg    = False                  # Flag to filter based on 

mnths      = [6,7,8]                # Months to filter on (these are the months to include data)

if loc.lower() == 'fl0':

    maxRMS     = 2.2                      # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = 0.5
    latFTS     = 40.04
    lonFTS     = -105.24

    pCols = [ [1.6, 3], [3, 5], [5, 7.5], [7.5, 10], [10,13], [13,17], [17,21] ]
    
else:
    print 'You need to input the maxRMS and minDOF!'
    exit() 

minSZA     = 0.0                   # Min SZA for filtering
maxSZA     = 90.0                   # Max SZA for filtering
maxCHI     = 2.0                    # Max CHI_y_2 value
maxTC      = 5.0E25                # Max Total column amount for filtering
minTC      = 0.0                 # Min Total column amount for filtering
sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
sclfctName = 'ppb'                 # Name of scale factor for labeling plots

#----------------------
# Date range to process
#----------------------
iyear      = 2010
imnth      = 1
iday       = 1
fyear      = 2016
fmnth      = 12
fday       = 31

#----------------------
# Directories
#----------------------
#retDir = [ '/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+v+'/' for g,v in zip(gasName,ver)] 
#ctlFile  = ['/data1/ebaumer/'+loc.lower()+'/'+g.lower()+'/'+'x.'+g.lower()+'/'+ ctlF[i] for i, g in enumerate(gasName)]

retDir = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'
#ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'

#------
# Files
#-----

if saveFlg:
    pltDir   = '/data/iortega/pbin/sondeFTIR/fig/'
    pltFile  = pltDir+loc.lower()+'_'+gasName.lower()+'_PrfVs.pdf'
     
