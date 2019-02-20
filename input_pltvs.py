#----------------------------------------------------------------------------------------
# Name:
#        pltts.py --> Plot time series of a series of gases 
#
# Purpose:
#        This is the input file for plte.py   --> Plot error for a serie of target gases defined below
#
#----------------------------------------------------------------------------------------
loc = 'tab'               # Name of station location

if loc == 'tab':

    #gasName    = ['c2h6', 'c2h6']
    #ver        = ['Current_NewSA_NewScale2', 'Current_v2']
    #ctlF       = ['sfit4.ctl', 'sfit4_v2.ctl']
    #ID         = ['v1', 'v2']     
    #ver_vs     = 'v1'
    #g4ratio    = 'H2O'

    gasName    = ['co', 'nh3']#, 'c2h6']
    ver        = ['Current_B3', 'Current_v2', 'Current_v2']
    ctlF       = ['sfit4_v3.ctl', 'sfit4_v2.ctl', 'sfit4_v2.ctl']
    ID         = ['co_B3', 'v2_nh3', 'v2_c2h6']     
    ver_vs     = 'v2_nh3'
    g4ratio    = 'H2O'
    
elif loc == 'mlo':
   
    #gasName    = ['c2h6', 'c2h6']
    #ver        = ['Current_newSA', 'Current_v2']
    #ctlF       = ['sfit4.ctl', 'sfit4_v2.ctl']
    #ID         = ['v1', 'v2']     
    #ver_vs     = 'v1'
    #g4ratio    = 'H2O'

    gasName    = ['co', 'co', 'co']
    ver        = ['Current', 'Current_v2', 'Current_v3']
    ctlF       = ['sfit4_h2o.ctl', 'sfit4_v2.ctl', 'sfit4_v3.ctl']
    ID         = ['v1', 'v2', 'v3']     
    ver_vs     = 'v1'
    g4ratio    = 'H2O'
   
elif loc == 'fl0':
    
    #--------------------
    #Below it is for NO2
    #--------------------
    #gasName    = ['no2', 'no2', 'no2', 'no2', 'no2']                                               # Name of gas
    #ver        = ['Current_v1', 'Current_v2', 'Current_v3', 'Current_v4', 'Current_v5']           # Name of retrieval version to process
    #ctlF       = ['sfit4_v1.ctl', 'sfit4_v2.ctl', 'sfit4_v3.ctl', 'sfit4_v4.ctl', 'sfit4_v5.ctl'] 
    #ID         = ['v1', 'v2', 'v3', 'v4', 'v5']
    #ver_vs     = 'v1'
    #g4ratio    = 'H2O'

    #gasName    = ['no2', 'no2']                                               # Name of gas
    #ver        = ['Current_v2', 'Current_v5']           # Name of retrieval version to process
    #ctlF       = ['sfit4_v2.ctl', 'sfit4_v5.ctl'] 
    #ID         = ['v2', 'v5']
    #ver_vs     = 'v2'
    #g4ratio    = 'H2O'

    #--------------------
    #Below it is for C2H6
    #--------------------
    #gasName    = ['c2h6', 'c2h6']#, 'c2h6']#, 'no2']                                               # Name of gas
    #ver        = ['Current', 'Current_v2']#, 'test']#, 'Current_v5']           # Name of retrieval version to process
    #ctlF       = ['sfit4.ctl', 'sfit4_v2.ctl']#, 'sfit4_v2.ctl' ]#, 'sfit4_v5.ctl'] 
    #ID         = ['v1', 'v2']#, 'test']#, 'v5']
    #ver_vs     = 'v1'
    #g4ratio    = 'H2O'

    #--------------------
    #Below it is for H2CO
    #-------------------- 
    #gasName    = ['h2co', 'h2co', 'h2co', 'h2co', 'h2co']                                               # Name of gas
    #ver        = ['Current_WP_v2', 'Current_WP_v3', 'Current_WP_v4',  'Current_WP_v5', 'Current_WP_v6']           # Name of retrieval version to process
    #ctlF       = ['sfit4_v2.ctl', 'sfit4_v3.ctl', 'sfit4_v4.ctl', 'sfit4_v5.ctl', 'sfit4_v6.ctl'] 
    #ID         = ['v2', 'v3', 'v4', 'v5', 'v6']
    #ver_vs     = 'v3'
    #g4ratio    = 'CH4'

    #--------------------
    #Below it is for CCL4
    #-------------------- 
    #gasName    = ['ccl4', 'ccl4', 'ccl4', 'ccl4']                                               # Name of gas
    #ver        = ['Current_v1', 'Current_v2', 'Current_v3', 'Current_v4']           # Name of retrieval version to process
    #ctlF       = ['sfit4_v1.ctl', 'sfit4_v2.ctl', 'sfit4_v3.ctl', 'sfit4_v4.ctl'] 
    #ID         = ['v1', 'v2', 'v3', 'v4']
    #ver_vs     = 'v1'
    #g4ratio    = 'H2O'

    gasName    = ['hcooh', 'hcooh']                                               # Name of gas
    ver        = ['Current_v1', 'Current_v2']           # Name of retrieval version to process
    ctlF       = ['sfit4.ctl', 'sfit4_v2.ctl'] 
    ID         = ['v1', 'v2']
    ver_vs     = 'v1'
    g4ratio    = 'H2O' 

else:
    print 'You need to input a location!'
    exit() 

#------
# Flags
#------
saveFlg    = False                  # Flag to either save data to pdf file (saveFlg=True) or plot to screen (saveFlg=False)
errorFlg   = True                  # Flag to process error data
fltrFlg    = True                   # Flag to filter the data
byYrFlg    = False                  # Flag to create plots for each individual year in date range
szaFlg     = True                   # Flag to filter based on min and max SZA
dofFlg     = True                   # Flag to filter based on min DOFs
pcNegFlg   = True                   # Flag to filter profiles with negative partial columns
tcNegFlg   = True                  # Flagsag to filter profiles with negative total columns
tcMMFlg    = False                  # Flag to filter based on min and max total column amount
cnvrgFlg   = True                   # Flag to filter profiles that did not converge
rmsFlg     = True                   # Flag to filter based on max RMS
chiFlg     = False                  # Flag to filter based on max CHI_2_Y
mnthFlg    = False                  # Flag to filter based on 

mnths      = [6,7,8]                # Months to filter on (these are the months to include data)

if loc == 'tab':
    maxRMS     = [5.0, 5.0, 5.0, 5.0, 5.0]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.0, 0.0, 0.0, 0.0, 0.0]                    # Min DOFs for filtering

elif loc == 'mlo': 
    maxRMS     = [2.0, 2.0, 2.0, 2.0, 2.0]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0.0, 0.0, 0.0, 0.0, 0.0]                    # Min DOFs for filtering

elif loc == 'fl0':
    #maxRMS     = [2.2, 0.3, 0.8, 5 , 1.0, 0.6, 5]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    #minDOF     = [0.5, 0.9, 0.8, 0 , 2,0, 0.9, 0]

    maxRMS     = [0.8, 1.0, 2.0, 2.0, 2.0]                       # Max Fit RMS to filter data. Data is filtered according to <= maxrms
    minDOF     = [0, 0, 0, 0, 0, 0, 0]

else:
    print 'You need to input the maxRMS and minDOF!'
    exit() 

minSZA     = 0.0                   # Min SZA for filtering
maxSZA     = 90.0                   # Max SZA for filtering
maxCHI     = 2.0                    # Max CHI_y_2 value
maxTC      = 5.0E25                # Max Total column amount for filtering
minTC      = 0.0                 # Min Total column amount for filtering
sclfct     = 1.0E9                  # Scale factor to apply to vmr plots (ppmv=1.0E6, ppbv=1.0E9, etc)
sclfctName = 'ppbv'                 # Name of scale factor for labeling plots

#----------------------
# Date range to process
#----------------------
iyear      = 2017
imnth      = 8
iday       = 19
fyear      = 2017
fmnth      = 8
fday       = 23

#----------------------------
# Partial Columns Bounds [km] 
# [lower bound, upper bound]
# To turn off set to False
#----------------------------
pCols = [[0.0,8.0]]
# Directories
#------------ -
#retDir = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+ver+'/'  
retDir = '/data1/ebaumer/'+loc[0].lower()+'/'+gasName[0].lower()+'/'+ver[0]+'/'  

#------
# Files
#-----
#ctlFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/'+gasName.lower()+'/'+'x.'+gasName.lower()+'/'+ctlF
ctlFile  = '/data1/ebaumer/'+loc.lower()+'/'+gasName[0].lower()+'/'+'x.'+gasName[0].lower()+'/'+ctlF[0]
#pltFile  = '/Volumes/data1/ebaumer/'+loc.lower()+'/' + 'Plots/' + loc + '_' + gasName + '_' + ver + '.pdf'

if saveFlg: 
    pltFile  = '/data/iortega/results/'+loc.lower()+'/fig/'+loc.lower()+'_'+gasName[0].lower()+'_vs.pdf'
else:
    pltFile  = ''

       
