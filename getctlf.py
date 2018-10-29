#!/usr/bin/python
##! /usr/local/python-2.7/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        getctlf.py
#
# Purpose:
#
# Get information in the ctl file
#
# Usage (Example):
#      mkListFile.py -i sfit4.ctl -p sfit4_out.ctl -sda -ce -l mlo
#
# Version History:
#       Created, June, 2016  Ivan Ortega (iortega@ucar.edu)
#       Version history stored in git repository
#----------------------------------------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import os
import getopt
import dataOutplts as dc
import getopt
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import myfunctions as mf
import numpy as np
#import sondeclass as rs
from matplotlib.backends.backend_pdf import PdfPages
import dataOutClass as dc
import myfunctions as mf
import fnmatch

def usage():
    ''' Prints to screen standard program usage'''
    print 'getctlf.py -i <file> -p <fileout> -s -c -l -?'
    print '  -i <file>      : Flag to specify the ctl input file'
    print '  -p <fileout>   : Flag to specify the output file printing the sorted ctl file'
    print '  -s <str>       : Flag for showing Sa matrix: d = plot Sa, a = plot Sa and apriori from WACCM (NEEDS -l)'
    print '  -c <str>       : Flag to calculates a second Sa matrix: s = scaled, c = constant; and the scale fractional factor'
    print '  -l <str>       : Flag for location in case a priori is plotted, options for now are (mlo/fl0/tab)'
    print '  -?             : Show all flags'
    print ' Example: getctlf.py -i sfit4.ctl -p sfit4_out.ctl -cs0.3 : To printout sfit4_out.ctl and show the scaled diagonal of Sa with 0.3'

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

def readRefPrf(fname='', parms=''):
        ''' Reads in reference profile, an input file for sfit4 (raytrace) '''
        refPrf = {}
            
        #parms = ['ALTITUDE','PRESSURE','TEMPERATURE']
        
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
   
        #------------------------
        # Convert to numpy arrays
        # and sort based on date
        #------------------------
        for k in refPrf:
            refPrf[k] = np.asarray(refPrf[k])

        return refPrf

def sortDict(DataDict,keyval):
    ''' Sort all values of dictionary based on values of one key'''
    base = DataDict[keyval]
    for k in DataDict:
        DataDict[k] = [y for (x,y) in sorted(zip(base,DataDict[k]))]
    return DataDict

def main(argv):

    savePDF    = False
    PrintFile  = False
    PlotSa     = False
    PlotSa2    = False
    PlotApr    = False
    loc = ''
    pltFile    = 'Sa_ctl.pdf'
    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:p:s::c:l:?')

    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit()

    #-----------------------------
    # Parse command line arguments
    #-----------------------------
    for opt, arg in opts:
        # Check input file flag and path
        if opt == '-i':

            ctlFile = arg

            ckFile(arg,exit=True)

        #
        elif opt == '-p':
            if not arg or arg.startswith('-'):
                usage()
                sys.exit() 
            PrintFile = True
            Fileout = arg

        #
        elif opt == '-s':
            flgs = list(arg)
            for f in flgs:
                if   f.lower() == 'd': PlotSa   = True
                elif f.lower() == 'a': PlotApr  = True
                else: print '{} not an option for -f ... ignored'.format(f)

        elif opt == '-c':
            if not arg or arg.startswith('-'):
                usage()
                sys.exit() 
            flgs = list(arg)
            PlotSa2 = True
            for f in flgs:
                if   f.lower() == 'e': Satype   = 'exp'
                elif arg[0].lower() == 's': 
                    Satype   = 'scaled'
                    SaScl    = float(arg[1:5])
                elif f.lower() == 'c': 
                    Satype   = 'constant'
                    SaScl    = float(arg[1:5])
                else: print '{} not an option for -f ... ignored'.format(f)
            

        elif opt == '-l':
            print arg
            if not arg or arg.startswith('-'):
                usage()
                sys.exit() 
            loc = arg

        # Show all command line flags
        elif opt == '-?':
            usage()
            sys.exit()

        else:
            print 'Unhandled option: ' + opt
            sys.exit()

    if savePDF: pdfsav = PdfPages(pltFile)
   
    #-----------------------------------------------
    #READING THE CTL FILE AND DEFINE PARAMETERS
    #-----------------------------------------------
    (PrimaryGas, ctlinfo) = dc.readCtlF(ctlFile)

    #-----------------------------------------------
    #READING THE CTL FILE IN BLOCKS  
    #-----------------------------------------------

    if PrintFile:
        # Create an order of keys to write spectral db file
        blocknames = ['file.in*', 'file.out*', 'gas.profile*', 'gas.column*', 'fw.', 'rt.*', 'kb.*', 'band*', 'sp*', 'out*']

        FalseVar = []

        with open(Fileout, 'wb') as fopen:
            fopen.write('{}''\n'.format('#-------------------------------------------------------------------------------------' ))
            fopen.write('{}''\n'.format('#---------------------------CTL FILE GENERATED WITH getctl.py-------------------------' ))
            fopen.write('{}''\n''\n'.format('#-------------------------------------------------------------------------------------' ))
            
            print '\n'
            print '#---------------------------CTL FILE GENERATED WITH getctl.py-------------------------'
            print 'Printing file....'
            
            for i, b in enumerate(blocknames):
                dvars=filter(lambda k: fnmatch.fnmatch(k,b),ctlinfo.keys())
                dvars = sorted(dvars, reverse=False)       

                for row in dvars:# for k in sorted(dvars)]):
                    npnt = len(ctlinfo[row])

                    if npnt == 1:
                        if ctlinfo[row][0] == 'F': continue     
                        else:
                            print '{:30}{}{}'.format(row,'=',ctlinfo[row][0])
                            fopen.write('{:40}{}{:}''\n'.format(row,'= ',ctlinfo[row][0]))
                    
                    if (npnt > 1) & (npnt < 10):            
                        fopen.write('{:40}{}'.format(row,'= ')) 
                        strformat = [' {'+str(i)+':<5}' for i in range(npnt)]
                        strformat = ''.join(strformat).lstrip().rstrip()+ '\n'
                        fopen.write(strformat.format( *[p for p in ctlinfo[row]]))

                    if npnt >= 10:
                        fopen.write('{:40}{}''\n'.format(row,'= '))
                        

                        rowsprf = int(np.ceil(npnt/6.0))
                        k = 0
                        for p in range(rowsprf):
                            stuff = ctlinfo[row][k:k+6]
                            strformat = [' {'+str(i)+':<7}' for i in range(len(stuff))]
                            strformat = ''.join(strformat).lstrip().rstrip()+ '\n'
                            print stuff
                            k+= 6 
                            fopen.write(strformat.format( *[p for p in stuff]))
                        fopen.write('\n'.format())
                    
                fopen.write('\n'.format())

            fopen.write('{}''\n'.format('#-------------------------------------------------------------------------------------' ))
            fopen.write('{}''\n'.format('#---------------------------VARIABLES NOT USED-------------------------' ))
            fopen.write('{}''\n''\n'.format('#-------------------------------------------------------------------------------------' ))
            
            print '\n'
            print '#---------------------------VARIABLES NOT USED-------------------------'

            for i, b in enumerate(blocknames):
                dvars=filter(lambda k: fnmatch.fnmatch(k,b),ctlinfo.keys())    

                for row in sorted(dvars, reverse=False):# for k in sorted(dvars)]):
                    npnt = len(ctlinfo[row])

                    if npnt == 1:
                        if ctlinfo[row][0] == 'F':     
                            print '{:30}{}{}'.format(row,'=',ctlinfo[row][0])
                            fopen.write('{:40}{}{:}''\n'.format(row,'= ',ctlinfo[row][0]))

        print 'End printing file....''\n'
 

    stfile = ctlinfo['file.in.stalayers'][0]
    PrimaryGas = ctlinfo['gas.profile.list'][0]     
    Sa = ctlinfo['gas.profile.'+str(PrimaryGas)+'.sigma']
    Sa =np.asarray(Sa)
    scl = ctlinfo['gas.profile.'+str(PrimaryGas)+'.scale']

    #-----------------------------------------------
    #READING THE STATION LAYER FILE
    #-----------------------------------------------
    ckFile(stfile,exit=True)
    #stfile = file(stfile, 'r')
    cols, indexToName = mf.getColumns(stfile, headerrow=2, delim=' ', header=True)
    midpoint = np.asarray(cols['midpnt'][0:-1]).astype(np.float)
    thick    = np.asarray(cols['thick'][0:-1]).astype(np.float)
    level    = np.asarray(cols['level']).astype(np.float)

    print 'Sa in Ctl file:'
    print Sa
    print '\n'


    if PlotSa2:
        #-----------------------------------------------
        #CONSTRUCT A NEW Sa
        #-----------------------------------------------
        
        Sa2 = np.zeros(len(Sa))                             #CONSTANT
        if Satype == 'constant':  Sa2[:] = SaScl

        #if Satype == 'exp':  Sa2 = 2.0*np.exp(-1.0* np.divide(midpoint, 9.0 ) )  #EXP DECREASE
        
        if Satype == 'scaled':
            Sa2 = SaScl/np.sqrt(thick)                            #scaled by the the root of the interlayer thickness
            Sa2 = np.array(Sa2)

            #----------Test For OCS (scale the HIPPO+ACE by the layer thickness)
            #Sa2 = Sa/np.sqrt(thick)
            #----------
        
        print 'Sa2 calculated:'
        print Sa2

        #----------FOR CO
        # # inds = np.where((midpoint >= 20.) & (midpoint < 60.) )[0]
        # # print 'Sa2 calculated:'
        # # Sa2[inds]= Sa2[inds]*0.75 #- MLO
        # # print Sa2

        #inds = np.where((midpoint >= 15.))[0]
        #print 'Sa2 calculated:'
        #Sa2[inds]= Sa2[inds]*0.5 #- MLO
        #print Sa2

        # # #inds = np.where(midpoint >= 50.)[0]
        # # #Sa2[inds]=0.1
        # # #print Sa2
        #----------
    
    if PlotApr:
        #-----------------------------------------------
        #READING WACCM PROFILE
        #-----------------------------------------------
        if PlotApr:
            if not loc: 
                print 'location is needed!'
                loc = raw_input('Input a location for apriori (mlo/tab/fl0) >>> ')
 
            WACCMfile = '/data/Campaign/'+loc.upper()+'/waccm/reference.prf.REFC1.3'
            ckFile(WACCMfile,exit=True)
            prfwaccm = readRefPrf(fname=WACCMfile, parms = ['ALTITUDE','PRESSURE','TEMPERATURE', PrimaryGas.upper()])
            aprprf =  np.asarray(prfwaccm[PrimaryGas.upper()][0]) * scl 
            z = np.asarray(prfwaccm['ALTITUDE'][0])



    if PlotSa:
        #-----------------------------------------------
        #plots
        #-----------------------------------------------
        fig, ax = plt.subplots(1, figsize=(7,8), sharex=True)
        
        ax.plot(Sa, midpoint, color='b', linestyle='-', marker ='o', markersize=6, label='Sa (file)')
        if PlotSa2: 
            ax.plot(Sa2, midpoint, color='g', linestyle='-', marker ='o', markersize=6, label='Sa (new)')
            ax.legend(prop={'size':14})
        ax.xaxis.set_tick_params(which='major',labelsize=12)
        ax.xaxis.set_tick_params(which='minor',labelbottom='off')
        #ax[0].set_xlim((dt.date(iyear,imnth, iday), dt.date(fyear,fmnth,fday)))
        ax.set_ylabel('Altitude [km]', fontsize=16)
        ax.set_xlabel('gas.profile.'+str(PrimaryGas)+'.sigma', fontsize=16, color='b')
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.grid(True)

        for tl in ax.get_xticklabels():
            tl.set_color('b')
        
        if PlotApr:
            ax2 = ax.twiny()
            ax2.plot(aprprf*1e9, z, color='r', linestyle='-', marker ='o', markersize=6, label = 'a priori')
            ax2.xaxis.set_tick_params(which='major',labelsize=12)
            ax2.xaxis.set_tick_params(which='minor',labelbottom='off')
            ax2.tick_params(axis='both', which='major', labelsize=14)
            ax2.set_xlabel('A priori [ppb] - WACCM)', fontsize=16, color='r')


            for tl in ax2.get_xticklabels():
                tl.set_color('r')

        fig.subplots_adjust(left = 0.12, bottom=0.1, top=0.9, right = 0.95)

        if savePDF: 
            pdfsav.savefig(fig,dpi=200)
            pdfsav.close() 
        else:           
            plt.show(block= False)
            user_input = raw_input('Press any key to exit >>> ')
            sys.exit()   
         
    

if __name__ == "__main__":
    main(sys.argv[1:])