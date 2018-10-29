#! /usr/bin/python2.7
#----------------------------------------------------------------------------------------
# Name:
#        pltsondeFTS.py
#
# Purpose:
#
# Plot sonde vs FTS
#
# Input:
#       input_pltsonde.py
#----------------------------------------------------------------------------------------

#---------------
# Import modules
#---------------
import sys
import os
import getopt
import dataOutplts as dc
import time
import datetime as dt
import matplotlib.dates as md
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator
from matplotlib.ticker import FormatStrFormatter, MultipleLocator,AutoMinorLocator,ScalarFormatter
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
import myfunctions as mf
import numpy as np
import classSondeFTS as rs
from matplotlib.backends.backend_pdf import PdfPages #to save multiple pages in 1 pdf...
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as colorbar
import matplotlib.colors as colors
from matplotlib.dates import DateFormatter, MonthLocator, YearLocator, DayLocator, WeekdayLocator, MONDAY
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from matplotlib.ticker import MaxNLocator

import pylab as P

from scipy import interpolate

from math import acos, asin, atan2, cos, hypot
from math import degrees, pi as PI, radians, sin, tan


def usage():
    ''' Prints to screen standard program usage'''
    print 'pltsonde.py -i <inputfile> -?'
    print '  -i <file> : Run pltsonde.py with specified input file'
    print '  -?        : Show all flags'

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

def getidver(version):

    labstr = version.strip().split('/')
    if labstr[1] == 'Current_ERA': idver = 'ERA-d'
    if labstr[1] == 'Current_ERA_v66': idver = 'ERA-6' 
    if labstr[1] == 'Current_WACCM':  idver = 'WACCM'
    if labstr[1] == 'Current_NCEP': idver = 'NCEP-d'

    return idver

def main(argv):

    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:?')

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

            pltInputs = {}

            ckFile(arg,exit=True)

            try:
                execfile(arg, pltInputs)
            except IOError as errmsg:
                print errmsg + ' : ' + arg
                sys.exit()

            if '__builtins__' in pltInputs:
                del pltInputs['__builtins__']               


        # Show all command line flags
        elif opt == '-?':
            usage()
            sys.exit()

        else:
            print 'Unhandled option: ' + opt
            sys.exit()


    sondeFlg = True
    FTSFlg   = True

    if pltInputs['saveFlg']: pdfsav = PdfPages(pltInputs['pltFile'])
    
    #---------------------------------
    #
    #     -- Read RadioSonde --
    #
    #---------------------------------
    if sondeFlg:
        s  = rs.SondeClass(pltInputs['sondeDir'], pltInputs['loc'], iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
        fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'], fleoutFlg=pltInputs['fleoutFlg'])

        s.readfilessonde()
        s.sondePrf()

        sonde = s.sonde

        sondealt  = np.asarray(sonde['Alt_km_a'])
        sondeprf  = np.asarray(sonde['H2Omr_ppmv_a'])


    if FTSFlg:

        #---------------------------------
        #
        #     -- Read FTS --
        #
        #---------------------------------
        
        fts = rs.FTSClass( pltInputs['gasName'], pltInputs['retDir'],  pltInputs['ctlF'], pltInputs['loc'], pltInputs['ver'], iyear=pltInputs['iyear'],imnth=pltInputs['imnth'],iday=pltInputs['iday'],
        fyear=pltInputs['fyear'],fmnth=pltInputs['fmnth'],fday=pltInputs['fday'], incr=1)

        fts.ReadFTS(fltrFlg=pltInputs['fltrFlg'],allGas=False,sclfct=pltInputs['sclfct'],sclname=pltInputs['sclfctName'],mnthFltr=pltInputs["mnths"],mnthFltFlg=pltInputs["mnthFlg"],
               errFlg=pltInputs['errorFlg'],minSZA=pltInputs['minSZA'],maxSZA=pltInputs['maxSZA'],maxRMS=pltInputs['maxRMS'],minTC=pltInputs['minTC'],maxTC=pltInputs['maxTC'],
               minDOF=pltInputs['minDOF'],maxCHI=pltInputs['maxCHI'],dofFlg=pltInputs['dofFlg'],rmsFlg=pltInputs['rmsFlg'],tcFlg=pltInputs['tcNegFlg'],
               pcFlg=pltInputs['pcNegFlg'],szaFlg=pltInputs['szaFlg'],cnvrgFlg=pltInputs['cnvrgFlg'],chiFlg=pltInputs['chiFlg'],tcMMflg=pltInputs['tcMMFlg'], 
               pCols=pltInputs['pCols'], smthFlg=pltInputs['smthFlg'], doiplt=pltInputs['doi'], loc=pltInputs['loc'], latFTS =pltInputs['latFTS'], lonFTS =pltInputs['lonFTS'],
               altCG= pltInputs['altCG'], diffT= pltInputs['diffT'], dfValue= pltInputs['dfValue'])

        # if pltInputs['altCG']:
        
        #     szarad = [radians(s) for s in fts.sza[pltInputs['ver'][0]]]

        #     distFTS = [pltInputs['altCG']*tan(s) for s in   szarad]

        #     latFTS2  = []
        #     lonFTS2  = []

        #     for i, a in enumerate(fts.saa[pltInputs['ver'][0]]):
        #         lat2, lon2 = mf.destination(pltInputs['latFTS'], pltInputs['lonFTS'], distFTS[i], a, radius=6371.008771415)
        #         latFTS2.append(lat2)
        #         lonFTS2.append(lon2)

        #     latFTS2 = np.asarray(latFTS2)
        #     lonFTS2 = np.asarray(lonFTS2)

        # Dist  = [mf.haversine(pltInputs['lonFTS'], pltInputs['latFTS'], lo, la) for (lo, la) in zip(lonFTS2, latFTS2) ]

        # print 'Mean horizontal positions of the center of gravity = {0:.2f} $\pm$ {1:.2f}'.format(np.mean(Dist), np.std(Dist))

    
            #listdates = list(set(listdates))
    
    if FTSFlg & sondeFlg: 


        #-------------------------------------------------------
        # HARD CODED DEFINITIONS
        #-------------------------------------------------------
        maxalt    = 25. #km
        dfValue   = str(pltInputs['dfValue'])
        diffT     = pltInputs['diffT']
        pCols     = pltInputs['pCols']
        loc       = pltInputs['loc']
        fleoutFlg = pltInputs['fleoutFlg']
        timeFlg   = pltInputs['timeFlg']
        if pltInputs['saveFlg']: pltDir    = pltInputs['pltDir']



        Qpercent = 95      

        clr = mf.clrplt()

        #-------------------------------------------------------
        #FINDING COINCIDENT DATES
        #-------------------------------------------------------
        doy_sonde = dc.toYearFraction(sonde['date'])

        #---------------------------
        # Define Delta time,  plot Profiles, and calculate columns/VRM(weighted) for both FTS and sonde
        #---------------------------

        prfmean_day            = {}
        prferr_day             = {}
        prfSTD_day             = {}
        aPrf_day               = {}
        Airmassmean            = {}
        NobsFTS                = {}

        latmean                = {}
        lonmean                = {}

        avkSCFday              = {}
        avkSCFVMRday           = {}
        rPrfMolmean_day        = {}
        aPrfMol_day            = {}
        prferrMol_day          = {}

        sondePrf               = {}
        sondePrfsd             = {}
        sondeairmass           = {}
        sondePrfMol            = {}
        sondePrfMolsd          = {}
        sondePrflat            = {}
        sondePrflon            = {}
        sondedates             = {}

        Prfsonde_interp        = {}
        Prfsonde_sd_interp     = {}
        PrfsondeMol_interp     = {}
        PrfsondeMol_sd_interp  = {}
        sondeairmass_a_interp  = {}
        Prfsondelat_interp     = {}
        Prfsondelon_interp     = {}
        sondealt2              = {}
        sondedt2               = {}

        PrfDiff                = {}
        PrfDiffApr             = {}
        PrfDiffRel             = {}
        PrfDiffAprRel          = {}
        PrfDist                = {}

        FTSvmrP                = {}
        FTSsumP                = {}
        FTSvmrsdP              = {}
        FTSsumsdP              = {}
        FTSdates               = {}

        FTSAprvmrP             = {}
        FTSAprsumP             = {}

        

        sondelatP              = {}
        sondelonP              = {}
        distP                  = {}
        
        for v in pltInputs['ver']:
            listdates = [ dt.date(singDate.year, singDate.month, singDate.day) for singDate in fts.dates[v]  ]
            doy_fts   = dc.toYearFraction(listdates)
            

            intrsctVals = np.intersect1d(doy_fts, doy_sonde, assume_unique=False)
        
            inds1       = np.nonzero( np.in1d( doy_sonde, intrsctVals, assume_unique=False ) )[0]
            inds2       = np.nonzero( np.in1d( doy_fts, intrsctVals, assume_unique=False ) )[0]

            print 'Total Number of coincident dates between FTS and sondes = ' +str(len(intrsctVals))+'\n'
            #---------------------------------
            # Re - Defining variables (FTIR) using Max/Min altitude same sonde dates
            #---------------------------------
            indsalt              = np.where(fts.alt[v] <= maxalt)[0]
            fts.alt[v]           = fts.alt[v][indsalt]
            fts.rPrfVMR[v]       = fts.rPrfVMR[v][:, indsalt];                              fts.rPrfVMR[v][inds2, :]
            fts.aPrfVMR[v]       = fts.aPrfVMR[v][:, indsalt];                              fts.aPrfVMR[v][inds2, :]
            fts.rPrfMol[v]       = fts.rPrfMol[v][:, indsalt];                              fts.rPrfMol[v][inds2, :]
            fts.aPrfMol[v]       = fts.aPrfMol[v][:, indsalt];                              fts.aPrfMol[v][inds2, :]
            fts.avkSCFav[v]      = fts.avkSCFav[v][indsalt[0]:, indsalt[0]:]            
            fts.avkVMRav[v]      = fts.avkVMRav[v][indsalt[0]:, indsalt[0]:]
            fts.avkSCF[v]        = fts.avkSCF[v][:, indsalt[0]:, indsalt[0]: ];   fts.avkSCF[v][inds2, :, : ]
            fts.avkVMR[v]        = fts.avkVMR[v][:, indsalt[0]:, indsalt[0]: ];   fts.avkVMR[v][inds2, :, : ]
            fts.Airmass[v]       = fts.Airmass[v][:, indsalt];                              fts.Airmass[v][inds2, :]

            fts.dates[v]          = fts.dates[v][inds2]
            fts.sza[v]            = fts.sza[v][inds2]
            fts.saa[v]            = fts.saa[v][inds2]
            
            if pltInputs['errorFlg']:
                fts.tot_errvmr[v]    = fts.tot_errvmr[v][:, indsalt];             fts.tot_errvmr[v][inds2, :]
                fts.rand_errvmr[v]   = fts.rand_errvmr[v][:, indsalt];            fts.rand_errvmr[v][inds2, :]
                fts.sys_errvmr[v]    = fts.sys_errvmr[v][:, indsalt];             fts.sys_errvmr[v][inds2, :]


            sondedt            = np.asarray(sonde['dt'])[inds1]
            sondedate          = np.asarray(sonde['date'])[inds1]
            sondealt           = np.asarray(sonde['Alt_km_a'])[inds1] 
            sondeH2OPrf_a      = np.asarray(sonde['H2Omr_ppmv_a'])[inds1]
            sondeH2OPrfsd_a    = np.asarray(sonde['H2Osd_ppmv_a'])[inds1]
            sondealt_a         = np.asarray(sonde['Alt_km_a'])[inds1]
            sondeairmass_a     = np.asarray(sonde['airmass_a'])[inds1]

            if fleoutFlg: 
                sondelat_a     = np.asarray(sonde['latFleout_a'])[inds1]
                sondelon_a     = np.asarray(sonde['lonFleout_a'])[inds1]


            szarad = [radians(s) for s in fts.sza[v]]

            distFTS = [pltInputs['altCG']*tan(s) for s in   szarad]

            latFTS2  = []
            lonFTS2  = []

            for i, a in enumerate(fts.saa[v]):
                lat2, lon2 = mf.destination(pltInputs['latFTS'], pltInputs['lonFTS'], distFTS[i], a, radius=6371.008771415)
                latFTS2.append(lat2)
                lonFTS2.append(lon2)

            latFTS2 = np.asarray(latFTS2)
            lonFTS2 = np.asarray(lonFTS2)

            Dist  = [mf.haversine(pltInputs['lonFTS'], pltInputs['latFTS'], lo, la) for (lo, la) in zip(lonFTS2, latFTS2) ]

            print 'Mean horizontal positions of the center of gravity = {0:.2f} $\pm$ {1:.2f}'.format(np.mean(Dist), np.std(Dist))


            for t, df in enumerate(diffT):

                dfStr = str(df)

                for d, da in enumerate(sondedate):
            
                    deltaT    = sondedt[d] -   fts.dates[v] 
                    deltaTmin =  np.asarray( [k.total_seconds() for k in deltaT])/60.

                    if pltInputs['timeFlg'].lower() == 'int':

                        if t == 0: inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= df])
                        else: inds = np.array([x for x, dff in enumerate(deltaTmin) if( (abs(dff) >= diffT[t-1]) and ( abs(dff) < df))])

                    elif pltInputs['timeFlg'].lower() == 'inc':
                        inds = np.array([x for x, dff in enumerate(deltaTmin) if abs(dff) <= df])

                    if len(inds) >= 1:

                        NobsFTS.setdefault(dfStr+v, []).append(len(inds))

                        prfmean_day.setdefault(dfStr+v, []).append(np.average(fts.rPrfVMR[v][inds],axis=0, weights=fts.tot_errvmr[v][inds]) )
                        prferr_day.setdefault(dfStr+v, []).append(np.sqrt (np.sum(fts.tot_errvmr[v][inds]**2, axis=0))) 
                        prfSTD_day.setdefault(dfStr+v, []).append(np.std(fts.rPrfVMR[v][inds], axis=0))
                        aPrf_day.setdefault(dfStr+v, []).append(np.mean(fts.aPrfVMR[v][inds], axis=0))

                        Airmassmean.setdefault(dfStr+v, []).append(np.mean(fts.Airmass[v][inds], axis=0))
                        avkSCFday.setdefault(dfStr+v, []).append(np.mean(fts.avkSCF[v][inds],axis=0))
                        avkSCFVMRday.setdefault(dfStr+v, []).append(np.mean(fts.avkVMR[v][inds],axis=0))

                        day = np.asarray([dt.date(dd.year,dd.month,dd.day) for dd in fts.dates[v][inds]])
                        uniqueDay = list(set(day))          # Find a list of unique days   
                        FTSdates.setdefault(dfStr+v, []).append(uniqueDay)

                        #--------------------------------------------------------
                        # Get location of FTS based on altitude profile weighted from sone, azimuth (bearing) and SZA
                        #--------------------------------------------------------


                        
                        latmean.setdefault(dfStr+v, []).append(np.mean(latFTS2[inds]))
                        lonmean.setdefault(dfStr+v, []).append(np.mean(lonFTS2[inds]))


                        sondeH2OPrf_i      = sondeH2OPrf_a[d]
                        sondeH2OPrfsd_i    = sondeH2OPrfsd_a[d]
                        sondealt_i          = sondealt_a[d]
                        sondeairmass_i      = sondeairmass_a[d]

                        sondePrf.setdefault(dfStr+v, []).append(sondeH2OPrf_i)
                        sondePrfsd.setdefault(dfStr+v, []).append(sondeH2OPrfsd_i)

                        sondedates.setdefault(dfStr+v, []).append(da)
                        sondedt2.setdefault(dfStr+v, []).append(sondedt[d])
                        sondealt2.setdefault(dfStr+v, []).append(sondealt_i)

                        if fleoutFlg: 
                            sondelat_i            = np.asarray(sondelat_a[d])
                            sondelon_i            = np.asarray(sondelon_a[d])

                    else: continue

                    Prfsonde_interp_i           = interpolate.interp1d(sondealt_i, sondeH2OPrf_i, axis=0, fill_value=sondeH2OPrf_i[0], bounds_error=False, kind='linear')(fts.alt[v] )
                    Prfsonde_sd_interp_i        = interpolate.interp1d(sondealt_i, sondeH2OPrfsd_i, axis=0, fill_value=sondeH2OPrfsd_i[0], bounds_error=False, kind='linear')(fts.alt[v] )
                    sondeairmass_a_interp.setdefault(dfStr+v, []).append(interpolate.interp1d(sondealt_i, sondeairmass_i, axis=0, fill_value=sondeairmass_i[0], bounds_error=False, kind='linear')(fts.alt[v] ))
                    
                    if fleoutFlg:
                        #sondePrflon_interp_i            = interpolate.interp1d(sondealt_a, sondelon_a, axis=0,  bounds_error=False, kind='linear')(alt)
                        #sondePrflat_interp_i            = interpolate.interp1d(sondealt_a, sondelat_a, axis=0,  bounds_error=False, kind='linear')(alt)

                        sondePrflon_interp_i            = interpolate.interp1d(sondealt_i, sondelon_i, axis=0, fill_value=sondelon_i[0], bounds_error=False, kind='linear')(fts.alt[v])
                        sondePrflat_interp_i            = interpolate.interp1d(sondealt_i, sondelat_i, axis=0, fill_value=sondelat_i[0], bounds_error=False, kind='linear')(fts.alt[v])

                        Prfsondelon_interp.setdefault(dfStr+v, []).append(sondePrflon_interp_i)
                        Prfsondelat_interp.setdefault(dfStr+v, []).append(sondePrflat_interp_i)

                    if pltInputs['smthFlg']:
                        #---------------------------------
                        # Smoothing using FTIR AK and apriori
                        #---------------------------------
                        Prfsonde_interp.setdefault(dfStr+v, []).append(np.mean(fts.aPrfVMR[v][inds], axis=0)/1e3 + np.dot(np.mean(fts.avkVMR[v][inds],axis=0), (Prfsonde_interp_i -  np.mean(fts.aPrfVMR[v][inds], axis=0)/1e3)))
                        Prfsonde_sd_interp.setdefault(dfStr+v, []).append(np.mean(fts.aPrfVMR[v][inds], axis=0)/1e3 + np.dot(np.mean(fts.avkVMR[v][inds],axis=0), (Prfsonde_sd_interp_i -  np.mean(fts.aPrfVMR[v][inds], axis=0)/1e3)))

                        
                    else:
                        Prfsonde_interp.setdefault(dfStr+v, []).append(Prfsonde_interp_i)
                        Prfsonde_sd_interp.setdefault(dfStr+v, []).append(Prfsonde_sd_interp_i)
            
                prfmean_day[dfStr+v]           = np.asarray(prfmean_day[dfStr+v])
                prferr_day[dfStr+v]            = np.asarray(prferr_day[dfStr+v])
                prfSTD_day[dfStr+v]            = np.asarray(prfSTD_day[dfStr+v])
                aPrf_day[dfStr+v]              = np.asarray(aPrf_day[dfStr+v])
                Airmassmean[dfStr+v]           = np.asarray(Airmassmean[dfStr+v])
                avkSCFday[dfStr+v]             = np.asarray(avkSCFday[dfStr+v])
                FTSdates[dfStr+v]              = np.asarray(FTSdates[dfStr+v])
                NobsFTS[dfStr+v]               = np.asarray(NobsFTS[dfStr+v])
                latmean[dfStr+v]               = np.asarray(latmean[dfStr+v])
                lonmean[dfStr+v]               = np.asarray(lonmean[dfStr+v])

                Prfsonde_interp[dfStr+v]       = np.asarray(Prfsonde_interp[dfStr+v])
                Prfsonde_sd_interp[dfStr+v]    = np.asarray(Prfsonde_sd_interp[dfStr+v])
                sondeairmass_a_interp[dfStr+v] = np.asarray(sondeairmass_a_interp[dfStr+v])
                sondedates[dfStr+v]            = np.asarray(sondedates[dfStr+v])
                sondealt2[dfStr+v]             = np.asarray(sondealt2[dfStr+v])
                sondePrf[dfStr+v]              = np.asarray(sondePrf[dfStr+v])
                sondePrfsd[dfStr+v]            = np.asarray(sondePrfsd[dfStr+v])
                
                PrfDiff[dfStr+v]               = np.asarray(prfmean_day[dfStr+v] - Prfsonde_interp[dfStr+v])
                PrfDiffApr[dfStr+v]            = np.asarray(aPrf_day[dfStr+v] - Prfsonde_interp[dfStr+v])

                PrfDiffRel[dfStr+v]            = np.true_divide(PrfDiff[dfStr+v], Prfsonde_interp[dfStr+v])*100.
                PrfDiffAprRel[dfStr+v]         = np.true_divide(PrfDiffApr[dfStr+v], Prfsonde_interp[dfStr+v])*100.

                if fleoutFlg:
                    Prfsondelon_interp[dfStr+v]    = np.asarray(Prfsondelon_interp[dfStr+v])
                    Prfsondelat_interp[dfStr+v]    = np.asarray(Prfsondelat_interp[dfStr+v])
                    PrfDist[dfStr+v]  = [ [mf.haversine(np.nanmean(lonmean[dfStr+v]), np.nanmean(latmean[dfStr+v]), lo2, la2 ) for (lo2, la2) in  zip(lo,la)] for (lo, la) in zip(Prfsondelon_interp[dfStr+v], Prfsondelat_interp[dfStr+v]) ]
                    PrfDist[dfStr+v]    = np.asarray(PrfDist[dfStr+v])           
        
       
        #---------------------------------
        #Figure: selected profiles (raw grid)
        #---------------------------------
        if pltInputs['timeFlg'].lower() == 'inc':
           
            if pltInputs['doi']:
                doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in pltInputs['doi']]
                doidt = np.asarray(doidt) 

                for v in pltInputs['ver']:

                    fig, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
                    ndoi = 0

                    
                    for d, da in enumerate(sondedates[dfValue+v]):

                        deltadoi = doidt - da

                        if ndoi < len(doidt):

                            if deltadoi[ndoi] == dt.timedelta(0):

                                if loc.lower() == 'mlo': indsFPH =  np.where(sondealt2[dfValue+v][d] >= 3.0)[0]
                                elif loc.lower() == 'fl0': indsFPH =  np.where(sondealt2[dfValue+v][d] >= 1.0)[0]


                                if ndoi < len(doidt):

                                    if ndoi<=4:

                                        if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=16)


                                        ax0[0, ndoi].plot(prfmean_day[dfValue+v][d]/1e3,fts.alt[v], color='b',  linewidth=2.0, label='FTIR', zorder=5)
                                        ax0[0, ndoi].scatter(prfmean_day[dfValue+v][d]/1e3,fts.alt[v],facecolors='white', s=35, color='b', zorder=6)
                                        ax0[0, ndoi].fill_betweenx(fts.alt[v],prfmean_day[dfValue+v][d]/1e3-prferr_day[dfValue+v][d]/1e3,prfmean_day[dfValue+v][d]/1e3+prferr_day[dfValue+v][d]/1e3,alpha=0.25,color='blue')

                                        ax0[0, ndoi].set_ylim(1, 15)
                                        if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                        if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                           
                                        ax0[0, ndoi].grid(True,which='both', alpha=0.35)
                                        ax0[0, ndoi].tick_params(which='both',labelsize=16)
                                        ax0[0, ndoi].set_title(da, fontsize=16)   

                                        ax0[0, ndoi].plot(aPrf_day[dfValue+v][d]/1e3, fts.alt[v], '-', color='gray', label='a priori', linewidth=2.0, zorder=3)
                                        ax0[0, ndoi].scatter(aPrf_day[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='gray', zorder=4)
                                           
                                        ax0[0, ndoi].plot(sondePrf[dfValue+v][d][indsFPH]/1e3, sondealt2[dfValue+v][d][indsFPH],color='k', label='FPH', linewidth=2.0, zorder=1)
                                        #ax0[0, ndoi].scatter(sondePrf[dfValue+v][d]/1e3, sondealt2[dfValue+v][d],facecolors='white', s=35, color='k', zorder=2)
                                        #ax0[0, ndoi].fill_betweenx(sondealt2[dfValue+v][d][indsFPH],sondePrf[dfValue+v][d][indsFPH]/1e3-sondePrf[dfValue+v][d][indsFPH]/1e3*0.05,sondePrf[dfValue+v][d][indsFPH]/1e3+sondePrf[dfValue+v][d][indsFPH]/1e3*0.05,alpha=0.25,color='k')
                                        ax0[0, ndoi].fill_betweenx(sondealt2[dfValue+v][d][indsFPH],sondePrf[dfValue+v][d][indsFPH]/1e3-sondePrfsd[dfValue+v][d][indsFPH]/1e3,sondePrf[dfValue+v][d][indsFPH]/1e3+sondePrfsd[dfValue+v][d][indsFPH]/1e3,alpha=0.25,color='k')


                                        ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue+v][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=16)

                                        ax0[0, ndoi].set_xlim(xmin=0)
                                        if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                                        if aPrf_day[dfValue+v][d][-1]/1e3 <= 3: ax0[0, ndoi].xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
                                        else: ax0[0, ndoi].xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

                                    if ndoi>=5:

                                        ax0[1, (ndoi-5)].plot(prfmean_day[dfValue+v][d]/1e3,fts.alt[v], color='b',  linewidth=2.0, label='Mean (Retrieved)', zorder=5)
                                        ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue+v][d]/1e3,fts.alt[v],facecolors='white', s=35, color='b', zorder=6)
                                        ax0[1, (ndoi-5)].fill_betweenx(fts.alt[v],prfmean_day[dfValue+v][d]/1e3-prferr_day[dfValue+v][d]/1e3,prfmean_day[dfValue+v][d]/1e3+prferr_day[dfValue+v][d]/1e3,alpha=0.25,color='blue')

                                        if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=16)

                                        if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                        if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)
     
                                        if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=16)
                                        ax0[1, (ndoi-5)].grid(True,which='both', alpha=0.35)
                                        ax0[1, (ndoi-5)].tick_params(which='both',labelsize=16)
                                        ax0[1, (ndoi-5)].set_title(da, fontsize=16)

                                        ax0[1, (ndoi-5)].plot(aPrf_day[dfValue+v][d]/1e3, fts.alt[v], '-', color='gray', label='a priori', linewidth=2.0, zorder=3)
                                        ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='gray', zorder=2)
                                           
                                        #ax0[1, (ndoi-5)].plot(sondePrf[dfValue+v][d]/1e3, sondealt2[dfValue+v][d],color='k', linewidth=2.0, zorder=1)
                                        ax0[1, (ndoi-5)].plot(sondePrf[dfValue+v][d][indsFPH]/1e3, sondealt2[dfValue+v][d][indsFPH],color='k', label='FPH', linewidth=2.0, zorder=1)
                                        ##ax0[1, (ndoi-5)].scatter(sondePrf[dfValue+v][d]/1e3, sondealt2[dfValue+v][d],facecolors='white', s=35, color='k',zorder=2)
                                        #ax0[1, (ndoi-5)].fill_betweenx(sondealt2[dfValue+v][d],sondePrf[dfValue+v][d]/1e3-sondePrf[dfValue+v][d]/1e3*0.05,sondePrf[dfValue+v][d]/1e3+sondePrf[dfValue+v][d]/1e3*0.05,alpha=0.25,color='k')
                                        #ax0[1, (ndoi-5)].fill_betweenx(sondealt2[dfValue+v][d][indsFPH],sondePrf[dfValue+v][d][indsFPH]/1e3-sondePrf[dfValue+v][d][indsFPH]/1e3*0.05,sondePrf[dfValue+v][d][indsFPH]/1e3+sondePrf[dfValue+v][d][indsFPH]/1e3*0.05,alpha=0.25,color='k')
                                        ax0[1, (ndoi-5)].fill_betweenx(sondealt2[dfValue+v][d][indsFPH],sondePrf[dfValue+v][d][indsFPH]/1e3-sondePrfsd[dfValue+v][d][indsFPH]/1e3,sondePrf[dfValue+v][d][indsFPH]/1e3+sondePrfsd[dfValue+v][d][indsFPH]/1e3,alpha=0.25,color='k')

                                        ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue+v][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=16)

                                        ax0[1, (ndoi-5)].set_xlim(xmin=0)

                                        if aPrf_day[dfValue+v][d][-1]/1e3 <= 3: ax0[1, (ndoi-5)].xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
                                        else: ax0[1, (ndoi-5)].xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

                                    ndoi+=1

                    fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 

                    idver = getidver(v)
                                    
                    if pltInputs['saveFlg']: 
                        pdfsav.savefig(fig,dpi=200)
                        if pltInputs['smthFlg']: 
                            plt.savefig(pltDir+'Selected_Prf_'+loc.upper()+'_'+idver+'_smooth.pdf', bbox_inches='tight')
                        else:
                            plt.savefig(pltDir+'Selected_Prf_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
                    else:       
                        plt.show(block=False)
                        user_input = raw_input('Press any key to exit >>> ')
                        sys.exit() 


            #---------------------------------
            #Figure: selected profiles (smoothed or not and same grid)
            #---------------------------------           
            if pltInputs['doi']:
                doidt = [dt.date(int(d[0:4]), int(d[4:6]), int(d[6:8])) for d in pltInputs['doi']]
                doidt = np.asarray(doidt) 

                dfValue = str(pltInputs['dfValue'])

                for v in pltInputs['ver']:

                    fig, ax0 = plt.subplots(2, 5, figsize=(15,10), sharey=True, sharex=False)
                    ndoi = 0

                    for d, da in enumerate(sondedates[dfValue+v]):

                        deltadoi = doidt - da

                        if ndoi < len(doidt):

                            if deltadoi[ndoi] == dt.timedelta(0):

                                if ndoi < len(doidt):

                                    if ndoi<=4:

                                        if ndoi ==0: ax0[0, ndoi].set_ylabel('Altitude [km]', fontsize=14)

                                        ax0[0, ndoi].plot(prfmean_day[dfValue+v][d]/1e3,fts.alt[v], color='b',  linewidth=2.0, label='FTIR', zorder=5)
                                        ax0[0, ndoi].scatter(prfmean_day[dfValue+v][d]/1e3,fts.alt[v],facecolors='white', s=35, color='b', zorder=6)
                                        ax0[0, ndoi].fill_betweenx(fts.alt[v],prfmean_day[dfValue+v][d]/1e3-prferr_day[dfValue+v][d]/1e3,prfmean_day[dfValue+v][d]/1e3+prferr_day[dfValue+v][d]/1e3,alpha=0.25,color='blue')

                                        ax0[0, ndoi].set_ylim(1, 15)
                                        if loc.lower() == 'fl0': ax0[0, ndoi].set_ylim(1, 15)
                                        if loc.lower() == 'mlo': ax0[0, ndoi].set_ylim(3, 15)
                                           
                                        ax0[0, ndoi].grid(True,which='both', alpha=0.35)
                                        ax0[0, ndoi].tick_params(which='both',labelsize=14)
                                        ax0[0, ndoi].set_title(da, fontsize=14)   

                                        ax0[0, ndoi].plot(aPrf_day[dfValue+v][d]/1e3, fts.alt[v], '-', color='gray', label='a priori', linewidth=2.0, zorder=3)
                                        ax0[0, ndoi].scatter(aPrf_day[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='gray', zorder=4)
                                           
                                        ax0[0, ndoi].plot(Prfsonde_interp[dfValue+v][d]/1e3, fts.alt[v],color='k', label='CFH', linewidth=2.0, zorder=1)
                                        ax0[0, ndoi].scatter(Prfsonde_interp[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='k', zorder=2)

                                        ax0[0, ndoi].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue+v][d]), va='center',transform=ax0[0, ndoi].transAxes,fontsize=14)

                                        ax0[0, ndoi].set_xlim(xmin=0)
                                        if ndoi == 0: ax0[0, ndoi].legend(prop={'size':11})

                                        if prfmean_day[dfValue+v][d][-1]/1e3 <= 3: ax0[0, ndoi].xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
                                        else: ax0[0, ndoi].xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))


                                    if ndoi>=5:

                                        ax0[1, (ndoi-5)].plot(prfmean_day[dfValue+v][d]/1e3,fts.alt[v], color='b',  linewidth=2.0, label='Mean (Retrieved)', zorder=5)
                                        ax0[1, (ndoi-5)].scatter(prfmean_day[dfValue+v][d]/1e3,fts.alt[v],facecolors='white', s=35, color='b', zorder=6)
                                        ax0[1, (ndoi-5)].fill_betweenx(fts.alt[v],prfmean_day[dfValue+v][d]/1e3-prferr_day[dfValue+v][d]/1e3,prfmean_day[dfValue+v][d]/1e3+prferr_day[dfValue+v][d]/1e3,alpha=0.25,color='blue')

                                        if ndoi ==5: ax0[1, (ndoi-5)].set_ylabel('Altitude [km]', fontsize=14)

                                        if loc.lower() == 'fl0': ax0[1, (ndoi-5)].set_ylim(1, 15)
                                        if loc.lower() == 'mlo': ax0[1, (ndoi-5)].set_ylim(3, 15)
        
                                        if ndoi ==7 :ax0[1, (ndoi-5)].set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                                        ax0[1, (ndoi-5)].grid(True,which='both', alpha=0.35)
                                        ax0[1, (ndoi-5)].tick_params(which='both',labelsize=14)
                                        ax0[1, (ndoi-5)].set_title(da, fontsize=14)

                                        ax0[1, (ndoi-5)].plot(aPrf_day[dfValue+v][d]/1e3, fts.alt[v], '-', color='gray', label='a priori', linewidth=2.0, zorder=3)
                                        ax0[1, (ndoi-5)].scatter(aPrf_day[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='gray', zorder=4)
                                           
                                        ax0[1, (ndoi-5)].plot(Prfsonde_interp[dfValue+v][d]/1e3, fts.alt[v],color='k', linewidth=2.0, zorder=1)
                                        ax0[1, (ndoi-5)].scatter(Prfsonde_interp[dfValue+v][d]/1e3, fts.alt[v],facecolors='white', s=35, color='k', zorder=2)

                                        ax0[1, (ndoi-5)].text(0.05, 0.95, 'N = {}'.format(NobsFTS[dfValue+v][d]), va='center',transform=ax0[1, (ndoi-5)].transAxes,fontsize=14)

                                        ax0[1, (ndoi-5)].set_xlim(xmin=0)

                                        if prfmean_day[dfValue+v][d][-1]/1e3 <= 3: ax0[1, (ndoi-5)].xaxis.set_major_locator(MaxNLocator(nbins=4, integer=False))
                                        else: ax0[1, (ndoi-5)].xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))

                                    ndoi+=1

                    fig.subplots_adjust(bottom=0.075,top=0.95, left=0.05, right=0.95) 

                    idver = getidver(v)
                                    
                    if pltInputs['saveFlg']: 
                        pdfsav.savefig(fig,dpi=200)
                        if pltInputs['smthFlg']: 
                            plt.savefig(pltDir+'Selected_Prf_Interp_'+loc.upper()+'_'+idver+'_smooth.pdf', bbox_inches='tight')
                        else:
                            plt.savefig(pltDir+'Selected_Prf_Interp_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
                    else:       
                        plt.show(block=False)

        #------------------------------------
        #HISTOGRAM OF RESIDUALS BEFORE FILTERING
        #------------------------------------
        for v in pltInputs['ver']:
            
            fig = plt.figure(figsize=(10,7.5))

            for p, pcol in enumerate(pCols[0:-1]):

                gs1 = gridspec.GridSpec(1, 3)

                if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
                if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                
                if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
                if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

                ax1 = plt.subplot(gs1[0:1, :])

                indsH = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <=pcol[1])  )[0]
                x = PrfDiff[dfValue+v][:,indsH]
                x = x[~np.isnan(x)]
            
                n, bins, patches = ax1.hist(x/1e3, bins=20, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
                P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
                mu = np.nanmean(x/1e3)
                me = np.nanmedian(x/1e3)
                sigma = np.nanstd(x/1e3)
                #add a line showing the expected distribution
                y = P.normpdf( bins, mu, sigma)
                l = ax1.plot(bins, y, 'k--', linewidth=3)

                ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent), color='k', linestyle='--')
                ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent)*-1., color='k', linestyle='--')

                ax1.axvline(x=mu, color='blue', linestyle='--')
                ax1.axvline(x=me, color='green', linestyle='--')

                ax1.axvline(x=mu + sigma, color='r', linestyle='--')
                ax1.axvline(x=mu - sigma, color='r', linestyle='--')


                ax1.grid(True,which='both', alpha=0.35)
                ax1.tick_params(which='both',labelsize=12)#, labelbottom='off')
                #ax1.set_ylabel('Probability', fontsize=12)
                #ax1.set_title(str(pcol[0])+' - '+str(pcol[1])+' km',horizontalalignment='center', verticalalignment='baseline', fontsize=14)
                #ax1.set_xlim(-1, 1)
                ax1.text(0.95, 0.92, str(pcol[0])+'-'+str(pcol[1])+' km', va='center',transform=ax1.transAxes,fontsize=13, ha='right')

                if (p == 4) or (p == 5): 
                    ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                    ax1.tick_params(labelbottom='on')
                if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)
               
                # print '\nAltitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
                # #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
                # print 'Mean Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, sigma)
                # print 'Median Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, sigma)
                # print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100.)
                # print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100.)

                # print 'Percentile 95 in Delta T {} min : {}'.format(dfValue+v, np.percentile(np.abs(x/1e3), Qpercent))

            #plt.suptitle(v, fontsize=16)

            idver = getidver(v)

            if pltInputs['saveFlg']: 
                pdfsav.savefig(fig,dpi=200)
                if pltInputs['smthFlg']: 
                    plt.savefig(pltDir+'Hist1_'+loc.upper()+'_'+idver+'_smooth.pdf', bbox_inches='tight')
                else:
                    plt.savefig(pltDir+'Hist1_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
            else:       
                plt.show(block=False)

        
        #------------------------------------
        #FILTERING OUTLIERS BASED ON PERCENTILE
        #------------------------------------
        Prct        = {}

        for v in pltInputs['ver']:

            for df in diffT:

                df = str(df)

                for pn, pcol in enumerate(pCols):

                    indsH = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <=pcol[1])  )[0]

                    if df == dfValue: 
                        print 'Percentile 95 in Delta T {} min and Pcol {} : {}'.format(df, pcol, np.percentile(np.abs(PrfDiff[df+v][:,indsH]), Qpercent))
                        Prct.setdefault(pcol[0], []).append(np.percentile(np.abs(PrfDiff[df+v][:,indsH]), Qpercent))

                    if loc.lower() == 'fl0': iBad = np.where( (np.abs(PrfDiff[df+v][:,indsH]) >= np.percentile(np.abs(PrfDiff[df+v][:,indsH]), Qpercent)) | (PrfDist[df+v][:,indsH] > 300.0) )
                    if loc.lower() == 'mlo': iBad = np.where( (np.abs(PrfDiff[df+v][:,indsH]) >= np.percentile(np.abs(PrfDiff[df+v][:,indsH]), Qpercent)))

                    iBad = np.asarray(iBad)
                    
                    if  df == 30: print 'Percent of bad points in Delta T = {} minutes: {}'.format(df+v, float(iBad.shape[1])/float((PrfDiffRel[df+v][:,indsH].shape[0] * PrfDiffRel[df+v][:,indsH].shape[1]))*100.)
                    #print PrfDiffRel[df]
                    prfmean_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(prfmean_day[df], indsBad, axis=0)
                    
                    prferr_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(prferr_day[df], indsBad, axis=0)
                    prfSTD_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(prfSTD_day[df], indsBad, axis=0)
                    aPrf_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                        = np.delete(aPrf_day[df], indsBad, axis=0)
                    Airmassmean[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                     = np.delete(Airmassmean[df], indsBad, axis=0)
                    avkSCFday[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                       = np.delete(avkSCFday[df], indsBad, axis=0)
                    #rPrfMolmean_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                 = np.delete(rPrfMolmean_day[df], indsBad, axis=0)
                    #prferrMol_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(prferrMol_day[df], indsBad, axis=0)
                    #aPrfMol_day[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                     = np.delete(aPrfMol_day[df], indsBad, axis=0)

                    Prfsonde_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                 = np.delete(Prfsonde_interp[df], indsBad, axis=0)
                    Prfsonde_sd_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsonde_sd_interp[df], indsBad, axis=0)
                    #PrfsondeMol_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(PrfsondeMol_interp[df], indsBad, axis=0)
                    #PrfsondeMol_sd_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#           = np.delete(PrfsondeMol_sd_interp[df], indsBad, axis=0)
                    sondeairmass_a_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#           = np.delete(sondeairmass_a_interp[df], indsBad, axis=0)

                    PrfDiff[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                         = np.delete(PrfDiff[df], indsBad, axis=0)
                    PrfDiffApr[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(PrfDiffApr[df], indsBad, axis=0)
                    PrfDiffRel[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                      = np.delete(PrfDiffRel[df], indsBad, axis=0)

                    PrfDiffAprRel[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                   = np.delete(PrfDiffAprRel[df], indsBad, axis=0)

                    if fleoutFlg:
                        Prfsondelon_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsondelon_interp[df], indsBad, axis=0)
                        Prfsondelat_interp[df+v][iBad[0], indsH[iBad[1]]] = np.nan#              = np.delete(Prfsondelat_interp[df], indsBad, axis=0)
                        PrfDist[df+v][iBad[0], indsH[iBad[1]]] = np.nan#                         = np.delete(PrfDist[df], indsBad, axis=0)

                    #if fleoutFlg:
                        #sondelatP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelat_interp[df][:,indsH]), axis=1)
                        #sondelonP[str(df)+'_'+str(pn)]    = np.ma.mean(np.ma.masked_invalid(Prfsondelon_interp[df][:,indsH]), axis=1)
                        #distP[str(df)+'_'+str(pn)]        = np.ma.mean(np.ma.masked_invalid(PrfDist[df][:,indsH]), axis=1)
        
        #----------------------------------------
        #
        #----------------------------------------
        # ind = np.arange(len(pltInputs['ver']))

        # fig, ax = plt.subplots(figsize=(9, 6))

        # xlabel = [getidver(vi) for vi in pltInputs['ver']]

        # ax.bar(ind-0.2, Prct[dfValue], 0.2, align='center', color = 'r', ecolor = 'k', label='ERA-I-6')
    

        # ax.legend(prop={'size':11}, loc=2)
        
        # ax.yaxis.grid(True, alpha=0.5)
        # ax.set_ylabel('Difference [x10$^3$] (Percentile 95%)', fontsize=14)
        # ax.set_xlabel('Layer [km]', fontsize=14)    
        # ax.set_xticks(ind)
        # ax.set_xticklabels(xlabel, rotation=45)
        # ax.tick_params(labelsize=14)
        # #ax.set_xlabel('Site')
        # #ax.axvline(0, color='black', lw=1)

        # #plt.gca().invert_yaxis()

        # plt.show(block=False)

        #-------------------
        #Bar Plot with Bias
        #-------------------



        for df in diffT:

            df = str(df)

            if df == dfValue:

                fig, ax  = plt.subplots(figsize=(9,6), sharex=True)

                xlabel = [getidver(vi) for vi in pltInputs['ver']]

                ind = np.arange(len(pltInputs['ver'] ))                         
                width = 1. / (len(pCols[0:-2]) + 1)

                bar_groups = []
                labels     = []
                for pn, pcol in enumerate(pCols[0:-2]):

                    Prct_i = np.asarray(Prct[pcol[0]], dtype=np.float32)
                    
                    labels.append(str(pcol[0])+'-'+str(pcol[1])+' km')                    
                    bars = ax.bar(ind+pn*width, Prct_i/1e3, width, color=clr[pn % len(clr)])
                    bar_groups.append(bars)

                ax.set_xlim(-width,len(ind)+width*0.4)
                ax.set_ylabel('Difference in VMR [x10$^3$ ppm] - Percentile 95%', fontsize=14)
                ax.set_xlabel('A priori source', fontsize=14)
                #ax.set_ylim(0,1.5)
                ax.legend([b[0] for b in bar_groups], labels, fontsize=11, ncol=len(labels), loc='upper center', bbox_to_anchor=(0.5, 1.05))
                #ax.text(0.01, 0.95, '(a)', va='center', ha='left', transform=ax.transAxes,fontsize=14)

                ax.set_xticks(ind+width*2)
                ax.tick_params(which='both',labelsize=14)
                #ax.set_ylim(-0.5,0.3)
                ax.set_xticklabels(xlabel, rotation=0)

                fig.subplots_adjust(bottom=0.125,top=0.975, left=0.15, right=0.95)
                
                if pltInputs['saveFlg']: 
                    pdfsav.savefig(fig,dpi=200)
                    if pltInputs['smthFlg']: 
                        plt.savefig(pltDir+'PercentileBars_'+loc.upper()+'_'+df+'_smooth.pdf', bbox_inches='tight')
                    else:
                        plt.savefig(pltDir+'PercentileBars_'+loc.upper()+'_'+df+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)

        #user_input = raw_input('Press any key to exit >>> ')
        #sys.exit()           # Exit program    


        #------------------------------------
        #HISTOGRAM OF RESIDUALS AFTER FILTERING
        #------------------------------------
        for v in pltInputs['ver']:
            
            fig = plt.figure(figsize=(10,7.5))

            for p, pcol in enumerate(pCols[0:-1]):

                gs1 = gridspec.GridSpec(1, 3)

                if p == 0: gs1.update(left=0.075, right=0.5, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                if p == 2: gs1.update(left=0.075, right=0.5, top = 0.64, bottom=0.38, wspace=0.05, hspace=0.08)
                if p == 4: gs1.update(left=0.075, right=0.5, top = 0.33, bottom=0.07, wspace=0.05, hspace=0.08)
                
                if p == 1: gs1.update(left=0.555, right=0.98, top=0.95, bottom=0.69,  wspace=0.05,  hspace=0.08)
                if p == 3: gs1.update(left=0.555, right=0.98, top=0.64, bottom=0.38, wspace=0.05,   hspace=0.08)
                if p == 5: gs1.update(left=0.555, right=0.98, top=0.33, bottom=0.07, wspace=0.05,   hspace=0.08)

                ax1 = plt.subplot(gs1[0:1, :])

                indsH = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <=pcol[1])  )[0]
                x = PrfDiff[dfValue+v][:,indsH]
                x = x[~np.isnan(x)]
            
                n, bins, patches = ax1.hist(x/1e3, bins=20, normed=True, alpha=0.8 ,histtype='step', fill=True, color='k')  #color='green'
                P.setp(patches, 'facecolor', 'g', 'alpha', 0.8)
                mu = np.nanmean(x/1e3)
                me = np.nanmedian(x/1e3)
                sigma = np.nanstd(x/1e3)
                #add a line showing the expected distribution
                y = P.normpdf( bins, mu, sigma)
                l = ax1.plot(bins, y, 'k--', linewidth=3)

                #ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent), color='k', linestyle='--')
                #ax1.axvline(x=np.percentile(np.abs(x/1e3), Qpercent)*-1., color='k', linestyle='--')

                ax1.axvline(x=mu, color='blue', linestyle='--')
                ax1.axvline(x=me, color='green', linestyle='--')

                ax1.axvline(x=mu + sigma, color='r', linestyle='--')
                ax1.axvline(x=mu - sigma, color='r', linestyle='--')


                ax1.grid(True,which='both')
                ax1.tick_params(which='both',labelsize=12)#, labelbottom='off')
                #ax1.set_ylabel('Probability', fontsize=12)
                #ax1.set_title(str(pcol[0])+' - '+str(pcol[1])+' km',horizontalalignment='center', verticalalignment='baseline', fontsize=14)
                #ax1.set_xlim(-1, 1)
                ax1.text(0.95, 0.92, str(pcol[0])+'-'+str(pcol[1])+' km', va='center',transform=ax1.transAxes,fontsize=13, ha='right')

                if (p == 4) or (p == 5): 
                    ax1.set_xlabel('VMR [x10$^3$ ppm]', fontsize=14)
                    ax1.tick_params(labelbottom='on')
                if (p == 0) or (p == 2) or (p == 4): ax1.set_ylabel('Probability', fontsize=14)
               
                print '\nAltitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
                #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
                print 'Mean Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, sigma)
                print 'Median Histogram (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, sigma)
                print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100.)
                print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100., sigma/np.nanmean(Prfsonde_interp[dfValue+v][:, indsH]/1e3) * 100.)

                print 'Percentile 95 in Delta T {} min : {}'.format(dfValue+v, np.percentile(np.abs(x/1e3), Qpercent))

            #plt.suptitle(v, fontsize=16)

            idver = getidver(v)

            if pltInputs['saveFlg']: 
                pdfsav.savefig(fig,dpi=200)
                if pltInputs['smthFlg']: 
                    plt.savefig(pltDir+'Hist2_'+loc.upper()+'_'+idver+'_smooth.pdf', bbox_inches='tight')
                else:
                    plt.savefig(pltDir+'Hist2_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
            else:       
                plt.show(block=False)


        #------------------------------------
        #Fig: Stats with Coincident Interval
        #------------------------------------
        if pltInputs['timeFlg'].lower() == 'inc':
            
            for v in pltInputs['ver']:

                fig1, ax1  = plt.subplots(2, figsize=(8,7), sharex=True)

                for pn, pcol in enumerate(pCols):

                    stdDF        = []
                    Npnts        = []
                    Npnts2       = []
                    biasDF       = []
                    rmseDF       = []

                    for df in diffT:

                        df = str(df)

                        inds = np.where( (fts.alt[v] > pcol[0]) & (fts.alt[v] <=pcol[1])  )[0]

                        x = PrfDiffRel[df+v][:, inds]
                        x = np.asarray(x)
                        x = x.flatten()

                        biasDF.append(np.nansum( x) / len(x))

                        ftsPrf = [i[inds] for i in prfmean_day[df+v]]
                        sonPrf = [i[inds] for i in Prfsonde_interp[df+v]]

                        st = prfSTD_day[df+v][:, inds]
                        st = np.asarray(st)
                        st = st.flatten()

                        stdDF.append(float(np.nanmean(st)/np.nanmean(ftsPrf)))


                        ftsPrf = np.asarray(ftsPrf)
                        ftsPrf = ftsPrf.flatten()

                        sonPrf = np.asarray(sonPrf)
                        sonPrf = sonPrf.flatten()

                        #------------------------------
                        # Root Mean Square Error (RMSE)
                        #------------------------------
                        ss_res = np.nansum( x**2)
                        rmseDF.append(np.sqrt( ss_res / len(x) ))

                        Npnts.append(float(np.sum(NobsFTS[df+v])))
                        Npnts2.append(float(len(NobsFTS[df+v])))

                    Npnts  = np.asarray(Npnts)
                    Npnts2 = np.asarray(Npnts2)
                    stdDF  = np.asarray(stdDF)
                    #with open('/data/iortega/Manuscripts/Fig/'+fileName, 'wb') as fopen:
                    #    fopen.write('{0:.4f}\t{1:.4f}\t{2:.4f}\n'.format(Npnts2, Npnts, np.asarray(stdDF)*100.))

                    if pn == 0:

                        ax1[0].plot(diffT, Npnts2,  color='k',linewidth=2.0, zorder=3)
                        ax1[0].scatter(diffT, Npnts2, facecolors='white', color='k',s=60, label='Number of dates', zorder=4)
                        
                        ax1[0].plot(diffT, Npnts,  color='blue',linewidth=2.0, zorder=1)
                        ax1[0].scatter(diffT, Npnts, facecolors='white', color='blue',s=60, label='Number of profiles', zorder=2)
                        

                        print 'Number of profiles = {}'.format(Npnts)
                        print 'Number of dates = {}'.format(Npnts2)
                        print 'diffT = {}'.format(diffT)

                    ax1[1].plot(diffT, np.asarray(stdDF)*100.,   color=clr[pn], linewidth=2.0, zorder=1)
                    ax1[1].scatter(diffT, np.asarray(stdDF)*100., color=clr[pn], s=60, facecolors='white', label= str(pcol[0])+'-'+str(pcol[1])+' km', zorder=2)
                    
                
                ax1[0].grid(True, alpha=0.35) 
                ax1[0].tick_params(which='both',labelsize=16)
                ax1[0].set_ylabel('Number of observations', fontsize=16)
                ax1[0].set_ylim(bottom=0)
                ax1[0].legend(prop={'size':11}, loc=2)

                ax1[1].grid(True, alpha=0.5)
                ax1[1].set_ylabel('Variability [%]', fontsize=16)
                ax1[1].tick_params(which='both',labelsize=16)
                ax1[1].set_ylim(bottom=0)
                ax1[1].legend(prop={'size':10.5}, loc=2)
                ax1[1].set_xlim(left=0, right=np.max(diffT) + 10)
                ax1[1].set_xlabel('$\Delta$t [min]', fontsize=18)

                if loc.lower() == 'fl0': ptitle = 'BLD'
                elif loc.lower() == 'mlo': ptitle = 'MLO'

                plt.suptitle(ptitle, fontsize=16)

                fig1.subplots_adjust(bottom=0.075,top=0.95, left=0.1, right=0.95)

                idver = getidver(v)


                if pltInputs['saveFlg']: 
                    pdfsav.savefig(fig1,dpi=200)
                    plt.savefig(pltDir+'Stats_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
               
                else:       
                    plt.show(block=False)

       
        #-------------------------------------------------
        # Profile position diference at specific Delta T
        #-------------------------------------------------
        if fleoutFlg:
           
            if pltInputs['timeFlg'].lower() == 'int':

                for v in pltInputs['ver']:
            
                    clmap     = 'jet'
                    fig, ax   =plt.subplots(figsize=(7, 7))
                    cm        = plt.get_cmap(clmap)
                  
                    norm = colors.Normalize(vmin=-15, vmax=30)

                    #times    =  [15, 30, 60, 120, 240]
                    times    =  [30, 60, 90, 120, 150]
                    markers  =  ['o', '>', 'D', '*', 's', '<']

                    for ti, t in enumerate(times):

                        if (ti == 0) or (ti ==3):

                            inds = np.where(np.asarray(diffT) == t)[0]

                            print inds.shape
                            print inds

                            df = str(diffT[inds[0]])

                            #PrfDiff_i = np.asarray(PrfDiffRel[diffT[inds]])
                    
                            #------------------------------
                            # Root Mean Square Error (RMSE)
                            #------------------------------
                            #SS_res = np.nansum( (Prfsonde_interp[diffT[inds]] - prfmean_day[diffT[inds]])**2, axis=0 )
                            #SS_res = np.asarray(SS_res)                
                            #rmse = np.sqrt( SS_res / len(Prfsonde_interp[diffT[inds]]) )

                            #------------------------------
                            # Bias
                            #------------------------------
                            biasCalc = np.nansum(prfmean_day[df+v] - Prfsonde_interp[df+v], axis=0 ) /Prfsonde_interp[df+v].shape[0] 
                            PrfDiff_i = biasCalc/np.nanmean(Prfsonde_interp[df+v], axis=0) * 100.


                            biasCalc_e = np.nansum(prfSTD_day[df+v], axis=0 ) /Prfsonde_interp[df+v].shape[0] 
                            PrfDiff_i_e = biasCalc_e/np.nanmean(Prfsonde_interp[df+v], axis=0) * 100.

                            Nx = prfmean_day[df+v].shape[0]
                            Ny = prfmean_day[df+v].shape[1]

                            Nx = int(np.sum(NobsFTS[df+v]))

                            PrfDist_i = np.asarray(PrfDist[df+v])

                            PrfDistMean  = np.nanmean(PrfDist_i, axis=0)
                            PrfDistStd   = np.nanstd(PrfDist_i, axis=0)
                            PrfDisMax    = np.nanmax(PrfDist_i, axis=0)
                            PrfDisMin    = np.nanmin(PrfDist_i, axis=0)
                            
                            #PrfDifftMean = np.nanmean(PrfDiff_i, axis=0)
                            #PrfDifftStd  = np.nanmean(PrfDiff_i, axis=0)

                            print '\nMean Relative difference for Delta T = {} minutes'.format(t)

                            for pn, pcol in enumerate(pCols):
                                
                                indsH = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <pcol[1])  )[0]
                                print 'Mean Horizontal DIfference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmean(PrfDistMean[indsH]), np.nanmean(PrfDistStd[indsH]), Nx)
                                print 'Median Horizontal DIfference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmedian(PrfDistMean[indsH]), np.nanmean(PrfDistStd[indsH]), Nx)
                                print 'Mean Bias Difference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}'.format(pcol,  np.nanmean(PrfDiff_i[indsH]), np.nanstd(PrfDiff_i[indsH]), Nx)
                                print 'Median Bias Difference in Layer {0:} = {1:.1f} +/- {2:.1f}, N= {3:}\n'.format(pcol,  np.nanmedian(PrfDiff_i[indsH]), np.nanstd(PrfDiff_i[indsH]), Nx)

                            #ax.plot(PrfDistMean, alt, color='gray', alpha=0.5, zorder=1)
                            if ti == 0: lab = '0-'+str(t) + ' min (N = {})'.format(Nx)
                            else: lab = str(times[ti-1])+'-'+str(t) + ' min (N = {})'.format(Nx)
                            #sc1 = ax.scatter(PrfDistMean, alt, c=PrfDiff_i, s=50, cmap=cm, norm=norm, zorder=2, marker=markers[ti], label=lab)
                            #sc1 = ax.scatter(PrfDistMean, alt, s=50,  marker=markers[ti], label=lab)

                            ax.plot(PrfDistMean, fts.alt[v],   color=clr[ti], linewidth=2, zorder=1)
                            #ax.errorbar(PrfDistMean, fts.alt[v],xerr=PrfDistStd, markersize=60, linestyle='-', color=clr[ti], ecolor=clr[ti])
                            ax.scatter(PrfDistMean, fts.alt[v], color=clr[ti], s=60, label= lab, facecolors='white', zorder=2)

                            #ax.fill_betweenx(fts.alt[v],PrfDistMean-PrfDistStd,PrfDistMean+PrfDistStd, alpha=0.25, color=clr[ti])


                            #ax.plot(PrfDisMax, fts.alt[v],   color=clr[ti], linewidth=2,linestyle='--', zorder=1)
                            #ax.plot(PrfDisMin, fts.alt[v],   color=clr[ti], linewidth=2,linestyle='--', zorder=1)
                            #ax.errorbar(PrfDistMean, fts.alt[v],xerr=PrfDistStd, markersize=60, linestyle='-', color=clr[ti], ecolor=clr[ti])

                            #print np.transpose([PrfDiff_i, alt])

                        #for pi in range(Nx):
                        #    #ax.plot(PrfDist_i[pi, :], alt, color='gray', alpha=0.5, zorder=1)
                        #    sc1 = ax.scatter(PrfDist_i[pi, :], alt, c=PrfDiff_i[pi, :], s=30, cmap=cm, norm=norm, zorder=2)
                        
                    ax.set_ylabel('Altitude [km]', fontsize=16)
                    ax.set_xlabel('Spatial mismatch [km]', fontsize=16)
                    ax.grid(True, alpha=0.35)
                    ax.legend(prop={'size':12}, loc=2)
                    ax.set_ylim(0, 20)
                    ax.set_xlim(xmin=0)

                    #cax  = fig.add_axes([0.86, 0.1, 0.03, 0.8])
                    
                    #cbar = fig.colorbar(sc1, cax=cax,format='%3i')
                    #cbar.set_label('Bias [%]', fontsize=13)
                    #ax.set_title('Averaging Kernels Scale Factor')
                    
                    ax.tick_params(which='both',labelsize=16)       

                    #fig.subplots_adjust(right=0.82) 

                    idver = getidver(v)

                    if pltInputs['saveFlg']: 
                        pdfsav.savefig(fig,dpi=200)
                    
                        plt.savefig(pltDir+'Dist_Prf_'+loc.upper()+'_'+idver+'.pdf', bbox_inches='tight')
                        
                    else:           
                        plt.show(block=False)
                        #user_input = raw_input('Press any key to exit >>> ')
                        #sys.exit()           # Exit program    

        #---------------------------------
        # CALCULATION OF BIAS (MEDIAN OF DIFFERENCES), PRECISION (STDV OF RESIDUALS), CORRELATION AND HISTOGRAM OF BIAS
        # PARTIAL COLUMNS
        #---------------------------------
        bias         = {}
        bias_e       = {}

        prec         = {}
        prec_e       = {}

        bias_perc    = {}
        bias_perc_e  = {}

        prec_perc    = {}
        prec_perc_e  = {}

        biasApr      = {}
        bias_eApr    = {}

        precApr      = {}
        prec_eApr    = {}

        slope        = {}
        slope_e      = {}
        
        intercept    = {}
        intercept_e  = {}
        
        rvalue       = {}

        slopeApr     = {}
        slope_eApr   = {}
        
        interceptApr ={}
        intercept_eApr = {}
        rvalueApr    = {}
        
        for v in pltInputs['ver']: 

            for df in diffT:

                df = str(df)

                #---------------------------------
                # Figure: HISTOGRAMS FOR DIFFERENT ALTITUDE LAYERS
                # link: http://viceroy.eeb.uconn.edu/EstimateS/EstimateSPages/EstSUsersGuide/References/WaltherAndMoore2005.pdf
                #---------------------------------
                
                for p, pcol in enumerate(pCols[0:-1]):

                    inds = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <pcol[1])  )[0]
                    #--------------
                    #Bias and precision in Retrieval
                    #--------------
                    bias_n = np.sum(prfmean_day[df+v][:, inds] - Prfsonde_interp[df+v][:, inds], axis=1 )/float(len(inds))#Prfsonde_interp[df][:, inds].shape[1] 
                    #biasCalc = prfmean_day[df][:, inds] - Prfsonde_interp[df][:, inds]
                    bias_n = bias_n[~np.isnan(bias_n)]/1e3

                    mu = np.nanmean(bias_n)
                    me = np.nanmedian(bias_n)
                    sigma = np.std(bias_n)
                    stdE = sigma/np.sqrt(len(bias_n))      #http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf

                    prec_n = sigma/np.sqrt(len(bias_n)) * 2.

                    bias.setdefault(df+v, []).append(me)
                    bias_e.setdefault(df+v, []).append(stdE)

                    bias_perc.setdefault(df+v, []).append(me/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    bias_perc_e.setdefault(df+v, []).append(stdE/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)

                    prec.setdefault(df+v, []).append(prec_n)
                    #prec_e.setdefault(df+v, []).append(stdE * 0.71)
                    prec_e.setdefault(df+v, []).append(stdE * 0.0)

                    prec_perc.setdefault(df+v, []).append(prec_n/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    #prec_perc_e.setdefault(df+v, []).append((stdE * 0.71)/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    prec_perc_e.setdefault(df+v, []).append((stdE * 0.0)/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)

                    #--------------
                    #Bias and precision in Apriori
                    #--------------
                    bias_nApr = np.sum(aPrf_day[df+v][:, inds] - Prfsonde_interp[df+v][:, inds], axis=1 )/float(len(inds))#Prfsonde_interp[df][:, inds].shape[1] 
                    bias_nApr = bias_nApr[~np.isnan(bias_nApr)]/1e3

                    muApr = np.nanmean(bias_nApr)
                    meApr = np.nanmedian(bias_nApr)
                    sigmaApr = np.std(bias_nApr)
                    stdEApr = sigmaApr/np.sqrt(len(bias_nApr)) * 2.  #http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf
                    #add a line showing the expected distribution



                    biasApr.setdefault(df+v, []).append(meApr)
                    bias_eApr.setdefault(df+v, []).append(sigmaApr/np.sqrt(len(bias_nApr)))

                    precApr.setdefault(df+v, []).append(sigmaApr/np.sqrt(len(bias_nApr)) * 2.)
                    #prec_eApr.setdefault(df+v, []).append(sigmaApr/np.sqrt(len(bias_nApr)) * 0.71)
                    prec_eApr.setdefault(df+v, []).append(sigmaApr/np.sqrt(len(bias_nApr)) * 0.0)

                    #--------------
                    #Orthogonal Regression in Retrieval
                    #--------------
                    xx = Prfsonde_interp[df+v][:, inds]
                    xx = xx[~np.isnan(xx)]/1e3

                    xx_e = Prfsonde_sd_interp[df+v][:, inds]
                    xx_e = xx_e[~np.isnan(xx_e)]/1e3

                    yy = prfmean_day[df+v][:, inds]
                    yy = yy[~np.isnan(yy)]/1e3

                    yy_e = prferr_day[df+v][:, inds]
                    yy_e = yy_e[~np.isnan(yy_e)]/1e3

                    odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx_e, yerr=yy_e, InError=True)
                    slope.setdefault(df+v, []).append(float(odr[0]))
                    intercept.setdefault(df+v, []).append(float(odr[1]))

                    slope_e.setdefault(df+v, []).append(float(odrErr[0]))
                    intercept_e.setdefault(df+v, []).append(float(odrErr[1]))

                    slopelr, interceptlr, r_valueln, p_valuelr, std_errlr = stats.linregress(xx, yy)
                    rvalue.setdefault(df+v, []).append(float(r_valueln))
                    
                    #--------------
                    #Orthogonal Regression in Apriori
                    #--------------
                    yyApr = aPrf_day[df+v][:, inds]
                    yyApr = yyApr[~np.isnan(yyApr)]/1e3

                    odrApr, odrErrApr  = mf.orthoregress(xx, yyApr, xerr=xx_e, yerr=yyApr*0.,InError=True)
                    slopeApr.setdefault(df+v, []).append(float(odrApr[0]))
                    interceptApr.setdefault(df+v, []).append(float(odrApr[1]))

                    slope_eApr.setdefault(df+v, []).append(float(odrErrApr[0]))
                    intercept_eApr.setdefault(df+v, []).append(float(odrErrApr[1]))

                    slope2Apr, intercept2Apr, r_value2Apr, p_value2Apr, std_err2Apr = stats.linregress(xx, yyApr)
                    
                    rvalueApr.setdefault(df+v, []).append(float(r_value2Apr))

                    #--------------
                    #Print in Terminal 
                    #--------------

                    if df == str(dfValue):

                        print '\nBias and Correlation - version: {} and time {}'.format(v, df)
                        print 'Number of Profiles: {}'.format(prfmean_day[df+v][:, inds].shape)
                    
                        print 'Altitude Layer: '+str(pcol[0])+' - '+str(pcol[1])+' km'
                        #print 'bias (vmr)        = {0:.2f}'.format(biasCalc2/1e3)
                        print 'Mean Bias (vmr)   = {0:.3f} +/- {1:.3f}'.format(mu, stdE)
                        print 'Median Bias (vmr)   = {0:.3f} +/- {1:.3f}'.format(me, stdE)
                        print 'Percent wrst Mean Sonde = {0:.3f} +/- {1:.3f}'.format(mu/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100., stdE/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                        print 'Percent wrst Median Sonde = {0:.3f} +/- {1:.3f}'.format(me/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100., stdE/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)

                        print '\nPrecision (vmr)   = {0:.3f} +/- {1:.3f}'.format(prec_n, stdE * 0.71)
                        print 'Precision Percent wrst Sonde = {0:.3f} +/- {1:.3f}'.format(prec_n/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100., (stdE * 0.71)/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                        
                        print '\nSlope: {0:.2f} +/- {1:.2f}'.format(float(odr[0]), float(odrErr[0]))
                        print 'Intercept = {0:.3f} +/- {1:.3f}'.format(float(odr[1]), float(odrErr[1]))
                        print 'R value = {0:.2f}'.format(float(r_valueln))
                        print 'Slope Apriori: {0:.2f} +/- {1:.2f}'.format(float(odrApr[0]), float(odrErrApr[0]))
                        print 'Intercept Apriori = {0:.3f} +/- {1:.3f}'.format(float(odrApr[1]), float(odrErrApr[1]))
                        print 'R value Apriori = {0:.2f}'.format(float(r_value2Apr)) 

        
        #---------------------------------
        # CALCULATION OF BIAS (MEDIAN OF DIFFERENCES), PRECISION (STDV OF RESIDUALS), CORRELATION AND HISTOGRAM OF BIAS
        # PARTIAL COLUMNS
        #---------------------------------
        bias2         = {}
        bias_e2       = {}

        prec2         = {}
        prec_e2       = {}

        bias_perc2    = {}
        bias_perc_e2  = {}

        prec_perc2    = {}
        prec_perc_e2  = {}

        slope2        = {}
        slope_e2      = {}
        
        intercept2    = {}
        intercept_e2  = {}
        
        rvalue2       = {}

        for v in pltInputs['ver']: 

            for p, pcol in enumerate(pCols[0:-1]):

                pcolstr = str(pcol[0])

                for df in diffT:

                    df = str(df)

                    inds = np.where( (fts.alt[v] >= pcol[0]) & (fts.alt[v] <pcol[1])  )[0]
                    #--------------
                    #Bias and precision in Retrieval
                    #--------------
                    bias_n = np.sum(prfmean_day[df+v][:, inds] - Prfsonde_interp[df+v][:, inds], axis=1 )/float(len(inds))
                    bias_n = bias_n[~np.isnan(bias_n)]/1e3

                    mu = np.nanmean(bias_n)
                    me = np.nanmedian(bias_n)
                    sigma = np.std(bias_n)
                    stdE = sigma/np.sqrt(len(bias_n))      #http://www.mas.ncl.ac.uk/~njnsm/medfac/docs/se&ci.pdf
                    prec_n = sigma/np.sqrt(len(bias_n)) * 2.

                    bias2.setdefault(pcolstr+v, []).append(me)
                    bias_e2.setdefault(pcolstr+v, []).append(stdE)

                    bias_perc2.setdefault(pcolstr+v, []).append(me/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    bias_perc_e2.setdefault(pcolstr+v, []).append(stdE/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)

                    prec2.setdefault(pcolstr+v, []).append(prec_n)
                    #prec_e2.setdefault(pcolstr+v, []).append(stdE * 0.71)
                    prec_e2.setdefault(pcolstr+v, []).append(stdE * 0.0)

                    prec_perc2.setdefault(pcolstr+v, []).append(prec_n/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    #prec_perc_e2.setdefault(pcolstr+v, []).append((stdE * 0.71)/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)
                    prec_perc_e2.setdefault(pcolstr+v, []).append((stdE * 0.0)/np.nanmean(Prfsonde_interp[df+v][:, inds]/1e3) * 100.)

                    #--------------
                    #Orthogonal Regression in Retrieval
                    #--------------
                    xx = Prfsonde_interp[df+v][:, inds]
                    xx = xx[~np.isnan(xx)]/1e3

                    xx_e = Prfsonde_sd_interp[df+v][:, inds]
                    xx_e = xx_e[~np.isnan(xx_e)]/1e3

                    yy = prfmean_day[df+v][:, inds]
                    yy = yy[~np.isnan(yy)]/1e3

                    yy_e = prferr_day[df+v][:, inds]
                    yy_e = yy_e[~np.isnan(yy_e)]/1e3

                    odr, odrErr  = mf.orthoregress(xx, yy, xerr=xx_e, yerr=yy_e, InError=True)
                    slope2.setdefault(pcolstr+v, []).append(float(odr[0]))
                    intercept2.setdefault(pcolstr+v, []).append(float(odr[1]))

                    slope_e2.setdefault(pcolstr+v, []).append(float(odrErr[0]))
                    intercept_e2.setdefault(pcolstr+v, []).append(float(odrErr[1]))

                    slope2lr, intercept2lr, r_value2lr, p_value2lr, std_err2lr = stats.linregress(xx, yy)
                    rvalue2.setdefault(pcolstr+v, []).append(float(r_value2lr))
                    
                   
        #-------------------
        #Bar Plot with slope, intercept and Rvalue
        #-------------------
        colors2 = 'rgbcmyk'

        for df in diffT:

            df = str(df)

            if df == dfValue:

                fig, (ax, ax2, ax3)  = plt.subplots(3, 1, figsize=(7.5,9), sharex=True)

                ind = np.arange(len(pCols[0:-1]))                        
                width = 1. / (len(pltInputs['ver']) + 1)

                bar_groups = []
                labels     = []
                for c, v in enumerate(pltInputs['ver']):
                    
                    labstr = v.strip().split('/')
                    if labstr[1] == 'Current_ERA': labels.append('ERA-d')
                    if labstr[1] == 'Current_ERA_v66': labels.append('ERA-6')
                    if labstr[1] == 'Current_WACCM': labels.append('WACCM')
                    if labstr[1] == 'Current_NCEP': labels.append('NCEP-d')
                    
                    bars = ax.bar(ind+c*width, slope[df+v], width, yerr=slope_e[df+v], ecolor='k', color=colors2[c % len(colors2)])
                    bar_groups.append(bars)

                ax.set_xlim(-width,len(ind)+width*0.4)
                ax.set_ylabel('Slope', fontsize=16)
                ax.set_ylim(0,1.5)
                ax.axhline(y=1.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
                ax.legend([b[0] for b in bar_groups], labels, fontsize=12, ncol=len(labels), loc='upper center', bbox_to_anchor=(0.5, 1.05))
                #ax.text(0.01, 0.95, '(a)', va='center', ha='left', transform=ax.transAxes,fontsize=14)

                ax.set_xticks(ind+width*2)
                ax.tick_params(which='both',labelsize=16)

                for c, v in enumerate(pltInputs['ver']):
                    intercept_v = np.asarray(intercept[df+v], dtype='float')
                    intercept_v[-3:] = intercept_v[-3:]*10.

                    intercept_ve = np.asarray(intercept_e[df+v], dtype='float')
                    intercept_ve[-3:] = intercept_ve[-3:]*10.

                    ax2.bar(ind+c*width, intercept_v, width, yerr=intercept_ve, ecolor='k', color=colors2[c % len(colors2)])

                ax2.set_xlim(-width,len(ind)+width*0.4)
                ax2.set_ylabel('Intercept [x10$^3$ ppm$_v$]', fontsize=16)
                ax2.set_xticks(ind+width*2)
                ax2.tick_params(which='both',labelsize=16)
                if loc.lower() =='fl0': 
                    ax2.set_ylim(-0.8,0.45)
                    ax2.text(0.56, 0.8, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)
                    ax2.text(0.71, 0.8, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)
                    ax2.text(0.87, 0.7, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)

                elif loc.lower() =='mlo': 
                    ax2.set_ylim(-0.15,0.15)
                    ax2.text(0.6, 0.94, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)
                    ax2.text(0.71, 0.62, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)
                    ax2.text(0.87, 0.62, 'x10', va='center', ha='left', transform=ax2.transAxes,fontsize=14)
                    


                ax2.axhline(y=0.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
                #ax2.text(0.01, 0.95, '(b)', va='center', ha='left', transform=ax2.transAxes,fontsize=14)

                for c, v in enumerate(pltInputs['ver']):
                    ax3.bar(ind+c*width, rvalue[df+v], width, color=colors2[c % len(colors2)])

                ax3.set_xlim(-width,len(ind)+width*0.4)
                ax3.set_ylabel('r-value', fontsize=16)
                xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
                ax3.set_xticks(ind+width*2)
                xtickNames = ax3.set_xticklabels(xTickMarks)
                plt.setp(xtickNames, rotation=0, fontsize=11)
                ax3.tick_params(which='both',labelsize=16)
                ax3.set_xlabel('Layer [km]', fontsize=16)
                ax3.set_ylim(0,1)
                #ax3.text(0.01, 0.95, '(c)', va='center', ha='left', transform=ax3.transAxes,fontsize=14)

                if loc.lower() == 'fl0': ptitle = 'BLD'
                elif loc.lower() == 'mlo': ptitle = 'MLO'

                plt.suptitle(ptitle, fontsize=16)

                fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)
                
                if pltInputs['saveFlg']: 
                    pdfsav.savefig(fig,dpi=200)
                    if pltInputs['smthFlg']: 
                        plt.savefig(pltDir+'StatsBars_'+loc.upper()+'_'+df+'_smooth.pdf', bbox_inches='tight')
                    else:
                        plt.savefig(pltDir+'StatsBars_'+loc.upper()+'_'+df+'.pdf', bbox_inches='tight')
                    #plt.savefig('/data/iortega/Manuscripts/Fig/StatsBars_'+self.loc.upper()+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)

        #-------------------
        #Bar Plot with Bias
        #-------------------
        for df in diffT:

            df = str(df)

            if df == dfValue:

                fig, (ax, ax2, ax3, ax4)  = plt.subplots(4, 1, figsize=(7.5,9), sharex=True)

                ind = np.arange(len(pCols[0:-1]))                        
                width = 1. / (len(pltInputs['ver']) + 1)

                bar_groups = []
                labels     = []
                for c, v in enumerate(pltInputs['ver']):
                    
                    labstr = v.strip().split('/')
                    if labstr[1] == 'Current_ERA': labels.append('ERA-d')
                    if labstr[1] == 'Current_ERA_v66': labels.append('ERA-6')
                    if labstr[1] == 'Current_WACCM': labels.append('WACCM')
                    if labstr[1] == 'Current_NCEP': labels.append('NCEP-d')

                    bias_v = np.asarray(bias[df+v], dtype='float')
                    bias_v[-2:] = bias_v[-2:]*10.

                    bias_ve = np.asarray(bias_e[df+v], dtype='float')
                    bias_ve[-2:] = bias_ve[-2:]*10.
                    
                    bars = ax.bar(ind+c*width, bias_v, width, yerr=bias_ve, ecolor='k', color=colors2[c % len(colors2)])
                    bar_groups.append(bars)

                
                ax.set_xlim(-width,len(ind)+width*0.4)
                ax.set_ylabel('Bias [x10$^3$ ppm$_v$]', fontsize=16)
                #ax.set_ylim(0,1.5)
                ax.axhline(y=0.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
                ax.legend([b[0] for b in bar_groups], labels, fontsize=12, ncol=len(labels), loc='upper center', bbox_to_anchor=(0.5, 1.05))
                #ax.text(0.01, 0.95, '(a)', va='center', ha='left', transform=ax.transAxes,fontsize=14)

                ax.set_xticks(ind+width*2)
                ax.tick_params(which='both',labelsize=16)
                if loc.lower() == 'mlo':   
                    ax.set_ylim(-0.3,0.2)
                    ax.text(0.69, 0.5, 'x10', va='center', ha='left', transform=ax.transAxes,fontsize=14)
                    ax.text(0.85, 0.5, 'x10', va='center', ha='left', transform=ax.transAxes,fontsize=14)

                elif loc.lower() == 'fl0': 
                    ax.set_ylim(-0.5, 0.3)
                    ax.text(0.69, 0.5, 'x10', va='center', ha='left', transform=ax.transAxes,fontsize=14)
                    ax.text(0.85, 0.5, 'x10', va='center', ha='left', transform=ax.transAxes,fontsize=14)
                #ax.annotate('x10', xy=(0.65, 0.12), size=16, xycoords='axes fraction', xytext=(0.69, 0.45), textcoords='axes fraction', arrowprops=dict(facecolor='black', shrink=0.05))

                for c, v in enumerate(pltInputs['ver']):
                    ax2.bar(ind+c*width, bias_perc[df+v], width, yerr=bias_perc_e[df+v], ecolor='k', color=colors2[c % len(colors2)])

                ax2.set_xlim(-width,len(ind)+width*0.4)
                ax2.set_ylabel('Bias [%]', fontsize=16)
                ax2.set_xticks(ind+width*2)
                ax2.tick_params(which='both',labelsize=16)
                ax2.set_ylim(-15, 30)

                ax2.axhline(y=0.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
                #ax2.text(0.01, 0.95, '(b)', va='center', ha='left', transform=ax2.transAxes,fontsize=14)

                for c, v in enumerate(pltInputs['ver']):

                    prec_v = np.asarray(prec[df+v], dtype='float')
                    prec_v[-2:] = prec_v[-2:]*10.

                    prec_ve = np.asarray(prec_e[df+v], dtype='float')
                    prec_ve[-2:] = prec_ve[-2:]*10.
                    
                    ax3.bar(ind+c*width, prec_v, width, yerr=prec_ve, ecolor='k', color=colors2[c % len(colors2)])

                ax3.set_xlim(-width,len(ind)+width*0.4)
                ax3.set_ylabel('Precision [x10$^3$ ppm$_v$]', fontsize=16)
                ax3.set_xticks(ind+width*2)
                ax3.tick_params(which='both',labelsize=16)
                ax3.axhline(y=0.0, linestyle='--', linewidth=1.5, color='k', alpha=0.5)
                #ax2.text(0.01, 0.95, '(b)', va='center', ha='left', transform=ax2.transAxes,fontsize=14)

                if loc.lower() == 'mlo': 
                    ax3.text(0.69, 0.15, 'x10', va='center', ha='left', transform=ax3.transAxes,fontsize=14)
                    ax3.text(0.85, 0.075, 'x10', va='center', ha='left', transform=ax3.transAxes,fontsize=14)
                elif loc.lower() == 'fl0':
                    ax3.text(0.69, 0.3, 'x10', va='center', ha='left', transform=ax3.transAxes,fontsize=14)
                    ax3.text(0.85, 0.15, 'x10', va='center', ha='left', transform=ax3.transAxes,fontsize=14)



                for c, v in enumerate(pltInputs['ver']):
                    ax4.bar(ind+c*width, prec_perc[df+v], width, yerr=prec_perc_e[df+v], ecolor='k', color=colors2[c % len(colors2)])

                ax4.set_xlim(-width,len(ind)+width*0.4)
                ax4.set_ylabel('Precision [%]', fontsize=16)
                xTickMarks = [str(pcol[0])+'-'+str(pcol[1]) for pcol in pCols]
                ax4.set_xticks(ind+width*2)
                xtickNames = ax4.set_xticklabels(xTickMarks)
                plt.setp(xtickNames, rotation=0, fontsize=12)
                ax4.tick_params(which='both',labelsize=16)
                ax4.set_xlabel('Layer [km]', fontsize=16)
                #ax3.text(0.01, 0.95, '(c)', va='center', ha='left', transform=ax3.transAxes,fontsize=14)

                if loc.lower() == 'fl0': ptitle = 'BLD'
                elif loc.lower() == 'mlo': ptitle = 'MLO'

                plt.suptitle(ptitle, fontsize=16)

                fig.subplots_adjust(bottom=0.075,top=0.95, left=0.15, right=0.95)
                
                if pltInputs['saveFlg']: 
                    pdfsav.savefig(fig,dpi=200)
                    if pltInputs['smthFlg']: 
                        plt.savefig(pltDir+'StatsBars2_'+loc.upper()+'_'+df+'_smooth.pdf', bbox_inches='tight')
                    else:
                        plt.savefig(pltDir+'StatsBars2_'+loc.upper()+'_'+df+'.pdf', bbox_inches='tight')
                else:
                    plt.show(block=False)

        #-------------------
        #
        #-------------------
        voi = pltInputs['ver'][0]#'Current_ERA'
        voistr = voi.strip().split('/')[1]
        #if voi[1] == 'Current_ERA': labels.append('ERA-d')

        fig1, ax1  = plt.subplots(4, figsize=(8,10), sharex=True)
        
        for c, pcol in enumerate(pCols[0:-1]):

            pcolstr = str(pcol[0])

            ax1[0].plot(diffT, slope2[pcolstr+voi],   color=colors2[c % len(colors2)], linewidth=2.0, zorder=1)
            ax1[0].scatter(diffT, slope2[pcolstr+voi], color=colors2[c % len(colors2)], s=60, facecolors='white', label= str(pcol[0])+'-'+str(pcol[1])+' km', zorder=2)

            ax1[1].plot(diffT, intercept2[pcolstr+voi],   color=colors2[c % len(colors2)], linewidth=2.0, zorder=1)
            ax1[1].scatter(diffT, intercept2[pcolstr+voi], color=colors2[c % len(colors2)], s=60, facecolors='white', zorder=2)

            ax1[2].plot(diffT, bias_perc2[pcolstr+voi], linewidth=2.0, zorder=1, color=colors2[c % len(colors2)])
            ax1[2].scatter(diffT, bias_perc2[pcolstr+voi], facecolors='white', color='k',s=60, zorder=2)

            ax1[3].plot(diffT, prec_perc2[pcolstr+voi],   color=colors2[c % len(colors2)], linewidth=2.0, zorder=1)
            ax1[3].scatter(diffT, prec_perc2[pcolstr+voi], color=colors2[c % len(colors2)], s=60, facecolors='white', zorder=2)

        ax1[0].grid(True) 
        ax1[0].tick_params(which='both',labelsize=16)
        ax1[0].set_ylabel('Slope', fontsize=16)
        #ax1[0].set_ylim(bottom=0)
        ax1[0].legend(prop={'size':12}, loc=2)

        ax1[1].grid(True) 
        ax1[1].tick_params(which='both',labelsize=16)
        ax1[1].set_ylabel('Intercept [x10$^3$ ppm$_v$]', fontsize=16)
        #ax1[0].set_ylim(bottom=0)
        #ax1[1].legend(prop={'size':11}, loc=2)

        ax1[2].grid(True) 
        ax1[2].tick_params(which='both',labelsize=16)
        ax1[2].set_ylabel('Bias [x10$^3$ ppm$_v$]', fontsize=16)
        #ax1[0].set_ylim(bottom=0)
        #ax1[2].legend(prop={'size':11}, loc=2)

        ax1[3].grid(True)
        ax1[3].set_ylabel('Precision [x10$^3$ ppm$_v$]', fontsize=16)
        ax1[3].tick_params(which='both',labelsize=16)
        #ax1[3].set_ylim(bottom=0)
        ax1[3].set_xlim(left=0, right=np.max(diffT) + 10)
        ax1[3].set_xlabel('$\Delta$t [min]', fontsize=16)

        if loc.lower() == 'fl0': ptitle = 'BLD'
        elif loc.lower() == 'mlo': ptitle = 'MLO'

        plt.suptitle(ptitle, fontsize=16)

        fig1.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)

        idver = getidver(v)

        if pltInputs['saveFlg']: 
            pdfsav.savefig(fig1,dpi=200)
            if pltInputs['smthFlg']: 
                plt.savefig(pltDir+'StatsTime_'+voistr+'_'+loc.upper()+'_smooth.pdf', bbox_inches='tight')
            else:
                plt.savefig(pltDir+'StatsTime_'+voistr+'_'+loc.upper()+'.pdf', bbox_inches='tight')

        else:       
            plt.show(block=False)

        #-------------------
        #
        #-------------------
        for pi, pcol in enumerate(pCols[0:-1]):

            if pi == 0:

                pcolstr = str(pcol[0])

                fig1, ax1  = plt.subplots(4, figsize=(8,10), sharex=True)

                labels     = []
                for c, v in enumerate(pltInputs['ver']):
                    
                    labstr = v.strip().split('/')
                    if labstr[1] == 'Current_ERA': labels = 'ERA-d'
                    if labstr[1] == 'Current_ERA_v66': labels = 'ERA-6'
                    if labstr[1] == 'Current_WACCM': labels = 'WACCM'
                    if labstr[1] == 'Current_NCEP': labels = 'NCEP-d'

                    ax1[0].plot(diffT, slope2[pcolstr+v],   color=colors2[c % len(colors2)], label=labels, linewidth=2.0, zorder=1)
                    ax1[0].scatter(diffT, slope2[pcolstr+v], color=colors2[c % len(colors2)], s=60, facecolors='white', zorder=2)

                    ax1[1].plot(diffT, intercept2[pcolstr+v],   color=colors2[c % len(colors2)], linewidth=2.0, zorder=1)
                    ax1[1].scatter(diffT, intercept2[pcolstr+v], color=colors2[c % len(colors2)], s=60, facecolors='white',  zorder=2)

                    ax1[2].plot(diffT, bias2[pcolstr+v], linewidth=2.0, zorder=1, color=colors2[c % len(colors2)], )
                    ax1[2].scatter(diffT, bias2[pcolstr+v], facecolors='white', color='k',s=60, zorder=2)

                    ax1[3].plot(diffT, prec2[pcolstr+v],   color=colors2[c % len(colors2)], linewidth=2.0, zorder=1)
                    ax1[3].scatter(diffT, prec2[pcolstr+v], color=colors2[c % len(colors2)], s=60, facecolors='white', zorder=2)

                ax1[0].grid(True) 
                ax1[0].tick_params(which='both',labelsize=16)
                ax1[0].set_ylabel('Slope', fontsize=16)
                #ax1[0].set_ylim(bottom=0)
                ax1[0].legend(prop={'size':12}, ncol=len(pltInputs['ver']), loc='upper center', bbox_to_anchor=(0.5, 1.05))

                ax1[1].grid(True) 
                ax1[1].tick_params(which='both',labelsize=16)
                ax1[1].set_ylabel('Intercept [x10$^3$ ppm$_v$]', fontsize=16)
                #ax1[0].set_ylim(bottom=0)
                #ax1[1].legend(prop={'size':11}, loc=2)

                ax1[2].grid(True) 
                ax1[2].tick_params(which='both',labelsize=16)
                ax1[2].set_ylabel('Bias [x10$^3$ ppm$_v$]', fontsize=16)
                #ax1[0].set_ylim(bottom=0)
                #ax1[2].legend(prop={'size':11}, loc=2)

                ax1[3].grid(True)
                ax1[3].set_ylabel('Precision [x10$^3$ ppm$_v$]', fontsize=16)
                ax1[3].tick_params(which='both',labelsize=16)
                ax1[3].set_ylim(bottom=0)
                ax1[3].set_xlim(left=0, right=np.max(diffT) + 10)
                ax1[3].set_xlabel('$\Delta$t [min]', fontsize=16)

                if loc.lower() == 'fl0': ptitle = 'BLD'
                elif loc.lower() == 'mlo': ptitle = 'MLO'

                plt.suptitle(ptitle, fontsize=16)

                fig1.subplots_adjust(bottom=0.075,top=0.95, left=0.12, right=0.95)

                idver = getidver(v)

                if pltInputs['saveFlg']: 
                    pdfsav.savefig(fig1,dpi=200)
                    if pltInputs['smthFlg']: 
                        plt.savefig(pltDir+'StatspCols_'+pcolstr+'_'+loc.upper()+'_smooth.pdf', bbox_inches='tight')
                    else:
                        plt.savefig(pltDir+'StatspCols_'+pcolstr+'_'+loc.upper()+'.pdf', bbox_inches='tight')
        
                else:       
                    plt.show(block=False)
     
    if pltInputs['saveFlg']: pdfsav.close()

    print('\nFinished Plots.......\n')
       
 #    #--------------------------------
 #    # Pause so user can look at plots
 #    #--------------------------------
    if not pltInputs['saveFlg']:
        user_input = raw_input('Press any key to exit >>> ')
        sys.exit()           # Exit program       

if __name__ == "__main__":
    main(sys.argv[1:])