#! /usr/bin/python

#----------------------------------------------------------------------------------------
# Name:
#        pullsondes.py
#
# Purpose:
#       This program pulls water vapor sondes from NOAA
#           
#
# Notes:
#       1) Command line arguments
#            -d   YYYY       : Specify the start and stop date to get data.
#                                         Default is previous utc day
#
# Usage:
#     pullsondes.py 
#
# Examples:
#    ./pullsondes.py 
#
# Version History:
#  1.0     Created, November, 2013  Ivan Ortega
#
#----------------------------------------------------------------------------------------


                            #-------------------------#
                            # Import Standard modules #
                            #-------------------------#

import datetime as dt
import sys
import os
import itertools
import getopt
import shutil
import subprocess as sp
import logging
import glob

                            #-------------------------#
                            # Define helper functions #
                            #-------------------------#
                            
def usage():
    ''' Prints to screen standard program usage'''
    print 'pullsondes.py [-d YYYY -s fl0]'                            

def subProcRun( sysCall, logF=False, shellFlg=False ):
    '''This runs a system command and directs the stdout and stderr'''
    rtn = sp.Popen( sysCall, stdout=sp.PIPE, stderr=sp.PIPE, shell=shellFlg )
    stdoutInfo, stderrInfo = rtn.communicate()

    if logF:
        if stdoutInfo: logF.info( stdoutInfo )
        if stderrInfo: logF.error( stderrInfo )
               
    return (stdoutInfo,stderrInfo)

def chMod(PrntPath):
    for dirpath, dirnames, filenames in os.walk(PrntPath):
        try:    os.chmod(dirpath,0o777)
        except: pass        
        for filename in filenames:
            path = os.path.join(dirpath, filename)
            try:    os.chmod(path, 0o777)   
            except: pass    

def ckDirMk(dirName,logFlg=False):
    ''' '''
    if not ( os.path.exists(dirName) ):
        os.makedirs( dirName, mode=0777 )
        if logFlg: logFlg.info( 'Created folder {}'.format(dirName))
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

def ckDir(dirName,logFlg=False,exitFlg=False):
    ''' '''
    if not os.path.exists( dirName ):
        print 'Input Directory %s does not exist' % (dirName)
        if logFlg: logFlg.error('Directory %s does not exist' % dirName)
        if exitFlg: sys.exit()
        return False
    else:
        return True   


def pull(web='', idstr='', ldir=''):
   
   #try:
        cmnd = ['wget', '-nc','-P'+ldir, web +'*'+idstr+'*.hdf']
        (stdoutInfo,stderrInfo) = subProcRun(cmnd)
    #except: print 'A'



                            #----------------------------#
                            #                            #
                            #        --- Main---         #
                            #                            #
                            #----------------------------#
   

def main():

    
    station    = 'kiruna'
    idstn      = 'kir'
    idname     = 'ocs'
    website    = 'ftp://ftp.cpc.ncep.noaa.gov/ndacc/station/'+station.lower()+'/hdf/ftir/'
    locDir     = '/data1/projects/ocs' 

    ckDir(locDir, exitFlg=True)
   
    if not( locDir.endswith('/') ):
        locDir = locDir + '/'

    finaldir   = locDir+idstn.lower()

    ckDirMk(finaldir)


    pull(web=website, idstr=idname, ldir=finaldir)


    
if __name__ == "__main__":
    main()
