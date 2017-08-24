#! /usr/bin/python
#----------------------------------------------------------------------------------------
# Name:
#        FTIRwebParser.py
#
# Purpose:
#       Interface with OPUS software to run scans and read results
#
#
#
# Notes:
#       
#
#
# License:
#    Copyright (c) 2013-2014 NDACC/IRWG
#    This file is part of sfit4.
#
#    sfit4 is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    any later version.
#
#    sfit4 is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with sfit4.  If not, see <http://www.gnu.org/licenses/>
#
#----------------------------------------------------------------------------------------


#---------------
# Import modules
#---------------		
import sys 		
import urllib2
import os
import socket
import select
import datetime       as     dt
from time             import sleep
from trackerUtils     import *
from remoteData       import FTIRdataClient
from bs4              import BeautifulSoup  as bs


class FTIRwebParser(FTIRdataClient):
    
    def __init__(self,EWS_IP="192.168.1.102",TCP_IP_in="192.168.1.100",TCP_Port_in=5555,BufferSize_in=4024):
        
        self.TCPserverIP = TCP_IP_in
        self.webAddress  = EWS_IP
        self.data        = {}
        
        FTIRdataClient.__init__(self,TCP_IP=self.TCPserverIP,TCP_Port=TCP_Port_in,BufferSize=BufferSize_in)
        
    
    def parseDiagnostics(self):
        
        #----------------------------
        # Check if address is visible
        #----------------------------      
        if not self.checkWebAddress(): return False                
        
        diagPage = "http://" + self.webAddress + "/brow_diag.htm"  
        webpg    = urllib2.urlopen(diagPage)
        data     = []
        
        #-----------------------------
        # Open page with Beautifulsoup
        #-----------------------------
        soup = bs(webpg)
        
        #--------------
        # Extract Table
        #--------------
        table = soup.find('table')
        
        #---------
        # Get rows
        #---------
        rows = table.find_all('tr')
        
        #------
        # Parse
        #------
        for row in rows:
            cols = row.find_all('td')
            cols = [element.text.strip() for element in cols]
            data.append([element for element in cols if element])
            
        #---------------------------
        # Add elements to dictionary
        #---------------------------
        for row in data:
            self.data["BRUKER_"+('_').join(row[0].upper().split())] = row[1].upper()   
            
            
    def parseCommandList(self):
        
        #----------------------------
        # Check if address is visible
        #----------------------------      
        if not self.checkWebAddress(): return False          
        
        cmndPage = "http://" + self.webAddress + "/config/cmdlist.htm"
        webpg    = urllib2.urlopen(cmndPage)
        data     = []
        getList  = ['APT','OF1','OPF','DTC','VAC','VLV','LSR','PGN','GNS','LPF','SG2','BMS','HFW','LFW','LWN','NSS','RDY',]
        
        #-----------------------------
        # Open page with Beautifulsoup
        #-----------------------------
        soup     = bs(webpg)
        
        #--------------
        # Extract Table
        #--------------
        table = soup.find('table')
        
        #---------
        # Get rows
        #---------
        rows = table.find_all('tr')
        
        #------
        # Parse
        #------
        for row in rows[1:]:
            cols = row.find_all('td')
            cols = [element.text.strip() for element in cols]
            data.append([element for element in cols if element])
            
        #---------------------------
        # Add elements to dictionary
        #---------------------------
        # Grab only certain elements
        #---------------------------
        for row in data:
            if row[0] in getList: self.data["BRUKER_"+row[0].upper()] = row[2].upper()        
            
            
    def parseMeasStatus(self):
        
        #----------------------------
        # Check if address is visible
        #----------------------------      
        if not self.checkWebAddress(): return False          

        measPage = "http://" + self.webAddress + "/brow_stat.htm"
        webpg    = urllib2.urlopen(measPage)
        data     = []
        
        #-----------------------------
        # Open page with Beautifulsoup
        #-----------------------------
        soup     = bs(webpg)
        
        #--------------
        # Extract Table
        #--------------
        table = soup.find('table')
        
        #---------
        # Get rows
        #---------
        rows = table.find_all('tr')
        
        #------
        # Parse
        #------
        for row in rows:
            cols = row.find_all('td')
            cols = [element.text.strip() for element in cols]
            data.append([element for element in cols if element])
            
        #---------------------------
        # Add elements to dictionary
        #---------------------------
        for row in data:
            self.data["BRUKER_"+('_').join(row[0].upper().split())] = row[1].upper()        
            
    def checkWebAddress(self):
        
        #------------------------
        # Check if address exists
        #------------------------           
        try:
            urllib2.urlopen("http://" + self.webAddress)
            return True
        except:
            print "Unable to open web address: {}".format(self.webAddress)
            return False
            
            
    def sendData(self):
        #--------------------------
        # Write data to data server
        #--------------------------
        if not self.data:
            print "No data to write from web parser..."
            return False
        
        for key in self.data:
            message = key + " " + self.data[key]
            self.setParam(message)
            
        #---------------------------------------------
        # Send time stamp of connection to data server
        #---------------------------------------------
        crntDateTime = dt.datetime.utcnow()
        dateTstr = "{0:}{1:02}{2:02}_{3:02}{4:02}{5:02}".format(crntDateTime.year,crntDateTime.month,crntDateTime.day,crntDateTime.hour,
                                                                crntDateTime.minute,crntDateTime.second)
        self.setParam("BRUKER_TS " + dateTstr)           


    def getDIOdata(self,crntDataDir):
        ''' This method parses the data from the DIO board output '''
        
        #--------------------------------
        # Open log file and get last line
        #--------------------------------
        with open(crntDataDir,'r') as fopen:
            fopen.seek(-1024,os.SEEK_END)
            lastline = fopen.readlines()[-1].decode()
        
        #--------------------------------------
        # Determine if we got valid data or not
        #--------------------------------------
        if any([not lastline,lastline[0]=="#"]): return False

        parts = lastline.strip().split()
        
        #-----------
        # Parse data
        #-----------
        # Date/Time
        #----------
        yr   = int(parts[0][:4])
        mnth = int(parts[0][4:6])
        day  = int(parts[0][6:])
        hr   = int(parts[1][:2])
        mn   = int(parts[1][3:5])
        sec  = int(parts[1][6:])
        tstamp = dt.datetime(yr,mnt,day,hr,mn,sec)
        self.data['HARDWARE_TSTAMP'] = tstamp
        
        self.data['LN2_DEWAR_PRESSURE']            = parts[2]
        self.data['OPTIC_BENCH_BASEPLATE_T']       = parts[3]
        self.data['BEAMSPLITTER_BODY_T']           = parts[4]
        self.data['INSB_BODY_T']                   = parts[5]
        self.data['MCT_BODY_T']                    = parts[6]
        self.data['BRUKER_OPTIC_RH']               = parts[7]
        self.data['E_RADIANCE']                    = parts[8]
        self.data['W_RADIANCE']                    = parts[9]
        self.data['E_RADIANCES']                   = parts[10]
        self.data['W_RADIANCES']                   = parts[11]
        self.data['OUTSIDE_RH']                    = parts[12]
        self.data['LASER_PS_T']                    = parts[13]
        self.data['WIND_SPEED']                    = parts[14]
        self.data['WIND_DIR_E_OF_N']               = parts[15]
        self.data['LASER_T']                       = parts[16]
        self.data['OUTSIDE_T']                     = parts[17]
        self.data['DEC_B_ELEC_CONTROLS']           = parts[18][2]
        self.data['MID-IR_COOLER']                 = parts[19][1]
        self.data['LN2_FILL']                      = parts[19][2]
        self.data['VACUUM_VALVE']                  = parts[19][3]
        self.data['VACUUM_PUMP']                   = parts[20][1]
        self.data['HATCH_RELAY']                   = parts[20][2]
        self.data['VACUUM_PUMP_CONTROL']           = parts[20][3]
        self.data['SOLAR_SEEKER_ON_RELAY']         = parts[21][0]
        self.data['SOLAR_SEEKER_OFF_RELAY']        = parts[21][1]
        self.data['DYN_MIRROR_PWR']                = parts[21][2]
        self.data['DEC_A_BRUKER']                  = parts[21][3]
        self.data['HATCH_POSITION']                = parts[22][0]
        self.data['28V_SOLAR_SEEKER_PWR']          = parts[22][2]
        self.data['DOY']                           = parts[25]
        self.data['DAYMIN']                        = parts[26]
        self.data['DAYFRAC']                       = parts[27]
        self.data['SZENITH']                       = parts[28]
        self.data['SAZIMUTH']                      = parts[29]


if __name__ == "__main__":
    
    #--------------------------------
    # Retrieve command line arguments
    #--------------------------------
    try:
        ctlFile  = sys.argv[1]
    except:
        ctlFile  = "/home/tabftir/remote/ops/TAB_Defaults.input" 

    #----------------------------
    # Check existance of ctl file
    #----------------------------
    ckFile(ctlFile,exitFlg=True)
    
    #-------------------------
    # Import ctlFile variables
    #-------------------------
    ctlFvars = mainInputParse(ctlFile)
    
    #-----------------------------
    # Initialize EWS parser object
    #-----------------------------
    ewsData = FTIRwebParser(EWS_IP=ctlFvars["Bruker_IP"],TCP_IP_in=ctlFvars["FTS_DataServ_IP"],
                            TCP_Port_in=int(ctlFvars["FTS_DataServ_PORT"]),BufferSize_in=int(ctlFvars["FTS_DATASERV_BSIZE"]))

    #-----------------------------------------
    # Start indefinite loop to update TCP data
    # server with Bruker EWS information
    #-----------------------------------------
    while True:
        
        #------------------
        # Parse diagnostics
        #------------------
        try:
            ewsData.parseDiagnostics()
        except:
            print "Unable to parse Bruker EWS Diagnostics page!"
    
        #-------------------
        # Parse command list
        #-------------------
        try:
            ewsData.parseCommandList()        
        except:
            print "Unable to parse Bruker EWS Command page!"
    
        #-------------------------
        # Parse measurement status
        #-------------------------
        try:
            ewsData.parseMeasStatus()     
        except:
            print "Unable to parse Bruker EWS Measurement Status page!"
        
        #------------------------
        # Send data to TCP server
        #------------------------
        try:
            ewsData.sendData()
        except:
            print "Unable to send data from web parser program to TCP data server!"
        
        #--------------------------
        # Hold for updated interval
        #--------------------------
        sleep(int(ctlFvars["Bruker_UpdateInt"]))
        

    Contact GitHub API Training Shop Blog About 

    © 2016 GitHub, Inc. Terms Privacy Security Status Help 

