#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2022-09-11'


# -- Built-in modules -- #
import os
# -- Third-part modules -- #
import datetime as dt
import numpy as np
from netCDF4 import Dataset
# -- Proprietary modules -- #
from read_ccip_AutoROI_final import collocation

def check_filesize(OBSID,SATELLITE, check_count, save_path):
        
        if check_count == 0:
            ofile = OBSID + '_' + SATELLITE + '_CCIp.out'
        else:
            ofile = OBSID + '_' + SATELLITE + str(check_count) + '_CCIp.out'
            
        # output file
        CompleteName = os.path.join(save_path, ofile)
        
        # check for file size
        try:
            file_size = os.path.getsize(CompleteName)
        except:
            file_size = 0
            print('New file just created')
        if file_size > 10e8: # if larger than 1GB start new file
            check_count += 1
            print(file_size)
        
        return ofile
    
#---------------------------------
''' Define these parameters'''
OBSID = "AEM-AWI" # Name of campaign
SATELLITE = "ERS2" # Satelitte name (either ENV, CS2, ERS1 or ERS2)
# name of reference data file
obsfile =  os.path.dirname(os.getcwd()) + "/FINAL/" + OBSID + "/final/ESACCIplus-SEAICE-RRDP2+-SIT-" + OBSID + ".dat" 
HS   = 'NH' # Hemisphere of observations (either Northen Hemisphere (NH) or Southern Hemisphere (SH))
# -------------------------------------------

## Read reference observations
ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')
# find relevant years and months in obsfile
ObsDate = ObsData['date']
formating = "%Y-%m-%dT%H:%M:%S"
ObsDate_object = [dt.datetime.strptime(Date, formating) for Date in ObsDate]
ObsYears = [ObsDate.year for ObsDate in ObsDate_object]
ObsMonth = [ObsDate.month for ObsDate in ObsDate_object]

# Location of output files
save_path =  os.path.dirname(os.getcwd()) + '/satellite/Satellite_subsets/' + OBSID

# if satellite is either ERS1 or ERS2
if 'ERS' in SATELLITE:
    # Satellite data directory
    sat_dir = os.path.dirname(os.getcwd()) + "/satellite/ERS_data/"+HS + "/"
    files = os.listdir(sat_dir)
    file = [f for f in files if SATELLITE in f][0]
    ifile = os.path.join(sat_dir, file)
    CCIpFile = Dataset(ifile,'r')
    # time unit of ERS files: days since 1950-01-01
    time = [dt.date(1950, 1, 1) + dt.timedelta(days = int(t)) for t in CCIpFile.variables['time'][:]] ;
    
    count = 0;
    check_count = 0; ## to check for file size
    ## Loop through all relevant dates in satellite file
    for t,i in zip(time, range(len(time))):
        if t.year in ObsYears:
            if t.month in ObsMonth:
                ofile = check_filesize(OBSID,SATELLITE, check_count, save_path)
                # Extract satellite subset
                collocation(OBSID, obsfile, ifile, ofile, count, i)
                count += 1


else:  # is satellite is ENV or CS2
    # Satellite data directory
    sat_dir = os.path.dirname(os.getcwd()) + "/satelitte/ENV_CS2_data/"+SATELLITE.replace('-','')+"_data/" + HS + "/"
    
    ## Loop through all satelitte files, with relevant dates
    years = os.listdir(sat_dir)
    count = 0;
    check_count = 0; ## to check for file size
    for year in years:
        ## Read only file if date is relevant
        if int(year) in np.unique(ObsYears):
            print(year)
            months = os.listdir(sat_dir + year)
            for month in months:
                ## Find months in obs in year
                obsMonth_in_year = [ObsDate.month for ObsDate in ObsDate_object if ObsDate.year==int(year)]
                if int(month) in np.unique(obsMonth_in_year):
                    files = os.listdir(sat_dir + year + '/' + month)
                    for file in files:
                        print(file)
                        day = file[-14:-12]
                        
                        ofile = check_filesize(OBSID,SATELLITE, check_count, save_path)
                        
                        ifile = os.path.join(sat_dir + year + '/' + month, file)
                        # Extract satellite subset
                        collocation(OBSID, obsfile, ifile, ofile, count)
                        count += 1


                