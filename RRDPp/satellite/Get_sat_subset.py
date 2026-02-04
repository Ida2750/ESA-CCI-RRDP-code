#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2022-09-11'


# -- Built-in modules -- #
import os
import sys
# -- Third-part modules -- #
import datetime as dt
import numpy as np
from netCDF4 import Dataset
import xarray as xr
# -- Proprietary modules -- #
from read_ccip_AutoROI_final import collocation
sys.path.append(os.path.dirname(os.getcwd()) + '/code') 
import Functions

def check_filesize(OBSID,SATELLITE, check_count, save_path, HS=''):
        
        if check_count == 0:
            ofile = OBSID + '_' + SATELLITE + f'-{HS}_CCIp.out'
        else:
            ofile = OBSID + '_' + SATELLITE + str(check_count) + f'-{HS}_CCIp.out'
            
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
        
        return ofile, check_count

#---------------------------------
''' Define these parameters'''
if len(sys.argv) != 4:
    print("Usage: python script.py <HS> <SATELLITE> <OBSID>")
    sys.exit(1)

HS = sys.argv[1]
SATELLITE = sys.argv[2]
OBSID = sys.argv[3]

print(f"Running script with HS={HS}, SATELLITE={SATELLITE}, OBSID={OBSID}")

if HS=='SH':
    inpath = os.path.dirname(os.getcwd()) + f"/FINAL/Antarctic/{OBSID}/final"
else:
    inpath = os.path.dirname(os.getcwd()) + f"/FINAL/{OBSID}/final"
# get files
files = os.listdir(inpath)
print(inpath)
file = [f for f in files if f.endswith('.nc')]
print(file)
if len(file)>1:
    if OBSID=='MOSAIC': # run first 0 then 1
        #files = [f"{inpath}/{file[0]}",f"{inpath}/{file[1]}"] 
        #ds = xr.open_mfdataset(files, combine='by_coords', parallel=True)
        infile = f"{inpath}/{file[1]}"
        Functions.netcdf_to_txt(infile, infile.replace('.nc', '.dat')) # convert to .dat file
        obsfile = infile.replace('.nc', '.dat')
    if OBSID=='Nansen_legacy': # run first 0 then 1
        infile = f"{inpath}/{file[1]}"
        Functions.netcdf_to_txt(infile, infile.replace('.nc', '.dat')) # convert to .dat file
        obsfile = infile.replace('.nc', '.dat')
    if OBSID=='NICE':
        infile = f"{inpath}/{file[1]}"
        Functions.netcdf_to_txt(infile, infile.replace('.nc', '.dat')) # convert to .dat file
        obsfile = infile.replace('.nc', '.dat')

else:
    infile = f"{inpath}/{file[0]}" 
    Functions.netcdf_to_txt(infile, infile.replace('.nc', '.dat')) # convert to .dat file
    obsfile = infile.replace('.nc', '.dat')

print(obsfile)

if 'SID' in file[0]: #OBSID=='SCICEX' or OBSID=='BGEP' or OBSID=='TRANSDRIFT' or OBSID=='AWI-ULS' or OBSID=='NPEO' or OBSID=='NPI':
    var = ['SID']
else:
    var = [ 'SIT', 'SD', 'FRB'] # list of variables relevant for campaign (either ['SID'] or [ 'SIT', 'SD', 'FRB'])
#---------------------------------

## Read reference observations
ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')


# find relevant years and months in obsfile
ObsDate = ObsData['date']
formating = "%Y-%m-%dT%H:%M:%S"
ObsDate_object = [dt.datetime.strptime(Date, formating) for Date in ObsDate]
ObsYears = [ObsDate.year for ObsDate in ObsDate_object]
ObsMonth = [ObsDate.month for ObsDate in ObsDate_object]

# Location of output files
#save_path =  os.path.dirname(os.getcwd()) + '/satellite/Satellite_subsets/' + OBSID
save_path =  f'/dmidata/projects/cmems2/C3S/RRDPp/satellite/Satellite_subsets/{OBSID}'
if not os.path.exists(save_path):os.makedirs(save_path)
print('------------------')
print(save_path)

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
    
    count = 0
    check_count = 0; ## to check for file size
    ## Loop through all relevant dates in satellite file
    for t,i in zip(time, range(len(time))):
        if t.year in ObsYears:
            if t.month in ObsMonth:
                ofile, check_count = check_filesize(OBSID,SATELLITE, check_count, save_path)
                # Extract satellite subset
                collocation(OBSID, obsfile, ifile, ofile, count, i)
                count += 1


else:  # if satellite is ENV or CS2
    # Satellite data directory
    # sat_dir = os.path.dirname(os.getcwd()) + "/satellite/ENV_CS2_data/"+SATELLITE.replace('-','')+"_data/" + HS + "/"
    sat_dir = f"/dmidata/projects/cmems2/C3S/RRDPp/satellite/ENV_CS2_data/"+SATELLITE.replace('-','')+"_data/" + HS + "/"
    ## Loop through all satelitte files, with relevant dates
    years = os.listdir(sat_dir)
    count = 0
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
                    files = sorted(os.listdir(sat_dir + year + '/' + month))
                    for file in files:
                        #print(file)
                        day = file[-14:-12]
                        
                        if HS=='SH':
                            ofile, check_count = check_filesize(OBSID,SATELLITE, check_count, save_path, HS=HS)
                        else:
                            ofile, check_count = check_filesize(OBSID,SATELLITE, check_count, save_path)
                        
                        ifile = os.path.join(sat_dir + year + '/' + month, file)
                        # Extract satellite subset
                        collocation(OBSID, obsfile, ifile, ofile, count, save_path)
                        count += 1


                
