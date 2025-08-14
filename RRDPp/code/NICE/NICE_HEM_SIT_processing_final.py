# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:24:45 2023

@author: Ida Olsen

Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Ice mass buoys from the The CRREL- Dartmouth Mass Balance Buoys Program
includes Warren snow depths and densities
"""
# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '1'
__date__ = '2023-06-12'

# -- Built-in modules -- #
import os.path
import sys
import glob
import re
import datetime as dt

# -- Third-part modules -- #
import numpy as np
from numpy.lib import recfunctions as rfn
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
#import PyPDF2
import pandas as pd
import csv

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions


#%% Support functions
def Get_ppflag(self, SD, lnSD, lnSIT):
    """
    Assign pre-processing flag, where
    the data gets pp-flag=2 if any position used in the average
    was obtained with a temporal discrepancy of more than 12 hours
    """
    # find valid SD entries
    pp_non_nan = np.array(self.pp_flag)[~np.isnan(SD)]
    # print(pp_non_nan)
    pp_flag = np.zeros(len(lnSD))
    summ = np.cumsum(lnSD).astype(int)
    try:
        pp_flag[0] = np.nanmax(pp_non_nan[:summ[0]])
        # print(np.max(pp_non_nan[:summ[0]]))
    except: # no data in datapoint
        pp_flag[0] = np.nan
    # find max pre-processing flag within the interval used for the average
    for i in range(1, len(summ)):
        # print(summ[i])
        try:
            pp_flag[i] = np.nanmax(pp_non_nan[summ[i-1]:summ[i]])
        except: # no data in datapoint
            pp_flag[i] = np.nan
    self.pp_flag = pp_flag
    print(pp_flag)
    return self.pp_flag


#%% Main

# Information
name = 'NICE'
save_path = os.path.dirname(os.path.dirname(os.getcwd())) + f'/{name}/'
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + f'/FINAL/{name}/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
ofile = f'ESACCIplus-SEAICE-RRDP2+-SIT-{name}-HEM_v2.nc'
gridres = 25000  # grid resolution
dtint = 30  # days per mean
## saving locations
ofile = os.path.join(save_path_data, ofile)


saveplot = os.path.join(os.path.dirname(
    os.path.dirname(os.getcwd())), f'FINAL/{name}/fig/')
# create directory if they do not exist
if not os.path.exists(saveplot): os.makedirs(saveplot)
files = glob.glob(f'/dmidata/users/ilo/projects/RRDPp/RawData/{name}/HEM_v2/NICE2015-HEM-V2/2*.txt')

datalen = 0
count = 0  # used to locate first file
# Prepare empty lists to accumulate all files' data
all_lat = []
all_lon = []
all_SIT = []
all_time = []

#defines variables in output file
count+=1
dataOut = Functions.Final_Data(Type='SIT', count_head=count)
dataOut.pp_flag = 0
for ifile in files:
    Snow=False #default as most files does not contain snowdepth measurements
    file = os.path.basename(ifile)
    print(file)
    
    data = np.genfromtxt(ifile, dtype=None, delimiter='\t', names=True, encoding=None)

    latitude = data['Lat']
    longitude = data['Lon']
    SIT = data['Thick_m']
    # date is provided as year, month, day + seconds since start of day
    date = data['Date']
     
    # replace faulty values if applicable
    SIT[SIT < 0] = np.nan
    SIT[SIT > 10] = np.nan

    datalen += len(SIT[np.isfinite(SIT)])

    # convert date to datetime
    t = [dt.datetime.strptime(d, "%Y-%m-%dT%H:%M:%S.%fZ") for d in date]

    # Append current file's data to master lists
    all_lat.append(latitude)
    all_lon.append(longitude)
    all_SIT.append(SIT)
    all_time.append(t)

# Combine lists into final arrays
latitude = np.concatenate(all_lat)
longitude = np.concatenate(all_lon)
SIT = np.concatenate(all_SIT)
t = np.concatenate([np.array(times) for times in all_time])

## Define fixed uncertainty
SIT_unc = np.ones(len(latitude)) * 0.10

# Define obsID
dataOut.obsID = f'N-ICE2015_HEMv2'

# Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
G = EASEgrid.Gridded()
G.SetHemisphere('N')
G.CreateGrids(gridres)
(index_i, index_j) = G.LatLonToIdx(latitude, longitude)

# Takes the time for each grid cell into account and calculate averages
avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, var1, var2, dataOut.QFT, dataOut.QFS, dataOut.QFG = G.GridData(
    dtint, latitude, longitude, t, SIT=SIT, SIT_unc=SIT_unc, dtype='AEM')

dataOut.unc_flag = np.ones(len(avgSIT)).astype(int)*2
dataOut.unc_flag[avgSIT>3] = 3
    
if len(time) > 0:
    Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS='NH')
    Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]', saveplot)
    if Snow:
        Functions.scatter(dataOut.obsID, t, SD, time, avgSD, 'SD [m]', saveplot)


# Correlates HEM data with Warren snow depth and snow density
for ll in range(np.size(avgSD,0)):
    (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
    dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
    (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
    w_density=int((wswe/w_SD)*1000)
    dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)

#Change names to correct format names
dataOut.obsID = [dataOut.obsID]*len(lat)
## pp flag + unc flag
dataOut.pp_flag = [dataOut.pp_flag]*len(lat)
dataOut.lat_final = lat
dataOut.lon_final = lon
for ll in range(np.size(time,0)):
    dataOut.date_final = np.append(dataOut.date_final,dt.datetime.strftime(time[ll],"%Y-%m-%dT%H:%M:%S"))
dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
dataOut.SD_final = avgSD
dataOut.SD_std = stdSD
dataOut.SD_ln = lnSD
dataOut.SD_unc = uncSD
dataOut.SIT_final = avgSIT
dataOut.SIT_std = stdSIT
dataOut.SIT_ln = lnSIT
dataOut.SIT_unc = uncSIT

# fill empty arrays with NaN values
dataOut.Check_Output()
    
if count>1:
    subset = dataOut.Create_NC_file(ofile,primary='SIT')
    df = Functions.Append_to_NC(df, subset)
else:
    df = dataOut.Create_NC_file(ofile,primary='SIT', datasource='Airborne Electromagnetic Measurement, N-ICE2015, DOI:https://doi.org/10.21334/NPOLAR.2016.70352512', key_variables='Ice thickness: sea ice thickness + snow depth')

print(datalen)
# Save data to NetCDF
Functions.save_NC_file(df, ofile, primary='SIT')

    
