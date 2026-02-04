# -*- coding: utf-8 -*-
"""
Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by upward looking sonars in the Laptev Sea from the Russian-German TRANSDRIFT program
includes Warren snow depths and densities

# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2023-06-12'
"""

# -- Built-in modules -- #
import os.path
import sys

# -- Third-part modules -- #
import numpy as np
import datetime as dt

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions
import EASEgrid_correct as EASEgrid
from Warren import SnowDepth, SWE

# %% Main
dtint = 30
gridres = 25000
# Reads input data
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/TRANSDRIFT/final/'
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/TRANSDRIFT/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd(
))) + '/FINAL/TRANSDRIFT/final/ESACCIplus-SEAICE-RRDP2+-SID-TRANSDRIFT.nc'
directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.getcwd())))) + '/RRDPp/RawData/TRANSDRIFT/datasets'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

total_obs_valid_SID_ULS = 0
total_obs_valid_SID_ACDP = 0
count = 0 # count for header
for filename in os.listdir(directory):
    if filename.endswith('.tab'):
        #try:
            ifile = os.path.join(directory, filename)
            file = open(ifile, 'r')
    
            # Find number of headerlines
            lookup = '*/'
            #with open(ifile,'rb') as myFile:
            for num, line in enumerate(file, 1):
                # print(line.strip())
                if line.strip() == lookup:
                    skip_header = num
                    break
    
            # Reads observation data from ASCII-file
            data = np.genfromtxt(file, dtype=None, skip_header=skip_header, names=[
                                 'date', 'lat', 'lon', 'SID'], delimiter='\t')
    
            #defines variables in output file
            count += 1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
            dataOut.pp_flag = 0
            dataOut.unc_flag = 2
    
            ## newer files
            if filename.startswith('ULS_1893'):
                dataOut.obsID = filename[0:13]
            elif filename.startswith('ULS_Taymyr'):
                dataOut.obsID = filename[0:15]
            else:
                dataOut.obsID = filename.split('_')[0]
    
            # Load data
            date = data['date']
            dates = [dat.decode('utf8') for dat in date]
            latitude = data['lat']  # degrees
            longitude = data['lon']  # degrees
            SID = data['SID']  # meters
    
            # uncertainty of 10 cm approx.
            if filename.startswith('ULS_1893') or filename.startswith('ULS_Taymyr'):
                SID_Unc = np.array([0.05 for sid in SID])
                
                total_obs_valid_SID_ACDP += len(SID[np.isfinite(SID)])
                
                # set non existing data to nan
                SID[SID == -99999.0] = np.nan
            else: # ACDP's - these have a higher uncertainty than the ULS's
                SID_Unc = np.array([0.96 for sid in SID])
    
                # set non existing data to nan
                SID[SID == -99999.0] = np.nan

                total_obs_valid_SID_ULS += len(SID[np.isfinite(SID)])

print(total_obs_valid_SID_ULS)
print(total_obs_valid_SID_ACDP)
print(total_obs_valid_SID_ULS + total_obs_valid_SID_ACDP)
