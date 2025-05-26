#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 25 km gridded means of input files containing measurements of SIT and SD from 
airborne campaigns using the EM-bird, and from the Snow and Ice Mass Balance Array both
obtained during the MOSAIC campaign
includes Warren snow depths and densities
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2025-04-30'

# -- Built-in modules -- #
import os.path
import datetime as dt
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions

#%%
# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean
name = 'Nansen_legacy'

# Define locations of input and output files
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/{name}/final/'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) +  f'/FINAL/{name}/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/{name}/final/ESACCIplus-SEAICE-RRDP2+-SIT-{name}.nc'
#directory = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) + f'/RRDPp/RawData/{name}'
directory = f'/dmidata/projects/cmems2/C3S/RRDPp/RawData/{name}'

# create directories if they do not exist
if not os.path.exists(save_path_data): os.makedirs(save_path_data)
if not os.path.exists(saveplot): os.makedirs(saveplot)

count = 0  # used to locate first file

# Access HEM data
directories = sorted(os.listdir(f'{directory}/HEM'))
print(directories)

datalen = 0
for dirr in directories:
    
    files = sorted(os.listdir(f'{directory}/HEM/{dirr}'))
    print(files)
    for ifile in files:
        Snow=False #default as most files does not contain snowdepth measurements
        file = os.path.join(f'{directory}/HEM/{dirr}/{ifile}')
        # Access all SIT files in AEM_AWI
        if file.endswith('.csv'):
            
            #defines variables in output file
            count+=1
            dataOut = Functions.Final_Data(Type='SIT', count_head=count)
            dataOut.pp_flag = 0
                
            data = np.genfromtxt(file, dtype=None, delimiter=',', names=True, encoding=None)

            print(data.dtype.names)

            try:
                latitude = data['latitude_degree']
                longitude = data['longitude_degree']
                SIT = data['total_thickness_m']
                # date is provided as year, month, day + seconds since start of day
                year = data['Year']
                month = data['month']
                day = data['day']
                seconds = data['seconds_since_the_day_start']
            except:
                latitude = data['Latitude_degree']
                longitude = data['Longitude_degree']
                SIT = data['Total_thickness_m']
                # date is provided as year, month, day + seconds since start of day
                year = data['Year']
                month = data['Month']
                day = data['Day']
                seconds = data['secs_since_the_day_start']      

            ## Define fixed uncertainty
            SIT_unc = np.ones(len(latitude)) * 0.10

            # replace faulty values if applicable
            SIT[SIT < 0] = np.nan
            SIT[SIT > 10] = np.nan

            datalen += len(SIT[np.isfinite(SIT)])

            # combine time information
            t = np.array([
                dt.datetime(y, m, d) + dt.timedelta(seconds=s)
                for y, m, d, s in zip(year, month, day, seconds)
            ])

            # Define obsID
            if 'JC' in dirr:
                flight = ifile[9:17]
            else:
                flight = ifile[7:15]

            dataOut.obsID = f'{name}_' + flight

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


            # Correlates SB-AWI buoy data with Warren snow depth and snow density
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
                df = dataOut.Create_NC_file(ofile,primary='SIT', datasource='Airborne Electromagnetic Measurement, Nansen Legacy, DOI: 10.21334/npolar.2023.1a9cc2df + 10.21334/npolar.2023.c1cfd5dd', key_variables='Ice thickness: sea ice thickness + snow depth')

    print(datalen)
    # Save data to NetCDF
    Functions.save_NC_file(df, ofile, primary='SIT')

            
