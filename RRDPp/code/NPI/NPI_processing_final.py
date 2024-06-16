#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data from NPI monthly mean SID values for the period 1990-2018
No data errors provided!

No information is provided regarding number of datapoints, therefore 
the calculated uncertainty is an upper bound and is likely an overestimation!
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2021-08-12'

# -- Built-in modules -- #
import os.path
import datetime as dt
from datetime import timedelta as td
import sys

# -- Third-part modules -- #
import numpy as np
import datetime as dt
import netCDF4 as nc

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions
from Warren import SnowDepth, SWE

# %% Main
dtint = 30
gridres = 25000
# Reads input data
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPI/final/'
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPI/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd(
))) + '/FINAL/NPI/final/ESACCIplus-SEAICE-RRDP2+-SID-NPI-V3.dat'
directory = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/RawData/NPI_ULS_FS/ULS_ncfiles_incl_README/'  

count=0
for filename in os.listdir(directory):
    if filename.endswith(".nc"):
        # print(filename)
        fn=directory+filename
        ds = nc.Dataset(fn)

        #defines variables in output file
        count += 1
        dataOut = Functions.Final_Data(Type='FRB', count_head=count)
        dataOut.pp_flag = 0
        dataOut.unc_flag = 1
        
        # for var in ds.variables.values():
        #     count1+=1

        time = ds['TIME'][:] #time since 1950-01-01T00:00:00Z
        ID_mooring=ds['PLATFORM'][:]
        #[np.argsort(time,axis=0)]
        ID=ID_mooring[0].decode('utf8')+ID_mooring[1].decode('utf8')+ID_mooring[2].decode('utf8')
        dataOut.lat_final = np.concatenate([ds['LATITUDE'][:] for t in time])
        dataOut.lon_final = np.concatenate([ds['LONGITUDE'][:] for t in time])
        dataOut.SID_final = ds['DRAFT'][:][np.argsort(time,axis=0)]
        
        dataOut.SID_ln = [np.nan for s in dataOut.SID_final]
        dataOut.SID_std = [np.nan for s in dataOut.SID_final]
        
        entries_usefull=ds['DATA_COVERAGE'][:]
        entries_usefull=entries_usefull[0][np.argsort(time,axis=0)] #decimal
          
        date0=dt.date(1950,1,1)
        dates=[date0+td(seconds=t) for t in time[np.argsort(time,axis=0)]] #list comprehension more efficient memory wise
        
        for date in dates:
            if date.year < 2005:
                dataOut.SID_unc = np.append(dataOut.SID_unc, 0.20)  # 20 cm prior to 2005
            elif date.year >= 2005:
                dataOut.SID_unc = np.append(dataOut.SID_unc, 0.10)  # 10 cm post 2005


                
        # Correlate NPI data with Warren snow depth and snow density
        for ll in range(np.size(dataOut.SID_final, 0)):
            (w_SD, w_SD_epsilon) = SnowDepth(
                dataOut.lat_final[ll], dataOut.lon_final[ll], dates[ll].month)
            dataOut.w_SD_final = np.append(dataOut.w_SD_final, w_SD)
            (wswe, wswe_epsilon) = SWE(
                dataOut.lat_final[ll], dataOut.lon_final[ll], dates[ll].month)
            w_density = int((wswe/w_SD)*1000)
            dataOut.w_density_final = np.append(
                dataOut.w_density_final, w_density)

        dataOut.obsID = 'NPI_' + ID # + '_' + str(dates[0]).replace('-','')
        dataOut.date_final=np.array([dt.datetime.strftime(date,"%Y-%m-%dT%H:%M:%S") for date in dates])
        # fill empty arrays with NaN values
        dataOut.Check_Output()

        # print data to output file
        dataOut.Print_to_output(ofile, primary='SID')

# Sort final data based on date
Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH', primary='SID')

        
