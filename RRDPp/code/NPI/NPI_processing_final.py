#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data is from Norwegian Polar Institue (NPI). https://doi.org/10.21334/npolar.2022.b94cb848 
providing monthly mean SID values for the period 1990-2018.

No information is provided regarding number of datapoints, therefore 
the calculated uncertainty is an upper bound uncertainty
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2024-08-03'

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
dtint = 30 # days per mean
gridres = 25000 # grid resolution
# Reads input data
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPI/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPI/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd(
))) + '/FINAL/NPI/final/ESACCIplus-SEAICE-RRDP2+-SID-NPI.nc'


directory = os.path.dirname(os.path.dirname(os.path.dirname(
    os.getcwd()))) + '/RRDPp/RawData/NPI_ULS_FS/ULS_ncfiles_incl_README/'  

import xarray as xr

count=0
for filename in os.listdir(directory):
    if filename.endswith(".nc"):
        print(filename)
        # print(filename)
        fn=directory+filename
        #ds = xr.open_dataset(fn)
        ds = nc.Dataset(fn)
        
        #print(ds)
        #defines variables in output file
        count += 1
        dataOut = Functions.Final_Data(Type='SID', count_head=count)

        
        # for var in ds.variables.values():
        #     count1+=1

        time = ds['TIME'][:] #time since 1950-01-01T00:00:00Z
        #print(time)
        ID_mooring=ds['PLATFORM'][:]
        #[np.argsort(time,axis=0)]
        ID=ID_mooring[0].decode('utf8')+ID_mooring[1].decode('utf8')+ID_mooring[2].decode('utf8')
        dataOut.lat_final = np.concatenate([ds['LATITUDE'][:] for t in time])
        dataOut.lon_final = np.concatenate([ds['LONGITUDE'][:] for t in time])
        dataOut.SID_final = ds['DRAFT'][:][np.argsort(time,axis=0)]
        
        #dataOut.SID_ln = [np.nan for s in dataOut.SID_final]
        dataOut.SID_std = [np.nan for s in dataOut.SID_final]
        
        entries_usefull=ds['DATA_COVERAGE'][:]
        entries_usefull=entries_usefull[0][np.argsort(time,axis=0)] #decimal
        
          
        date0=dt.date(1950,1,1)
        dates=[date0+td(seconds=t) for t in time[np.argsort(time,axis=0)]] #list comprehension more efficient memory wise
        #date0 = []
        for date in dates:
            # source: https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-29470-7/MediaObjects/41467_2022_29470_MOESM1_ESM.pdf
            if 'ES300' in filename:
                dataOut.SID_ln = np.append(dataOut.SID_ln, 1e4)
                if date.year < 1992:
                    dataOut.SID_unc = np.append(dataOut.SID_unc, 2.7e-3)

                elif date.year < 2003 and date.year>=1992:
                    dataOut.SID_unc = np.append(dataOut.SID_unc, 1.9e-3)

                elif date.year < 2007:
                    dataOut.SID_unc = np.append(dataOut.SID_unc, 1.3e-3)
                

            elif 'IPS' in filename: # newer instuments used IPS4 or IPS5
                dataOut.SID_unc = np.append(dataOut.SID_unc, 0.08e-3)  # 10 cm post 2005
                dataOut.SID_ln = np.append(dataOut.SID_ln, 1e6)

        dataOut.pp_flag = [1]*len(dataOut.SID_final)
        dataOut.unc_flag = [1]*len(dataOut.SID_final)
        dataOut.QFT = [np.nan]*len(dataOut.SID_final)
        dataOut.QFS = [np.nan]*len(dataOut.SID_final)
        dataOut.QFG = [np.nan]*len(dataOut.SID_final)
                
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

        dataOut.obsID = ['NPI_' + ID]*len(dataOut.SID_final) # + '_' + str(dates[0]).replace('-','')
        #print(dataOut.obsID)
        dataOut.time = [np.datetime64(d) for d in dates]

        dataOut.date_final=np.array([dt.datetime.strftime(date,"%Y-%m-%dT%H:%M:%S") for date in dates])
        # fill empty arrays with NaN values
        dataOut.Check_Output()

        # print data to output file
        #dataOut.Print_to_output(ofile, primary='SID')
        
        
        if count>1:
            subset = dataOut.Create_NC_file(ofile, primary='SID')
            df = Functions.Append_to_NC(df, subset)
        else:
            df = dataOut.Create_NC_file(ofile, primary='SID', datasource='Monthly mean sea ice draft from the Fram Strait Arctic Outflow Observatory since 1990, doi: 10.21334/npolar.2021.5b717274', key_variables='Sea Ice Draft', comment='The number of observations in SID_ln is an approximated number, see section 3.2 in CCI SIT RRDP associated publication')
            

# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SID')
#Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH', primary='SID')

        