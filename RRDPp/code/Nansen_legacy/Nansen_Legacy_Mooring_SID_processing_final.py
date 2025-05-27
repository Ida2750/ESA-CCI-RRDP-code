# -*- coding: utf-8 -*-
"""
Uses EASE-grid to produce 25 km grid mean values. Uses input files measured by upward looking sonars in the Laptev Sea from the Nansen Legacy program
includes Warren snow depths and densities

# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2025-04-30'
"""

# -- Built-in modules -- #
import os.path
import sys

# -- Third-part modules -- #
import numpy as np
import xarray as xr
import datetime as dt

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions
import EASEgrid_correct as EASEgrid
from Warren import SnowDepth, SWE

# %% Main
dtint = 30
gridres = 25000
name = 'Nansen_legacy'

# Reads input data
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + f'/FINAL/{name}/final/'
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + f'/FINAL/{name}/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd(
))) + f'/FINAL/{name}/final/ESACCIplus-SEAICE-RRDP2+-SID-{name}.nc'
directory = os.path.dirname(os.path.dirname(os.path.dirname(
    os.getcwd()))) + f'/RRDPp/RawData/{name}'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

# read Moorings
directory = f'{directory}/Mooring'
count = 0 # count for header

filenames = sorted(os.listdir(directory))
filenames = [filename for filename in filenames if filename.endswith('.nc')]
print(filenames)

def robust_std(data):
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median))
    return mad * 1.4826  # scaling factor for normal distribution

datalen = 0
for filename in filenames:

    ifile = os.path.join(directory, filename)
    data = xr.open_dataset(ifile)
    print(data)

    #defines variables in output file
    count += 1
    dataOut = Functions.Final_Data(Type='SID', count_head=count)
    dataOut.pp_flag = 0
    dataOut.unc_flag = 2

    ## define obsID
    dataOut.obsID = f'{name}_' + filename[:12]
    print(dataOut.obsID)

    # Load data
    dates = data['TIME'].to_numpy().astype(str)
    # round date to seconds (contains nanoseconds)
    dates = [d.split('.')[0] for d in dates]
    latitude = data['LATITUDE'].to_numpy()  # degrees
    longitude = data['LONGITUDE'].to_numpy()  # degrees
    SID = data['SEA_ICE_DRAFT_MEDIAN'].to_numpy()  # meters

    # uncertainty of 5 cm approx.
    SID_Unc = np.array([np.nan for sid in SID])

    # set negative to nan
    SID[SID<0] = np.nan
    SID[SID>8] = np.nan
    
    datalen += len(SID[np.isfinite(SID)])
    # calculates average, number of data, std and unc of each month
    t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                    for s in dates])

    t2 = np.array(
        [(time-dt.datetime(1970, 1, 1)).total_seconds() for time in t])

    months = [time.month for time in t]
    mondiff = np.where(~(np.diff(months) == 0))[0]  # index of when month changes
    # add start and end indexes on month
    if mondiff[0]!=0:
        index = np.insert(np.append(mondiff, len(months)-1), 0, 0)
    
    # compute monthly SID values
    dataOut.SID_final=np.array([np.nanmedian(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    #dataOut.SID_std=np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
    dataOut.SID_std = np.array([robust_std(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])

    dataOut.SID_ln = np.array([len(SID[index[i]:index[i+1]])
                                for i in range(len(index)-1)])
    # compute uncertainty
    for i in range(len(index)-1):
        start = index[i]
        end = index[i+1]
        dataOut.SID_unc = np.append(
            dataOut.SID_unc, 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2)))

    #find median date within month (tells about which part of the month majority measurements are from)
    avgDates = np.array([np.median(t2[index[i]:index[i+1]])
                        for i in range(len(index)-1)])

    #change date to right format
    dates_final = np.array(
        [dt.datetime.fromtimestamp(int(sec)) for sec in avgDates])
    dataOut.date_final = np.array([dt.datetime.strftime(
        date, "%Y-%m-%dT%H:%M:%S") for date in dates_final])

            
    ######### QUALITY FLAGS ############
    dataOut.QFT = [] # temporal
    dataOut.QFS = [] # spatial
    dataOut.QFG = [] # global threshold

    days = [time.day for time in t]
    for i in range(len(index)-1):
        # find number of days
        unique = len(np.unique(days[index[i]:index[i+1]]))
        if np.any(SID[index[i]:index[i+1]]>6):
            dataOut.QFG.append(1)
        else:
            dataOut.QFG.append(0)
        if unique==1:
            dataOut.QFT.append(3)
            dataOut.QFS.append(3)
        elif unique<=5:
            dataOut.QFT.append(2)
            dataOut.QFS.append(3)
        elif unique<15:
            dataOut.QFT.append(1)
            dataOut.QFS.append(3)
        elif unique>=15:
            dataOut.QFT.append(0)
            dataOut.QFS.append(0)
    
    dataOut.QFT = np.array(dataOut.QFT)
    dataOut.QFS = np.array(dataOut.QFS)
    dataOut.QFG = np.array(dataOut.QFG)
    ######### QUALITY FLAGS ############    
    
    # get latitude and longtitude information
    dataOut.lat_final = [latitude] * len(dataOut.SID_final)
    dataOut.lon_final = [longitude] * len(dataOut.SID_final)

    # plot data before and after processing
    if len(dates_final) > 0:
        Functions.scatter(dataOut.obsID, t, SID, dates_final,
                            dataOut.SID_final, 'SID [m]', saveplot)

    # Correlate {name} data with Warren snow depth and snow density
    for ll in range(np.size(dataOut.SID_final, 0)):
        (w_SD, w_SD_epsilon) = SnowDepth(
            dataOut.lat_final[ll], dataOut.lon_final[ll], dates_final[ll].month)
        dataOut.w_SD_final = np.append(dataOut.w_SD_final, w_SD)
        (wswe, wswe_epsilon) = SWE(
            dataOut.lat_final[ll], dataOut.lon_final[ll], dates_final[ll].month)
        w_density = int((wswe/w_SD)*1000)
        dataOut.w_density_final = np.append(
            dataOut.w_density_final, w_density)

    dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
    dataOut.pp_flag = [dataOut.pp_flag]*len(dataOut.SID_final)
    dataOut.unc_flag = [2]*len(dataOut.SID_final)
    dataOut.obsID = [dataOut.obsID]*len(dataOut.SID_final)
    
    # fill empty arrays with NaN values
    dataOut.Check_Output()

    if count>1:
        subset = dataOut.Create_NC_file(ofile, primary='SID')
        df = Functions.Append_to_NC(df, subset)
    else:
        df = dataOut.Create_NC_file(ofile, primary='SID', datasource='Sea ice draft from mooring observations in the northwestern Barents Sea from 2018 onward, Nansen Legacy, DOI: 10.21334/npolar.2024.c19e8a7d', key_variables='Sea Ice Draft')
        

# Sort final data based on date
print(datalen)
Functions.save_NC_file(df, ofile, primary='SID')
