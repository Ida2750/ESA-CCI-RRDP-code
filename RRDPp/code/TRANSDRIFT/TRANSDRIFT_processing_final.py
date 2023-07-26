# -*- coding: utf-8 -*-
Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by upward looking sonars in the Laptev Sea from the Russian-German TRANSDRIFT program
includes Warren snow depths and densities
"""
# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
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
))) + '/FINAL/ULS_RUSSIA/final/ESACCIplus-SEAICE-RRDP2+-SID-TRANSDRIFTdat'
directory = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/RawData/TRANSDRIFT/datasets'


count = 0 # count for header
for filename in os.listdir(directory):
    if filename.endswith('.tab'):

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
        dataOut = Functions.Final_Data(Type='FRB', count_head=count)
        dataOut.pp_flag = 0
        dataOut.unc_flag = 2

        ## newer files
        if filename.startswith('ULS_1893'):
            dataOut.obsID = filename[0:13]
        else:
           dataOut.obsID = filename[0:15]

        # Load data
        date = data['date']
        dates = [dat.decode('utf8') for dat in date]
        latitude = data['lat']  # degrees
        longitude = data['lon']  # degrees
        SID = data['SID']  # meters

        # uncertainty of 10 cm approx.
        SID_Unc = np.array([0.10 for sid in SID])

        # set non existing data to nan
        SID[SID == -99999.0] = np.nan

        # calculates average, number of data, std and unc of each month
        t = np.array([dt.datetime.strptime(s, "%Y-%m-%d")
                      for s in dates])

        t2 = np.array(
            [(time-dt.datetime(1970, 1, 1)).total_seconds() for time in t])

        months = [time.month for time in t]
        mondiff = np.where(~(np.diff(months) == 0))[0]  # index of when month changes
        # add start and end indexes on month
        index = np.insert(np.append(mondiff, len(months)-1), 0, 0)
        
        # compute monthly SID values
        dataOut.SID_final = np.array(
            [np.nanmean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
        dataOut.SID_std = np.array(
            [np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
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
        
        # get latitude and longtitude information
        dataOut.lat_final = latitude[:len(dataOut.SID_final)]
        dataOut.lon_final = longitude[:len(dataOut.SID_final)]

        # plot data before and after processing
        if len(dates_final) > 0:
            Functions.scatter(dataOut.obsID, t, SID, dates_final,
                              dataOut.SID_final, 'SID [m]', saveplot)

        # Correlate TRANSDRIFT data with Warren snow depth and snow density
        for ll in range(np.size(dataOut.SID_final, 0)):
            (w_SD, w_SD_epsilon) = SnowDepth(
                dataOut.lat_final[ll], dataOut.lon_final[ll], dates_final[ll].month)
            dataOut.w_SD_final = np.append(dataOut.w_SD_final, w_SD)
            (wswe, wswe_epsilon) = SWE(
                dataOut.lat_final[ll], dataOut.lon_final[ll], dates_final[ll].month)
            w_density = int((wswe/w_SD)*1000)
            dataOut.w_density_final = np.append(
                dataOut.w_density_final, w_density)

        # fill empty arrays with NaN values
        dataOut.Check_Output()

        # print data to output file
        dataOut.Print_to_output(ofile, primary='SID')

# Sort final data based on date
Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH', primary='SID')
