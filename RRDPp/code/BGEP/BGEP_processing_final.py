# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 13:35:37 2021

@author: Olsen
"""

"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured during the Beaufort Gyre Exploration Project
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2021-10-22'

# -- Built-in modules -- #
import os.path
import datetime as dt
import glob
import re
import sys

# -- Third-part modules -- #
import numpy as np
import pandas as pd

# -- Proprietary modules -- #
from Warren import SnowDepth, SWE
sys.path.append(os.path.dirname(os.getcwd()))
import Functions

dtint = 30  # days
gridres = 25000  # m

directory = os.path.dirname(os.path.dirname(os.getcwd())) + '/RawData/BPR_ULS_BGEP'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/fig/'
save_path_data= os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/final/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/BGEP/final/ESACCIplus-SEAICE-RRDP2+-SID-BGEP-V3.dat'


files = glob.glob(save_path_data+'*')
# for f in files:
#     os.remove(f)

lat=0
lon=0

count = 0
for dir in os.listdir(directory):
    print(dir)
    dir_data=os.path.join(directory, dir)
    if os.path.isdir(dir_data) and not dir_data.endswith('no_use'):        
        for ifile in os.listdir(dir_data):
            file=os.path.join(dir_data, ifile)
            #defines variables in output file
            count+=1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
            dataOut.pp_flag = 0
            with open(file,'rb') as myFile:
                for num, line in enumerate(myFile, 1):
                    print(line.strip())
                    if num == 1:
                        line = line.strip()
                        line = line.split()
                        print(line)
                        dataOut.obsID = 'BGEP_Mooring' + line[3].decode('utf8')[:-1]
                        dataOut.lat = float(line[4].decode('utf8')) + float(line[5].decode('utf8'))/60 # convert to deicmal degrees
                        if line[9].decode('utf8') == 'W':
                            dataOut.lon = - float(line[7].decode('utf8')) - float(line[8].decode('utf8'))/60
                        elif line[9].decode('utf8') == 'E':
                            dataOut.lon = float(line[7].decode('utf8')) + float(line[8].decode('utf8'))/60
                        break
 
            # Reads observation data from ASCII-file
            df = pd.read_table(file, skiprows=1, sep="\s+", dtype = 'float32', converters = {'%date': str})
                        
            # Load data
            dates = df['%date'].to_numpy()
            time = df['time(UTC)'].to_numpy() #UTC - duration of each time used on each measurement
            SID = df['draft(m)'].to_numpy() #m
            SID_Unc = np.array([0.10 for sid in SID]) # uncertainty of 10 cm approx.
            print('Getting months')
            
            months = [int(date[4:6]) for date in dates]
            days = [int(date[6:]) for date in dates]
            mondiff=np.where(~(np.diff(months) == 0))[0] #index of when month changes
            index=np.insert(np.append(mondiff,len(months)-1),0,0) #add start and end indexes on month
            time_in = [dt.datetime(int(y[:4]),int(m),int(d)) for y,m,d in zip(dates, months, days)]
            
            # Avg_draft
            dataOut.SID_final = np.array([np.mean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.SID_std = np.array([np.std(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.SID_ln = np.array([len(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            
            uncSID = []
            for i in range(len(index)-1): # calculate uncertainty 
                start = index[i]
                end = index[i+1]
                uncSID = np.append(uncSID, 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2)))
            
            iddd = np.array([index[i] for i in range(len(index)-1)])
            # Avg date
            mean_day = np.array([np.median(days[index[i]:index[i+1]]) for i in range(len(index)-1)])
            
            print('Created SID variables')


            # Correlates BGEP data with Warren snow depth and snow density
            mon_uni = [months[ind] for ind in index[1:]]
            for month in mon_uni:
                (w_SD,w_SD_epsilon) = SnowDepth(lat,lon,month)
                dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
                (wswe,wswe_epsilon) = SWE(lat,lon,month)
                w_density=int((wswe/w_SD)*1000)
                dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)

            print('Warren climatology created')
            months_final = ['0'+str(mon_uni[i]) if mon_uni[i]<10 else str(mon_uni[i]) for i in range(len(mon_uni))]
            days_final = ['0'+str(int(d)) if int(d)<10 else str(int(d)) for d in mean_day]
            years = [dates[index[1:]][i][:4] for i in range(len(index[1:]))]
            dataOut.date_final=[years[i] +'-'+ months_final[i] +'-'+ days_final[i] + 'T00:00:00' for i in range(len(years))]
            dataOut.SID_unc = uncSID
            dataOut.lat_final = [dataOut.lat for i in range(len(dataOut.date_final))]
            dataOut.lon_final = [dataOut.lon for i in range(len(dataOut.date_final))]

            time_out = [dt.datetime(int(y),int(m),int(d)) for y,m,d in zip(years, months_final, mean_day)]
            if len(time_out) > 0:
                # Functions.plot(dataOut.lat, dataOut.lon, dataOut.obsID, time_out,saveplot, HS='NH')
                Functions.scatter(dataOut.obsID + '_' + str(time_in[0].year), time_in[::1000], SID[::1000], time_out, dataOut.SID_final, 'SID [m]', saveplot)
            
            # fill empty arrays with NaN values
            dataOut.Check_Output()
                
            # print data to output file
            dataOut.Print_to_output(ofile, primary='SID')
    
Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH', primary='SID')
