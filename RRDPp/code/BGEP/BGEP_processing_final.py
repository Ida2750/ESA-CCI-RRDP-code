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
__date__ = '2024-08-12'

# -- Built-in modules -- #
import os.path
import datetime as dt
import glob
import re
import sys

# -- Third-part modules -- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import Functions

def robust_std(data):
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median))
    return mad * 1.4826  # scaling factor for normal distribution

dtint = 30  # days
gridres = 25000  # m

#directory = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) + '/RRDPp/RawData/BPR_ULS_BGEP'
directory = '/dmidata/projects/cmems2/C3S/RRDPp/RawData/BPR_ULS_BGEP'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/fig/'
save_path_data= os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/final/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/BGEP/final/ESACCIplus-SEAICE-RRDP2+-SID-BGEP_test.nc'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

files = glob.glob(save_path_data+'*')

lat=0
lon=0

count = 0
datalen = 0
for dir in ['2010-2011']: #2017-2021', '2021-2022', '2022-2023']: #sorted(os.listdir(directory)):
    print(dir)
    dir_data=os.path.join(directory, dir)
    if os.path.isdir(dir_data) and not dir_data.endswith('no_use'):        
        for ifile in os.listdir(dir_data)[:1]:
            file=os.path.join(dir_data, ifile)
            #defines variables in output file
            count+=1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
            dataOut.pp_flag = 0
            with open(file,'rb') as myFile:
                for num, line in enumerate(myFile, 1):
                    if num == 1:
                        line = line.strip()
                        line = line.split()
                        print(line)
                        if dir=='2022-2023': # different format
                            try:
                                dataOut.obsID = 'BGEP_Mooring' + 'B'
                                dataOut.lat = float(line[3].decode('utf8')) + float(line[4].decode('utf8'))/60 # convert to deicmal degrees
                                if line[8].decode('utf8') == 'W':
                                    dataOut.lon = - float(line[6].decode('utf8')) - float(line[7].decode('utf8'))/60
                                elif line[8].decode('utf8') == 'E':
                                    dataOut.lon = float(line[6].decode('utf8')) + float(line[7].decode('utf8'))/60
                                break
                            except:
                                dataOut.obsID = 'BGEP_Mooring' + line[3].decode('utf8')[:-1]
                                dataOut.lat = float(line[4].decode('utf8')) + float(line[5].decode('utf8'))/60 # convert to deicmal degrees
                                if line[9].decode('utf8') == 'W':
                                    dataOut.lon = - float(line[7].decode('utf8')) - float(line[8].decode('utf8'))/60
                                elif line[9].decode('utf8') == 'E':
                                    dataOut.lon = float(line[7].decode('utf8')) + float(line[8].decode('utf8'))/60
                                break  
                        else: # all other folders
                            dataOut.obsID = 'BGEP_Mooring' + line[3].decode('utf8')[:-1]
                            dataOut.lat = float(line[4].decode('utf8')) + float(line[5].decode('utf8'))/60 # convert to deicmal degrees
                            if line[9].decode('utf8') == 'W':
                                dataOut.lon = - float(line[7].decode('utf8')) - float(line[8].decode('utf8'))/60
                            elif line[9].decode('utf8') == 'E':
                                dataOut.lon = float(line[7].decode('utf8')) + float(line[8].decode('utf8'))/60
                            break
 

                        
            # Load data
            if dir == '2022-2023': # different format
                # Reads observation data from ASCII-file
                dff = pd.read_table(file, skiprows=2, names=['date', 'UTC', 'draft(m)'],sep="\s+", dtype = 'float32', converters = {'date': str})
                dates = dff['date'].to_numpy()
                time = dff['UTC'].to_numpy() #UTC - duration of each time used on each measurement
            else:
                # Reads observation data from ASCII-file
                dff = pd.read_table(file, skiprows=1, sep="\s+", dtype = 'float32', converters = {'%date': str})
                dates = dff['%date'].to_numpy()
                time = dff['time(UTC)'].to_numpy() #UTC - duration of each time used on each measurement
            
            SID = dff['draft(m)'].to_numpy() #m
            SID[SID>8] = np.nan
            SID[SID<0] = np.nan
            datalen += len(SID[np.isfinite(SID)])
            SID_Unc = np.array([0.10 for sid in SID]) # uncertainty of 10 cm approx.
            #print('Getting months')
            
            months = [int(date[4:6]) for date in dates]
            days = [int(date[6:]) for date in dates]
            mondiff=np.where(~(np.diff(months) == 0))[0] #index of when month changes
            index=np.insert(np.append(mondiff,len(months)-1),0,0) #add start and end indexes on month
            time_in = [dt.datetime(int(y[:4]),int(m),int(d)) for y,m,d in zip(dates, months, days)]
            
            # Avg_draft
            #dataOut.SID_final = np.array([np.nanmean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            #dataOut.SID_std = np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            
            fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(12, 12), constrained_layout=True)
            ax = ax.flatten()
            for i in range(len(index)-1)[:12]:
                median_val = np.nanmedian(SID[index[i]:index[i+1]])
                mean_val   = np.nanmean(SID[index[i]:index[i+1]])
                ax[i].hist(SID[index[i]:index[i+1]], bins=50, label=f'Median: {median_val:.2f}, Mean: {mean_val:.2f}')
                ax[i].set_xlabel('SIT [m]')
                ax[i].set_ylabel('count')
                ax[i].legend()
            fig.suptitle('BGEP: Example of monthly distribution of data')
            plt.savefig(f'sid_test_BGEP.png')   
            plt.close()
            
            dataOut.SID_final=np.array([np.nanmedian(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.SID_std = np.array([robust_std(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])

            dataOut.SID_ln = np.array([len(SID[np.isfinite(SID)][index[i]:index[i+1]]) for i in range(len(index)-1)])
            
            ######### QUALITY FLAGS ############
            dataOut.QFT = [] # temporal
            dataOut.QFS = [] # spatial
            dataOut.QFG = [] # global threshold

            for i in range(len(index)-1):
                # find number of days
                unique = len(np.unique(days[index[i]:index[i+1]]))
                if np.any(SID[index[i]:index[i+1]]>8):
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
                df = dataOut.Create_NC_file(ofile, primary='SID', datasource='Beaufort Gyre Exploration Project, Mooring Data: https://www2.whoi.edu/site/beaufortgyre/data/mooring-data/', key_variables='Sea Ice Draft')
                

print(datalen)
# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SID')