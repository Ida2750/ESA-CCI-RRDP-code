#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 25 km gridded means of input files containing measurements of SIT from 
airborne campaigns using the EM-bird from the Alfred Wegener institut 
includes Warren snow depths and densities
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2022-01-25'

# -- Built-in modules -- #
import os.path
import datetime as dt

import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from Warren import SnowDepth, SWE
import EASEgrid as EASEgrid
from Get_Dates_From_Header import Get_date as gd
from Get_New_AEM_Data import Get_New_AEM_data as gnad

# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean

# Define locations of input and output files
save_path_data = '/media/s174020/My Passport2/ESACCIplus/RRDPp/FINAL/OIB/final/'
ofile = '/media/s174020/My Passport2/ESACCIplus/RRDPp/FINAL/AEM_AWI/final/ESACCIplus-SEAICE-RRDP2+-SIT-AEM-AWI.dat'
directory = '/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/AEM_AWI'

count_head = 0  # used to locate header

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

for dir in directories:
    # dir = '2019'
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all SIT files in AEM_AWI
            if (file.endswith('.tab') or file.endswith('.dat') or file.endswith('.nc')) and 'free' not in ifile:
                print(file)
                count_head += 1
                
                #defines variables in output file
                flag = 0
                obsID = [];
                date_final = []; 
                lat_final = []; # decimal degrees
                lon_final = []; # decimal degrees
                SD_final = [];  # snow depth [m]
                SD_std = [];  # standard deviation of snow depth
                SD_ln = 0   # number of SD measurements in each gridcell
                SD_unc = []  # SD uncertainty
                SIT_final = [];  # Sea Ice Thickness [m]
                SIT_std = [];  # sea Ice thickness standard deviation
                freeboard_final = [];  # freeboard [m]
                freeboard_std = [] # freeboard standard deviation [m]
                Frb_ln = 0  # number of FRB measurements in each gridcell
                FRB_unc = []  # FRB uncertainty
                sur_temp_final = [];  # surface temperature [celcius]
                air_temp_final = [];  # air temperature [celcius]
                w_SD_final = [];  # snow depth from climatology (Warren)
                w_density_final = [];  # snow density from climatology (Warren)
                pp_flag = 'no'  # pre-processing flag
                unc_flag = 'yes'  # Uncertainty flag
                
                if file.endswith('.nc') or dir == '2019': # accessing newer files (has different format)
                    # print('Entering script')
                    obsID, latitude, longitude, SIT, t = gnad(file)
                    # print('Finished script')
                else:  # acessing older files (prior to 2012)
                    if dir == 'Feb_Mar_2004':  # Interpolated dates
                        names = ['DateTime', 'distance', 'Latitude', 'Longitude','Sample ID','EsEs','Height [m]','time_diff']
                        header = -1
                        pp_flag = 'yes'
                    else:  # find number of headerlines and names of variables in file
                        lookup = b'*/'  # notation for when header is done
                        count = 0
                        header = -999
                        with open(file,'rb') as myFile:
                            for num, line in enumerate(myFile, 1):
                                # print(line.strip())
                                count+=1
                                # break loop when lookup is found
                                if lookup in line:
                                    header = count  # note number of headerlines
                                if count == header+1:
                                    headerline = line.decode('utf8').replace("[m]","")  #remove [m]
                                    names = headerline.split("\t")
                                    names = [name.strip() for name in names]
                                    break
    
                    # Read data
                    if dir == 'Mar_Apr_2003':  # Problems with gaps in height data
                        # names = names[]
                        data = np.genfromtxt(file, dtype=None, names=names, skip_header=header+1, usecols=np.arange(0,6)) 
                    else:
                        data = np.genfromtxt(file, dtype=None, names=names, skip_header=header+1) 
    
                    latitude = data['Latitude']
                    longitude = data['Longitude']
                    SIT = data['EsEs']
                    # set nonexistin data to nan
                    SIT[SIT==-99999.0] = np.nan
    
                    try: # get date
                        date = data['DateTime']
                        # Convert date format into python date-time format
                        dates = date.astype(str)
                        if dir =='Apr_2009_SIZONet':
                            t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M")
                                         for s in dates])
                        elif dir == 'Mar_Apr_2003' or dir == 'May_2004':
                            t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%f")
                                         for s in dates])
                        else:
                            t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                                         for s in dates])
                    except ValueError:  # no date provided in file
                        t = gd(file, ifile)  # Calls Get_Date file
                        flag = 1
                        pp_flag = 'yes'
                        # print(t)
                    
                    # Define obsID
                    flight = ifile.replace('sea-ice_thickness.tab', '')
                    flight = ifile.replace('ice_thickness.tab', '')
                    flight = flight.replace('HEM_', '')
                    if 'GreenICE' in flight:  # avoid very long obsID names
                        flight.replace('GreenICE', 'GICE')
                    elif 'PoleAirship' in flight:  # avoid very long obsID names
                        flight.replace('PoleAirship', 'PAS')
                    elif 'Hendricks' in flight:
                        flight = 'STJ06_01'  # Storfjorden flight name
                    date_start = str(t[0].date()).replace('-', '')
                    obsID = 'AEM_AWI_' + flight + date_start        
                
                # Create EASEgrid and return grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # Take the time for each grid cell into account and calculate 25 km averages
                SIT_unc = np.ones(len(latitude)) * 0.10
                SIT_unc[np.isnan(SIT)] = np.nan
                (lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT) = G.GridData(dtint, latitude,
                                                                             longitude, SIT, SIT_unc, t)
               
                # Sort data according to date
                index_array = np.zeros(np.shape(time)) 
                time_sorted  =  []
                avgSIT_sorted = []
                stdSIT_sorted = []
                lnSIT_sorted  = []
                uncSIT_sorted = []
                lat_sorted   =  []
                lon_sorted   =  []
                
                for kk in range(len(time)):
                    index_array = np.argsort(time)
                    time_sorted = np.append(time_sorted,time[index_array[kk]])
                    avgSIT_sorted = np.append(avgSIT_sorted,avgSIT[index_array[kk]])
                    stdSIT_sorted = np.append(stdSIT_sorted,stdSIT[index_array[kk]])
                    lnSIT_sorted = np.append(lnSIT_sorted,lnSIT[index_array[kk]])
                    uncSIT_sorted = np.append(uncSIT_sorted,uncSIT[index_array[kk]])
                    lat_sorted = np.append(lat_sorted,lat[index_array[kk]])
                    lon_sorted = np.append(lon_sorted,lon[index_array[kk]])
                
                # Correlate AEM_AWI data with Warren snow depth and snow density
                w_SD_final      = []
                w_density_final = []
                for ll in range(np.size(avgSIT_sorted,0)):
                    (w_SD,w_SD_epsilon) = SnowDepth(lat_sorted[ll],lon_sorted[ll],time_sorted[ll].month)
                    w_SD_final = np.append(w_SD_final,w_SD)
                    (wswe,wswe_epsilon) = SWE(lat_sorted[ll],lon_sorted[ll],time_sorted[ll].month)
                    w_density=int((wswe/w_SD)*1000)
                    w_density_final = np.append(w_density_final,w_density)
                
                # Change names to correct format names
                lat_final=lat_sorted
                lon_final=lon_sorted
                
                for ll in range(np.size(avgSIT_sorted,0)):
                    if flag == 1:
                         date_final=np.append(date_final,dt.datetime.strftime(t[0],"%Y-%m-%dT%H:%M:%S"))
                    else:
                        date_final=np.append(date_final,dt.datetime.strftime(time_sorted[ll],"%Y-%m-%dT%H:%M:%S"))        
                SIT_final = avgSIT_sorted
                SIT_std = stdSIT_sorted
                SIT_ln = lnSIT_sorted
                SIT_unc = uncSIT_sorted
 
                # Add nan values to nonexisting data
                if len(SD_final)==0:
                    SD_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_std)==0:
                    SD_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_unc)==0:
                    SD_unc = np.zeros(np.shape(lat_final)) * np.nan
                if SD_ln == 0:
                    SD_ln = np.zeros(np.shape(lat_final))
                if len(SIT_final)==0:
                    SIT_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SIT_std)==0:
                    SIT_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(freeboard_final)==0:
                    freeboard_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(freeboard_std)==0:
                    freeboard_std = np.zeros(np.shape(lat_final)) * np.nan
                if Frb_ln == 0:
                    Frb_ln = np.zeros(np.shape(lat_final))
                if len(FRB_unc) == 0:
                    FRB_unc = np.zeros(np.shape(lat_final)) * np.nan
                if len(sur_temp_final)==0:
                    sur_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(air_temp_final)==0:
                    air_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                
                # print data to output file
                if count_head == 1: # make empty output file with header
                    output = open(ofile,'w')
                    print('{:^30s} {:^20s} {:^8s} {:^8s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^6s} {:^7s} {:^7s} {:^6s} {:^4s} {:^7s} {:^7s}'.format('obsID', 'date', 'lat', 'lon', 'SD', 'SDstd', 'SDln', 'SDunc', 'SIT', 'SITstd', 'SITln','SITunc', 'Frb', 'FRBstd', 'FRBln', 'FRBunc', 'Tsur', 'Tair', 'wsd', 'w\u03C1', 'pp-flag', 'unc-flag'),file=output)
                output = open(ofile,'a')
                for ll in range(len(lat_final)):
                    if SIT_ln[ll] != 0: # print output whenever SIT data exist
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], freeboard_final[ll],freeboard_std[ll], Frb_ln[ll], FRB_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag),file=output)
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], freeboard_final[ll],freeboard_std[ll], Frb_ln[ll], FRB_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag))
                
                
                output.close()
                    
