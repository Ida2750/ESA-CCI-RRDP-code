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
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from Get_Dates_From_Header import Get_date as gd
from Get_New_AEM_Data import Get_New_AEM_data as gnad

from Warren import SnowDepth, SWE
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean

# Define locations of input and output files
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/final/'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/AEM-AWI/final/ESACCIplus-SEAICE-RRDP2+-SIT-AEM-AWI-V3.dat'
directory = os.path.dirname(os.path.dirname(os.getcwd())) +'/RawData/AEM_AWI'

count = 0  # used to locate header

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

for dir in directories:
    # dir = '2013'
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all SIT files in AEM_AWI
            if (file.endswith('.tab') or file.endswith('.dat') or file.endswith('.nc')) and 'free' not in ifile:
                print(file)
                
                #defines variables in output file
                count+=1
                dataOut = Functions.Final_Data(Type='SIT', count_head=count)
                dataOut.pp_flag = 0
                
                if file.endswith('.nc') or dir == '2019': # accessing newer files (has different format)
                    # print('Entering script')
                    dataOut.obsID, latitude, longitude, SIT, t = gnad(file)
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
                    SIT[SIT < 0] = np.nan
    
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
                        dataOut.pp_flag = 1
                        # print(t)
                    

                    # Define obsID
                    flight = ifile.replace('sea-ice_thickness.tab', '')
                    flight = flight.replace('ice_thickness.tab', '')
                    flight = flight.replace('HEM_', '')
                    if 'GreenICE' in flight:  # avoid very long obsID names
                        flight = flight.replace('GreenICE', 'GICE')
                    elif 'PoleAirship' in flight:  # avoid very long obsID names
                        flight = flight.replace('PoleAirship', 'PAS')
                    elif 'Hendricks' in flight:
                        flight = 'STJ06_01'  # Storfjorden flight name
                    dataOut.obsID = 'AEM_AWI_' + flight[:-1]       
                
                ## Define fixed uncertainty
                SIT_unc = np.ones(len(latitude)) * 0.10
                dataOut.unc_flag = 2
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i, index_j) = G.LatLonToIdx(latitude, longitude)
    
                # Takes the time for each grid cell into account and calculate averages
                avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc = G.GridData(
                    dtint, latitude, longitude, t, SIT=SIT, SIT_unc=SIT_unc)
        
                if len(time) > 0:
                    Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS='NH')
                    Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]', saveplot)


                # Correlates SB-AWI buoy data with Warren snow depth and snow density
                for ll in range(np.size(avgSD,0)):
                    (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
                    dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
                    (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
                    w_density=int((wswe/w_SD)*1000)
                    dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)
    
                #Change names to correct format names
                dataOut.lat_final = lat
                dataOut.lon_final = lon
                for ll in range(np.size(time,0)):
                    dataOut.date_final = np.append(dataOut.date_final,dt.datetime.strftime(time[ll],"%Y-%m-%dT%H:%M:%S"))
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
                    
                # print data to output file
                dataOut.Print_to_output(ofile, primary='SIT')
            
Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH')

                    