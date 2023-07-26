#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 16:38:42 2022

@author: s174020
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 25 km gridded means of input files containing measurements of FRB from 
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
from Warren import SnowDepth, SWE
from Get_Dates_From_Header import Get_date as gd
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean

# Define locations of input and output files
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/fig/'
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/final/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/final/ESACCIplus-SEAICE-RRDP2+-FRB-AEM-AWI-V3.dat'
directory = os.path.dirname(os.path.dirname(os.getcwd())) + '/RawData/AEM_AWI'

count = 0  # used to locate header

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

for dir in directories:
    # dir = 'Feb_Mar_2004'
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all SIT files in AEM_AWI
            if file.endswith('.tab')  and 'free' in ifile:
                print(file)
                
                #defines variables in output file
                count+=1
                dataOut = Functions.Final_Data(Type='FRB', count_head=count)
                dataOut.pp_flag = 0

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
                FRB = data['Freeboard']
                # set nonexistin data to nan
                FRB[FRB==-99999.0] = np.nan

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
                    dataOut.pp_flag = 1
                
                # Define obsID
                flight = ifile.replace('_freeboard.tab', '')
                flight = flight.replace('HEM_', '')
                if 'GreenICE' in flight:  # avoid very long obsID names
                    flight = flight.replace('GreenICE', 'GICE')
                flight = flight.split('_')[0]
                dataOut.obsID = 'AEM_AWI_' + flight      
                
                # Create EASEgrid and return grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # Take the time for each grid cell into account and calculate 25 km averages
                FRB_unc = np.ones(len(latitude)) * 0.1  # uncertainty of approx. 10 cm 
                dataOut.unc_flag = 2
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i, index_j) = G.LatLonToIdx(latitude, longitude)
    
                # Takes the time for each grid cell into account and calculate averages
                avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, uncFRB = G.GridData(
                    dtint, latitude, longitude, t, FRB=FRB, FRB_unc=FRB_unc)
    
                if len(time) > 0:
                    Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS='NH')
                    Functions.scatter(dataOut.obsID, t, FRB, time, avgFRB, 'FRB [m]', saveplot)


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
                dataOut.FRB_final = avgFRB
                dataOut.FRB_std = stdFRB
                dataOut.FRB_ln = lnFRB
                dataOut.FRB_unc = uncFRB

                
                # fill empty arrays with NaN values
                dataOut.Check_Output()
                    
                # print data to output file
                dataOut.Print_to_output(ofile, primary='FRB')
        
Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH')

                    