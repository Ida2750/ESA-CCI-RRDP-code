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
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2024-08-04'

# -- Built-in modules -- #
import os.path
import datetime as dt
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Get_Dates_From_Header import Get_date as gd
from Get_New_AEM_Data import Get_New_AEM_data as gnad
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions

#%%
# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean

# Define locations of input and output files
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/AEM-AWI/final/ESACCIplus-SEAICE-RRDP2+-SIT-AEM-AWI-V3.nc'
directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) +'/RRDPp/RawData/AEM_AWI'

count = 0  # used to locate first file

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

# count data
total_obs_valid_SIT = 0
total_obs_valid_SIT_aircraft = 0
total_obs_valid_SIT_heli = 0

for dir in directories:
    # dir = '2013'
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            aircraft=0
            file = os.path.join(dir_data, ifile)
            # Access all SIT files in AEM_AWI
            if (file.endswith('.tab') or file.endswith('.dat') or file.endswith('.nc')) and 'free' not in ifile:
                #print(file)
                
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
                    else:  # find number of headerlines and names of variables in file
                        lookup = b'*/'  # notation for when header is done
                        counthead = 0
                        header = -999
                        with open(file,'rb') as myFile:
                            for num, line in enumerate(myFile, 1):
                                if 'Aircraft' in line.decode('utf8'):
                                    print(line)
                                    aircraft = 1
                                #print(line)
                                counthead+=1
                                # break loop when lookup is found
                                if lookup in line:
                                    header = counthead  # note number of headerlines
                                if counthead == header+1:
                                    headerline = line.decode('utf8').replace("[m]","")  #remove [m]
                                    names = headerline.split("\t")
                                    names = [name.strip() for name in names]
                                    break
    
                    # Read data
                    if dir == 'Mar_Apr_2003':  # Problems with gaps in height data
                        # names = names[]
                        data = np.genfromtxt(file, dtype=None, names=names, skip_header=header+1, usecols=np.arange(0,6), encoding=None) 
                    else:
                        data = np.genfromtxt(file, dtype=None, names=names, skip_header=header+1, encoding=None) 
    
                    latitude = data['Latitude']
                    longitude = data['Longitude']
                    SIT = data['EsEs']
                    # set nonexistin data to nan
                    SIT[SIT==-99999.0] = np.nan
                    SIT[SIT < 0] = np.nan
    
                    total_obs_valid_SIT += len(SIT[np.isfinite(SIT)])
                    if aircraft==1:
                        total_obs_valid_SIT_aircraft += len(SIT[np.isfinite(SIT)])
                    else:
                        total_obs_valid_SIT_heli += len(SIT[np.isfinite(SIT)])

print(total_obs_valid_SIT_aircraft)
print(total_obs_valid_SIT_heli)

                    
