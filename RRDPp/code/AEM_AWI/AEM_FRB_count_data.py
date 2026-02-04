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
__contact__ = ['ilo@dmi.dk']
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
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/AEM-AWI/fig/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/AEM-AWI/final/ESACCIplus-SEAICE-RRDP2+-FRB-AEM-AWI.nc'
directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) +'/RRDPp/RawData/AEM_AWI'


count = 0  # used to locate header

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

# count data
total_obs_valid_FRB = 0
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
                counthead = 0
                header = -999
                with open(file,'rb') as myFile:
                    for num, line in enumerate(myFile, 1):
                        # print(line.strip())
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
                FRB = data['Freeboard']
                # set nonexistin data to nan
                FRB[FRB==-99999.0] = np.nan
                
total_obs_valid_FRB += len(FRB[np.isfinite(FRB)])

print(total_obs_valid_FRB)
                    
