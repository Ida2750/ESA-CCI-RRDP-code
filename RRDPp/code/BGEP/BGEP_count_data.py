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

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import Functions

dtint = 30  # days
gridres = 25000  # m

directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/BPR_ULS_BGEP'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/fig/'
save_path_data= os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/BGEP/final/'
ofile = os.path.dirname(os.path.dirname(os.getcwd())) +'/FINAL/BGEP/final/ESACCIplus-SEAICE-RRDP2+-SID-BGEP.nc'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

files = glob.glob(save_path_data+'*')
# for f in files:
#     os.remove(f)

lat=0
lon=0

total_obs_valid_SID = 0
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
            dff = pd.read_table(file, skiprows=1, sep="\s+", dtype = 'float32', converters = {'%date': str})
                        
            # Load data
            dates = dff['%date'].to_numpy()
            time = dff['time(UTC)'].to_numpy() #UTC - duration of each time used on each measurement
            SID = dff['draft(m)'].to_numpy() #m


            total_obs_valid_SID += len(SID[np.isfinite(SID)])

print(total_obs_valid_SID)