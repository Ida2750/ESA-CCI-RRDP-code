# -*- coding: utf-8 -*-
"""

Creates 25 km distance means of input files measured by airborne altimetry from the Operation IceBridge by NASA
includes Warren snow depths and densities

"""
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2024-08-02'

# -- Built-in modules -- #
import os.path
import glob
import datetime as dt
import re
import sys

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions

dtint = 30  # days per mean
hemisphere = 'SH'
if hemisphere == 'NH':
    # determine gridsize and temporal limits
    gridres = 25000  # grid resolution
    # Reads input data
    save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/OIB/final/'
    saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/OIB/fig/'
    ofile = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/OIB/final/ESACCIplus-SEAICE-RRDP2+-SIT-OIB-{hemisphere}.nc'
    directory1 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/OIB/IDCS4'
    directory2 = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/OIB/QuickLook'
elif hemisphere == 'SH':
    # determine gridsize and temporal limits
    gridres = 50000  # grid resolution
    # Reads input data
    save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/OIB/final/'
    saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/OIB/fig/'
    ofile = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/Antarctic/OIB/final/ESACCIplus-SEAICE-RRDP2+-SIT-OIB-{hemisphere}.nc'
    directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/Antarctic/OIB'

if not os.path.exists(save_path_data):os.makedirs(save_path_data)
count = 0  # used to locate header

# Access data
if hemisphere =='NH':
    directories = glob.glob(f'{directory1}/*') + glob.glob(f'{directory2}/*')
else:
    directories = glob.glob(f'{directory}/*')
    
directories.sort()  # sort to ensure starting with the earliest date

total_obs_valid_SIT = 0
total_obs_valid_SD = 0
total_obs_valid_FRB = 0

for dir_data in directories:
    if os.path.isdir(dir_data):        
        for ifile in os.listdir(dir_data):
            file = os.path.join(dir_data, ifile)
            if file.endswith('txt'):
                
                # Reads observation data from ASCII-file
                data = np.genfromtxt(file,dtype=None,names=True,delimiter=',',usecols=np.arange(0,50), encoding=None) #usecols: specifies which collumns to read
                
                #defines variables in output file
                count+=1
                dataOut = Functions.Final_Data(Type='FRB', count_head=count)
                dataOut.pp_flag = 0
                dataOut.unc_flag = 0
                
                #define obsID
                ID_numbers = re.findall(r'\d+', ifile)  # extract numbers from title
                ID_date = [int(num) for num in ID_numbers if int(num) > 100]  # Get date
                
                # Determine if data is from QL or IDCS4
                if 'QuickLook' in dir_data: 
                    dataOut.obsID = 'OIB_QL_' + str(ID_date[0])
                elif 'IDCS4' in dir_data:    
                    dataOut.obsID = 'OIB_IDCS4_'  + str(ID_date[0])
                else:
                    # Specify that data is from IDCS4   
                    dataOut.obsID = 'OIB_IDCS4_'  + str(ID_date[0])
                
                # Read data
                latitude = data['lat']
                longitude = data['lon']
                date = data['date']
                time = data['elapsed']
                SIT = data['thickness']
                SIT_unc = data['thickness_unc']
                meanFrb = data['mean_fb']
                frbUnc = data['fb_unc']
                OIBSnow = data['snow_depth']
                OIBSnowUnc = data['snow_depth_unc']
                roughness = data['surface_roughness']
                
                #convert lon to be between -180 and 180
                for ii in range(len(longitude)):
                    if longitude[ii] > 180:
                        longitude[ii] = longitude[ii]-360
                    else:
                        longitude[ii] = longitude[ii]
                
                # Set non-existing data to nan
                SIT[SIT==-99999.0] = np.nan
                SIT_unc[SIT_unc==-99999.0] = np.nan
                meanFrb[meanFrb==-99999.0] = np.nan
                frbUnc[frbUnc==-99999.0] = np.nan
                OIBSnow[OIBSnow==-99999.0] = np.nan
                OIBSnowUnc[OIBSnowUnc==-99999.0] = np.nan
                roughness[roughness==-99999.0] = np.nan
                
                total_obs_valid_SIT += len(SIT[np.isfinite(SIT)])
                total_obs_valid_SD += len(OIBSnow[np.isfinite(OIBSnow)])
                total_obs_valid_FRB += len(meanFrb[np.isfinite(meanFrb)])
                

print(total_obs_valid_SIT)
print(total_obs_valid_SD)
print(total_obs_valid_FRB)