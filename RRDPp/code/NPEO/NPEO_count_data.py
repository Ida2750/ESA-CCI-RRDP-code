# -*- coding: utf-8 -*-
"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured at the North Pole Enviromental Observatory
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

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import Functions


def Get_lat_lon(file):
    """
    Get lat and lon information from file header

    Parameters
    ----------
    file : TYPE
        inputfilename with absolute path

    Returns
    -------
    lat : float
        Mooring latitude
    lon : float
        Mooring longitude
    header : int
        number of headerlines

    """
    lat = 0
    lon = 0
    #to find end of header
    lookup = b'#'
    #to find location of lat and lon
    latlon = b'Position'
    count = 0
    with open(file, 'rb') as myFile:
        for num, line in enumerate(myFile, 1):
            #print(line.strip())
            count += 1
            #either if latlon in line or if the lat/lon is split over two lines
            if latlon in line or (lat != 0 and lon == 0):
                if 'N' in line.decode('latin-1') and 'E' in line.decode('latin-1'):
                    #print(line.decode('latin-1'))
                    #print('-----------')
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = float(latlon_num[3])+lonmin/60
                elif 'N' in line.decode('latin-1') and 'W' in line.decode('latin-1'):
                    ##if coordinates given to the west then the longitude must be -
                    #print(line.decode('latin-1'))
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = -(float(latlon_num[3])+lonmin/60)

                # if lat/lon is split over two lines
                elif 'North' in line.decode('latin-1'):
                    lat_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(lat_num[1]+'.'+lat_num[2])
                    # convert to decimal degrees
                    lat = float(lat_num[0])+latmin/60

                elif 'East' in line.decode('latin-1'):
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = float(lon_num[0])+lonmin/60
                elif 'West' in line.decode('latin-1'):
                    print(line.decode('latin-1'))
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = -float(lon_num[0])+lonmin/60
            #break loop when lookup is found
            if lookup in line:
                header = count  # note number of headerlines
                break

    return lat, lon, header


dtint = 30  # days
gridres = 25000  # m

directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.getcwd())))) + '/RRDPp/RawData/NPEO/data'

# saving locations
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPEO/final/'
ofile = save_path_data + 'ESACCIplus-SEAICE-RRDP2-SID-NPEO.nc'
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPEO/fig/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

total_obs_valid_SID = 0
count = 0
for dir in os.listdir(directory):
    ## the 2001 dir is different from the rest
    if '2001' in dir:
        dir_data = os.path.join(directory, dir+'/ULS')
    else:
        dir_data = os.path.join(directory, dir+'/'+dir+'/ULS')
    #print(dir_data)
    for ifile in os.listdir(dir_data):
        print(ifile)
        #only the files containing min are of interest to us
        if 'min' in ifile:
            file = os.path.join(dir_data, ifile)
            lat, lon, header = Get_lat_lon(file)

            #defines variables in output file
            count += 1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
            dataOut.pp_flag = 0
            dataOut.unc_flag = 2

            # # Reads observation data from ASCII-file
            data = np.genfromtxt(file, dtype=None, skip_header=header, names=[
                                 'YEAR', 'MO', 'DAY', 'HR', 'MIN', 'SEC', 'TIME', 'DEPTH', 'DRAFT'], encoding='latin-1')

            years = data['YEAR']
            months = data['MO']
            days = data['DAY']
            hr = data['HR']
            min = data['MIN']
            sec = data['SEC']
            dates = [dt.datetime(years[i], months[i], days[i],
                                 hr[i], min[i], sec[i]) for i in range(len(years))]
            dataOut.obsID = 'NPEO_' + ifile.split('_')[2]
            SID = data['DRAFT'].astype(float)

            # uncertainty of 5 cm approx.
            SID_Unc = np.array([0.05 for sid in SID])
            # set non existing data to nan
            SID[SID == -999.000] = np.nan

            total_obs_valid_SID += len(SID[np.isfinite(SID)])

print(total_obs_valid_SID)