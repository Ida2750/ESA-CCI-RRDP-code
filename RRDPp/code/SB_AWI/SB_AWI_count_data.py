# -*- coding: utf-8 -*-
"""
Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Alfred Wegener Institute Snow Depth Buoys.
includes Warren snow depths and densities
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '1'
__date__ = '2023-06-20'

# -- Built-in modules -- #
import os.path
import pdb
import sys
import datetime as dt

# -- Third-part modules -- #
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions
from PDF_read_initial_SD import pdf_read_initial


def write_info_table(self, date, SD_init, SIT_init, num_sh):
    """
    Parameters
    ----------
    information: time span of mass balance measurements, 
    of temperature measurements, 
    of meteorological data
    
    pre-processing flag
    uncertainty flag
    
    variables in file
    -------
    None.
    """

    date_start = str(date[0].date())
    print(date_start)
    date_end = str(date[-1].date())
    
    num_sh = len(num_sh)
    ofile = 'SB_'+hemisphere+'_info_table.dat'
    output = open(ofile, 'a')
    
    print('{:^s}  & {:^s} & {:^s} & {:^s} & {:^6.2f} & {:^6.2f} & {:^1.0f} & {:^1.0f} & {:^1.0f} \\'.format(
        self.obsID.replace('_', '-'), hemisphere + 'H', date_start, date_end, SD_init, SIT_init, num_sh, np.nanmean(self.pp_flag), np.nanmean(self.unc_flag)), file=output)
    output.close()

    return None

#%% Main
# Arctic
# dtint=30 #days
# gridres=25000 #m
# hemisphere = 'NH'

# # Antarctic
dtint = 30  # days
gridres = 50000  # m
hemisphere = 'SH'

if hemisphere == 'SH':
    directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/SDB_AWI/Antarctic'
    save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/SB_AWI/final/'
    if not os.path.exists(save_path_data):os.makedirs(save_path_data)
    saveplot = os.path.join(os.path.dirname(
        os.path.dirname(os.getcwd())), 'Final/Antarctic/SB_AWI/fig/')
else:
    directory =  os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/SDB_AWI/Arctic'
    save_path_data =  os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/SB_AWI/final/'
    if not os.path.exists(save_path_data):os.makedirs(save_path_data)
    saveplot = os.path.join(os.path.dirname(
        os.path.dirname(os.getcwd())), 'Final/SB_AWI/fig/')
ofile = f'ESACCIplus-SEAICE-RRDP2+-SD-SB-AWI-{hemisphere}.nc'

# count data
total_obs_valid_SD = 0

for dir in os.listdir(directory):
    count = 0
    print(os.listdir(directory))
    if dir.startswith('bouy'):
        dir_dir = os.path.join(directory, dir)
        for dir2 in os.listdir(dir_dir):
            dir_data = os.path.join(dir_dir, dir2)
            print(dir_data)
            if os.path.isdir(dir_data):  # only enter if is directory
                files = sorted(os.listdir(dir_data), reverse=True)
                for ifile in files:
                    # reverse order to get initial info first
                    # break
                    if 'deployment' in ifile:
                        [sd1_i, sd2_i, sd3_i, sd4_i, SD_init, SIT_init] = pdf_read_initial(
                            os.path.join(dir_data, ifile))
                        snows = [sd1_i, sd2_i, sd3_i, sd4_i]
                        num_sensors = [i for i, s in zip(range(len(snows)), snows) if ~np.isnan(s)]
                        
                    elif ifile.endswith('proc.csv'):
                        file = os.path.join(dir_data, ifile)
                        ofile = os.path.join(save_path_data, ofile)

                        #defines variables in output file
                        count+=1
                        dataOut = Functions.Final_Data(Type='SD', count_head=count)

                        # read data from one ASCII-file, content of 1 month
                        data2 = open(file, 'r')
                        sb = np.genfromtxt(data2, skip_header=1, dtype=None, usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
                                           delimiter=',', names=('date', 'lat', 'lon', 'sd1', 'sd2', 'sd3', 'sd4', 'patm', 'tair', 'tbody'), encoding=None)

                        date = sb['date']
                        latitude = sb['lat']
                        longitude = sb['lon']
                        longitude = np.array(
                            [l if l <= 180 else l-360 for l in longitude])
                        # remember to add initial SD
                        snow1 = (sb['sd1'] + sd1_i)
                        snow2 = (sb['sd2'] + sd2_i)
                        snow3 = (sb['sd3'] + sd3_i)
                        snow4 = (sb['sd4'] + sd4_i)
                        if any(snow1 < 0) or any(snow2 < 0) or any(snow3 < 0) or any(snow4 < 0):
                            break
                            print('SD lower than 0')
                        Tair = (sb['tair'])

                        #define obsID
                        dataOut.obsID = 'SB_AWI_'+ifile[:6]

                        # Calculate average snow depth of snow1-snow4
                        snow1234 = np.row_stack((snow1, snow2, snow3, snow4))
                        SD = np.nanmean(snow1234, axis=0)
                        index = ~np.isnan(SD)

                        total_obs_valid_SD += len(SD[np.isfinite(SD)])

print(total_obs_valid_SD)