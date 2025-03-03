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
# hemisphere = 'N'

# # Antarctic
dtint = 30  # days
gridres = 25000  # m
hemisphere = 'NH'

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


                        # Reads date format into datetime format
                        date_string = date.astype(str)

                        t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                                      for s in date_string])

                        # SD uncertainty based on Nicholaus et al. 2017
                        SD_unc = np.array([0.02 for s in SD])
                        dataOut.unc_flag = 3
                        dataOut.pp_flag = 0
                        
                        write_info_table(dataOut, t, SD_init, SIT_init, num_sensors)
                        
                        # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
                        G = EASEgrid.Gridded()
                        if hemisphere=='SH':
                            G.SetHemisphere('S')
                        else:
                            G.SetHemisphere('N')
                        G.CreateGrids(gridres)
                        (index_i, index_j) = G.LatLonToIdx(latitude[index], longitude[index])
            
                        # Takes the time for each grid cell into account and calculate averages
                        avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, avgTair, avgTsurf = G.GridData(
                            dtint, latitude[index], longitude[index], t[index], SD=SD[index], SD_unc=SD_unc[index], VAR1=Tair[index])

            
                        if len(time) > 0:
                            Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS=hemisphere)
                            Functions.scatter(dataOut.obsID, t[index], SD[index], time, avgSD, 'SD [m]', saveplot)

                        if hemisphere=='NH':
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
                        dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
                        dataOut.SD_final = avgSD
                        dataOut.SD_std = stdSD
                        dataOut.SD_ln = lnSD
                        dataOut.SD_unc = uncSD
                        dataOut.SIT_final = avgSIT
                        dataOut.SIT_std = stdSIT
                        dataOut.SIT_ln = lnSIT
                        dataOut.SIT_unc = uncSIT
                        dataOut.air_temp_final = avgTair
                        dataOut.unc_flag = [dataOut.unc_flag]*len(lat)
                        dataOut.pp_flag = [dataOut.pp_flag]*len(lat)
                        
                        dataOut.obsID = [dataOut.obsID]*len(dataOut.SD_final)
                        
                        # fill empty arrays with NaN values
                        dataOut.Check_Output()
                            
                        if count>1:
                            subset = dataOut.Create_NC_file(ofile, primary='SD')
                            df = Functions.Append_to_NC(df, subset)
                        else:
                            df = dataOut.Create_NC_file(ofile, primary='SD', datasource='Autonomous measurements of sea ice properties: https://data.meereisportal.de/relaunch/buoy.php?lang=en&active-tab1=method&active-tab2=buoy&ice-type=buoy&region=all&buoytype=SB&buoystate=all&expedition=all&showMaps=y&dateRepeat=n&timeline=', key_variables='Snow depth')

# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SD')
