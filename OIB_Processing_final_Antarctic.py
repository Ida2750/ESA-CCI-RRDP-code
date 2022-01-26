#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creates 25 km distance means of input files measured by airborne altimetry from the Operation IceBridge by NASA
includes Warren snow depths and densities

Pre-processing flag is used to indicate whether out of the ordinary
pre-processing has been done. E.g. interpolation of location and/or dates.

Uncertainty flag is used to indicate whether the uncertainties are linked to
approximations and/or assumptions. E.g. if uncertainties are only provided for level ice
"""
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2022-01-24'

# -- Built-in modules -- #
import os.path
import datetime as dt
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from Warren import SnowDepth, SWE
import EASEgrid_OIB as EASEgrid

# determine gridsize and temporal limits
gridres = 50000  # grid resolution
dtint = 30  # days per mean

# Reads input data
save_path_data = '/media/s174020/My Passport2/ESACCIplus/RRDPp/FINAL/Antarctic/OIB/final/'
ofile = '/media/s174020/My Passport2/ESACCIplus/RRDPp/FINAL/Antarctic/OIB/final/ESACCIplus-SEAICE-RRDP2+-SIT-OIB-IDCS4-ANTARCTIC.dat'
directory = '/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/Antarctic/OIB'

count = 0  # used to locate header

# Access data
directories = os.listdir(directory)
directories.sort()  # sort to ensure starting with the earliest date

for dir in directories:
    print(dir)
    dir_data = os.path.join(directory, dir)
    if os.path.isdir(dir_data):        
        for ifile in os.listdir(dir_data):
            file = os.path.join(dir_data, ifile)
            if file.endswith('txt'):
                count += 1
                
                # Reads observation data from ASCII-file
                data = np.genfromtxt(file,dtype=None,names=True,delimiter=',',usecols=np.arange(0,50)) #usecols: specifies which collumns to read
                
                #defines variables in output file
                obsID = [];
                date_final = []; 
                lat_final = []; # decimal degrees
                lon_final = []; # decimal degrees
                SD_final = [];  # snow depth [m]
                SD_std = [];  # standard deviation of snow depth
                SIT_final = [];  # Sea Ice Thickness [m]
                SIT_std = [];  # sea Ice thickness standard deviation
                freeboard_final = [];  # freeboard [m]
                freeboard_std = [] # freeboard standard deviation [m]
                sur_temp_final = [];  # surface temperature [celcius]
                air_temp_final = [];  # air temperature [celcius]
                w_SD_final = [];  # snow depth from climatology (Warren)
                w_density_final = [];  # snow density from climatology (Warren)
                pp_flag = 'no'  # pre-processing flag
                unc_flag = 'no'  # Uncertainty flag
                
                #define obsID
                ID_numbers = re.findall(r'\d+', ifile)  # extract numbers from title
                ID_date = [int(num) for num in ID_numbers if int(num) > 100]  # Get date
                
                # Specify that data is from IDCS4   
                obsID = 'OIB_IDCS4_'  + str(ID_date[0])
                
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
                
                ## Convert date format into python date-time format    
                date_string = date.astype(str)
                t = np.array([dt.datetime.strptime(s, "%Y%m%d")
                                 for s in date_string])
                dates = [t[i] + dt.timedelta(seconds=time[i]) for i in range(len(time))] 
                
                
                # Create EASEgrid and return grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('S')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # Take the time for each grid cell into account and calculate 25 km averages
                (avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFrb,
                 stdFrb, lnFrb, uncFrb) = G.GridData(dtint, latitude, longitude, OIBSnow, OIBSnowUnc,
                                                     SIT, SIT_unc, dates, meanFrb, frbUnc)
               
                # Sort data according to date
                index_array = np.zeros(np.shape(time)) 
                time_sorted  =  []
                avgSD_sorted =  []
                stdSD_sorted =  []
                lnSD_sorted  =  []
                uncSD_sorted =  []
                avgSIT_sorted = []
                stdSIT_sorted = []
                lnSIT_sorted  = []
                uncSIT_sorted = []
                avgFrb_sorted = []
                stdFrb_sorted = []
                lnFrb_sorted =  []
                uncFrb_sorted = []
                lat_sorted   =  []
                lon_sorted   =  []
                
                for kk in range(len(time)):
                    index_array = np.argsort(time)
                    time_sorted = np.append(time_sorted,time[index_array[kk]])
                    avgSD_sorted = np.append(avgSD_sorted,avgSD[index_array[kk]])
                    stdSD_sorted = np.append(stdSD_sorted,stdSD[index_array[kk]])
                    lnSD_sorted = np.append(lnSD_sorted,lnSD[index_array[kk]])
                    uncSD_sorted = np.append(uncSD_sorted,uncSD[index_array[kk]])
                    avgSIT_sorted = np.append(avgSIT_sorted,avgSIT[index_array[kk]])
                    stdSIT_sorted = np.append(stdSIT_sorted,stdSIT[index_array[kk]])
                    lnSIT_sorted = np.append(lnSIT_sorted,lnSIT[index_array[kk]])
                    uncSIT_sorted = np.append(uncSIT_sorted,uncSIT[index_array[kk]])
                    avgFrb_sorted = np.append(avgFrb_sorted,avgFrb[index_array[kk]])
                    stdFrb_sorted = np.append(stdFrb_sorted,stdFrb[index_array[kk]])
                    lnFrb_sorted = np.append(lnFrb_sorted,lnFrb[index_array[kk]])
                    uncFrb_sorted = np.append(uncFrb_sorted,uncFrb[index_array[kk]])
                    lat_sorted = np.append(lat_sorted,lat[index_array[kk]])
                    lon_sorted = np.append(lon_sorted,lon[index_array[kk]])
                
                
                # Change names to correct format names
                lat_final=lat_sorted
                lon_final=lon_sorted
                for ll in range(np.size(avgFrb_sorted,0)):
                    date_final=np.append(date_final,dt.datetime.strftime(time_sorted[ll],"%Y-%m-%dT%H:%M:%S"))
                SD_final = avgSD_sorted
                SD_std = stdSD_sorted
                SD_ln = lnSD_sorted
                SD_unc = uncSD_sorted
                
                SIT_final = avgSIT_sorted
                SIT_std = stdSIT_sorted
                SIT_ln = lnSIT_sorted
                SIT_unc = uncSIT_sorted
                
                freeboard_final = avgFrb_sorted
                freeboard_std = stdFrb_sorted
                Frb_ln = lnFrb_sorted
                Frb_unc = uncFrb_sorted
                
                
                # Add nan values to nonexisting data
                if len(SD_final)==0:
                    SD_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_std)==0:
                    SD_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(SIT_final)==0:
                    SIT_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SIT_std)==0:
                    SIT_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(freeboard_final)==0:
                    freeboard_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(freeboard_std)==0:
                    freeboard_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(sur_temp_final)==0:
                    sur_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(air_temp_final)==0:
                    air_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(w_SD_final)==0:
                    w_SD_final=np.zeros(np.shape(lat_final)) * np.nan
                if len(w_density_final)==0:
                    w_density_final=np.zeros(np.shape(lat_final)) * np.nan   
    
                if count == 1: # make empty output file with header
                    output = open(ofile,'w')
                    # print header
                    print('{:^20s} {:^20s} {:^8s} {:^8s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^7s} {:^4s} {:^7s} {:^7s}'.format('obsID', 'date', 'lat', 'lon', 'SD', 'SDstd', 'SDln', 'SDunc', 'SIT', 'SITstd', 'SITln','SITunc', 'Frb', 'FRBstd', 'FRBln', 'FRBunc', 'Tsur', 'Tair', 'wsd', 'w\u03C1', 'pp-flag', 'unc-flag'),file=output)
                
                for ll in range(len(lat_final)):
                    if Frb_ln[ll] != 0 or SIT_ln[ll] != 0 or SD_ln[ll] != 0: # print output whenever Frb, SD or SIT data exist
                        print('{:^20s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], freeboard_final[ll],freeboard_std[ll], Frb_ln[ll], Frb_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag),file=output)
                        print('{:^20s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], freeboard_final[ll],freeboard_std[ll], Frb_ln[ll], Frb_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag))
                
                
output.close()
                    