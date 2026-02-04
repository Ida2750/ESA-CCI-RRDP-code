#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Creates 25 km gridded means of input files containing measurements of SIT and SD from 
airborne campaigns using the EM-bird, and from the Snow and Ice Mass Balance Array both
obtained during the MOSAIC campaign
includes Warren snow depths and densities
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2025-04-30'

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
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions

#%%
# specify gridsize and temporal limits
gridres = 25000  # grid resolution
dtint = 30  # days per mean
name = 'MOSAIC'

# Define locations of input and output files
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/{name}/final/'

if not os.path.exists(save_path_data):os.makedirs(save_path_data)
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) +  f'/FINAL/{name}/fig/'

#directory = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) + f'/RRDPp/RawData/{name}'
directory = f'/dmidata/projects/cmems2/C3S/RRDPp/RawData/{name}'

# create directories if they do not exist
if not os.path.exists(save_path_data): os.makedirs(save_path_data)
if not os.path.exists(saveplot): os.makedirs(saveplot)

count = 0  # used to locate first file

# Access HEM data
directories = sorted(os.listdir(directory))


for dirr in directories:
    datalen_SD = 0
    datalen_SIT = 0
    ofile = os.path.dirname(os.path.dirname(os.getcwd())) + f'/FINAL/{name}/final/ESACCIplus-SEAICE-RRDP2+-SIT-{name}-{dirr}.nc'
    files = sorted(os.listdir(f'{directory}/{dirr}'))
    #print(files)
    for ifile in files:
        Snow=False #default as most files does not contain snowdepth measurements
        file = os.path.join(f'{directory}/{dirr}/{ifile}')
        # Access all SIT files in AEM_AWI
        if (file.endswith('.tab') or file.endswith('.dat') or file.endswith('.nc')) and 'free' not in ifile:
            
            #defines variables in output file
            count+=1
            dataOut = Functions.Final_Data(Type='SIT', count_head=count)
            dataOut.pp_flag = 0
            
            lookup = b'*/'  # notation for when header is done
            counthead = 0
            header = -999
            with open(file,'rb') as myFile:
                for num, line in enumerate(myFile, 1):
                    counthead+=1
                    # break loop when lookup is found
                    if lookup in line:
                        header = counthead  # note number of headerlines
                    if counthead == header+1:
                        headerline = line.decode('utf8').replace("[m]","")  #remove [m]
                        names = headerline.split("\t")
                        names = [name.strip() for name in names]
                        break
                
            data = np.genfromtxt(file, dtype=None, delimiter='\t', names=names, skip_header=header+1, encoding=None)
            latitude = data['Latitude']
            longitude = data['Longitude']
            if dirr=='SIMBA':
                Snow=True
                SIT = data['Ice_thick']
                SD = data['Snow_thick']
                SD[SD<0] = np.nan
                SD[SD>2] = np.nan
                # source https://doi.pangaea.de/10.1594/PANGAEA.938244
                SIT_unc = np.ones(len(latitude)) * 0.02
                SD_unc = np.ones(len(latitude)) * 0.02

            elif dirr=='HEM':
                QF = data['QF_Reliability']
                #print(np.unique(QF))
                FIL_MOD = data['Filter_Moderate_filter']
                FIL_STR = data['Filter_Strict_filter']
                SIT = data['EsEs__ice_and_snow_Electromagneti']
                # filter data based on flags
                filter_comb = (QF<=2) & (FIL_MOD==1) & (FIL_STR==1) 
                #print(f' included: {len(filter_comb[filter_comb])}')
                #print(f' removed: {len(filter_comb[~filter_comb])}')
                SIT = SIT[filter_comb]
                latitude = latitude[filter_comb]
                longitude = longitude[filter_comb]

                ## Define fixed uncertainty
                SIT_unc = np.ones(len(latitude)) * 0.10

            # replace faulty values if applicable
            SIT[SIT < 0] = np.nan
            SIT[SIT > 10] = np.nan

            try:
                datalen_SD += len(SD[np.isfinite(SD)])
            except:
                pass
            datalen_SIT += len(SIT[np.isfinite(SIT)])

            date = data['DateTime']
            # Convert date format into python date-time format
            dates = date.astype(str)
            if dirr=='HEM':
                t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                                    for s in dates])
                t = t[filter_comb]
            elif dirr=='SIMBA':
                t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M")
                    for s in dates])
            # Define obsID
            flight = ifile[:8]
            dataOut.obsID = f'{name}_' + flight

            # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
            G = EASEgrid.Gridded()
            G.SetHemisphere('N')
            G.CreateGrids(gridres)
            (index_i, index_j) = G.LatLonToIdx(latitude, longitude)

            # Takes the time for each grid cell into account and calculate averages
            if Snow: # SIMBA
                avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, var1, var2, dataOut.QFT, dataOut.QFS, dataOut.QFG = G.GridData(
                    dtint, latitude, longitude, t, SD=SD, SD_unc=SD_unc, SIT=SIT, SIT_unc=SIT_unc, dtype='buoy')
            else:
                avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, var1, var2, dataOut.QFT, dataOut.QFS, dataOut.QFG = G.GridData(
                    dtint, latitude, longitude, t, SIT=SIT, SIT_unc=SIT_unc, dtype='AEM')


            #print(np.nanmax(uncSIT))  
            dataOut.unc_flag = np.ones(len(avgSIT)).astype(int)*2
            if dirr=='HEM':
                dataOut.unc_flag[avgSIT>3] = 3
                
            if len(time) > 0:
                Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS='NH')
                Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]', saveplot)
                if Snow:
                    Functions.scatter(dataOut.obsID, t, SD, time, avgSD, 'SD [m]', saveplot)


            # Correlates SB-AWI buoy data with Warren snow depth and snow density
            for ll in range(np.size(avgSD,0)):
                (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
                dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
                (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
                w_density=int((wswe/w_SD)*1000)
                dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)

            #Change names to correct format names
            dataOut.obsID = [dataOut.obsID]*len(lat)
            ## pp flag + unc flag
            dataOut.pp_flag = [dataOut.pp_flag]*len(lat)
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
            
            # fill empty arrays with NaN values
            dataOut.Check_Output()
                
            if count>1:
                subset = dataOut.Create_NC_file(ofile,primary='SIT')
                df = Functions.Append_to_NC(df, subset)
            else:
                if dirr=='HEM':
                    df = dataOut.Create_NC_file(ofile,primary='SIT', datasource='Airborne Electromagnetic Measurement, MOSAIC, https://doi.pangaea.de/10.1594/PANGAEA.934578', key_variables='Ice thickness: sea ice thickness + snow depth')
                elif dirr=='SIMBA':
                    print('SIMBA')
                    df = dataOut.Create_NC_file(ofile,primary='SIT', datasource='Snow and Ice Mass Balance Array (SIMBA), MOSAIC, https://doi.pangaea.de/10.1594/PANGAEA.938244', key_variables='sea ice thickness and snow depth')
    
    # Save data to NetCDF
    print(f'sd len {datalen_SD}')
    print(f'sit len {datalen_SIT}')
    Functions.save_NC_file(df, ofile, primary='SIT')

            
