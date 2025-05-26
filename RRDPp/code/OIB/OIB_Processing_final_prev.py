# -*- coding: utf-8 -*-
"""
Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by NASA's Operation Icebridge
includes Warren snow depths and densities

Choose whether data from the southern or the northern hemisphere is examined
"""
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '1'
__date__ = '2021-11-15'

# -- Built-in modules -- #
import os.path
import datetime as dt
import re
import sys

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from Warren import SnowDepth, SWE
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

# SET THIS PARAMETER
hemisphere = 'SH'

if hemisphere == 'NH':
    # determine gridsize and temporal limits
    gridres = 25000  # grid resolution
    # Reads input data
    save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/OIB/final/'
    saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/OIB/fig/'
    ofile = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/OIB/final/ESACCIplus-SEAICE-RRDP2+-SIT-OIB-QL.dat'
    directory = os.path.dirname(os.path.dirname(os.getcwd())) + '/RawData/OIB/QuickLook'
elif hemisphere == 'SH':
    # determine gridsize and temporal limits
    gridres = 50000  # grid resolution
    # Reads input data
    save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/OIB/final/'
    saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/OIB/fig/'
    ofile = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/OIB/final/ESACCIplus-SEAICE-RRDP2+-SIT-OIB.dat'
    directory = os.path.dirname(os.path.dirname(os.getcwd())) + '/RawData/Antarctic/OIB'

dtint = 30  # days per mean
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
                
                # Reads observation data from ASCII-file
                data = np.genfromtxt(file,dtype=None,names=True,delimiter=',',usecols=np.arange(0,50)) #usecols: specifies which collumns to read
                
                #defines variables in output file
                count+=1
                dataOut = Functions.Final_Data(Type='FRB', count_head=count)
                dataOut.pp_flag = 0
                dataOut.unc_flag = 0
                
                #define obsID
                ID_numbers = re.findall(r'\d+', ifile)  # extract numbers from title
                ID_date = [int(num) for num in ID_numbers if int(num) > 100]  # Get date
                
                # Determine if data is from QL or IDCS4
                if 'QuickLook' in directory: 
                    dataOut.obsID = 'OIB_QL_' + str(ID_date[0])
                elif 'IDCS4' in directory:    
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
                
                ## Convert date format into python date-time format    
                date_string = date.astype(str)
                t = np.array([dt.datetime.strptime(s, "%Y%m%d")
                                 for s in date_string])
                dates = [t[i] + dt.timedelta(seconds=time[i]) for i in range(len(time))] 
                
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
                G = EASEgrid.Gridded()
                G.SetHemisphere(hemisphere[:-1])
                G.CreateGrids(gridres)
                (index_i, index_j) = G.LatLonToIdx(latitude, longitude)
    
                # Take the time for each grid cell into account and calculate 25 km averages
                (avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB,
                 stdFRB, lnFRB, uncFRB) = G.GridData(dtint, latitude, longitude, dates, SD=OIBSnow, SD_unc=OIBSnowUnc,
                                                     SIT=SIT, SIT_unc=SIT_unc, FRB=meanFrb, FRB_unc=frbUnc)
    
                if len(time) > 0:
                    Functions.plot(latitude, longitude, dataOut.obsID, time,saveplot, HS=hemisphere)
                    Functions.scatter(dataOut.obsID, dates, meanFrb, time, avgFRB, 'FRB [m]', saveplot)
                    Functions.scatter(dataOut.obsID, dates, SIT, time, avgSIT, 'SIT [m]', saveplot)
                    Functions.scatter(dataOut.obsID, dates, OIBSnow, time, avgSD, 'SD [m]', saveplot)


                if hemisphere == 'NH':
                    # Correlates OIB data with Warren snow depth and snow density
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
                    
                # print data to output file
                dataOut.Print_to_output(ofile, primary='FRB')
        
Functions.sort_final_data(ofile, saveplot=saveplot, HS=hemisphere)
