# -*- coding: utf-8 -*-

"""
Creates 25 km gridded means of input files containing measurements of SIT and SD from 
input files measured by ice breakers from the Antarctic Sea Ice Processes and Climate (ASSIST)
Uses EASE-grid to produce 25 km grid mean values.
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2020-07-13'

# -- Built-in modules -- #
import os.path
import datetime as dt
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

#%% Functions

def compute_SD_SIT(conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T,SIT_T, SD_T):
    SD = [];
    SIT = [];
    for kk in range(len(conc_tot)):
        # print(kk)
        if conc_tot[kk]>0 and ~np.isnan(conc_tot[kk]):    
    
            SIT_eff_P = np.multiply(np.divide(cc_P[kk],conc_tot[kk]),SIT_P[kk])
            SIT_eff_S = np.multiply(np.divide(cc_S[kk],conc_tot[kk]),SIT_S[kk])
            SIT_eff_T = np.multiply(np.divide(cc_T[kk],conc_tot[kk]),SIT_T[kk])
    
            SD_eff_P = np.multiply(np.divide(cc_P[kk],conc_tot[kk]),SD_P[kk])
            SD_eff_S = np.multiply(np.divide(cc_S[kk],conc_tot[kk]),SD_S[kk])
            SD_eff_T = np.multiply(np.divide(cc_T[kk],conc_tot[kk]),SD_T[kk])
            
            # Append non nan value if:
            # 1. The total concentration is the same as the sum of the partial concentrations
            # 2. As a minimum one of the SIT/SD_eff are non nan
            try:
                assert ~np.isnan(cc_P[kk]) or ~np.isnan(cc_S[kk]) or ~np.isnan(cc_T[kk])
                assert np.nansum([cc_P[kk],cc_S[kk],cc_T[kk]]) == conc_tot[kk]
                assert ~np.isnan(SIT_eff_P) or ~np.isnan(SIT_eff_S) or ~np.isnan(SIT_eff_T)
                SIT = np.append(SIT,np.nansum([SIT_eff_P,SIT_eff_S,SIT_eff_T]))
            except AssertionError:
                SIT = np.append(SIT, np.nan)
            try:
                assert ~np.isnan(cc_P[kk]) or ~np.isnan(cc_S[kk]) or ~np.isnan(cc_T[kk])
                assert np.nansum([cc_P[kk],cc_S[kk],cc_T[kk]]) == conc_tot[kk]
                assert ~np.isnan(SD_eff_P) or ~np.isnan(SD_eff_S) or ~np.isnan(SD_eff_T)
                SD = np.append(SD,np.nansum([SD_eff_P,SD_eff_S,SD_eff_T]))
            except AssertionError:
                SD = np.append(SD, np.nan)
     
        elif conc_tot[kk]==0: # Everything is water      
            SIT = np.append(SIT,0)
            SD = np.append(SD,0)
        else:
            SIT = np.append(SIT,np.nan)
            SD = np.append(SD,np.nan)
    return SD,SIT

def Get_unc(SD, SIT):
    SD_unc = []
    SIT_unc = []
    for sit, sd in zip(SIT,SD):
        if np.isnan(sit):
            SIT_unc.append(np.nan)
        elif sit<=0.10: # m
            SIT_unc.append(sit*0.5)
        elif sit<=0.3: # m
            SIT_unc.append(sit*0.3)
        elif sit>0.3:
            SIT_unc.append(sit*0.2)    
        if np.isnan(sd):
            SD_unc.append(np.nan)
        elif sd<=0.10: # m
            SD_unc.append(sd*0.5)
        elif sd<=0.3: # m
            SD_unc.append(sd*0.3)
        elif sd>0.3:
            SD_unc.append(sd*0.2)
    return SD_unc, SIT_unc 

#%% Main

parrent = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))
# Reads input data
directory = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) + '/RRDPp/RawData/ASSIST'
# saving location
save_path_data=parrent + '/RRDPp/FINAL/ASSIST/final/'
saveplot = parrent + '/RRDPp/FINAL/ASSIST/fig/'
ofile ='ESACCIplus-SEAICE-RRDP2+-SIT-ASSIST-t2.nc'
ofile = os.path.join(save_path_data,ofile)

gridres = 25000  # grid resolution
dtint = 30; # days per mean

# Access data
directories = os.listdir(directory)
directories = [d for d in directories if d.startswith('2')]
# sort to ensure starting with the earliest date
# numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
# index = np.argsort(numbers)
# directories = [directories[ind] for ind in index]

# count for header
#total_obs = 0
total_obs_valid_SIT = 0
total_obs_valid_SD = 0
total_obs_out_SIT = 0
total_obs_out_SD = 0
count = 0
for dir in directories:
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all ASSIST data
            if file.endswith('.csv') and dir:
                
                count +=1
                #defines variables in output file
                dataOut = Functions.Final_Data(Type='SIT', count_head=count)

                # Reads observation data from ASCII-file
                dtype = [object if i<2 else float for i in range(70)]
                data = np.genfromtxt(file,skip_header=0,names=True, dtype=dtype,delimiter=',',usecols=np.arange(0,70), encoding=None)
                
                # Load data
                date = data['Date']
                # faulty date measurement
                if date[0].decode('utf8').startswith(''):
                    data = np.delete(data, 0, 0)
                    date = data['Date']
                latitude = data['LAT']  # degrees
                longitude = data['LON']  # degrees
                conc_tot = data['TC']  # total concentration precentage value between 0 and 10
                
                # Define obsID
                obsName = ifile.split('.')[0]
                file_num = obsName.split('-')[1]
                dataOut.obsID = 'ASSIST_' + 'obs' + file_num 
                
                # Ice type:Primary
                cc_P = data['PPC'] #precentage value between 0 and 10
                SIT_P = data['PZ'] # SIT [cm]
                SD_P = data['PSH'] # SD [cm]
                # Ice type:Secondary
                cc_S = np.array(data['SPC']) #precentage secondary ice concentration value between 0 and 10
                SIT_S = data['SZ'] # SIT [cm]
                SD_S = data['SSH'] # SD [cm]
                # Ice type:Tertiary
                cc_T = np.array(data['TPC']) #precentage value between 0 and 10
                SIT_T = data['TZ'] # SIT [m]
                SD_T = data['TSH'] # SD [m]
                
                
                #change errorcode for conc from -1 to nan
                conc_tot[conc_tot==-1] = np.nan
                cc_P[cc_P==-1] = np.nan
                cc_S[cc_S==-1] = np.nan
                cc_T[cc_T==-1] = np.nan
                
                SIT_P[SIT_P==-1] = np.nan
                SIT_S[SIT_S==-1] = np.nan
                SIT_T[SIT_T==-1] = np.nan
                
                SD_P[SD_P==-1] = np.nan
                SD_S[SD_S==-1] = np.nan
                SD_T[SD_T==-1] = np.nan
                
                # Get weighted average
                cc_tot = np.add(np.add((cc_P),(cc_S)),(cc_T))
                
                SD, SIT = Functions.compute_SD_SIT(conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T,SIT_T, SD_T)
                
                # Change measurements from cm to m
                SIT = SIT/100
                SD = SD/100
                total_obs_valid_SIT += len(SIT[np.isfinite(SIT)])
                total_obs_valid_SD += len(SD[np.isfinite(SD)])
                #print(len(SIT[np.isfinite(SIT)]))
                
                t = np.array([dt.datetime.strptime(s.decode('utf8'), "%Y-%m-%d %H:%M:%S UTC")
                                 for s in date])
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # # Define uncertainties
                # SD_unc = np.ones(len(latitude)) * 0.20 # precision of observations (approximate uncertainty)
                # SD_unc[np.isnan(SD)] = np.nan
                # SIT_unc = np.ones(len(latitude)) * 0.20
                # SIT_unc[np.isnan(SIT)] = np.nan
                
                SD_unc, SIT_unc = Functions.Get_unc(SD, SIT)
                # convert to numpy arrays
                SD_unc = np.array(SD_unc)
                SIT_unc = np.array(SIT_unc)
                # Takes the time for each grid cell into account and calculate averages
                # (avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB,
                #  lnFRB, FRB_Unc,var1, var2) = G.GridData(dtint, latitude, longitude, t, SD, SD_unc, SIT, SIT_unc, FRB=[], FRB_unc=[])
                
                # total_obs_out_SIT += len(avgSIT[np.isfinite(avgSIT)])
                # total_obs_out_SD += len(avgSD[np.isfinite(avgSD)])

                
print(total_obs_valid_SD)
print(total_obs_valid_SIT)
# print(total_obs_out_SIT)
# print(total_obs_out_SD)
                