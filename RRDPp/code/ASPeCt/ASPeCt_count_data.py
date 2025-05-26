# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:46:03 2024

@author: Ida Olsen

Uses EASE-grid to produce 50 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files of visual observations from ships (ice breakers) in the Southern Hemisphere
includes Warren snow depths and densities
"""
# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2024-08-04'

# -- Build-in modules -- #
import os.path
import sys
import datetime as dt

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions

#%% Functions

def OBSID_names(SITln, obsID):
    OBSID_out = [obsID[el-1] for el in np.cumsum(SITln).astype(int)]
    return OBSID_out

def Get_obsID(file):

    obsID = []
    u=0
    #uses u as index to seperate files
    for line in file:
        if line.strip().startswith('%') and  'ASPECT' in line:
            ID = line.strip().replace('%','')
            u=1
        elif u>0:
            obsID.append(ID)
    return obsID

def Make_Gridded_Product(SD, SIT):
        
    # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation 
    G = EASEgrid.Gridded()
    G.SetHemisphere('S')
    G.CreateGrids(gridres)
    (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
    
    SD_unc, SIT_unc = Functions.Get_unc(SD, SIT)
    # Takes the time for each grid cell into account and calculate averages
    avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, var1, var2 = G.GridData(dtint, latitude, longitude, dates, SD=SD, SD_unc=SD_unc, SIT=SIT, SIT_unc=SIT_unc)
    
    if len(time)>0:
        try:
            Functions.plot(lat, lon, dataOut.obsID, time,saveplot, HS='SH')
            Functions.scatter(dataOut.obsID, dates, SD, time, avgSD, 'SD [m]', saveplot)
            Functions.scatter(dataOut.obsID, dates, SIT, time, avgSIT, 'SIT [m]',saveplot)
        except:
            Functions.plot(lat, lon, 'ASPeCt', time,saveplot, HS='SH')
            Functions.scatter('ASPeCt', dates, SD, time, avgSD, 'SD [m]', saveplot)
            Functions.scatter('ASPeCt', dates, SIT, time, avgSIT, 'SIT [m]',saveplot)
    
    #Change names to correct format names
    dataOut.obsID = [dataOut.obsID]*len(lat)
    ## pp flag + unc flag
    dataOut.pp_flag = [dataOut.pp_flag]*len(lat)
    dataOut.unc_flag = [dataOut.unc_flag]*len(lat)
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

#%% MAIN

total_obs_valid_SIT = 0
total_obs_valid_SD = 0
#total_obs_out_SIT = 0
#total_obs_out_SD = 0

# Reads input data
gridres = 50000  # grid resolution
dtint = 30  # days per mean
ifile = 'ASPECT_Allvoys_obs_mindist0.txt'
ifile2 = 'ASPeCt-ASSIST__standardized__ShipBasedSeaIceObservations__SH__UHAM-ICDC_v2.0_fv0.01.txt'

file = open(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/Antarctic/ASPeCt/Original_data/'+ifile, 'r')
save_path_data = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/ASPeCt/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
ofile = save_path_data + 'ESACCIplus-SEAICE-RRDP2+-SIT-ASPeCt-new2.nc'
saveplot =  os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/Antarctic/ASPeCt/fig/'

names =['Year','Mon','DoY','Lon','Lat','Tcc', 'ccP','ICP','ZiP','ZrP','A%rP','aHtP','SCP','SzP','ccS','ICS','ZiS','ZrS','A%rS','aHtS','SCS','SzS','ccT','ICT','ZiT','ZrT','A%rT','aHtT','SCT','SzT']
dtype = [float for el in names]
#     if i<2:
obsID =  Get_obsID(file)
file = open(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/Antarctic/ASPeCt/Original_data/'+ifile, 'r')
data = np.genfromtxt(file, skip_header=22, names=names, comments='%', encoding=None, dtype=dtype)

# seperate observations
index = [i+1 for i in range(len(obsID)-1) if obsID[i]!=obsID[i+1]]
index.append(len(obsID))
index.insert(0,0)

count = 0
for in0, in1 in zip(index[:-1], index[1:]):
    
    count +=1
    #defines variables in output file
    dataOut = Functions.Final_Data(Type='SIT', count_head=count)
    
    datasubset = data[in0:in1]
    dataOut.obsID = obsID[in0]
    
    dates = [dt.datetime(int(y), 1, 1) + dt.timedelta(days=d) for y,m,d in zip(datasubset['Year'], datasubset['Mon'], datasubset['DoY'])]
    
    # use only subset of dataset to avoid duplicates from the second dataset
    index = [True if d<dt.datetime(2002, 5, 8) else False for d in dates]
    
    datasubset = datasubset[index]
    dates = np.array(dates)[index]
    
    # Load datasubset
    latitude = datasubset['Lat']  # degrees
    longitude = datasubset['Lon']  # degrees
    longitude = np.array([l if l<=180 else l-360 for l in longitude])
    ice_conc_tot = datasubset['Tcc']  # total concentration precentage value between 0 and 10

    # Ice type:Primary
    cc_P = np.array(datasubset['ccP']) #precentage value between 0 and 10
    SIT_P = datasubset['ZrP'] # SIT [cm]
    SD_P = datasubset['SzP'] # SD [cm]
    # Ice type:Secondary
    cc_S = np.array(datasubset['ccS']) #precentage secondary ice concentration value between 0 and 10
    SIT_S = datasubset['ZrS'] # SIT [cm]
    SD_S = datasubset['SzS'] # SD [cm]
    # Ice type:Tertiary
    cc_T = np.array(datasubset['ccT']) #precentage value between 0 and 10
    SIT_T = datasubset['ZrT'] # SIT [m]
    SD_T = datasubset['SzT'] # SD [m]

    #Change missing data to nan
    SD_P[SD_P==-9.99] = np.nan
    SD_S[SD_S==-9.99] = np.nan
    SD_T[SD_T==-9.99] = np.nan

    SIT_P[SIT_P==-9.99] = np.nan
    SIT_S[SIT_S==-9.99] = np.nan
    SIT_T[SIT_T==-9.99] = np.nan

    SD, SIT = Functions.compute_SD_SIT(ice_conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T,SIT_T, SD_T)
    if len(SIT)>0 or len(SD)>0:
        total_obs_valid_SIT += len(SIT[np.isfinite(SIT)])
        total_obs_valid_SD += len(SD[np.isfinite(SD)])

    #Make_Gridded_Product(SD, SIT)


#%%
## Read information from newer ASPeCt data file (use data from september 2005, when Worby data ends)
file = open(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) +  '/RRDPp/RawData/Antarctic/ASPeCt/Original_data/'+ifile2, 'r')

#data = np.genfromtxt(file, names=True, encoding=None, delimiter=',')
dtypes = [object if i<1 or i==3 else float for i in range(len(data.dtype.names))]
#file = open(parrent + '/RawData/Antarctic/ASPeCt/Original_data/'+ifile2, 'r')
data = np.genfromtxt(file, names=True, encoding=None, delimiter=',', dtype=dtypes)

dataOut = Functions.Final_Data(Type='SIT', count_head=-1)
#print(data['identifier'])
dataOut.obsID = [idd.decode('utf8') for idd in data['identifier']]

# Load data
date = data[data.dtype.names[0]]
dates = np.array([dt.datetime.strptime(s.decode('utf8'), "%Y-%m-%dT%H:%M:%SZ")
                                  for s in date])
latitude = data['latitude_degrees']  # degrees
longitude = data['longitude_degrees']  # degrees
longitude = np.array([l if l<=180 else l-360 for l in longitude])
ice_conc_tot = data['SIC_total_percent']  # total concentration precentage value between 0 and 10

# Ice type:Primary
cc_P = np.array(data['SIC_primary_percent']) #precentage value between 0 and 10
SIT_P = data['SIT_primary_m'] # SIT [cm]
SD_P = data['snow_depth_primary_m'] # SD [cm]
# Ice type:Secondary
cc_S = np.array(data['SIC_secondary_pct']) #precentage secondary ice concentration value between 0 and 10
SIT_S = data['SIT_secondary_m'] # SIT [cm]
SD_S = data['snow_depth_secondary_m'] # SD [cm]
# Ice type:Tertiary
cc_T = np.array(data['SIC_tertiary_pct']) #precentage value between 0 and 10
SIT_T = data['SIT_tertiary_m'] # SIT [m]
SD_T = data['snow_depth_tertiary_m'] # SD [m]

#Change missing data to nan
ice_conc_tot[ice_conc_tot==-9.] = np.nan
cc_P[cc_P==-9.] = np.nan
cc_S[cc_S==-9.] = np.nan
cc_T[cc_T==-9.] = np.nan

SD_P[SD_P==-9.9] = np.nan
SD_S[SD_S==-9.9] = np.nan
SD_T[SD_T==-9.9] = np.nan

SIT_P[SIT_P==-9.9] = np.nan
SIT_S[SIT_S==-9.9] = np.nan
SIT_T[SIT_T==-9.9] = np.nan

SD, SIT = Functions.compute_SD_SIT(ice_conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T,SIT_T, SD_T)
total_obs_valid_SIT += len(SIT[np.isfinite(SIT)])
total_obs_valid_SD += len(SD[np.isfinite(SD)])

#print(dates[0])
print(total_obs_valid_SD)
print(total_obs_valid_SIT)

print(dates[-1])