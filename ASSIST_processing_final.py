# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 12:44:44 2020

@author: Olsen
@author: s174020

Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by ice breakers from the Antarctic Sea Ice Processes and Climate
Does not includes Warren snow depth and densities
"""
import os.path
#import math
import pdb
import numpy as np
import sys
import datetime as dt
import re
import EASEgrid_ASSIST as EASEgrid

# Reads input data
directory = '/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/ASSIST'
save_path_data='/media/s174020/My Passport2/ESACCIplus/RRDPp/FINAL/ASSIST/final/'

gridres = 25000 #grid resolution
dtint = 30; # days per mean

ifile = '/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/ASSIST/observations-92.csv'
ofile ='ESACCIplus-SEAICE-RRDP2+-SIT-ASSIST.dat'

# saving locations
ofile = os.path.join(save_path_data,ofile)

count_head = 0  # used to locate header

# Access data
directories = os.listdir(directory)
# sort to ensure starting with the earliest date
numbers = [re.findall(r'[0-9]+', directory) for directory in directories]
index = np.argsort(numbers)
directories = [directories[ind] for ind in index]

for dir in directories:
    # dir = '2021'
    dir_data = os.path.join(directory, dir)                
    if os.path.isdir(dir_data):          
        files = os.listdir(dir_data)
        files.sort()
        for ifile in files:
            file = os.path.join(dir_data, ifile)
            # Access all SIT files in AEM_AWI
            if file.endswith('.csv') and dir:
                print(file)
                count_head += 1
                
                #defines variables in output file
                flag = 0
                obsID = [];
                date_final = []; 
                lat_final = []; # decimal degrees
                lon_final = []; # decimal degrees
                SD_final = [];  # snow depth [m]
                SD_std = [];  # standard deviation of snow depth
                SD_ln = []   # number of SD measurements in each gridcell
                SD_unc = []  # SD uncertainty
                SIT_final = [];  # Sea Ice Thickness [m]
                SIT_std = [];  # sea Ice thickness standard deviation
                SIT_ln =[]
                SIT_unc = []
                FRB_final = [];  # freeboard [m]
                FRB_std = [] # freeboard standard deviation [m]
                FRB_ln = []  # number of FRB measurements in each gridcell
                FRB_unc = []  # FRB uncertainty
                sur_temp_final = [];  # surface temperature [celcius]
                air_temp_final = [];  # air temperature [celcius]
                w_SD_final = [];  # snow depth from climatology (Warren)
                w_density_final = [];  # snow density from climatology (Warren)
                pp_flag = 'yes'  # pre-processing flag
                unc_flag = 'yes'  # Uncertainty flag


                # Reads observation data from ASCII-file
                data = np.genfromtxt(file,dtype=None,skip_header=0,names=True,delimiter=',',usecols=np.arange(0,110))
                
                # Load data
                date = data['Date']
                latitude = data['LAT']  # degrees
                longitude = data['LON']  # degrees
                ice_conc_tot = data['TC']  # total concentration precentage value between 0 and 10
                
                # Define obsID
                obsName = ifile.split('.')[0]
                file_num = obsName.split('-')[1]
                try:
                    date0 = date[0].decode('utf8')
                except:
                    date0 = date[0]
                IDdate = date0.split(' ')[0].replace('-','')
                obsID = 'ASSIST_' + 'obs' + file_num + '_' + IDdate
                
                # Ice type:Primary
                icc_P = np.array(data['PPC']) #precentage value between 0 and 10
                SIT_prim = data['PZ'] #[cm]
                SD_P = data['PSH'] #[cm]
                # Ice type:Secondary
                icc_S = np.array(data['SPC']) #precentage secondary ice concentration value between 0 and 10
                SIT_sec = data['SZ'] #[cm]
                SD_S = data['SSH'] #[cm]
                # Ice type:Tertiary
                icc_T = np.array(data['TPC']) #precentage value between 0 and 10
                SIT_tert = data['TZ'] #[m]
                SD_T = data['TSH'] #[m]
                
                # define variables
                conc_tot = [];
                cc_P = [];
                cc_S = [];
                cc_T = [];
                SIT_P = [];
                SIT_S = [];
                SIT_T = [];
                
                SD_prim = [];
                SD_sec = [];
                SD_tert = [];
                
                SD = [];
                SIT = [];
                
                for num in range(len(ice_conc_tot)):
                
                    conc_tot = np.append(conc_tot,float(ice_conc_tot[num]))    
                    #for num in icc_P: 
                    cc_P = np.append(cc_P,float(icc_P[num]))
                    cc_S = np.append(cc_S,float(icc_S[num]))
                    cc_T = np.append(cc_T,float(icc_T[num]))
                    SIT_P = np.append(SIT_P,float(SIT_prim[num]))
                    SIT_S = np.append(SIT_S,float(SIT_sec[num]))
                    SIT_T = np.append(SIT_T,float(SIT_tert[num]))
                    #convet integers into float number
                    SD_prim = np.append(SD_prim,float(SD_P[num]))
                    SD_sec = np.append(SD_sec,float(SD_S[num]))
                    SD_tert = np.append(SD_tert,float(SD_T[num]))
                
                #change errorcode for conc from -1 to nan
                conc_tot[conc_tot==-1] = np.nan
                cc_P[cc_P==-1] = np.nan
                cc_S[cc_S==-1] = np.nan
                cc_T[cc_T==-1] = np.nan
                
                SIT_P[SIT_P==-1] = np.nan
                SIT_S[SIT_S==-1] = np.nan
                SIT_T[SIT_T==-1] = np.nan
                
                SD_prim[SD_prim==-1] = np.nan
                SD_sec[SD_sec==-1] = np.nan
                SD_tert[SD_tert==-1] = np.nan
                
                #weighted average - ice thickness
                cc_P = np.array(cc_P)
                cc_S = np.array(cc_S)
                cc_T = np.array(cc_T)
                cc_tot = np.add(np.add((cc_P),(cc_S)),(cc_T))
                
                for kk in range(len(cc_tot)):
                    print(kk)
                    if cc_tot[kk]>0:    
                        #calculate weight
                        SIT_eff_P = np.multiply(np.divide(cc_P[kk],cc_tot[kk]),SIT_P[kk])
                        SIT_eff_S = np.multiply(np.divide(cc_S[kk],cc_tot[kk]),SIT_S[kk])
                        SIT_eff_T = np.multiply(np.divide(cc_T[kk],cc_tot[kk]),SIT_T[kk])
                        #estimated total ice thickness
                        SIT = np.append(SIT,SIT_eff_P+SIT_eff_S+SIT_eff_T)
                    
                        #average snow thickness? Wighted average
                        SD_eff_P = np.multiply(np.divide(cc_P[kk],cc_tot[kk]),SD_prim[kk])
                        SD_eff_S = np.multiply(np.divide(cc_S[kk],cc_tot[kk]),SD_sec[kk])
                        SD_eff_T = np.multiply(np.divide(cc_T[kk],cc_tot[kk]),SD_tert[kk])
                        
                        SD = np.append(SD,SD_eff_P+SD_eff_S+SD_eff_T)
                
                    else:       
                        SIT = np.append(SIT,0)
                        SD = np.append(SD,0)
                
                #change measurements from cm to m
                SIT = SIT/100
                SD = SD/100
                
                # Reads date format into datetime format
                date = date.astype(str)
                
                if date[0].startswith(''):
                    date = np.delete(date,0)
                    latitude = np.delete(latitude,0)
                    longitude = np.delete(longitude,0)
                    SD = np.delete(SD,0)
                    SIT = np.delete(SIT,0)
                
                t = np.array([dt.datetime.strptime(s, "%Y-%m-%d %H:%M:%S UTC")
                                 for s in date])
                
                # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation 
                G = EASEgrid.Gridded()
                G.SetHemisphere('N')
                G.CreateGrids(gridres)
                (index_i,index_j) = G.LatLonToIdx(latitude,longitude)
                
                # Define uncertainties
                SD_unc = np.ones(len(latitude)) * 0.20 # precision of observations (approximate uncertainty)
                SD_unc[np.isnan(SD)] = np.nan
                SIT_unc = np.ones(len(latitude)) * 0.20
                SIT_unc[np.isnan(SIT)] = np.nan
                # Takes the time for each grid cell into account and calculate averages
                (avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB,
                 lnFRB, FRB_Unc,) = G.GridData(dtint, latitude, longitude, SD, SD_unc, SIT, SIT_unc, t, Frb=[], Frb_unc=[])
                
                
                # Sort data according to date
                index_array = np.zeros(np.shape(time))
                time_sorted  = []
                avgSD_sorted = []
                stdSD_sorted = []
                lnSD_sorted  = []
                uncSD_sorted = []
                avgSIT_sorted = []
                stdSIT_sorted = []
                lnSIT_sorted  = []
                uncSIT_sorted = []
                lat_sorted   = []
                lon_sorted   = []
                
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
                    lat_sorted = np.append(lat_sorted,lat[index_array[kk]])
                    lon_sorted = np.append(lon_sorted,lon[index_array[kk]])
                
                #Correlate obs data with Warren snow depth and snow density
                from Warren import SnowDepth, SWE
                w_SD_final      = []
                w_density_final = []
                for ll in range(np.size(avgSD_sorted,0)):
                    (w_SD,w_SD_epsilon) = SnowDepth(lat_sorted[ll],lon_sorted[ll],time_sorted[ll].month)
                    w_SD_final = np.append(w_SD_final,w_SD)
                    (wswe,wswe_epsilon) = SWE(lat_sorted[ll],lon_sorted[ll],time_sorted[ll].month)
                    w_density=int((wswe/w_SD)*1000)
                    w_density_final = np.append(w_density_final,w_density)
                
                #Change names to correct format names
                lat_final = lat_sorted
                lon_final = lon_sorted
                for ll in range(np.size(avgSD_sorted,0)):
                    date_final = np.append(date_final,dt.datetime.strftime(time_sorted[ll],"%Y-%m-%dT%H:%M:%S"))
                SD_final = avgSD_sorted
                SD_std = stdSD_sorted
                SD_ln = lnSD_sorted
                SD_unc = uncSD_sorted
                SIT_final = avgSIT_sorted
                SIT_std = stdSIT_sorted
                SIT_ln = lnSIT_sorted
                SIT_unc = uncSIT_sorted
                
                # Add nan values to nonexisting data
                if len(SD_final)==0:
                    SD_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_std)==0:
                    SD_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_unc)==0:
                    SD_unc = np.zeros(np.shape(lat_final)) * np.nan
                if len(SD_ln) == 0:
                    SD_ln = np.zeros(np.shape(lat_final))
                if len(SIT_final)==0:
                    SIT_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(SIT_std)==0:
                    SIT_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(SIT_ln)==0:
                    SIT_ln = np.zeros(np.shape(lat_final))
                if len(SIT_unc)==0:
                    SIT_unc = np.zeros(np.shape(lat_final)) * np.nan
                if len(FRB_final)==0:
                    FRB_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(FRB_std)==0:
                    FRB_std = np.zeros(np.shape(lat_final)) * np.nan
                if len(FRB_ln) == 0:
                    FRB_ln = np.zeros(np.shape(lat_final))
                if len(FRB_unc) == 0:
                    FRB_unc = np.zeros(np.shape(lat_final)) * np.nan
                if len(sur_temp_final)==0:
                    sur_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                if len(air_temp_final)==0:
                    air_temp_final = np.zeros(np.shape(lat_final)) * np.nan
                    
                # print data to output file
                if count_head == 1: # make empty output file with header
                    output = open(ofile,'w')
                    print('{:^30s} {:^20s} {:^8s} {:^8s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^6s} {:^7s} {:^7s} {:^6s} {:^4s} {:^7s} {:^7s}'.format('obsID', 'date', 'lat', 'lon', 'SD', 'SDstd', 'SDln', 'SDunc', 'SIT', 'SITstd', 'SITln','SITunc', 'Frb', 'FRBstd', 'FRBln', 'FRBunc', 'Tsur', 'Tair', 'wsd', 'w\u03C1', 'pp-flag', 'unc-flag'),file=output)
                output = open(ofile,'a')
                for ll in range(np.size(SD_final,0)):
                    if SIT_ln[ll] != 0 or SD_ln[ll] != 0:
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], FRB_final[ll],FRB_std[ll], FRB_ln[ll], FRB_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag),file=output)
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7s} {:^7s}'.format(obsID,(date_final[ll]),lat_final[ll],lon_final[ll],SD_final[ll],SD_std[ll],SD_ln[ll], SD_unc[ll], SIT_final[ll],SIT_std[ll],SIT_ln[ll], SIT_unc[ll], FRB_final[ll],FRB_std[ll], FRB_ln[ll], FRB_unc[ll], sur_temp_final[ll],air_temp_final[ll],w_SD_final[ll]/100.,w_density_final[ll], pp_flag, unc_flag))
                
                output.close()


