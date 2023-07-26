# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:23:56 2023

@author: Ida Olsen
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2023-06-27'


# -- Built-in modules -- #
import os
import datetime as dt
import re
import sys

# -- Third-part modules -- #
import pandas as pd
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions
from Warren import SnowDepth, SWE


savepath = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/NP/final/'
saveplot = os.path.dirname(os.path.dirname(os.getcwd())) + '/FINAL/NP/fig/'
path = os.path.dirname(os.path.dirname(os.getcwd())) + '/RawData/NP/'
files = [fil for fil in os.listdir(path) if fil.endswith('xlsx') and fil.startswith('NP')]
files_pos = [fil for fil in os.listdir(path + 'pos/') if fil.endswith('.csv')]
ofile = savepath + 'ESACCIplus-SEAICE-RRDP2+-SIT-NP.dat'

def lat_lon2(file):
    data = np.genfromtxt(os.path.join(path + '/pos', file), names=True, delimiter=',')
    lat = data['latitude']
    lon = data['longitude']
    day = data['date']
    index = list(np.where(np.diff(day)==1)[0]+1)
    index.insert(0,0)
    index.append(len(day))
    for ind1, ind2 in zip(index[:-1], index[1:]):
        dataOut.lat_final.append(np.nanmedian(lat[ind1:ind2]))
        dataOut.lon_final.append(np.nanmedian(lon[ind1:ind2]))

def Get_unc(self):
    ## very crude estimation of unc
    sit_unc = [0.10]
    sd_unc = [0.10]
    frb_unc = [0.10]
    
    self.FRB_unc = [1/el*np.sqrt(np.nansum((np.array(frb_unc*el))**2)) if el>0 else np.nan for el in self.FRB_ln]
    self.SIT_unc = [1/el*np.sqrt(np.nansum((np.array(sit_unc*el))**2)) if el>0 else np.nan for el in self.SIT_ln]
    self.SD_unc = [1/el*np.sqrt(np.nansum((np.array(sd_unc*el))**2)) if el>0 else np.nan for el in self.SD_ln]

count=0
for fil, fil_pos in zip(files, files_pos):
    assert fil[:4]==fil_pos[:4]
    df = pd.read_excel( path + fil, header=1, skiprows=[37,38,39], sheet_name='Data_sheet')


    # print(df)
    count += 1
    dataOut = Functions.Final_Data(Type='SIT', count_head=count)
    dataOut.unc_flag = 3
    dataOut.pp_flag = 0
    listt = np.linspace(1,50,50)
    datoer = pd.read_excel( path + fil, skiprows=listt, sheet_name='Data_sheet').keys()
    dates = [el for el in datoer if 'Unnamed' not in str(el)]
    dataOut.date_final = [dt.datetime.strftime(dat,"%Y-%m-%dT%H:%M:%S") for dat in dates]
    dataOut.obsID = 'NP-' + re.findall('\d+', fil)[0]
    
    for key in df.keys():
        if 'Ice_thick' in key and not 'rel' in key or 'ice_thick' in key:
            dataOut.SIT_final.append(np.nanmean(df[key])/100)
            dataOut.SIT_ln.append(len(df[key][~np.isnan(df[key])]))
            dataOut.SIT_std.append(np.nanstd(df[key])/100)
        if 'Freeboard' in key:
            dataOut.FRB_final.append(np.nanmean(df[key])/100)
            dataOut.FRB_ln.append(len(df[key][~np.isnan(df[key])]))
            dataOut.FRB_std.append(np.nanstd(df[key])/100)
        if 'Sn_Depth' in key or 'Sn_depth' in key:    
            dataOut.SD_final.append(np.nanmean(df[key])/100)
            dataOut.SD_ln.append(len(df[key][~np.isnan(df[key])]))
            dataOut.SD_std.append(np.nanstd(df[key])/100)
    
    ## Get positions
    lat_lon2(fil_pos)

    # Correlates data with Warren snow depth and snow density
    for ll in range(np.size(dataOut.SIT_final,0)):
        (w_SD,w_SD_epsilon) = SnowDepth(dataOut.lat_final[ll],dataOut.lon_final[ll],dates[ll].month)
        dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
        (wswe,wswe_epsilon) = SWE(dataOut.lat_final[ll],dataOut.lon_final[ll],dates[ll].month)
        w_density=int((wswe/w_SD)*1000)
        dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)
        
    # fill empty arrays with NaN values
    dataOut.Check_Output()
    if len(dataOut.FRB_final)<len(dataOut.SIT_final):
        num = len(dataOut.SIT_final)-len(dataOut.FRB_final)
        dataOut.FRB_final = np.append(dataOut.FRB_final, [np.nan for el in range(num)])
        dataOut.FRB_std = np.append(dataOut.FRB_std, [np.nan for el in range(num)])
        dataOut.FRB_ln = np.append(dataOut.FRB_ln, [0 for el in range(num)])
    if len(dataOut.SD_final)<len(dataOut.SIT_final):
        num = len(dataOut.SIT_final)-len(dataOut.SD_final)
        dataOut.SD_final = np.append(dataOut.SD_final, [np.nan for el in range(num)])
        dataOut.SD_std = np.append(dataOut.SD_std, [np.nan for el in range(num)])
        dataOut.SD_ln = np.append(dataOut.SD_ln, [0 for el in range(num)])
    
    Get_unc(dataOut)
    # print data to output file
    dataOut.Print_to_output(ofile, primary='SIT')

Functions.sort_final_data(ofile, saveplot=saveplot, HS='NH', primary='SIT')
