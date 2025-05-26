#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Get data from newer AEM files (all files post 2012)
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2022-01-26'

import os
import numpy as np
import datetime as dt
import netCDF4 as nc
from datetime import date, timedelta, datetime

def Get_New_AEM_data(ifile):
    #ifile="/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/AEM_AWI/2012/HEM_PAM12_20120403T143731_20120403T151350.nc"
    
    if ifile.endswith('.nc') or ifile.endswith('.dat')  or ifile.endswith('.tab'):
                #ifile="/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/AEM_AWI/2012/HEM_PAM12_20120403T143731_20120403T151350.nc"
                if ifile.endswith('.nc'):
                    data = nc.Dataset(ifile)
                    year=data['YEAR'][:]
                    month=data['MONTH'][:]
                    day=data['DAY'][:]
                    time=data['TIME'][:]
                    lat=data['LATITUDE'][:]
                    lon=data['LONGITUDE'][:]
                    SIT=data['THICKNESS'][:]
                    obsID = 'AEM_AWI_' + ifile[-40:-26]
                else:
                    print('entered')
                    print(os.path.basename(ifile))
                    if os.path.basename(ifile).startswith('P'): # for newer files with snow data
                        ObsNames=['DateTime','LATITUDE','LONGITUDE','Sea_icesnow thick_m','Sea_icesnow thick unc_m', 'Snow_thick_m',	'Snow thick unc [m]',	'Snow thick std dev [±]'] #,	'Snow thick min [m]'	,'Snow thick max [m]',	'NOBS [#]',	'Snow freeboard [m]',	'Snow freeboard unc [m]',	'Snow freeboard std dev [±]',	'Snow freeboard min [m]',	'Snow freeboard max [m]',	'Surf temp [°C]',	'Surf temp std dev [±]',	'Surf temp min [°C]',	'Surf temp max [°C]',	'EsEs [m]',	'Sea ice freeboard [m]',	'Density ice [kg/m**3]',	'Density ice unc [kg/m**3]',	'Sea ice age [a]',	'Ice type',	'QF1',	'QF2',	'QF3',	'QF4',	'QF5',	'QF6']
                        dtypes = [object, float, float, float, float, float, float, float]
                        data = np.genfromtxt(ifile, names=ObsNames, dtype=dtypes, skip_header=49, usecols=(0,1,2,3,4,5,6,7))
                        date=data['DateTime']
                        date = [d.decode('utf8') for d in date]
                        lat=data['LATITUDE']
                        lon=data['LONGITUDE']
                        SIT=data['Sea_icesnow_thick_m']
                        SIT_unc=data['Sea_icesnow_thick_unc_m']
                        SD=data['Snow_thick_m']
                        SD_unc=data['Snow_thick_unc_m']
                        # surf temp?
                        
                        obsID = 'AEM_AWI_' + ifile[-25:-15]
                    else: # for newer files with no snow data
                        ObsNames=['YEAR','MONTH','DAY','TIME','??','LATITUDE','LONGITUDE','???','THICKNESS','????']
                        data = np.genfromtxt(ifile, names=ObsNames)
                        year=data['YEAR'].astype(int)
                        month=data['MONTH'].astype(int)
                        day=data['DAY'].astype(int)
                        time=data['TIME']
                        lat=data['LATITUDE']
                        lon=data['LONGITUDE']
                        SIT=data['THICKNESS']
                        obsID = 'AEM_AWI_' + ifile[-21:-13]                        
                
                if not os.path.basename(ifile).startswith('P'):
                    # create datetime array in correct format
                    offset = []
                    date = []
        
                    for i in range(len(time)):
                        # seconds = np.round(time[i])  # given in decimal seconds
                        seconds = time[i]
                        start = datetime(year[1], month[1], day[1], 0, 0, 0)
                        
                        delta = timedelta(seconds=np.float64(seconds))
                        offset = start + delta
                        
                        date.append(offset.strftime('%Y-%m-%dT%H:%M:%S'))
                try:
                    t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                                 for s in date])
                except:
                    t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S.%f")
                                 for s in date])
                # sort t as measurements are not in timewise order
                index = t.argsort()
                t.sort()
                
                if os.path.basename(ifile).startswith('P'):
                    return [obsID, lat[index], lon[index], SIT[index], SIT_unc[index], SD[index], SD_unc[index], t]
                else:
                    return [obsID, lat[index], lon[index], SIT[index], t]
