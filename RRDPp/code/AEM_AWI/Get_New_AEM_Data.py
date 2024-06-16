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

import numpy as np
import datetime as dt
import netCDF4 as nc
def Get_New_AEM_data(ifile):
    #ifile="/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/AEM_AWI/2012/HEM_PAM12_20120403T143731_20120403T151350.nc"
    
    if ifile.endswith('.nc') or ifile.endswith('.dat'):
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
                    
                # create datetime array in correct format
                from datetime import date, timedelta, datetime
                offset = []
                date = []
    
                for i in range(len(time)):
                    # seconds = np.round(time[i])  # given in decimal seconds
                    seconds = time[i]
                    start = datetime(year[1], month[1], day[1], 0, 0, 0)
                    
                    delta = timedelta(seconds=np.float64(seconds))
                    offset = start + delta
                    
                    date.append(offset.strftime('%Y-%m-%dT%H:%M:%S'))
                
                t = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S")
                             for s in date])
                # sort t as measurements are not in timewise order
                index = t.argsort()
                t.sort()
                
                
                return [obsID, lat[index], lon[index], SIT[index], t]