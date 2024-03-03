# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 21:40:19 2022

@author: Ida Olsen (s174020)
"""


import numpy as np
import datetime as dt
import sys
import scipy.spatial as ss 
import pyproj
from pyproj import CRS
import os
import math

# path = "C:/Users/Ida Olsen/Documents/work/RRDPp/satelitte/Collocationv3p3/OIB"
# CSfile = os.path.join(path, "OIB_ENV_CCIp.out")

def collocation_part2(obsfile, CSfile, count, ofile, var, HS='NH'):

    # Define standad projection for the Lambert Azimuthal Equeal Area projection
    if HS == 'NH':
        Pollat0  = 90.
        Pollon0  = 0.
    elif HS == 'SH':
        Pollat0  = -90.
        Pollon0  = 0.
    ell = 'WGS84'
    PolProj=pyproj.Proj(proj='laea', ellps=ell, datum='WGS84', lat_0=Pollat0, lon_0=Pollon0, units='m')
    daylim = dt.timedelta(days=15)

    
    CSnames = ['time','lat','lon','sit','sitUnc','sif','sifUnc','sd','sdUnc']
    CSdata = np.genfromtxt(CSfile, dtype=None, names=CSnames, usecols=(0,1,2,3,4,5,6,7,8))
    print('data has been loaded')
    
    CSinptime = np.empty(np.shape(CSdata['time'][:]), dtype=dt.datetime)
    
    for jj in range(len(CSdata['time'][:])):
        CSinptime[jj] = dt.datetime.fromtimestamp(CSdata['time'][jj])
     
                                           
    #%% calculate draft
    CSdraft = CSdata['sit'] - CSdata['sif']
    CSdraftunc = CSdata['sitUnc'] + CSdata['sifUnc']
    print('draft computed')
    
    ## convert from 0 to 360 to -180 to 180 for projection
    CSdata['lon'] = np.array([CSlon-360 if CSlon>180 else CSlon for CSlon in CSdata['lon']]).flatten()
    CSx,CSy = PolProj(CSdata['lon'],CSdata['lat'])
    print('projection succesfull')
    # Reads input data from observation file (RRDP), ASCII format
    file = open(obsfile, 'r')
    #%%

    ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')
    
    obsID =  ObsData['obsID']
    obsDate = ObsData['date']
    ## lat/lon are 
    obslat = ObsData['lat']
    obslon = ObsData['lon']
    
    if 'SID' in var:    
        # obsSID = ObsData['SID']
        obsSID = [abs(sid) for sid in ObsData['SID']]
        obsSID_std = ObsData['SIDstd']
        obsSID_ln = ObsData['SIDln']
        try:
            obsSID_unc = ObsData['SIDunc']
        except:
            print('Missing unc information')
            obsSID_unc =  np.array([np.nan] *len(obslat)).flatten()
    if 'SIT' in var:
        obsSIT = ObsData['SIT']
        obsSIT_std = ObsData['SITstd']
        obsSIT_ln = ObsData['SITln']
        try:
            obsSIT_unc = ObsData['SITunc']
        except:
            obsSIT_unc = np.array([np.nan] *len(obslat)).flatten()
            print('Missing unc information')
    if 'FRB' in var:
        obsFRB = ObsData['FRB']
        obsFRB_std = ObsData['FRBstd']
        obsFRB_ln = ObsData['FRBln']
        obsFRB_unc = ObsData['FRBunc']
        
    if 'SD' in var:
        obsSD = ObsData['SD']
        obsSD_std = ObsData['SDstd']
        obsSD_ln = ObsData['SDln']
        try:
            obsSD_unc = ObsData['SDunc']
        except:
            obsSD_unc = np.array([np.nan] *len(obslat)).flatten()
            print('Missing unc information')
    
    obstime = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") 
                     for s in obsDate])
    
    obsx,obsy = PolProj(obslon,obslat)
    print('projection of observation data succesfull')
    if HS == 'SH':
        dist = 50000/2
    else:
        dist = 25000/2;  # search radius
    
    tree = ss.cKDTree(list(zip(CSx, CSy)))
    KDstruct = tree.query_ball_point(list(zip(obsx, obsy)), dist, p=2)
    print('KDstruct made')
    
    # create a 1D numpy-array of length obslat with elements NaN 
    CSsit = np.zeros(np.shape(obslat)) * np.nan
    CSsif = np.zeros(np.shape(obslat)) * np.nan
    CSsid = np.zeros(np.shape(obslat)) * np.nan
    CSsd = np.zeros(np.shape(obslat)) * np.nan
    ## std
    CSsitstd = np.zeros(np.shape(obslat)) * np.nan
    CSsifstd = np.zeros(np.shape(obslat)) * np.nan
    CSsidstd = np.zeros(np.shape(obslat)) * np.nan
    CSsdstd = np.zeros(np.shape(obslat)) * np.nan
    ## ln
    CSsitln = np.zeros(np.shape(obslat)) * np.nan
    CSsifln = np.zeros(np.shape(obslat)) * np.nan
    CSsidln = np.zeros(np.shape(obslat)) * np.nan
    CSsdln = np.zeros(np.shape(obslat)) * np.nan
    ## unc
    CSsitunc = np.zeros(np.shape(obslat)) * np.nan
    CSsifunc = np.zeros(np.shape(obslat)) * np.nan
    CSsidunc = np.zeros(np.shape(obslat)) * np.nan
    CSsdunc = np.zeros(np.shape(obslat)) * np.nan
    
    print(ofile)
    print(count)
    if count == 0 and not 'QL' in CSfile:
        output = open(ofile,'w')
        
    else:
        output = open(ofile,'a')
    
    time_string = []
    # The length of KDstruct equals the number of observations
    for ii in range(len(KDstruct)):
        sit = CSdata['sit'][KDstruct[ii]]
        sif = CSdata['sif'][KDstruct[ii]]
        sid = CSdraft[KDstruct[ii]]
        sd = CSdata['sd'][KDstruct[ii]]
        # unc
        sit_unc = CSdata['sitUnc'][KDstruct[ii]]
        sif_unc = CSdata['sifUnc'][KDstruct[ii]]
        sid_unc = CSdraftunc[KDstruct[ii]]
        sd_unc = CSdata['sdUnc'][KDstruct[ii]]
        
        #print np.mean(free)
        CStime = CSinptime[KDstruct[ii]]
        time =np.zeros(np.shape(CStime))
        for ij in range(len(time)):
            time[ij] = (CStime[ij]-obstime[ii]).days
        CSsit[ii] = np.nanmean(sit[np.abs(time)<15])
        CSsif[ii] = np.nanmean(sif[np.abs(time)<15])
        CSsd[ii]  = np.nanmean(sd[np.abs(time)<15])
        CSsid[ii]  = np.nanmean(sid[np.abs(time)<15])
	        
        ## uncertainties + number of measurements
        CSsifln[ii] = len(sif[np.isnan(sif)==0]) 
        CSsifstd[ii] = np.std(sif[np.abs(time)<15])
        CSsifunc[ii] = 1/CSsifln[ii]*np.sqrt(np.nansum(sif_unc[np.abs(time)<15]**2))
        
        CSsitln[ii] = len(sit[np.isnan(sit)==0]) 
        CSsitstd[ii] = np.std(sit[np.abs(time)<15])
        CSsitunc[ii] = 1/CSsitln[ii]*np.sqrt(np.nansum(sit_unc[np.abs(time)<15]**2))
        
        CSsdln[ii] = len(sd[np.isnan(sd)==0]) 
        CSsdstd[ii] = np.std(sd[np.abs(time)<15])
        CSsdunc[ii] = 1/CSsdln[ii]*np.sqrt(np.nansum(sd_unc[np.abs(time)<15]**2))
        
        CSsidln[ii] = len(sid[np.isnan(sid)==0]) 
        CSsidstd[ii] = np.std(sid[np.abs(time)<15])
        CSsidunc[ii] = 1/CSsidln[ii]*np.sqrt(np.nansum(sid_unc[np.abs(time)<15]**2))
    
        time_string = dt.datetime.strftime(obstime[ii],"%Y-%m-%dT%H:%M:%S")
        
        if not (math.isnan(CSsit[ii]) and math.isnan(CSsif[ii])):
            if 'SID' in var:
                print(time_string,obslat[ii],obslon[ii],obsSID[ii],CSsid[ii],obsSID_std[ii],
                      obsSID_ln[ii], obsSID_unc[ii], CSsidstd[ii], CSsidln[ii],CSsidunc[ii], file=output)
                print(time_string,obslat[ii],obslon[ii],obsSID[ii],CSsid[ii],obsSID_std[ii],
                      obsSID_ln[ii], obsSID_unc[ii], CSsidstd[ii], CSsidln[ii], CSsidunc[ii])
            else:
                print(time_string,obslat[ii],obslon[ii],obsSD[ii],obsSIT[ii],obsFRB[ii],
                      CSsd[ii],CSsit[ii],CSsif[ii],obsSD_std[ii],obsSD_ln[ii],obsSD_unc[ii],
                      obsSIT_std[ii], obsSIT_ln[ii],obsSIT_unc[ii],obsFRB_std[ii],obsFRB_ln[ii],
                      obsFRB_unc[ii], CSsdstd[ii], CSsdln[ii], CSsdunc[ii],
                      CSsitstd[ii], CSsitln[ii], CSsitunc[ii], CSsifstd[ii], CSsifln[ii],CSsifunc[ii],file=output)
                
    output.close()
    
    
