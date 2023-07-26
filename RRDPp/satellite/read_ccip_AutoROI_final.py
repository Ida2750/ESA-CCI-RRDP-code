#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Note: In upcoming versions, the orbit data will be provided in
      self describing netcdf data format

Adapted for use of ESA CCI phase 2, SIT RRDP, Henriette Skourup, DTU Space, February 2016 
"""

# -- File info -- #
__author__ = 'Henriette Skorup'
__contributors__ = 'Ida Olsen'
__contact__ = ['hsk@dtu.dk']
__version__ = '0'
__date__ = '2016-02-01'


# -- Built-in modules -- #
import os.path
# -- Third-part modules -- #
from netCDF4 import Dataset
import numpy as np
import datetime as dt
# -- Proprietary modules -- #


class collocation:
  def __init__(self,campaign, obsfile, ifile, ofile, count, i=0):
    """

    Parameters
    ----------
    campaign : string
        Name of Campaign
    obsfile : string
        Filename of reference data file (observations)
    ifile : string
        Filename of input file from satellitte
    ofile : string
        Filename of output file
    count : integer
        Used to identify the file number (for count=0 a new file is created)
    i : integer
        Used to access the correct data in the ERS files

    Returns
    -------
    None.

    """
    self.save_path = os.path.dirname(os.getcwd()) + '/satellite/Satellite_subsets/' + campaign
    self.campaign = campaign
    self.CCIpFile = Dataset(ifile,'r')
    self.lat = self.CCIpFile.variables['lat'][:].flatten()
    lon = self.CCIpFile.variables['lon'][:].flatten()
    self.lon = np.array([l+360 if l<0 else l for l in lon])
    # time unit: seconds since 1970-01-01 for ENV and CS2; days since 1950-01-01 for ERS1/2
    self.time = self.CCIpFile.variables['time'][:] 
    # Reads input data from observation file, ASCII format
    ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')  
    obslat = ObsData['lat']
    obslon = ObsData['lon']
    
    ## Convert lon to 0 to 360
    obslon = np.array([lon+360 if lon<0 else lon for lon in obslon]).flatten()
    ## set longitude limits
    self.lonmin = np.amin(obslon)-5.0 
    self.lonmax = np.amax(obslon)+5.0
    ## set latitude limits
    self.latmin = np.amin(obslat)-0.5
    self.latmax = np.amax(obslat)+0.5
    
    # Path to where data is saved
    CompleteName = os.path.join( self.save_path, ofile)
    ##create file
    if count==0:
        self.output=open(CompleteName,'w')
    # append to file
    else:
        self.output=open(CompleteName,'a')
    
    if 'ERS' in ifile:
        collocation_ERS(self, i)
    else:
        collocation_ENV_CS2(self)
        

def collocation_ERS(self, i):
        sif = self.CCIpFile.variables['radar_freeboard_corr'][:][i].flatten()
        try:
            sifUnc = self.CCIpFile.variables['radar_freeboard_corr_unc'][:][i].flatten()
        except:
            sifUnc = [np.nan for el in sif]
       
        t = dt.date(1950, 1, 1) + dt.timedelta(days = int(self.time[i]))
        # convert time into seconds since 1970
        time = (t - dt.date(1970, 1, 1)).total_seconds()
        # print information to file
        for ii in range(np.size(self.lon)):
            if self.lat[ii]>self.latmin and self.lat[ii]<self.latmax and self.lon[ii]>self.lonmin and self.lon[ii]<self.lonmax: 
                    print(time,self.lat[ii],self.lon[ii],sif[ii],sifUnc[ii], file=self.output)   
               
        self.output.close()   
            
def collocation_ENV_CS2(self):

    sit = self.CCIpFile.variables['sea_ice_thickness'][:]
    sitUnc = self.CCIpFile.variables['sea_ice_thickness_uncertainty'][:]
    sif = self.CCIpFile.variables['sea_ice_freeboard'][:]
    sifUnc = self.CCIpFile.variables['sea_ice_freeboard_uncertainty'][:]
    SnowDepth = self.CCIpFile.variables['snow_depth'][:]
    SnowDepthUnc = self.CCIpFile.variables['snow_depth_uncertainty'][:]
    RegionCode = self.CCIpFile.variables['region_code'][:]
    SIDensity = self.CCIpFile.variables['sea_ice_density'][:]
    SIType  = self.CCIpFile.variables['sea_ice_type'][:]
    SnowDensity = self.CCIpFile.variables['snow_density'][:]

    for ii in range(np.size(self.time)):
        if self.lat[ii]>self.latmin and self.lat[ii]<self.latmax and self.lon[ii]>self.lonmin and self.lon[ii]<self.lonmax: 
                print(self.time[ii],self.lat[ii],self.lon[ii],sit[ii],sitUnc[ii],sif[ii],sifUnc[ii],SnowDepth[ii],SnowDepthUnc[ii],
                    RegionCode[ii],SIDensity[ii],SIType[ii],SnowDensity[ii], file=self.output)       
            
    self.output.close()   
