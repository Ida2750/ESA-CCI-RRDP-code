# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Henriette Skorup'
__contributors__ = 'Ida Olsen'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2022-07-11'


# -- Built-in modules -- #

# -- Third-part modules -- #
import numpy as np
import datetime as dt
import scipy.spatial as ss 
import pyproj
import math
import warnings
# -- Proprietary modules -- #

def ERS_data(ObsData, KDstruct, CSdata, CSinptime, CSx, CSy, output):
    obsDate = ObsData['date']
    obstime = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") 
                     for s in obsDate])
    
    # create a 1D numpy-array of length obslat with elements NaN 
    CSsif = np.zeros(np.shape(ObsData['lat'])) * np.nan

    time_string = []
    # The length of KDstruct equals the number of observations
    for ii in range(len(KDstruct)):
        sif = CSdata['sif'][KDstruct[ii]]
        CStime = CSinptime[KDstruct[ii]]
        time =np.zeros(np.shape(CStime))
    
        for ij in range(len(time)):
            time[ij] = (CStime[ij]-obstime[ii]).days
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if not np.isnan(np.nanmean(sif[np.abs(time)<15])):
                CSsif[ii] = np.nanmean(sif[np.abs(time)<15])
            
                time_string = dt.datetime.strftime(obstime[ii],"%Y-%m-%dT%H:%M:%S")      
                print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],
                      CSsif[ii],file=output)
                    
    output.close()

def collocation_part2(obsfile, CSfile, count, ofile, var, HS='NH'):
    """

    Parameters
    ----------
    obsfile : String
        Name of reference data file.
    CSfile : String
        Name of satellite data file.
    count : Integer
        Counter to know whether to write or append to output file.
    ofile : String
        Name of output file
    var : List of strings
        A list with relevant variables
    HS : String, optional
        Hemisphere, either NH for Northern or SH for Southern. 
        The default is 'NH'.

    Returns
    -------
    None.

    """
    # Defines standad projection for the Lambert Azimuthal Equeal Area projection
    if HS == 'NH':
        Pollat0  = 90.
        Pollon0  = 0.
    elif HS == 'SH':
        Pollat0  = -90.
        Pollon0  = 0.
    ell = 'WGS84'
    PolProj=pyproj.Proj(proj='laea', ellps=ell, datum='WGS84', lat_0=Pollat0, lon_0=Pollon0, units='m')

    ## Loads satellite data
    if 'ENV' in CSfile or 'CS2' in CSfile:
        CSnames = ['time','lat','lon','sit','sitUnc','sif','sifUnc','sd','sdUnc']
        CSdata = np.genfromtxt(CSfile, dtype=None, names=CSnames, usecols=(0,1,2,3,4,5,6,7,8))
        
        # Calculates sea ice draft
        CSdraft = CSdata['sit'] - CSdata['sif']
        CSdraftunc = CSdata['sitUnc'] + CSdata['sifUnc']
        print('draft computed')
        
    elif 'ERS' in CSfile:
        CSnames = ['time','lat','lon','sif','sifUnc']
        dtypes = [float for n in CSnames]
        CSdata = np.genfromtxt(CSfile, dtype=dtypes, names=CSnames,encoding='utf8') 
    print('data has been loaded')
    
    CSinptime = np.empty(np.shape(CSdata['time'][:]), dtype=dt.datetime)
    
    for jj in range(len(CSdata['time'][:])):
        CSinptime[jj] = dt.datetime.fromtimestamp(CSdata['time'][jj])
     
    # Converts from 0 to 360 to -180 to 180 for projection
    CSdata['lon'] = np.array([CSlon-360 if CSlon>180 else CSlon for CSlon in CSdata['lon']]).flatten()
    
    # Makes projection of satellite data
    CSx,CSy = PolProj(CSdata['lon'],CSdata['lat'])
    print('projection succesfull')
    
    # Reads input data from observation file (RRDP), ASCII format
    ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')
    obsDate = ObsData['date']

    obstime = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") 
                     for s in obsDate])
    
    obsx,obsy = PolProj(ObsData['lon'],ObsData['lat'])
    print('projection of observation data succesfull')
    
    if HS == 'SH':
        dist = 50000/2
    else:
        dist = 25000/2;  # search radius
        
    tree = ss.cKDTree(list(zip(CSx, CSy)))
    KDstruct = tree.query_ball_point(list(zip(obsx, obsy)), dist, p=2)
    print('KDstruct made')
    
    # create a 1D numpy-arrays of length obslat with elements NaN 
    ## variables
    CSsit = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsif = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsid = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsd = np.zeros(np.shape(ObsData['lat'])) * np.nan
    ## std
    CSsitstd = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsifstd = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsidstd = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsdstd = np.zeros(np.shape(ObsData['lat'])) * np.nan
    ## ln
    CSsitln = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsifln = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsidln = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsdln = np.zeros(np.shape(ObsData['lat'])) * np.nan
    ## unc
    CSsitunc = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsifunc = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsidunc = np.zeros(np.shape(ObsData['lat'])) * np.nan
    CSsdunc = np.zeros(np.shape(ObsData['lat'])) * np.nan
    
    
    # Make output file
    print(ofile)
    if count == 0 and not 'QL' in CSfile:
        output = open(ofile,'w')    
    else:
        output = open(ofile,'a')
    
    if 'ERS' in CSfile:
        ERS_data(ObsData, KDstruct, CSdata, CSinptime, CSx, CSy, output)
    else: # ENV or CS2
        ##%% Loops through elements of the KDstruct and calculates the mean
        # of elements within +/-15 days of the time of the reference data  
        time_string = []
        # The length of KDstruct equals the number of reference data observations
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
            
            CStime = CSinptime[KDstruct[ii]]
            time =np.zeros(np.shape(CStime))
            for ij in range(len(time)):
                # find time differences
                time[ij] = (CStime[ij]-obstime[ii]).days
            # Get average value of satellite observations within time difference
            CSsit[ii] = np.nanmean(sit[np.abs(time)<15])
            CSsif[ii] = np.nanmean(sif[np.abs(time)<15])
            CSsd[ii]  = np.nanmean(sd[np.abs(time)<15])
            CSsid[ii]  = np.nanmean(sid[np.abs(time)<15])
	        
            # freeboard
            # number of measurements
            CSsifln[ii] = len(sif[np.isnan(sif)==0]) 
            # standard deviation of co-located satellite observations
            CSsifstd[ii] = np.std(sif[np.abs(time)<15])
            # Uncertainty of co-located satellite observations
            CSsifunc[ii] = 1/CSsifln[ii]*np.sqrt(np.nansum(sif_unc[np.abs(time)<15]**2))
            
            # sea ice thickness
            CSsitln[ii] = len(sit[np.isnan(sit)==0]) 
            CSsitstd[ii] = np.std(sit[np.abs(time)<15])
            CSsitunc[ii] = 1/CSsitln[ii]*np.sqrt(np.nansum(sit_unc[np.abs(time)<15]**2))
            
            # snow depth
            CSsdln[ii] = len(sd[np.isnan(sd)==0]) 
            CSsdstd[ii] = np.std(sd[np.abs(time)<15])
            CSsdunc[ii] = 1/CSsdln[ii]*np.sqrt(np.nansum(sd_unc[np.abs(time)<15]**2))
            
            # sea ice draft
            CSsidln[ii] = len(sid[np.isnan(sid)==0]) 
            CSsidstd[ii] = np.std(sid[np.abs(time)<15])
            CSsidunc[ii] = 1/CSsidln[ii]*np.sqrt(np.nansum(sid_unc[np.abs(time)<15]**2))
            
            # formats datetime object into a string
            time_string = dt.datetime.strftime(obstime[ii],"%Y-%m-%dT%H:%M:%S")
            
            # prints co-located values to output file
            if not (math.isnan(CSsit[ii]) and math.isnan(CSsif[ii]) and math.isnan(CSsd[ii])):
                if 'SID' in var:
                    print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],ObsData['SID'][ii],
                          CSsid[ii],ObsData['SIDstd'][ii],ObsData['SIDln'][ii], ObsData['SIDunc'][ii],
                          CSsidstd[ii], CSsidln[ii],CSsidunc[ii], file=output)
                else:
                    print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],ObsData['SD'][ii],
                          ObsData['SIT'][ii],ObsData['FRB'][ii],CSsd[ii],CSsit[ii],CSsif[ii],
                          ObsData['SDstd'][ii],ObsData['SDln'][ii],ObsData['SDunc'][ii],
                          ObsData['SITstd'][ii],ObsData['SITln'][ii],ObsData['SITunc'][ii],
                          ObsData['FRBstd'][ii],ObsData['FRBln'][ii],ObsData['FRBunc'][ii],
                          CSsdstd[ii], CSsdln[ii], CSsdunc[ii], CSsitstd[ii], CSsitln[ii], 
                          CSsitunc[ii], CSsifstd[ii], CSsifln[ii],CSsifunc[ii],file=output)
                    
        output.close()
    
    
