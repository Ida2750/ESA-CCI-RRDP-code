# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Henriette Skorup'
__contributors__ = 'Ida Olsen'
__contact__ = ['ilo@dmi.dk']
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

def robust_std(data):
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median))
    return mad * 1.4826  # scaling factor for normal distribution

def ERS_data(ObsData, KDstruct, CSdata, CSinptime, CSx, CSy, output, days_half):

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
            if not np.isnan(np.nanmean(sif[np.abs(time)<=days_half])):
                CSsif[ii] = np.nanmean(sif[np.abs(time)<=days_half])
            
                time_string = dt.datetime.strftime(obstime[ii],"%Y-%m-%dT%H:%M:%S")      
                print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],
                      CSsif[ii],file=output)
                    
    output.close()

def collocation_part2(obsfile, CSfile, count, ofile, var, HS='NH', days=30, resolution=25000):
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

    ccc = 0

    # Defines standad projection for the Lambert Azimuthal Equeal Area projection
    if HS == 'NH':
        Pollat0  = 90.
        Pollon0  = 0.
        dist = resolution/2
    elif HS == 'SH':
        Pollat0  = -90.
        Pollon0  = 0.
        dist = (resolution*2)/2
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
        
        #apply limits
        CSdraft[CSdraft>8] = np.nan
        CSdata['sit'][CSdata['sit']>10] = np.nan
        CSdata['sif'][CSdata['sif']>2] = np.nan
        #CSdata['sif'][(CSdata['sif'] + CSdata['sd'])>2] = np.nan
        CSdata['sd'][CSdata['sd']>2] = np.nan

        CSdraft[CSdraft<0] = np.nan
        CSdata['sit'][CSdata['sit']<0] = np.nan
        CSdata['sif'][CSdata['sif']<0] = np.nan
        #CSdata['sif'][(CSdata['sif'] + CSdata['sd'])<0] = np.nan
        CSdata['sd'][CSdata['sd']<0] = np.nan

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
    try:
        ObsData = np.genfromtxt(obsfile,dtype=None,names=True, encoding='Latin1')
    except:
        ObsData = np.genfromtxt(obsfile.replace('SIT', 'SD'),dtype=None,names=True, encoding='Latin1')
    obsDate = ObsData['date']

    obstime = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") 
                     for s in obsDate])
    
    obsx,obsy = PolProj(ObsData['lon'],ObsData['lat'])
    print('projection of observation data succesfull')
    

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
    #print(ofile)
    # temporal resolution
    days_half = days/2

    if count == 0 and not 'QL' in CSfile:
        output = open(ofile,'w')    
    else:
        output = open(ofile,'a')
    
    if 'ERS' in CSfile:
        ERS_data(ObsData, KDstruct, CSdata, CSinptime, CSx, CSy, output, days_half)
    else: # ENV or CS2
        ##%% Loops through elements of the KDstruct and calculates the mean
        # of elements within +/-15 days of the time of the reference data  
        time_string = []
        # The length of KDstruct equals the number of reference data observations
        print(len(ObsData))
        print(len(KDstruct))
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
            CSsit[ii] = np.nanmedian(sit[np.abs(time)<=days_half])
            CSsif[ii] = np.nanmedian(sif[np.abs(time)<=days_half])
            CSsd[ii]  = np.nanmedian(sd[np.abs(time)<=days_half])
            CSsid[ii]  = np.nanmedian(sid[np.abs(time)<=days_half])
	        
            # freeboard
            # number of measurements
            CSsifln[ii] = len(sif[np.abs(time)<=days_half][np.isnan(sif[np.abs(time)<=days_half])==0]) 
            # standard deviation of co-located satellite observations
            CSsifstd[ii] = robust_std(sif[np.abs(time)<=days_half])
            # Uncertainty of co-located satellite observations
            CSsifunc[ii] = 1/CSsifln[ii]*np.sqrt(np.nansum(sif_unc[np.abs(time)<=days_half]**2))
            
            # sea ice thickness
            CSsitln[ii] =  len(sit[np.abs(time)<=days_half][np.isnan(sit[np.abs(time)<=days_half])==0]) 
            CSsitstd[ii] = robust_std(sit[np.abs(time)<=days_half])
            CSsitunc[ii] = 1/CSsitln[ii]*np.sqrt(np.nansum(sit_unc[np.abs(time)<=days_half]**2))
            
            # snow depth
            CSsdln[ii] =  len(sd[np.abs(time)<=days_half][np.isnan(sid[np.abs(time)<=days_half])==0]) 
            CSsdstd[ii] = robust_std(sd[np.abs(time)<=days_half])
            CSsdunc[ii] = 1/CSsdln[ii]*np.sqrt(np.nansum(sd_unc[np.abs(time)<=days_half]**2))
            
            # sea ice draft
            CSsidln[ii] =  len(sid[np.abs(time)<=days_half][np.isnan(sid[np.abs(time)<=days_half])==0]) 
            CSsidstd[ii] = robust_std(sid[np.abs(time)<=days_half])
            CSsidunc[ii] = 1/CSsidln[ii]*np.sqrt(np.nansum(sid_unc[np.abs(time)<=days_half]**2))
            
            # formats datetime object into a string
            time_string = dt.datetime.strftime(obstime[ii],"%Y-%m-%dT%H:%M:%S")
            
            # prints co-located values to output file
            if not (math.isnan(CSsit[ii]) and math.isnan(CSsif[ii]) and math.isnan(CSsd[ii])):
                if 'SID' in var:
                    print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],ObsData['SID'][ii],
                          CSsid[ii],ObsData['SIDstd'][ii],ObsData['SIDln'][ii], ObsData['SIDunc'][ii],
                          ObsData['QFT'][ii],ObsData['QFS'][ii],ObsData['QFG'][ii],
                          CSsidstd[ii], CSsidln[ii],CSsidunc[ii], file=output)
                else:
                    print(time_string,ObsData['lat'][ii],ObsData['lon'][ii],ObsData['SD'][ii],
                          ObsData['SIT'][ii],ObsData['FRB'][ii],CSsd[ii],CSsit[ii],CSsif[ii],
                          ObsData['SDstd'][ii],ObsData['SDln'][ii],ObsData['SDunc'][ii],
                          ObsData['SITstd'][ii],ObsData['SITln'][ii],ObsData['SITunc'][ii],
                          ObsData['FRBstd'][ii],ObsData['FRBln'][ii],ObsData['FRBunc'][ii],
                          ObsData['QFT'][ii],ObsData['QFS'][ii],ObsData['QFG'][ii],
                          CSsdstd[ii], CSsdln[ii], CSsdunc[ii], CSsitstd[ii], CSsitln[ii], 
                          CSsitunc[ii], CSsifstd[ii], CSsifln[ii],CSsifunc[ii], ii, file=output)
                    
                    lat = CSdata["lat"][KDstruct[ii]]
                    lon = CSdata["lon"][KDstruct[ii]]
                    if ObsData['obsID'][ii]=='OIB_QL_20190422':
                        ccc +=1
                        # Make sure all array contents are fully printed
                        np.set_printoptions(threshold=np.inf)
                        with open(f'output_4_{ccc}.txt', 'w') as f100:
                            f100.write(f'obsline {ObsData[ii]}\n')
                            f100.write(f'satlat: {lat[np.abs(time)<=2]}\n')
                            f100.write(f'satlon: {lon[np.abs(time)<=2]}\n')
                            f100.write(f'sattime: {CStime[np.abs(time)<=2]}\n')
                        with open(f'output_16_{ccc}.txt', 'w') as f100:
                            f100.write(f'obsline {ObsData[ii]}\n')
                            f100.write(f'satlat: {lat[np.abs(time)<=8]}\n')
                            f100.write(f'satlon: {lon[np.abs(time)<=8]}\n')
                            f100.write(f'sattime: {CStime[np.abs(time)<=8]}\n')
                        with open(f'output_30_{ccc}.txt', 'w') as f100:
                            f100.write(f'obsline {ObsData[ii]}\n')
                            f100.write(f'satlat: {lat[np.abs(time)<=15]}\n')
                            f100.write(f'satlon: {lon[np.abs(time)<=15]}\n')
                            f100.write(f'sattime: {CStime[np.abs(time)<=15]}\n')
        output.close()
    
    
