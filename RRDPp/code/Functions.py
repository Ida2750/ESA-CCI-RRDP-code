# -*- coding: utf-8 -*-
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2023-06-27'
#%% Importing libaries and data
import os

# Third-party imports
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import datetime as dt
import xarray as xr

#%% Functions

# =============================================================================
# #%% FINAL DATA CLASS
# =============================================================================


class Final_Data:
    def __init__(self, Type='SIT', count_head=0):

        self.obsID = []
        self.date_final = []
        self.lat_final = []  # decimal degrees
        self.lon_final = []  # decimal degrees
        self.SD_final = []  # snow depth [m]
        self.SD_std = []  # standard deviation of snow depth
        self.SD_ln = []   # number of SD measurements in each gridcell
        self.SD_unc = []  # SD uncertainty
        self.SID_final = []  # snow depth [m]
        self.SID_std = []  # standard deviation of snow depth
        self.SID_ln = []   # number of SID measurements in each gridcell
        self.SID_unc = []  # SID uncertainty
        self.SIT_final = []  # Sea Ice Thickness [m]
        self.SIT_std = []  # sea Ice thickness standard deviation
        self.SIT_ln = []
        self.SIT_unc = []
        self.FRB_final = []  # freeboard [m]
        self.FRB_std = []  # freeboard standard deviation [m]
        self.FRB_ln = []  # number of FRB measurements in each gridcell
        self.FRB_unc = []  # FRB uncertainty
        self.sur_temp_final = []  # surface temperature [celcius]
        self.air_temp_final = []  # air temperature [celcius]
        self.w_SD_final = []  # snow depth from climatology (Warren)
        self.w_density_final = []  # snow density from climatology (Warren)
        self.pp_flag = 2  # pre-processing flag
        self.unc_flag = 1  # Uncertainty flag

        # header count
        self.count_head = count_head

    def Check_Output(self, *args, **kwargs):

        # Add nan values to nonexisting data
        if len(self.SD_final) == 0:
            self.SD_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SD_std) == 0:
            self.SD_std = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SD_unc) == 0:
            self.SD_unc = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SD_ln) == 0:
            self.SD_ln = np.zeros(np.shape(self.lat_final))
        if len(self.SID_final) == 0:
            self.SID_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SID_std) == 0:
            self.SID_std = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SID_unc) == 0:
            self.SID_unc = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SID_ln) == 0:
            self.SID_ln = np.zeros(np.shape(self.lat_final))
        if len(self.SIT_final) == 0:
            self.SIT_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SIT_std) == 0:
            self.SIT_std = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.SIT_ln) == 0:
            self.SIT_ln = np.zeros(np.shape(self.lat_final))
        if len(self.SIT_unc) == 0:
            self.SIT_unc = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.FRB_final) == 0:
            self.FRB_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.FRB_std) == 0:
            self.FRB_std = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.FRB_ln) == 0:
            self.FRB_ln = np.zeros(np.shape(self.lat_final))
        if len(self.FRB_unc) == 0:
            self.FRB_unc = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.w_density_final) == 0:
            self.w_density_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.w_SD_final) == 0:
            self.w_SD_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.sur_temp_final) == 0:
            self.sur_temp_final = np.zeros(np.shape(self.lat_final)) * np.nan
        if len(self.air_temp_final) == 0:
            self.air_temp_final = np.zeros(np.shape(self.lat_final)) * np.nan
        
        # if any Warren SD are negative set them to NaN
        self.w_SD_final[self.w_SD_final<0]=np.nan
        # if any Warren rho are negative or larger than 900 set them to NaN
        self.w_density_final[self.w_density_final<0]=np.nan
        self.w_density_final[self.w_density_final<0]=np.nan
        

    def Print_to_output(self, ofile, primary='SIT'):
        # print header on output file
        if self.count_head == 1 and primary=='SID': 
            output = open(ofile, 'w')
            print('{:^25s} {:^20s} {:^8s} {:^8s} {:^7s} {:^9s} {:^8s} {:^7s} {:^7s} {:^6s} {:^8s} {:^8s}'.format(
                'obsID', 'date', 'lat', 'lon', 'SID', 'SIDstd', 'SIDln', 'SIDunc', 'wSD', 'w-rho', 'pp-flag', 'unc-flag'), file=output)
            output.close()

        elif self.count_head == 1: # SIT, SD, FRB file
            output = open(ofile, 'w')
            print('{:^30s} {:^20s} {:^8s} {:^8s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^7s} {:^7s} {:^7s} {:^6s} {:^6s} {:^7s} {:^7s} {:^6s} {:^4s} {:^7s} {:^7s}'.format(
                'obsID', 'date', 'lat', 'lon', 'SD', 'SDstd', 'SDln', 'SDunc', 'SIT', 'SITstd', 'SITln', 'SITunc', 'FRB', 'FRBstd', 'FRBln', 'FRBunc', 'Tsur', 'Tair', 'wSD', 'w-rho', 'pp-flag', 'unc-flag'), file=output)
            output.close()
        
        # Append processed data to file
        output = open(ofile, 'a')
        for ll in range(np.size(self.date_final, 0)):
            if ((primary == 'SIT') | (primary == 'SD') | (primary == 'FRB')):
                if ((self.SIT_ln[ll] != 0) | (self.SD_ln[ll] != 0) | (self.FRB_ln[ll] != 0)):  # if we have non nan data for SIT, SD or FRB
                    if type(self.obsID)==str: # and  type(self.pp_flag)==str:
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7d} {:^7d}'.format(self.obsID, (self.date_final[ll]), self.lat_final[ll], self.lon_final[ll], self.SD_final[ll], self.SD_std[ll], self.SD_ln[
                              ll], self.SD_unc[ll], self.SIT_final[ll], self.SIT_std[ll], self.SIT_ln[ll], self.SIT_unc[ll], self.FRB_final[ll], self.FRB_std[ll], self.FRB_ln[ll], self.FRB_unc[ll], self.sur_temp_final[ll], self.air_temp_final[ll], self.w_SD_final[ll]/100., self.w_density_final[ll], self.pp_flag, self.unc_flag[ll]), file=output)
                        # print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7d} {:^7d}'.format(self.obsID, (self.date_final[ll]), self.lat_final[ll], self.lon_final[ll], self.SD_final[ll], self.SD_std[ll], self.SD_ln[
                        #       ll], self.SD_unc[ll], self.SIT_final[ll], self.SIT_std[ll], self.SIT_ln[ll], self.SIT_unc[ll], self.FRB_final[ll], self.FRB_std[ll], self.FRB_ln[ll], self.FRB_unc[ll], self.sur_temp_final[ll], self.air_temp_final[ll], self.w_SD_final[ll]/100., self.w_density_final[ll], self.pp_flag, self.unc_flag))
                    elif type(self.pp_flag)==np.ndarray:
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7.0f} {:^7d}'.format(self.obsID, (self.date_final[ll]), self.lat_final[ll], self.lon_final[ll], self.SD_final[ll], self.SD_std[ll], self.SD_ln[
                              ll], self.SD_unc[ll], self.SIT_final[ll], self.SIT_std[ll], self.SIT_ln[ll], self.SIT_unc[ll], self.FRB_final[ll], self.FRB_std[ll], self.FRB_ln[ll], self.FRB_unc[ll], self.sur_temp_final[ll], self.air_temp_final[ll], self.w_SD_final[ll]/100., self.w_density_final[ll], self.pp_flag[ll], self.unc_flag), file=output)
                    else:
                        print(self.obsID[ll])
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7d} {:^7d}'.format(self.obsID[ll], (self.date_final[ll]), self.lat_final[ll], self.lon_final[ll], self.SD_final[ll], self.SD_std[ll], self.SD_ln[
                              ll], self.SD_unc[ll], self.SIT_final[ll], self.SIT_std[ll], self.SIT_ln[ll], self.SIT_unc[ll], self.FRB_final[ll], self.FRB_std[ll], self.FRB_ln[ll], self.FRB_unc[ll], self.sur_temp_final[ll], self.air_temp_final[ll], self.w_SD_final[ll]/100., self.w_density_final[ll], self.pp_flag, self.unc_flag), file=output)
            elif primary == 'SID':
                if self.SID_ln[ll] != 0:  # if we have non nan data for SID
                    try:
                        print('{:^25s} {:^20s} {:^7.3f} {:^8.3f} {:^7.3f} {:^9.3f} {:^8.0f} {:^7.3f} {:^7.3f} {:^6.0f} {:^6d} {:^6d}'.format(self.obsID, self.date_final[ll], self.lat_final[ll], self.lon_final[ll], self.SID_final[ll], self.SID_std[ll], self.SID_ln[ll], self.SID_unc[ll], self.w_SD_final[ll]/100.,self.w_density_final[ll], self.pp_flag, self.unc_flag),file=output)
                    except: # varying flags
                        print('{:^25s} {:^20s} {:^7.3f} {:^8.3f} {:^7.3f} {:^9.3f} {:^8.0f} {:^7.3f} {:^7.3f} {:^6.0f} {:^6d} {:^6d}'.format(self.obsID, self.date_final[ll], self.lat_final[ll], self.lon_final[ll], self.SID_final[ll], self.SID_std[ll], self.SID_ln[ll], self.SID_unc[ll], self.w_SD_final[ll]/100.,self.w_density_final[ll], self.pp_flag, self.unc_flag[ll]),file=output)

            else:
                print('WRONG VARIABLE NAME')

        output.close()
        
    
    def Create_NC_file(self, primary='SIT', datasource=''):
        """
        Create netCDF file with processed information from all files belonging to the datasource
        
        Parameters
        ----------
        primary : string, optional
            Primary variable identifier. The default is 'SIT'.
        datasource : string, optional
            Name of the datasource e.g. Norwegian Polar Institute SID Fram Srait. The default is ''.

        Returns
        -------
        ds : xarray dataset.

        """
             
        if primary=='SID':
            # Create dataset
            ds = xr.Dataset(
                # variables
                data_vars = dict(
                    obsID = (['time'], self.obsID),
                    date = (['time'], self.date_final),
                    lat = (['time'], self.lat_final),
                    lon = (['time'], self.lon_final),
                    SID = (['time'], self.SID_final),
                    SIDstd = (['time'], self.SID_std),
                    SIDln = (['time'], self.SID_ln),
                    SIDunc = (['time'], self.SID_unc),
                    wSD = (['time'], self.w_SD_final),
                    wrho = (['time'], self.w_density_final),
                    ppflag = (['time'], self.pp_flag),
                    uncflag = (['time'], self.unc_flag),
                    ),
                # coordinates
                coords ={
                    "time": self.time,
                    #"lat": self.lat_final,
                    #"lon": self.lon_final,
                },
                attrs = dict(title = 'Sea ice thickness reference measurements (ESA CCI SIT RRDP)',
                             institution = 'Technical University of Denmark (DTU)',
                             grid = 'Equal-Area Scalable Earth Grid in version 2 (EASE2) from the National Snow and Ice Data Center (NSIDC)',
                             Datasource = datasource,
                             source = os.path.basename(__file__),
                             history = 'created: ' + dt.datetime.today().strftime('%Y-%m-%dT%H:%M')
                             )
            )
        else:
            # Create dataset
            ds = xr.Dataset(
                # variables
                data_vars = dict(
                    obsID = (['time'], self.obsID),
                    date = (['time'], self.date_final),
                    lat = (['time'], self.lat_final),
                    lon = (['time'], self.lon_final),
                    
                    SIT = (['time'], self.SIT_final),
                    SITstd = (['time'], self.SIT_std),
                    SITln = (['time'], self.SIT_ln),
                    SITunc = (['time'], self.SIT_unc),

                    SD = (['time'], self.SD_final),
                    SDstd = (['time'], self.SD_std),
                    SDln = (['time'], self.SD_ln),
                    SDunc = (['time'], self.SD_unc),

                    FRB = (['time'], self.FRB_final),
                    FRBstd = (['time'], self.FRB_std),
                    FRBln = (['time'], self.FRB_ln),
                    FRBunc = (['time'], self.FRB_unc),
                    
                    wSD = (['time'], self.w_SD_final),
                    wrho = (['time'], self.w_density_final),
                    ppflag = (['time'], self.pp_flag),
                    uncflag = (['time'], self.unc_flag),
                    ),
                # coordinates
                coords ={
                    "date": self.date_final,

                },
                # attributes
                attrs = dict(title = 'Sea ice thickness reference measurements (ESA CCI SIT RRDP)',
                             institution = 'Technical University of Denmark (DTU)',
                             grid = 'Equal-Area Scalable Earth Grid in version 2 (EASE2) from the National Snow and Ice Data Center (NSIDC)',
                             Datasource = datasource,
                             source = os.path.basename(__file__),
                             history = 'created: ' + dt.datetime.today().strftime('%Y-%m-%dT%H:%M')
                             )
            )
        return ds


def save_NC_file(ds, ofile, primary):
    """
    Allocates attributes of xarray dataset variables and saves dataset as netCDF file

    Parameters
    ----------
    ds : xarray dataset object
        Dataset containing the data  to be saved
    ofile : string
        Name of output file
    primary : string
        Name of primary varibale either SID, SD, SIT or FRB

    Returns
    -------
    None.

    """
    
    if primary=='SID':
        # Variable attributes to SID data
       
        ds['SID'].attrs['standard_name'] = 'Sea Ice Draft'
        ds['SID'].attrs['long_name'] = 'Gridded, monthly mean sea ice draft'
        #ds['SID'].attrs['_FillValue'] = '-999'
        #ds['SID'].attrs['missing_value'] = '-999'
        ds['SID'].attrs['Unit'] = 'meters (m)'
    
       
        ds['SIDstd'].attrs['standard_name'] = 'Sea Ice Draft standard deviation'
        ds['SIDstd'].attrs['long_name'] = 'standard deviation of monthly mean sea ice draft'
        #ds['SIDstd'].attrs['_FillValue'] = '-999'
        #ds['SIDstd'].attrs['missing_value'] = '-999'
        ds['SIDstd'].attrs['Unit'] = 'meters (m)'
    
       
        ds['SIDln'].attrs['standard_name'] = 'Number of Sea Ice Draft OBS'
        ds['SIDln'].attrs['long_name'] = 'Number of observations included in Sea Ice Draft estimate'
        #ds['SIDln'].attrs['_FillValue'] = '-999'
        #ds['SIDln'].attrs['missing_value'] = '-999'
    
       
        ds['SIDunc'].attrs['standard_name'] = 'Sea Ice Draft uncertainty'
        ds['SIDunc'].attrs['long_name'] = 'Approximate uncertainty of Sea Ice Draft estimate'
        #ds['SIDunc'].attrs['_FillValue'] = '-999'
        #ds['SIDunc'].attrs['missing_value'] = '-999'
        ds['SIDunc'].attrs['Unit'] = 'meters (m)'
    
    else:
        # variable attributes to SIT, FRB and SD data 
        
        ## Sea Ice thicknessv variable attributes
        ds['SIT'].attrs['standard_name'] = 'Sea Ice Thickness'
        ds['SIT'].attrs['long_name'] = 'Gridded, monthly mean Sea Ice Thickness'
        #ds['SIT'].attrs['_FillValue'] = '-999'
        #ds['SIT'].attrs['missing_value'] = '-999'
        ds['SIT'].attrs['Unit'] = 'meters (m)'
    
       
        ds['SITstd'].attrs['standard_name'] = 'Sea Ice Thickness standard deviation'
        ds['SITstd'].attrs['long_name'] = 'standard deviation of monthly mean Sea Ice Thickness'
        #ds['SITstd'].attrs['_FillValue'] = '-999'
        #ds['SITstd'].attrs['missing_value'] = '-999'
        ds['SITstd'].attrs['Unit'] = 'meters (m)'
    
       
        ds['SITln'].attrs['standard_name'] = 'Number of Sea Ice Thickness OBS'
        ds['SITln'].attrs['long_name'] = 'Number of observations included in Sea Ice Thickness estimate'
        #ds['SITln'].attrs['_FillValue'] = '-999'
        #ds['SITln'].attrs['missing_value'] = '-999'
    
       
        ds['SITunc'].attrs['standard_name'] = 'Sea Ice Thickness uncertainty'
        ds['SITunc'].attrs['long_name'] = 'Approximate uncertainty of Sea Ice Thickness estimate'
        #ds['SITunc'].attrs['_FillValue'] = '-999'
        #ds['SITunc'].attrs['missing_value'] = '-999'
        ds['SITunc'].attrs['Unit'] = 'meters (m)'
        

        ## Snow depth variable attributes
        ds['SD'].attrs['standard_name'] = 'Snow Depth'
        ds['SD'].attrs['long_name'] = 'Gridded, monthly mean Snow Depth'
        #ds['SD'].attrs['_FillValue'] = '-999'
        #ds['SD'].attrs['missing_value'] = '-999'
        ds['SD'].attrs['Unit'] = 'meters (m)'
    
        ds['SDstd'].attrs['standard_name'] = 'Snow Depth standard deviation'
        ds['SDstd'].attrs['long_name'] = 'standard deviation of monthly mean Snow Depth'
        #ds['SDstd'].attrs['_FillValue'] = '-999'
        #ds['SDstd'].attrs['missing_value'] = '-999'
        ds['SDstd'].attrs['Unit'] = 'meters (m)'
    
        ds['SDln'].attrs['standard_name'] = 'Number of Snow Depth OBS'
        ds['SDln'].attrs['long_name'] = 'Number of observations included in Snow Depth estimate'
        #ds['SDln'].attrs['_FillValue'] = '-999'
        #ds['SDln'].attrs['missing_value'] = '-999'
    
        ds['SDunc'].attrs['standard_name'] = 'Snow Depth uncertainty'
        ds['SDunc'].attrs['long_name'] = 'Approximate uncertainty of Snow Depth estimate'
        #ds['SDunc'].attrs['_FillValue'] = '-999'
        #ds['SDunc'].attrs['missing_value'] = '-999'
        ds['SDunc'].attrs['Unit'] = 'meters (m)'

        ## Freeboard variable attributes
        ds['FRB'].attrs['standard_name'] = 'Freeboard'
        ds['FRB'].attrs['long_name'] = 'Gridded, monthly mean Freeboard (FRB)'
        ds['FRB'].attrs['Additonal_INFO'] =  ' Total FRB e.g. snow depth and sea ice thickness above local sea level'
        #ds['FRB'].attrs['_FillValue'] = '-999'
        #ds['FRB'].attrs['missing_value'] = '-999'
        ds['FRB'].attrs['Unit'] = 'meters (m)'
    
        ds['FRBstd'].attrs['standard_name'] = 'Freeboard standard deviation'
        ds['FRBstd'].attrs['long_name'] = 'standard deviation of monthly mean Freeboard'
        #ds['FRBstd'].attrs['_FillValue'] = '-999'
        #ds['FRBstd'].attrs['missing_value'] = '-999'
        ds['FRBstd'].attrs['Unit'] = 'meters (m)'
    
        ds['FRBln'].attrs['standard_name'] = 'Number of Freeboard OBS'
        ds['FRBln'].attrs['long_name'] = 'Number of observations included in Freeboard estimate'
        #ds['FRBln'].attrs['_FillValue'] = '-999'
        #ds['FRBln'].attrs['missing_value'] = '-999'
    
        ds['FRBunc'].attrs['standard_name'] = 'Freeboard uncertainty'
        ds['FRBunc'].attrs['long_name'] = 'Approximate uncertainty of Freeboard estimate'
        #ds['FRBunc'].attrs['_FillValue'] = '-999'
        #ds['FRBunc'].attrs['missing_value'] = '-999'
        ds['FRBunc'].attrs['Unit'] = 'meters (m)'
        


    ds['date'].attrs['standard_name'] = 'Date'
    #ds['date'].attrs['long_name'] = 'Gridded, monthly mean sea ice draft'
    #ds['SID'].attrs['_FillValue'] = '-999'
    #ds['SID'].attrs['missing_value'] = '-999'
    ds['date'].attrs['Unit'] = 'format: %Y%m%dT%H%M%s'
        
   
    ds['lat'].attrs['standard_name'] = 'lat'
    ds['lat'].attrs['long_name'] = 'lattitude'
    #ds['SID'].attrs['_FillValue'] = '-999'
    #ds['SID'].attrs['missing_value'] = '-999'
    ds['lat'].attrs['Unit'] = 'degrees north'

   
    ds['lon'].attrs['standard_name'] = 'lon'
    ds['lon'].attrs['long_name'] = 'longitude'
    #ds['SID'].attrs['_FillValue'] = '-999'
    #ds['SID'].attrs['missing_value'] = '-999'
    ds['lon'].attrs['Unit'] = 'degrees east'

   
    ds['wSD'].attrs['standard_name'] = 'Warren (1999) Snow Depth'
    ds['wSD'].attrs['long_name'] = 'Snow Depth for NH from Warren climatology (1999) doi:https://doi.org/10.1175/1520-0442(1999)012%3C1814:SDOASI%3E2.0.CO;2'
    #ds['wSD'].attrs['_FillValue'] = '-999'
    #ds['wSD'].attrs['missing_value'] = '-999'
    ds['wSD'].attrs['Unit'] = 'meters (m)'

   
    ds['wrho'].attrs['standard_name'] = 'Warren (1999) Snow Density'
    ds['wrho'].attrs['long_name'] = 'Snow Density for NH from Warren climatology (1999) doi:https://doi.org/10.1175/1520-0442(1999)012%3C1814:SDOASI%3E2.0.CO;2'
    #ds['wrho'].attrs['_FillValue'] = '-999'
    #ds['wrho'].attrs['missing_value'] = '-999'
    ds['wrho'].attrs['Unit'] = 'kilograms per meters cubed (kg/m3)'
    
   
    ds['ppflag'].attrs['standard_name'] = 'Pre-proceessing flag'
    ds['ppflag'].attrs['long_name'] = 'Pre-proceessing flag, indicating the amount of pre-processing done to data'
    ds['ppflag'].attrs['Values'] = '0: No pre-processing, 1: Very minor pre-processing, 2: Minor pre-processing, 3: Major pre-processing'
    #ds['ppflag'].attrs['_FillValue'] = '-999'
    #ds['ppflag'].attrs['missing_value'] = '-999'

   
    ds['uncflag'].attrs['standard_name'] = 'Uncertainty flag'
    ds['uncflag'].attrs['long_name'] = 'Uncertainty flag, indicating the confidence in the uncertainty estimates with 0 being highest confidence and 3 being the lowest'
    ds['uncflag'].attrs['Values'] = '0: Individual uncertainties, 1: Some degree of distinction in uncertainties, 2: Same uncertainty for all data, 3: Same uncertainty and incorrect assumptions'
    #ds['uncflag'].attrs['_FillValue'] = '-999'
    #ds['uncflag'].attrs['missing_value'] = '-999'
    
   
    # Save NC file
    for var in ds:
        ds[var].encoding.update(dict(zlib=True, complevel=6))
    ds.to_netcdf(ofile, format="NETCDF4", mode="w")
    ds.close()
    
    print(f'File {ofile} created')
    
def Append_to_NC(ds, subset):
    """
    Concatenate xarray datasets along the time dimension

    Parameters
    ----------
    ds : xarray dataset
        Primary dataset to combine data to.
    subset : xarray dataset
        Subset of data from file to add to the combined dataset

    Returns
    -------
    ds :  xarray dataset
        Combined dataset

    """
    
    ds = xr.concat([ds,subset], dim="time")   
    return ds
# =============================================================================
# # %% Plotting functions
# =============================================================================


def plot(latitude, longitude, obsID, date, saveplot, HS='NH', obstype=''):
    """
    creates geographical plots of data after processing
    
    Parameters
    ----------
    latitude : numpy array
    longitude : numpy array
    obsID : string
        Observation Identifier
    date : list of datetime object
        list of date range of observations
    Returns
    -------
    None.
    """

    plt.figure(figsize=(6, 6))
    if HS == 'SH':
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -50, -90], ccrs.PlateCarree())
    else:
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    ax.gridlines()
    try:
        plt.title(obstype + obsID + ': ' +
                  str(date[0].date()) + ' - ' + str(date[-1].date()))
    except:
        plt.title(obstype + str(date[0][:10]) + ' - ' + str(date[-1][:10]))
    plt.scatter(longitude, latitude,
                s=10, c='k', alpha=0.9,
                transform=ccrs.PlateCarree(),)
    plt.savefig(saveplot + '/' + obsID + '.png')
    plt.show()

    return None


def scatter(obsID, date_pre, var_pre, date_post, var_post, varname, saveplot):
    """
    creates scatterplot of data before and after processing

    Parameters
    ----------
    obsID : string
        Observation Identifier.
    date_pre : list of datetime objects
        Dates of observations before processing
    var_pre : numpy array
        values of relevant variable (SID; SIT..) before processing
    date_post :  list of datetime objects
        Dates of observations after processing.
    var_post : numpy array
        values of relevant variable (SID; SIT..) after processing
    varname : string
        Name of variable e.g. 'SID', 'SIT', 'SD' or 'FRB'
    saveplot : string
        Folder to where the plots should be saved

    Returns
    -------
    None.

    """

    plt.figure(figsize=(6, 6))
    plt.scatter(date_pre, var_pre, s=10, label='Original data')
    plt.scatter(date_post, var_post, s=10, label='Processed data')
    plt.xlabel('Time')
    plt.ylabel(varname)
    if 'SD' in varname:
        plt.ylim(-0.1, 2)
    elif 'SID' in varname:
        plt.ylim(-1, 50)
    else:
        plt.ylim(-1, 10)
    plt.legend()
    plt.grid()
    plt.xticks(rotation=45)
    plt.title(obsID)
    plt.savefig(saveplot + '/' + obsID + varname +
                '_scatter.png', bbox_inches='tight')
    plt.show()

    return None

# =============================================================================
# # Computing SIT, SD and UNC for ASSIST and ASPeCt
# =============================================================================


def compute_SD_SIT(conc_tot, cc_P, SIT_P, SD_P, cc_S, SIT_S, SD_S, cc_T, SIT_T, SD_T):
    """
    Computes the combined SD and SIT for ASSIST and ASPeCt data based on 
    partial concentrations and thicknesses.

    Parameters
    ----------
    conc_tot : numpy array of floats
        Total combined concentration
    cc_P : numpy array of floats
        Partial concentration of primary SIT
    SIT_P : numpy array of floats
        Partial thickness of primary SIT.
    SD_P : numpy array of floats
        Partial thickness of primary SD.
    cc_S : numpy array of floats
        Partial concentration of seconday SIT.
    SIT_S : numpy array of floats
        Partial thickness of seconday SIT.
    SD_S : numpy array of floats
        Partial thickness of seconday SIT.
    cc_T : numpy array of floats
        Partial concentration of tertiary SIT.
    SIT_T : numpy array of floats
        Partial thickness of tertiary SIT.
    SD_T : numpy array of floats
        Partial thickness of tertiary SIT.

    Returns
    -------
    SD : numpy array of floats
        combined snow depth
    SIT : numpy array of floats
        combined sea ice thickness

    """
    SD = []
    SIT = []
    for kk in range(len(conc_tot)):
        # enter if concentration is more than 0 and is defined (non nan)
        if conc_tot[kk] > 0 and ~np.isnan(conc_tot[kk]):
            # compute the amount of the total ice which is primary, seconday and tertiary
            SIT_eff_P = np.multiply(
                np.divide(cc_P[kk], conc_tot[kk]), SIT_P[kk])
            SIT_eff_S = np.multiply(
                np.divide(cc_S[kk], conc_tot[kk]), SIT_S[kk])
            SIT_eff_T = np.multiply(
                np.divide(cc_T[kk], conc_tot[kk]), SIT_T[kk])
            # compute the amount of the total snow which is primary, seconday and tertiary
            SD_eff_P = np.multiply(np.divide(cc_P[kk], conc_tot[kk]), SD_P[kk])
            SD_eff_S = np.multiply(np.divide(cc_S[kk], conc_tot[kk]), SD_S[kk])
            SD_eff_T = np.multiply(np.divide(cc_T[kk], conc_tot[kk]), SD_T[kk])

            print(SIT_P)
            # Append non nan value if:
            # 1. The total concentration is the same as the sum of the partial concentrations
            # 2. As a minimum one of the SIT/SD_eff are non nan
            try:
                assert ~np.isnan(cc_P[kk]) or ~np.isnan(
                    cc_S[kk]) or ~np.isnan(cc_T[kk])
                assert np.nansum(
                    [cc_P[kk], cc_S[kk], cc_T[kk]]) == conc_tot[kk]
                assert ~np.isnan(SIT_eff_P) or ~np.isnan(
                    SIT_eff_S) or ~np.isnan(SIT_eff_T)
                SIT = np.append(SIT, np.nansum(
                    [SIT_eff_P, SIT_eff_S, SIT_eff_T]))
            except AssertionError:
                SIT = np.append(SIT, np.nan)
            try:
                assert ~np.isnan(cc_P[kk]) or ~np.isnan(
                    cc_S[kk]) or ~np.isnan(cc_T[kk])
                assert np.nansum(
                    [cc_P[kk], cc_S[kk], cc_T[kk]]) == conc_tot[kk]
                assert ~np.isnan(SD_eff_P) or ~np.isnan(
                    SD_eff_S) or ~np.isnan(SD_eff_T)
                SD = np.append(SD, np.nansum([SD_eff_P, SD_eff_S, SD_eff_T]))
            except AssertionError:
                SD = np.append(SD, np.nan)

        elif conc_tot[kk] == 0:  # Everything is water
            SIT = np.append(SIT, 0)
            SD = np.append(SD, 0)
        else: # if the total concentration is not defined
            SIT = np.append(SIT, np.nan)
            SD = np.append(SD, np.nan)
    return SD, SIT


def Get_unc(SD, SIT):
    """
    Append uncertainty to ASSIST and ASPeCt data according to the sea ice
    thickness. See the belonging publication for the origin of the uncertainty
    estimates

    Parameters
    ----------
    SD : numpy array
        Observed snow depths
    SIT : numpy array
        Observaed sea ice thicknesses.

    Returns
    -------
    SD_unc : list
        uncertainty of individual SD measurements
    SIT_unc : list
        uncertainty of individual SIT measurements.

    """
    SD_unc = []
    SIT_unc = []
    for sit, sd in zip(SIT, SD):
        if np.isnan(sit):
            SIT_unc.append(np.nan)
        elif sit <= 0.10:  # m
            SIT_unc.append(sit*0.5)
        elif sit <= 0.3:  # m
            SIT_unc.append(sit*0.3)
        elif sit > 0.3:
            SIT_unc.append(sit*0.2)
        if np.isnan(sd):
            SD_unc.append(np.nan)
        elif sd <= 0.10:  # m
            SD_unc.append(sd*0.5)
        elif sd <= 0.3:  # m
            SD_unc.append(sd*0.3)
        elif sd > 0.3:
            SD_unc.append(sd*0.2)
    return SD_unc, SIT_unc

# =============================================================================
# # Sorting function - to have final data sorted by date
# =============================================================================


def sort_final_data(ofile, saveplot, HS='NH', primary='SIT'):
    """
    Sort final files according to date
        
    Parameters
    ----------
    ofile : string
        complete path to input file to be sorted according to date
    saveplot : string
        Complete path of saving location of final files
    HS : string, optional
        Hemisphere identifier, either NH or SH. The default is 'NH'.
    primary : string, optional
        identifier of primariy variable (SD, FRB, SIT or SID). The default is 'SIT'.

    Returns
    -------
    None.

    """
    if primary == 'SIT':
        inputfile = open(ofile, 'r')
        dtype = [object, object, float, float,
                 float, float, float, float,
                 float, float, float, float,
                 float, float, float, float,
                 float, float, float, float,
                 int, int]
        data = np.genfromtxt(inputfile, dtype=dtype, names=True, encoding=None)
        dates = [dt.datetime.strptime(dat.decode(
            "utf-8"), '%Y-%m-%dT%H:%M:%S') for dat in data['date']]
        sorting = np.argsort(dates)

        data = data[sorting]

        print('Average SD: ', np.nanmean(data['SD']))
        print('Average STD SD: ', np.nanmean(data['SDstd']))
        print('Average Unc SD: ', np.nanmean(data['SDunc']))

        print('Average SIT: ', np.nanmean(data['SIT']))
        print('Average STD SIT: ', np.nanmean(data['SITstd']))
        print('Average Unc SIT: ', np.nanmean(data['SITunc']))

        print('Average FRB: ', np.nanmean(data['FRB']))
        print('Average STD FRB: ', np.nanmean(data['FRBstd']))
        print('Average Unc FRB: ', np.nanmean(data['FRBunc']))

        # print('Average SID: ', np.nanmean(data['SID']))
        # print('Average STD SID: ', np.nanmean(data['SIDstd']))
        # print('Average Unc SID: ', np.nanmean(data['SIDunc']))

        plot(data['lat'], data['lon'], 'All', dates, saveplot, HS=HS)

        for var in ['SD', 'SIT', 'FRB']:
            plt.figure(figsize=(6, 6))
            if var == 'SD':
                plt.xlim(-0.1, 1)
                bins = np.linspace(-2, 2, 200)
            else:
                plt.xlim(-0.1, 6)
                bins = np.linspace(-2, 10, 200)
            plt.hist(data[var], bins=bins)
            plt.grid()

            plt.xlabel(var + ' [m]')
            plt.ylabel('count')
            plt.show()

        inputfile.close()

        ofile = ofile[:-4] + 'sorted.dat'
        output = open(ofile, 'w')
        print('{:^20s} {:^20s} {:^8s} {:^8s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^7s} {:^6s} {:^6s} {:^6s}'.format(
            'obsID', 'date', 'lat', 'lon', 'SD', 'SDstd', 'SDln', 'SDunc', 'SIT', 'SITstd', 'SITln', 'SITunc', 'FRB', 'FRBstd', 'FRBln', 'FRBunc', 'Tsur', 'Tair', 'wSD', 'w-rho', 'pp-flag', 'unc-flag'), file=output)

        for ll in range(len(data['lat'])):
            print('{:^20s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^7.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.0f} {:^7d} {:^7d}'.format(
                data['obsID'][ll].decode('utf8'),
                data['date'][ll].decode('utf8'),
                data['lat'][ll],
                data['lon'][ll],
                data['SD'][ll],
                data['SDstd'][ll],
                data['SDln'][ll],
                data['SDunc'][ll],
                data['SIT'][ll],
                data['SITstd'][ll],
                data['SITln'][ll],
                data['SITunc'][ll],
                data['FRB'][ll],
                data['FRBstd'][ll],
                data['FRBln'][ll],
                data['FRBunc'][ll],
                data['Tsur'][ll],
                data['Tair'][ll],
                data['wSD'][ll],
                data['wrho'][ll],
                data['ppflag'][ll],
                data['uncflag'][ll]
            ), file=output)
        output.close()

    elif primary == 'SID':

        inputfile = open(ofile, 'r')
        dtype = [object, object, float, float, float,
                 float, float, float, float, float, int, int]
        
        data = np.genfromtxt(inputfile, dtype=dtype, names=True)
        dates = [dt.datetime.strptime(dat.decode(
            "utf-8"), '%Y-%m-%dT%H:%M:%S') for dat in data['date']]
        sorting = np.argsort(dates)

        data = data[sorting]

        print('Average SID: ', np.mean(data['SID']))
        print('Average STD SID: ', np.mean(data['SIDstd']))
        print('Average Unc SID: ', np.mean(data['SIDunc']))

        plot(data['lat'], data['lon'], 'All', dates, saveplot, HS=HS)

        for var in ['SID']:
            plt.figure(figsize=(6, 6))
            plt.xlim(-0.1, 6)
            bins = np.linspace(-2, 10, 50)
            plt.hist(data[var], bins=bins)
            plt.grid()

            plt.xlabel(var + ' [m]')
            plt.ylabel('count')
            plt.show()

        inputfile.close()

        ofile = ofile[:-4] + 'sorted.dat'
        output = open(ofile, 'w')
        print('{:^25s} {:^20s} {:^8s} {:^8s} {:^7s} {:^9s} {:^8s} {:^7s} {:^7s} {:^6s} {:^8s} {:^8s}'.format(
            'obsID', 'date', 'lat', 'lon', 'SID', 'SIDstd', 'SIDln', 'SIDunc', 'wSD', 'w-rho', 'pp-flag', 'unc-flag'), file=output)

        for ll in range(len(data['lat'])):
            print('{:^25s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^9.3f} {:^8.0f} {:^7.3f} {:^7.3f} {:^6.0f} {:^8d} {:^8d}'.format(
                data['obsID'][ll].decode('utf8'),
                data['date'][ll].decode('utf8'),
                data['lat'][ll],
                data['lon'][ll],
                data['SID'][ll],
                data['SIDstd'][ll],
                data['SIDln'][ll],
                data['SIDunc'][ll],
                data['wSD'][ll],
                data['wrho'][ll],
                data['ppflag'][ll],
                data['uncflag'][ll]),
                file=output)
