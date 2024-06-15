# -*- coding: utf-8 -*-
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2023-06-27'
#%% Importing libaries and data

# Third-party imports
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import datetime as dt

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
                    if type(self.obsID)==str and  type(self.pp_flag)==str:
                        print('{:^30s} {:^20s} {:^8.3f} {:^8.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^6.0f} {:^7.3f} {:^7.3f} {:^7.3f} {:^7.3f} {:^4.0f} {:^7d} {:^7d}'.format(self.obsID, (self.date_final[ll]), self.lat_final[ll], self.lon_final[ll], self.SD_final[ll], self.SD_std[ll], self.SD_ln[
                              ll], self.SD_unc[ll], self.SIT_final[ll], self.SIT_std[ll], self.SIT_ln[ll], self.SIT_unc[ll], self.FRB_final[ll], self.FRB_std[ll], self.FRB_ln[ll], self.FRB_unc[ll], self.sur_temp_final[ll], self.air_temp_final[ll], self.w_SD_final[ll]/100., self.w_density_final[ll], self.pp_flag, self.unc_flag), file=output)
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
