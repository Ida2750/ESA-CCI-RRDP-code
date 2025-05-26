# -*- coding: utf-8 -*-
"""
Created on Mon Jun 12 13:24:45 2023

@author: Ida Olsen

Uses EASE-grid to produce 25 km grid mean values. Mean values are obtained either by the distance limit or the time limit of 30 days.
Uses input files measured by Ice mass buoys from the The CRREL- Dartmouth Mass Balance Buoys Program
includes Warren snow depths and densities
"""
# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '1'
__date__ = '2023-06-12'

# -- Built-in modules -- #
import os.path
import sys
import re
import datetime as dt

# -- Third-part modules -- #
import numpy as np
from numpy.lib import recfunctions as rfn
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
#import PyPDF2
import pandas as pd
import csv

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import EASEgrid_correct as EASEgrid
import Functions


#%% Support functions


class IMB:
    def __init__(self, obsID, data):
        self.data = data
        self.obsID = obsID
        self.date = []
        self.lat = []
        self.lon = []
        self.SD = []
        self.SIT = []
        self.FRB = []
        self.unc_flag = 3
        self.pp_flag = []

def write_info_table(self, initial_SD, initial_SIT, files, MassDate, TempDate, MetDate, SD, SIT):
    """
    Parameters
    ----------
    ID
    date
    location?
    number of sd sensors available
    initial sd and initial sit
    
    pre-processing flag
    uncertainty flag
    
    variables in file
    -------
    None.
    """
    #################
    pdf_file = [file for file in files if file.endswith('.pdf')][0]
    pdfFileObj = open(os.path.join(dataPath, pdf_file), 'rb')
      
    # creating a pdf reader object
    pdfReader = PyPDF2.PdfReader(pdfFileObj)
      
    # creating a page object - get first page
    pageObj = pdfReader.pages[0]
    
    # extract text from page
    text = pageObj.extract_text()
    
    # substitue special characters
    text = re.sub('[^A-Za-z,]+', ' ', text)
    
    # search for position
    search1 = re.search('Autonomous Ice Mass Balance Buoy Observations, Buoy', text)
    position = text[search1.span()[1]:search1.span()[1]+35]
    position = position.split(',')[1]
    
    # closing the pdf file object
    pdfFileObj.close()
    #######################
    intervals = []
    for d in [MassDate, TempDate, MetDate]:
        try:
            time = [dt.datetime.strptime(dd, '%Y-%m-%d %H:%M') for dd in d]
        except: 
            time = [dt.datetime.strptime(dd, '%m/%d/%y %H:%M') for dd in d]
        if len(intervals)==0: ## append dates of mass information
            date_start = str(time[0].date())
            date_end = str(time[-1].date())
            print(date_start)
        intervals.append(np.nanmean(np.diff(time)).total_seconds()/3600)
                         
    s_SD = 'O'
    s_SIT = 'O'
    if any(~np.isnan(SD)):
        s_SD = 'X'
    if any(~np.isnan(SIT)):
        s_SIT = 'X'
        
    ofile = 'IMB_info_table_Henriette.dat'
    output = open(ofile, 'a')
    #  & {:^6.2f} & {:^6.2f}

    try:
        print('{:^s} & {:^s} & {:^s} & {:^s} {:^6.2f} & {:^6.2f}'.format(
            self.obsID, position, date_start, date_end, initial_SD, initial_SIT), file=output) 
        # print('{:^s} & {:^s} & {:^s} & {:^s} & {:^6.2f} & {:^6.2f} & {:^6.2f} & {:^s} & {:^s} &  {:^6.2f} & {:^6.2f}'.format(
        #     self.obsID, position, date_start, date_end, intervals[0], intervals[1], intervals[2], s_SD, s_SIT, SD[~np.isnan(SD)][0], SIT[~np.isnan(SIT)][0]),file=output)  #, np.nanmean(self.pp_flag)# np.nanmean(self.unc_flag)), file=output)
    except:
        # print('{:^s} & {:^s} & {:^s} & {:^s} & {:^6.2f} & {:^6.2f} & {:^6.2f} & {:^s} & {:^s} &  {:^6.2f} & {:^6.2f}'.format(
        #     self.obsID, position, date_start, date_end, intervals[0], intervals[1], intervals[2], s_SD, s_SIT, SD[0], SIT[0]),file=output)  #, np.nanmean(self.pp_flag)# np.nanmean(self.unc_flag)), file=output)
        print('{:^s} & {:^s} & {:^s} & {:^s} {:^6.2f} & {:^6.2f}'.format(
            self.obsID, position, date_start, date_end, SD[0], SIT[0]),file=output) 
    output.close()

    return None

def interpolate_positions(self, MassDate, year):
    # pre-processing flag! There are several days betwen day of
    # position measurements and day of massbalance data aquisition
    if year>=2017:
        Massdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M:%S') -
                    dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    else:
        Massdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M') -
                    dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    if year>=2017:
        Posdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M:%S') -
                dt.datetime(1970, 1, 1)).total_seconds() for dat in self.data['Date']]
    else:
        try:
            if len(self.data['Date'][0])>10:
                Posdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M') -
                            dt.datetime(1970, 1, 1)).total_seconds() for dat in self.data['Date']]
            else: # to time provided
                Posdates = [(dt.datetime.strptime(dat, '%Y-%m-%d') -
                            dt.datetime(1970, 1, 1)).total_seconds() for dat in self.data['Date']]        
        except:  # different format
            Posdates = [(dt.datetime.strptime(dat, '%m/%d/%y %H:%M') -
                        dt.datetime(1970, 1, 1)).total_seconds() for dat in self.data['Date']]

    # print(Posdates)
    # convert to total seconds since 1970
    self.deltaT = []
    for date in Massdates:
        ## find location with smallest possible time difference
        index = np.argmin(abs(np.array(date)-np.array(Posdates)))
        # print(np.min(abs(np.array(date)-np.array(Posdates))))
        if np.min(abs(np.array(date)-np.array(Posdates))) > 12*3600:  # more difference than 0.5 day!
            self.pp_flag.append(2)
        else:
            self.pp_flag.append(1)
        self.deltaT.append(np.min(abs(np.array(date)-np.array(Posdates))))
        self.lat.append(self.data['Latitude'][index])
        self.lon.append(self.data['Longitude'][index])


    # # some fauly (nan) latitudes and longitudes exist - these must be removed
    # #self.lat = self.[np.isnan(self.lat)]
    #     #print(self.lat)
    #     if any(np.isnan(self.lat)):

    #         print(self.data['Latitude'])
    #         print(self.data['Latitude'][np.isnan(self.data['Latitude'])])
    #         break
        


def Find_temp_file_overlap(MassDate, TempDate):
    # pre-processing flag! There are several days betwen day of
    # position measurements and day of massbalance data aquisition
    if year>=2017:
        Massdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M:%S') -
                    dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    else:
        Massdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M') -
                    dt.datetime(1970, 1, 1)).total_seconds() for dat in MassDate]
    if year>=2017:
        Tempdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M:%S') -
                dt.datetime(1970, 1, 1)).total_seconds() for dat in TempDate]
    else:
        try:
            Tempdates = [(dt.datetime.strptime(dat, '%Y-%m-%d %H:%M') -
                        dt.datetime(1970, 1, 1)).total_seconds() for dat in TempDate]
        except:  # different format
            Tempdates = [(dt.datetime.strptime(dat, '%m/%d/%y %H:%M') -
                        dt.datetime(1970, 1, 1)).total_seconds() for dat in TempDate]

    xy, x_ind, y_ind = np.intersect1d(
        Massdates, Tempdates, return_indices=True)

    return x_ind, y_ind

def Get_ppflag(self, SD, lnSD, lnSIT):
    """
    Assign pre-processing flag, where
    the data gets pp-flag=2 if any position used in the average
    was obtained with a temporal discrepancy of more than 12 hours
    """
    # find valid SD entries
    pp_non_nan = np.array(self.pp_flag)[~np.isnan(SD)]
    # print(pp_non_nan)
    pp_flag = np.zeros(len(lnSD))
    summ = np.cumsum(lnSD).astype(int)
    try:
        pp_flag[0] = np.nanmax(pp_non_nan[:summ[0]])
        # print(np.max(pp_non_nan[:summ[0]]))
    except: # no data in datapoint
        pp_flag[0] = np.nan
    # find max pre-processing flag within the interval used for the average
    for i in range(1, len(summ)):
        # print(summ[i])
        try:
            pp_flag[i] = np.nanmax(pp_non_nan[summ[i-1]:summ[i]])
        except: # no data in datapoint
            pp_flag[i] = np.nan
    self.pp_flag = pp_flag
    print(pp_flag)
    return self.pp_flag


#%% Main


# Information
save_path = os.path.dirname(os.path.dirname(os.getcwd())) + '/IMB/'
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/IMB/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
ofile = 'ESACCIplus-SEAICE-RRDP2+-SIT-IMB-countlast.nc'
gridres = 25000  # grid resolution
dtint = 30  # days per mean
## saving locations
ofile = os.path.join(save_path_data, ofile)


saveplot = os.path.join(os.path.dirname(
    os.path.dirname(os.getcwd())), 'FINAL/IMB/fig/')
# create directory if they do not exist
if not os.path.exists(saveplot): os.makedirs(saveplot)
#dataDir = os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd()))) + '/RRDPp/RawData/IMB_CRREL/'
dataDir = '/dmidata/projects/cmems2/C3S/RRDPp/RawData/IMB_CRREL/'

count = 0
datalen_SD = 0
datalen_SIT = 0
for dirr in os.listdir(dataDir):
    # Enter Snow files
    if dirr.startswith('Snow'):
            path = os.path.join(dataDir, dirr)
            # Enter identifier (e.g. 2002A)
            for dirrr in [ '2018A','2018B','2018C','2018D','2018E']: # sorted(os.listdir(path)):
                    year = int(re.findall(r'\d+\.?\d*', dirrr)[0])
                    print('---------------------------')
                    print(' Printing data for buoy: ', dirrr)
                    print('---------------------------')
                    
                    #defines variables in output file
                    dataOut = Functions.Final_Data(Type='SIT', count_head=count+1)
        
                    path2 = os.path.join(path, dirrr)
                    dataPath = path2 + '/' + os.listdir(path2)[0] + '/data/'

                    print(dataPath)
                    files = os.listdir(dataPath)
                    files_csv = [f for f in files if f.endswith('.csv')]
                    pos_file = [file for file in files_csv if 'Position' in file][0]
                    MassBalance_file = [
                        file for file in files_csv if 'Mass_Balance' in file][0]
                    Temp_file = [file for file in files_csv if 'Temp' in file][0]
                    Met_file = [file for file in files_csv if 'Meteo' in file][0]
                    try:
                        Meta_file = [file for file in files_csv if 'Metadata' in file][0] 
                    except: # newer data not awailable initial conditions are added
                        pass 

                    ## Read position data
                    file = open(os.path.join(dataPath, pos_file), 'r')
                    # Reads observation data from ASCII-file
                    Pos_data = np.genfromtxt(
                        file, dtype=None, names=True, delimiter=',', encoding=None)

                    if year>=2017: # differnet format of newer files
                        # 'buoy_name', 'Date_mmddyyyy_hhmm', 'Latitude_Degrees', 'Longitude_Degrees' 
                        Pos_data.dtype.names = ('bouy_id', 'Date', 'Latitude', 'Longitude')
                    # remove invalid latitude and longitude entries (-9999)
                    index = [False if la <= 0 or lo < -180 else True for la,
                             lo in zip(Pos_data['Latitude'], Pos_data['Longitude'])]
                    Pos_data = Pos_data[index]
                    # remove nan values
                    index = np.isfinite(Pos_data['Latitude']) & np.isfinite(Pos_data['Longitude'])
                    Pos_data = Pos_data[index]
                    #print(Pos_data)
                    
                    try:
                        # look for initial snowdepth
                        file = open(os.path.join(dataPath, Meta_file), 'r')
                        Meta_data = csv.reader(file,  delimiter=',')
                        for row in Meta_data:
                            
                            if 'Initial Snow Depth' in ', '.join(row):
                                try:
                                    initial_SD = float(re.findall(r'\d+\.\d+|\d+', ', '.join(row))[0])
                                except: # if not defined add 0
                                    initial_SD = np.nan
                            if 'Initial Ice' in ', '.join(row):
                                try:
                                    initial_SIT = float(re.findall(r'\d+\.\d+|\d+', ', '.join(row))[0])
                                except: # not defined
                                    initial_SIT = np.nan
                    except: 
                        initial_SD = np.nan
                        initial_SIT = np.nan
        
                    ## Read Mass Balance data
                    file = open(os.path.join(dataPath, MassBalance_file), 'r')
                    Mass_data = np.genfromtxt(
                        file, dtype=None, names=True, delimiter=',', encoding=None)

                    ## Read Ice Temperature data
                    file = open(os.path.join(dataPath, Temp_file), 'r')
                    Temp_data = np.genfromtxt(
                        file, dtype=None, names=True, delimiter=',', encoding=None)
        
                    ## Read Air Temperature data
                    file = open(os.path.join(dataPath, Met_file), 'r')
                    Met_data = np.genfromtxt(
                        file, dtype=None, names=True, delimiter=',', encoding=None)
                    if year>=2017: # different format of newer data
                        Mass_data.dtype.names = ('bouy_id', 'Date', 'TOIP', 'BOIP', 'Snow_Depth', 'Ice_Thickness', 'offset')
                        Met_data.dtype.names = ('bouy_id', 'Date', 'Air_Temp', 'Air_Pressure')
                        Temp_data = rfn.rename_fields(Temp_data, {'Date_mmddyyyy_hhmm': 'Date'})
                        Temp_data = rfn.rename_fields(Temp_data, {'Temp_0cm': '0'})
                    # Interpolate positions from position file to the mass data
                    IMB1 = IMB(dirrr, Pos_data)
                    if Mass_data.size > 1:
                        count+=1
                       
                        interpolate_positions(IMB1, Mass_data['Date'], year)

                        #if np.isfinite(initial_SD):
                        SD = Mass_data['Snow_Depth'].astype(float) # + initial_SD  # meters
                        #else:
                        #    SD = Mass_data['Snow_Depth'].astype(float)
                        SIT = Mass_data['Ice_Thickness'].astype(
                            float)  # m calculated by TOP-BOP
            
                        ## Temperature and mass balance data are recorded with time
                        # intervals - signifying that if dates are correct they can be
                        # used simultaniously

                        Mass_index, Temp_index = Find_temp_file_overlap(
                            Mass_data['Date'], Temp_data['Date'])
                           
                        nan_elements = np.delete([i for i in range(len(SD))], Mass_index)
                        sur_temp = Temp_data['0'][Temp_index].astype(
                            float)  # degrees celcius at ice surface

                        for el in nan_elements:
                            sur_temp = np.insert(sur_temp, el, np.nan)
            
                        Mass_index, Met_index = Find_temp_file_overlap(
                            Mass_data['Date'], Met_data['Date'])
                        nan_elements = np.delete([i for i in range(len(SD))], Mass_index)
                        air_temp = Met_data['Air_Temp'][Met_index].astype(
                            float)  # degrees celcius
                        for el in nan_elements:
                            air_temp = np.insert(air_temp, el, np.nan)
            
                        #set '-9999' values to np.nan
                        SD[SD < 0] = np.nan
                        SD[SD > 2] = np.nan
                        SIT[SIT < 0] = np.nan
                        SIT[SIT > 10] = np.nan
                        air_temp[air_temp < -900] = np.nan
                        sur_temp[sur_temp < -900] = np.nan
            

                        datalen_SD += len(SD[np.isfinite(SD)])
                        datalen_SIT += len(SIT[np.isfinite(SIT)])

                        # Changes date format into date time format
                        if year>=2017:
                            t = np.array([dt.datetime.strptime(s, "%Y-%m-%d %H:%M:%S")
                              for s in Mass_data['Date']])
                        else:
                            t = np.array([dt.datetime.strptime(s, "%Y-%m-%d %H:%M")
                              for s in Mass_data['Date']])
            
                        # uncertainty estimate
                        unc = 0.01  # m
                        SD_unc = np.array([unc for i in range(len(SD))])
                        SIT_unc = np.array([unc for i in range(len(SIT))])
                        dataOut.unc_flag = 3 
                        
                        # Create EASEgrid and returns grid cell indicies (index_i, index_j) for each observation
                        G = EASEgrid.Gridded()
                        G.SetHemisphere('N')
                        G.CreateGrids(gridres)
                        (index_i, index_j) = G.LatLonToIdx(IMB1.lat, IMB1.lon)
                        
                        IMB1.lat = np.array(IMB1.lat)
                        #print(IMB1.lon)
                        # Takes the time for each grid cell into account and calculate averages
                        avgSD, stdSD, lnSD, uncSD, lat, lon, time, avgSIT, stdSIT, lnSIT, uncSIT, avgFRB, stdFRB, lnFRB, FRB_Unc, avgTair, avgTsurf, dataOut.QFT, dataOut.QFS, dataOut.QFG = G.GridData(
                            dtint, IMB1.lat, IMB1.lon, t, SD=SD, SD_unc=SD_unc, SIT=SIT, SIT_unc=SIT_unc, VAR1=air_temp, VAR2=sur_temp, dtype='buoy')
                        if len(time) > 0:
                            dataOut.obsID = 'IMB' + dirrr
                            Functions.plot(IMB1.lat, IMB1.lon, dataOut.obsID, time,saveplot, HS='NH')
                            Functions.scatter(dataOut.obsID, t, SD, time, avgSD, 'SD [m]', saveplot)
                            Functions.scatter(dataOut.obsID, t, SIT, time, avgSIT, 'SIT [m]',saveplot)
                        ## Assign pp-flags
                        dataOut.pp_flag = Get_ppflag(IMB1, SD, lnSD, lnSIT)
                        
                        if year<=2016: # there are no pdf files in newer data
                            write_info_table(dataOut, initial_SD, initial_SIT, files, Mass_data['Date'],Temp_data['Date'],Met_data['Date'], SD, SIT)
                        
                        # ## assign uncertainties
                        # SD_unc, SIT_unc = Get_uncertainties(lnSD, lnSIT)
            
                        # Correlates IMB buoy data with Warren snow depth and snow density
                        for ll in range(np.size(avgSD,0)):
                            (w_SD,w_SD_epsilon) = SnowDepth(lat[ll],lon[ll],time[ll].month)
                            dataOut.w_SD_final = np.append(dataOut.w_SD_final,w_SD)
                            (wswe,wswe_epsilon) = SWE(lat[ll],lon[ll],time[ll].month)
                            w_density=int((wswe/w_SD)*1000)
                            dataOut.w_density_final = np.append(dataOut.w_density_final,w_density)
            
                        #Change names to correct format names
                        dataOut.lat_final = lat
                        dataOut.lon_final = lon
                        for ll in range(np.size(time,0)):
                            dataOut.date_final = np.append(dataOut.date_final,dt.datetime.strftime(time[ll],"%Y-%m-%dT%H:%M:%S"))
                        dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
                        dataOut.SD_final = avgSD
                        dataOut.SD_std = stdSD
                        dataOut.SD_ln = lnSD
                        dataOut.SD_unc = uncSD
                        dataOut.SIT_final = avgSIT
                        dataOut.SIT_std = stdSIT
                        dataOut.SIT_ln = lnSIT
                        dataOut.SIT_unc = uncSIT
                        dataOut.air_temp_final = avgTair
                        dataOut.sur_temp_final = avgTsurf
                        dataOut.unc_flag = [dataOut.unc_flag]*len(lat)
                        dataOut.obsID = [dataOut.obsID]*len(lat)
                        
                        # fill empty arrays with NaN values
                        dataOut.Check_Output()
                        
                        if count>1:
                            subset = dataOut.Create_NC_file(ofile, primary='SIT')
                            df = Functions.Append_to_NC(df, subset)
                        else:
                            df = dataOut.Create_NC_file(ofile, primary='SIT', datasource=' Monitoring the mass balance, motion, and thickness of Arctic sea ice, http://imb-crrel-dartmouth.org', key_variables='Sea ice thickness and snow depth')

print(f'sd len {datalen_SD}')
print(f'sit len {datalen_SIT}')
# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SIT')
