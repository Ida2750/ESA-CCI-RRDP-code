# -*- coding: utf-8 -*-
"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured at the North Pole Enviromental Observatory
"""

# -- File info -- #

__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2024-08-12'

# -- Built-in modules -- #
import os.path
import datetime as dt
import glob
import re
import sys

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
from Warren import SnowDepth, SWE
import Functions


def Get_lat_lon(file):
    """
    Get lat and lon information from file header

    Parameters
    ----------
    file : TYPE
        inputfilename with absolute path

    Returns
    -------
    lat : float
        Mooring latitude
    lon : float
        Mooring longitude
    header : int
        number of headerlines

    """
    lat = 0
    lon = 0
    #to find end of header
    lookup = b'#'
    #to find location of lat and lon
    latlon = b'Position'
    count = 0
    with open(file, 'rb') as myFile:
        for num, line in enumerate(myFile, 1):
            #print(line.strip())
            count += 1
            #either if latlon in line or if the lat/lon is split over two lines
            if latlon in line or (lat != 0 and lon == 0):
                if 'N' in line.decode('latin-1') and 'E' in line.decode('latin-1'):
                    #print(line.decode('latin-1'))
                    #print('-----------')
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = float(latlon_num[3])+lonmin/60
                elif 'N' in line.decode('latin-1') and 'W' in line.decode('latin-1'):
                    ##if coordinates given to the west then the longitude must be -
                    #print(line.decode('latin-1'))
                    latlon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(latlon_num[1]+'.'+latlon_num[2])
                    # convert to decimal degrees
                    lat = float(latlon_num[0])+latmin/60
                    lonmin = float(latlon_num[4]+'.'+latlon_num[5])
                    # convert to decimal degrees
                    lon = -(float(latlon_num[3])+lonmin/60)

                # if lat/lon is split over two lines
                elif 'North' in line.decode('latin-1'):
                    lat_num = re.findall('[0-9]+', line.decode('latin-1'))
                    latmin = float(lat_num[1]+'.'+lat_num[2])
                    # convert to decimal degrees
                    lat = float(lat_num[0])+latmin/60

                elif 'East' in line.decode('latin-1'):
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = float(lon_num[0])+lonmin/60
                elif 'West' in line.decode('latin-1'):
                    print(line.decode('latin-1'))
                    lon_num = re.findall('[0-9]+', line.decode('latin-1'))
                    lonmin = float(lon_num[1]+'.'+lon_num[2])
                    # convert to decimal degrees
                    lon = -float(lon_num[0])+lonmin/60
            #break loop when lookup is found
            if lookup in line:
                header = count  # note number of headerlines
                break

    return lat, lon, header

def robust_std(data):
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median))
    return mad * 1.4826  # scaling factor for normal distribution

dtint = 30  # days
gridres = 25000  # m

# directory = os.path.dirname(os.path.dirname(os.path.dirname(
#     os.getcwd()))) + '/RRDPp/RawData/NPEO/data'
directory = '/dmidata/projects/cmems2/C3S/RRDPp/RawData/NPEO/data'

# saving locations
save_path_data = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPEO/final/'
ofile = save_path_data + 'ESACCIplus-SEAICE-RRDP2-SID-NPEO.nc'
saveplot = os.path.dirname(os.path.dirname(
    os.getcwd())) + '/FINAL/NPEO/fig/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
if not os.path.exists(saveplot):os.makedirs(saveplot)

count = 0
for dir in os.listdir(directory):
    ## the 2001 dir is different from the rest
    if '2001' in dir:
        dir_data = os.path.join(directory, dir+'/ULS')
    else:
        dir_data = os.path.join(directory, dir+'/'+dir+'/ULS')
    #print(dir_data)
    for ifile in os.listdir(dir_data):
        print(ifile)
        #only the files containing min are of interest to us
        if 'min' in ifile:
            file = os.path.join(dir_data, ifile)
            lat, lon, header = Get_lat_lon(file)

            #defines variables in output file
            count += 1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
            dataOut.pp_flag = 0
            dataOut.unc_flag = 2

            # # Reads observation data from ASCII-file
            data = np.genfromtxt(file, dtype=None, skip_header=header, names=[
                                 'YEAR', 'MO', 'DAY', 'HR', 'MIN', 'SEC', 'TIME', 'DEPTH', 'DRAFT'], encoding='latin-1')

            years = data['YEAR']
            months = data['MO']
            days = data['DAY']
            hr = data['HR']
            min = data['MIN']
            sec = data['SEC']
            dates = [dt.datetime(years[i], months[i], days[i],
                                 hr[i], min[i], sec[i]) for i in range(len(years))]
            dataOut.obsID = 'NPEO_' + ifile.split('_')[2]
            SID = data['DRAFT'].astype(float)

            # uncertainty of 5 cm approx.
            SID_Unc = np.array([0.05 for sid in SID])
            # set non existing data to nan
            SID[SID == -999.000] = np.nan
            SID[SID < 0] = np.nan
            SID[SID > 8] = np.nan

            t = np.array([(date-dt.datetime(1970, 1, 1)).total_seconds()
                         for date in dates])

            # calculates average, number of data, std and unc of each month
            mondiff = np.where(~(np.diff(months) == 0))[
                0]  # index of when month changes
            # add start and end indexes on month
            index = np.insert(np.append(mondiff, len(months)-1), 0, 0)

            # dataOut.SID_final = np.array(
            #     [np.nanmean(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            # dataOut.SID_std = np.array(
            #     [np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.SID_final=np.array([np.nanmedian(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            #dataOut.SID_std=np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.SID_std = np.array([robust_std(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
            dataOut.unc_flag = [3 if sid>2 else 2 for sid in dataOut.SID_final]
            dataOut.SID_ln = np.array(
                [len(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])

            for i in range(len(index)-1):
                start = index[i]
                end = index[i+1]
                dataOut.SID_unc = np.append(
                    dataOut.SID_unc, 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2)))

            #find median date within month (tells about which part of the month majority measurements are from)
            avgDates = np.array([np.median(t[index[i]:index[i+1]])
                                for i in range(len(index)-1)])

            ######### QUALITY FLAGS ############
            dataOut.QFT = [] # temporal
            dataOut.QFS = [] # spatial
            dataOut.QFG = [] # global threshold

            days = [time.day for time in dates]
            for i in range(len(index)-1):
                # find number of days
                unique = len(np.unique(days[index[i]:index[i+1]]))
                if np.any(SID[index[i]:index[i+1]]>6):
                    dataOut.QFG.append(1)
                else:
                    dataOut.QFG.append(0)
                if unique==1:
                    dataOut.QFT.append(3)
                    dataOut.QFS.append(3)
                elif unique<=5:
                    dataOut.QFT.append(2)
                    dataOut.QFS.append(3)
                elif unique<15:
                    dataOut.QFT.append(1)
                    dataOut.QFS.append(3)
                elif unique>=15:
                    dataOut.QFT.append(0)
                    dataOut.QFS.append(0)
            
            dataOut.QFT = np.array(dataOut.QFT)
            dataOut.QFS = np.array(dataOut.QFS)
            dataOut.QFG = np.array(dataOut.QFG)
            ######### QUALITY FLAGS ############  

            #change date to right format
            dates_final = np.array(
                [dt.datetime.fromtimestamp(int(sec)) for sec in avgDates])
            dataOut.date_final = np.array([dt.datetime.strftime(
                date, "%Y-%m-%dT%H:%M:%S") for date in dates_final])

            # Correlate NPEO data with Warren snow depth and snow density
            for ll in range(np.size(dataOut.SID_final, 0)):
                (w_SD, w_SD_epsilon) = SnowDepth(
                    lat, lon, dates_final[ll].month)
                dataOut.w_SD_final = np.append(dataOut.w_SD_final, w_SD)
                (wswe, wswe_epsilon) = SWE(
                    lat, lon, dates_final[ll].month)
                w_density = int((wswe/w_SD)*1000)
                dataOut.w_density_final = np.append(
                    dataOut.w_density_final, w_density)

            dataOut.lat_final = [lat for el in dataOut.date_final]
            dataOut.lon_final = [lon for el in dataOut.date_final]

            dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
            dataOut.pp_flag = [dataOut.pp_flag]*len(dataOut.SID_final)
            dataOut.obsID = [dataOut.obsID]*len(dataOut.SID_final)

            #time_in = [dt.datetime(int(y[:4]),int(m),int(d)) for y,m,d in zip(dates, months, days)]
            Functions.scatter(dataOut.obsID[0], dates, SID, dataOut.time, dataOut.SID_final, 'SID [m]', saveplot)
            
            # fill empty arrays with NaN values
            dataOut.Check_Output()

            if count>1:
                subset = dataOut.Create_NC_file(ofile, primary='SID')
                df = Functions.Append_to_NC(df, subset)
            else:
                df = dataOut.Create_NC_file(ofile, primary='SID', datasource='North Pole Environmental Observatory (NPEO) Oceanographic Mooring Data: https://doi.org/10.5065/D6P84921', key_variables='Sea Ice Draft')
                

# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SID')       
