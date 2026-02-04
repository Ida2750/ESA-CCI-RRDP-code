# -*- coding: utf-8 -*-
"""
Calculates 25 km gridded, monthly means of input files containing 
measurements of SID from submarines 

ANALOG files are recorded based on roll numbers, e.g. the analog roll number 
and the segment of the roll used. Therefore these should be combined
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
import re
import sys

# -- Third-part modules -- #
import numpy as np
from haversine import haversine, inverse_haversine, Unit
import math
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import EASEgrid_correct as EASEgrid
import Functions
from Warren import SnowDepth, SWE

#%% Functions
def drft_files(scicex_class, dir_data, ifile):
    file=os.path.join(dir_data, ifile)
    
    ### OBTAIN LAT; LON AND DATE INFORMATION FROM FILE HEADER
    file_end = -1
    with open(file, 'r') as the_file:
        c = 0
        marker = 0
        third_of_month = []
        year = []
        month =  []
        day =  []
        latBeg = []
        latEnd = []
        lonBeg = []
        lonEnd = []
        for line in the_file.readlines() : 
            #print(line)
            c+=1
            if re.search('Beginning Latitude', line) or re.search('beginning latitude', line):
                latBeg = float(re.findall('\d+\.\d+', line)[0]) # starting lattitude
            if re.search('Ending Latitude', line) or re.search('ending latitude', line):
                latEnd = float(re.findall('\d+\.\d+', line)[0]) # ending lattitude
            if re.search('Beginning Longitude', line) or re.search('beginning longitude', line):
                neg = re.findall('-', line)
                lonBeg = float(re.findall('\d+\.\d+', line)[0]) # starting longitude
                if '-' in neg:
                    lonBeg = -lonBeg
            if re.search('Ending Longitude', line) or re.search('ending longitude', line):
                neg = re.findall('-', line)
                lonEnd = float(re.findall('\d+\.\d+', line)[0]) # ending longitude
                if '-' in neg:
                    lonEnd = -lonEnd
            ## DATE INFORMATION
            if re.search('Year', line):
                year = int(re.findall('[0-9]+',line)[0]) # starting lattitude
            if re.search('Month', line) and not re.search('Third of Month', line):
                month = re.findall(r':(\w+)', line.strip().replace(' ', ''))[0] # ending lattitude
                month = mtn(month)
            if re.search('Third of Month', line):
                third_of_month = re.findall('[0-9]+',line)[0] # starting longitude
                day = third_of_mtn(third_of_month)
                s.third_of_month = day
            elif re.search('Day', line):
                day = int(re.findall('[0-9]+',line)[0]) # starting longitude
            if re.search('START DATE/TIME', line) or re.search('start date/time', line):
                # 70476   80000
                try:
                    dates = re.findall('[0-9]+',line)[0]
                    dateStart = dt.datetime.strptime(dates,"%d%m%Y%H%M%S")
                except:
                    dates = re.findall('[0-9]+',line)
                    if len(dates[0])<6:
                        dates[0] = '0' + dates[0]
                    if len(dates[1])<6:
                        dates[1] = '0' + dates[1]
                    dates = dates[0]+dates[1]
                    dateStart = dt.datetime.strptime(dates,"%d%m%y%H%M%S")
            if re.search('END DATE/TIME', line) or re.search('end date/time', line):
                # 70476  103822
                try:
                    dates = re.findall('[0-9]+',line)[0]
                    dateEnd = dt.datetime.strptime(dates,"%d%m%Y%H%M%S")
                except:
                    dates = re.findall('[0-9]+',line)
                    if len(dates[0])<6:
                        dates[0] = '0' + dates[0]
                    if len(dates[1])<6:
                        dates[1] = '0' + dates[1]
                    dates = dates[0]+dates[1]
                    dateEnd = dt.datetime.strptime(dates,"%d%m%y%H%M%S")
                    print(dateEnd)
            try:
                if float(line.strip()) and marker==0:
                    marker=1
                    #re.search('Begin', line): ### DOES NOT CONTAIN BEGIN
                    header_end = c-1
            except:
                pass
            if re.search('!!!End', line.strip().replace(' ', '')):
                file_end = c-1
            elif 'End of Data' in line or 'end of data' in line:
                file_end = c-1
        #### WRITE SEGMENT DESCRIPTION TO FILE
        if re.search('SEGMENT DESCRIPTION', line):
            ## write segment description to file along with filename and date
            d=1
            ## if third of month
            pp_flag_explanation = ' this flag is used as the date is approximated by the third'
            
    
    
        #### WRITE PP-flag + comment in description file
        s.pp_flag = 3 # reason approximate dates and locations (linear interpolation)

        
        if file_end==-1:
            file_end=c
        data = np.genfromtxt(file,dtype=None,skip_header=header_end, skip_footer=c-file_end, delimiter=',')
        try: # format of data strange with "" used to mark the data rows
            #dist = [float(d[0].decode('utf8').replace('"','')) for d in data]
            SID = [float(d.decode('utf8').replace('"','')) for d in data] # sea ice draft
        except:
            #dist = [d[0] for d in data]
            SID = [d for d in data]
        start = (latBeg, lonBeg)
        end = (latEnd, lonEnd)  
        # print(np.array(SID) + 0.29)
    
            
        ## method for calculating lat, lon
        """ Initital steps
            1. make a linearly spaced latitude array between start and end
            2. 
        """

        # lets say we wanted to figure out a series of points between the two given with linear interpolation
        if latEnd!=latBeg:
            latitude = np.linspace(latBeg, latEnd, len(SID))  # ten points
            longitude = (lonBeg - lonEnd)/(latBeg - latEnd)*(latitude - latEnd) + lonEnd
        else:
            latEnd = latEnd+0.001
            latitude = np.linspace(latBeg, latEnd, len(SID))  # ten points
            longitude = (lonBeg - lonEnd)/(latBeg - latEnd)*(latitude - latEnd) + lonEnd
        # if lon<-180:
        #     lon = lon+360
        assert np.round(latitude[-1],1),np.round(longitude[-1],1) == (latEnd, lonEnd)
        print('starting lat, lon', latitude[0],longitude[0],latBeg, lonBeg)
        print('ending lat, lon', latitude[-1],longitude[-1],latEnd, lonEnd)
        
        longitude = [l if l<=180 else l-360 for l in longitude]
    
        """ Increment time with 1 second to assure that Easegrid2 accepts the individual datapoints 
            Assign center date of the third of the month recorded in the file subsequently"""
    
        if s.third_of_month=='1':
            t = np.array([dt.datetime(year, month, 1, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(SID))])
        elif s.third_of_month=='2':
            t = np.array([dt.datetime(year, month, 11, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(SID))])
        elif s.third_of_month=='3':
            t = np.array([dt.datetime(year, month, 21, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(SID))])
        else:
            try:
                t = np.array([dt.datetime(year, month, day, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(SID))])
            except:
                ddt = ((dateEnd-dateStart).total_seconds())/(len(latitude)-1)
                t = np.array([dateStart + dt.timedelta(seconds=ddt*i) for i in range(len(SID))])
                s.pp_flag = 1 # reason exact dates but approximate locations (linear interpolation)

        # import cartopy.crs as ccrs
        # #import numpy as np
        # import matplotlib.pyplot as plt
        # import matplotlib.path as mpath
        # import cartopy.feature as cfeature

        # ax = plt.axes(projection=ccrs.NorthPolarStereo())
        # ax.coastlines()
        # # ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
        # # ax.add_feature(cfeature.OCEAN)        
        # ax.add_feature(cfeature.LAND)
        # ax.gridlines()

        # plot=plt.scatter(longitude, latitude,
        #           s=10, c=SID,alpha=0.5,
        #           transform=ccrs.PlateCarree(),)
        # plt.show()
        
        """ append data to class"""
        s.lat.append(latitude)
        s.lon.append(longitude)
        s.SID.append(SID)
        s.date.append(t)


class SCICEX_cruise:
  def __init__(self, obsID):
    self.obsID = obsID
    self.date = []
    self.lat = []
    self.lon = []
    self.SID = []
    self.third_of_month = []
    self.unc_flag = 2
    self.pp_flag = 2

def mtn(x):
    ## convert month abbreviation to number
    months = {
        'jan': 1,
        'feb': 2,
        'mar': 3,
        'apr':4,
         'may':5,
         'jun':6,
         'jul':7,
         'aug':8,
         'sep':9,
         'oct':10,
         'nov':11,
         'dec':12
        }
    a = x.strip()[:3].lower()
    try:
        ez = months[a]
        return ez
    except:
        raise ValueError('Not a month')

def third_of_mtn(x):
    ## convert third of month to a date
    third_of_months = {
        '1': 5, # date span 1-10
        '2': 15, # date span 11-20
        '3': 25, # date span 21-31
        }
    try:
        ez = third_of_months[x]
        return ez
    except:
        raise ValueError('Not a valid number')

def initial_bearing(pointA, pointB):

    if (type(pointA) != tuple) or (type(pointB) != tuple):
        raise TypeError("Only tuples are supported as arguments")

    lat1 = math.radians(pointA[0])
    lat2 = math.radians(pointB[0])

    
    dl = math.radians(pointB[1] - pointA[1])
    X = math.cos(lat2)*math.sin(dl)
    print(X)
    
    Y = math.cos(lat1)*math.sin(lat2)-math.sin(lat1)*math.sin(lat2)*math.cos(dl)

    initial_bearing = math.atan2(X, Y)
    
    bearing = (np.rad2deg(initial_bearing) + 360) % 360
    print(bearing)

    return np.deg2rad(bearing)

def calc_bearing(lat1, long1, lat2, long2):
  # Convert latitude and longitude to radians
  lat1 = math.radians(lat1)
  long1 = math.radians(long1)
  lat2 = math.radians(lat2)
  long2 = math.radians(long2)
  
  # Calculate the bearing
  bearing = math.atan2(
      math.sin(long2 - long1) * math.cos(lat2),
      math.cos(lat1) * math.sin(lat2) - math.sin(lat1) * math.cos(lat2) * math.cos(long2 - long1)
  )
  
  # Convert the bearing to degrees
  bearing = math.degrees(bearing)
  
  # Make sure the bearing is positive
  bearing = (bearing + 360) % 360
  
  return np.deg2rad(bearing)

def sort(directory, directories):
    combined = []
    relevant = [l for l in directories if os.path.isdir(os.path.join(directory,l))]
    print(relevant)
    relevant = [l for l in relevant if bool(re.search(r'\d', l))]
    print(relevant)
    # relevant = [l for l in relevant if os.path.isdir(os.path.join(directory,l))]
    
    for r in relevant:
        numbers = ''.join([l for l in r if l.isdigit()])
        if len(numbers)==2:
            numbers = '19' + numbers
        elif len(numbers)==3:
            numbers = numbers[-2:]
            numbers = '19' + numbers
        combined.append(numbers)
    print(combined)
    order = np.argsort(combined)
        
    return np.array(relevant)[order]
#%% Main


dtint = 30 # days
gridres = 25000 # resolution of gridded product in meters

campaign = 'SCICEX'
save_directory=save_directory = os.path.dirname(os.path.dirname(os.getcwd())) + '/Final/'
save_path_data= save_directory + 'SCICEX/final/'
if not os.path.exists(save_path_data):os.makedirs(save_path_data)
ofile='ESACCIplus-SEAICE-RRDP2+-SID-SCICEX.nc'
ofile =os.path.join(save_path_data,ofile)
output=open(ofile,'w')
# Host directory of the data
directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/SCICEX/data'
saveplot = os.path.join(os.path.dirname(
    os.path.dirname(os.getcwd())), 'Final/SCICEX/fig/')
if not os.path.exists(saveplot):os.makedirs(saveplot)

# loop throguh directories

directories = sort(directory, os.listdir(directory))
count=0
total_obs_valid_SID = 0
for dirr in directories:
    try:
        dirr_check = int(re.findall('\d+', dirr)[0])
    except:
        print('contains no numbers ', dirr)
    # dirr_check = dirr.translate({ord(i): None for i in 'abcdefghijklmnopqrstuABCDEFGHIJKLMNOPQRS'})
    if os.path.isdir(os.path.join(directory,dirr)): # only enter if is directory
    #if 'scicex' in dirr or 'UK' in dirr or 'L2' in dirr or 'grayling' in dirr: #filenaming convention
        dir_data=os.path.join(directory,dirr)
        if 'scicex' in dirr or 'UK' in dirr or 'L2' in dirr or 'grayling' in dirr: #filenaming convention
            if 'scicex' in dirr:
                name = dirr.replace('scicex', '19') + 'S'
            if 'UK' in dirr:
                name = dirr.replace('UK', '19') + 'U'
            if 'L2' in dirr:
                name = dirr.replace('L2', '19') + 'L'
            if 'grayling' in dirr:
                name = dirr.replace('grayling', '19') + 'G'
            name = name.replace('_', '')
        else:
            name = dirr
        s = SCICEX_cruise('SCICEX-' + name)
        if os.path.isdir(dir_data): # check if file is a directory
            # Get files from data folder
            for ifile, cc in zip(os.listdir(dir_data), range(len(os.listdir(dir_data)))):
                # check if it is a data file
                if 'drft' in ifile and not os.path.isdir(ifile):
                    
                    print(ifile)
                    drft_files(s, dir_data, ifile)
                elif ifile.endswith('_-uwa.series') or ifile.endswith('.Ser') or ifile.endswith('.txt') and not os.path.isdir(ifile):
                    print(ifile)
                    file=os.path.join(dir_data, ifile)
                    
                    ### OBTAIN LAT; LON AND DATE INFORMATION FROM FILE HEADER
                    with open(file, 'r') as the_file:
                        c = 0
                        for line in the_file.readlines() :
                            if ifile.endswith('_-uwa.series') or ifile.endswith('.Ser'):
                                c+=1
                                #print(line)
                                if re.search('Beginning Latitude', line):
                                    latBeg = float(re.findall('\d+\.\d+', line)[0]) # starting lattitude
                                if re.search('Ending Latitude', line):
                                    latEnd = float(re.findall('\d+\.\d+', line)[0]) # ending lattitude
                                if re.search('Beginning Longitude', line):
                                    lonBeg = float(re.findall('\d+\.\d+', line)[0]) # starting longitude
                                if re.search('Ending Longitude', line):
                                    lonEnd = float(re.findall('\d+\.\d+', line)[0]) # ending longitude 
                                ## DATE INFORMATION
                                if re.search('Year', line):
                                    year = int(re.findall('[0-9]+',line)[0]) # starting lattitude
                                if re.search('Month', line) and not re.search('Third of Month', line):
                                    month = re.findall(r':(\w+)', line.strip().replace(' ', ''))[0] # ending lattitude
                                    month = mtn(month)
                                if re.search('Third of Month', line):
                                    print('third of month')
                                    third_of_month = re.findall('[0-9]+',line)[0] # starting longitude
                                    day = third_of_mtn(third_of_month)
                                    s.third_of_month = third_of_month
                                if re.search('Begin', line):
                                    header_end = c
                                if re.search('!!!End', line.strip().replace(' ', '')):
                                    file_end = c-1
                            elif ifile.endswith('.txt'): # if files are newer (2011, 2014)
                                c+=1
                                #print(line)
                                if re.search('Beg lat', line):
                                    latBeg = float(re.findall('\d+\.\d+', line)[0]) # starting lattitude
                                if re.search('End lat', line):
                                    latEnd = float(re.findall('\d+\.\d+', line)[0]) # ending lattitude
                                if re.search('Beg lon', line):
                                    lonBeg = float(re.findall('\d+\.\d+', line)[0]) # starting longitude
                                if re.search('End lon', line):
                                    lonEnd = float(re.findall('\d+\.\d+', line)[0]) # ending longitude 
                                ## DATE INFORMATION                            
                                if re.search('Begin Time', line):
                                    dateStart = dt.datetime.strptime(line.replace('Begin Time:', '').strip(),"%m/%d/%Y %H:%M:%S")
                                if re.search('End Time', line):
                                    dateEnd = dt.datetime.strptime(line.replace('End Time:', '').strip(),"%m/%d/%Y %H:%M:%S")
                                ## find header information
                                if 'distance(m)' in line.strip():
                                    header_end = c+1
                                    
                                # if re.search('!!!End', line.strip().replace(' ', '')):
                                file_end = c-1
                    
                    if lonBeg>180:
                        lonBeg = lonBeg-360
                    if lonEnd>180:
                        lonEnd = lonEnd-360
                    
                    data = np.genfromtxt(file,dtype=None,skip_header=header_end, skip_footer=c-file_end, delimiter=',')
                    try: # format of data strange with "" used to mark the data rows
                        dist = [float(d[0].decode('utf8').replace('"','')) for d in data]
                        SID = [float(d[1].decode('utf8').replace('"','')) for d in data] # sea ice draft
                    except:
                        dist = [d[0] for d in data]
                        SID = [d[1] for d in data]
                    start = (latBeg, lonBeg)
                    end = (latEnd, lonEnd)  

                    
                    ## method for calculating lat, lon
                    """ Initital steps
                        1. Calculate initial bearing between start and end points
                        2. Calculate the coordinates located the informed distance away 
                        (e.g. from the file)
                        3. use these coordinates as new starting point for calculating bearing
                        and to use for the inverse haversine formula
                        Subsequent steps
                        4. Calculate delta distance
                        5. Update bearing
                        6. Update coordinates
                        7. continue until all points have been processed
                    """
                    d_pre = 0 # initial distance is 0
                    latitude = []
                    longitude = []
                    
                    #latitude.append(latBeg)
                    #longitude.append(lonBeg)
                    for d in dist:
                        # calculate distance between points in file
                        delta_d = d - d_pre
                        # calculate the bearing between the coordinates
                        bearing = calc_bearing(latBeg, lonBeg, latEnd, lonEnd)
                        # print('current bearing: ', bearing)
                        # calculate lat, lon from the haversine formula
                        # - great circle distance on a sphere:, https://en.wikipedia.org/wiki/Haversine_formula
                        lat, lon = inverse_haversine(start, delta_d, bearing, unit='m')
                        latitude.append(lat)
                        longitude.append(lon)
                        # update starting coordinate
                        latBeg, lonBeg = lat, lon
                        start = lat, lon
                        # update previous distance
                        d_pre = d
                    
                    if lon<-180:
                        lon = lon+360
                    assert np.round(lat,1),np.round(lon,1) == (latEnd, lonEnd)
                    print(np.round(lat,1),np.round(lon,1))
                    # print('ending lat, lon', lat,lon,latEnd, lonEnd)

                    """ Increment time with 1 second to assure that Easegrid2 accepts the individual datapoints 
                        Assign center date of the third of the month recorded in the file subsequently"""

                    if ifile.endswith('_-uwa.series') or ifile.endswith('.Ser'):
                        if s.third_of_month=='1':
                            t = np.array([dt.datetime(year, month, 1, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(dist))])
                        elif s.third_of_month=='2':
                            t = np.array([dt.datetime(year, month, 11, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(dist))])
                        elif s.third_of_month=='3':
                            t = np.array([dt.datetime(year, month, 21, 0,0,0) + dt.timedelta(seconds=i) for i in range(len(dist))])

                    elif ifile.endswith('.txt'):
                        dtt = dateEnd-dateStart
                        dt_inv = dtt/len(SID)
                        t = [dateStart + dt_inv*tt for tt in range(len(SID))]
                    """ append data to class"""
                    s.lat.append(latitude)
                    s.lon.append(longitude)
                    s.SID.append(SID)
                    s.date.append(t)
                    # s.obsID.append([dirr for i in range(len(latitude))])

                    """ name the cruise - using SCICEX + year + specifier + segment number """
                    if len(ifile.split('-'))>3:
                        obsID = 'SCICEX-' + ifile.split('-')[0] + '-' + ifile.split('-')[2].replace('_', '')
                    elif len(ifile.split('-'))==1:
                        obsID = 'SCICEX-' + ifile.split('_')[0]
                    else:
                        obsID = 'SCICEX-' + ifile.split('-')[0] + '-' + ifile.split('-')[1].replace('_', '')
                    print(ifile)


            
        """ -------------------------------------------- """
        """ Main part of script does processing to fit 25km grid
        and 1 month using EaseGrid2"""
        """ -------------------------------------------- """
        test_drft = any(['drft' in f for f in os.listdir(dir_data)]) # check if any drft files exist in folder
        test_series  =any(['_-uwa.series' in f for f in os.listdir(dir_data)])
        test_ser  =any(['.Ser' in f for f in os.listdir(dir_data)])
        test_txt = any(['.txt' in f for f in os.listdir(dir_data)])
        if test_series or test_ser or test_txt or test_drft:
            latitude = np.concatenate((s.lat))
            longitude = np.concatenate((s.lon))
            SID = np.concatenate((s.SID))
            t = np.concatenate((s.date))
            # obsID = np.concatenate((s.obsID))
            
            #defines variables in output file
            count+=1
            dataOut = Functions.Final_Data(Type='SID', count_head=count)
    
        
            # set non existing data to nan
            SID = np.array(SID)
            SID[SID==-99999.0] = np.nan
            
            total_obs_valid_SID += len(SID[np.isfinite(SID)])
        
        
print(total_obs_valid_SID)
        
        
        
            
            
    
        

            