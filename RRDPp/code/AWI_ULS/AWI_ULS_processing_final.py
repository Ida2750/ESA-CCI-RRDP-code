# -*- coding: utf-8 -*-

"""
Calculates monthly means of input files containing measurements of SID from stationary moorings 
measured at the Weddel Sea Antarctica by Alfred Wegener Institute

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
import sys
import re

# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
sys.path.append(os.path.dirname(os.getcwd()))
import Functions
#%% Functions
def get_coordinates(filename):
    """
    get latitude, longitude and skipheader information from file
    """
    with open(filename, 'r') as f:
        content = f.read()
        values1 = content.split("\n")
        values = content.split("\t")
        index = [i for val, i in zip(values, range(len(values))) if len(
            re.findall(r'\bLA[A-Z]*', val)) > 0]
        lat, lon = re.findall(r"[-+]?(?:\d*\.*\d+)", values[index[0]])
        index = [i for val, i in zip(values1, range(len(values1))) if len(
            re.findall(r'\bDate/Time', val)) > 0]
    return lat, lon, index[-1]

def write_info_table(obsID, timespan, quality):
    """
    Parameters
    ----------
    obsID : String
        DESCRIPTION.
    timespan : TYPE
        DESCRIPTION.
    quality : Integer
        Describes the correction applied to the raw SID
    Returns
    -------
    None.
    """
    ofile = 'AWI-ULS_info_table_'+str(Bias)+'.dat'
    output = open(ofile, 'a')
    ## print data to output
    if quality == 1:
        unc = 3
        best = 'Raw Draft'
    elif quality == 2:
        best = 'Modelled Draft'
        unc = 3
    elif quality == 3:
        best = 'Model correction'
        unc = 1
    elif quality == 4:
        best = 'Zero line correction'
        unc = 1
    print('{:^s} & {:^s} & {:^d} & {:^d} & {:^s} \\'.format(
        obsID, timespan, 0, unc, best), file=output)
    output.close()

    return None

def Bias_correction(SID):
    """
    Parameters
    ----------
    SID : Numpy Array
        SID data array.

    Returns
    -------
    SID : Numpy Array
        SID data array after correction for Bias

    """
    for i in range(len(SID)):
        if SID[i]>0.42:
            if SID[i] <= 1.05:
                SID[i] = SID[i]-0.42
            elif SID[i] <= 1.15:
                SID[i] = SID[i] -0.45
            elif SID[i] <= 1.25:
                SID[i] = SID[i] -0.48
            elif SID[i] <= 1.35:
                SID[i] = SID[i] -0.52
            elif SID[i] <= 1.45:
                SID[i] = SID[i] -0.55
            elif SID[i] <= 1.55:
                SID[i] = SID[i] -0.58
            elif SID[i] <= 1.65:
                SID[i] = SID[i] -0.62
            elif SID[i] <= 1.75:
                SID[i] = SID[i] -0.65
            elif SID[i] > 1.75:
                SID[i] = SID[i] -0.68
        elif SID[i]<0:
            SID[i] = np.nan
    return SID

def robust_std(data):
    median = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - median))
    return mad * 1.4826  # scaling factor for normal distribution
#%% Main

#define output variables
Bias=True
# raw dat directory
#directory = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.getcwd())))) + '/RRDPp/RawData/Antarctic/AWI_ULS/Behrendt_2013/datasets/'
directory = '/dmidata/projects/cmems2/C3S/RRDPp/RawData/Antarctic/AWI_ULS/Behrendt_2013/datasets'
# saving path for output files and figures
savepath = os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'FINAL/Antarctic/AWI_ULS/final/')
saveplot = os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'FINAL/Antarctic/AWI_ULS/fig/')
# name of output file
ofile = 'ESACCIplus-SEAICE-RRDP2+-SID-AWI-ULS_Bias_'+str(Bias)+'-robusttest.nc'
ofile =  os.path.join(savepath, ofile)
if not os.path.exists(savepath):os.makedirs(savepath)
if not os.path.exists(saveplot):os.makedirs(saveplot)


count = 0 # used to locate header
"""
Main Loop through raw data which takes the monthly average and writes output data to ofile
""" 

for ifile in os.listdir(directory):
    if ifile.endswith('.tab'): 
        print(ifile)
        file=os.path.join(directory, ifile)
        # initiate object from final data class
        count += 1
        dataOut = Functions.Final_Data(Type='SID', count_head=count)
        # Assign pre-processing flag
        dataOut.pp_flag = 0
        
        # get latitude, longitude information from file header
        lat, lon, header_count = get_coordinates(file)
        # Read data from ASCII-file
        try:
            dtype = [object, float, float, float, float,float, float, float, float]
            data = np.genfromtxt(file,dtype=dtype,skip_header=header_count,names=True,delimiter='\t')
        except: # fewer arrays in file 229-8 - something off with the formatting
             dtype = [object, float, float]
             Names = ['DateTime', 'Draft_m_Upward_looking_sonar_ULS', 'Draft_m_Calculated']
             data = np.genfromtxt(file,skip_header=header_count+1,names=Names, usecols=(0,3,4), dtype=dtype) #,delimiter='\t')

        date = data['DateTime']
        months = [];
        dates = [];
        # Create datetime array, seconds array and month array
        try:
            dates = [dt.datetime.strptime(dat.decode("utf-8"), '%Y-%m-%dT%H:%M') for dat in date]
        except: # off formatting in 229-8 file
            dates = [dt.datetime.strptime(dat.decode("utf-8"), '%Y-%m-%dT%H:%M:%S') for dat in date]
        # Define observation identifier (obsID)
        splitted = re.split('- |_|-|!', ifile)
        try:
            int(splitted[1])
            dataOut.obsID= splitted[0] + '-' + splitted[1] #'AWI-' + str(dates[0]).split()[0].replace('-','')
        except:
            dataOut.obsID= splitted[0]
        
        try:
            months=[int(dt.datetime.strptime(dat.decode("utf-8"), '%Y-%m-%dT%H:%M').strftime("%m")) for dat in date]
        except: # off formatting in 229-8 file
           months=[int(dt.datetime.strptime(dat.decode("utf-8"), '%Y-%m-%dT%H:%M:%S').strftime("%m")) for dat in date]
 
         # seconds array
        t2=np.array([(date-dt.datetime(1970,1,1)).total_seconds() for date in dates])

        #load SID DATA
        try:
            SID_raw=data['Draft_m_Upward_looking_sonar_ULS'] #meters
            unc_flag = [3 for SID in SID_raw]
        except ValueError:
            SID_raw=[]
        try:
            SID_calculated=data['Draft_m_Calculated'] #meters
            unc_flag = [1 for SID in SID_calculated]
        except ValueError:
            SID_calculated=[]
        try:    
            #Draft [m] (zero line correction)	Flag (zero line correction)	Draft [m] (model correction)	Flag (model correction)
            SID_zero_corr=data['Draft_m_zero_line_correction']#meters
            unc_flag = [1 for SID in SID_zero_corr]
        except ValueError:
            SID_zero_corr=[]
        try:
            SID_model_corr = data['Draft_m_model_correction']
            unc_flag = [1 for SID in SID_model_corr]
        except ValueError:
             SID_model_corr=[]
        
        ## get the best version of SID avalible from file
        # and define uncertainty
        if len(SID_zero_corr)!=0 or len(SID_calculated)!=0:
            if len(SID_zero_corr)!=0:
                SID = SID_zero_corr
            else:
                SID = SID_calculated              
            if '206-4' in ifile or '227-3' in ifile:
                unc_flag = [3 for S in SID]
                
            # define uncertainty
            SID_Unc = np.zeros(SID.shape)
            for i, dat in zip(range(len(dates)), dates):
                if dat.month < 6 or dat.month > 10:
                    SID_Unc[i] = 0.05
                else:
                    SID_Unc[i] = 0.12
        elif len(SID_model_corr)!=0:
            SID = SID_model_corr
            # define uncertainty
            SID_Unc = np.array([0.23 for S in SID])
            

        ## Bias correction based on table4 A. Behrendt el. al 2013
        if Bias==True:
            SID = Bias_correction(SID)

        SID[SID<=0] = np.nan
        SID[SID>8] = np.nan

        # Find indexes of new months
        mondiff=np.where(~(np.diff(months) == 0))[0] #index of when month changes
        index=np.insert(np.append(mondiff,len(months)-1),0,0) #add start and end indexes on month
        
        # Get average value of variables for each month        
        # loop over variables, but assure that the length takes only the non nan value elements
        dataOut.SID_final=np.array([np.nanmedian(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
        #dataOut.SID_std=np.array([np.nanstd(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])
        dataOut.SID_std = np.array([robust_std(SID[index[i]:index[i+1]]) for i in range(len(index)-1)])

        dataOut.SID_ln=np.array([len(SID[index[i]:index[i+1]][~np.isnan(SID[index[i]:index[i+1]])]) for i in range(len(index)-1)])
        dataOut.unc_flag = np.array([unc_flag[index[i]] for i in range(len(index)-1)])
        
        # Get uncertainty of average
        dataOut.SID_unc = np.zeros(dataOut.SID_final.shape)
        for i in range(len(index)-1):
            start = index[i]
            end = index[i+1]
            dataOut.SID_unc[i] = 1/dataOut.SID_ln[i] * np.sqrt(np.nansum(SID_Unc[start:end]**2))

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
        
        # Get median data of each month
        avgDates=np.array([np.median(t2[index[i]:index[i+1]]) for i in range(len(index)-1)])
        dates_final=np.array([dt.datetime.fromtimestamp(int(sec)) for sec in avgDates])
        dataOut.date_final=np.array([dt.datetime.strftime(date,"%Y-%m-%dT%H:%M:%S") for date in dates_final])
        
        #define values of output variables
        dataOut.lat_final = [float(lat) for num in dataOut.SID_final]
        dataOut.lon_final = [float(lon) for num in dataOut.SID_final]
        
        Functions.scatter(dataOut.obsID, dates, SID, dates_final, dataOut.SID_final, 'SID [m]', saveplot)
        
        dataOut.time = [np.datetime64(d) for d in dataOut.date_final]
        dataOut.pp_flag = [dataOut.pp_flag]*len(dataOut.SID_final)
        dataOut.obsID = [dataOut.obsID]*len(dataOut.SID_final)
        
        # fill empty arrays with NaN values
        dataOut.Check_Output()

        if count>1:
            subset = dataOut.Create_NC_file(ofile, primary='SID')
            df = Functions.Append_to_NC(df, subset)
        else:
            df = dataOut.Create_NC_file(ofile, primary='SID', datasource='Sea ice draft measured by upward looking sonars in the Weddell Sea (Antarctica): https://doi.org/10.1594/PANGAEA.785565', key_variables='Sea Ice Draft')
            

# Sort final data based on date
Functions.save_NC_file(df, ofile, primary='SID')         
