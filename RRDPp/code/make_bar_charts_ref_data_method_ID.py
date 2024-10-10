# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:44:25 2023

@author: Ida Olsen
"""

"""
Title: Sea Ice Data Analysis
Author: Ida Olsen
Description: This script performs analysis on Sea Ice data, including Sea Ice Thickness (SIT), Snow Depth (SD), Freeboard (FRB), and Sea Ice Drift (SID).
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['ilo@dmi.dk']
__version__ = '1'
__date__ = '2024-10-10'

# -- Built-in modules -- #
import os
import glob
import datetime as dt
import xarray as xr

# -- Third-part modules -- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Data:
    def __init__(self):
        """
        Initialize Data class with empty lists to store data.
        """
        self.name = []
        self.var = ''
        self.date = []
        self.lat = []
        self.lon = []
        self.obs_sit = []
        self.obs_sd = []
        self.obs_sid = []
        self.obs_sif = []

def get_name(ifile):
    ## Get name
    filename = ifile.split("\\")[-1]
                           #[-1]
    print(filename)
    # name = ''.join(filter(str.isalnum, ifile))
    name = filename.replace("ESACCIplus-SEAICE-RRDP2+-SIT-", "")
    name = name.replace("ESACCIplus-SEAICE-RRDP2+-SID-", "")
    name = name.replace("ESACCIplus-SEAICE-RRDP2+-FRB-", "")
    name = name.replace("ESACCIplus-SEAICE-RRDP2+-SD-", "")
    #name = name.replace("-IDCS4", "")
    #name = name.replace("-QL", "")
    name = name.replace("I-ANT", "I-SH") 
    name = name.replace("sorted", "")
    name = name.replace(".dat", "")
    name = name.replace(".nc", "")
    
    print(name)
    
    return name

def load_data_sid(self, ifile):
    name = get_name(ifile)
    # Load data for Sea Ice Drift (SID)
    #obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
    #print(ifile)
    obs_data = xr.open_dataset(ifile)
    #crs_env = np.arange(0, len(obs_data))
    obs_date = obs_data['time'].to_numpy()#[crs_env]
    obs_lat = obs_data['lat'].to_numpy()#[crs_env]
    obs_lon = obs_data['lon'].to_numpy()#[crs_env]
    obs_sid = obs_data['SID'].to_numpy()#[crs_env]
    print(obs_sid)
    #print(obs_sid)
    m = ~np.isnan(obs_sid)

    if any(m):
        self.lat.append(obs_lat[m])
        self.lon.append(obs_lon[m])
        self.obs_sid.append(obs_sid[m])
        self.date.append(obs_date[m])
        self.name.append(name)

def load_data_sit(self, ifile):
    
    name = get_name(ifile)
    # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
    #obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
    #print(ifile)
    obs_data = xr.open_dataset(ifile)
    #print(obs_data)
    #crs_env = np.arange(0, len(obs_data))
    obs_date = obs_data['time'].to_numpy()#[crs_env]
    obs_lat = obs_data['lat'].to_numpy()#[crs_env]
    obs_lon = obs_data['lon'].to_numpy()#[crs_env]
    obs_sif = obs_data['FRB'].to_numpy()#[crs_env]
    obs_sd = obs_data['SD'].to_numpy()#[crs_env]
    obs_sit = obs_data['SIT'].to_numpy()#[crs_env]
    
    m = ~np.isnan(obs_sd)
    n = ~np.isnan(obs_sit)
    o = ~np.isnan(obs_sif)
    
    if any(n) and self.var=='SIT':
        self.obs_sit.append(obs_sit[n])
        self.date.append(obs_date[n])
        self.name.append(name)

    if any(m) and self.var=='SD':
        self.obs_sd.append(obs_sd[m])
        self.date.append(obs_date[m])
        self.name.append(name)

    if any(o) and self.var=='SIF':
        self.obs_sif.append(obs_sif[o])
        self.date.append(obs_date[o])
        self.name.append(name)
    

def get_data(HS):
    directory = "C:/Users/Ida Olsen/Documents/work/ESA-CCI-RRDP-code/RRDPp/FINAL/"
    if HS != 'NH':
        directory += "Antarctic"
    
    files = glob.glob(f'{directory}/*/*/ESACCI*.nc')
    #files = [f for f in os.listdir(directory) if f.startswith('ESACCI') ]
    #print(files)
    d_sid = Data()
    d_sid.var = 'SID'
    d_sit = Data()
    d_sit.var = 'SIT'
    d_sd = Data()
    d_sd.var = 'SD'
    d_sif = Data()
    d_sif.var = 'SIF'
    
    
    for ifile in files:
        if 'SID' in ifile:
            load_data_sid(d_sid, ifile)

        elif 'SIT' in ifile or 'SD' in ifile or 'FRB' in ifile:
            
            [load_data_sit(d, ifile) for d in [d_sit, d_sd, d_sif]]

    return [d_sit, d_sd, d_sif, d_sid]


def plot_distribution_bar(data, title, x_label,HS):
    """
    Function to plot a stacked bar chart for the given data.

    Parameters:
    data (dict): A dictionary containing data to be plotted.
    title (str): Title of the plot.
    x_label (str): Label for the x-axis.

    Returns:
    None
    """
    from matplotlib.ticker import MultipleLocator
    from itertools import cycle, islice

    plt.rcParams.update({'font.size': 14})
    # colors =[u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd']
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ## colorblind colors 
    my_colors = ['#006BA4', '#FF800E', '#ABABAB',
                      '#595959', '#5F9ED1', '#C85200',
                      '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']
    #my_colors = list(islice(cycle(colors), None, len(data)))
    #print(my_colors)
    if HS == 'NH':
        ylim = [0,6000]
        c = my_colors
        
    else:
        ylim = [0,2000]
        data.pop('Submarine')
        c = [my_colors[0]] +  my_colors[2:]
        print(c)
    df = pd.DataFrame(data).sort_index()

    d = df.plot(kind='bar',
                figsize=(10, 8),
                fontsize=14,
                title=title,
                rot=45,
                grid=True,
                ylim=(ylim),
                edgecolor='white',
                width=0.8,
                linewidth=1.5,
                xlabel='Months',
                #legend=False,
                color=c,
                ylabel='Observation count',
                ).legend(loc='upper right')
    
    plt.show()


def identify_observation_type(n):
    if 'AEM' in n or 'OIB' in n:
        ID = 'Airborne'
    elif n=='NP':
        ID = 'DS'
    elif 'SCICEX' in n:
        ID = 'Submarine'
    elif 'ASSIST' in n or 'ASPeCt' in n:
        ID = 'Ship'
    elif 'IMB' in n or 'SB' in n:
        ID = 'Drifting Buoy'
    else:
        ID = 'Mooring'
        print(n)
    
    return ID
    
    
def Avoid_Duplicates(self):
    # identify number of observations
    l_out = 0
    for el in self:
        index = np.where(el.name=='OIB')[0]
        l = len(np.concatenate((el.date[index])))
        if l>l_out:
            l_out = l
            var_count = el.var
            
            
                
            
            # Most FRB measurements.
            
        

def histograms_time(self, HS):
    """
    Function to calculate and plot the monthly and yearly distribution of observations.

    Parameters:
    d_sit (Data): Data instance for Sea Ice Thickness (SIT) observations.
    d_sd (Data): Data instance for Snow Depth (SD) observations.
    d_sif (Data): Data instance for Freeboard (FRB) observations.
    d_sid (Data): Data instance for Sea Ice Drift (SID) observations.

    Returns:
    None
    """
    # if self.var=='SID':
    #     categories = {'ULS': {}, 'Submarine': {}}
    #else:
    categories = {'Mooring': {m: 0 for m in range(1, 13)},
                  'Submarine': {m: 0 for m in range(1, 13)},
                  'Airborne': {m: 0 for m in range(1, 13)},
                  'Drifting Buoy': {m: 0 for m in range(1, 13)},
                  'Ship': {m: 0 for m in range(1, 13)}}
    
    names = {}
    rel_var = {}
    for el in self:
        for n in el.name:
            names.update({n:0})
            rel_var.update({n:''})
    
    for el in self:
        for n, i in zip(el.name, range(len(el.name))):
            if names[n] < len(el.date[i]):
                
                names[n] = len(el.date[i])
                print(n, names[n])
                rel_var[n] = el.var
                
                
    
    
    for el in self:
        ## avoid duplicate counts
        for n, i in zip(el.name, range(len(el.name))):
            if rel_var[n]==el.var:

                # Extract dates, years, and months from the Data instances
                dates = [dt.datetime.strptime(str(dd)[:26], '%Y-%m-%dT%H:%M:%S.%f') for dd in el.date[i]]
                months = [d.month for d in dates]
                
                ID = identify_observation_type(n)
                if ID!='DS':
                    for m in months:
                        categories[ID][m] += 1

       
    # # Plot the monthly distribution
    plot_distribution_bar(categories, 'Monthly distribution of observations, ' + HS, 'Months', HS)

# Set Hemisphere
HS = 'NH'
NH_data = get_data('NH')
SH_data = get_data('SH')

# Call the histograms_time function with your
# for el in NH_data:
histograms_time(NH_data, 'NH')
histograms_time(SH_data, 'SH')