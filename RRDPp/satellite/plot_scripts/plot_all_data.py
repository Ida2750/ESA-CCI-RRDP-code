# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2023-07-25'


# -- Built-in modules -- #
# import sys
import os
# -- Third-part modules -- #
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
# -- Proprietary modules -- #
import polar_plots as pp
from plotting_functions import Histograms, scatter, Data, stats

plt.style.use('tableau-colorblind10')
## colorblind colors 
colors = ['#006BA4', '#FF800E', '#ABABAB',
                  '#595959', '#5F9ED1', '#C85200',
                  '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']
# '006BA4'	Cerulean/Blue
# 'FF800E'	Pumpkin/Orange
# 'ABABAB'	Dark Gray/Gray
# '595959'	Mortar/Grey
# '5F9ED1'	Picton Blue/Blue
# 'C85200'	Tenne (Tawny)/Orange
# '898989'	Suva Grey/Grey
# 'A2C8EC'	Sail/Blue
# 'FFBC79'	Macaroni And Cheese/Orange

# Choose a satellite (CS2 or ENV)
sat = 'ENV'
ofile = os.getcwd() + '/' + sat + '_statistics.dat'
output=open(ofile,'w')
directory = os.path.dirname(os.getcwd()) + "/Final_files/"
files = os.listdir(directory)

if sat=='ENV':
    cc = ['orange', 'g']
else:
    cc = ['blue', 'r']

# start a data class for each variable type
d_SID = Data(sat)
d_SID.var = 'SID'
d_SIT = Data(sat)
d_SIT.var = 'SIT'
d_SD = Data(sat)
d_SD.var = 'SD'
d_SIF = Data(sat)
d_SIF.var = 'SIF'

# loop through all files and append data to data class
for ifile, i in zip(files, range(len(files))):
    if ifile.startswith('ESACCI') and sat in ifile.replace('-',''):
        print(ifile)
        ## Get name
        name = ''.join(filter(str.isalnum, ifile))
        name = name.replace("ESACCIplusSEAICERRDP2SIT", "")
        name = name.replace("ESACCIplusSEAICERRDP2SID", "")
        name = name.replace("ESACCIplusSEAICERRDP2FRB", "")
        name = name.replace("ESACCIplusSEAICERRDP2SD", "")
        name = name.replace("CCIpv3p0rc2dat", "")
        name = name.replace("IDCS4", "")
        name = name.replace("AWI", "")
        name = name.replace(sat, "")
        name = name.replace("ANTARCTIC", "-SH") 
        name = name.replace("BANT", "B-SH") 
        name = name.replace("sorted", "")   
        print(name)
        if 'OIB' in name:
            c=colors[2] # 'grey'
        elif 'ASSIST' in name or 'ASPeCt' in name:
            c=colors[1] # 'orange'
        elif 'IMB' in name:
            name = name.replace('B', 'B-CRREL')
            c="#832db6" # purple
        elif 'AEM' in name:
            name = name.replace('M', 'M-AWI')
            c=colors[0] # blue
        elif 'SB' in name or 'SB-SH' in name:
            name = name.replace('B', 'B-AWI')
            c=colors[4] # light blue
        ## SID 
        elif 'SCICEX' in name:
            c=colors[5] # dark orange
        elif 'BGEP' in name:
            c=colors[6] # grey
        elif 'NPI' in name:
            name = name + '-FS'
            c=colors[0] # blye
        elif 'TRANSDRIFT' in name:
            c=colors[1] # orange
        elif 'NP' in name or 'ULS' in name:
            c='k' # black
        else:
            c=colors[0]
        
        
        if 'SID' in ifile:
            ObsNames=['date','lat','lon','obsSID', 'satSID', 'obsSIDstd','obsSIDln','obsSIDunc',
                      'satSIDstd','satSIDln', 'satSIDunc']
        
            ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None,names=ObsNames)
        
            CRS_ENV = np.arange(0,len(ObsData))
            obsDate = ObsData['date'][CRS_ENV]
            obslat  = ObsData['lat'][CRS_ENV]
            obslon  = ObsData['lon'][CRS_ENV]
            obsSID  = ObsData['obsSID'][CRS_ENV]
            satSID  = ObsData['satSID'][CRS_ENV]                
            m=~np.isnan(obsSID)

            
            if any(m):
                print(d_SID.var)
                d_SID.lat.append(obslat[m])
                d_SID.lon.append(obslon[m])
                d_SID.obsSID.append(obsSID[m])
                d_SID.satSID.append(satSID[m])
                d_SID.name.append(name)
                d_SID.c.append(c)

        
        
        elif 'SIT' in ifile or 'SD' in ifile or 'FRB' in ifile:
            ObsNames=['date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT',
                      'satSIF', "obsSD_std","obsSD_ln", "obsSD_unc","obsSIT_std","obsSIT_ln", "obsSIT_unc","obsFRB_std",
                      "obsFRB_ln", "obsFRB_unc", "satSD_std","satSD_ln", "satSD_unc","satSIT_std","satSIT_ln", "satSIT_unc","satFRB_std",
                      "satFRB_ln", "satFRB_unc"]
            ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None,names=ObsNames)
        
        #try:
            CRS_ENV = np.arange(0,len(ObsData))
            obsDate = ObsData['date'][CRS_ENV]
            obslat  = ObsData['lat'][CRS_ENV]
            obslon  = ObsData['lon'][CRS_ENV]
            obsSIF  = ObsData['obsSIF'][CRS_ENV]
            obsSD   = ObsData['obsSD'][CRS_ENV]
            obsSIT  = ObsData['obsSIT'][CRS_ENV]
            satSIF  = ObsData['satSIF'][CRS_ENV]
            satSIT  = ObsData['satSIT'][CRS_ENV]
            satSD   = ObsData['satSD'][CRS_ENV]
            if 'OIB' in name and 'SH' in name or 'AEM' in name: # convert sea ice frb to total frb
                satSIF = satSIF + satSD
            elif 'OIB' in ifile:
                obsSIF = obsSIF - obsSD
            if 'AEM-AWI' in ifile:
                satSIT = satSIT + satSD # AEM measures SIT+SD
            
            m=~np.isnan(obsSD)
            n=~np.isnan(obsSIT)
            o=~np.isnan(obsSIF)
            
            if any(n):
                d_SIT.lat.append(obslat[n])
                d_SIT.lon.append(obslon[n])
                d_SIT.obsSIT.append(obsSIT[n])
                d_SIT.satSIT.append(satSIT[n])
                d_SIT.name.append(name)
                # d_SIT.var.append('SIT')
                d_SIT.c.append(c)
            
                

            if any(m):
                d_SD.lat.append(obslat[m])
                d_SD.lon.append(obslon[m])
                d_SD.obsSD.append(obsSD[m])
                d_SD.satSD.append(satSD[m])
                d_SD.name.append(name)
                # d_SD.var.append('SD')
                d_SD.c.append(c)

            if any(o):
                d_SIF.lat.append(obslat[o])
                d_SIF.lon.append(obslon[o])
                d_SIF.obsSIF.append(obsSIF[o])
                d_SIF.satSIF.append(satSIF[o])
                d_SIF.name.append(name)
                # d_SIF.var.append('SIF')
                d_SIF.c.append(c)


## sort data for visualisation purposes
d_SIT.sort()
d_SD.sort()
d_SID.sort()
d_SIF.sort()
# calculate statistics of data
stats(d_SID,ofile)
stats(d_SIT,ofile)
stats(d_SD,ofile)
stats(d_SIF,ofile)

#%% Make histograms
Histograms(d_SID)
Histograms(d_SIT)
Histograms(d_SIF)
Histograms(d_SD)

#%% Make scatterplots
scatter(d_SID)
scatter(d_SIT)
scatter(d_SIF)
scatter(d_SD)
#%% Polar plots
bool_list = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIF.name])
savename =  'FRB_data_'+sat+'.png'
title = 'FRB Validation Data for ' + sat
pp.plot_all(d_SIF.lat[bool_list], d_SIF.lon[bool_list], d_SIF.obsSIF[bool_list], title=title, ylabel=d_SIF.name[bool_list], s=5, c=d_SIF.c[bool_list], savename=savename)


bool_list = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIT.name])
savename =  'SIT_data_'+sat+'.png'
title = 'SIT Validation Data for ' + sat
pp.plot_all(d_SIT.lat[bool_list], d_SIT.lon[bool_list], d_SIT.obsSIT[bool_list], title=title, ylabel=d_SIT.name[bool_list], s=5, c=d_SIT.c[bool_list], savename=savename)


bool_list = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SD.name])
savename = 'SD_data_'+sat+'.png'
title = 'SD Validation Data for '  + sat
pp.plot_all(d_SD.lat[bool_list], d_SD.lon[bool_list], d_SD.obsSD[bool_list], title=title, ylabel=d_SD.name[bool_list], c=d_SD.c[bool_list], savename=savename)

bool_list = np.invert(['ASPeCt'==name or 'SH' in name or 'ULS' in name for name in d_SID.name])
savename = 'SID_data_'+sat+'.png'
s = 20
title = 'SID Validation Data for '  + sat
pp.plot_all(d_SID.lat[bool_list], d_SID.lon[bool_list], d_SID.obsSID[bool_list], title=title, ylabel=d_SID.name[bool_list], c=d_SID.c[bool_list], savename=savename)

#%%
## ANTARCTIC
bool_list_SIT = ['ASPeCt'==name or 'SH' in name for name in d_SIT.name]
bool_list_SD = ['ASPeCt'==name or 'SH' in name for name in d_SD.name]
bool_list_SID = ['ASPeCt'==name or 'SH' in name or 'AWI-ULS' in name for name in d_SID.name]
bool_list_SIF = ['ASPeCt'==name or 'SH' in name for name in d_SIF.name]
savename =  'SH Validation_data_'+sat+'.png'
title = 'Validation Data for ' + sat

lats = np.concatenate([d_SIT.lat[bool_list_SIT], d_SD.lat[bool_list_SD], d_SID.lat[bool_list_SID], d_SIF.lat[bool_list_SIF]])
lons = np.concatenate([d_SIT.lon[bool_list_SIT], d_SD.lon[bool_list_SD], d_SID.lon[bool_list_SID], d_SIF.lon[bool_list_SIF]])
obs = np.concatenate([d_SIT.obsSIT[bool_list_SIT], d_SD.obsSD[bool_list_SD], d_SID.obsSID[bool_list_SID], d_SIF.obsSIF[bool_list_SIF]])
names = np.concatenate([d_SIT.name[bool_list_SIT], d_SD.name[bool_list_SD], d_SID.name[bool_list_SID], d_SIF.name[bool_list_SIF]])
colors = np.concatenate([d_SIT.c[bool_list_SIT], d_SD.c[bool_list_SD], d_SID.c[bool_list_SID], d_SIF.c[bool_list_SIF]])
pp.plot_all(lats, lons, obs, title=title, ylabel=names, s=5, c=colors, savename=savename, NP=False)
