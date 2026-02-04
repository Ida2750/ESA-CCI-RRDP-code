# -*- coding: utf-8 -*-
"""
Created on Fri May 23 13:01:15 2025

@author: rmfha
"""

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.dates as mdates
import cartopy.feature as cfeature
import cartopy
from matplotlib import colors
import scipy.signal as scipy
import matplotlib.font_manager as fm
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
from sklearn.linear_model import LinearRegression
from scipy.interpolate import griddata
import numpy as np
from datetime import datetime, timedelta
import copy
from math import radians, cos, sin, asin, sqrt
import os
import numpy as np
import proplot as pplt
import h5py
import pandas as pd
import netCDF4
import sys
# import dates
from matplotlib.colors import LogNorm
from scipy import signal
import cartopy.crs as ccrs


#%%% Buoy time-series ... flagging before monthly

fig,ax = pplt.subplots(ncols=2, nrows=2, sharex=True, refwidth=5, refheight=1, sharey=False)
fig.patch.set_facecolor('white')

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


axs = ax[0]
QFT_val = 4

fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SD-SB_AWI-NH-CS2-NH-CCIp-v3p0-rc2.dat'
df_AWI_SD_NH = pd.read_csv(fp, sep = " ")
#fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-SB_AWI-NH-CS2-NH-CCIp-v3p0-rc2.dat'
#df_AWI_SD_NH_v2 = pd.read_csv(fp, sep = " ")
#df_AWI_SD_NH = pd.concat([df_AWI_SD_NH_v1, df_AWI_SD_NH_v2])
df_AWI_SD_NH = df_AWI_SD_NH[df_AWI_SD_NH['QFT']<=QFT_val]
df_AWI_SD_NH['datetime'] = pd.to_datetime(df_AWI_SD_NH['date'],format="%Y-%m-%dT%H:%M:%S") 
df_AWI_SD_NH_filt = df_AWI_SD_NH.resample('M', on='datetime').mean()
df_AWI_SD_NH_filt['datetime'] = df_AWI_SD_NH_filt.index

df_sat, val1, val2 = df_AWI_SD_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])

df = df_AWI_SD_NH_filt
axs.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
axs.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
axs.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='SB-AWI, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), title ="QFT$\leq${}".format(QFT_val))

ax1 = axs.panel('top', space=0, width=1)
#fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-IMB-ENV-NH-CCIp-v3p0-rc2.dat'
#df_IMB_NH_v1 = pd.read_csv(fp, sep = " ")
fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-IMB-CS2-NH-CCIp-v3p0-rc2.dat'
#df_IMB_NH_v2 = pd.read_csv(fp, sep = " ")
df_IMB_NH = pd.read_csv(fp, sep = " ")
#df_IMB_NH = pd.concat([df_IMB_NH_v1, df_IMB_NH_v2])
df_IMB_NH_alm = df_IMB_NH
df_IMB_NH = df_IMB_NH[df_IMB_NH['QFT']<=QFT_val]
df_IMB_NH['datetime'] = pd.to_datetime(df_IMB_NH['date'],format="%Y-%m-%dT%H:%M:%S") 

df = df_IMB_NH
start_date = '2017-08-01 00:00:00'
end_date = '2018-06-01 00:00:00'
mask = (df['datetime'] > start_date) & (df['datetime'] <= end_date)

df_IMB_NH[mask] = np.nan

df_IMB_NH_filt = df_IMB_NH.resample('M', on='datetime').mean()
df_IMB_NH_filt['datetime'] = df_IMB_NH_filt.index

df_sat, val1, val2 = df_IMB_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


df = df_IMB_NH_filt
ax1.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
ax1.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
ax1.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='IMB, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE))

ax2 = axs.panel('top', space=0, width=1)
fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-MOSAIC-SIMBA-CS2-NH-CCIp-v3p0-rc2.dat'
df_MOSAIC_NH = pd.read_csv(fp, sep = " ")
df_MOSAIC_NH = df_MOSAIC_NH[df_MOSAIC_NH['QFT']<=QFT_val]
df_MOSAIC_NH['datetime'] = pd.to_datetime(df_MOSAIC_NH['date'],format="%Y-%m-%dT%H:%M:%S") 
df_MOSAIC_NH_filt = df_MOSAIC_NH.resample('M', on='datetime').mean()
df_MOSAIC_NH_filt['datetime'] = df_MOSAIC_NH_filt.index

df_sat, val1, val2 = df_MOSAIC_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


df = df_MOSAIC_NH_filt
ax2.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
ax2.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
ax2.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='MOSAIC SIMBA, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE))


axs = ax[2]
df = df_IMB_NH_filt

df_sat, val1, val2 = df_IMB_NH_filt, 'obsSIT', 'satSIT'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


axs.errorbar(x = df['datetime'], y=df['obsSIT'], yerr=df['obsSIT_std'], zorder=10)
axs.errorbar(x = df['datetime'], y=df['satSIT'], yerr=df['satSIT_std'])
axs.format(ylabel='sea ice\nthickness (m)', ultitle='IMB, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), ylim=(-0.1, 4.3), title ="QFT$\leq${}".format(QFT_val))

ax3 = axs.panel('top', space=0, width=1)
df = df_MOSAIC_NH_filt

df_sat, val1, val2 = df_MOSAIC_NH_filt, 'obsSIT', 'satSIT'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


ax3.errorbar(x = df['datetime'], y=df['obsSIT'], yerr=df['obsSIT_std'], zorder=10)
ax3.errorbar(x = df['datetime'], y=df['satSIT'], yerr=df['satSIT_std'])
ax3.format(ylabel='sea ice\n thickness (m)',  ultitle='MOSAIC SIMBA, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), ylim=(-0.1, 4.3))


axs = ax[1]
QFT_val = 1

fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SD-SB_AWI-NH-CS2-NH-CCIp-v3p0-rc2.dat'
df_AWI_SD_NH = pd.read_csv(fp, sep = " ")
#fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-SB_AWI-NH-CS2-NH-CCIp-v3p0-rc2.dat'
#df_AWI_SD_NH_v2 = pd.read_csv(fp, sep = " ")
#df_AWI_SD_NH = pd.concat([df_AWI_SD_NH_v1, df_AWI_SD_NH_v2])
df_AWI_SD_NH = df_AWI_SD_NH[df_AWI_SD_NH['QFT']<=QFT_val]
df_AWI_SD_NH['datetime'] = pd.to_datetime(df_AWI_SD_NH['date'],format="%Y-%m-%dT%H:%M:%S") 
df_AWI_SD_NH_filt = df_AWI_SD_NH.resample('M', on='datetime').mean()
df_AWI_SD_NH_filt['datetime'] = df_AWI_SD_NH_filt.index

df_sat, val1, val2 = df_AWI_SD_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])

df = df_AWI_SD_NH_filt
axs.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
axs.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
axs.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='SB-AWI, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), title ="QFT$\leq${}".format(QFT_val))

ax1 = axs.panel('top', space=0, width=1)
#fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-IMB-ENV-NH-CCIp-v3p0-rc2.dat'
#df_IMB_NH_v1 = pd.read_csv(fp, sep = " ")
fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-IMB-CS2-NH-CCIp-v3p0-rc2.dat'
#df_IMB_NH_v2 = pd.read_csv(fp, sep = " ")
df_IMB_NH = pd.read_csv(fp, sep = " ")
#df_IMB_NH = pd.concat([df_IMB_NH_v1, df_IMB_NH_v2])

df_IMB_NH = df_IMB_NH[df_IMB_NH['QFT']<=QFT_val]
df_IMB_NH['datetime'] = pd.to_datetime(df_IMB_NH['date'],format="%Y-%m-%dT%H:%M:%S") 

df = df_IMB_NH
start_date = '2017-08-01 00:00:00'
end_date = '2018-06-01 00:00:00'
mask = (df['datetime'] > start_date) & (df['datetime'] <= end_date)

df_IMB_NH[mask] = np.nan

df_IMB_NH_filt = df_IMB_NH.resample('M', on='datetime').mean()
df_IMB_NH_filt['datetime'] = df_IMB_NH_filt.index

df_sat, val1, val2 = df_IMB_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


df = df_IMB_NH_filt
ax1.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
ax1.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
ax1.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='IMB, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE))

ax2 = axs.panel('top', space=0, width=1)
fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-MOSAIC-SIMBA-CS2-NH-CCIp-v3p0-rc2.dat'
df_MOSAIC_NH = pd.read_csv(fp, sep = " ")
df_MOSAIC_NH = df_MOSAIC_NH[df_MOSAIC_NH['QFT']<=QFT_val]
df_MOSAIC_NH['datetime'] = pd.to_datetime(df_MOSAIC_NH['date'],format="%Y-%m-%dT%H:%M:%S") 
df_MOSAIC_NH_filt = df_MOSAIC_NH.resample('M', on='datetime').mean()
df_MOSAIC_NH_filt['datetime'] = df_MOSAIC_NH_filt.index

df_sat, val1, val2 = df_MOSAIC_NH_filt, 'obsSD', 'satSD'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


df = df_MOSAIC_NH_filt
ax2.errorbar(x = df['datetime'], y=df['obsSD'], yerr=df['obsSD_std'], zorder=10)
ax2.errorbar(x = df['datetime'], y=df['satSD'], yerr=df['satSD_std'])
ax2.format(ylabel='snow depth (m)', ylim=(-0.1,1.1), ultitle='MOSAIC SIMBA, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE))


axs = ax[3]
df = df_IMB_NH_filt

df_sat, val1, val2 = df_IMB_NH_filt, 'obsSIT', 'satSIT'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


axs.errorbar(x = df['datetime'], y=df['obsSIT'], yerr=df['obsSIT_std'], zorder=10)
axs.errorbar(x = df['datetime'], y=df['satSIT'], yerr=df['satSIT_std'])
axs.format(ylabel='sea ice\nthickness (m)', ultitle='IMB, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), ylim=(-0.1, 4.3), title ="QFT$\leq${}".format(QFT_val))

ax3 = axs.panel('top', space=0, width=1)
df = df_MOSAIC_NH_filt

df_sat, val1, val2 = df_MOSAIC_NH_filt, 'obsSIT', 'satSIT'
df_sat = df_sat[(df_sat[val1].notna())&(df_AWI_SD_NH_filt[val2].notna())]
len_N = len(df_sat)
corr = df_sat[val1].corr(df_sat[val2])
bias = np.mean(df_sat[val1]-df_sat[val2])
RMSE = rmse(df_sat[val1], df_sat[val2])


ax3.errorbar(x = df['datetime'], y=df['obsSIT'], yerr=df['obsSIT_std'], label="Reference measurements", zorder=10)
ax3.errorbar(x = df['datetime'], y=df['satSIT'], yerr=df['satSIT_std'], label="Satellite measurements")
ax3.format(ylabel='sea ice\n thickness (m)',  ultitle='MOSAIC SIMBA, N = {}, R = {:.2f}, bias = {:.2f} m, RMSE = {:.2f} m'.format(len_N, corr, bias, RMSE), ylim=(-0.1, 4.3))




fig.legend(loc='b')

fig.format(abc='(a)',abcloc='l')

fig.save(r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\representation_error_study_file_sharing\timeseries-buoys.png', dpi=300)

#%%

df = df_IMB_NH_alm
df_IMB_NH_alm['datetime'] = pd.to_datetime(df_IMB_NH_alm['date'],format="%Y-%m-%dT%H:%M:%S") 
start_date = '2010-08-01 00:00:00'
end_date = '2011-06-01 00:00:00'
mask = (df['datetime'] > start_date) & (df['datetime'] < end_date)

df_IMB_NH_alm_filt = df_IMB_NH_alm[mask] 



#%%

fig,ax = pplt.subplots(ncols=1, nrows=1, sharex=True, refwidth=3, refheight=3, sharey=False, proj='npstere')
fig.patch.set_facecolor('white')
resol = '10m'
land = cartopy.feature.NaturalEarthFeature('physical', 'land',
                                           scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])



axs = ax[0]
axs.add_feature(cfeature.LAND, facecolor='lightgrey')
axs.coastlines(resolution=resol, color='k')

QFT_val = 4

fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-SB_AWI-NH-CS2-NH-CCIp-v3p0-rc2.dat'
df_AWI_SD_NH = pd.read_csv(fp, sep = " ")
df_AWI_SD_NH = df_AWI_SD_NH[df_AWI_SD_NH['QFT']<=QFT_val]
df_AWI_SD_NH['datetime'] = pd.to_datetime(df_AWI_SD_NH['date'],format="%Y-%m-%dT%H:%M:%S") 

fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-IMB-CS2-NH-CCIp-v3p0-rc2.dat'
#df_IMB_NH_v2 = pd.read_csv(fp, sep = " ")
df_IMB_NH = pd.read_csv(fp, sep = " ")
df_IMB_NH = df_IMB_NH[df_IMB_NH['QFT']<=QFT_val]

#df_IMB_NH_alm = df_IMB_NH_alm[df_IMB_NH_alm['QFT']<=QFT_val]

fp = r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\buoy_timeseries_study\ESACCIplus-SEAICE-RRDP2+-SIT-MOSAIC-SIMBA-CS2-NH-CCIp-v3p0-rc2.dat'
df_MOSAIC_NH = pd.read_csv(fp, sep = " ")
df_MOSAIC_NH = df_MOSAIC_NH[df_MOSAIC_NH['QFT']<=QFT_val]

df = df_IMB_NH
val = 'obsSIT'
axs.scatter(df[df[val].notna()]['lon'],df[df[val].notna()]['lat'], s=1)
#cb = axs.scatter(df_IMB_NH_alm_filt['lon'],df_IMB_NH_alm_filt['lat'], c=df_IMB_NH_alm_filt['obsSIT'], vmin=0., vmax=4, s=1, extend="max",cmap='blues')
axs.scatter(df_IMB_NH_alm_filt[df_IMB_NH_alm_filt['QFT']<=QFT_val]['lon'],df_IMB_NH_alm_filt[df_IMB_NH_alm_filt['QFT']<=QFT_val]['lat'], c='r', s=1)
axs.format(boundinglat=60)

fig.save(r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\representation_error_study_file_sharing\buoy-spatial-comb_QFT={}_{}_{}.png'.format(QFT_val, 'IMB', val), dpi=300)

#%%

fig,ax = pplt.subplots(ncols=2, nrows=1, sharex=True, refwidth=3, refheight=3, sharey=False, proj='npstere')
fig.patch.set_facecolor('white')
resol = '10m'
land = cartopy.feature.NaturalEarthFeature('physical', 'land',
                                           scale=resol, edgecolor='k', facecolor=cfeature.COLORS['land'])



axs = ax[0]
axs.add_feature(cfeature.LAND, facecolor='lightgrey')
axs.coastlines(resolution=resol, color='k')


axs.scatter(df_IMB_NH_alm_filt['lon'],df_IMB_NH_alm_filt['lat'], c=df_IMB_NH_alm_filt['obsSIT'], vmin=0., vmax=3, s=1,cmap='blues')
axs.format(boundinglat=60)


axs = ax[1]
axs.add_feature(cfeature.LAND, facecolor='lightgrey')
axs.coastlines(resolution=resol, color='k')


cb = axs.scatter(df_IMB_NH_alm_filt[df_IMB_NH_alm_filt['QFT']<=1]['lon'],df_IMB_NH_alm_filt[df_IMB_NH_alm_filt['QFT']<=1]['lat'], c=df_IMB_NH_alm_filt[df_IMB_NH_alm_filt['QFT']<=1]['obsSIT'], vmin=0., vmax=4, s=1, extend="max",cmap='blues')
axs.colorbar(cb, label='sea ice thickness (m)',loc='r')
axs.format(boundinglat=60)
fig.format()

fig.save(r'C:\Users\rmfha\OneDrive - Danmarks Tekniske Universitet\DTU\Projects\CCIp\representation_error_study_file_sharing\buoy-spatial-distr-IMB-check.png', dpi=300)


