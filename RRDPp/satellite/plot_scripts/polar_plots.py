# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:00:19 2020

@author: Olsen
"""
# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2023-07-25'

# -- Built-in modules -- #
import os
# -- Third-part modules -- #
from matplotlib.lines import Line2D
import matplotlib as mpl
import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
# -- Proprietary modules -- #


def plot(lat,lon,z, title, ylabel='SIT [m]', clim=False, savename=False, NP=True, saveloc=False, s=8):

    plt.figure(figsize=(7, 7))
    if NP==True:
        ax = plt.axes(projection=ccrs.NorthPolarStereo())
        ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
    else: # we are plotting for the South Pole
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180,180,-90,-55],ccrs.PlateCarree())  
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.gridlines(draw_labels=False)

    plot=plt.scatter(lon, lat,
              s=s, c=z,alpha=1,cmap='bwr',
              transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    cmap = mpl.cm.cool
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    cb.ax.set_ylabel(ylabel);
    cb.ax.tick_params(labelsize=12) 
       
    if clim==False:
        plt.clim(np.nanmin(z),np.nanmax(z))
    else:
        plt.clim(clim)
    if savename:
        if saveloc:
            save = os.path.join(saveloc,savename)
            print(save)
            plt.savefig(os.path.join(saveloc,savename), bbox_inches='tight')
            plt.show()
        else:
            plt.savefig(savename, bbox_inches='tight')
    
        