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
              s=s, c=z, alpha=1,cmap='bwr',
              transform=ccrs.PlateCarree(),)
    plt.title(title,fontsize=16,fontweight='bold')
    
    cmap = mpl.cm.cool
    ax_cb = plt.axes([0.92, 0.25, 0.015, 0.5])
    #cb = plt.colorbar(plot, cax=ax_cb, orientation='vertical')
    #cb.ax.set_ylabel(ylabel);
    #cb.ax.tick_params(labelsize=12) 
       
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


def plot_all(lat, lon, z, title, c, ylabel='SIT [m]', clim=False,
             savename=False, NP=True, saveloc=False, s=18,
             fig=None, ax=None, sat=None):
    
    if fig is None or ax is None:
        fig = plt.figure(figsize=(35, 35))
        if NP:
            ax = plt.axes(projection=ccrs.NorthPolarStereo())
            ax.set_extent([-180, 180, 65, 90], ccrs.PlateCarree())
        else:
            ax = plt.axes(projection=ccrs.SouthPolarStereo())
            ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())

    ax.coastlines()
    ax.add_feature(cfeature.LAND)
    
    gl = ax.gridlines(draw_labels=False, color='grey')
    gl.xlabels_top = True
    gl.ylabels_left = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 48, 'color': 'black', 'weight': 'bold'}

    legend_elements = []
    for la, lo, label, cc in zip(np.flip(lat), np.flip(lon), np.flip(ylabel), np.flip(c)):
        # Size override
        if label == 'SCICEX':
            s = 200
        elif label == 'NPI-FS':
            s = 500

        ax.scatter(lo, la, s=s, c=cc, alpha=1, label=label,
                   transform=ccrs.PlateCarree())

        legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                      label=label, markerfacecolor=cc, markersize=40))

    # Title on the side
    ax.text(-0.03, 0.5, title, va='bottom', ha='center', rotation='vertical',
            rotation_mode='anchor', fontsize=60, transform=ax.transAxes)

    # Legend
    if any(label == 'SCICEX' for label in ylabel):
        ax.legend(handles=legend_elements, fontsize=48)
    else:
        ax.legend(markerscale=4, fontsize=48)

    # Subplot label in bottom-right
    subplot_label = None
    if sat == 'CS2':
        subplot_label = '(a)'
    elif sat == 'ENV':
        subplot_label = '(d)'

    if subplot_label:
        ax.text(0.90, 0.05, subplot_label, transform=ax.transAxes,
                fontsize=80, fontweight='bold', ha='right', va='bottom')

    # Save if requested
    if savename:
        if saveloc:
            full_path = os.path.join(saveloc, savename)
            print(full_path)
            plt.savefig(full_path, bbox_inches='tight')
        else:
            plt.savefig(savename, bbox_inches='tight')
        #plt.show()
        plt.close()

    return fig, ax

# def plot_all(lat,lon,z, title, c, ylabel='SIT [m]', clim=False, savename=False, NP=True, saveloc=False, s=18):
    
#     plt.figure(figsize=(30, 30))
#     if NP==True:
#         ax = plt.axes(projection=ccrs.NorthPolarStereo())
#         ax.set_extent([-180,180,65,90],ccrs.PlateCarree())
#     else: # we are plotting for the South Pole
#         ax = plt.axes(projection=ccrs.SouthPolarStereo())
#         ax.set_extent([-180,180,-90,-55],ccrs.PlateCarree())  
#     ax.coastlines()
#     #ax.add_feature(cfeature.OCEAN)        
#     ax.add_feature(cfeature.LAND)
#     gl = ax.gridlines(draw_labels=False, color='grey')

#     gl.xlabels_top = True
#     gl.ylabels_left = False
#     gl.ylabels_right = False
#     # gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
#     gl.xformatter = LONGITUDE_FORMATTER
#     # gl.yformatter = LATITUDE_FORMATTER
#     #gl.xlabel_style = {'size': 15, 'color': 'gray'}
#     gl.xlabel_style = {'size': 48, 'color': 'black', 'weight': 'bold'}
    
#     legend_elements = []    
#     for la, lo, label, cc in zip(np.flip(lat), np.flip(lon), np.flip(ylabel), np.flip(c)):
#         #s=22
#         if label=='SCICEX':
#             s=200
#         elif label=='NPI-FS':
#             s=500

#         ax.scatter(lo, la,
#                   s=s, c=cc,alpha=1, label=label,
#                   transform=ccrs.PlateCarree(),)
        
#         legend_elements.append(Line2D([0], [0], marker='o', color='w', label=label, markerfacecolor=cc, markersize=40))
        
#     #ax.set_ylabel(title,fontsize=16)
#     ax.text(-0.03, 0.50, title, va='bottom', ha='center',
#         rotation='vertical', rotation_mode='anchor', fontsize=60,
#         transform=ax.transAxes)
    
        
#     if any(ylabel=='SCICEX'):
#         legend = plt.legend(handles=legend_elements, fontsize=48)
#     else:
#         legend = plt.legend(markerscale=4, fontsize=48)
    

#     if savename:
#         if saveloc:
#             save = os.path.join(saveloc,savename)
#             print(save)
#             plt.savefig(os.path.join(saveloc,savename), bbox_inches='tight')
#             plt.show()
#         else:
#             plt.savefig(savename, bbox_inches='tight')
#             plt.show()
    
        