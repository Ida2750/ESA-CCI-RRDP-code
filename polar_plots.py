# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 11:00:19 2020

@author: Olsen
"""
def plot(lat,lon,z, cmin, cmax, extent=[-180, 180, 65, 90]):
    import cartopy.crs as ccrs
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.path as mpath
    import cartopy.feature as cfeature
    import  matplotlib
    

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(121, projection=ccrs.NorthPolarStereo())
    # make geographical background and axes
    ax.coastlines()
    ax.add_feature(cfeature.OCEAN)        
    ax.add_feature(cfeature.LAND)
    ax.set_extent(extent, ccrs.PlateCarree())
    ax.gridlines()
    # make circular plot
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)
    # plot data
    plot = ax.scatter(lon, lat, c=z, s=1, alpha=0.5,
             transform=ccrs.PlateCarree(),)
    # create colorbar
    cbar = plt.colorbar(plot)
    plot.set_clim(cmin, cmax)
    cbar.ax.tick_params(labelsize=11) 
