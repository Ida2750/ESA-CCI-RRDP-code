# -*- coding: utf-8 -*-
"""
Created on Thu May  8 11:34:07 2025

@author: Ida Olsen
"""
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from scipy.stats import pearsonr

import pandas as pd
import matplotlib.patches as mpatches

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os

# === Config ===
path = "/dmidata/users/ilo/projects/RRDPp/satellite/Final_files"

name = 'OIB-NH'
var = 'SIT'  # Only affects column choice for SID vs non-SID
sat = 'CS2'
HS = 'NH'

valid_variables = ['SIT', 'SD', 'SIF']
days_list = [ 30, 16 ,4] #, 16] #, 30] #, 30]
res_list = [25000, 15000, 5000]  # assumed resolution in meters
filt = True

ObsNames_SID = [
    'date','lat','lon','obsSID','satSID',
    "obsSID_std","obsSID_ln","obsSID_unc", "QFT", "QFS", "QFG",
    "satSID_std","satSID_ln","satSID_unc",
]

ObsNames2 = [
    'date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT','satSIF',
    "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
    "obsFRB_std","obsFRB_ln","obsFRB_unc", "QFT", "QFS", "QFG",
    "satSD_std","satSD_ln","satSD_unc","satSIT_std","satSIT_ln","satSIT_unc",
    "satFRB_std","satFRB_ln","satFRB_unc", "index"
]

ObsNames_new = [
    'date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT','satSIF',
    "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
    "obsFRB_std","obsFRB_ln","obsFRB_unc", "QFT", "QFS", "QFG",
    "satSD_std","satSD_ln","satSD_unc","satSIT_std","satSIT_ln","satSIT_unc",
    "satFRB_std","satFRB_ln","satFRB_unc", "index"
]
# === Load Datasets ===
datasets = []
dataset_labels = []

# define gridcells to include
ref_grid = f'ESACCIplus-SEAICE-RRDP2+-{var}-{name}-{sat}-{HS}-4-5000-CCIp-v3p0-rc2.dat'
filepath = os.path.join(path, ref_grid)
data_ref = np.genfromtxt(filepath, dtype=None, skip_header=1, names=ObsNames_new, encoding='utf-8')

index_ref = data_ref['index']

for days in days_list:
    for res in res_list:
        filename = f'ESACCIplus-SEAICE-RRDP2+-{var}-{name}-{sat}-{HS}-{days}-{res}-CCIp-v3p0-rc2.dat'
        filepath = os.path.join(path, filename)
        if not os.path.exists(filepath):
            print(f"File not found: {filename} — skipping.")
            continue

        names = ObsNames_SID if var == 'SID' else ObsNames2
        data = np.genfromtxt(filepath, dtype=None, skip_header=1, names=names, encoding='utf-8')
        data['obsSIT'][data['obsSIT']>10] = np.nan
        mask = np.isin(data['index'], index_ref)
        if filt:
            data = data[mask]
        if var == 'SID':
            data['satSID'][data['satSID'] < 0] = np.nan
        else:
            data['satSIT'][data['satSIT'] < 0] = np.nan
            #data['satSIF'] = data['satSIF'] + data['satSD']
            data['obsSIF'] = data['obsSIF'] - data['obsSD']
            data['satSIF'][data['satSIF'] < 0] = np.nan
            data['obsSIF'][data['obsSIF'] < 0] = np.nan 

        datasets.append(data)
        dataset_labels.append(f"{days} days / {res//1000} km")

# === Plot Loop ===
colors = plt.cm.viridis(np.linspace(0, 1, len(datasets)))  # Unique color per dataset

for plot_var in valid_variables:
    obs_col = f'obs{plot_var}'
    sat_col = f'sat{plot_var}'

    
    #plt.figure(figsize=(12, 9))
    fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 12), constrained_layout=True)
    ax = ax.flatten()
    for i, data in enumerate(datasets):
        if obs_col not in data.dtype.names or sat_col not in data.dtype.names:
            print(f"Skipping {plot_var} for dataset {i} — missing column.")
            continue

        obs = data[obs_col]
        sat = data[sat_col]

        mask = ~np.isnan(obs) & ~np.isnan(sat)
        obs_clean = obs[mask]
        sat_clean = sat[mask]
        if plot_var=='SIF':
            plot_var = 'FRB'
        obs_len_clean = data[f'obs{plot_var}_ln'][mask]
        sat_len_clean = data[f'sat{plot_var}_ln'][mask]

        if len(obs_clean) == 0 or len(sat_clean) == 0:
            print(f"No valid data for {plot_var} in dataset {i}")
            continue

        median_diff = np.median(sat_clean - obs_clean)
        corr, _ = pearsonr(obs_clean, sat_clean)

        footprint_diameter_CS2 = (1650+300)/2
        #footprint_diameter_OIB = 40/2

        area_OIB = 40**2*obs_len_clean #3.14*(footprint_diameter_OIB/2)**2 *obs_len_clean
        area_CS2 = 3.14*(footprint_diameter_CS2/2)**2 *sat_len_clean
        # source: https://earth.esa.int/eogateway/documents/20142/37627/CryoSat-Footprints-ESA-Aresys.pdf

        percentage_total_area_OIB = np.median(area_OIB)/(3.14*(25000/2)**2)*100
        percentage_total_area_CS2 = np.median(area_CS2)/(3.14*(25000/2)**2)*100

        # , area % {percentage_total_area_CS2:.2f}
        # , area % {percentage_total_area_OIB:.2f}

        # plt.scatter(obs_clean, sat_clean, s=10, alpha=0.5, color=colors[i],
        #             label=f'{dataset_labels[i]}\nΔmedian: {median_diff:.2f}, R: {corr:.2f}, Median len SAT: {np.mean(sat_len_clean):.0f}, Median len OBS: {np.mean(obs_len_clean):.0f}')

        h = ax[i].hist2d(obs_clean, sat_clean, bins=(30,30), cmap=plt.cm.Blues)

        # # make a fake artist for legend
        # patch = mpatches.Patch(color='blue', alpha=0.5,
        #     label=f'{dataset_labels[i]}\nΔmedian: {median_diff:.2f}, R: {corr:.2f}, '
        #         f'Median len SAT: {np.mean(sat_len_clean):.0f}, Median len OBS: {np.mean(obs_len_clean):.0f}')
        # leg = ax[i].legend(handles=[patch], fontsize=11, loc='upper left')
        # Reference 1:1
        all_obs_vals = np.concatenate([
            data[obs_col][~np.isnan(data[obs_col])] for data in datasets if obs_col in data.dtype.names
        ])
        if len(all_obs_vals) > 0:
            min_val, max_val = np.nanmin(all_obs_vals), np.nanmax(all_obs_vals)
            ax[i].plot([min_val, max_val], [min_val, max_val], 'k--', lw=1)

        from matplotlib.lines import Line2D

        custom_legend = [Line2D([0], [0], color='w', markerfacecolor='blue', marker='s',
                                markersize=10, label=f'{dataset_labels[i]}\nΔmedian: {median_diff:.2f}, \nR: {corr:.2f}, \nMedian len SAT: {np.mean(sat_len_clean):.0f}, \nMedian len OBS: {np.mean(obs_len_clean):.0f}')]
        ax[i].legend(handles=custom_legend, fontsize=11, loc='upper left')
        ax[i].grid(True)

        #plt.show()
        ax[i].set_xlim(0,10)
        ax[i].set_ylim(0,10)


    fig.supxlabel(f'Observed {plot_var}', fontsize=14)
    fig.supylabel(f'Satellite {plot_var}', fontsize=14)
    fig.suptitle(f'Observed vs Satellite {plot_var} (All Day/Res Combinations): Number of gridcells: {len(obs_clean)}', fontsize=16)
    plt.tight_layout()
    plt.savefig(f'scatter_{name}_{plot_var}_all_filter_{filt}.png') 
