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

# df_ref = pd.DataFrame(data_ref)
# #df_ref = df_ref.dropna(subset=['obsSIF', 'satSIF'])
# df_ref = df_ref.drop_duplicates(subset=df_ref.columns[:6])
# columns_to_check = df_ref.columns[:6]
# keys1 = np.core.defchararray.add(
#     np.char.add(data_ref['date'].astype(str), data_ref['lat'].astype(str)),
#     np.char.add(data_ref['lon'].astype(str), data_ref['obsSIF'].astype(str))
# )
# unique1 = np.unique(data_ref)

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
        # df = pd.DataFrame(data)
        # #df = df.dropna(subset=['obsSIF', 'satSIF'])
        # df = df.drop_duplicates(subset=df_ref.columns[:6])
        # matching_rows = df[df[columns_to_check].apply(tuple, axis=1).isin(
        # df_ref[columns_to_check].apply(tuple, axis=1)
        # )]
        #dtype = data.dtype
        #data = np.unique(data, axis=0)

        # keys2 = np.core.defchararray.add(
        #     np.char.add(data['date'].astype(str), data['lat'].astype(str)),
        #     np.char.add(data['lon'].astype(str), data['obsSIF'].astype(str))
        # )

        # keys1 = set(zip(data_ref['date'], data_ref['lat'], data_ref['lon'], data_ref['obsSIF']))
        # mask = np.array([ (t, la, lo, obs) in keys1 for t, la, lo, obs in zip(data['date'], data['lat'], data['lon'], data['obsSIF']) ])
        # matched_data2 = data[mask]

        # Step 1: Remove duplicates from both
        # unique2 = np.unique(data)

        # # Step 2: Convert to sets of tuples for easy comparison
        # set1 = set(map(tuple, unique1))
        # set2 = set(map(tuple, unique2))

        # # Step 3: Intersection to find matching records
        # matched_set = set2.intersection(set1)

        # # Step 4: Convert matched result back to structured array
        # matched_data2 = np.array(list(matched_set), dtype=dtype)

        # Result
        

        # # Apply mask
        # if filt:
        #     #data = data[mask]
        #     #data = np.copy(matched_data2)
        #     data = matching_rows.to_records(index=False)
        #     #print(data)
        #     #print(len(matching_rows))  # Now this should never exceed len(ref_coords_time)
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

    plt.figure(figsize=(12, 9))

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

        plt.scatter(obs_clean, sat_clean, s=10, alpha=0.5, color=colors[i],
                    label=f'{dataset_labels[i]}\nΔmedian: {median_diff:.2f}, R: {corr:.2f}, Median len SAT: {np.mean(sat_len_clean):.0f}, Median len OBS: {np.mean(obs_len_clean):.0f}')

    # Reference 1:1
    all_obs_vals = np.concatenate([
        data[obs_col][~np.isnan(data[obs_col])] for data in datasets if obs_col in data.dtype.names
    ])
    if len(all_obs_vals) > 0:
        min_val, max_val = np.nanmin(all_obs_vals), np.nanmax(all_obs_vals)
        plt.plot([min_val, max_val], [min_val, max_val], 'k--', lw=1)

    plt.xlabel(f'Observed {plot_var}', fontsize=14)
    plt.ylabel(f'Satellite {plot_var}', fontsize=14)
    plt.title(f'Observed vs Satellite {plot_var} (All Day/Res Combinations): Number of gridcells: {len(obs_clean)}', fontsize=16)
    plt.legend(fontsize=11, loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    #plt.show()
    plt.xlim(0,10)
    plt.ylim(0,10)
    plt.savefig(f'scatter_{name}_{plot_var}_all_filter_{filt}.png')  # Optional saving


# # File path and filenames
# path = "/dmidata/users/ilo/projects/RRDPp/satellite/Final_files"
# # file3 = 'ESACCIplus-SEAICE-RRDP2+-SIT-OIB-CS2-CCIp-v3p0-rc2.DAT'
# # file2 = 'ESACCIplus-SEAICE-RRDP2+-SIT-OIB-CS2-CCIp-v3p0-rc2.DAT'
# # file1 = 'ESACCIplus-SEAICE-RRDP2+-SIT-OIB-CS2-CCIp-v3p0-rc2.DAT'

# name = 'BGEP'
# var = 'SID'
# sat = 'CS2'
# HS = 'NH'

# file3 = f'ESACCIplus-SEAICE-RRDP2+-{var}-{name}-{sat}-{HS}-CCIp-v3p0-rc2.dat'
# file2 = f'ESACCIplus-SEAICE-RRDP2+-{var}-{name}-{sat}-{HS}-CCIp-v3p0-rc2.dat'
# file1 = f'ESACCIplus-SEAICE-RRDP2+-{var}-{name}-{sat}-{HS}-CCIp-v3p0-rc2.dat'

# valid_variables = ['SID']
# # path2 = "C:/Users/Ida Olsen/Documents/work/ESA-CCI-RRDP-code/RRDPp/FINAL/OIB/final"
# # file = "ESACCIplus-SEAICE-RRDP2+-SIT-OIB-NH.nc"
# # flags_data = xr.open_dataset(f'{path2}/{file}')
# # QFT = flags_data['QFT'].to_numpy()
# # QFS = flags_data['QFS'].to_numpy()
# # QFG = flags_data['QFG'].to_numpy()

# # Column names

# ObsNames_SID = [
#     'date','lat','lon','obsSID','satSID',
#     "obsSID_std","obsSID_ln","obsSID_unc", "QFT", "QFS", "QFG","satSID_std","satSID_ln","satSID_unc",
# ]

# ObsNames2 = [
#     'date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT','satSIF',
#     "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
#     "obsFRB_std","obsFRB_ln","obsFRB_unc", "QFT", "QFS", "QFG","satSD_std","satSD_ln","satSD_unc",
#     "satSIT_std","satSIT_ln","satSIT_unc","satFRB_std","satFRB_ln","satFRB_unc"
# ]

# #footprint_OIB =
# #footprint_CS2 = 

# # Load data
# if var=='SID':
#     data1 = np.genfromtxt(f'{path}/{file1}', dtype=None, skip_header=1, names=ObsNames_SID, encoding='utf-8')
#     data2 = np.genfromtxt(f'{path}/{file2}', dtype=None, skip_header=1, names=ObsNames_SID, encoding='utf-8')
#     data3 = np.genfromtxt(f'{path}/{file3}', dtype=None, skip_header=1, names=ObsNames_SID, encoding='utf-8')

#     data1['satSID'][data1['satSID']<0] = np.nan
#     data2['satSID'][data2['satSID']<0] = np.nan
#     data3['satSID'][data3['satSID']<0] = np.nan
# else:
#     data1 = np.genfromtxt(f'{path}/{file1}', dtype=None, skip_header=1, names=ObsNames2, encoding='utf-8')
#     data2 = np.genfromtxt(f'{path}/{file2}', dtype=None, skip_header=1, names=ObsNames2, encoding='utf-8')
#     data3 = np.genfromtxt(f'{path}/{file3}', dtype=None, skip_header=1, names=ObsNames2, encoding='utf-8')

#     data1['satSIT'][data1['satSIT']<0] = np.nan
#     data2['satSIT'][data2['satSIT']<0] = np.nan
#     data3['satSIT'][data3['satSIT']<0] = np.nan

#     data1['satSIF'] = data1['satSIF'] + data1['satSD'] 
#     data2['satSIF'] = data2['satSIF'] + data2['satSD']
#     data3['satSIF'] = data3['satSIF'] + data3['satSD']

#     data1['satSIF'][data1['satSIF']<0] = np.nan
#     data2['satSIF'][data2['satSIF']<0] = np.nan
#     data3['satSIF'][data3['satSIF']<0] = np.nan

# # List of datasets and names
# datasets = [data1, data2, data3]
# dataset_names = ['30 days', '14 days', '6 days']

# # Loop over each dataset and each variable
# for i, data in enumerate(datasets):
#     print(f"\n=== Statistics for {dataset_names[i]} ===")
#     for var in valid_variables:
#         obs = data[f'obs{var}']
#         sat = data[f'sat{var}']
        
#         # Compute metrics
#         satAvg = np.nanmedian(sat)
#         satStd = np.nanstd(sat)
#         obsAvg = np.nanmedian(obs)
#         obsStd = np.nanstd(obs)
#         bias = np.nanmedian(obs - sat)

#         # Print results
#         print(f"\nVariable: {var}")
#         print(f"  Satellite Median: {satAvg:.3f}")
#         print(f"  Satellite Std Dev: {satStd:.3f}")
#         print(f"  Observed Median: {obsAvg:.3f}")
#         print(f"  Observed Std Dev: {obsStd:.3f}")
#         print(f"  Bias (Obs - Sat): {bias:.3f}")

# # Plot SIT
# for vv in valid_variables:
#     plt.figure(figsize=(10,10))
#     plt.title(f"Observed vs Satellite {vv}")
#     index = np.isfinite(data1[f'obs{vv}']) & np.isfinite(data1[f'sat{vv}']) 
#     corr = np.round(pearsonr(data1[f'obs{vv}'][index], data1[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data1[f'obs{vv}'], data1[f'sat{vv}'], label=f'30 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data2[f'obs{vv}']) & np.isfinite(data2[f'sat{vv}']) 
#     corr = np.round(pearsonr(data2[f'obs{vv}'][index], data2[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data2[f'obs{vv}'], data2[f'sat{vv}'], label=f'14 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data3[f'obs{vv}']) & np.isfinite(data3[f'sat{vv}']) 
#     corr = np.round(pearsonr(data3[f'obs{vv}'][index], data3[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data3[f'obs{vv}'], data3[f'sat{vv}'], label=f'6 days R:{corr}', alpha=0.6, s=1)
#     plt.xlabel(f"Observed {vv}")
#     plt.ylabel(f"Satellite {vv}")
#     plt.legend()
#     plt.grid(True)
#     plt.show()

# #%%
# ##############################################
# # Global limits

# data1['obsSIT'][data1['obsSIT']>8] = np.nan
# data1['satSIT'][data1['satSIT']>8] = np.nan
# data2['obsSIT'][data2['obsSIT']>8] = np.nan
# data2['satSIT'][data2['satSIT']>8] = np.nan
# data3['obsSIT'][data3['obsSIT']>8] = np.nan
# data3['satSIT'][data3['satSIT']>8] = np.nan

# data1['obsSIF'][data1['obsSIF']>2] = np.nan
# data1['satSIF'][data1['satSIF']>2] = np.nan
# data2['obsSIF'][data2['obsSIF']>2] = np.nan
# data2['satSIF'][data2['satSIF']>2] = np.nan
# data3['obsSIF'][data3['obsSIF']>2] = np.nan
# data3['satSIF'][data3['satSIF']>2] = np.nan

# data1['obsSD'][data1['obsSD']>2] = np.nan
# data1['satSD'][data1['satSD']>2] = np.nan
# data2['obsSD'][data2['obsSD']>2] = np.nan
# data2['satSD'][data2['satSD']>2] = np.nan
# data3['obsSD'][data3['obsSD']>2] = np.nan
# data3['satSD'][data3['satSD']>2] = np.nan


# # Plot SIT
# for vv in valid_variables:
#     plt.figure(figsize=(10,10))
#     plt.title(f"Observed vs Satellite {vv}")
#     index = np.isfinite(data1[f'obs{vv}']) & np.isfinite(data1[f'sat{vv}']) 
#     corr = np.round(pearsonr(data1[f'obs{vv}'][index], data1[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data1[f'obs{vv}'], data1[f'sat{vv}'], label=f'30 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data2[f'obs{vv}']) & np.isfinite(data2[f'sat{vv}']) 
#     corr = np.round(pearsonr(data2[f'obs{vv}'][index], data2[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data2[f'obs{vv}'], data2[f'sat{vv}'], label=f'14 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data3[f'obs{vv}']) & np.isfinite(data3[f'sat{vv}']) 
#     corr = np.round(pearsonr(data3[f'obs{vv}'][index], data3[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data3[f'obs{vv}'], data3[f'sat{vv}'], label=f'6 days R:{corr}', alpha=0.6, s=1)
#     plt.xlabel(f"Observed {vv}")
#     plt.ylabel(f"Satellite {vv}")
#     plt.legend()
#     plt.grid(True)
#     plt.show()



# #%%
# ##############################################
# # number of observations # has more or less no impact
# obsnum1 = (data1['obsSIT_ln']<20) | (data1['satSIT_ln']<20)
# obsnum2 = (data2['obsSIT_ln']<20) | (data2['satSIT_ln']<20)
# obsnum3 = (data3['obsSIT_ln']<20) | (data3['satSIT_ln']<20)

# data1['obsSIT'][obsnum1] = np.nan
# data1['satSIT'][obsnum1] = np.nan
# data2['obsSIT'][obsnum2] = np.nan
# data2['satSIT'][obsnum2] = np.nan
# data3['obsSIT'][obsnum3] = np.nan
# data3['satSIT'][obsnum3] = np.nan

# obsnum1 = (data1['obsFRB_ln']<20) | (data1['satFRB_ln']<20)
# obsnum2 = (data2['obsFRB_ln']<20) | (data2['satFRB_ln']<20)
# obsnum3 = (data3['obsFRB_ln']<20) | (data3['satFRB_ln']<20)

# data1['obsSIF'][obsnum1] = np.nan
# data1['satSIF'][obsnum1] = np.nan
# data2['obsSIF'][obsnum2] = np.nan
# data2['satSIF'][obsnum2] = np.nan
# data3['obsSIF'][obsnum3] = np.nan
# data3['satSIF'][obsnum3] = np.nan

# obsnum1 = (data1['obsSD_ln']<20) | (data1['satSD_ln']<20)
# obsnum2 = (data2['obsSD_ln']<20) | (data2['satSD_ln']<20)
# obsnum3 = (data3['obsSD_ln']<20) | (data3['satSD_ln']<20)

# data1['obsSD'][obsnum1] = np.nan
# data1['satSD'][obsnum1] = np.nan
# data2['obsSD'][obsnum2] = np.nan
# data2['satSD'][obsnum2] = np.nan
# data3['obsSD'][obsnum3] = np.nan
# data3['satSD'][obsnum3] = np.nan

# # Plot SIT
# for vv in valid_variables:
#     plt.figure(figsize=(10,10))
#     plt.title(f"Observed vs Satellite {vv}")
#     index = np.isfinite(data1[f'obs{vv}']) & np.isfinite(data1[f'sat{vv}']) 
#     corr = np.round(pearsonr(data1[f'obs{vv}'][index], data1[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data1[f'obs{vv}'], data1[f'sat{vv}'], label=f'30 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data2[f'obs{vv}']) & np.isfinite(data2[f'sat{vv}']) 
#     corr = np.round(pearsonr(data2[f'obs{vv}'][index], data2[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data2[f'obs{vv}'], data2[f'sat{vv}'], label=f'14 days R:{corr}', alpha=0.6, s=1)
#     index = np.isfinite(data3[f'obs{vv}']) & np.isfinite(data3[f'sat{vv}']) 
#     corr = np.round(pearsonr(data3[f'obs{vv}'][index], data3[f'sat{vv}'][index]), 2)[0]
#     plt.scatter(data3[f'obs{vv}'], data3[f'sat{vv}'], label=f'6 days R:{corr}', alpha=0.6, s=1)
#     plt.xlabel(f"Observed {vv}")
#     plt.ylabel(f"Satellite {vv}")
#     plt.legend()
#     plt.grid(True)
#     plt.show()
