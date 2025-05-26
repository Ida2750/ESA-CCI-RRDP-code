import os
import xarray as xr
import glob
import numpy as np
from numpy.lib import recfunctions as rfn


def compute_qfs(obs_val, q1, q2, q3):
    #print(obs_val)
    if obs_val > q3:
        return 0
    elif obs_val > q2:
        return 1
    elif obs_val > q1:
        return 2
    else:
        return 3

def remove_qfg(data):
    if 'QFG' in data.dtype.names:
        data = data[[name for name in data.dtype.names if name != 'QFG']]
    return data

def update_qfs_array(data, var_name):
    if var_name not in data.dtype.names or 'QFS' not in data.dtype.names:
        return data
    values = data[var_name]
    #print(values)
    q1, q2, q3 = np.percentile(values[~np.isnan(values)], [25, 50, 75])
    #print(q1, q2, q3)
    qfs = np.array([compute_qfs(val, q1, q2, q3) for val in values])
    data['QFS'] = qfs
    return data

def update_qfs_xr(ds, var_name):
    if var_name not in ds or 'QFS' not in ds:
        return ds
    #print(ds)
    values = ds[var_name].values
    q1, q2, q3 = np.percentile(values[~np.isnan(values)], [25, 50, 75])
    qfs = np.array([compute_qfs(val, q1, q2, q3) for val in values])
    ds['QFS'] = xr.DataArray(qfs, coords=ds[var_name].coords, dims=ds[var_name].dims)
    return ds

# ---------------------- Process .dat and .nc files in FINAL ------------------------
basepath = '/dmidata/users/ilo/projects/RRDPp/FINAL'
output_basepath = os.path.join(basepath, 'combined_final')
os.makedirs(output_basepath, exist_ok=True)

OBSID = os.listdir(basepath) + ['ASPeCt', 'AWI_ULS']
for obsid in OBSID:
    if obsid == 'combined_final':
        continue  # skip output folder
    files_dat = glob.glob(f'{basepath}/{obsid}/final/*.dat') + glob.glob(f'{basepath}/Antarctic/{obsid}/final/*.dat')
    files_nc = glob.glob(f'{basepath}/{obsid}/final/*.nc') + glob.glob(f'{basepath}/Antarctic/{obsid}/final/*.nc')

    for file in files_dat:
        data = np.genfromtxt(file, names=True, dtype=None, encoding=None)

        data = remove_qfg(data)


        if '-IMB' in file or 'SB-AWI' in file or 'ASPeCt' in file or 'ASSIST' in file or 'SB_AWI' in file:
            # update QFS to 3 if it is not already
            if 'QFS' in data.dtype.names:
                data['QFS'][data['QFS'] != 3] = 3

        if '-IMB' in file or 'SCICEX' in file:
            years = np.array([int(str(d)[:4]) if str(d)[:4].isdigit() else np.nan for d in data['date']])
            if '-IMB' in file:
                data = data[years != 2017]
            if 'SCICEX' in file:
                data = data[years != 2014]

        if 'BGEP' in file or ('Nansen' in file and 'SID' in file):
            data['uncflag'][data['uncflag'] == 1] = 2

        if 'SCICEX' in file or 'AEM-AWI' in file:
            if 'ppflag' in data.dtype.names:
                data['ppflag'][data['ppflag'] == 1] = 2
        
        if 'OIB' in file:
            var = 'FRBln'
        elif 'SCICEX' in file:
            var = 'SIDln'
        elif 'AEM-AWI' in file or 'HEM' in file:
            var = 'SITln'
        else:
            var = None  # No QFS update for these

        if var is not None:
            data = update_qfs_array(data, var)
        
        outfile = os.path.join(output_basepath, os.path.basename(file))
        np.savetxt(outfile, data, header=' '.join(data.dtype.names), fmt='%s', comments='')

    for file in files_nc:
        ds = xr.open_dataset(file)

        # Remove QFG if it exists
        if 'QFG' in ds:
            ds = ds.drop_vars('QFG')

        # Force QFS to 3 for specific datasets
        if '-IMB' in file or 'SB-AWI' in file or 'ASPeCt' in file or 'ASSIST' in file or 'SB_AWI' in file:
            if 'QFS' in ds:
                ds['QFS'] = ds['QFS'].where(ds['QFS'] == 3, other=3)

        # Filter out faulty years for IMB (2017) and SCICEX (2014)
        if '-IMB' in file or 'SCICEX' in file:
            #if 'date' in ds:
            # Ensure datetime64 dtype
            if not np.issubdtype(ds['time'].dtype, np.datetime64):
                ds['time'] = xr.decode_cf(ds).date

            years = ds['time'].dt.year

            if '-IMB' in file:
                ds = ds.where(years != 2017, drop=True)
            if 'SCICEX' in file:
                ds = ds.where(years != 2014, drop=True)
            
            print(ds)

        # Update unc from 1 to 2 where needed
        if 'BGEP' in file or ('Nansen' in file and 'SID' in file):
            ds['uncflag'] = ds['uncflag'].where(ds['uncflag'] != 1, other=2)

        # Update ppflag from 1 to 2 where applicable
        if 'SCICEX' in file or 'AEM-AWI' in file:
            if 'ppflag' in ds:
                ds['ppflag'] = ds['ppflag'].where(ds['ppflag'] != 1, other=2)

        if 'QFG' in ds:
            ds = ds.drop_vars('QFG')
        if 'OIB' in file:
            var = 'FRBln'
        elif 'SCICEX' in file:
            var = 'SIDln'
        elif 'AEM-AWI' in file or 'HEM' in file:
            var = 'SITln'
        else:
            var = None  # No QFS update for these

        if var is not None:
            ds = update_qfs_xr(ds, var)

        outfile = os.path.join(output_basepath, os.path.basename(file))
        ds.to_netcdf(outfile)

#%%
# ---------------------- Process collocated satellite .dat files ------------------------
sat_path = '/dmidata/users/ilo/projects/RRDPp/satellite/Final_files'
sat_outpath = os.path.join(sat_path, 'combined_final')
os.makedirs(sat_outpath, exist_ok=True)

files_dat = glob.glob(f'{sat_path}/*.dat')

ObsNames = [
    'date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT','satSIF',
    "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
    "obsFRB_std","obsFRB_ln","obsFRB_unc", "QFT", "QFS", "QFG",
    "satSD_std","satSD_ln","satSD_unc","satSIT_std","satSIT_ln","satSIT_unc",
    "satFRB_std","satFRB_ln","satFRB_unc", "index"
]

ObsNames_SID = [
    'date','lat','lon','obsSID','satSID',
    "obsSID_std","obsSID_ln","obsSID_unc", "QFT", "QFS", "QFG",
    "satSID_std","satSID_ln","satSID_unc",
]

for file in files_dat:
    try:
        if 'SID' in file:
            data = np.genfromtxt(file, dtype=None, skip_header=1, names=ObsNames_SID, encoding='utf-8')
        else:
            data = np.genfromtxt(file, dtype=None, skip_header=1, names=ObsNames, encoding='utf-8')

        # Skip if file is empty or has no rows
        if data.size == 0:
            print(f"Skipping empty file: {file}")
            continue
    except Exception as e:
        print(f"Error reading file {file}: {e}")
        continue

    # Force QFS = 3 for specific sources
    if '-IMB' in file or 'SB-AWI' in file or 'ASPeCt' in file or 'ASSIST' in file or 'SB_AWI' in file:
        if 'QFS' in data.dtype.names:
            data['QFS'][data['QFS'] != 3] = 3

    # Remove faulty years
    if '-IMB' in file or 'SCICEX' in file:
        if 'date' in data.dtype.names:
            years = np.array([int(str(d)[:4]) if str(d)[:4].isdigit() else np.nan for d in data['date']])
            if '-IMB' in file:
                data = data[years != 2017]
            if 'SCICEX' in file:
                data = data[years != 2014]
                
    # Always remove QFG
    data = remove_qfg(data)

    # Conditionally update QFS
    if 'OIB' in file:
        var = 'obsFRB_ln'
    elif 'SCICEX' in file:
        var = 'obsSID_ln'
    elif 'AEM-AWI' in file or 'HEM' in file:
        var = 'obsSIT_ln'
    else:
        var = None

    if var is not None:
        data = update_qfs_array(data, var)

    outfile = os.path.join(sat_outpath, os.path.basename(file))
    np.savetxt(outfile, data, header=' '.join(data.dtype.names), fmt='%s', comments='')
