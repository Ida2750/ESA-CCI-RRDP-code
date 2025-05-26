import glob
import os
import numpy as np

basepath = '/dmidata/users/ilo/projects/RRDPp/FINAL'
OBSID = os.listdir(basepath) + os.listdir(f'{basepath}/Antarctic')
for obsid in OBSID:
    if obsid == 'combined_final':
        continue  # skip output folder

    #/dmidata/users/ilo/projects/RRDPp/FINAL/Antarctic/AWI_ULS/final
    files_dat =  glob.glob(f'{basepath}/Antarctic/{obsid}/final/*.dat')#glob.glob(f'{basepath}/{obsid}/final/*.dat') + \


    for file in files_dat:
        try:
            data = np.genfromtxt(file, names=True, dtype=None, encoding=None)
            #print(data.dtype.names)
            if data.size == 0:
                continue  # skip empty files

            for var in ['Tair', 'Tsur', 'SIT', 'SD', 'FRB', 'SID']:
                if var in data.dtype.names:
                    non_nan_count = np.count_nonzero(~np.isnan(data[var]))
                    print(f"{file} - {var}: {non_nan_count} non-NaN values")
                    #break  # stop after first matching variable
        except Exception as e:
            print(f"Error reading {file}: {e}")