# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['ilo@dmi.dk']
__version__ = '0'
__date__ = '2022-07-11'


# -- Built-in modules -- #
import os
import sys
# -- Third-part modules -- #
from multiprocessing import Pool, cpu_count
import numpy as np
# -- Proprietary modules -- #
from collocation_p2 import collocation_part2


#---------------------------------
''' Define these parameters'''
if len(sys.argv) != 4:
    print("Usage: python script.py <HS> <SATELLITE> <OBSID>")
    sys.exit(1)

HEMISPHERE = sys.argv[1]
SATELLITE = sys.argv[2]
OBSID = sys.argv[3]

print(f"Running script with HS={HEMISPHERE}, SATELLITE={SATELLITE}, OBSID={OBSID}")

if OBSID=='SCICEX' or OBSID=='BGEP' or OBSID=='TRANSDRIFT' or OBSID=='AWI-ULS' or OBSID=='NPEO' or OBSID=='NPI':
    var = ['SID']
    name = f'ESACCIplus-SEAICE-RRDP2+-SID-{OBSID}'
else:
    var = [ 'SIT', 'SD', 'FRB'] # list of variables relevant for campaign (either ['SID'] or [ 'SIT', 'SD', 'FRB'])
    name = 'ESACCIplus-SEAICE-RRDP2+-SIT-' + OBSID # name of file in the FINAL data folder
    if OBSID=='SB_AWI' or OBSID=='OIB':
         name = f'ESACCIplus-SEAICE-RRDP2+-SIT-{OBSID}-{HEMISPHERE}' # contains data for both NH and SH
#---------------------------------

## Defines path to satellite files
parrent_dir = os.path.dirname(os.getcwd())
# parrent_dir = '/dmidata/projects/cmems2/C3S/RRDPp'
# CSdir = parrent_dir + '/satellite/Satellite_subsets/' + OBSID
CSdir = '/dmidata/projects/cmems2/C3S/RRDPp/satellite/Satellite_subsets/' + OBSID
## Defines obsfile
if HEMISPHERE == 'NH':
    obsdir = parrent_dir + '/FINAL/' + OBSID  + '/final/'
elif HEMISPHERE == 'SH':
    obsdir = parrent_dir + 'FINAL/Antarctic/' + OBSID  + '/final/'
obsfile = os.path.join(obsdir, name + '.dat')

## defines output file
ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', name  + '-' + SATELLITE + '-' + HEMISPHERE + '-CCIp-v3p0-rc2.dat')
if not os.path.exists(os.path.join(parrent_dir, 'satellite', 'Final_files')): os.makedirs(os.path.join(parrent_dir, 'satellite', 'Final_files'))

# loops over satellite files
count = 1

# for days in [4, 16, 30]:
#     for resolution in [5000, 15000, 25000]:
#         ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', f'{name}-{SATELLITE}-{HEMISPHERE}-{days}-{resolution}-CCIp-v3p0-rc2.dat')
#         for CSfile in sorted(os.listdir(CSdir)):
#             if CSfile.startswith(OBSID)  and CSfile.__contains__(SATELLITE) and CSfile.endswith('_CCIp.out'):
#                     CSfile = os.path.join(CSdir, CSfile)
#                     collocation_part2(obsfile, CSfile, count, ofile, var, HS=HEMISPHERE, days=days, resolution=resolution)
#                     count += 1

def combine_data(CSdir, OBSID, SATELLITE):
    all_data = []

    # Loop through and read each matching file
    for filename in sorted(os.listdir(CSdir)):
        if filename.startswith(OBSID) and SATELLITE in filename and filename.endswith('_CCIp.out'):
            print(filename)
            filepath = os.path.join(CSdir, filename)
            try:
                data = np.loadtxt(filepath)  # or use np.genfromtxt if there's missing data
                all_data.append(data)
            except Exception as e:
                print(f"Skipping {filename} due to error: {e}")

    # Combine all arrays into one big array
    if all_data:
        combined_array = np.vstack(all_data)
        print(f"Combined shape: {combined_array.shape}")
    else:
        print("No valid data files found.")

    output_path = os.path.join(CSdir, f"{OBSID}_{SATELLITE}_combined.out")

    np.savetxt(output_path, combined_array, fmt="%.6f")
    print(f"Saved combined data to: {output_path}")

def process_days_resolution(args):
    days, resolution = args
    count = 0

    ofile = os.path.join(
        parrent_dir, 'satellite', 'Final_files',
        f'{name}-{SATELLITE}-{HEMISPHERE}-{days}-{resolution}-test-CCIp-v3p0-rc2.dat'
    )

    CSfile = f"{OBSID}_{SATELLITE}_combined.out"
    CSfile_path = os.path.join(CSdir, CSfile)
    if not os.path.exists(CSfile_path):
        combine_data(CSdir, OBSID, SATELLITE)
    print('data combined')
    # for CSfile_name in sorted(os.listdir(CSdir)):
    #     if (CSfile_name.startswith(OBSID) and 
    #         SATELLITE in CSfile_name and 
    #         CSfile_name.endswith('_CCIp.out')):

    #CSfile_path = os.path.join(CSdir, CSfile)
    collocation_part2(
        obsfile,
        CSfile_path,
        count,
        ofile,
        var,
        HS=HEMISPHERE,
        days=days,
        resolution=resolution
    )
    count += 1

# Prepare all (days, resolution) combinations
param_combinations = [(d, r) for d in [30] for r in [25000]]

# Run parallel processes
if __name__ == '__main__':
    num_processes = 1  # Set this to however many parallel runs you want
    with Pool(processes=num_processes) as pool:
        pool.map(process_days_resolution, param_combinations)