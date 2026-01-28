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
import glob
# -- Third-part modules -- #
import numpy as np

# -- Proprietary modules -- #
from collocation_p2 import collocation_part2
sys.path.append(os.path.dirname(os.getcwd()) + '/code') 
import Functions

#---------------------------------
''' Define these parameters'''
if len(sys.argv) != 4:
    print("Usage: python script.py <HS> <SATELLITE> <OBSID>")
    sys.exit(1)

HEMISPHERE = sys.argv[1]
SATELLITE = sys.argv[2]
OBSID = sys.argv[3]

print(f"Running script with HS={HEMISPHERE}, SATELLITE={SATELLITE}, OBSID={OBSID}")

# if OBSID=='SCICEX' or OBSID=='BGEP' or OBSID=='TRANSDRIFT' or OBSID=='AWI-ULS' or OBSID=='NPEO' or OBSID=='NPI':
#     var = ['SID']
#     name = f'ESACCIplus-SEAICE-RRDP2+-SID-{OBSID}'
# else:
#     var = [ 'SIT', 'SD', 'FRB'] # list of variables relevant for campaign (either ['SID'] or [ 'SIT', 'SD', 'FRB'])
#     name = 'ESACCIplus-SEAICE-RRDP2+-SIT-' + OBSID # name of file in the FINAL data folder
#     if OBSID=='SB_AWI' or OBSID=='OIB':
#          name = f'ESACCIplus-SEAICE-RRDP2+-SIT-{OBSID}-{HEMISPHERE}' # contains data for both NH and SH
#---------------------------------

## Defines path to satellite files
parrent_dir = os.path.dirname(os.getcwd())
# parrent_dir = '/dmidata/projects/cmems2/C3S/RRDPp'
# CSdir = parrent_dir + '/satellite/Satellite_subsets/' + OBSID
CSdir = '/dmidata/projects/cmems2/C3S/RRDPp/satellite/Satellite_subsets/' + OBSID

# Define obsdir based on hemisphere and OBSID
if HEMISPHERE == 'NH':
    obsdir = os.path.join(parrent_dir, 'FINAL', OBSID, 'final')
elif HEMISPHERE == 'SH':
    obsdir = os.path.join(parrent_dir, 'FINAL', 'Antarctic', OBSID, 'final')

# Pattern to match observation files — adjust if needed
pattern = os.path.join(obsdir, f'ESACCIplus-SEAICE-RRDP2+-*{OBSID}*.nc')

# Find all matching .nc files
nc_files = glob.glob(pattern)

if len(nc_files)==0:
    OBSID_filename = OBSID.replace('_', '-')
    pattern = os.path.join(obsdir, f'ESACCIplus-SEAICE-RRDP2+-*{OBSID_filename}*.nc')
    nc_files = glob.glob(pattern)

dat_files = []
names = []
vars = []
for nc_file in nc_files:
    dat_file = nc_file.replace('.nc', '.dat')
    Functions.netcdf_to_txt(nc_file, dat_file)  # Your conversion function
    dat_files.append(dat_file)
    if 'SID' in dat_file:
        names.append(f'ESACCIplus-SEAICE-RRDP2+-SID-{OBSID}')
        vars.append(['SID'])
    elif OBSID=='SB_AWI' or OBSID=='OIB':
        names.append(f'ESACCIplus-SEAICE-RRDP2+-SIT-{OBSID}-{HEMISPHERE}')
        vars.append([ 'SIT', 'SD', 'FRB'])
    elif OBSID=='NICE' and 'SD' in os.path.basename(nc_file):
        names.append(f'ESACCIplus-SEAICE-RRDP2+-SD-{OBSID}')
        vars.append([ 'SIT', 'SD', 'FRB'])
    else:
        names.append(f'ESACCIplus-SEAICE-RRDP2+-SIT-{OBSID}')
        vars.append([ 'SIT', 'SD', 'FRB'])

#print(dat_files)
#dat_files = [dat_files[1]]
print(dat_files)
for obsfile, name, var in zip(dat_files, names, vars):

    ## defines output file
    if OBSID=='MOSAIC':
        if 'HEM' in obsfile:
            ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', name  + '-HEM-' + SATELLITE + '-' + HEMISPHERE + '-CCIp-v3p0-rc2.dat')
        elif 'SIMBA' in obsfile:
            ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', name  + '-SIMBA-' + SATELLITE + '-' + HEMISPHERE + '-CCIp-v3p0-rc2.dat')
    else:
        ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', name  + '-' + SATELLITE + '-' + HEMISPHERE + '-CCIp-v3p0-rc2.dat')
    print(ofile)
    if not os.path.exists(os.path.join(parrent_dir, 'satellite', 'Final_files')): os.makedirs(os.path.join(parrent_dir, 'satellite', 'Final_files'))
    if os.path.exists(ofile): os.remove(ofile)

    def combine_data(CSdir, OBSID, SATELLITE, HS='NH'):
        all_data = []

        # Loop through and read each matching file
        for filename in sorted(os.listdir(CSdir)):
            if HS=='NH':
                if filename.startswith(OBSID) and 'SH' not in filename and SATELLITE in filename and filename.endswith('_CCIp.out'):
                    print(filename)
                    filepath = os.path.join(CSdir, filename)
                    try:
                        data = np.loadtxt(filepath)  # or use np.genfromtxt if there's missing data
                        all_data.append(data)
                    except Exception as e:
                        print(f"Skipping {filename} due to error: {e}")
            elif HS=='SH':
                if filename.startswith(OBSID) and 'SH' in filename and SATELLITE in filename and filename.endswith('_CCIp.out'):
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

        if HS=='NH':
           output_path = os.path.join(CSdir, f"{OBSID}_{SATELLITE}_combined.out") 
        elif HS=='SH':
            output_path = os.path.join(CSdir, f"{OBSID}_{SATELLITE}_{HS}_combined.out")

        np.savetxt(output_path, combined_array, fmt="%.6f")
        print(f"Saved combined data to: {output_path}")

    # loops over satellite files
    count = 1
    if HEMISPHERE=='NH':
        CSfile = f"{OBSID}_{SATELLITE}_combined.out"
    elif HEMISPHERE=='SH':
        CSfile = f"{OBSID}_{SATELLITE}_{HEMISPHERE}_combined.out"
    CSfile_path = os.path.join(CSdir, CSfile)
    if not os.path.exists(CSfile_path):
        combine_data(CSdir, OBSID, SATELLITE, HS=HEMISPHERE)
    print('data combined')
    collocation_part2(obsfile, CSfile_path, count, ofile, var, HS=HEMISPHERE)
    count += 1

    # # loops over satellite files
    # count = 1
    # for CSfile in sorted(os.listdir(CSdir)):
    #     if CSfile.startswith(OBSID)  and CSfile.__contains__(SATELLITE) and CSfile.endswith('_CCIp.out'):
    #             #print(CSfile)
    #             CSfile = os.path.join(CSdir, CSfile)
    #             collocation_part2(obsfile, CSfile, count, ofile, var, HS=HEMISPHERE)
    #             count += 1
