# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = 'Henriette Skorup'
__contact__ = ['s174020@student.dtu.dk']
__version__ = '0'
__date__ = '2022-07-11'


# -- Built-in modules -- #
import os
# -- Third-part modules -- #

# -- Proprietary modules -- #
from collocation_p2 import collocation_part2


#---------------------------------
''' Define these parameters'''
OBSID      = 'TRANSDRIFT' # Name of campaign
SATELLITE = 'CS2' # Satelitte name (either ENV, CS2, ERS1 or ERS2)
HEMISPHERE = 'NH' # Hemisphere of observations (either Northen Hemisphere (NH) or Southern Hemisphere (SH))
var        =  ['SID'] #[ 'SIT', 'SD', 'FRB'] # list of variables relevant for campaign (either ['SID'] or [ 'SIT', 'SD', 'FRB'])
name       = 'ESACCIplus-SEAICE-RRDP2+-SID-' + OBSID # name of file in the FINAL data folder
#---------------------------------

## Defines path to satellite files
parrent_dir = os.path.dirname(os.getcwd())
CSdir = parrent_dir + '/satellite/Satellite_subsets/' + OBSID

## Defines obsfile
if HEMISPHERE == 'NH':
    obsdir = parrent_dir + '/FINAL/' + OBSID  + '/final/'
elif HEMISPHERE == 'SH':
    obsdir = parrent_dir + 'FINAL/Antarctic/' + OBSID  + '/final/'
obsfile = os.path.join(obsdir, name + '.dat')

## defines output file
ofile = os.path.join(parrent_dir, 'satellite', 'Final_files', name  + '-' + SATELLITE + '-CCIp-v3p0-rc2.dat')

# loops over satellite files
count = 1
for CSfile in os.listdir(CSdir):
    if CSfile.startswith(OBSID)  and CSfile.__contains__(SATELLITE) and CSfile.endswith('_CCIp.out'):
            print(CSfile)
            CSfile = os.path.join(CSdir, CSfile)
            collocation_part2(obsfile, CSfile, count, ofile, var, HS=HEMISPHERE)
            count += 1