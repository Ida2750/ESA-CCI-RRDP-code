# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 17:24:09 2024

@author: Ida Olsen
"""

from glob import glob
import Functions

datapath='C:/Users/Ida Olsen/Documents/work/RRDPp/satellite/Final_files/final'

files = glob(f'{datapath}/A*/*.dat')
for ifile in files:
    if 'SID' in ifile:
        key_variables='Sea Ice Draft'
        primary = 'SID'
        
    if 'IMB' in ifile:    
        primary = 'SIT'
        datasource = ' Monitoring the mass balance, motion, and thickness of Arctic sea ice, http://imb-crrel-dartmouth.org'
        key_variables='Sea ice thickness and snow depth'
    elif 'SB-AWI' in ifile:
        primary = 'SD'
        datasource = 'Autonomous measurements of sea ice properties: https://data.meereisportal.de/relaunch/buoy.php?lang=en&active-tab1=method&active-tab2=buoy&ice-type=buoy&region=all&buoytype=SB&buoystate=all&expedition=all&showMaps=y&dateRepeat=n&timeline='
        key_variables = 'Snow depth'
    elif 'OIB' in ifile and 'SH' not in ifile:
        primary = 'FRB'
        datasource = 'IceBridge Sea Ice Freeboard, Snow Depth, and Thickness Quick Look, Version 1, doi: 10.5067/GRIXZ91DE0L9 + IceBridge L4 Sea Ice Freeboard, Snow Depth, and Thickness, Version 1, doi:10.5067/G519SHCKWQV6'
        key_variables = 'Sea ice thickness, sea ice freeboard and snow depth'
    elif 'OIB' in ifile and 'SH' in ifile:
        primary = 'FRB'
        datasource = 'IceBridge L4 Sea Ice Freeboard, Snow Depth, and Thickness, Version 1, doi:10.5067/G519SHCKWQV6'
        key_variables = 'Freeboard: sea ice freeboard + snow depth'
    elif 'AEM-AWI' in ifile and 'SIT' in ifile:
        primary = 'SIT'
        datasource = 'Airborne Electromagnetic Measurement, Alfred Wegener Institute for Polar and Marine Research, doi:10.2312/polfor.2016.011.'
        key_variables = 'Ice thickness: sea ice thickness + snow depth'
    elif 'AEM-AWI' in ifile and 'FRB' in ifile:
        primary = 'FRB'
        datasource = 'Airborne Electromagnetic Measurement, Alfred Wegener Institute for Polar and Marine Research, doi:10.2312/polfor.2016.011.'
        key_variables = 'Freeboard: sea ice freeboard + snow depth'
    elif 'ASSIST' in ifile:
        primary = 'SIT'
        datasource = 'ASSIST: https://icewatch.met.no/'
        key_variables = 'Sea ice thickness and snow depth'
    elif 'ASPeCt' in ifile:
        primary = 'SIT'
        datasource = 'ASPeCt, 1980-2005: Thickness distribution of Antarctic sea ice doi:10.1029/2007JC004254, 2005-2019: ESA-CCI_Phase2_Standardized_Manual_Visual_Ship-Based_SeaIceObservations_v02, doi:10.26050/WDCC/ESACCIPSMVSBSIOV2'
        key_variables = 'Sea ice thickness and snow depth'
    elif 'AWI-ULS' in ifile:
        datasource = 'Sea ice draft measured by upward looking sonars in the Weddell Sea (Antarctica): https://doi.org/10.1594/PANGAEA.785565'
    elif 'BGEP' in ifile:
        datasource='Beaufort Gyre Exploration Project, Mooring Data: https://www2.whoi.edu/site/beaufortgyre/data/mooring-data/'
    elif 'NPEO' in ifile:
        datasource = 'North Pole Environmental Observatory (NPEO) Oceanographic Mooring Data: https://doi.org/10.5065/D6P84921'
    elif 'NPI' in ifile:
        datasource = 'Monthly sea ice thickness distribution in Fram Strait. Norwegian Polar Institue. https://doi.org/10.21334/npolar.2022.b94cb848'
    elif 'SCICEX' in ifile:
        datasource = 'SCICEX: Submarine Arctic Science Program, Version 1: https://doi.org/10.7265/N5930R3Z + Submarine Upward Looking Sonar Ice Draft Profile Data and Statistics, Version 1: https://doi.org/10.7265/N54Q7RW'
    elif 'TRANSDRIFT' in ifile:
        datasource = 'Daily mean sea ice draft from moored Upward-Looking Sonars in the Laptev Sea between 2013 and 2015: https://doi.pangaea.de/10.1594/PANGAEA.899275 + Daily mean sea ice draft from moored upward-looking Acoustic Doppler Current Profilers (ADCPs) in the Laptev Sea from 2003 to 2016: https://doi.pangaea.de/10.1594/PANGAEA.912927'
                
    Functions.txt_to_netcdf(datapath, ifile, primary, datasource, key_variables)
    

#%%

# from glob import glob
# import Functions

# datapath='C:/Users/Ida Olsen/Documents/work/ESA-CCI-RRDP-code/RRDPp/FINAL/*/final'
# #datapath='C:/Users/Ida Olsen/Documents/work/RRDPp/FINAL/*/*/final'

# files = glob(f'{datapath}/*FRB*.nc')
# for ifile in files:
#     print(ifile)
#     #try:
#     Functions.get_count(ifile)
#     #except:
#     #    pass