## RRDPp Code Folder Overview

Each folder within the `RRDPp/code` directory contains scripts for converting raw observations into gridded monthly data:
- **25 km resolution** for the Northern Hemisphere (NH)
- **50 km resolution** for the Southern Hemisphere (SH)

For further details, see the associated publication:

*A first approach towards dual-hemisphere sea ice reference measurements from multiple data sources repurposed for evaluation and product intercomparison of satellite altimetry*

## Description of Scripts

- **EASEgrid_correct.py**  
  Creates temporal and spatial averages of input data based on specified temporal and spatial windows.

- **Functions.py**  
  Collection of supporting functions.

- **Warren.py**  
  Assigns Warren climatology snow depths based on time of year.

- **count_file_length.py**  
  Counts the number of observations in output files.

- **fix_outputfiles.py**  
  Updates relevant information in processed files.

- **run_final_files.sh**  
  Creates final files for all defined campaigns (as specified in the `OBSID` variable).

- **txt_to_nc.py**  
  Converts text files to NetCDF format.
