# Project's Title:
Code for sea ice thickness reference measurements


# Project Description:
The Climate Change Initiative sea ice thickness Round Robin Data Package (CCI SIT RRDP) is a set of space delimitted ASCII files which contain a range of relevant reference sea ice variables, including sea ice thickness, snow depth, draft and freeboard, along with snow depth and density from Warren climatology (1999) for Arctic. Sea ice surface temperature and air temperature is also provided when available. The data has been gridded into a 25 by 25 km product for Arctic and 50 by 50 km for Antarctic using EASE-grid 2.0 and contains monthly averages of the given gridpoint. Reference data from autonomous buoys, moorings, submarines, ships (ice breakers), airborne campaigns (radar/laser altimetry and electromagnetic induction system) and in-situ measurements is provided in the final dataset. This dataset is available from DTU DATA ([DOI:10.11583/DTU.23735679](https://doi.org/10.11583/DTU.24787341)) and described in the publication [preprint] [DOI:10.5194/essd-2024-234](https://doi.org/10.5194/essd-2024-234).

This repository provides the code used for processing from raw reference data to the final gridded product and the code used for co-locating reference data to satellite measurements from Envisat, CryoSat-2, ERS-1 and ERS-2. For examples of use, further information and inspiration we refer to the belonging paper [DOI:10.5194/essd-2024-234](https://doi.org/10.5194/essd-2024-234).

# How to Install and Run the Project:
The project was made using python version 3.13.3, along with the packages: numpy, matplotlib, cartopy, PyPDF2, netcdf4, pandas and more (see associated rrdp.yml file to copy the environment structure) \
To use the project download relevant reference data and Satellite data (see table 3 in [preprint] [DOI:10.5194/essd-2024-234](https://doi.org/10.5194/essd-2024-234))

# How to Use the Project:
Currently the *../RRDPp/FINAL*, *../RRDPp/satellite/Final_files* and *../RRDPp/RawData* folders are left intentionally blank.

Data for the  *../RRDPp/FINAL* and *../RRDPp/satellite/Final_files* folders are readily available for use from DTU DATA ([DOI:10.11583/DTU.23735679](https://doi.org/10.11583/DTU.24787341))

The *../RRDPp/code* folder contains the scripts for the individual campaigns which processes the raw data to the final gridded data.
A folder with the naming convention of the *../RRDPp/code* folder should be made in the *../RRDPp/FINAL* folder.

In the ../RRDPp/RawData folder the user should locate raw data in folders following the naming convention of the *../RRDPp/code* folder.
Individual links to raw data are available from the belonging publication in table 3 [preprint] [DOI:10.5194/essd-2024-234](https://doi.org/10.5194/essd-2024-234)

The *../RRDPp/satellite* folder contains scripts that are related to co-locating data from CryoSat-2, ENVISAT, ERS-1 and ERS-2 to data in the *../RRDPp/FINAL* folder.

To do the co-location reference data from the CCI SIT RRDP ([DOI:10.11583/DTU.23735679](https://doi.org/10.11583/DTU.24787341)) and satelitte data from table 3 [preprint] [DOI:10.5194/essd-2024-234](https://doi.org/10.5194/essd-2024-234) must be downloaded.

To change the spatial and temporal resolutions from the default (30 days and 25km (northern hemisphere) and 50km (southern hemisphere)) one must adjust the parameters gridres = 25000  # grid resolution and
dtint = 30  # days per mean, which are defined in the start of the processing scripts. For collocation the script collocate.py must be adjusted by defining the days and resolution parameters in the function collocation_part2.

# Example use (run BGEP)
1. Download relevant data
2. Create environment using conda env create -f rrdp.yml
4. Run python RRDPp/code/BGEP_processing_final.py (adjust dtint = 30 and gridres = 25000 if necessary) -> creates output file and figures in *../RRDPp/FINAL/BGEP/final* and *../RRDPp/FINAL/BGEP/fig* repectively
6. Run ./RRDPp/satellitte/collocation.sh (adjust OBSID=('BEGP'), SATELLITE=('CS2' 'ENV') and HS=('NH')) and adjust days and resolution in collocation.py if necessary)
7. Run python RRDPp/code/fix_outputfiles.py

# Credits
This project was made Henriette Skourup and Ida Lundtorp Olsen. Data used in the project were gathered by several external organisations. 
For data availability, credits and acknowledgements, please see the publication linked to this project and remember to acknowledge and 
cite relevant data sources when using data from this project

# License
This project is made available under the CC BY 4.0 license meaning that you are free to:
Share, copy and redistribute the material. Adapt, remix , transform and build upon the material.

When using code, or data from this project please include the citation:
Olsen, I. L. and Skourup, H.: Sea ice thickness reference measurements (ESA CCI SIT RRDP), Dataset, https://doi.org/10.11583/DTU.24787341, 2024.

and the acknowledgement:
This datapackage and code was made by Ida Lundtorp Olsen from the Danish Meteorological Institute and Henriette Skourup from the Technical University of Denmark

You are also required to provide a citation to the raw data sources, when using this data. 

# Contact
Ida Lundtorp Olsen
ilo@dmi.dk or ida.lundtorp@outlook.dk

