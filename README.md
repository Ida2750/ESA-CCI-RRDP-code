# Project's Title:
Dual hemisphere reference dataset for sea ice thickness, snow depth, draft and freeboard, Climate Change Initiative sea ice thickness Round Robin Data Package (CCI+ SIT RRDP)


# Project Description:
The Climate Change Initiative sea ice thickness Round Robin Data Package (CCI+ SIT RRDP) is a set of space delimitted ASCII files which contain a range of relevant reference sea ice variables, including sea ice thickness, snow depth, draft and freeboard, along with snow depth and density from Warren climatology (1999) for Arctic. Sea ice surface temperature and air temperature is also provided when available. The data has been gridded into a 25 by 25 km product for Arctic and 50 by 50 km for Antarctic using EASE-grid 2.0 and contains monthly averages of the given gridpoint. Reference data from autonomous buoys, moorings, submarines, ships (ice breakers), airborne campaigns (radar/laser altimetry and electromagnetic induction system) and in-situ measurements is provided in the final dataset. This dataset is available from DTU DATA (DOI:10.11583/DTU.23735679) and is described in the belonging publication with DOI: **INSERT DOI**.\

This repository provides the code used for processing from raw referende data to the final gridded product and the code used for co-locating reference data to satellite measurements from Envisat, CryoSat-2, ERS-1 and ERS-2. For examples of use, further information and inspiration we refer to the belonging paper: *"Dual-hemisphere sea ice thickness referencemeasurements from multiple data sources for evaluationand product inter-comparison of satellite altimetry"* with DOI: **INSERT DOI**
We highly recomend potential users of the repository to read the paper.

# How to Install and Run the Project:
The project was made using python version 3.9.13, along with the packages: numpy, matplotlib, cartopy, PyPDF2, netcdf4, pandas and datetime.
Furhtermore NASA's EaseGrid-2.0 was used and modified to make the gridded final data. This script is found in the code folder under the name: EASEgrid_correct

To use the project extract the data from the zip RRDPp zip folder. Download relevant reference data and Satellite data (see table 3 in **INSERT DOI**). See the How to Use the Project section for more information.

# The structure of the project:
The structure of the project is illustrated on the diagrams below, with the blue boxes showning folders and pink boxes showing scripts. Black arrows indicate folder/subfolder connections, whereas the green lines show data/code dependencies. The content of the satellite folder is shown individually.
![image](https://github.com/Idalundtorp/ESACCI-/assets/70795109/f52f888f-4e12-42a5-947f-6852c2dbf021)

![image](https://github.com/Idalundtorp/ESACCI-/assets/70795109/b0277fd5-89d6-4725-821b-92b5cd421c4d)


# How to Use the Project:
Currently the *../RRDPp/FINAL*, *../RRDPp/satellite/Final_files* and *../RRDPp/RawData* folders are left intentionally blank.

Data for the  *../RRDPp/FINAL* and *../RRDPp/satellite/Final_files* folders are readily available for use from DTU DATA
DOI: 10.11583/DTU.23735679 \

The *../RRDPp/code* folder contains the scripts for the individual campaigns which processes the raw data to the final gridded data.\
A folder with the naming convention of the *../RRDPp/code* folder should be made in the *../RRDPp/FINAL* folder.\

In the ../RRDPp/RawData folder the user should locate raw data in folders following the naming convention of the *../RRDPp/code* folder.\ 
Individual links to raw data are available from the belonging publication in table 3 **INSERT DOI** \

The *../RRDPp/satellite* folder contains scripts that are related to co-locating data from CryoSat-2, ENVISAT, ERS-1 and ERS-2 to data in the *../RRDPp/FINAL* folder.\

To do the co-location reference data from the CCI+ SIT RRDP DOI: 10.11583/DTU.23735679 and satelitte data from able 3 **INSERT DOI** must be downloaded.

For an overview of the satellite folder see the diagram above.\


# Credits
This project was made Henriette Skourup and Ida Olsen. Data used in the project were gathered by several external organisations. \
For data availability, credits and acknowledgements, please see the publication linked to this project and remember to acknowledge and \
cite relevant data sources when using data from this project\

# License
This project is made available under the CC BY 4.0 license meaning that you are free to:\
Share, copy and redistribute the material. Adapt, remix , transform and build upon the material.\

When using code, or data from this project please include the citation:\
**INSERT CITATION**

and the acknowledgement:\
This datapackage and code was made by Ida Olsen and Henriette Skourup from the Technical University of Denmark\

You are also required to provide a citation to the raw data sources, when using this data. \
