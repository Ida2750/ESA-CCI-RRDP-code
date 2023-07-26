# Project's Title:
Dual hemisphere reference dataset for sea ice thickness, snow depth, draft and freeboard, 
Climate Change Initiative sea ice thickness Round Robin Data Package (CCI+ SIT RRDP)


# Project Description:
The Climate Change Initiative sea ice thickness Round Robin Data Package (CCI+ SIT RRDP) is a set of space delimitted ASCII files which contain a 
range of relevant reference sea ice variables, including sea ice thickness, snow depth, draft and freeboard, along with snow depth and density from 
Warren climatology (1999) for Arctic. Sea ice surface temperature and air temperature is also provided when available.
The data has been gridded into a 25 by 25 km product for Arctic and 50 by 50 km for Antarctic using EASE-grid 2.0 and contains monthly averages 
of the given gridpoint.
Reference data from autonomous buoys, moorings, submarines, ships (ice breakers), airborne campaigns (radar/laser altimetry and 
electromagnetic induction system) and in-situ measurements is provided in the final dataset. 
Along with the final dataset is provided code to process data from raw files into the final files. 

For further information please read the paper: " Dual-hemisphere sea ice thickness reference
measurements from multiple data sources for evaluation
and product inter-comparison of satellite altimetry" with DOI: **INSERT DOI**

# How to Install and Run the Project:
The project was made using python version 3.9.13, along with the packages:
numpy, EaseGrid-2.0 (link), matplotlib, cartopy, PyPDF2, pandas and datetime.
The version of EaseGrid-2.0 used in this project is found in the code folder under the name: EASEgrid_correct

To use the project extract the data from the zip RRDPp zip folder. 
Download relevant reference data and Satellite data (see table 3 in **INSERT DOI**)
See the How to Use the Project section for more information.

# The structure of the project:

![image](https://github.com/Idalundtorp/ESACCI-/assets/70795109/f52f888f-4e12-42a5-947f-6852c2dbf021)

![image](https://github.com/Idalundtorp/ESACCI-/assets/70795109/b0277fd5-89d6-4725-821b-92b5cd421c4d)


# How to Use the Project:
In the ../RRDPp/FINAL folder you find the gridded reference data, which is available for use.
The ../RRDPp/code folder contains the scripts for the individual campaigns which processes the raw data to the final gridded data.
In the ../RRDPp/RawData folder the user should locate raw data in folders following the naming convention of the ../RRDPp/code folder. 
Individual links to raw data are available from the belonging publication in table 3 INSERT DOI
The ../RRDPp/satellite folder contains scripts that are related to co-locating data from CryoSat-2, ENVISAT, ERS-1 and ERS-2 to data in the ../RRDPp/FINAL folder.
For an overview of the structure of this folder see the manual.


# Credits
This project was made Henriette Skourup and Ida Olsen. Data used in the project were gathered by several external organisations. 
For data availability, credits and acknowledgements, please see the publication linked to this project and remember to acknowledge and 
cite relevant data sources when using data from this project

# License
This project is made available under the CC BY 4.0 license meaning that you are free to:
Share, copy and redistribute the material. Adapt, remix , transform and build upon the material.

When using code, or data from this project please include the citation:
## INSERT CITATION

and the acknowledgement:
This datapackage and code was made by Ida Olsen and Henriette Skourup from the Technical University of Denmark

You are also required to provide a citation to the raw data sources, when using this data! 

-------------------------------------------------------------------------------------------------------------------------------------------------------
# Information about reference data final files.
-------------------------------------------------------------------------------------------------------------------------------------------------------

For files with the primary variable being SID.

Header includes:
obsID, date, lat, lon, SID, SIDstd, SIDln, SIDunc, wSD, w-rho, pp-flag, unc-flag

SID = sea ice draft in units of meters (m)
SIDstd = standard deviation of sea ice draft measurements used for each grid value (m)
SIDln = number of observations in each snow depth grid value
SIDunc = estimated uncertainty for each grid value (m).
wSD = snow depth from Warren climatology in units of (m)
w-rho = estimated density of snow from Warren climatology in units of kg/m^3
pp-flag = has values 0 to 3, with 0 meaning that no out of the ordinary processing was required of the raw data
unc-flag = has values 0 to 3, with 0 meaning that individual uncertainties were available for all raw data
For more information about pre processing and uncertainty flags see the belonging publication **INSERT DOI**


## Reference data sources SID files

ARCTIC sources:
BGEP: Stationary moored upward-looking sonar - Beufort Gyre Exploration Project
TRANSDRIFT: Stationary moored upward-looking sonar - Russian-German TRANSDRIFT project
NPEO:  Stationary moored upward-looking sonar - North Pole Enviromental Observatory
NPI:  Stationary moored upward-looking sonar - from the Fram Strait Arctic Outflow Observatory by the Norwegian Polar Institute 
SCICEX: Submarine-mounted upward-looking sonar - Submarine Arctic Science Program

ANTARCTIC sources:
AWI-ULS: Stationary moored upward-looking sonar - Alfred Wegener Institute

-------------------------------------------------------------------------------------------------------------------------------------------------------
For files with the primary variable being sea ice thickness (SIT), snow depth (SD) or freeboard (FRB).

Header includes:
obsID, date, lat, lon, SD, SDstd, SDln, SDunc, SIT, SITstd, SITln, SITunc, FRB, FRBstd, FRBln, FRBunc, Tsur, Tair, wSD, w-rho, pp-flag, unc-flag

SIT = sea ice thickness (m)
SD = snow depth (m)
FRB = freeboard (m). This is the sea ice freeboard for NP and AEM-AWI, whereas OIB measures the total freeboard.

The resulting variables follow the structure of the files with primary variable being SID.

-------------------------------------------------------------------------------------------------------------------------------------------------------
## Reference data sources SIT, SD and FRB files

ARCTIC sources:
OIB: Airborne laser and radar altimetry and snow radar - NASA Operation IceBridge
AEM-AWI: Airborne electromagnetic measurements - Alfred Wegener Institute
ASSIST: Visual observations from ice breakers - Arctic Shipborne Sea Ice Standardization Tool 
SB-AWI:  Snow depth buoys - Alfred Wegener Institute
IMB-CRREL: Ice mass balance buoys - Cold Regions Research and Engineering Laboratory
NP: Drifting research stations - North Pole Drifting Stations

ANTARCTIC sources:
OIB: Airborne laser and radar altimetry
SB-AWI: Snow depth buoys - Alfred Wegener Institute
ASPeCt: Visual observations from ice breakers - Antarctic Sea ice Processes and Climate
