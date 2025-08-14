# merge netcdf files

import xarray as xr

# Open both NetCDF files
ds1 = xr.open_dataset("/dmidata/users/ilo/projects/RRDPp/FINAL/NICE/final/ESACCIplus-SEAICE-RRDP2+-SIT-NICE-HEM-p1_v2.nc")
ds2 = xr.open_dataset("/dmidata/users/ilo/projects/RRDPp/FINAL/NICE/final/ESACCIplus-SEAICE-RRDP2+-SIT-NICE-HEM-p2_v2.nc")

# Merge along a dimension (e.g., time)
combined = xr.concat([ds1, ds2], dim="time")

# Update the global attribute
combined.attrs["Datasource"] = (
    "Airborne Electromagnetic Measurement, N-ICE2015, "
    "DOI:https://doi.org/10.21334/NPOLAR.2016.70352512 + https://doi.org/10.21334/NPOLAR.2016.AA3A5232"
)

# Save the combined dataset
combined.to_netcdf("/dmidata/users/ilo/projects/RRDPp/FINAL/NICE/final/ESACCIplus-SEAICE-RRDP2+-SIT-NICE-HEM.nc")
