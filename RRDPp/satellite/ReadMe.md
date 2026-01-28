## Collocation Code (CCI-SIT-RRDP and Satellite Observations)

This code is used to collocate files from the CCI-SIT-RRDP with satellite observations.

### Workflow

1. Run `collocate_data.sh`
   
   1.1 Runs `Get_sat_subset.py`  
   Extracts relevant satellite subsets by calling `read_ccip_AutoROI_final.py`.

   1.2 Runs `collocation.py`  
   Collocates the satellite subsets with observations by calling `collocation_p2.py`.

2. Run `fix_outputfiles.py` from  
   `/dmidata/users/ilo/projects/RRDPp/code`  
   Adjusts output files where necessary (both before and after collocation).
