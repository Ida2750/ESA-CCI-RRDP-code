#!/bin/bash

OBSID=( 'OIB') #'SCICEX' 'TRANSDRIFT' 'NPI' 'NPEO' 'BGEP' 'AWI-ULS' 'AEM-AWI' 'ASPeCt' 'MOSAIC' 'Nansen_legacy') # 'ASSIST' 'OIB' 'IMB' 'SB_AWI')
SATELLITE=('CS2')
HS=('NH')

for hs in "${HS[@]}"; do
    for sat in "${SATELLITE[@]}"; do
        for id in "${OBSID[@]}"; do
            #python -W ignore Get_sat_subset.py "$hs" "$sat" "$id"
            python -W ignore representation_error_study.py "$hs" "$sat" "$id"
        done
    done
done

