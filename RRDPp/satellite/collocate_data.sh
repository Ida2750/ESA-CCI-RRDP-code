#!/bin/bash

OBSID=('OIB') #'AWI_ULS') # 'ASPeCt' 'IMB') # 'MOSAIC' 'IMB' 'AEM-AWI' ASPeCt' )
SATELLITE=('CS2')
HS=('NH') #'NH')

for hs in "${HS[@]}"; do
    for sat in "${SATELLITE[@]}"; do
        for id in "${OBSID[@]}"; do
            #python -W ignore Get_sat_subset.py "$hs" "$sat" "$id"
            python -W ignore collocation.py "$hs" "$sat" "$id"
        done
    done
done

