#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 12:04:04 2022

@author: s174020

Get dates

DATE/TIME START: 2004-02-05T12:10:00 * DATE/TIME END: 2004-02-05T12:27:00 
"""

# file = '/media/s174020/My Passport2/ESACCIplus/RRDPp/RawData/AEM_AWI/Jul_Aug_2004/HEM_ARK20_01_sea-ice_thickness.tab'
# ifile = 'HEM_ARK20_01_sea-ice_thickness.tab'


def Get_date(file, ifile):
    import re
    import numpy as np
    import datetime as dt
    lookup_start = '\tDATE/TIME START:'
    #lookup_end = 'DATE/TIME END'
    lookup = '*/'
    count = 0
    header = -999
    numbers = []
    with open(file,'rb') as myFile:
        for num, line in enumerate(myFile, 1):
            count+=1
            # break loop when lookup is found
            if lookup_start in line.decode('utf8'):  # find date information in header
                numbers = re.findall(r'[0-9]+', line.decode('utf8'))
                print(line)
                date_start = numbers[0] + numbers[1] + numbers[2] + numbers[3] + numbers[4] + numbers[5]
                date_end = numbers[6] + numbers[7] + numbers[8] + numbers[9] + numbers[10] + numbers[11]
    
            if lookup in line.decode('utf8'):
                header = count  # note number of headerlines
            if count == header+1:
                headerline = line.decode('utf8').replace("[m]","")
                names = headerline.split()
                break
    data = np.genfromtxt(file, dtype=None, names=names, skip_header=header+1) 
    
    latitude = data['Latitude']
    
    if date_start == date_end:
        t = dt.datetime.strptime(date_start, "%Y%m%d%H%M%S")
        date = [t + dt.timedelta(seconds=i) for i in range(len(latitude))]  # add progressively 1 second to ensure that date changes
    
    return date