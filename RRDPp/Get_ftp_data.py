# -*- coding: utf-8 -*-
"""
Script to get data from ftp server
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2022-07-17'

# -- Built-in modules -- #
import os
# -- Third-part modules -- #
import ftplib
# -- Proprietary modules -- #


HOSTNAME = "anon-ftp.ceda.ac.uk"
USERNAME = "anonymous"
PASSWORD = "" 
DOWNLOAD_DIR = os.getcwd() +"/satelitte/ENV_CS2_data"
HS = 'NH'
sat = 'ENV'

YEARS_ENV = [str(2002 + i) for i in range(11)]
YEARS_CS2 = [str(2010 + i) for i in range(12)]
MONTHS = [str(0) + str(1 + i) if i < 9 else str(i + 1) for i in range(12) ]

def get_along_track(YEARS, MONTHS):
    for YEAR in YEARS:
        os.chdir(DOWNLOAD_DIR + "/"+sat+"_data/" + HS)
        # make directory of year and enter
        try:
            os.mkdir(YEAR)
        except:
            pass
        # os.chdir(YEAR)
        for MONTH in MONTHS:
            os.chdir(DOWNLOAD_DIR + "/"+sat+"_data/" + HS + "/" + YEAR)
            # make directory of month and enter
            try:
                os.mkdir(MONTH)
            except:
                pass
            os.chdir(MONTH)
          
            # Connect FTP Server
            ftp_server = ftplib.FTP(HOSTNAME, USERNAME)
             
            # force UTF-8 encoding
            ftp_server.encoding = "utf-8"

            # Change dir of ftp server
            if sat=='CS2':
                satelitte = "cryosat2"
            else:
                satelitte = "envisat"
            data_dir = "neodc/esacci/sea_ice/data/sea_ice_thickness/L2P/"+satelitte+"/v2.0/"+HS+"/"
            ftp_server.cwd(data_dir)
            ftp_server.cwd(YEAR)
            try:
                ftp_server.cwd(MONTH)
                
                # Enter File Name with Extension
                DAYS_L = [str(num) for num in range(1,32) if num>10]
                DAYS_S = ["0" + str(num) for num in range(1,32) if num<10]
                DAYS = DAYS_S + DAYS_L
                
        
                if sat=='CS2':
                    filenames = ["ESACCI-SEAICE-L2P-SITHICK-SIRAL_CRYOSAT2-"+ HS +"-" + YEAR + MONTH + DAY + "-fv2.0.nc" for DAY in DAYS]
                else:
                    filenames = ["ESACCI-SEAICE-L2P-SITHICK-RA2_ENVISAT-"+ HS +"-" + YEAR + MONTH + DAY + "-fv2.0.nc" for DAY in DAYS]
                    
                # Write file in binary mode
                for filename in filenames:
                    print(filename)
                    try:
                        with open(filename, "wb") as file:
                            # Command for Downloading the file "RETR filename"
                            ftp_server.retrbinary(f"RETR {filename}", file.write)
                    except:
                        print('No file called:', filename)
                        # remove empty files
                        os.remove(filename)
            except:
                pass 
            # Close the Connection
            ftp_server.quit()


get_along_track(YEARS_ENV, MONTHS)

