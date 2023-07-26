#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DISCLAIMER: This is a modified version of NASA's Ease-Grid 2.0

The only modification made compared to the original scripts is to include the
computation of additional variables with the grid.

The method of computation follows the original script.
"""

# -- File info -- #
__author__ = 'NASA'
__contributors__ = 'Ida Olsen, Henriette Skorup'
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date_modified__ = '2022-07-13'

class Gridded:
    def __init__(self):
        import pdb;
        import pyproj
        import numpy as np

        self.str_proj = '+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
        self.EASE_proj = pyproj.Proj(self.str_proj)
        self.spatial_resolution = "25.0 km grid spacing" 

        return


    def SetHemisphere(self,str_hemisphere):
        import numpy as np
        import pyproj

        if str_hemisphere == 'S':
            self.hemisphere = 'South'
            self.geospatial_lat_max = np.float32(-16.62393)
            self.geospatial_lat_min = np.float32(-90.0)
            self.str_proj = '+proj=laea +lat_0=-90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
            self.EASE_proj = pyproj.Proj(self.str_proj)
        if str_hemisphere == 'N':
            self.hemisphere = 'North'
            self.geospatial_lat_min = np.float32(16.62393)
            self.geospatial_lat_max = np.float32(90.0)
            self.str_proj = '+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84 +datum=WGS84 +units=m'
            self.EASE_proj = pyproj.Proj(self.str_proj)



    def CreateGrids(self, resolution):
        import numpy as np
    #    import pyproj
    #    import Warren

        # Default for 100 km EASE grid!???
        self.max_x = 5400000.0
        self.min_x = -5400000.0

        self.max_y = 5400000.0
        self.min_y = -5400000.0


        self.resolution = resolution

        self.xc = np.arange(self.min_x,self.max_x, self.resolution) + self.resolution/2
        self.yc = np.arange(self.max_y,self.min_y, -self.resolution) - self.resolution/2

        lon = np.zeros([len(self.yc),len(self.xc)])
        lat = np.zeros([len(self.yc),len(self.xc)])

        for i in range(len(self.xc)):
            for j in range(len(self.yc)):
                [lon[j,i], lat[j,i]] = self.EASE_proj(self.xc[i], self.yc[j],inverse=True)
                

        self.lat = lat
        self.lon = lon

        sd_list = list()
        sd_unc_list = list()
        sit_list = list()
        sit_unc_list = list()
        Frb_list = list()
        Frb_unc_list = list()
        lat_list = list()
        lon_list = list()
        tn_list = list()
        td_list = list()

        #Create 2D empty list of lists corresponding to the grid-size
        for i in range(len(self.xc)):
            sd_list.append(list())
            sd_unc_list.append(list())
            sit_list.append(list())
            sit_unc_list.append(list()) 
            Frb_list.append(list())
            Frb_unc_list.append(list())
            lon_list.append(list())
            lat_list.append(list())
            tn_list.append(list())
            td_list.append(list())
            for j in range(len(self.yc)):
                sd_list[i].append(list())
                sd_unc_list[i].append(list())
                sit_list[i].append(list())
                sit_unc_list[i].append(list()) 
                Frb_list[i].append(list())
                Frb_unc_list[i].append(list())
                lat_list[i].append(list())
                lon_list[i].append(list())
                td_list[i].append(list())
                tn_list[i].append(list())
        self.sd_list = sd_list
        self.sd_unc_list = sd_unc_list
        self.sit_list = sit_list
        self.sit_unc_list = sit_unc_list
        self.Frb_list = Frb_list
        self.Frb_unc_list = Frb_unc_list
        self.lon_list = lon_list
        self.lat_list = lat_list
        self.td_list = td_list
        self.tn_list = tn_list
        
        self.n = []
        return



    def LatLonToIdx(self, vec_lat, vec_lon):
        import numpy as np
        import pyproj
        EASE_proj = pyproj.Proj(self.str_proj)
        
        dX = self.xc[0]-self.xc[1]
        dY = self.yc[0]-self.yc[1]

        vec_x, vec_y = EASE_proj(vec_lon, vec_lat)

        
        i = np.round((self.xc[0] - vec_x) / dX)
        j = np.round((self.yc[0] - vec_y) / dY)

        return i,j



    def GridData(self, tlim, vec_lat, vec_lon, t, SD=[], SD_unc=[], SIT=[], SIT_unc=[], FRB=[], FRB_unc=[]):
        #import pdb;
        import numpy as np
        import warnings
        import datetime as dt
        from datetime import timedelta
        import numpy.ma as ma
        
        if len(SD) == 0:
            SD = np.array(vec_lat) * np.nan

        if len(SD_unc) == 0:
            SD_unc = np.array(vec_lat) * np.nan
        
        if len(SIT) == 0:
            SIT = np.array(vec_lat) * np.nan

        if len(SIT_unc) == 0:
            SIT_unc = np.array(vec_lat) * np.nan

        if len(FRB) == 0:
            FRB = np.array(vec_lat) * np.nan

        if len(FRB_unc) == 0:
            FRB_unc = np.array(vec_lat) * np.nan            

        vec_i, vec_j = self.LatLonToIdx(vec_lat, vec_lon)

        # Tranform geographic lat/lon to Lambert Azimuthal Equal Area
        xx,yy = self.EASE_proj(vec_lon,vec_lat)
        # Transforms the 'laea' back to geographic lat/lon
        #lon,lat = self.EASE_proj(xx,yy,inverse=True)

        # Adds observations with same indicies to a list of list of lists
        # Removes data if time is not changing over successive measurements 
        count=0
        last = dt.datetime(2000,1,1,0,0,0)
        dt0 = timedelta(seconds=0)
        for i in range(len(t)):
            dt = t[i]-last
            if dt != dt0:
                count += 1
                # fill observations into grid
                self.sd_list[int(vec_i[i])][int(vec_j[i])].append(SD[i])
                self.sd_unc_list[int(vec_i[i])][int(vec_j[i])].append(SD_unc[i])
                self.sit_list[int(vec_i[i])][int(vec_j[i])].append(SIT[i])
                self.sit_unc_list[int(vec_i[i])][int(vec_j[i])].append(SIT_unc[i])
                self.Frb_list[int(vec_i[i])][int(vec_j[i])].append(FRB[i])
                self.Frb_unc_list[int(vec_i[i])][int(vec_j[i])].append(FRB_unc[i])
                self.td_list[int(vec_i[i])][int(vec_j[i])].append(t[i])
                self.lat_list[int(vec_i[i])][int(vec_j[i])].append(xx[i])
                self.lon_list[int(vec_i[i])][int(vec_j[i])].append(yy[i])
            last = t[i]
        print('Number of observations after removing dt=0: ',count)

        # Average observations for each grid cell of maximum time interval (tlim)               
        # output = open('test.dat','a')
        outxx   =[]
        outyy   =[]
        outtime  =[]
        # SD
        outavgsd =[]
        outstdsd =[]
        outlnsd  =[]
        outUncSD = []
        
        # FRB
        outFrb = []
        outstdFrb =[]
        outlnFrb =[]
        outUncFrb = []
        
        # SIT
        outSIT = []
        outstdSIT = []
        outlnSIT = []
        outUncSIT = []
        
        ## longitude
        for i in range(len(self.xc)):
            ## lattitude
            for j in range(len(self.yc)):
                ## enter if there are any data
                if self.sd_list[i][j] or self.sit_list[i][j] or self.Frb_list[i][j]: 
                    t0 = self.td_list[i][j][0]
                    tn = 1 
                    count = 0
                    
                    lat = []
                    lon = []
                    dt0 = []
                    
                    sd =  []
                    sd_unc = []
                    sit = []
                    sit_unc = []
                    Frb = []
                    frb_unc = []

                    for n in range(len(self.td_list[i][j])):
                        ## Add data if the time limit is less than 30 days
                        if (self.td_list[i][j][n] - t0).days < tlim:
                            sd=np.append(sd,self.sd_list[i][j][n])
                            sd_unc=np.append(sd_unc,self.sd_unc_list[i][j][n])
                            
                            sit=np.append(sit,self.sit_list[i][j][n])
                            sit_unc=np.append(sit_unc,self.sit_unc_list[i][j][n])
                            
                            Frb=np.append(Frb,self.Frb_list[i][j][n])
                            frb_unc=np.append(frb_unc,self.Frb_unc_list[i][j][n])
                            
                            lat=np.append(lat,self.lat_list[i][j][n])
                            lon=np.append(lon,self.lon_list[i][j][n])
                            dt0=np.append(dt0,(self.td_list[i][j][n] - t0).total_seconds())
                            count +=1
                            flag = 0
                        else:
                            ## Compute datapoint
                            count += 1
                            flag = 1
                            # supress warning - mean of empty numpy array
                            with warnings.catch_warnings():
                                warnings.simplefilter("ignore", category=RuntimeWarning)
                                avgSD = np.nanmean(sd)
                                stdSD = np.std(sd[np.isnan(sd)==0])
                                lnSD = len(sd[np.isnan(sd)==0])
                                avgSIT = np.nanmean(sit)
                                stdSIT = np.std(sit[np.isnan(sit)==0])
                                lnSIT = len(sit[np.isnan(sit)==0])
                                stdFrb = np.std(Frb[np.isnan(Frb)==0])
                                lnFrb = len(Frb[np.isnan(Frb)==0])
                                avgFrb = np.nanmean(Frb)
                                avgLat = np.mean(lat)
                                avgLon = np.mean(lon)
                                avgDT = np.mean(dt0)
                                avgTime = t0+timedelta(seconds=avgDT)  
                            # Calculation of uncertainty
                            if lnSD > 0:
                                sigma_SD = 1/lnSD*np.sqrt(np.nansum(sd_unc**2))
                            else:
                                sigma_SD = np.nan
                            if lnSIT > 0:
                                sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                            else:
                                sigma_SIT = np.nan
                            if lnFrb > 0:
                                sigma_Frb = 1/lnFrb*np.sqrt(np.nansum(frb_unc**2))
                            else:
                                sigma_Frb = np.nan			
				
                            # Append to output
                            outavgsd = np.append(outavgsd,avgSD) 
                            outstdsd = np.append(outstdsd,stdSD) 
                            outlnsd = np.append(outlnsd,lnSD) 
                            outFrb = np.append(outFrb,avgFrb) 
                            outstdFrb = np.append(outstdFrb,stdFrb)
                            outlnFrb = np.append(outlnFrb,lnFrb)
                            outSIT = np.append(outSIT,avgSIT)
                            outstdSIT = np.append(outstdSIT,stdSIT)
                            outlnSIT = np.append(outlnSIT,lnSIT)
                            outxx = np.append(outxx,avgLat) 
                            outyy = np.append(outyy,avgLon) 
                            outtime = np.append(outtime,avgTime)
                            # uncertainty
                            outUncSD = np.append(outUncSD,sigma_SD)
                            outUncSIT = np.append(outUncSIT,sigma_SIT)
                            outUncFrb = np.append(outUncFrb,sigma_Frb)

                            ## update time and append data
                            t0 = self.td_list[i][j][n]
                            self.n.append(n)
                            
                            sd =  []
                            sd_unc = []
                            sit = []
                            sit_unc = []
                            lat = []
                            lon = []
                            Frb = []
                            frb_unc = []
                            dt0 = []
                            
                            sd = np.append(sd,self.sd_list[i][j][n])
                            sd_unc = np.append(sd_unc,self.sd_unc_list[i][j][n])
                            Frb = np.append(Frb,self.Frb_list[i][j][n])
                            frb_unc = np.append(frb_unc,self.Frb_unc_list[i][j][n])
                            sit = np.append(sit,self.sit_list[i][j][n])
                            sit_unc = np.append(sit_unc,self.sit_unc_list[i][j][n])
                            lat = np.append(lat,self.lat_list[i][j][n])
                            lon = np.append(lon,self.lon_list[i][j][n])
                            dt0 = np.append(dt0,(self.td_list[i][j][n] - t0).total_seconds())

                            tn +=1

                            if count == len(self.td_list[i][j]): #done to include the last element which is otherwise excluded??
                                # supress warning - mean of empty numpy array
                                with warnings.catch_warnings():
                                    warnings.simplefilter("ignore", category=RuntimeWarning)
                                    avgSD = np.nanmean(sd)
                                    stdSD = np.std(sd[np.isnan(sd)==0])
                                    lnSD = len(sd[np.isnan(sd)==0])
                                    avgSIT = np.nanmean(sit)
                                    stdSIT = np.std(sit[np.isnan(sit)==0])
                                    lnSIT = len(sit[np.isnan(sit)==0])
                                    stdFrb = np.std(Frb[np.isnan(Frb)==0])
                                    lnFrb = len(Frb[np.isnan(Frb)==0])
                                    avgFrb = np.nanmean(Frb)
                                    avgLat = np.mean(lat)
                                    avgLon = np.mean(lon)
                                    avgDT = np.mean(dt0)
                                    avgTime = t0+timedelta(seconds=avgDT)                 

                                # Calculation of uncertainty
                                if lnSD > 0:
                                    sigma_SD = 1/lnSD*np.sqrt(np.nansum(sd_unc**2))
                                else:
                                    sigma_SD = np.nan
                                if lnSIT > 0:
                                    sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                                else:
                                    sigma_SIT = np.nan
                                if lnFrb > 0:
                                    sigma_Frb = 1/lnFrb*np.sqrt(np.nansum(frb_unc**2))
                                else:
                                    sigma_Frb = np.nan

                                outavgsd = np.append(outavgsd,avgSD) 
                                outstdsd = np.append(outstdsd,stdSD) 
                                outlnsd = np.append(outlnsd,lnSD)  
                                outFrb = np.append(outFrb,avgFrb)
                                outstdFrb = np.append(outstdFrb,stdFrb)
                                outlnFrb = np.append(outlnFrb,lnFrb)
                                outSIT = np.append(outSIT,avgSIT)
                                outstdSIT = np.append(outstdSIT,stdSIT)
                                outlnSIT = np.append(outlnSIT,lnSIT)
                                outxx = np.append(outxx,avgLat) 
                                outyy = np.append(outyy,avgLon) 
                                outtime = np.append(outtime,avgTime)
                                # uncertainty
                                outUncSD = np.append(outUncSD,sigma_SD)
                                outUncSIT = np.append(outUncSIT,sigma_SIT)
                                outUncFrb = np.append(outUncFrb,sigma_Frb)                                              
                    if flag == 0:
                        # supress warning - mean of empty numpy array
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=RuntimeWarning)
                            avgSD = np.nanmean(sd)
                            stdSD = np.std(sd[np.isnan(sd)==0])
                            lnSD = len(sd[np.isnan(sd)==0])
                            avgSIT = np.nanmean(sit)
                            stdSIT = np.std(sit[np.isnan(sit)==0])
                            lnSIT = len(sit[np.isnan(sit)==0])
                            stdFrb = np.std(Frb[np.isnan(Frb)==0])
                            lnFrb = len(Frb[np.isnan(Frb)==0])
                            avgFrb = np.nanmean(Frb)
                            avgLat = np.mean(lat)
                            avgLon = np.mean(lon)
                            avgDT = np.mean(dt0)
                            avgTime = t0+timedelta(seconds=avgDT)             

                        # Calculation of uncertainty
                        if lnSD > 0:
                            sigma_SD = 1/lnSD*np.sqrt(np.nansum(sd_unc**2))
                        else:
                            sigma_SD = np.nan
                        if lnSIT > 0:
                            sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                        else:
                            sigma_SIT = np.nan
                        if lnFrb > 0:
                            sigma_Frb = 1/lnFrb*np.sqrt(np.nansum(frb_unc**2))
                        else:
                            sigma_Frb = np.nan
        
                        outavgsd = np.append(outavgsd,avgSD) 
                        outstdsd = np.append(outstdsd,stdSD) 
                        outlnsd = np.append(outlnsd,lnSD) 
                        outFrb = np.append(outFrb,avgFrb) 
                        outstdFrb = np.append(outstdFrb,stdFrb)
                        outlnFrb = np.append(outlnFrb,lnFrb)
                        outSIT = np.append(outSIT,avgSIT)
                        outstdSIT = np.append(outstdSIT,stdSIT)
                        outlnSIT = np.append(outlnSIT,lnSIT)
                        outxx = np.append(outxx,avgLat) 
                        outyy = np.append(outyy,avgLon) 
                        outtime = np.append(outtime,avgTime)
                        # uncertainty
                        outUncSD = np.append(outUncSD,sigma_SD)
                        outUncSIT = np.append(outUncSIT,sigma_SIT)
                        outUncFrb = np.append(outUncFrb,sigma_Frb)
                          
       
        outlon,outlat = self.EASE_proj(outxx,outyy,inverse=True)

        return outavgsd, outstdsd, outlnsd, outUncSD, outlat, outlon, outtime, outSIT, outstdSIT, outlnSIT, outUncSIT, outFrb, outstdFrb, outlnFrb, outUncFrb

         

    def convert_to_partial_year(self,vec_time=[]):

        import calendar
        from datetime import datetime
        from datetime import timedelta

        if calendar.isleap(vec_time.year):
            dayinyear = 366
        else:
            dayinyear = 365

        d_dayinyear = timedelta(days=dayinyear)
        day_one = datetime(vec_time.year,1,1)
        d = vec_time - day_one

        dec_year = vec_time.year + d.total_seconds()/d_dayinyear.total_seconds()

        return dec_year





            

