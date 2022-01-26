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

        # Default for 25 km EASE grid! https://nsidc.org/ease/ease-grid-projection-gt
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
        SIT_list = list()
        sit_unc_list = list()
        tsur_list = list()
        lat_list = list()
        lon_list = list()
        tn_list = list()
        td_list = list()

        
        #Create 2D empty list of lists corresponding to the grid-size
        for i in range(len(self.xc)):
            sd_list.append(list())
            SIT_list.append(list())
            sit_unc_list.append(list())
            tsur_list.append(list())
            lon_list.append(list())
            lat_list.append(list())
            tn_list.append(list())
            td_list.append(list())
            for j in range(len(self.yc)):
                sd_list[i].append(list())
                SIT_list[i].append(list())
                sit_unc_list[i].append(list())
                tsur_list[i].append(list())
                lat_list[i].append(list())
                lon_list[i].append(list())
                td_list[i].append(list())
                tn_list[i].append(list())
        self.sd_list = sd_list
        self.SIT_list = SIT_list
        self.sit_unc_list = sit_unc_list
        self.tsur_list = tsur_list
        self.lon_list = lon_list
        self.lat_list = lat_list
        self.td_list = tn_list
        self.tn_list = tn_list
      
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



    def GridData(self, tlim, vec_lat, vec_lon, SIT=[], SIT_unc=[], t=[], snow_depth=[], SD_unc=[],  sur_temp=[]):
        #import pdb;
        import numpy as np
        import datetime as dt
        from datetime import timedelta
        import numpy.ma as ma
        
        if len(snow_depth) == 0:
            snow_depth = np.array(vec_lat) * np.nan

        if len(SD_unc) == 0:
            SD_unc = np.array(vec_lat) * np.nan
        
        if len(SIT) == 0:
            SIT = np.array(vec_lat) * np.nan

        if len(SIT_unc) == 0:
            SIT_unc = np.array(vec_lat) * np.nan

        if len(t) == 0:
            t = np.array(vec_lat) * np.nan

        if len(sur_temp) == 0:
            sur_temp = np.array(vec_lat) * np.nan
            

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
        for i in range(len(snow_depth)):
            dt = t[i]-last
            if dt != dt0:
                count += 1
                self.sd_list[int(vec_i[i])][int(vec_j[i])].append(snow_depth[i])
                self.SIT_list[int(vec_i[i])][int(vec_j[i])].append(SIT[i])
                self.sit_unc_list[int(vec_i[i])][int(vec_j[i])].append(SIT_unc[i])
                self.tsur_list[int(vec_i[i])][int(vec_j[i])].append(sur_temp[i])
                self.td_list[int(vec_i[i])][int(vec_j[i])].append(t[i])
                self.lat_list[int(vec_i[i])][int(vec_j[i])].append(xx[i])
                self.lon_list[int(vec_i[i])][int(vec_j[i])].append(yy[i])
            last = t[i]
            #
        print('Number of observations after removing dt=0: ',count)

        # Average observations for each grid cell of maximum time interval (tlim)               
        output=open('test.dat','a')
        outxx   =[]
        outyy   =[]
        outtime  =[]
        outavgsd =[]
        outstdsd =[]
        outlnsd  =[]
        outtsur = []
        outstdfrb =[]
        outlnfrb =[]
        outSIT = []
        outstdSIT = []
        outlnSIT = []
        outUncSIT = []
        for i in range(len(self.xc)):
            for j in range(len(self.yc)):
                if self.sd_list[i][j]: 
                    t0 = self.td_list[i][j][0]
                    tn = 1 
                    count = 0
		    
                    sd =  []
                    sit = []
                    sit_unc = []
                    lat = []
                    lon = []
                    tsur = []
                    dt0 = []	
	    
                    for n in range(len(self.td_list[i][j])):
                        if (self.td_list[i][j][n] - t0).days < tlim:
                            sd=np.append(sd,self.sd_list[i][j][n])
                            sit=np.append(sit,self.SIT_list[i][j][n])
                            sit_unc = np.append(sit_unc,self.sit_unc_list[i][j][n])
                            tsur=np.append(tsur,self.tsur_list[i][j][n])
                            lat=np.append(lat,self.lat_list[i][j][n])
                            lon=np.append(lon,self.lon_list[i][j][n])
                            dt0=np.append(dt0,(self.td_list[i][j][n] - t0).total_seconds())
                            count +=1
                            flag=0
                            print("%6i %6i %6i %6i %6.2f %6.2f %6.2f %6.2f %12.6f %12.0f" % (i,j,tn,(self.td_list[i][j][n] - t0).days,self.sd_list[i][j][n],self.SIT_list[i][j][n],self.tsur_list[i][j][n],self.lat_list[i][j][n],self.lon_list[i][j][n],(self.td_list[i][j][n] - t0).total_seconds()),output)
                        else:
                            count +=1
                            flag=1
                            avgSD=np.nanmean(sd)
                            avgSIT=np.nanmean(sit)
                            if np.isnan(avgSIT):
                                print('nan')
                            else:
                                avgSD=np.nanmean(sd)
                                stdSD=np.std(sd[np.isnan(sd)==0])
                                lnSD=len(sd[np.isnan(sd)==0])
                                avgSIT=np.nanmean(sit)
                                stdSIT=np.std(sit[np.isnan(sit)==0])
                                lnSIT=len(sit[np.isnan(sit)==0])
                                stdFrb=np.std(tsur[np.isnan(tsur)==0])
                                lnFrb=len(tsur[np.isnan(tsur)==0])
                                avgTsur=np.nanmean(tsur)
                                avgLat=np.mean(lat[np.isnan(lat)==0])
                                avgLon=np.mean(lon[np.isnan(lon)==0])
                                avgDT = np.mean(dt0[np.isnan(sit)==0])
                                avgTime = t0+timedelta(seconds=avgDT)
                                
                                # calculation of uncertainty
                                if lnSIT > 0:
                                    sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                                else:
                                    sigma_SIT = np.nan
    
                                #flag = 1
                                outavgsd=np.append(outavgsd,avgSD) 
                                outstdsd=np.append(outstdsd,stdSD) 
                                outlnsd=np.append(outlnsd,lnSD) 
                                outtsur=np.append(outtsur,avgTsur) 
                                outstdfrb=np.append(outstdfrb,stdFrb)
                                outlnfrb=np.append(outlnfrb,lnFrb)
                                outSIT=np.append(outSIT,avgSIT)
                                outstdSIT=np.append(outstdSIT,stdSIT)
                                outlnSIT = np.append(outlnSIT,lnSIT)
                                outUncSIT = np.append(outUncSIT,sigma_SIT)
                                outxx=np.append(outxx,avgLat) 
                                outyy=np.append(outyy,avgLon) 
                                outtime=np.append(outtime,avgTime)
    
                                print((i,j,tn,avgSD,stdSD,lnSD,avgTsur,stdFrb,lnFrb,avgSIT,stdSIT,lnSIT,avgLat,avgLon,avgDT),file=output)
                            t0 = self.td_list[i][j][n]
                            sd  = [] 
                            lat = []
                            lon = []
                            tsur = []
                            sit = []
                            dt0 = []
                            sd=np.append(sd,self.sd_list[i][j][n])
                            tsur=np.append(sd,self.tsur_list[i][j][n])
                            sit=np.append(sd,self.SIT_list[i][j][n])
                            sit_unc = np.append(sit_unc,self.sit_unc_list[i][j][n])
                            lat=np.append(lat,self.lat_list[i][j][n])
                            lon=np.append(lon,self.lon_list[i][j][n])
                            dt0=np.append(dt0,(self.td_list[i][j][n] - t0).total_seconds())

                            tn +=1
                            print("%6i %6i %6i %6i %6.2f %6.2f %6.2f %6.2f %12.6f %12.6f %12.0f" % (i,j,tn,(self.td_list[i][j][n] - t0).days,self.sd_list[i][j][n],self.SIT_list[i][j][n],self.tsur_list[i][j][n],self.lat_list[i][j][n],self.lon_list[i][j][n],(self.td_list[i][j][n] - t0).total_seconds()),output) 
                            if count == len(self.td_list[i][j]): #done to include the last element which is otherwise excluded??
                                avgSD=np.nanmean(sd)
                                avgSIT=np.nanmean(sit)
                                if np.isnan(avgSIT):
                                    print('nan')
                                else:
                                    stdSIT=np.std(sit[np.isnan(sit)==0])
                                    lnSIT=len(sit[np.isnan(sit)==0])
                                    stdSD=np.std(sd[np.isnan(sd)==0])
                                    lnSD=len(sd[np.isnan(sd)==0]) 
                                    avgTsur=np.nanmean(tsur)
                                    stdFrb=np.std(tsur[np.isnan(tsur)==0])
                                    lnFrb=len(tsur[np.isnan(tsur)==0])
                                    avgLat=np.mean(lat[np.isnan(lat)==0])
                                    avgLon=np.mean(lon[np.isnan(lon)==0])
                                    avgDT = np.mean(dt0[np.isnan(sit)==0])
                                    avgTime = t0+timedelta(seconds=avgDT)                 
                                    # calculation of uncertainty
                                    if lnSIT > 0:
                                        sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                                    else:
                                        sigma_SIT = np.nan
        
                                    #flag = 1
                                    outavgsd=np.append(outavgsd,avgSD) 
                                    outstdsd=np.append(outstdsd,stdSD) 
                                    outlnsd=np.append(outlnsd,lnSD) 
                                    outtsur=np.append(outtsur,avgTsur) 
                                    outstdfrb=np.append(outstdfrb,stdFrb)
                                    outlnfrb=np.append(outlnfrb,lnFrb)
                                    outSIT=np.append(outSIT,avgSIT)
                                    outstdSIT=np.append(outstdSIT,stdSIT)
                                    outlnSIT = np.append(outlnSIT,lnSIT)
                                    outUncSIT = np.append(outUncSIT,sigma_SIT)
                                    outxx=np.append(outxx,avgLat) 
                                    outyy=np.append(outyy,avgLon) 
                                    outtime=np.append(outtime,avgTime)
                                   
                                    print(i,j,tn,avgSD,stdSD,lnSD,avgTsur,stdFrb,lnFrb,avgSIT,stdSIT,lnSIT,avgLat,avgLon,avgDT,output)
                                    #print(i,j,tn,avgSD,stdSD,lnSD,avgTair,avgLat,avgLon,avgDT,output)

                    if flag == 0:
                        
                        avgSD=np.nanmean(sd)
                        avgSIT=np.nanmean(sit)
                        if np.isnan(avgSIT):
                            print('nan')
                        else:
                            stdSIT=np.std(sit[np.isnan(sit)==0])
                            lnSIT=len(sit[np.isnan(sit)==0])
                            stdSD=np.std(sd[np.isnan(sd)==0])
                            lnSD=len(sd[np.isnan(sd)==0])
                            avgTsur=np.nanmean(tsur)
                            stdFrb=np.std(tsur[np.isnan(tsur)==0])
                            lnFrb=len(tsur[np.isnan(tsur)==0])
                            avgLat=np.mean(lat[np.isnan(lat)==0])
                            avgLon=np.mean(lon[np.isnan(lon)==0])
                            avgDT = np.mean(dt0[np.isnan(sit)==0])
                            avgTime = t0+timedelta(seconds=avgDT)                 
                            # calculation of uncertainty
                            if lnSIT > 0:
                                sigma_SIT = 1/lnSIT*np.sqrt(np.nansum(sit_unc**2))
                                # sigma_SIT = 1/lnSIT*np.nansum(sit_unc) upper bound uncertainty
                            else:
                                sigma_SIT = np.nan
    
                            #flag = 1
                            outavgsd=np.append(outavgsd,avgSD) 
                            outstdsd=np.append(outstdsd,stdSD) 
                            outlnsd=np.append(outlnsd,lnSD) 
                            outtsur=np.append(outtsur,avgTsur) 
                            outstdfrb=np.append(outstdfrb,stdFrb)
                            outlnfrb=np.append(outlnfrb,lnFrb)
                            outSIT=np.append(outSIT,avgSIT)
                            outstdSIT=np.append(outstdSIT,stdSIT)
                            outlnSIT = np.append(outlnSIT,lnSIT)
                            outUncSIT = np.append(outUncSIT,sigma_SIT)
                            outxx=np.append(outxx,avgLat) 
                            outyy=np.append(outyy,avgLon) 
                            outtime=np.append(outtime,avgTime)
                    
                            print(i,j,tn,avgSD,stdSD,lnSD,avgTsur,avgSIT,avgLat,avgLon,avgDT,output)

                          
        output.close()
       
        outlon,outlat = self.EASE_proj(outxx,outyy,inverse=True)

        return outlat,outlon,outtime,outSIT,outstdSIT,outlnSIT, outUncSIT

         

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





            

