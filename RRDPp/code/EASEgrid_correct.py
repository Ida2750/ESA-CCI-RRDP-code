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

        # Default for 25 km EASE grid
        self.max_x = 9000000.0	#5400000.0
        self.min_x = -9000000.0	#-5400000.0

        self.max_y = 9000000.0	#5400000.0
        self.min_y = -9000000.0	#-5400000.0


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
        var1_list = list()
        var2_list = list()
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
            var1_list.append(list())
            var2_list.append(list())
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
                var1_list[i].append(list())
                var2_list[i].append(list())
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
        self.var1_list = var1_list
        self.var2_list = var2_list
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


    def robust_std(self,data):
        import numpy as np
        median = np.nanmedian(data)
        mad = np.nanmedian(np.abs(data - median))
        return mad * 1.4826  # scaling factor for normal distribution


    def GridData(self, tlim, vec_lat, vec_lon, t, SD=[], SD_unc=[], SIT=[], SIT_unc=[], FRB=[], FRB_unc=[], VAR1=[], VAR2=[], dtype='buoy'):
        #import pdb;
        import numpy as np
        import warnings
        import datetime as dt
        from datetime import timedelta
        import matplotlib.pyplot as plt
        import numpy.ma as ma
        
        plotcounter = 0
        name = 'SCICEX'
        fig, ax = plt.subplots(nrows=4, ncols=4, figsize=(12, 12), constrained_layout=True)
        ax = ax.flatten()

        # for spatial QF - find average frequency of measurements
        # freq = np.mean(np.diff(t)).total_seconds()
        # print(f'freq: {1/freq}')

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
        
        # e.g. air temperature
        if len(VAR1) == 0:
            VAR1 = np.array(vec_lat) * np.nan
        
        # e.g. surface temperature
        if len(VAR2) == 0:
            VAR2 = np.array(vec_lat) * np.nan


        # replace unc values with nan where data is nan
        SD_unc[np.isnan(SD)] = np.nan
        SIT_unc[np.isnan(SIT)] = np.nan
        FRB_unc[np.isnan(FRB)] = np.nan

        vec_i, vec_j = self.LatLonToIdx(vec_lat, vec_lon)
        #print(vec_i, vec_j)

        # Tranform geographic lat/lon to Lambert Azimuthal Equal Area
        xx,yy = self.EASE_proj(vec_lon,vec_lat)
        # Transforms the 'laea' back to geographic lat/lon
        #lon,lat = self.EASE_proj(xx,yy,inverse=True)

        # Adds observations with same indicies to a list of list of lists
        # Removes data if time is not changing over successive measurements 
        count=0
        last = dt.datetime(1950,1,1,0,0,0)
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
                self.var1_list[int(vec_i[i])][int(vec_j[i])].append(VAR1[i])
                self.var2_list[int(vec_i[i])][int(vec_j[i])].append(VAR2[i])
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
        
        # Additional variables e.g. surface temperature/air temperature
        outVAR1 = []
        outVAR2 = []

        # QF (quality flags)
        outQFT = []
        outQFS = []
        outQFG = []
        
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
                    var1 = []
                    var2 = []

                    for n in range(len(self.td_list[i][j])):
                        ## Add data if the time limit is less than 30 days
                        if (self.td_list[i][j][n] - t0).days < tlim:
                            sd=np.append(sd,self.sd_list[i][j][n])
                            sd_unc=np.append(sd_unc,self.sd_unc_list[i][j][n])
                            
                            sit=np.append(sit,self.sit_list[i][j][n])
                            sit_unc=np.append(sit_unc,self.sit_unc_list[i][j][n])
                            
                            Frb=np.append(Frb,self.Frb_list[i][j][n])
                            frb_unc=np.append(frb_unc,self.Frb_unc_list[i][j][n])
                            
                            var1 = np.append(var1,self.var1_list[i][j][n])
                            var2 = np.append(var2,self.var2_list[i][j][n])
                            
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
                                avgSD = np.nanmedian(sd)
                                stdSD = self.robust_std(sd[np.isnan(sd)==0])
                                #print(f'robust std: {stdSD}\n')
                                #print(f' std: {np.std(sd[np.isnan(sd)==0])}\n')
                                lnSD = len(sd[np.isnan(sd)==0])
                                avgSIT = np.nanmedian(sit)
                                stdSIT = self.robust_std(sit[np.isnan(sit)==0])
                                lnSIT = len(sit[np.isnan(sit)==0])
                                stdFrb = self.robust_std(Frb[np.isnan(Frb)==0])
                                lnFrb = len(Frb[np.isnan(Frb)==0])
                                avgFrb = np.nanmedian(Frb)
                                avgVAR1 = np.nanmedian(var1)
                                avgVar2 = np.nanmedian(var2)
                                avgLat = np.mean(lat)
                                avgLon = np.mean(lon)
                                avgDT = np.mean(dt0)
                                avgTime = t0+timedelta(seconds=avgDT)

                            # determine quality flags

                            ################# TEMPORAL FLAG START ##################
                            # compute time between first and last measuement
                            # timediff = sorted(dt0)[-1] - sorted(dt0)[0]
                            # compute number of days in gridcell
                            tny = [(t0+timedelta(seconds=t)).day for t in dt0]
                            # check number of unique days in gridcell
                            unique = len(np.unique(tny))

                            if unique==1:
                                QF_flag_temporal = 3
                            elif unique<=5:
                                QF_flag_temporal = 2
                            elif unique<15:
                                QF_flag_temporal = 1
                            elif unique>=15:
                                QF_flag_temporal = 0

                            ################# TEMPORAL FLAG END ##################
                            ################# SPATIAL FLAG START ##################
                            # footprint
                            if dtype == 'AEM':
                                footprint = 45 # meters source: https://epic.awi.de/id/eprint/50023/1/IceBird-2019-Winter-ICESat2DataAcquisitionReport.pdf
                            elif dtype == 'OIB':
                                footprint = 1 # meters source: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020RG000712
                            # elif dtype == 'ship':
                            #     footprint = 2000 # meters source: https://fiona.uni-hamburg.de/bdc0b40b/hutchings-ice-watch-manual-v4-1.pdf
                            elif dtype == 'buoy' or dtype == 'ship':
                                QF_flag_spatial = 3
                                footprint = np.nan
                            elif dtype == 'submarine': # diameter between 2.6 and 6 meters source: https://journals.ametsoc.org/view/journals/atot/24/11/jtech2097_1.pdf
                                footprint = (6+2.6)/2


                            if np.isfinite(footprint):
                                # compute area - upper limit of area assuming no overlap between measurements which is not true!
                                # area = pi* r²
                                area = 3.14*(footprint/2)**2 *lnSIT

                                if area/(self.resolution**2)*100 > 1: # what is a sufficiently high area covered?
                                    QF_flag_spatial = 0
                                    #print('good data')
                                    #print(f'lnSIT {lnSIT}')
                                else:
                                    QF_flag_spatial = 1

                            ################# SPATIAL FLAG END ##################
                            ################# GBL THRESHOLD FLAG START ################
                            QF_flag_threshold = 0 # default
                            if np.any(sit>8):
                                QF_flag_threshold = 1
                            if np.any(sd>2) and dtype!='submarine':
                                QF_flag_threshold = 1
                            if np.any(sd>6) and dtype=='submarine': # sd is a placeholder for SID
                                QF_flag_threshold = 1
                            if np.any(Frb>3):
                                QF_flag_threshold = 1
                            ################# GBL THRESHOLD FLAG END ##################

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
                            outVAR1 = np.append(outVAR1,avgVAR1)
                            outVAR2 = np.append(outVAR2,avgVar2)
                            outxx = np.append(outxx,avgLat) 
                            outyy = np.append(outyy,avgLon) 
                            outtime = np.append(outtime,avgTime)
                            # uncertainty
                            outUncSD = np.append(outUncSD,sigma_SD)
                            outUncSIT = np.append(outUncSIT,sigma_SIT)
                            outUncFrb = np.append(outUncFrb,sigma_Frb)
                            # QF
                            outQFT = np.append(outQFT,QF_flag_temporal)
                            outQFS = np.append(outQFS,QF_flag_spatial)
                            outQFG = np.append(outQFG,QF_flag_threshold)  

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
                            var1 = []
                            var2 = []
                            dt0 = []
                            
                            sd = np.append(sd,self.sd_list[i][j][n])
                            sd_unc = np.append(sd_unc,self.sd_unc_list[i][j][n])
                            Frb = np.append(Frb,self.Frb_list[i][j][n])
                            frb_unc = np.append(frb_unc,self.Frb_unc_list[i][j][n])
                            sit = np.append(sit,self.sit_list[i][j][n])
                            sit_unc = np.append(sit_unc,self.sit_unc_list[i][j][n])
                            var1 = np.append(var1,self.var1_list[i][j][n])
                            var2 = np.append(var2,self.var2_list[i][j][n])
                            lat = np.append(lat,self.lat_list[i][j][n])
                            lon = np.append(lon,self.lon_list[i][j][n])
                            dt0 = np.append(dt0,(self.td_list[i][j][n] - t0).total_seconds())

                            tn +=1

                            if count == len(self.td_list[i][j]): #done to include the last element which is otherwise excluded??
                                # supress warning - mean of empty numpy array
                                with warnings.catch_warnings():
                                    warnings.simplefilter("ignore", category=RuntimeWarning)
                                    avgSD = np.nanmedian(sd)
                                    stdSD = self.robust_std(sd[np.isnan(sd)==0])
                                    #print(f'robust std: {stdSD}\n')
                                    #print(f' std: {np.std(sd[np.isnan(sd)==0])}\n')
                                    lnSD = len(sd[np.isnan(sd)==0])
                                    avgSIT = np.nanmedian(sit)
                                    stdSIT = self.robust_std(sit[np.isnan(sit)==0])
                                    lnSIT = len(sit[np.isnan(sit)==0])
                                    stdFrb = self.robust_std(Frb[np.isnan(Frb)==0])
                                    lnFrb = len(Frb[np.isnan(Frb)==0])
                                    avgFrb = np.nanmedian(Frb)
                                    avgVAR1 = np.nanmedian(var1)
                                    avgVar2 = np.nanmedian(var2)
                                    avgLat = np.mean(lat)
                                    avgLon = np.mean(lon)
                                    avgDT = np.mean(dt0)
                                    avgTime = t0+timedelta(seconds=avgDT) 

                                    # ### plot distributions
                                    # plotcounter += 1
                                    # if any(np.isfinite(sd)):
                                    #     plt.figure()
                                    #     plt.hist(sd, bins=50)
                                    #     plt.title(f'Median: {np.nanmedian(sd)}, Mean: {np.nanmean(sd)}')
                                    #     plt.savefig(f'sd_test_{name}_{plotcounter}.png')   
                                    #     plt.close()

                                    # if any(np.isfinite(sit)):
                                    #     plt.figure()
                                    #     plt.hist(sit, bins=50)
                                    #     plt.title(f'Median: {np.nanmedian(sit)}, Mean: {np.nanmean(sit)}')
                                    #     plt.savefig(f'sit_test_{name}_{plotcounter}.png')   
                                    #     plt.close()
                                    ### ------------ ###



                                # determine quality flags

                                ################# TEMPORAL FLAG START ##################
                                # compute time between first and last measuement
                                # timediff = sorted(dt0)[-1] - sorted(dt0)[0]
                                # compute number of days in gridcell
                                tny = [(t0+timedelta(seconds=t)).day for t in dt0]
                                # check number of unique days in gridcell
                                unique = len(np.unique(tny))

                                if unique==1:
                                    QF_flag_temporal = 3
                                elif unique<=5:
                                    QF_flag_temporal = 2
                                elif unique<15:
                                    QF_flag_temporal = 1
                                elif unique>=15:
                                    QF_flag_temporal = 0

                                ################# TEMPORAL FLAG END ##################
                                ################# SPATIAL FLAG START ##################
                                # footprint
                                if dtype == 'AEM':
                                    footprint = 45 # meters source: https://epic.awi.de/id/eprint/50023/1/IceBird-2019-Winter-ICESat2DataAcquisitionReport.pdf
                                elif dtype == 'OIB':
                                    footprint = 1 # meters source: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020RG000712
                                # elif dtype == 'ship':
                                #     footprint = 2000 # meters source: https://fiona.uni-hamburg.de/bdc0b40b/hutchings-ice-watch-manual-v4-1.pdf
                                elif dtype == 'buoy' or dtype == 'ship':
                                    QF_flag_spatial = 3
                                    footprint = np.nan
                                elif dtype == 'submarine': # diameter between 2.6 and 6 meters source: https://journals.ametsoc.org/view/journals/atot/24/11/jtech2097_1.pdf
                                    footprint = (6+2.6)/2

                                if np.isfinite(footprint):
                                    # compute area - upper limit of area assuming no overlap between measurements which is not true!
                                    # area = pi* r²
                                    area = 3.14*(footprint/2)**2 *lnSIT

                                    if area/(self.resolution**2)*100 > 1: # what is a sufficiently high area covered?
                                        QF_flag_spatial = 0
                                        #print('good data')
                                        #print(f'lnSIT {lnSIT}')
                                    else:
                                        QF_flag_spatial = 1

                                ################# SPATIAL FLAG END ##################
                                ################# GBL THRESHOLD FLAG START ################
                                QF_flag_threshold = 0 # default
                                if np.any(sit>8):
                                    QF_flag_threshold = 1
                                if np.any(sd>2) and dtype!='submarine':
                                    QF_flag_threshold = 1
                                if np.any(sd>6) and dtype=='submarine': # sd is a placeholder for SID
                                    QF_flag_threshold = 1
                                if np.any(Frb>3):
                                    QF_flag_threshold = 1
                                ################# GBL THRESHOLD FLAG END ##################

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
                                outVAR1 = np.append(outVAR1,avgVAR1)
                                outVAR2 = np.append(outVAR2,avgVar2)
                                outxx = np.append(outxx,avgLat) 
                                outyy = np.append(outyy,avgLon) 
                                outtime = np.append(outtime,avgTime)
                                # uncertainty
                                outUncSD = np.append(outUncSD,sigma_SD)
                                outUncSIT = np.append(outUncSIT,sigma_SIT)
                                outUncFrb = np.append(outUncFrb,sigma_Frb) 
                                # QF
                                outQFT = np.append(outQFT,QF_flag_temporal)
                                outQFS = np.append(outQFS,QF_flag_spatial)
                                outQFG = np.append(outQFG,QF_flag_threshold)  

                    if flag == 0:
                        # supress warning - mean of empty numpy array
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=RuntimeWarning)
                            avgSD = np.nanmedian(sd)
                            stdSD = self.robust_std(sd[np.isnan(sd)==0])
                            #print(f'robust std: {stdSD}\n')
                            #print(f' std: {np.std(sd[np.isnan(sd)==0])}\n')
                            lnSD = len(sd[np.isnan(sd)==0])
                            avgSIT = np.nanmedian(sit)
                            stdSIT = self.robust_std(sit[np.isnan(sit)==0])
                            lnSIT = len(sit[np.isnan(sit)==0])
                            stdFrb = self.robust_std(Frb[np.isnan(Frb)==0])
                            lnFrb = len(Frb[np.isnan(Frb)==0])
                            avgFrb = np.nanmedian(Frb)
                            avgVAR1 = np.nanmedian(var1)
                            avgVar2 = np.nanmedian(var2)
                            avgLat = np.mean(lat)
                            avgLon = np.mean(lon)
                            avgDT = np.mean(dt0)
                            avgTime = t0+timedelta(seconds=avgDT)    


                            ### plot distributions
                            if any(np.isfinite(sit)):
                                if plotcounter<16:
                                    median_val = np.nanmedian(sit)
                                    mean_val   = np.nanmean(sit)
                                    ax[plotcounter].hist(sit, bins=50, label=f'Median: {median_val:.2f}, Mean: {mean_val:.2f}')
                                    ax[plotcounter].set_xlabel('SIT [m]')
                                    ax[plotcounter].set_ylabel('count')
                                    ax[plotcounter].legend()
                                
                                if plotcounter == 16:
                                    fig.suptitle(f'{name}: Example of gricell distribution of data')
                                    plt.savefig(f'sit_test_{name}.png')   
                                    plt.close()

                                plotcounter += 1


                            ### plot distributions
                            if any(np.isfinite(sd)):
                                if plotcounter<16:
                                    median_val = np.nanmedian(sd)
                                    mean_val   = np.nanmean(sd)
                                    ax[plotcounter].hist(sd, bins=50, label=f'Median: {median_val:.2f}, Mean: {mean_val:.2f}')
                                    ax[plotcounter].set_xlabel('SID [m]')
                                    ax[plotcounter].set_ylabel('count')
                                    ax[plotcounter].legend()
                                
                                if plotcounter == 16:
                                    fig.suptitle(f'{name}: Example of gricell distribution of data')
                                    plt.savefig(f'sid_test_{name}.png')   
                                    plt.close()

                                plotcounter += 1


                            # if any(np.isfinite(sd)):
                            #     plt.figure()
                            #     plt.hist(sd, bins=50)
                            #     plt.title(f'Median: {np.nanmedian(sd)}, Mean: {np.nanmean(sd)}')
                            #     plt.savefig(f'sd_test_{name}_{plotcounter}.png')   
                            #     plt.close()

                            # if any(np.isfinite(sit)):
                            #     plt.figure()
                            #     plt.hist(sit, bins=50)
                            #     plt.title(f'Median: {np.nanmedian(sit)}, Mean: {np.nanmean(sit)}')
                            #     plt.savefig(f'sit_test_{name}_{plotcounter}.png')   
                            #     plt.close()
                            ### ------------ ###

                        # determine quality flags
                        ################# TEMPORAL FLAG START ##################
                        # compute time between first and last measuement
                        # timediff = sorted(dt0)[-1] - sorted(dt0)[0]
                        # compute number of days in gridcell
                        tny = [(t0+timedelta(seconds=t)).day for t in dt0]
                        # check number of unique days in gridcell
                        unique = len(np.unique(tny))

                        if unique==1:
                            QF_flag_temporal = 3
                        elif unique<=5:
                            QF_flag_temporal = 2
                        elif unique<15:
                            QF_flag_temporal = 1
                        elif unique>=15:
                            QF_flag_temporal = 0

                        ################# TEMPORAL FLAG END ##################
                        ################# SPATIAL FLAG START ##################
                        # footprint
                        if dtype == 'AEM':
                            footprint = 45 # meters source: https://epic.awi.de/id/eprint/50023/1/IceBird-2019-Winter-ICESat2DataAcquisitionReport.pdf
                        elif dtype == 'OIB':
                            footprint = 1 # meters source: https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020RG000712
                        # elif dtype == 'ship':
                        #     footprint = 2000 # meters source: https://fiona.uni-hamburg.de/bdc0b40b/hutchings-ice-watch-manual-v4-1.pdf
                        elif dtype == 'buoy' or dtype == 'ship':
                            QF_flag_spatial = 3
                            footprint = np.nan
                        elif dtype == 'submarine': # diameter between 2.6 and 6 meters source: https://journals.ametsoc.org/view/journals/atot/24/11/jtech2097_1.pdf
                            footprint = (6+2.6)/2

                        if np.isfinite(footprint):
                            # compute area - upper limit of area assuming no overlap between measurements which is not true!
                            # area = pi* r²
                            area = 3.14*(footprint/2)**2 *lnSIT

                            if area/(self.resolution**2)*100 > 1: # what is a sufficiently high area covered?
                                QF_flag_spatial = 0
                                #print('good data')
                                #print(f'lnSIT {lnSIT}')
                            else:
                                QF_flag_spatial = 1

                        ################# SPATIAL FLAG END ##################
                        ################# GBL THRESHOLD FLAG START ################
                        QF_flag_threshold = 0 # default
                        if np.any(sit>8):
                            QF_flag_threshold = 1
                        if np.any(sd>2) and dtype!='submarine':
                            QF_flag_threshold = 1
                        if np.any(sd>6) and dtype=='submarine': # sd is a placeholder for SID
                            QF_flag_threshold = 1
                        if np.any(Frb>3):
                            QF_flag_threshold = 1
                        ################# GBL THRESHOLD FLAG END ##################

                        # if freq > 1/(5*60) and (lnSD>1440 or lnSIT>1440 or lnFrb>1440): # more than once every 5 minutes
                        #     QF_flag_spatial = 0
                        # elif (freq > 1/(0.5*3600) and freq < 1/(5*60)) and (lnSD>240 | lnSIT>240 | lnFrb>240):
                        #     QF_flag_spatial = 0
                        # elif freq < 1/(0.5*3600) and (lnSD>30 | lnSIT>30 | lnFrb>30): # once every 0.5 hour
                        #     QF_flag_spatial = 0
                        # else:
                        #     QF_flag_spatial = 1
                        # print(QF_flag_spatial)

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
                        outVAR1 = np.append(outVAR1,avgVAR1)
                        outVAR2 = np.append(outVAR2,avgVar2)
                        outxx = np.append(outxx,avgLat) 
                        outyy = np.append(outyy,avgLon) 
                        outtime = np.append(outtime,avgTime)
                        # uncertainty
                        outUncSD = np.append(outUncSD,sigma_SD)
                        outUncSIT = np.append(outUncSIT,sigma_SIT)
                        outUncFrb = np.append(outUncFrb,sigma_Frb)
                        # QF
                        outQFT = np.append(outQFT,QF_flag_temporal)
                        outQFS = np.append(outQFS,QF_flag_spatial)
                        outQFG = np.append(outQFG,QF_flag_threshold)
                          
       
        outlon,outlat = self.EASE_proj(outxx,outyy,inverse=True)

        return outavgsd, outstdsd, outlnsd, outUncSD, outlat, outlon, outtime, outSIT, outstdSIT, outlnSIT, outUncSIT, outFrb, outstdFrb, outlnFrb, outUncFrb, outVAR1, outVAR2, outQFT, outQFS, outQFG

         

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





            

