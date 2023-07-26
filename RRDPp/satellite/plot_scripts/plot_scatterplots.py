#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 11:28:06 2021

@author: hsko
"""

import numpy as np
import datetime as dt
import matplotlib.pyplot as plt

from sklearn.metrics import r2_score, mean_squared_error
import os,sys    

#--------- DEFINE THESE -------------------------
OBSID   = 'NP'
sat     = 'CS-2'
version = 'v3p0-rc2'
var     = 'SIT'
HS      = 'NH'

obsfile = "C:/Users/Ida Olsen/Documents/work/RRDPp/satelitte/Final_files/ESACCIplus-SEAICE-RRDP2+-"+var+"-"+OBSID+"-"+sat+"-CCIp-"+version+".dat"  #ESACCIplus-SEAICE-RRDP2+-SD-"+OBSID+"-" + sat + "-test-CCIp-v3p0-rc2.dat"
# obsfile = "C:/Users/Ida Olsen/Documents/work/RRDPp/satelitte/Final_files/ESACCIplus-SEAICE-RRDP2+-SIT-OIB-IDCS4-ANTARCTIC-ENV--CCIp-v3p0-rc2.dat"

# var = 'SD'
path = 'C:/Users/Ida Olsen/Documents/work/RRDPp/satelitte/Final_files'
obsfile = os.path.join(path, obsfile)

ofile = 'statistics.dat'
if sat=='ENV':
    col='r'
if sat=='CS-2':
    col='g'

# plotName=os.path.basename(obsfile)[:-4]

# Reads input data from observation file (RRDP), ASCII format
print('Reads obs-data',dt.datetime.now())
file = open(obsfile, 'r')

if var == 'SID':
    try:
        ObsNames=['date','lat','lon','obsSID','satSID','obsstd','obsln','satstd','satln']
        ObsData = np.genfromtxt(obsfile,dtype=None,names=ObsNames)
        
    except:
        ObsNames=['date','lat','lon','obsSID','satSID','obsstd','obsln','obsUnc','satstd','satln', 'satUnc']
        ObsData = np.genfromtxt(obsfile,dtype=None,names=ObsNames)
    
    obsDate = ObsData['date']
    obslat = ObsData['lat']
    obslon = ObsData['lon']
    obsVar = ObsData['obsSID']
    
    satVar = ObsData['satSID']
    try:
        obsUnc = ObsData['obsUnc']
        satUnc = ObsData['satUnc']
        obsStd = ObsData['obsstd']
        satStd = ObsData['satstd']
    except:
        print('no uncertainties')
else:
        
    ObsNames=['date','lat','lon','obsSD','obsSIT','obsFRB','satSD','satSIT',
          'satFRB', "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
          "obsFRB_std", "obsFRB_ln","obsFRB_unc", "satSD_std","satSD_ln","satSD_unc","satSIT_std",
          "satSIT_ln","satSIT_unc", "satFRB_std", "satFRB_ln", "satFRB_unc"]

    # ObsNames=['date','lat','lon','obsSD','obsSIT','obsSIF','satSIT',
    #           'satSIF', "obsSD_std","obsSD_ln", "obsSD_unc","obsSIT_std","obsSIT_ln", "obsSIT_unc","obsFRB_std",
    #           "obsFRB_ln", "obsFRB_unc", "satSIT_std","satSIT_ln", "satSIT_unc","satFRB_std",
    #           "satFRB_ln", "satFRB_unc"]
    ObsData = np.genfromtxt(obsfile,dtype=None,names=ObsNames)
    obsDate = ObsData['date']
    obslat = ObsData['lat']
    obslon = ObsData['lon']
    if var == 'SIT':
        obsVar = ObsData['obsSIT']
        satVar = ObsData['satSIT']
        obsUnc = ObsData['obsSIT_unc']
        satUnc = ObsData['satSIT_unc']
        obsStd = ObsData['obsSIT_std']
        satStd = ObsData['satSIT_std']
        
        if OBSID == 'AEM-AWI':
            satVar = satVar + ObsData['satSD']

    elif var == 'SD':
        obsVar = ObsData['obsSD']
        satVar = ObsData['satSD']   
        try:
            obsUnc = ObsData['obsSD_unc']
            satUnc = ObsData['satSD_unc']
            obsStd = ObsData['obsSD_std']
            satStd = ObsData['satSD_std']
        except:
            pass
    elif var == 'FRB':
        # obsVar = ObsData['obsFRB']
        # OBS IF OIB remember FRB = frb -SD
        if HS=='SH' and 'OIB' in OBSID: # total freeboard
            obsVar = ObsData['obsFRB']
            satVar = ObsData['satFRB'] + ObsData['satSD']
        else: # sea ice freeboard
            obsVar = ObsData['obsFRB'] - ObsData['obsSD']
            satVar = ObsData['satFRB']
            
            
        try:
            obsStd = ObsData['obsFRB_std']
            satStd = ObsData['satFRB_std']
            obsUnc = ObsData['obsFRB_unc']
            satUnc = ObsData['satFRB_unc']
        except:
            pass
        


unicode_obsdate = obsDate.view(np.chararray).decode('utf-8')

obstime = np.array([dt.datetime.strptime(s, "%Y-%m-%dT%H:%M:%S") 
                 for s in unicode_obsdate])

index = [True if ~np.isnan(ov) and ~np.isnan(sv) else False for ov,sv in zip(obsVar, satVar)]
# obsVarnan = obsVar[~np.isnan(obsVar)]
# satVarnan = satVar[~np.isnan(obsVar)] 

obsVarnan = obsVar[index]
satVarnan = satVar[index] 

#obsVarnan_unc = obsUnc[~np.isnan(obsVar)] + obsStd[~np.isnan(obsVar)]
#satVarnan_unc = satUnc[~np.isnan(obsVar)]+ satStd[~np.isnan(obsVar)]

obstimenan =  obstime[~np.isnan(satVar)]

# satVar2nan = satVarnan[~np.isnan(satVarnan)]
# obsVar2nan = obsVarnan[~np.isnan(satVarnan)]
# obstime2nan = obstimenan[~np.isnan(satVarnan)]

#obsVar2nan_unc = obsVarnan_unc[~np.isnan(satVarnan)]
#satVar2nan_unc = satVarnan_unc[~np.isnan(satVarnan)]

year = np.empty(np.shape(obstimenan))
for ii in range(len(obstimenan)):
    year[ii] = obstimenan[ii].year


model = np.polyfit(obsVarnan, satVarnan, 1)
model2 = np.polyfit(obsVarnan, obsVarnan, 1)


predict = np.poly1d(model)
predict2 = np.poly1d(model2)

rr=r2_score(satVarnan, predict(obsVarnan))
rr2=r2_score(satVarnan, predict2(obsVarnan))

rmse=np.sqrt(mean_squared_error(satVarnan, predict(obsVarnan)))
rmse2=np.sqrt(mean_squared_error(satVarnan, predict2(obsVarnan)))

x_lin_reg = range(-1, 8)
y_lin_reg = predict(x_lin_reg)

x_one2one = range(-1, 8)
y_one2one = predict2(x_lin_reg)

output=open(ofile,'a')
print('R-sqaured: ',rr)
print('RMSE: ',rmse)
print('R-sqaured: ',rr2)
print('RMSE: ',rmse2)
print('#points: ',len(satVarnan))
print('Linear fit: y=',model[0],'x + ',model[1])
print('%5s & %5s &  %5s & %5s & %5.3f & %5.3f & %6.3f & %6.3f & %8.2i & y=%6.3f*x + %6.3f' % (OBSID,sat,version,var,rr,rmse,rr2,rmse2,len(satVarnan),model[0],model[1]),file=output)
print(OBSID + '_' + sat + '_' + var + 'ANT.png',file=output)

if var=='SD':
    xmin,xmax=[-0.1,1.4]
    ymin,ymax=[-0.1,1.4]
if var=='SIT':
    xmin,xmax=[-0.5,7]
    ymin,ymax=[-0.5,7]
if var=='SID':
    xmin,xmax=[-0.5,6]
    ymin,ymax=[-0.5,6]
if var=='FRB':
    xmin,xmax=[-0.5,1.5]
    ymin,ymax=[-0.5,1.5]
    

fig=plt.figure(figsize=(7,7))
# plt.errorbar(obsVar2nan, satVar2nan, yerr=satVar2nan_unc, xerr=obsVar2nan_unc, fmt='.', color=col, alpha=0.7)
plt.scatter(obsVarnan,satVarnan, s=10, color=col,edgecolor=col)

plt.plot(x_lin_reg, y_lin_reg, c = 'gray')
plt.plot(x_one2one, y_one2one, c = 'black')

plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.tick_params(axis='both', labelsize=14)
plt.xlabel('%3s %3s (m)' % (OBSID.replace('-IDCS4', ''), var),fontsize=14)
plt.ylabel('%3s %3s (m)' % (sat, var),fontsize=14)
# plt.annotate("\n rmse: "+ str(round(rmse,3))+ '\n r² :'+ str(round(rr,3))+"\n rmse of 1/1 line: "+ str(round(rmse2,3))+ '\n r² of 1/1 line:'+ str(round(rr2,3)),(0,0))
plt.title('%3s vs %3s %5s' % (OBSID.replace('-IDCS4', ''), sat, version),fontsize=16)
plt.savefig(OBSID + '_' + sat + '_' + var + '.png', bbox_inches='tight')    

output.close()
file.close()