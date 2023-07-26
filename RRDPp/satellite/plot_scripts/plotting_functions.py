# -*- coding: utf-8 -*-

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@dtu.dk']
__version__ = '0'
__date__ = '2023-07-25'

# -- Built-in modules -- #
# -- Third-part modules -- #
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
# -- Proprietary modules -- #

#%% Creates statistics file
def stats(self, ofile):
    """
    Prints statistics regarding the fit of reference and satelitte data
    The statistics are printed to the file sat + satistics.dat

    Returns
    -------
    None.

    """
    if self.var=='SID':
        obsVar=self.obsSID
        satVar=self.satSID
    elif self.var=='SIT':
        obsVar=self.obsSIT
        satVar=self.satSIT
    elif self.var=='SD':
        obsVar=self.obsSD
        satVar=self.satSD
    elif self.var=='SIF':
        obsVar=self.obsSIF
        satVar=self.satSIF
    
    # loop through one campaign at a time
    for oV, sV, name in zip(obsVar, satVar, self.name):
        
        # use only data where both satelitte and reference data has non nan values
        index = [True if ~np.isnan(ov) and ~np.isnan(sv) else False for ov,sv in zip(oV, sV)]
        oV = oV[index]
        sV = sV[index]
        
        satAvg = np.nanmean(sV)
        satStd = np.nanstd(sV)
        obsAvg = np.nanmean(oV)
        obsStd = np.nanstd(oV)
        bias   = obsAvg - satAvg
        
        # best fit
        model  = np.polyfit(oV, sV, 1)
        # one to one line
        model2 = np.polyfit(oV, oV, 1)
    
        predict = np.poly1d(model)
        predict2 = np.poly1d(model2)
        
        # correlation coefficient
        r, _ = pearsonr(sV, predict(oV))
        r2, _ = pearsonr(sV, predict2(oV))
        # r2 score
        rr    = r2_score(sV, predict(oV))
        rr2   = r2_score(sV, predict2(oV))
        # root mean square error 
        rmse  = np.sqrt(mean_squared_error(sV, predict(oV)))
        rmse2 = np.sqrt(mean_squared_error(sV, predict2(oV)))
    
        output=open(ofile,'a')
        print('name: ',self.name)
        print('R: ',r)
        print('R-sqaured: ',rr)
        print('RMSE: ',rmse)
        print('R: ',r2)
        print('R-sqaured: ',rr2)
        print('RMSE: ',rmse2)
        print('#points: ',len(sV))
        print('Linear fit: y=',model[0],'x + ',model[1])
        print('%5s & %5s & %5.2f &  %5.2f & %5.2f &  %5.2f & %5.2f & %5.2f & %5.2f & %6.2f & %6.2f & %8i & y=%6.2f*x + %6.2f\\' % (name, self.var, obsAvg, satAvg, obsStd, satStd, bias, rr,rmse,rr2,rmse2,len(sV),model[0],model[1]),file=output)
        # print(OBSID + '_' + sat + '_' + var + 'ANT.png',file=output)
#%% Creates histograms
def Histograms(self):
    """
    Makes a combined histogram for each variable
    showing a comparison of the reference data and co-located
    satellite data

    Returns
    -------
    None.

    """
    
    # define variable
    if self.var=='SID':
        obsVar=self.obsSID
        satVar=self.satSID
        bins=np.linspace(-0.1, 3, 60)
    elif self.var=='SIT':
        obsVar=self.obsSIT
        satVar=self.satSIT
        bins = np.linspace(0, 8, 60)
    elif self.var=='SD':
        obsVar=self.obsSD
        satVar=self.satSD
        bins = np.linspace(0, 1, 60)
    elif self.var=='SIF':
        obsVar=self.obsSIF
        satVar=self.satSIF
        bins = np.linspace(-0.5, 1, 60)

    # Histograms
    fig, ax1 = plt.subplots(len(obsVar), figsize=(8,8), sharex=True, layout='constrained')
    
    for i, name in zip(range(len(obsVar)), self.name):
        ax1[i].set_ylabel(name, fontsize=14)
        ax1[i].set_yticks([])
        ax1[i].tick_params(axis="x", labelsize=14)
        ax2 = ax1[i].twinx()  # instantiate a second axes that shares the same x-axis
        if name=='ASPeCt' or 'SH' in name or name=='AWI-ULS':
            ax2.hist(obsVar[i], bins=bins,  rwidth=0.85, linewidth=0.8, 
             linestyle='-', edgecolor='k', alpha=0.7, stacked=True)
            ax2.hist(satVar[i], bins=bins,  rwidth=0.85, linewidth=0.8, 
             linestyle='-', edgecolor='k', alpha=0.7)
        else:
            ax2.hist(obsVar[i], bins=bins,  rwidth=0.85, alpha=0.7 ) #, linewidth=0.5, edgecolor='k', stacked=True)
            ax2.hist(satVar[i], bins=bins,  rwidth=0.85, alpha=0.7 ) #, linewidth=0.5, edgecolor='k')
        
        ax2.tick_params(axis="y", labelsize=14)
        ax2.axvline(x=np.nanmean(obsVar[i]), c='r', linewidth=3, label='obs mean')
        ax2.axvline(x=np.nanmean(satVar[i]), c='k', linewidth=3, label='sat mean')
        
    ax1[i].set_xlabel('[m]', fontsize=14)
    # Make custom legend
    legend_elements = [mpatches.Patch(color='tab:blue', alpha=0.7, label='obs ' + self.var + ' [m]'),
                      mpatches.Patch(color='tab:orange', alpha=0.7,label=self.sat+ ' '+ self.var + ' [m]'),
                      mpatches.Patch(facecolor='tab:blue', alpha=0.7, edgecolor='k', linewidth=1.2, 
                       linestyle='-',  label='obs SH. ' + self.var + ' [m]'),
                      mpatches.Patch(facecolor='tab:orange', alpha=0.7, edgecolor='k', linewidth=1.2, 
                       linestyle='-',  label='obs SH. ' + self.var + ' [m]'),
                     mlines.Line2D([], [],color='r', linewidth=3,label='obs mean [m]'),
                   mlines.Line2D([], [],color='k', linewidth=3,label='sat mean [m]')]
    if self.sat=='ENV':
        fig.legend(handles=legend_elements,prop={"size":14}, bbox_to_anchor=(1.35, 0.63), bbox_transform=fig.transFigure)
    else:
        fig.legend(handles=legend_elements,prop={"size":14}, bbox_to_anchor=(1.3, 0.63), bbox_transform=fig.transFigure)


    plt.savefig(self.var + '_hist_'+self.sat+'.png', bbox_inches='tight')
    plt.show()
    
#%% Make scatterplots of all data available for each variable
def scatter(self):
    """
    Makes a combined scatterplots for each variable
    showing a comparison of the reference data and co-located
    satellite data

    Returns
    -------
    None.

    """
    if self.var=='SID':
        obsVar=self.obsSID
        satVar=self.satSID
    elif self.var=='SIT':
        obsVar=self.obsSIT
        satVar=self.satSIT
    elif self.var=='SD':
        obsVar=self.obsSD
        satVar=self.satSD
    elif self.var=='SIF':
        obsVar=self.obsSIF
        satVar=self.satSIF

    plt.figure(figsize=(8,8))
    for i, name in zip(range(len(obsVar)), self.name):
        if self.name[i]=='OIB':
            obsVar[i] = obsVar[i][::10]
            satVar[i] = satVar[i][::10]
        plt.scatter(obsVar[i], satVar[i], color=self.c[i], label=self.name[i], alpha=0.4 + 1/(i+2))
    
    if self.var!='SD':
        plt.plot([-0.10, np.nanmax(np.concatenate(satVar))], [0, np.nanmax(np.concatenate(satVar))], c='k')
        plt.xlim([-0.10,np.nanmax(np.concatenate(satVar))])
        plt.ylim([-0.10,np.nanmax(np.concatenate(satVar))])
    else:
        plt.plot([0, 1], [0, 1], c='k')
        plt.xlim([-0.10,1])
        plt.ylim([-0.10,1])        
    plt.legend()
    plt.xlabel('Observation ' + self.var + ' [m]')
    plt.ylabel('Satelitte ' + self.var + ' [m]')
    plt.show()
    
#%% Data class to be created for each variable
class Data:
    def __init__(self,sat):
      self.name = []
      self.sat = sat
      self.var = ''
      self.c = []
      self.lat = []
      self.lon = []
      self.obsSIT = []
      self.obsSD = []
      self.obsSID = []
      self.obsSIF = []
      self.satSIT = []
      self.satSD = []
      self.satSID = []
      self.satSIF = []
    
    def sort(self):
        self.name = [name if name!='ULS' else 'AWI-ULS' for name in self.name]
        ## find Antarctic campaigns and put these last
        order = [i for i in range(len(self.name))]
        
        if 'ASPeCt' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name)))if name=='ASPeCt'][0]
            order.remove(ind)
            order.append(ind)
        if 'SB-AWI-SH' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name)))if name=='SB-AWI-SH'][0]
            order.remove(int(ind))
            order.append(int(ind))
        if 'AWI-ULS' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name)))if name=='AWI-ULS'][0]
            order.remove(int(ind))
            order.append(int(ind))
        
        self.name = np.array(self.name)[order]
        self.lat = np.array(self.lat)[order]
        self.lon = np.array(self.lon)[order]
        self.c = np.array(self.c)[order]
        
        if self.var=='SIT':
            self.obsSIT = np.array(self.obsSIT)[order]
            self.satSIT = np.array(self.satSIT)[order]   
        if self.var=='SD':
            self.obsSD = np.array(self.obsSD)[order]
            self.satSD = np.array(self.satSD)[order]
        if self.var=='SID':
            self.obsSID = np.array(self.obsSID)[order]
            self.satSID = np.array(self.satSID)[order]        
        if self.var=='SIF':
            self.obsSIF = np.array(self.obsSIF)[order]
            self.satSIF = np.array(self.satSIF)[order]   