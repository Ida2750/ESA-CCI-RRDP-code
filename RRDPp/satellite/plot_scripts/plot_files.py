# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 15:31:46 2022

@author: Ida Olsen

"""

# Built-in packages
import sys
import os
import datetime as dt

# Append custom plot scripts directory to system path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'plot_scripts'))

# Third-party packages
import numpy as np
import io
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from sklearn.metrics import r2_score, mean_squared_error
from scipy.stats import pearsonr
import scipy.odr as spodr
from sklearn.linear_model import LinearRegression
import warnings
from netCDF4 import Dataset
import xarray as xr
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

# Set matplotlib parameters for better visualization
mpl.rcParams['font.size'] = 48.0
mpl.rcParams['lines.linewidth'] = 5

# Ignore warnings
warnings.filterwarnings("ignore")

# Homemade modules
import polar_plots as pp

# %% Weighted linear fit
def odr_fit(x_o, y_o, sx, sy):

    # Check if the standard deviations sx have any non-NaN values
    if len(sx[~np.isnan(sx)]) == 0:
        # If sx has no valid entries, set both sx and sy to ones (default weight)
        sx = np.ones(len(x_o))
        sy = np.ones(len(y_o))
        print('entered')  # Indicate that the default values are being used
    else:
        # Replace NaN values in sx and sy with their respective means
        sx[np.isnan(sx)] = np.nanmean(sx)
        sy[np.isnan(sy)] = np.nanmean(sy)

    # Add a tiny offset to the weights if they are zero to avoid division by zero issues
    sx[sx == 0] = 1e-3
    sy[sy == 0] = 1e-3

    # Define a linear function for the fitting process
    def linear_func(coefs, x):
        w, b = coefs  # Unpack coefficients into slope (w) and intercept (b)
        return w * x + b  # Linear equation y = wx + b

    # Create a model for fitting (linear).
    linear_model = spodr.Model(linear_func)

    # Create a RealData object using the provided data.
    # This object incorporates weighting based on standard deviations and uncertainties.
    data = spodr.RealData(x_o, y_o, sx=sx, sy=sy)

    # Set up Orthogonal Distance Regression (ODR) with the model and data.
    # beta0 provides initial guesses for the model parameters.
    odr = spodr.ODR(data, linear_model, beta0=[1., 1.])  # Initial guess: w=1, b=1

    # Run the regression analysis.
    out = odr.run()

    # Extract the slope (w) and intercept (b) from the output
    w, b = out.beta

    # Make predictions based on the fitted model
    modelPredictions_odr = linear_func(out.beta, x_o)

    # Calculate the absolute errors between the model predictions and actual values
    absError_odr = modelPredictions_odr - y_o
    # Uncomment the following line to print the mean squared error of the absolute errors
    # print(np.mean(np.square(absError_odr)))
    
    # Calculate squared errors
    SE_odr = (absError_odr) ** 2  # Squared errors
    MSE_odr = np.nanmean(SE_odr)   # Mean squared errors, ignoring NaNs

    # Calculate Root Mean Squared Error (RMSE)
    RMSE_odr = np.sqrt(MSE_odr)

    # Calculate the R-squared value for the fit
    Rsquared_odr = 1.0 - (np.var(absError_odr) / np.var(y_o))

    # Return the slope and intercept of the linear fit
    return w, b

#%% Histograms per month/year
    
def Histograms_time(d_SIT, d_SD, d_SIF, d_SID):
    import pandas as pd
    import numpy as np
    import datetime as dt
    import matplotlib.pyplot as plt

    # Initialize a dictionary to hold monthly sums for each dataset
    d = {'SIT': {}, 'SD': {}, 'SIF': {}, 'SID': {}}

    # Iterate over each dataset
    for v in [d_SIT, d_SD, d_SIF, d_SID]:
        # Convert the date strings into datetime objects
        dates = [dt.datetime.strptime(dd.decode('utf8'), '%Y-%m-%dT%H:%M:%S') for dd in np.concatenate((v.date))]
        
        # Extract the year and month from the datetime objects
        years = [d.year for d in dates]  # List of years
        months = [d.month for d in dates]  # List of months
        
        # Initialize a dictionary to store the count of observations for each month
        sum_months = {month: 0 for month in range(1, 13)}  # Months 1-12 initialized to 0
        
        # Count the number of observations for each month
        for m in sum_months:
            sum_months[m] = (sum(np.array(months) == m))  # Count occurrences of each month
        
        # Store the monthly sums in the main dictionary under the appropriate key
        d[v.var] = sum_months  # Store the monthly count in the dictionary
    
    # Configure plotting parameters
    plt.rcParams.update({'font.size': 14})  # Update font size for plots
    
    # Create a DataFrame from the dictionary and plot it
    pd.DataFrame(d).plot(kind='bar',
                         figsize=(10, 7),  # Set figure size
                         fontsize=14,  # Set font size for plot labels
                         title='Monthly distribution of observations',  # Title of the plot
                         layout='constrained',  # Adjust layout for better spacing
                         xlabel='Months',  # X-axis label
                         ylabel='Number of observations')  # Y-axis label
    
    #plt.show()  # Display the plot
    
#%% Scatterplots
def scatter(self, save=False):
    import matplotlib.pyplot as plt

    # Select the observed and satellite variables based on the 'var' attribute
    if self.var == 'SID':
        obsVar = self.obsSID
        satVar = self.satSID
        obsErr = [s + y for s, y, n in zip(self.obsSIDunc, self.obsSIDstd, self.name)]
        satErr = [s + y for s, y, n in zip(self.satSIDunc, self.satSIDstd, self.name)]
    elif self.var == 'SIT':
        obsVar = self.obsSIT
        satVar = self.satSIT
        obsErr = [s + y for s, y in zip(self.obsSITunc, self.obsSITstd)]
        satErr = [s + y for s, y in zip(self.satSITunc, self.satSITstd)]
    elif self.var == 'SD':
        obsVar = self.obsSD
        satVar = self.satSD
        obsErr = [s + y for s, y in zip(self.obsSDunc, self.obsSDstd)]
        satErr = [s + y for s, y in zip(self.satSDunc, self.satSDstd)]
    elif self.var == 'FRB':
        obsVar = self.obsSIF
        satVar = self.satSIF
        obsErr = [s + y for s, y in zip(self.obsSIFunc, self.obsSIFstd)]
        satErr = [s + y for s, y in zip(self.satSIFunc, self.satSIFstd)]

    obsvar, satvar, obserr, saterr = [], [], [], []
    fig, ax = plt.subplots(figsize=(35, 35))
    ax.grid(alpha=0.5, linewidth=2)

    print(f'variable : {self.var}')
    print(f'***********************')

    for i, name in zip(range(len(obsVar)), self.name):
        if self.name[i] != 'ASPeCt' and 'SH' not in self.name[i] and 'AWI-ULS' not in self.name[i]:
            ax.scatter(obsVar[i], satVar[i], color=self.c[i], s=40, alpha=0.4)
            corr = np.round(pearsonr(obsVar[i][~np.isnan(satVar[i])], satVar[i][~np.isnan(satVar[i])]), 2)
            print(self.name[i])
            w, b = odr_fit(obsVar[i][~np.isnan(satVar[i])], satVar[i][~np.isnan(satVar[i])],
                           sx=obsErr[i][~np.isnan(satVar[i])], sy=satErr[i][~np.isnan(satVar[i])])
            print('RMSE', np.round(np.sqrt(mean_squared_error(satVar[i][~np.isnan(satVar[i])],
                                                               (obsVar[i][~np.isnan(satVar[i])]*w + b))), 2))
            ax.plot(range(-1, 8), [v*w + b for v in range(-1, 8)], color=self.c[i], linewidth=2,
                    label=f'{name}, R:{corr[0]}')
            print('R2', np.round(r2_score(satVar[i][~np.isnan(satVar[i])],
                                          (obsVar[i][~np.isnan(satVar[i])]*w + b)), 2))
            print(f' wb: {w} {b}')
            print(f'len: {len(satVar[i][~np.isnan(satVar[i])])}')
            obsvar.append(obsVar[i][~np.isnan(satVar[i])])
            satvar.append(satVar[i][~np.isnan(satVar[i])])
            obserr.append(obsErr[i][~np.isnan(satVar[i])])
            saterr.append(satErr[i][~np.isnan(satVar[i])])

    x = np.concatenate(obsvar)
    y = np.concatenate(satvar)
    sx = np.concatenate(obserr)
    sy = np.concatenate(saterr)
    corr = np.round(pearsonr(x[~np.isnan(y)], y[~np.isnan(y)]), 2)
    w, b = odr_fit(x, y, sx, sy)
    ax.plot(range(-1, 8), [v*w + b for v in range(-1, 8)], color='r', label=f'Combined, R:{corr[0]}')

    # Add label text based on satellite type (bottom right corner)
    if self.sat == 'CS2':
        ax.text(0.92, 0.05, '(b)', transform=ax.transAxes, fontsize=80, fontweight='bold', 
                va='bottom', ha='right')
    elif self.sat == 'ENV':
        ax.text(0.92, 0.05, '(e)', transform=ax.transAxes, fontsize=80, fontweight='bold', 
                va='bottom', ha='right')
        
    # Axis setup
    if self.var == 'SIT':
        ax.plot([-0.10, 8], [-0.10, 8], c='k', label='1 to 1 line')
        ax.set_xlim([-0.10, 7])
        ax.set_ylim([-0.10, 7])
    elif self.var == 'SID':
        ax.plot([-0.10, 5], [-0.10, 5], c='k', label='1 to 1 line')
        ax.set_xlim([-0.10, 4])
        ax.set_ylim([-0.10, 4])
    elif self.var == 'SD':
        ax.plot([-0.1, 1], [-0.1, 1], c='k')
        ax.set_xlim([-0.1, 1])
        ax.set_ylim([-0.1, 1])
    else:
        ax.plot([-0.1, 1], [-0.1, 1], c='k')
        ax.set_xlim([-0.10, 1])
        ax.set_ylim([-0.10, 1])

    leg = ax.legend(fontsize=50, loc='upper left')
    for line in leg.get_lines():
        line.set_linewidth(4.0)

    ax.set_xlabel('Reference ' + self.var + ' [m]', fontsize=50)
    ax.set_ylabel('Satellite ' + self.var + ' [m]', fontsize=50)

    if save:
        fig.savefig(f'scatter_{self.sat}_{self.var}.png', bbox_inches='tight')

    # Return the figure object for combining externally
    return fig, ax


#%% Get statistics
def stats(self, sat, ofile):
    # Determine the observed and satellite variables based on the selected variable type
    if self.var == 'SID':
        obsVar = self.obsSID
        satVar = self.satSID
        obsErr = [s + y for s, y in zip(self.obsSIDunc, self.obsSIDstd)]
        satErr = [s + y for s, y in zip(self.satSIDunc, self.satSIDstd)]
    elif self.var == 'SIT':
        obsVar = self.obsSIT
        satVar = self.satSIT
        obsErr = [s + y for s, y in zip(self.obsSITunc, self.obsSITstd)]
        satErr = [s + y for s, y in zip(self.satSITunc, self.satSITstd)]
    elif self.var == 'SD':
        obsVar = self.obsSD
        satVar = self.satSD
        obsErr = [s + y for s, y in zip(self.obsSDunc, self.obsSDstd)]
        satErr = [s + y for s, y in zip(self.satSDunc, self.satSDstd)]
    elif self.var == 'FRB':
        obsVar = self.obsSIF
        satVar = self.satSIF
        obsErr = [s + y for s, y in zip(self.obsSIFunc, self.obsSIFstd)]
        satErr = [s + y for s, y in zip(self.satSIFunc, self.satSIFstd)]
    
    # Loop through each dataset and compute statistics
    for oV, sV, oE, sE, name in zip(obsVar, satVar, obsErr, satErr, self.name):
        # Print the name of the dataset being analyzed
        print(self.name)
        
        # Create an index mask to filter out NaN values from observed and satellite data
        index = [True if ~np.isnan(ov) and ~np.isnan(sv) else False for ov, sv in zip(oV, sV)]
        
        # Use the index to filter the observed and satellite data along with their errors
        oV = oV[index]
        oE = oE[index]
        sV = sV[index]
        sE = sE[index]
        
        # Calculate statistical metrics
        satAvg = np.nanmean(sV)  # Mean of satellite values
        satStd = np.nanstd(sV)   # Standard deviation of satellite values
        obsAvg = np.nanmean(oV)  # Mean of observed values
        obsStd = np.nanstd(oV)   # Standard deviation of observed values
        bias = np.nanmean(oV - sV)  # Bias (mean difference) between observed and satellite
        
        # Perform Orthogonal Distance Regression (ODR) fitting
        w, b = odr_fit(oV, sV, sx=oE, sy=sE)  # w: slope, b: intercept
        
        # Generate model predictions using ODR
        modelPredictions_odr = w * oV + b
        absError_odr = modelPredictions_odr - sV  # Calculate absolute error
        SE_odr = (absError_odr) ** 2  # Squared errors
        MSE_odr = np.nanmean(SE_odr)  # Mean squared error
        RMSE_odr = np.sqrt(MSE_odr)  # Root Mean Squared Error
        
        # Calculate R-squared value for ODR
        Rsquared_odr = 1.0 - (np.var(absError_odr) / np.var(sV))

        # Perform regular linear regression
        model = np.polyfit(oV, sV, 1)  # Linear fit (slope and intercept)
        model2 = np.polyfit(oV, oV, 1)  # 1-to-1 line fit for comparison
        
        predict = np.poly1d(model)  # Create a polynomial object for predictions
        predict2 = np.poly1d(model2)  # Create a polynomial object for 1-to-1 line
        
        # Calculate correlation coefficients
        r, _ = pearsonr(sV, predict(oV))  # Pearson correlation for fitted model
        rr = r2_score(sV, predict(oV))  # R-squared for fitted model
        r2, _ = pearsonr(sV, predict2(oV))  # Pearson correlation for 1-to-1 fit
        

        rr2 = r2_score(sV, predict2(oV))  # R-squared for 1-to-1 fit
        rmse = np.sqrt(mean_squared_error(sV, predict(oV)))  # RMSE for fitted model

        rmse2 = np.sqrt(mean_squared_error(sV, predict2(oV)))  # RMSE for 1-to-1 fit
        
        # Prepare data for plotting the linear regression line
        x_lin_reg = range(-1, 8)
        y_lin_reg = predict(x_lin_reg)  # Predictions based on fitted model
        

        # Prepare data for plotting the 1-to-1 line
        x_one2one = range(-1, 8)
        y_one2one = predict2(x_lin_reg)  # Predictions for 1-to-1 line
        
        # Create output file path for statistics
        ofile = os.getcwd() + '/' + sat + '_statistics.dat'
        output = open(ofile, 'a')  # Open file to append results
        
        # Print statistics to console and file
        print('===========================')  # Separator for readability
        print('name: ', name)
        print(f'bias: {bias}')  # Print bias
        print('R-squared best fit: ', rr)  # Print R-squared for fitted model
        print('RMSE: ', rmse)  # Print RMSE for fitted model
        print('R-squared 1 to 1: ', rr2)  # Print R-squared for 1-to-1 line
        print('RMSE: ', rmse2)  # Print RMSE for 1-to-1 line
        print('R-squared ODR: ', Rsquared_odr)  # Print R-squared for ODR
        print('RMSE ODR: ', RMSE_odr)  # Print RMSE for ODR
        print('#points: ', len(sV))  # Print number of valid points
        print('Linear fit: y=', w, 'x + ', b)  # Print linear fit coefficients

        # Output formatted statistics to the file
        print('%5s & %5s & %5.2f &  %5.2f & %5.2f &  %5.2f & %5.2f &  %5.2f & %5.2f & %5.2f & %5.2f & %6.2f & %6.2f & %8i & y=%6.2f x + %6.2f\\' % 
              (name, self.var, obsAvg, satAvg, obsStd, satStd, bias, Rsquared_odr, RMSE_odr, rr, rmse, rr2, rmse2, len(sV), w, b), file=output)
        
        print('%5s & %5s & %5.2f &  %5.2f & %5.2f &  %5.2f & %5.2f &  %5.2f & %5.2f & %5.2f & %5.2f & %6.2f & %6.2f & %8i & y=%6.2f x + %6.2f\\' % 
              (name, self.var, obsAvg, satAvg, obsStd, satStd, bias, Rsquared_odr, RMSE_odr, rr, rmse, rr2, rmse2, len(sV), w, b))
        
        print(self.name)  # Print the name of the dataset again (could be redundant)

#%% Histogram distribution
def Histograms(self, fig=None, axes=None):
    # Define bin edges
    if self.var == 'SID':
        obsVar = self.obsSID
        satVar = self.satSID
        bins = np.arange(-0.1, 3, 0.10)
    elif self.var == 'SIT':
        obsVar = self.obsSIT
        satVar = self.satSIT
        bins = np.arange(0, 7, 0.10)
    elif self.var == 'SD':
        obsVar = self.obsSD
        satVar = self.satSD
        bins = np.arange(0, 1, 0.03)
    elif self.var == 'FRB':
        obsVar = self.obsSIF
        satVar = self.satSIF
        bins = np.arange(-0.5, 1, 0.03)

    # Create or use existing figure and axes
    if fig is None or axes is None:
        fig, ax1 = plt.subplots(len(obsVar), figsize=(35, 35), sharex=True)
    else:
        ax1 = axes

    for i, name in zip(range(len(obsVar)), self.name):
        threshold = 10  # Adjust this number as needed

        if len(name) > threshold and ':' in name:
            # Split roughly in the middle at the space before the unit
            parts = name.split(':')
            name = parts[0] + ':\n' + ' '.join(parts[1:])
        
        ax1[i].set_ylabel(name, fontsize=50)
        ax1[i].set_yticks([])
        ax1[i].tick_params(axis="x", labelsize=50)
        ax1[i].grid()

        ax2 = ax1[i].twinx()

        if name == 'ASPeCt' or 'SH' in name or name == 'AWI-ULS':
            ax2.hist(obsVar[i], bins=bins, rwidth=0.85, linewidth=4,
                     linestyle='-', edgecolor='k', alpha=0.7, stacked=True)
            ax2.hist(satVar[i], bins=bins, rwidth=0.85, linewidth=4,
                     linestyle='-', edgecolor='k', alpha=0.7)
        else:
            ax2.hist(obsVar[i], bins=bins, rwidth=0.85, alpha=0.7)
            ax2.hist(satVar[i], bins=bins, rwidth=0.85, alpha=0.7)

        ax2.tick_params(axis="y", labelsize=50)
        ax2.axvline(x=np.nanmean(obsVar[i]), c='r', linewidth=8, label='obs mean')
        ax2.axvline(x=np.nanmean(satVar[i]), c='k', linewidth=8, label='sat mean')
        ax2.grid()

    # Label for x-axis
    ax1[-1].set_xlabel(f'{self.var} [m]', fontsize=50)

    # Adjust satellite label for FRB
    if self.var == 'FRB':
        sat_label = 'ENV & CS2'
    else:
        sat_label = self.sat

    # Legend definition (not added here to avoid duplication in combined plots)
    legend_elements = [
        mpatches.Patch(color='tab:blue', alpha=0.7, label='obs ' + self.var + ' [m]'),
        mpatches.Patch(color='tab:orange', alpha=0.7, label=sat_label + ' ' + self.var + ' [m]'),
        mpatches.Patch(facecolor='tab:blue', alpha=0.7, edgecolor='k', linewidth=1.2,
                       linestyle='-', label='obs SH. ' + self.var + ' [m]'),
        mpatches.Patch(facecolor='tab:orange', alpha=0.7, edgecolor='k', linewidth=1.2,
                       linestyle='-', label=sat_label + ' SH. ' + self.var + ' [m]'),
        mlines.Line2D([], [], color='r', linewidth=3, label='obs mean [m]'),
        mlines.Line2D([], [], color='k', linewidth=3, label='sat mean [m]')
    ]

    # Super y-label
    fig.supylabel('count of observations in each bin', fontsize=65, position=(0.97, 0.5))

    # Add bottom-right subplot label
    label = None
    if self.sat == 'CS2':
        label = '(c)'
    elif self.sat == 'ENV':
        label = '(f)'

    if label:
        ax1[-1].text(0.92, 0.05, label, transform=ax1[-1].transAxes,
                     fontsize=80, fontweight='bold', ha='right', va='bottom')

    # Save if standalone
    if fig is None or axes is None:
        plt.savefig(f'{self.var}_hist_{self.sat}.png', bbox_inches='tight')
        print(f'{self.var}_hist_{self.sat}.png')
        #plt.show()
        plt.close()

    return fig, ax1

#%% Data class
class Data:
    # This class is designed to read and organize satellite data, 
    # assigning colors, names, and other relevant information.
    def __init__(self, sat):
        # Initialize the attributes of the class
        self.name = []  # List to store the names of the datasets
        self.var = ''  # Variable type (e.g., 'SIT' for Sea Ice Thickness, 'SD' for Snow Depth, etc.)
        self.sat = sat  # Satellite name passed during initialization
        self.c = []  # List to store colors associated with datasets
        
        # Initialize attributes related to hemispheres, dates, and coordinates
        self.HS = []  # Hemisphere (e.g., North/South)
        self.date = []  # Dates of observations
        self.lat = []  # Latitude coordinates
        self.lon = []  # Longitude coordinates
        
        # Observational and satellite data for different variables
        self.obsSIT = []  # Observed Sea Ice Thickness
        self.obsSITstd = []  # Standard deviation of observed Sea Ice Thickness
        self.obsSITunc = []  # Uncertainty of observed Sea Ice Thickness
        self.satSITstd = []  # Standard deviation of satellite Sea Ice Thickness
        self.satSITunc = []  # Uncertainty of satellite Sea Ice Thickness
        
        self.obsSD = []  # Observed Snow Depth
        self.obsSID = []  # Observed Sea Ice Draft
        self.obsSIDunc = []  # Uncertainty of observed Sea Ice Draft
        self.obsSIDstd = []  # Standard deviation of observed Sea Ice Draft
        self.obsSIF = []  # Observed Freeboard
        
        self.satSIT = []  # Satellite Sea Ice Thickness
        self.satSD = []  # Satellite Snow Depth
        
        self.obsSDstd = []  # Standard deviation of observed Snow Depth
        self.obsSDunc = []  # Uncertainty of observed Snow Depth
        self.satSDstd = []  # Standard deviation of satellite Snow Depth
        self.satSDunc = []  # Uncertainty of satellite Snow Depth
        
        self.satSID = []  # Satellite Sea Ice Draft
        self.satSIDunc = []  # Uncertainty of satellite Sea Ice Draft
        self.satSIDstd = []  # Standard deviation of satellite Sea Ice Draft
        self.satSIF = []  # Satellite Freeboard
        
        self.obsSIFunc = []  # Uncertainty of observed Freeboard
        self.obsSIFstd = []  # Standard deviation of observed Freeboard
        self.satSIFstd = []  # Standard deviation of satellite Freeboard
        self.satSIFunc = []  # Uncertainty of satellite Freeboard
    
    def sort(self):
        # Update dataset names, ensuring consistency in naming
        self.name = [name if name != 'ULS' else 'AWI-ULS' for name in self.name]

        # Create an order list to sort Antarctic campaigns last
        order = [i for i in range(len(self.name))]
        
        # Move specific datasets to the end of the list
        if 'ASPeCt' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name))) if name == 'ASPeCt'][0]
            order.remove(ind)  # Remove index of ASPeCt
            order.append(ind)  # Append ASPeCt index to the end
        
        if 'SB-AWI-SH' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name))) if name == 'SB-AWI-SH'][0]
            order.remove(int(ind))  # Remove index of SB-AWI-SH
            order.append(int(ind))  # Append SB-AWI-SH index to the end
        
        if 'AWI-ULS' in self.name:
            ind = [i for name, i in zip(self.name, range(len(self.name))) if name == 'AWI-ULS'][0]
            order.remove(int(ind))  # Remove index of AWI-ULS
            order.append(int(ind))  # Append AWI-ULS index to the end
        
        # Use the order list to sort various attributes
        self.name = np.array(self.name, dtype=object)[order]
        self.lat = np.array(self.lat, dtype=object)[order]
        self.lon = np.array(self.lon, dtype=object)[order]
        self.c = np.array(self.c, dtype=object)[order]
        self.HS = np.array(self.HS, dtype=object)[order]
        
        # Sort observational and satellite data based on the variable type
        if self.var == 'SIT':
            self.obsSIT = np.array(self.obsSIT, dtype=object)[order]
            self.satSIT = np.array(self.satSIT, dtype=object)[order]
            self.obsSITunc = np.array(self.obsSITunc, dtype=object)[order]
            self.satSITunc = np.array(self.satSITunc, dtype=object)[order]   
            self.obsSITstd = np.array(self.obsSITstd, dtype=object)[order]
            self.satSITstd = np.array(self.satSITstd, dtype=object)[order]  
        
        if self.var == 'SD':
            self.obsSD = np.array(self.obsSD, dtype=object)[order]
            self.satSD = np.array(self.satSD, dtype=object)[order]
            self.obsSDunc = np.array(self.obsSDunc, dtype=object)[order]
            self.satSDunc = np.array(self.satSDunc, dtype=object)[order]   
            self.obsSDstd = np.array(self.obsSDstd, dtype=object)[order]
            self.satSDstd = np.array(self.satSDstd, dtype=object)[order] 
        
        if self.var == 'SID':
            self.obsSID = np.array(self.obsSID, dtype=object)[order]
            self.satSID = np.array(self.satSID, dtype=object)[order]  
            self.obsSIDunc = np.array(self.obsSIDunc, dtype=object)[order]
            self.satSIDunc = np.array(self.satSIDunc, dtype=object)[order]   
            self.obsSIDstd = np.array(self.obsSIDstd, dtype=object)[order]
            self.satSIDstd = np.array(self.satSIDstd, dtype=object)[order]        
        
        if self.var == 'FRB':
            self.obsSIF = np.array(self.obsSIF, dtype=object)[order]
            self.satSIF = np.array(self.satSIF, dtype=object)[order]
            self.obsSIFunc = np.array(self.obsSIFunc, dtype=object)[order]
            self.satSIFunc = np.array(self.satSIFunc, dtype=object)[order]   
            self.obsSIFstd = np.array(self.obsSIFstd, dtype=object)[order]
            self.satSIFstd = np.array(self.satSIFstd, dtype=object)[order] 


#%% Get campaign name
#%% Get name of campaign
def Get_name(ifile):
    name = ''.join(filter(str.isalnum, ifile))
    name = name.replace("ESACCIplusSEAICERRDP2SIT", "")
    name = name.replace("ESACCIplusSEAICERRDP2SID", "")
    name = name.replace("ESACCIplusSEAICERRDP2FRB", "")
    name = name.replace("ESACCIplusSEAICERRDP2SD", "")
    name = name.replace("CCIpv3p0rc2dat", "")
    name = name.replace("CCIpv3p0rc2nc", "")
    name = name.replace("IDCS4", "")
    name = name.replace(sat, "")
    name = name.replace("ANTARCTIC", "-SH") 
    name = name.replace("AWISHSH", "AWI:SH") 
    name = name.replace("sorted", "") 
    name = name.replace("ULS", "-ULS") 
    name = name.replace("ENV", "-ENV")
    name = name.replace("CSIM", "C: SIM")
    name = name.replace("ASPeCtSH", "ASPeCt")
    name = name.replace("Nansenlegacy", "Nansen_legacy")
    if 'SIMBA' not in name:
        name = name.replace("IMB", "IMB-CRREL")
    name = name.replace("NPI", "NPI-FS")
    name = name.replace("NH", "")
    if 'ULS' not in name:
        name = name.replace("AWI", "-AWI")
        name = name.replace("AWIFRB", "AWI-FRB")
    name = name.replace("AEMcomb", "AEM: (AWI + MOSAIC)")
    return name

#%% assign color to name
def Get_color(name):

    ## colorblind colors 
    colors = ['#006BA4', '#FF800E', '#ABABAB',
                      '#595959', '#5F9ED1', '#C85200',
                      '#898989', '#A2C8EC', '#FFBC79', '#CFCFCF']
    if 'OIB' in name:
        c=colors[2] # 'grey'
    elif 'ASSIST' in name or 'ASPeCt' in name:
        #if ((sat=='CS2') & (HS=='SH')):
        #    c=colors[5] # 'orange'
        #else:
        c=colors[1] # 'orange'
    elif 'IMB' in name and 'SIMBA' not in name:
        name = name.replace('B', 'B-CRREL')
        c="#832db6" # purple
    elif 'AEM' in name:
        name = name.replace('M', 'M-AWI')
        c=colors[0] # blue
    elif 'SB' in name or 'SB-SH' in name:
        name = name.replace('B', 'B-AWI')
        c=colors[4] # light blue
    ## SID 
    elif 'SCICEX' in name:
        c=colors[5] # dark orange
    elif 'BGEP' in name:
        c=colors[6] # grey
    elif 'NPI' in name:
        name = name + '-FS'
        c='blue'
    elif 'TRANSDRIFT' in name:
        c='k' # orange
    elif 'ULS' in name:
        c='k' # black
    else:
        c='green'
        
    return c

def Append_data(d_SID, d_SIT, d_SD, d_SIF, ifile, name, c, directory, filt=True):
        ## append data
        if 'SID' in ifile:
            #print(name)
            #if 'AWI-ULS' not in name: 
                # define observation names
                ObsNames=['date','lat','lon','obsSID', 'satSID', 'obsSIDstd','obsSIDln','obsSIDunc',
                          'QFT', 'QFS', 'QFG', 'satSIDstd','satSIDln', 'satSIDunc', 'index']
            
                # Read data
                ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None,names=True) #,skip_header=1,names=ObsNames)
            
                CRS_ENV = np.arange(0,len(ObsData))
                
                QFT = ObsData['QFT'][CRS_ENV]
    
                if filt:
                    index = QFT<3
                else:
                    index = (QFT<4) | (np.isnan(QFT))
                    if name=='SCICEX':
                        years = np.array([dt.datetime.strptime(dd.decode('utf8'), '%Y-%m-%dT%H:%M:%S').year for dd in ObsData['date'][CRS_ENV]])
                        index = years != 2014
                    if 'AWI-ULS' in name:
                        SID = ObsData['obsSID'][CRS_ENV][index]
                        index = SID!=0
                        
                    
                obsDate = ObsData['date'][CRS_ENV][index]
                obslat  = ObsData['lat'][CRS_ENV][index]
                obslon  = ObsData['lon'][CRS_ENV][index]
                obsSID  = ObsData['obsSID'][CRS_ENV][index]
                obsSIDunc = ObsData['obsSID_unc'][CRS_ENV][index]
                obsSIDstd = ObsData['obsSID_std'][CRS_ENV][index]
                obsSIDln = ObsData['obsSID_ln'][CRS_ENV][index]
                
                satSID  = ObsData['satSID'][CRS_ENV][index]
                satSIDunc = ObsData['satSID_unc'][CRS_ENV][index]
                satSIDstd = ObsData['satSID_std'][CRS_ENV][index]
               
                m=np.logical_and(~np.isnan(obsSID), obsSIDln>0)
    
                # assign data to Data object
                if any(m):
                    #print(d_SID.var)
                    d_SID.lat.append(obslat[m])
                    d_SID.lon.append(obslon[m])
                    d_SID.obsSID.append(obsSID[m])
                    d_SID.obsSIDunc.append(obsSIDunc[m])
                    d_SID.obsSIDstd.append(obsSIDstd[m])
                    d_SID.satSID.append(satSID[m])
                    d_SID.satSIDunc.append(satSIDunc[m])
                    d_SID.satSIDstd.append(satSIDstd[m])
                    d_SID.name.append(name)
                    d_SID.date.append(obsDate[m])
                    d_SID.c.append(c)
                    d_SID.HS.append(HS)

        
        
        elif 'SIT' in ifile or 'SD' in ifile or 'FRB' in ifile:
            
            # define observation names
            ObsNames=['date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT',
                      'satSIF', "obsSD_std","obsSD_ln", "obsSD_unc","obsSIT_std","obsSIT_ln", "obsSIT_unc","obsFRB_std",
                      "obsFRB_ln", "obsFRB_unc", 'QFT', 'QFS', "satSD_std","satSD_ln", "satSD_unc","satSIT_std","satSIT_ln", "satSIT_unc","satFRB_std",
                      "satFRB_ln", "satFRB_unc", 'index']
            # Read data
            try:
                ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None, names=True) #skip_header=1, names=ObsNames)
            except:
                ObsNames2 = [
                    'date','lat','lon','obsSD','obsSIT','obsSIF','satSD','satSIT','satSIF',
                    "obsSD_std","obsSD_ln","obsSD_unc","obsSIT_std","obsSIT_ln","obsSIT_unc",
                    "obsFRB_std","obsFRB_ln","obsFRB_unc", "QFT", "QFS","satSD_std","satSD_ln","satSD_unc",
                    "satSIT_std","satSIT_ln","satSIT_unc","satFRB_std","satFRB_ln","satFRB_unc"
                ]
                ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None, names=True) #skip_header=1, names=ObsNames2)
            
            #print(ObsData.dtype.names)
            # identify copy rows
            #ObsData = np.unique(ObsData, axis=0)
 
            CRS_ENV = np.arange(0,len(ObsData))
            
            QFS = ObsData['QFS'][CRS_ENV]
            QFT = ObsData['QFT'][CRS_ENV]
            #lim = QFS<3
            
            #lim = 0
            if 'SB_AWI' in ifile or 'ASPeCt' in ifile:
                QFS = np.ones(len(QFS))*3
            if filt:
                index = (QFS<3) & (QFT<3)
            else:
                index = QFS<4
                if name=='IMB-CRREL':
                    years = np.array([dt.datetime.strptime(dd.decode('utf8'), '%Y-%m-%dT%H:%M:%S').year for dd in ObsData['date'][CRS_ENV]])
                    index = years != 2017

            obsDate = ObsData['date'][CRS_ENV][index]
            obslat  = ObsData['lat'][CRS_ENV][index]
            obslon  = ObsData['lon'][CRS_ENV][index]          
            obsSIF  = ObsData['obsSIF'][CRS_ENV][index]
            obsSD   = ObsData['obsSD'][CRS_ENV][index]
            obsSIT  = ObsData['obsSIT'][CRS_ENV][index]
            satSIF  = ObsData['satSIF'][CRS_ENV][index]
            satSIT  = ObsData['satSIT'][CRS_ENV][index]
            satSD   = ObsData['satSD'][CRS_ENV][index]
            obsSITunc = ObsData['obsSIT_unc'][CRS_ENV][index]
            obsSITstd = ObsData['obsSIT_std'][CRS_ENV][index]
            obsSITln = ObsData['obsSIT_ln'][CRS_ENV][index]

            #print(obsSITln)
            obsSIFunc = ObsData['obsFRB_unc'][CRS_ENV][index]
            obsSIFstd = ObsData['obsFRB_std'][CRS_ENV][index]
            obsSDunc = ObsData['obsSD_unc'][CRS_ENV][index]
            obsSDstd = ObsData['obsSD_std'][CRS_ENV][index]
            satSITunc = ObsData['satSIT_unc'][CRS_ENV][index]
            satSITstd = ObsData['satSIT_std'][CRS_ENV][index]
            satSIFunc = ObsData['satFRB_unc'][CRS_ENV][index]
            satSIFstd = ObsData['satFRB_std'][CRS_ENV][index]
            satSDunc = ObsData['satSD_unc'][CRS_ENV][index]
            satSDstd = ObsData['satSD_std'][CRS_ENV][index]

            if 'OIB' in name and 'SH' in name or 'AEM' in name: # convert sea ice frb to total frb
                satSIF = satSIF + satSD
            elif 'OIB' in ifile:
                obsSIF = obsSIF - obsSD
                #satSIF = satSIF + satSD
            if 'AEM-AWI' in ifile or 'HEM' in ifile:
                satSIT = satSIT + satSD # AEM measures SIT+SD
            
            m=~np.isnan(obsSD)
            n=np.logical_and(~np.isnan(obsSIT), obsSITln>0)
            o=~np.isnan(obsSIF)
            
            # assign non nan data to Data object
            if any(n) and name!='NP' and 'CS2' not in name:
                # if name == 'AEM-AWI' or name =='MOSAICHEM':
                #     name = 'AEM'
                
                d_SIT.lat.append(obslat[n])
                d_SIT.lon.append(obslon[n])
                d_SIT.obsSIT.append(obsSIT[n])
                d_SIT.satSIT.append(satSIT[n])
                d_SIT.obsSITunc.append(obsSITunc[n])
                d_SIT.satSITunc.append(satSITunc[n])
                d_SIT.obsSITstd.append(obsSITstd[n])
                d_SIT.satSITstd.append(satSITstd[n])
                d_SIT.name.append(name)
                d_SIT.date.append(obsDate[n])
                d_SIT.HS.append(HS)
                d_SIT.c.append(c)
            
                
            if any(m) and name!='NP' and 'CS2' not in name:
                d_SD.lat.append(obslat[m])
                d_SD.lon.append(obslon[m])
                d_SD.obsSD.append(obsSD[m])
                d_SD.satSD.append(satSD[m])
                d_SD.obsSDunc.append(obsSDunc[m])
                
                
                
                
                d_SD.satSDunc.append(satSDunc[m])
                d_SD.obsSDstd.append(obsSDstd[m])
                d_SD.satSDstd.append(satSDstd[m])
                d_SD.name.append(name)
                d_SD.date.append(obsDate[m])
                d_SD.HS.append(HS)
                d_SD.c.append(c)

            if any(o) and name!='NP':
                # change name and color
                if name=='OIB':
                    name='OIB-ENV'
                elif name=='OIBCS2':
                    name='OIB-CS2'
                    c='dimgray'
                elif name=='OIBSHSH':
                    name='OIB-ENV-SH'
                

                d_SIF.lat.append(obslat[o])
                d_SIF.lon.append(obslon[o])
                d_SIF.obsSIF.append(obsSIF[o])
                d_SIF.satSIF.append(satSIF[o])
                d_SIF.obsSIFunc.append(obsSIFunc[o])
                d_SIF.satSIFunc.append(satSIFunc[o])
                d_SIF.obsSIFstd.append(obsSIFstd[o])
                d_SIF.satSIFstd.append(satSIFstd[o])
                d_SIF.name.append(name)
                d_SIF.date.append(obsDate[o])
                d_SIF.HS.append(HS)
                d_SIF.c.append(c)

                
        return [d_SID, d_SIT, d_SD, d_SIF]

def Append_data_nc(d_SID, d_SIT, d_SD, d_SIF, ifile, name, c, directory):
    """
    Appends data from a .nc (NetCDF) file to corresponding Data objects.
    """

    # Open the NetCDF file
    dataset = xr.open_dataset(os.path.join(directory, ifile))
    
    #print(dataset)
    if 'SID' in ifile:

            # Extract variables from the NetCDF file
            obsSID = dataset['obsSID'].to_numpy()
            #print(obsSID)
            satSID = dataset['satSID'].to_numpy()
            obsSIDunc = dataset['obsSIDunc'].to_numpy()
            #print(obsSIDunc)
            obsSIDstd = dataset.variables['obsSIDstd'].to_numpy()
            satSIDunc = dataset.variables['satSIDunc'].to_numpy()
            satSIDstd = dataset.variables['satSIDstd'].to_numpy()
            obslat = dataset.variables['lat'].to_numpy()
            obslon = dataset.variables['lon'].to_numpy()
            obsDate = dataset.variables['date'].to_numpy()#.filled(fill_value=np.nan)
    
            m = ~np.isnan(obsSID)
    
            if any(m):
                d_SID.lat.append(obslat[m])
                d_SID.lon.append(obslon[m])
                d_SID.obsSID.append(obsSID[m])
                d_SID.obsSIDunc.append(obsSIDunc[m])
                d_SID.obsSIDstd.append(obsSIDstd[m])
                d_SID.satSID.append(satSID[m])
                d_SID.satSIDunc.append(satSIDunc[m])
                d_SID.satSIDstd.append(satSIDstd[m])
                d_SID.name.append(name)
                d_SID.date.append(obsDate[m])
                d_SID.c.append(c)
                d_SID.HS.append(HS)

    elif 'SIT' in ifile or 'SD' in ifile or 'FRB' in ifile:
        # Extract variables from the NetCDF file
        obsSIT = dataset['obsSIT'].to_numpy()
        obsSIT_ln = dataset['obsSITln'].to_numpy()
        print(obsSIT_ln)
        satSIT = dataset['satSIT'].to_numpy()
        obsSD = dataset['obsSD'].to_numpy()
        satSD = dataset['satSD'].to_numpy()
        obsSIF = dataset['obsFRB'].to_numpy()
        satSIF = dataset['satFRB'].to_numpy()
        obslat = dataset['lat'].to_numpy()
        obslon = dataset['lon'].to_numpy()
        obsDate = dataset['date'].to_numpy()#.filled(fill_value=np.nan)

        obsSIT_unc = dataset['obsSITunc'].to_numpy()
        satSIT_unc = dataset['satSITunc'].to_numpy()
        obsSD_unc = dataset['obsSDunc'].to_numpy()
        satSD_unc = dataset['satSDunc'].to_numpy()
        obsSIF_unc = dataset['obsFRBunc'].to_numpy()
        satSIF_unc = dataset['satFRBunc'].to_numpy()
        obsSIT_std = dataset['obsSITstd'].to_numpy()
        satSIT_std = dataset['satSITstd'].to_numpy()
        obsSD_std = dataset['obsSDstd'].to_numpy()
        satSD_std = dataset['satSDstd'].to_numpy()
        obsSIF_std = dataset['obsFRBstd'].to_numpy()
        satSIF_std = dataset['satFRBstd'].to_numpy()

        if 'OIB' in name and 'SH' in name or 'AEM' in name: # convert sea ice frb to total frb
            satSIF = satSIF + satSD
        elif 'OIB' in ifile:
            obsSIF = obsSIF - obsSD
        if 'AEM-AWI' in ifile:
            satSIT = satSIT + satSD # AEM measures SIT+SD
            
        m = ~np.isnan(obsSD)
        n = ~np.isnan(obsSIT)
        o = ~np.isnan(obsSIF)
        
        if any(n) and name != 'NP' and 'CS2' not in name:
            d_SIT.lat.append(obslat[n])
            d_SIT.lon.append(obslon[n])
            d_SIT.obsSIT.append(obsSIT[n])
            d_SIT.satSIT.append(satSIT[n])
            d_SIT.obsSITunc.append(obsSIT_unc[n])
            d_SIT.satSITunc.append(satSIT_unc[n])
            d_SIT.obsSITstd.append(obsSIT_std[n])
            d_SIT.satSITstd.append(satSIT_std[n])
            d_SIT.name.append(name)
            d_SIT.date.append(obsDate[n])
            d_SIT.HS.append(HS)
            d_SIT.c.append(c)

        if any(m) and name != 'NP' and 'CS2' not in name:
            d_SD.lat.append(obslat[m])
            d_SD.lon.append(obslon[m])
            d_SD.obsSD.append(obsSD[m])
            d_SD.satSD.append(satSD[m])
            d_SD.obsSDunc.append(obsSD_unc[m])
            d_SD.satSDunc.append(satSD_unc[m])
            d_SD.obsSDstd.append(obsSD_std[m])
            d_SD.satSDstd.append(satSD_std[m])
            d_SD.name.append(name)
            d_SD.date.append(obsDate[m])
            d_SD.HS.append(HS)
            d_SD.c.append(c)

        if any(o) and name != 'NP':
            if name == 'OIB':
                name = 'OIB-ENV'
            elif name == 'OIBCS2':
                name = 'OIB-CS2'
                c = 'dimgray'
            elif name == 'OIBSH':
                name = 'OIB-ENV-SH'

            d_SIF.lat.append(obslat[o])
            d_SIF.lon.append(obslon[o])
            d_SIF.obsSIF.append(obsSIF[o])
            d_SIF.satSIF.append(satSIF[o])
            d_SIF.obsSIFunc.append(obsSIF_unc[o])
            d_SIF.satSIFunc.append(satSIF_unc[o])
            d_SIF.obsSIFstd.append(obsSIF_std[o])
            d_SIF.satSIFstd.append(satSIF_std[o])
            d_SIF.name.append(name)
            d_SIF.date.append(obsDate[o])
            d_SIF.HS.append(HS)
            d_SIF.c.append(c)

    dataset.close()
    
    return d_SID, d_SIT, d_SD, d_SIF

def combine_existing_figures_as_images(map_fig, scatter_fig, hist_fig, savepath=None):
    """
    Combine three matplotlib figures horizontally by converting them to images.

    Parameters:
        map_fig: matplotlib Figure (e.g. map plot)
        scatter_fig: matplotlib Figure (e.g. scatter plot)
        hist_fig: matplotlib Figure (e.g. histogram plot)
        savepath: optional file path to save the combined figure

    Returns:
        combined matplotlib Figure object
    """
    def fig_to_array(fig):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
        buf.seek(0)
        img = Image.open(buf)
        return np.array(img)

    # Convert each figure to image arrays
    img_map = fig_to_array(map_fig)
    img_scatter = fig_to_array(scatter_fig)
    img_hist = fig_to_array(hist_fig)

    # Create a new figure with 1 row and 3 columns
    fig, axs = plt.subplots(1, 3, figsize=(36, 12))

    # Hide axes and show images
    for ax in axs:
        ax.axis('off')

    axs[0].imshow(img_map)
    axs[1].imshow(img_scatter)
    axs[2].imshow(img_hist)

    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, bbox_inches='tight')
        print(f"Combined figure saved to: {savepath}")
        plt.close(fig)
    else:
        plt.show()

    return fig

def combine_images_vertically(img1_path, img2_path, output_path, dpi=300):
    # Open the images
    img1 = Image.open(img1_path)
    img2 = Image.open(img2_path)

    # Ensure both are in RGB mode
    img1 = img1.convert("RGB")
    img2 = img2.convert("RGB")

    # Resize to match width if needed (optional)
    if img1.width != img2.width:
        new_width = max(img1.width, img2.width)
        img1 = img1.resize((new_width, int(img1.height * new_width / img1.width)))
        img2 = img2.resize((new_width, int(img2.height * new_width / img2.width)))

    # Create new image with combined height
    total_height = img1.height + img2.height
    combined_img = Image.new('RGB', (img1.width, total_height), (255, 255, 255))

    # Paste both images
    combined_img.paste(img1, (0, 0))
    combined_img.paste(img2, (0, img1.height))

    # Save with high resolution
    combined_img.save(output_path, dpi=(dpi, dpi))
    print(f"Saved high-resolution combined image to: {output_path}")

#%% MAIN
# Specify the directory containing the data files.
#directory = "C:/Users/Ida Olsen/Documents/work/RRDPp/satellite/Final_files/final/Arctic/"
directory = "C:/Users/Ida Olsen/Documents/work/RRDPp/satellite/Final_files/collocated_files-20250519/collocated_files"
# List all files in the specified directory.
files = os.listdir(directory)

sat = 'CS2'  # Set the satellite variable (either ENV for Envisat or CS2 for CryoSat2)
filt = False

# Initialize lists to hold data objects and boolean flags.
objects = []
bool_list = []

# Create empty data objects for different variable types.
d_SID = Data(sat)  # Sea Ice Draft
d_SID.var = 'SID'  # Set the variable type to 'SID'
d_SIT = Data(sat)  # Sea Ice Thickness
d_SIT.var = 'SIT'  # Set the variable type to 'SIT'
d_SD = Data(sat)   # Snow Depth
d_SD.var = 'SD'    # Set the variable type to 'SD'
d_SIF = Data(sat)  # Freeboard
d_SIF.var = 'FRB'  # Set the variable type to 'FRB'

# Loop through the files and append data to the respective objects.
for ifile, i in zip(files, range(len(files))):
    # Check if the file is a data file and matches the satellite criteria.
    if (ifile.endswith('.dat') and 
        (ifile.startswith('ESACCI') and (sat in ifile.replace('-', '')) or 
         (('CS2' in ifile) and ('ENV' in sat) and ('OIB' in ifile)))):
        # Get the name from the filename.
        name = Get_name(ifile)
        # Get the associated color for the dataset.
        c = Get_color(name)

        # Determine the hemisphere based on the name.
        HS = 'NH'  # Default to Northern Hemisphere
        if 'ASPeCt' in name or 'SH' in name or 'ULS' in name:
            HS = 'SH'  # Set to Southern Hemisphere if criteria are met
        
        # Append the data from the file to the respective Data objects.
        d_SID, d_SIT, d_SD, d_SIF = Append_data(d_SID, d_SIT, d_SD, d_SIF, ifile, name, c, directory=directory, filt=filt)


if sat=='CS2':
    d_SD.name[d_SD.name=='AEM: (AWI + MOSAIC)'] = 'AEM-AWI'
# Generate histograms for Sea Ice Thickness and Snow Depth.
fig_hist_SID, ax = Histograms(d_SID, sat)  # Uncomment to generate histograms for Sea Ice Draft
fig_hist_SIT, ax = Histograms(d_SIT, sat)  # Generate histogram for Sea Ice Thickness
if sat=='ENV':
    fig_hist_SIF, ax =Histograms(d_SIF, sat)  # Uncomment to generate histogram for Freeboard if satellite is ENV
fig_hist_SD, ax = Histograms(d_SD, sat)  # Generate histogram for Snow Depth

# Sort data based on names for all data types.
d_SIT.sort()
d_SD.sort()
d_SID.sort()
if sat == 'ENV':
     d_SIF.sort()

# Define the output file for statistics.
ofile = f'stats_{sat}_out.txt'
# Get statistics about the data and save them to the specified output file.
stats(d_SID, sat, ofile)  # Uncomment to get stats for Sea Ice Draft
# stats(d_SIT, sat, ofile)  # Get statistics for Sea Ice Thickness
# stats(d_SD, sat, ofile)  # Get statistics for Snow Depth
# stats(d_SIF, sat, ofile)  # Uncomment to get stats for Freeboard

# # Generate scatter plots for each data type.
fig_scatter_SID, ax = scatter(d_SID)  # Scatter plot for Sea Ice Draft
fig_scatter_SIT, ax = scatter(d_SIT)  # Scatter plot for Sea Ice Thickness
fig_scatter_SIF, ax = scatter(d_SIF)  # Scatter plot for Freeboard
fig_scatter_SD, ax = scatter(d_SD)   # Scatter plot for Snow Depth

#%% Polar plots
bool_list_FRB = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIF.name])
savename =  'FRB_data_'+sat+'.png'
title = 'FRB Validation Data' #' for ' + sat
if sat=='ENV':
    fig_geo_SIF, ax = pp.plot_all(d_SIF.lat[bool_list_FRB], d_SIF.lon[bool_list_FRB], d_SIF.obsSIF[bool_list_FRB], title=title, ylabel=d_SIF.name[bool_list_FRB], c=d_SIF.c[bool_list_FRB], savename=savename, s=40, sat=d_SIF.sat)

bool_list_SIT = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIT.name])
savename =  'SIT_data_'+sat+'.png'
title = 'SIT Validation Data for ' + sat
fig_geo_SIT, ax = pp.plot_all(d_SIT.lat[bool_list_SIT], d_SIT.lon[bool_list_SIT], d_SIT.obsSIT[bool_list_SIT], title=title, ylabel=d_SIT.name[bool_list_SIT], c=d_SIT.c[bool_list_SIT], savename=savename, s=40, sat=d_SIT.sat)


bool_list_SD = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SD.name])
savename = 'SD_data_'+sat+'.png'
title = 'SD Validation Data for '  + sat
fig_geo_SD, ax = pp.plot_all(d_SD.lat[bool_list_SD], d_SD.lon[bool_list_SD], d_SD.obsSD[bool_list_SD], title=title, ylabel=d_SD.name[bool_list_SD], c=d_SD.c[bool_list_SD], savename=savename, s=40, sat=d_SD.sat)

bool_list_SID = np.invert(['ASPeCt'==name or 'SH' in name or 'ULS' in name for name in d_SID.name])
savename = 'SID_data_'+sat+'.png'
s = 500
title = 'SID Validation Data for '  + sat
fig_geo_SID, ax = pp.plot_all(d_SID.lat[bool_list_SID], d_SID.lon[bool_list_SID], d_SID.obsSID[bool_list_SID], title=title,  c=d_SID.c[bool_list_SID], ylabel=d_SID.name[bool_list_SID], savename=savename, s=s, sat=d_SID.sat)

objects += [d_SD, d_SIT, d_SID]
bool_list += [bool_list_SD, bool_list_SIT, bool_list_SID]

#%%  ANTARCTIC
# bool_list_SIT = ['ASPeCt' in name or 'SH' in name for name in d_SIT.name]
# bool_list_SD =  ['SH' in name for name in d_SD.name] #'ASPeCt' in name or 
# bool_list_SID = ['ASPeCt' in name or 'SH' in name or 'AWI-ULS' in name for name in d_SID.name]
# bool_list_SIF = ['ASPeCt' in name or 'SH' in name for name in d_SIF.name]

# savename =  'SH Validation_data.png'
# title = "SH: Validation Data for CS2 & ENV"

# lats = np.concatenate([d_SIT.lat[bool_list_SIT], d_SD.lat[bool_list_SD], d_SID.lat[bool_list_SID], np.array(d_SIF.lat)[bool_list_SIF]])
# lons = np.concatenate([d_SIT.lon[bool_list_SIT], d_SD.lon[bool_list_SD], d_SID.lon[bool_list_SID], np.array(d_SIF.lon)[bool_list_SIF]])
# obs = np.concatenate([d_SIT.obsSIT[bool_list_SIT], d_SD.obsSD[bool_list_SD], d_SID.obsSID[bool_list_SID], np.array(d_SIF.obsSIF)[bool_list_SIF]])
# names = np.concatenate([d_SIT.name[bool_list_SIT], d_SD.name[bool_list_SD], d_SID.name[bool_list_SID], np.array(d_SIF.name)[bool_list_SIF]])
# names = np.array([name.replace('-SH','')+'-CS2' if ('ENV' not in name) else name.replace('-SH', '') for name in names])
# colors = np.concatenate([d_SIT.c[bool_list_SIT], d_SD.c[bool_list_SD], d_SID.c[bool_list_SID], np.array(d_SIF.c)[bool_list_SIF]])
# pp.plot_all(lats, lons, obs, title=title, ylabel=names, s=5, c=colors, savename=savename, NP=False)


#%% Combine plots

combine_existing_figures_as_images(fig_geo_SID, fig_scatter_SID, fig_hist_SID, savepath=f'{d_SID.var}_{sat}_combined_{filt}.png')
combine_existing_figures_as_images(fig_geo_SIT, fig_scatter_SIT, fig_hist_SIT, savepath=f'{d_SIT.var}_{sat}_combined_{filt}.png')
combine_existing_figures_as_images(fig_geo_SD, fig_scatter_SD, fig_hist_SD, savepath=f'{d_SD.var}_{sat}_combined_{filt}.png')
if sat=='ENV':
    combine_existing_figures_as_images(fig_geo_SIF, fig_scatter_SIF, fig_hist_SIF, savepath=f'{d_SIF.var}_{sat}_combined_{filt}.png')


# #%% Save final figure
# combine_images_vertically(f'{d_SIT.var}_CS2_combined_{filt}.png', f'{d_SIT.var}_ENV_combined_{filt}.png', f'{d_SIT.var}_combined_{filt}.png', dpi=300)
# combine_images_vertically(f'{d_SD.var}_CS2_combined_{filt}.png', f'{d_SD.var}_ENV_combined_{filt}.png', f'{d_SD.var}_combined_{filt}.png', dpi=300)
# combine_images_vertically(f'{d_SID.var}_CS2_combined_{filt}.png', f'{d_SID.var}_ENV_combined_{filt}.png', f'{d_SID.var}_combined_{filt}.png', dpi=300)