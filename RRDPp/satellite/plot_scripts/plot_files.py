# -*- coding: utf-8 -*-
"""
Created on Mon Nov 28 15:31:46 2022

@author: Ida Olsen
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 21 18:31:18 2021

@author: s174020

Make histograms of data to check distribution.
"""

# Built-in packages
import sys
import os
import datetime as dt

# Append custom plot scripts directory to system path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.getcwd())), 'plot_scripts'))

# Third-party packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
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
    
    plt.show()  # Display the plot
    
#%% Scatterplots
def scatter(self):
    # Select the observed and satellite variables based on the 'var' attribute
    if self.var == 'SID':
        obsVar = self.obsSID
        satVar = self.satSID
        # Calculate observation and satellite errors
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
    
    # Initialize lists to store non-NaN observations and errors
    obsvar = []
    satvar = []
    obserr = []
    saterr = []
    
    # Create a figure for the scatter plot
    fig = plt.figure(figsize=(30, 30))
    plt.grid(alpha=0.5, linewidth=2)
    
    # Print the variable being processed
    print(f'variable : {self.var}')
    print(f'***********************')
    
    # Loop through each observation variable and satellite variable
    for i, name in zip(range(len(obsVar)), self.name):
        # Optional: Skip specific satellite names (commented out)
        # if self.name[i] == 'OIB' and self.sat == 'CS2':
        #     obsVar[i] = obsVar[i][::5]
        #     satVar[i] = satVar[i][::5]
        
        # Only plot for NH data
        if self.name[i] != 'ASPeCt' and 'SH' not in self.name[i] and 'AWI-ULS' not in self.name[i]:
            # Create a scatter plot for observed vs. satellite values
            plt.scatter(obsVar[i], satVar[i], color=self.c[i], s=40, alpha=0.4)
            
            # Calculate the Pearson correlation coefficient
            corr = np.round(pearsonr(obsVar[i][~np.isnan(satVar[i])], satVar[i][~np.isnan(satVar[i])]), 2)
            print(self.name[i])
            
            # Fit a linear model to the observed and satellite data
            w, b = odr_fit(obsVar[i][~np.isnan(satVar[i])], satVar[i][~np.isnan(satVar[i])], 
                           sx=obsErr[i][~np.isnan(satVar[i])], sy=satErr[i][~np.isnan(satVar[i])])
            
            # Calculate and print the RMSE of the fit
            print('RMSE', np.round(np.sqrt(mean_squared_error(satVar[i][~np.isnan(satVar[i])], 
                             (obsVar[i][~np.isnan(satVar[i])]*w + b))), 2))
            
            # Plot the fitted line on the scatter plot
            plt.plot(range(-1, 8), range(-1, 8)*w + b, color=self.c[i], linewidth=2, 
                     label=f'{name}, R:{corr[0]}')
            
            # Print the R² value of the fit
            print('R2', np.round(r2_score(satVar[i][~np.isnan(satVar[i])], 
                                             (obsVar[i][~np.isnan(satVar[i])]*w + b)), 2))
            print(f' wb: {w} {b}')
            print(f'len: {len(satVar[i][~np.isnan(satVar[i])])}')
        
        # Collect non-NaN observations and errors for overall analysis
        if self.name[i] != 'ASPeCt' and 'SH' not in self.name[i] and 'AWI-ULS' not in self.name[i]:
            obsvar.append(obsVar[i][~np.isnan(satVar[i])])
            satvar.append(satVar[i][~np.isnan(satVar[i])])
            obserr.append(obsErr[i][~np.isnan(satVar[i])])
            saterr.append(satErr[i][~np.isnan(satVar[i])])
    
    # Concatenate all collected data for combined analysis
    x = np.concatenate(obsvar)
    y = np.concatenate(satvar)
    sx = np.concatenate(obserr)
    sy = np.concatenate(saterr)
    
    # Calculate the Pearson correlation for the combined dataset
    corr = np.round(pearsonr(x[~np.isnan(y)], y[~np.isnan(y)]), 2)
    
    # Fit a linear model to the combined dataset
    w, b = odr_fit(x, y, sx, sy)
    
    # Plot the fitted line for the combined dataset
    plt.plot(range(-1, 8), range(-1, 8)*w + b, color='r', label=f'Combined, R:{corr[0]}')
    
    # Plot the 1-to-1 line and set axes limits based on the variable type
    if self.var == 'SIT' or self.var == 'SID':
        plt.plot([-0.10, 8], 
                 [-0.10, 8], c='k', label='1 to 1 line')
        plt.xlim([-0.10, 7])
        plt.ylim([-0.10, 7])
    elif self.var == 'SD':
        plt.plot([-0.1, 1], [-0.1, 1], c='k')
        plt.xlim([-0.1, 1])
        plt.ylim([-0.1, 1])
    else:
        plt.plot([-0.1, 1], [-0.1, 1], c='k')
        plt.xlim([-0.10, 1])
        plt.ylim([-0.10, 1])

    # Get the legend object and adjust line widths for better visibility
    leg = plt.legend(fontsize=40, loc='upper left')
    for line in leg.get_lines():
        line.set_linewidth(4.0)

    # Set the axis labels
    plt.xlabel('Reference ' + self.var + ' [m]', fontsize=48)
    plt.ylabel('Satellite ' + self.var + ' [m]', fontsize=48)

    # Save the figure
    plt.savefig(f'scatter_{self.sat}_{self.var}.png', bbox_inches='tight')
    plt.show()
    plt.close()

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
def Histograms(self, sat):
    # Define bin edges for histograms based on the selected variable type
    if self.var == 'SID':
        obsVar = self.obsSID
        satVar = self.satSID
        bins = np.arange(-0.1, 3, 0.10)  # Bins for SID
    elif self.var == 'SIT':
        obsVar = self.obsSIT
        satVar = self.satSIT
        bins = np.arange(0, 7, 0.10)  # Bins for SIT
    elif self.var == 'SD':
        obsVar = self.obsSD
        satVar = self.satSD
        bins = np.arange(0, 1, 0.03)  # Bins for SD
    elif self.var == 'FRB':
        obsVar = self.obsSIF
        satVar = self.satSIF
        bins = np.arange(-0.5, 1, 0.03)  # Bins for FRB
        


    # Create subplots for histograms, sharing the x-axis
    fig, ax1 = plt.subplots(len(obsVar), figsize=(30, 30), sharex=True)

    # Loop through each observed variable to create individual histograms
    for i, name in zip(range(len(obsVar)), self.name):
        ax1[i].set_ylabel(name, fontsize=48)  # Set y-axis label for each subplot
        ax1[i].set_yticks([])  # Hide y-ticks for clarity
        ax1[i].tick_params(axis="x", labelsize=48)  # Set font size for x-ticks
        ax1[i].grid()  # Enable grid for better readability

        ax2 = ax1[i].twinx()  # Create a twin y-axis for overlaying histograms
        
        # Conditional plotting for specific datasets with additional formatting
        if name == 'ASPeCt' or 'SH' in name or name == 'AWI-ULS':
            # Plot histograms for observed and satellite data
            ax2.hist(obsVar[i], bins=bins, rwidth=0.85, linewidth=4, 
                     linestyle='-', edgecolor='k', alpha=0.7, stacked=True)
            ax2.hist(satVar[i], bins=bins, rwidth=0.85, linewidth=4, 
                     linestyle='-', edgecolor='k', alpha=0.7)
        else:
            ax2.hist(obsVar[i], bins=bins, rwidth=0.85, alpha=0.7)  # Normal histogram for observed data
            ax2.hist(satVar[i], bins=bins, rwidth=0.85, alpha=0.7)  # Normal histogram for satellite data
        
        ax2.tick_params(axis="y", labelsize=48)  # Set font size for y-ticks
        
        # Draw vertical lines indicating the means of observed and satellite data
        ax2.axvline(x=np.nanmean(obsVar[i]), c='r', linewidth=8, label='obs mean')  # Observed mean
        ax2.axvline(x=np.nanmean(satVar[i]), c='k', linewidth=8, label='sat mean')  # Satellite mean
        ax2.grid()  # Enable grid on the twin axis

    # Set x-axis label for the bottom subplot
    ax1[i].set_xlabel(f'{self.var} [m]', fontsize=48)

    # Adjust the satellite name for FRB variable
    if self.var == 'FRB':
        sat = 'ENV & CS2'
    else:
        sat = sat  # Retain the provided satellite name
    
    # Create legend elements for the plot
    legend_elements = [
        mpatches.Patch(color='tab:blue', alpha=0.7, label='obs ' + self.var + ' [m]'),
        mpatches.Patch(color='tab:orange', alpha=0.7, label=sat + ' ' + self.var + ' [m]'),
        mpatches.Patch(facecolor='tab:blue', alpha=0.7, edgecolor='k', linewidth=1.2, 
                       linestyle='-', label='obs SH. ' + self.var + ' [m]'),
        mpatches.Patch(facecolor='tab:orange', alpha=0.7, edgecolor='k', linewidth=1.2, 
                       linestyle='-', label=sat + ' SH. ' + self.var + ' [m]'),
        mlines.Line2D([], [], color='r', linewidth=3, label='obs mean [m]'),
        mlines.Line2D([], [], color='k', linewidth=3, label='sat mean [m]')
    ]
    
    # Add a super label to the figure
    fig.supylabel('count of observations in each bin', fontsize=60, position=(0.97, 0.5))
    
    # Save the figure to a PNG file
    plt.savefig(self.var + '_hist_' + sat + '.png', bbox_inches='tight')
    
    plt.show()  # Display the plot
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
    name = name.replace("BSH", "B-SH") 
    name = name.replace("sorted", "") 
    name = name.replace("ULS", "-ULS") 
    name = name.replace("ENV", "-ENV")
    name = name.replace("IMB", "IMB-CRREL")
    name = name.replace("NPI", "NPI-FS")
    if 'ULS' not in name:
        name = name.replace("AWI", "-AWI")
        
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
    elif 'IMB' in name:
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
        c=colors[0]
        
    return c

def Append_data(d_SID, d_SIT, d_SD, d_SIF, ifile, name, c,directory):
        ## append data
        if 'SID' in ifile:
            # define observation names
            ObsNames=['date','lat','lon','obsSID', 'satSID', 'obsSIDstd','obsSIDln','obsSIDunc',
                      'satSIDstd','satSIDln', 'satSIDunc']
        
            # Read data
            ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None,skip_header=1,names=ObsNames)
        
            CRS_ENV = np.arange(0,len(ObsData))
            obsDate = ObsData['date'][CRS_ENV]
            obslat  = ObsData['lat'][CRS_ENV]
            obslon  = ObsData['lon'][CRS_ENV]
            obsSID  = ObsData['obsSID'][CRS_ENV]
            obsSIDunc = ObsData['obsSIDunc'][CRS_ENV]
            obsSIDstd = ObsData['obsSIDstd'][CRS_ENV]
            satSID  = ObsData['satSID'][CRS_ENV]
            satSIDunc = ObsData['satSIDunc'][CRS_ENV]
            satSIDstd = ObsData['satSIDstd'][CRS_ENV]
           
            m=~np.isnan(obsSID)

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
                      "obsFRB_ln", "obsFRB_unc", "satSD_std","satSD_ln", "satSD_unc","satSIT_std","satSIT_ln", "satSIT_unc","satFRB_std",
                      "satFRB_ln", "satFRB_unc"]
            # Read data
            ObsData  = np.genfromtxt(os.path.join(directory, ifile),dtype=None, skip_header=1, names=ObsNames)
            
            # identify copy rows
            ObsData = np.unique(ObsData, axis=0)
 
            CRS_ENV = np.arange(0,len(ObsData))
            obsDate = ObsData['date'][CRS_ENV]
            obslat  = ObsData['lat'][CRS_ENV]
            obslon  = ObsData['lon'][CRS_ENV]         
            
            obsSIF  = ObsData['obsSIF'][CRS_ENV]
            obsSD   = ObsData['obsSD'][CRS_ENV]
            obsSIT  = ObsData['obsSIT'][CRS_ENV]
            satSIF  = ObsData['satSIF'][CRS_ENV]
            satSIT  = ObsData['satSIT'][CRS_ENV]
            satSD   = ObsData['satSD'][CRS_ENV]
            obsSITunc = ObsData['obsSIT_unc'][CRS_ENV]
            obsSITstd = ObsData['obsSIT_std'][CRS_ENV]
            obsSIFunc = ObsData['obsFRB_unc'][CRS_ENV]
            obsSIFstd = ObsData['obsFRB_std'][CRS_ENV]
            obsSDunc = ObsData['obsSD_unc'][CRS_ENV]
            obsSDstd = ObsData['obsSD_std'][CRS_ENV]
            satSITunc = ObsData['satSIT_unc'][CRS_ENV]
            satSITstd = ObsData['satSIT_std'][CRS_ENV]
            satSIFunc = ObsData['satFRB_unc'][CRS_ENV]
            satSIFstd = ObsData['satFRB_std'][CRS_ENV]
            satSDunc = ObsData['satSD_unc'][CRS_ENV]
            satSDstd = ObsData['satSD_std'][CRS_ENV]

            if 'OIB' in name and 'SH' in name or 'AEM' in name: # convert sea ice frb to total frb
                satSIF = satSIF + satSD
            elif 'OIB' in ifile:
                obsSIF = obsSIF - obsSD
            if 'AEM-AWI' in ifile:
                satSIT = satSIT + satSD # AEM measures SIT+SD
            
            m=~np.isnan(obsSD)
            n=~np.isnan(obsSIT)
            o=~np.isnan(obsSIF)
            
            # assign non nan data to Data object
            if any(n) and name!='NP' and 'CS2' not in name:
                
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
                elif name=='OIBSH':
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

#%% MAIN
# Specify the directory containing the data files.
directory = "C:/Users/Ida Olsen/Documents/work/RRDPp/satellite/Final_files/final/Arctic/"
# List all files in the specified directory.
files = os.listdir(directory)

sat = 'ENV'  # Set the satellite variable (either ENV for Envisat or CS2 for CryoSat2)

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
    if (ifile.endswith('.nc') and 
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
        d_SID, d_SIT, d_SD, d_SIF = Append_data_nc(d_SID, d_SIT, d_SD, d_SIF, ifile, name, c, directory=directory)

# Generate histograms for Sea Ice Thickness and Snow Depth.
Histograms(d_SID, sat)  # Uncomment to generate histograms for Sea Ice Draft
Histograms(d_SIT, sat)  # Generate histogram for Sea Ice Thickness
if sat=='ENV':
      Histograms(d_SIF, sat)  # Uncomment to generate histogram for Freeboard if satellite is ENV
Histograms(d_SD, sat)  # Generate histogram for Snow Depth

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
stats(d_SIT, sat, ofile)  # Get statistics for Sea Ice Thickness
stats(d_SD, sat, ofile)  # Get statistics for Snow Depth
stats(d_SIF, sat, ofile)  # Uncomment to get stats for Freeboard

# Generate scatter plots for each data type.
scatter(d_SID)  # Scatter plot for Sea Ice Draft
scatter(d_SIT)  # Scatter plot for Sea Ice Thickness
scatter(d_SIF)  # Scatter plot for Freeboard
scatter(d_SD)   # Scatter plot for Snow Depth

#%% Polar plots
#bool_list_FRB = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIF.name])
# savename =  'FRB_data_'+sat+'.png'
# title = 'FRB Validation Data' #' for ' + sat
# if sat=='ENV':
#     pp.plot_all(d_SIF.lat[bool_list_FRB], d_SIF.lon[bool_list_FRB], d_SIF.obsSIF[bool_list_FRB], title=title, ylabel=d_SIF.name[bool_list_FRB], c=d_SIF.c[bool_list_FRB], savename=savename, s=40)

# bool_list_SIT = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SIT.name])
# savename =  'SIT_data_'+sat+'.png'
# title = 'SIT Validation Data for ' + sat
# pp.plot_all(d_SIT.lat[bool_list_SIT], d_SIT.lon[bool_list_SIT], d_SIT.obsSIT[bool_list_SIT], title=title, ylabel=d_SIT.name[bool_list_SIT], c=d_SIT.c[bool_list_SIT], savename=savename, s=40)


# bool_list_SD = np.invert(['ASPeCt'==name or 'SH' in name for name in d_SD.name])
# savename = 'SD_data_'+sat+'.png'
# title = 'SD Validation Data for '  + sat
# pp.plot_all(d_SD.lat[bool_list_SD], d_SD.lon[bool_list_SD], d_SD.obsSD[bool_list_SD], title=title, ylabel=d_SD.name[bool_list_SD], c=d_SD.c[bool_list_SD], savename=savename, s=40)

# bool_list_SID = np.invert(['ASPeCt'==name or 'SH' in name or 'ULS' in name for name in d_SID.name])
# savename = 'SID_data_'+sat+'.png'
# s = 500
# title = 'SID Validation Data for '  + sat
# pp.plot_all(d_SID.lat[bool_list_SID], d_SID.lon[bool_list_SID], d_SID.obsSID[bool_list_SID], title=title,  c=d_SID.c[bool_list_SID], ylabel=d_SID.name[bool_list_SID], savename=savename, s=s)

# objects += [d_SD, d_SIT, d_SID]
# bool_list += [bool_list_SD, bool_list_SIT, bool_list_SID]

# #%%  ANTARCTIC
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
    