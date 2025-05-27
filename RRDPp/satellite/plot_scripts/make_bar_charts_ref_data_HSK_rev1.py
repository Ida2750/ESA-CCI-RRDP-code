"""
Title: Sea Ice Data Analysis
Author: Ida Olsen
Description: This script performs analysis on Sea Ice data, including Sea Ice Thickness (SIT), Snow Depth (SD), Freeboard (FRB), and Sea Ice Drift (SID).
Updated to make eventsplots of the data, Henriette Skourup, 2024
"""

# -- File info -- #
__author__ = 'Ida Olsen'
__contributors__ = ''
__contact__ = ['iblol@space.dtu.dk']
__version__ = '0'
__date__ = '2023-07-13'

# -- Built-in modules -- #
import os
import datetime as dt
from datetime import date, timedelta

# -- Third-part modules -- #
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class Data:
    def __init__(self):
        """
        Initialize Data class with empty lists to store data.
        """
        self.name = []
        self.var = ''
        self.date = []
        self.lat = []
        self.lon = []
        self.obs_sit = []
        self.obs_sd = []
        self.obs_sid = []
        self.obs_sif = []


def createList(r1, r2):
    return list(range(r1, r2+1))
     

def get_data(HS):
    # Loads data
#    directory = "C:/Users/Ida Olsen/Documents/work/RRDPp/FINAL/done_done_files/"
    directory = "C:/SICCI/Paper/DTUdata_nc_rev1/FINAL/"
    if HS != 'NH':
        directory += "Antarctic"
    
    files = os.listdir(directory)
 
    # Definition of reference observations
    d_oib = Data()
    d_oib.var = 'OIB'
    d_aem = Data()
    d_aem.var = 'AEM-AWI'
    d_hem = Data()
    d_hem.var = 'HEM-MOSAIC'
    d_hem_nl = Data()
    d_hem_nl.var = 'HEM-NL'
    d_assist = Data()
    d_assist.var = 'ASSIST'
    d_aspect = Data()
    d_aspect.var = 'ASPeCt'
    d_imb = Data()
    d_imb.var = 'IMB-CRREL'
    d_simba = Data()
    d_simba.var = 'SIMBA-MOSAIC'
    d_sb = Data()
    d_sb.var = 'SB-AWI'
    d_bgep = Data()
    d_bgep.var = 'BGEP'
    d_ulsawi = Data()
    d_ulsawi.var = 'AWI-ULS'
    d_trans = Data()
    d_trans.var = 'TRANSDRIFT'
    d_ulsnpi = Data()
    d_ulsnpi.var = 'NPI-FS'
    d_npeo = Data()
    d_npeo.var = 'NPEO'
    d_uls_nl = Data()
    d_uls_nl.var = 'ULS-NL'
    d_scicex = Data()
    d_scicex.var = 'SCICEX'

    for ifile in files:
        if ifile.startswith('ESACCI'):
            if 'OIB' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
#               print('SB ',ifile)
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sif = obs_data['FRB'][crs_env]
                k = ~np.isnan(obs_sif)
    
                if any(k):
                    d_oib.lat.append(obs_lat[k])
                    d_oib.lon.append(obs_lon[k])
                    d_oib.obs_sif.append(obs_sif[k])
                    d_oib.date.append(obs_date[k])

            elif 'SIT-AEM' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                l = ~np.isnan(obs_sit)
    
                if any(l):
                    d_aem.lat.append(obs_lat[l])
                    d_aem.lon.append(obs_lon[l])
                    d_aem.obs_sit.append(obs_sit[l])
                    d_aem.date.append(obs_date[l])
#                    print(d_aem.date)

            elif 'SIT-MOSAIC-HEM' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                l = ~np.isnan(obs_sit)
    
                if any(l):
                    d_hem.lat.append(obs_lat[l])
                    d_hem.lon.append(obs_lon[l])
                    d_hem.obs_sit.append(obs_sit[l])
                    d_hem.date.append(obs_date[l])
#                    print(d_aem.date)

            elif 'SIT-Nansen_legacy' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                l = ~np.isnan(obs_sit)
    
                if any(l):
                    d_hem_nl.lat.append(obs_lat[l])
                    d_hem_nl.lon.append(obs_lon[l])
                    d_hem_nl.obs_sit.append(obs_sit[l])
                    d_hem_nl.date.append(obs_date[l])
#                    print(d_aem.date)

            elif 'ASSIST' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                m = ~np.isnan(obs_sit)
    
                if any(m):
                    d_assist.lat.append(obs_lat[m])
                    d_assist.lon.append(obs_lon[m])
                    d_assist.obs_sit.append(obs_sit[m])
                    d_assist.date.append(obs_date[m])


            elif 'BGEP' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_bgep.lat.append(obs_lat[m])
                    d_bgep.lon.append(obs_lon[m])
                    d_bgep.obs_sid.append(obs_sid[m])
                    d_bgep.date.append(obs_date[m])

            elif 'NPI' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_ulsnpi.lat.append(obs_lat[m])
                    d_ulsnpi.lon.append(obs_lon[m])
                    d_ulsnpi.obs_sid.append(obs_sid[m])
                    d_ulsnpi.date.append(obs_date[m])

            elif 'TRANSDRIFT' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_trans.lat.append(obs_lat[m])
                    d_trans.lon.append(obs_lon[m])
                    d_trans.obs_sid.append(obs_sid[m])
                    d_trans.date.append(obs_date[m])

            elif 'NPEO' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_npeo.lat.append(obs_lat[m])
                    d_npeo.lon.append(obs_lon[m])
                    d_npeo.obs_sid.append(obs_sid[m])
                    d_npeo.date.append(obs_date[m])

            elif 'SID-Nansen_legacy' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_uls_nl.lat.append(obs_lat[m])
                    d_uls_nl.lon.append(obs_lon[m])
                    d_uls_nl.obs_sid.append(obs_sid[m])
                    d_uls_nl.date.append(obs_date[m])


            elif 'SCICEX' in ifile:
                # Load data for Sea Ice Drift (SID)
#                print('AEM ',ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                m = ~np.isnan(obs_sid)
    
                if any(m):
                    d_scicex.lat.append(obs_lat[m])
                    d_scicex.lon.append(obs_lon[m])
                    d_scicex.obs_sid.append(obs_sid[m])
                    d_scicex.date.append(obs_date[m])


            elif 'SB' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
 #               print('SB ',ifile)
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sd = obs_data['SD'][crs_env]
                m = ~np.isnan(obs_sd)
    
                if any(m):
                    d_sb.lat.append(obs_lat[m])
                    d_sb.lon.append(obs_lon[m])
                    d_sb.obs_sd.append(obs_sd[m])
                    d_sb.date.append(obs_date[m])

            elif 'CRREL-IMB' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
#                print('hej1')
#               print('SB ',ifile)
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                n = ~np.isnan(obs_sit)
    
                if any(n):
                    d_imb.lat.append(obs_lat[n])
                    d_imb.lon.append(obs_lon[n])
                    d_imb.obs_sit.append(obs_sit[n])
                    d_imb.date.append(obs_date[n])


            elif 'SIMBA' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                print('hej: ', ifile)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                n = ~np.isnan(obs_sit)
    
                if any(n):
                    d_simba.lat.append(obs_lat[n])
                    d_simba.lon.append(obs_lon[n])
                    d_simba.obs_sit.append(obs_sit[n])
                    d_simba.date.append(obs_date[n])



            elif 'ASPeCt' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
 #               print('SB ',ifile)
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sit = obs_data['SIT'][crs_env]
                n = ~np.isnan(obs_sit)
    
                if any(n):
                    d_aspect.lat.append(obs_lat[n])
                    d_aspect.lon.append(obs_lon[n])
                    d_aspect.obs_sit.append(obs_sit[n])
                    d_aspect.date.append(obs_date[n])

            elif 'AWI-ULS_Bias_False' in ifile: # or 'SD' in ifile or 'FRB' in ifile:
 #               print('SB ',ifile)
                # Load data for Sea Ice Thickness (SIT), Snow Depth (SD), and Freeboard (FRB)
                obs_data = np.genfromtxt(os.path.join(directory, ifile), dtype=None, names=True, encoding='Latin1')
                crs_env = np.arange(0, len(obs_data))
                obs_date = obs_data['date'][crs_env]
                obs_lat = obs_data['lat'][crs_env]
                obs_lon = obs_data['lon'][crs_env]
                obs_sid = obs_data['SID'][crs_env]
                n = ~np.isnan(obs_sid)
    
                if any(n):
                    d_ulsawi.lat.append(obs_lat[n])
                    d_ulsawi.lon.append(obs_lon[n])
                    d_ulsawi.obs_sid.append(obs_sid[n])
                    d_ulsawi.date.append(obs_date[n])


    return [d_oib, d_aem, d_hem, d_hem_nl, d_assist, d_bgep, d_trans, d_ulsnpi, d_npeo, d_uls_nl, d_scicex, d_imb, d_simba, d_sb, d_aspect, d_ulsawi]


def plot_distribution_bar(dataCDR,dataNH,dataSH, title, x_label):
#     """
#     Function to plot a stacked bar chart for the given data.

#     Parameters:
#     data (dict): A dictionary containing data to be plotted.
#     title (str): Title of the plot.
#     x_label (str): Label for the x-axis.
#
#    Returns:
#    None
#    """
#    print(dataCDR)
    # print(data)
    # print(len(data))
    # print(data.values())

#    for key in data:
# #        print(key)
#     # from matplotlib.ticker import MultipleLocator
#     # plt.rcParams.update({'font.size': 14})
    # if len(data)==2:
    #      df_NH = pd.DataFrame(data[0]).sort_index()
    #      df_SH = pd.DataFrame(data[1]).sort_index()
        
    #      print('sub-routine to plot NH data: ', df_NH)
    #      print('sub-routine to plot SH data: ', df_SH)
    
#     # fig, axes = plt.subplots(nrows=2, ncols=1, constrained_layout=True)
 
    
 
#     # for ax,df in zip(axes, [df_NH, df_SH]):
#     #     p = df.plot(ax=ax,kind='bar',
#     #                 figsize=(11, 7),
#     #                 fontsize=14,
#     #                 stacked=True,
#     #                 title=title,
#     #                 rot=45,
#     #                 legend=True,
#     #                 grid=True,
#     #                 # layout='constrained',
#     #                 xlabel=x_label,
#     #                 ylabel='Percentage of observations (%)'
#     #                 )
#     #     p.xaxis.set_major_locator(MultipleLocator(5))
#     # plt.show()

    # Create figure and subplots
#    fig, axs = plt.subplots(3, 1, figsize=(25, 5), gridspec_kw={'height_ratios': [1, 2.5, 1]})
    fig, axs = plt.subplots(3, 1, figsize=(18, 8), gridspec_kw={'height_ratios': [1.0, 3.5, 1.0]})
    plt.subplots_adjust(hspace=0.01)  # Adjust vertical spacing
    


# # #    mooring: RGB 0/107/164, #006BA4
# # #    Submarine: RGB 255/128/14, #FF800E  
# # #    Airborne: RGB 171/171/171, #ABABAB 
# # #    Drifting buoy: RGB 89/89/89, #595959 
# # #    Ship: RGB 95/158/209, #5F9ED1


    y_valuesCDR = list(dataCDR.keys())
    event_dataCDR = list(dataCDR.values())

    y_valuesNH = list(dataNH.keys())
    event_dataNH = list(dataNH.values())

    y_valuesSH = list(dataSH.keys())
    event_dataSH = list(dataSH.values())


    # Define colors for each event data
#    MyColors = ['red', 'green', 'blue']
    MyColorsCDR = ['#007200', '#007200', '#00b300', '#00b300']  # Example hex colors
    MyColorsNH = ['#5F9ED1', '#FF800E', '#595959', '#595959', '#595959', '#006BA4', '#006BA4', '#006BA4', '#006BA4', '#006BA4', '#ABABAB', '#ABABAB', '#ABABAB', '#ABABAB']  # Example hex colors
    MyColorsSH = ['#5F9ED1', '#595959', '#006BA4', '#ABABAB']  # Example hex colors

#d_assist, d_scicex, d_sb, d_simba, d_imb, d_npeo, d_trans, d_ulsnpi, d_bgep, d_aem, d_oib

    # Plot event data in each subplot
    axs[0].eventplot(event_dataCDR,lineoffsets=y_valuesCDR,orientation="horizontal", colors=MyColorsCDR, linelengths=0.2, linewidths=8, label='FDR4ALT')
#    axs[0].set_title('Event Plot 1')
    axs[0].set_ylim(-0.5,3.5)  # Set y-axis limits for subplot 1    plt.ylim(-0.5,9.5)    # SH

## #     ax.set_yticklabels(y_labels)
## Setting labels and title
##    plt.ylabel('NH RRDP')
    axs[0].set_ylabel('CDR')

#    # Remove x-axis labels
    axs[0].set_xticks([])
    axs[0].set_xticklabels([])
    axs[0].yaxis.set_label_coords(-0.077, .5)
##    axs[0].legend(loc="upper right")


    axs[1].eventplot(event_dataNH,lineoffsets=y_valuesNH,orientation="horizontal", colors=MyColorsNH, linelengths=0.2, linewidths=8)
##    axs[1].set_title('Event Plot 2')
    axs[1].set_ylim(-0.5,13.5)  # Set y-axis limits for subplot 1    plt.ylim(-0.5,9.5)    # SH
##    plt.ylim(-0.5,9.5)    # SH
## #     ax.set_yticklabels(y_labels)
## Setting labels and title
##    plt.ylabel('NH RRDP')
    axs[1].set_ylabel('NH RRDP')
    axs[1].set_xticklabels([])
#    axs[1].set_xticklabels([])
#    axs[1].yaxis.set_label_coords(-0.077, .5)

    axs[2].eventplot(event_dataSH,lineoffsets=y_valuesSH,orientation="horizontal", colors=MyColorsSH, linelengths=0.2, linewidths=8)
#    axs[2].set_title('Event Plot 3')
    axs[2].set_ylim(-0.5,3.5)  # Set y-axis limits for subplot 1    plt.ylim(-0.5,9.5)    # SH
    axs[2].set_ylabel('SH RRDP', ha='center')
    axs[2].yaxis.set_label_coords(-0.077, .5)


##    for key in data.keys():
##    plt.eventplot(event_data,lineoffsets=y_values,orientation="horizontal", colors=MyColorsNH, linelengths=0.2, linewidths=8)
##    plt.eventplot(event_data,lineoffsets=y_values,orientation="horizontal", colors=MyColorsSH, linelengths=0.2, linewidths=8)
##    plt.ylim(-0.5,9.5)    # NH
##    plt.ylim(-0.5,3.5)    # SH
    # Set global x-axis limits
    global_xlim = (1993, 2025)  # Set the limits as per your requirement
    for ax in axs:
        ax.set_xlim(global_xlim)

#    plt.xlim(1993,2021)
# #     ax.set_yticklabels(y_labels)
# Setting labels and title

# Set global x-ticks
    global_xticks = [1995, 2000, 2005, 2010, 2015,2020]  # Define your global x-ticks
    for ax in axs:
        ax.set_xticks(global_xticks)
    
    plt.xlabel('Year')
#    plt.ylabel('NH RRDP')
#    plt.title('Event Plot with Different Colors')

    # Show plot
    # Enable global grid
    for ax in axs:
        ax.grid(True)

    plt.tight_layout()
    plt.show()
    plt.savefig('Eventsplot_Figure_3_rev1_final_withLabel.png', bbox_inches='tight', dpi=300, transparent=True)
    
    
 
    


def generate_monthly_dates(start_year, end_year):
    dates = []
    for year in range(start_year, end_year + 1):
        for month in range(1, 13):
            dates.append(dt.date(year, month, 1))
    return dates



# def generate_year_month_range(start_year, start_month, end_year, end_month):
#     from datetime import datetime, timedelta
#     from itertools import product

#     start_date = datetime(start_year, start_month)
#     end_date = datetime(end_year, end_month)

#     all_months = [(start_date + timedelta(days=day)).strftime('%Y-%m') 
#                   for day in range((end_date - start_date).days + 1)]

#     return all_months

# # Example usage:
# start_year = 2020
# start_month = 1
# end_year = 2022
# end_month = 12

# year_month_range = generate_year_month_range(start_year, start_month, end_year, end_month)
# print(year_month_range)



def year_month_converter(dates): 


    years = [d.year for d in dates]
    months = [d.month for d in dates]
    months_dec = [(x-1) / 12. for x in months]
    
    # Sum each element in years with decimal months to create decimal-years
    year_month = [sum(i) for i in zip(years,months_dec)] 

    # using list comprehension + enumerate() to remove duplicated from list
    year_final = [i for n, i in enumerate(year_month) if i not in year_month[:n]]

    return year_final




def histograms_time(data, data_HS2 = []):
    """
    Function to calculate and plot the monthly and yearly distribution of observations.

    Parameters:
    d_sit (Data): Data instance for Sea Ice Thickness (SIT) observations.
    d_sd (Data): Data instance for Snow Depth (SD) observations.
    d_sif (Data): Data instance for Freeboard (FRB) observations.
    d_sid (Data): Data instance for Sea Ice Drift (SID) observations.

    Returns:
    None
    """
    
    dict_time_series_NH = {}
    dict_time_series_SH = {}
    dict_time_series_CDR = {}
    

    dict_time_series_CDR = {'ERS-1': {}, 'ERS-2': {}, '       ENVISAT': {}, 'CS-2': {}}


    #         # Defining the CDR coverage
    ERS1_start_year = 1993
    ERS1_start_month = 1
    ERS1_end_year = 1996   
    ERS1_end_month = 6
#    ERS1_monthly_dates = generate_year_month_range(ERS1_start_year, ERS1_start_month, ERS1_end_year, ERS1_end_month)
    ERS1_monthly_dates = generate_monthly_dates(ERS1_start_year, ERS1_end_year)
#    print(ERS1_monthly_dates)
    ERS1_year_dec_final = year_month_converter(ERS1_monthly_dates)

    ERS2_start_year = 1995
    ERS2_end_year = 2003   
    ERS2_monthly_dates = generate_monthly_dates(ERS2_start_year, ERS2_end_year)
    ERS2_year_dec_final = year_month_converter(ERS2_monthly_dates)

    ENV_start_year = 2002
    ENV_end_year = 2012   
    ENV_monthly_dates = generate_monthly_dates(ENV_start_year, ENV_end_year)
    ENV_year_dec_final = year_month_converter(ENV_monthly_dates)

    CS2_start_year = 2010
    CS2_end_year = 2025   
    CS2_monthly_dates = generate_monthly_dates(CS2_start_year, CS2_end_year)
    CS2_year_dec_final = year_month_converter(CS2_monthly_dates)

    dict_time_series_CDR['ERS-1'] = ERS1_year_dec_final
    dict_time_series_CDR['ERS-2'] = ERS2_year_dec_final
    dict_time_series_CDR['       ENVISAT'] = ENV_year_dec_final
    dict_time_series_CDR['CS-2'] = CS2_year_dec_final

#    print('this is satellite CDR: ',dict_time_series_CDR)
    # for timestamp, value in zip(timestamps, values):
    #     dict_time_series_CDR[timestamp] = 
    
    
    
    
    n=0
    
    for d, HS in zip([data, data_HS2], ['NH', 'SH']):
#        d_sit, d_sd, d_sif, d_sid = d
        if HS=='NH':
            d_oib, d_aem, d_hem, d_hem_nl, d_assist, d_bgep, d_trans, d_ulsnpi, d_npeo, d_uls_nl, d_scicex, d_imb, d_simba, d_sb, d_aspect, d_ulsawi = d
        

#        summ = 0

        # create dictonaries
#            dict_time_series_NH = {str(d_oib.var): {}, str(d_aem.var): {}, str(d_sb.var): {}}
            dict_time_series_NH = {str(d_assist.var): {}, str(d_scicex.var): {}, str(d_sb.var): {}, str(d_simba.var): {}, str(d_imb.var): {}, str(d_uls_nl.var): {}, str(d_npeo.var): {}, str(d_trans.var): {}, str(d_ulsnpi.var): {}, str(d_bgep.var): {}, str(d_hem_nl.var): {}, str(d_hem.var): {}, str(d_aem.var): {},  str(d_oib.var): {}}
#d_assist, d_scicex, d_sb, s_simba, d_imb, d_npeo, d_trans, d_ulsnpi, d_bgep, d_aem, d_oib

        # For all years in RRDP generate decimal year for each month
            start_year = 1993    # The full altimetry era
            end_year = 2025      # The full altimetry era
            monthly_dates = generate_monthly_dates(start_year, end_year)
            year_dec_final = year_month_converter(monthly_dates)


     #       print('OIB dates: ', d_oib.date)        
  #          print('SIMBA dates: ', d_simba.date)        

#        for v in [d_aem, d_sb]:
        
            for v in [d_assist, d_scicex, d_sb, d_simba, d_imb, d_uls_nl, d_npeo, d_trans, d_ulsnpi, d_bgep, d_hem_nl, d_hem, d_aem, d_oib]:

    # # print(d_oib)    
    # print(d_aem.date)    
    # # print(d_sb)    

                if len(v.date) > 0:
#                print("The list is not empty: ", v.date)    

           
            # Extract dates, years, and months from the Data instances
#            dates = [dt.datetime.strptime(dd, '%Y-%m-%dT%H:%M:%S') for dd in np.concatenate((v.date))]
                    obs_dates = [dt.datetime.strptime(dd, '%Y-%m-%dT%H:%M:%S') for dd in np.concatenate((v.date))]
                    obs_year_dec_final = year_month_converter(obs_dates)
                
#                print(v.var, obs_year_dec_final)
                
                    print(HS)
                    dict_time_series_NH[v.var] = obs_year_dec_final
# #                               
#                 # Appending lists to the dictionary
#                 if HS=='NH':
#                     dict_time_series_NH[v.var] = obs_year_dec_final
# #                 elif HS=='SH':    
# #                    print(n)
#                     dict_time_series_SH[v.var] = obs_year_dec_final


        elif HS=='SH':
            d_oib, d_aem, d_hem, d_hem_nl, d_assist, d_bgep, d_trans, d_ulsnpi, d_npeo, d_uls_nl, d_scicex, d_imb, d_simba, d_sb, d_aspect, d_ulsawi = d
#            dict_time_series_SH = {str(d_aspect.var): {}, str(d_sb.var): {}, str(d_ulsawi.var): {}, str(d_oib.var): {}}

#        summ = 0

        # create dictonaries
#            dict_time_series_NH = {str(d_oib.var): {}, str(d_aem.var): {}, str(d_sb.var): {}}
            dict_time_series_SH = {str(d_aspect.var): {}, str(d_sb.var): {}, str(d_ulsawi.var): {}, str(d_oib.var): {}}

        # For all years in RRDP generate decimal year for each month
            start_year = 1993    # The full altimetry era
            end_year = 2025      # The full altimetry era
            monthly_dates = generate_monthly_dates(start_year, end_year)
            year_dec_final = year_month_converter(monthly_dates)

        
            for v in [d_aspect, d_sb, d_ulsawi, d_oib]:


                if len(v.date) > 0:
#                print("The list is not empty: ", v.date)    

           
            # Extract dates, years, and months from the Data instances
#            dates = [dt.datetime.strptime(dd, '%Y-%m-%dT%H:%M:%S') for dd in np.concatenate((v.date))]
                    obs_dates = [dt.datetime.strptime(dd, '%Y-%m-%dT%H:%M:%S') for dd in np.concatenate((v.date))]
                    obs_year_dec_final = year_month_converter(obs_dates)
                
#                print(v.var, obs_year_dec_final)
                
                    print(HS)
                    dict_time_series_SH[v.var] = obs_year_dec_final

#        n=n+1

    print('Southern Hemisphere: ',dict_time_series_NH)        
    plot_distribution_bar(dict_time_series_CDR,dict_time_series_NH, dict_time_series_SH, f'Yearly distribution of observations, {HS}', 'Years')
                
        
        


  



# Set Hemisphere
NH_data = get_data('NH')
SH_data = get_data('SH')


# # Call the histograms_time function with your
histograms_time(NH_data, SH_data)