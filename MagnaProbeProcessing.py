# -*- coding: utf-8 -*-
"""
Update MagnaProbe Survey points with PPK position output from corrected .pos file from EmlidStudio

Inputs:
    Survey .csv file and .pos file
    
Notes:
    If it seems off it could be due to the time-offset (line 41) between UTC and GPST
    
@author: Thomas Van Der Weide
12/27/2023
"""

#! /usr/bin/env python
import pandas as pd
import numpy as np
from datetime import timedelta
import math
import matplotlib.pyplot as plt
import seaborn as sns
from geopy.distance import geodesic
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
import tkinter as tk
from tkinter import filedialog



def process(writeFlag, plotTitle, survey_pts_csv_fn, ppk_pos_fn):
    # Check if a file was selected
    if survey_pts_csv_fn:
        print("Using selected CSV File.")
    else:
        # Manually define CSV, POS, and output file locations
        survey_pts_csv_fn =  "P:/SnowDrones/Surveys/2024/2024-04-14_Utq/MagnaProbe/PlumberPoint_Shore2Tent.csv"
        print("Using hardcoded CSV File.")
    if ppk_pos_fn:
        print("Using selected .POS File.")
    else:
        ppk_pos_fn = "P:/SnowDrones/Surveys/2024/2024-04-14_Utq/GPSDATA/MagnaProbe/SEPT1050_Emlid_OPUS_Forward_35_GloON_15Deg.pos"
        print("Using hardcoded .POS File.")
    
    plotTitle = '2024-04-14_Utq'
    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, header = 0)
    header = 'Date GPST latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio'
    print('Loading: %s' % ppk_pos_fn)
    ppk_pos = pd.read_csv(ppk_pos_fn, comment='%', delim_whitespace=True, names=header.split(), parse_dates=[[0,1]])
    ## Extract the UTC time from the GPS NMEA String and format this as a GPST DateTime
    time_regex = r'(\d{6}\.\d{2})'
    survey_pts['rmcutc'] = survey_pts['RMCstring'].str.extract(time_regex) # Extract the time element using str.extract()
    survey_pts['Combined'] = survey_pts['rmcutcdate'].astype(str) + survey_pts['rmcutc'].astype(str)
    survey_pts['DateTime'] = pd.to_datetime(survey_pts['Combined'], format='%d%m%y%H%M%S.%f') - timedelta(seconds=18) # + timedelta(hours=9) 
    ## Format the POS file
    ppk_pos['DateTime'] = pd.to_datetime(ppk_pos['Date_GPST'], format='%Y-%m-%d %H:%M:%S.%f').dt.round('100ms')
    outDF = pd.merge(survey_pts, ppk_pos, on='DateTime', how='inner')
    
    ##Format Output to match original magnaProbe format
    #Convert DD to DDM and isolate the decimal portion
    outDF['latitude_a'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[1])
    outDF['latitude_b'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[0]).multiply(60)
    outDF['longitude_a'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[1])
    outDF['longitude_b'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[0]).multiply(60)
    outDF['latitudeDDDDD'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[0])
    outDF['longitudeDDDDD'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[0])
    #Isolate the timestamps
    outDF['month'] = [d.date().month for d in outDF['DateTime']]
    outDF['dayofmonth'] = [d.date().day for d in outDF['DateTime']]
    outDF['hourofday'] = [d.time().hour for d in outDF['DateTime']]
    outDF['minutes'] = [d.time().minute for d in outDF['DateTime']]
    outDF['seconds'] = [d.time().second for d in outDF['DateTime']]
    outDF['microseconds'] = [d.time().microsecond for d in outDF['DateTime']]
    #Create df to output
    out_df = outDF[['TIMESTAMP','RECORD','Counter', 'DepthCm', 'BattVolts','latitude_a','latitude_b','longitude_a','longitude_b', 'Q','ns', 'ggahdop', 'height(m)', 'DepthVolts','latitudeDDDDD','longitudeDDDDD', 'month', 'dayofmonth', 'hourofday', 'minutes', 'seconds', 'microseconds', 'latitude(deg)', 'longitude(deg)', 'sdn(m)','sde(m)','sdu(m)','sdne(m)', 'GGAstring','RMCstring']]
    out_df = out_df.rename(columns={'Q': 'fix_quality', 'ns': 'nmbr_satellites','ggahdop':'HDOP','height(m)':'altitudeB', 'GGAstring':'GGAstringOriginal','RMCstring':'RMCstringOriginal'})
    
    ##Write out new file
    if writeFlag:
        out_csv_fn = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if out_csv_fn:
            print("Writing out selected file: %s" % out_csv_fn)
        else:
            # Hardcode the output location
            out_csv_fn = 'P:/SnowDrones/Surveys/2024/2024-04-14_Utq/MagnaProbe/PlumberPoint_Shore2Tent_corrected2.csv'
            print("Writing out hardcoded file: %s" % out_csv_fn)
        out_df.to_csv(out_csv_fn, index=False)
        # out_df.to_feather(out_csv_fn.split(".")[0]+".feather")
    return plotTitle, ppk_pos, survey_pts, outDF


def HeightPlot(plotTitle, ppk_pos, survey_pts):
    ppk_pos["smoothed_height"] = ppk_pos['height(m)'].rolling(window=5, min_periods=1).mean()
    ## Remove extraneous points
    # filtered_df = survey_pts.iloc[:]
    filtered_df = survey_pts.iloc[16:-12]
    # filtered_df = survey_pts.iloc[100:155]
    filtered_df = filtered_df.reset_index(drop=True)
    #filtered_df = filtered_df[(filtered_df['DepthCm'] >= 0.03) & (filtered_df['DepthCm'] <= 118)]
    
    ## Use this to find where to start the ppk_pos file
    matching_rows = ppk_pos.loc[ppk_pos['DateTime'] == filtered_df['DateTime'][0]]
    matching_End_rows = ppk_pos.loc[ppk_pos['DateTime'] == filtered_df['DateTime'].iloc[-1]]
    first_match_index = matching_rows.index[0]
    last_match_index = matching_End_rows.index[0]
    ppk_pos = ppk_pos.iloc[first_match_index:last_match_index]
    
    matching_times = filtered_df["DateTime"].iloc[2:]
    matching_heights = ppk_pos.loc[ppk_pos['DateTime'].isin(matching_times), ['DateTime', 'height(m)']]

    ## Plot the GPS height data and SnowDepth vs Time
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot height vs. time on the first y-axis
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Height (m)', color='blue')
    ax1.plot(ppk_pos['DateTime'], ppk_pos['height(m)'], color='blue', label='GPS Height')
    ax1.fill_between(ppk_pos['DateTime'], ppk_pos['height(m)'] - 3*ppk_pos['sdu(m)'], ppk_pos['height(m)'] + 3*ppk_pos['sdu(m)'], color='lightblue', alpha=0.5, label='Uncertainty')
    ax1.scatter(matching_heights['DateTime'], matching_heights['height(m)'], color='red', marker='x', label='Recorded Depth', s=72)  # Adjust the size as needed
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.legend(loc='upper left')
    # Plot Depth vs. time on the second y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Depth (cm)', color='red')
    ax2.plot(filtered_df['DateTime'], filtered_df['DepthCm'], color='red', label='MagnaProbe Depth (cm)')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.legend(loc='upper right')
    
    # Plot shaded regions based on ppk_pos['Q']
    for i in range(len(ppk_pos)-1):
        if ppk_pos['Q'].iloc[i] == 5:
            ax1.axvspan(ppk_pos['DateTime'].iloc[i], ppk_pos['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='red', alpha=0.1)
        elif ppk_pos['Q'].iloc[i] == 2:
            ax1.axvspan(ppk_pos['DateTime'].iloc[i], ppk_pos['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='yellow', alpha=0.1)
        elif ppk_pos['Q'].iloc[i] == 1:
            ax1.axvspan(ppk_pos['DateTime'].iloc[i], ppk_pos['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='green', alpha=0.1)
    
    # Plot Details
    # Add a text box with the number of rows in 'filtered_df'
    num_rows = len(filtered_df)
    text = f'Number of Depth Measurements Shown: {num_rows}'
    ax1.text(0.51, 0.044, text, transform=ax1.transAxes, fontsize=12, verticalalignment='top')
    plt.title(plotTitle + ' MagnaProbe Depth and GPS Height')
    plt.xlabel('Time')
    plt.ylabel('Snow Depth (cm)')
    plt.grid(True)
    plt.show()
    
    
def heightComp(plotTitle, outDF):
    outDF['HeightDiff(m)'] = outDF['ggaaltitude'] - outDF['height(m)']
    # Set up a grid that has row heights in proportion to desired plot sizes
    fig = plt.figure(figsize=(10, 12))
    gs = gridspec.GridSpec(5, 3, figure=fig)
    
    # Create the first subplot on the top 2/3 of the figure
    ax1 = fig.add_subplot(gs[0:2, :])  # Spans first two rows
    
    # Plot 'ggaaltitude'
    ax1.plot(outDF['DateTime'], outDF['ggaaltitude'], label='GNGGA Altitude', marker='o', linestyle='-')
    
    # Plot 'height(m)'
    ax1.plot(outDF['DateTime'], outDF['height(m)'], label='Processed Height (m)', marker='x', linestyle='--')
    
    # Plot shaded regions based on 'Q' values
    for i in range(len(outDF)-1):
        if outDF['Q'].iloc[i] == 5:
            ax1.axvspan(outDF['DateTime'].iloc[i], outDF['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='red', alpha=0.1)
        elif outDF['Q'].iloc[i] == 2:
            ax1.axvspan(outDF['DateTime'].iloc[i], outDF['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='yellow', alpha=0.1)
        elif outDF['Q'].iloc[i] == 1:
            ax1.axvspan(outDF['DateTime'].iloc[i], outDF['DateTime'].iloc[i+1], ymin=0, ymax=0.05, color='green', alpha=0.1)
    
    # Adding titles and labels to the first subplot
    ax1.set_title(plotTitle + ' Original vs Processed Elevation over Time')
    ax1.set_xlabel('DateTime')
    ax1.set_ylabel('Altitude/Height')
    ax1.legend(loc='upper left')
    
    # Create the second subplot on the bottom 1/3 of the figure
    ax2 = fig.add_subplot(gs[2, :])  # Third row only
    
    # Plot height difference
    ax2.plot(outDF['DateTime'], outDF['HeightDiff(m)'], color='red', label='Height Diff (m)')
    ax2.set_ylabel('Height Diff (m)', color='red')
    ax2.set_ylabel('DateTime')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_title('Difference between GNGGA and Processed Altitudes')
    ax2.legend(loc='upper right')
    # Annotate number of entries and mean value on ax2
    num_entries = len(outDF)
    mean_value = outDF['HeightDiff(m)'].mean()
    annotation_text = f'n = {num_entries},\nMean = {mean_value:.2f} m'
    ax2.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10,
                 horizontalalignment='left', verticalalignment='top', color='black')
    
    # Location Accuracy Plot - fix the subplot indices
    ax4 = fig.add_subplot(gs[3, :])
    ax4.set_title('Processed MagnaProbe GPS Elevation Precision')
    ax3 = [fig.add_subplot(gs[3:5, i]) for i in range(3)]  # Create three subplots in the last row, each taking one cell of a three-column layout
    labels = ["Fix", "Float", "Single"]
    groups = outDF.groupby('Q')
    for i, (ID, group) in enumerate(groups):
        sns.boxplot(y=3*group['sdu(m)'].abs(), ax=ax3[i])
        ax3[i].set_xlabel(labels[i])
        ax3[i].set_ylim(3*group['sdu(m)'].abs().min()*0.85, 3*group['sdu(m)'].abs().max()*1.15)
        if i == 0:
            ax3[i].set_ylabel('3*SDU or Elevation Uncertainty (m)')
        else:
            ax3[i].set_ylabel('')
        ax3[i].text(0.1, 0.95, f"N= {len(group)}", transform=ax3[i].transAxes, size=10, weight='bold')
        ax3[i].text(0.1, 0.9, f"Mean= {3*round(group['sdu(m)'].abs().mean(), 2)} m", transform=ax3[i].transAxes, size=10, weight='bold')
    
    plt.tight_layout()
    plt.show()


def xyComp(plotTitle, outDF):
    # Convert latitude and longitude to decimal degrees
    outDF['Lat'] = outDF.apply(lambda row: dms_to_dd(row['ggailatitude'], row['ggan_s_ind']), axis=1)
    outDF['Lon'] = outDF.apply(lambda row: dms_to_dd(row['ggalongitude'], row['ggae_w_ind']), axis=1)
    # Apply the geodesic function to each row to calculate 'Diff(m)'
    outDF['xyDiff(m)'] = outDF.apply(lambda row: calculate_distance_geodesic(row['Lat'], row['Lon'], row['latitude(deg)'], row['longitude(deg)']), axis=1)
    
    # Create the figure layout
    fig = plt.figure(figsize=(10, 12))
    gs = gridspec.GridSpec(5, 3, figure=fig)  # Four rows, three column grid
    # Main XY Comparison Plot
    ax1 = fig.add_subplot(gs[0:2, :])
    ax1
    ax1.plot(outDF['Lat'], outDF['Lon'], label='GNGGA Pos', marker='o', linestyle='-', color='blue')
    colors = {5: 'red', 2: 'yellow', 1: 'green'}
    scatter = ax1.scatter(outDF['latitude(deg)'], outDF['longitude(deg)'], c=outDF['Q'].map(colors), label='Processed Pos', marker='x', s=100)
    legend_elements = [
        Line2D([0], [0], marker='x', color='g', label='GPS FIX', markerfacecolor='green', markersize=10),
        Line2D([0], [0], marker='x', color='y', label='GPS FLOAT', markerfacecolor='yellow', markersize=10),
        Line2D([0], [0], marker='x', color='r', label='GPS SINGLE', markerfacecolor='red', markersize=10),
        Line2D([0], [0], marker='o', color='blue', label='GNGGA Pos', markerfacecolor='blue', markersize=10)
    ]
    ax1.legend(handles=legend_elements)
    ax1.set_title(plotTitle + ' Original vs Processed Position')
    ax1.set_xlabel('Latitude')
    ax1.set_ylabel('Longitude')

    # XY Differences Plot
    ax2 = fig.add_subplot(gs[2, :])
    ax2.plot(outDF['Lat'], outDF['xyDiff(m)'], color='red', label='XYDiff (m)')
    ax2.set_ylabel('XYDiff (m)', color='red')
    ax2.set_ylabel('Latitude')
    ax2.legend(loc='upper right')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.set_title('Difference between GNGGA Pos and Processed Position')
    # Annotate number of entries and mean value on ax2
    num_entries = len(outDF)
    mean_value = outDF['xyDiff(m)'].mean()
    annotation_text = f'n = {num_entries},\nMean = {mean_value:.2f} m'
    ax2.annotate(annotation_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=10,
                 horizontalalignment='left', verticalalignment='top', color='black')

    # Location Accuracy Plot - fix the subplot indices
    ax4 = fig.add_subplot(gs[3, :])
    ax4.set_title('Processed MagnaProbe GPS Location Precision')
    ax3 = [fig.add_subplot(gs[3:5, i]) for i in range(3)]  # Create three subplots in the last row, each taking one cell of a three-column layout
    labels = ["Fix", "Float", "Single"]
    groups = outDF.groupby('Q')
    for i, (ID, group) in enumerate(groups):
        sns.boxplot(y=3*group['sdne(m)'].abs(), ax=ax3[i])
        ax3[i].set_xlabel(labels[i])
        ax3[i].set_ylim(3*group['sdne(m)'].abs().min()*0.85, 3*group['sdne(m)'].abs().max()*1.15)
        if i == 0:
            ax3[i].set_ylabel('3*SDNE or Location Uncertainty (m)')
        else:
            ax3[i].set_ylabel('')
        ax3[i].text(0.1, 0.95, f"N= {len(group)}", transform=ax3[i].transAxes, size=10, weight='bold')
        ax3[i].text(0.1, 0.9, f"Mean= {3*round(group['sdne(m)'].abs().mean(), 2)} m", transform=ax3[i].transAxes, size=10, weight='bold')
    plt.tight_layout()
    plt.show()

def dms_to_dd(dms, direction):
    # Split the DMS into degrees and minutes
    degrees = int(dms) // 100
    minutes = float(dms) % 100
    
    # Convert to decimal degrees
    decimal_degrees = degrees + (minutes / 60)
    
    # Adjust for direction
    if direction in ['S', 'W']:  # Assume S or W should be negative
        decimal_degrees = -decimal_degrees
    return decimal_degrees

# Function to calculate distance using geodesic method
def calculate_distance_geodesic(lat1, lon1, lat2, lon2):
    # Coordinates of the two points
    point1 = (lat1, lon1)
    point2 = (lat2, lon2)
    
    # Calculate distance using geodesic method
    distance = geodesic(point1, point2).meters
    return distance

    
if __name__ == "__main__":
    # Create a Tkinter root window
    root = tk.Tk()
    root.withdraw()  # Hide the root window
    # Prompt the user to select a file
    survey_pts_csv_fn = filedialog.askopenfilename(title="Select survey points CSV file", filetypes=[("CSV files", "*.csv")])
    ppk_pos_fn = filedialog.askopenfilename(title="Select rover .POS file", filetypes=[("POS files", "*.pos")])
    
    # Let user decided if they want to write the corrected magnaProbe file
    writeFlag = 0
    # Figure title ID
    plotTitle = '2024-04-14_Utq'
    
    plotTitle, ppk_pos, survey_pts, outDF = process(writeFlag, plotTitle, survey_pts_csv_fn, ppk_pos_fn)
    xyComp(plotTitle, outDF)
    heightComp(plotTitle, outDF)
    # HeightPlot(plotTitle, ppk_pos, survey_pts)

