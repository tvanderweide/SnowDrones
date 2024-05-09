# -*- coding: utf-8 -*-
"""
Update MagnaProbe with Septentrio GNSS Survey points to match previous format

Inputs:
    Magnaprobe data as a .csv file
    
    
@author: Thomas Van Der Weide
5/3/2024
"""

#! /usr/bin/env python

import pandas as pd
#import numpy as np
from datetime import timedelta
import math
import matplotlib.pyplot as plt
import seaborn as sns

# Function to extract latitude, longitude, and elevation
def extract_gga_info(gga_string):
    # Split the GGA string by commas
    parts = gga_string.split(',')

    # Extract latitude, longitude, and elevation
    latitude_deg, _, latitude_min = parts[2].rpartition('.')
    latitude = float(latitude_deg + '.' + latitude_min)
    if parts[3] == 'S':
        latitude *= -1

    longitude_deg, _, longitude_min = parts[4].rpartition('.')
    longitude = float(longitude_deg + '.' + longitude_min)
    if parts[5] == 'W':
        longitude *= -1
        
    ggahdop = float(parts[8])
    elevation = float(parts[9])
    Q = int(parts[6])
    ns = int(parts[7])

    return latitude, longitude, elevation, ggahdop, Q, ns
    
def process():
    #Define CSV, POS, and output file locations
    survey_pts_csv_fn =  "P:/SnowDrones/Surveys/2024/2024-04-16_Utq/BSU_MagnaProbe/MagnaProbe.csv"
    out_csv_fn = 'P:/SnowDrones/Surveys/2024/2024-04-16_Utq/BSU_MagnaProbe/MagnaProbe_corrected.csv'
    plotTitle = '2024-04-16_Utq'
    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, header = 0)
    ## Extract the UTC time from the GPS NMEA String and format this as a GPST DateTime
    time_regex = r'(\d{6}\.\d{2})'
    survey_pts['rmcutc'] = survey_pts['RMCstring'].str.extract(time_regex) # Extract the time element using str.extract()
    survey_pts['Combined'] = survey_pts['rmcutcdate'].astype(str) + survey_pts['rmcutc'].astype(str)
    survey_pts['DateTime'] = pd.to_datetime(survey_pts['Combined'], format='%d%m%y%H%M%S.%f') - timedelta(seconds=18) # + timedelta(hours=9) 
    outDF = survey_pts
    outDF[['latitude(deg)', 'longitude(deg)', 'height(m)','ggahdop', 'Q', 'ns']] = outDF['GGAstring'].apply(lambda x: pd.Series(extract_gga_info(x)))

    
    ##Format Output to match original magnaProbe format
    #Convert DMM to DD and isolate the decimal portion
    outDF['latitude_a'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[1])
    outDF['latitude_b'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[0]).multiply(60)
    outDF['longitude_a'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[1])
    outDF['longitude_b'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[0]).multiply(60)
    outDF['latitudeDDDDD'] = outDF['latitude(deg)'].apply(lambda x: math.modf(x)[0])
    outDF['longitudeDDDDD'] = outDF['longitude(deg)'].apply(lambda x: math.modf(x)[0])
    #Match timestamps to RMCUTC
    outDF['month'] = [d.date().month for d in outDF['DateTime']]
    outDF['dayofmonth'] = [d.date().day for d in outDF['DateTime']]
    outDF['hourofday'] = [d.time().hour for d in outDF['DateTime']]
    outDF['minutes'] = [d.time().minute for d in outDF['DateTime']]
    outDF['seconds'] = [d.time().second for d in outDF['DateTime']]
    outDF['microseconds'] = [d.time().microsecond for d in outDF['DateTime']]
    #Create df to output
    out_df = outDF[['TIMESTAMP','RECORD','Counter', 'DepthCm', 'BattVolts','latitude_a','latitude_b','longitude_a','longitude_b', 'Q','ns', 'ggahdop', 'height(m)', 'DepthVolts','latitudeDDDDD','longitudeDDDDD', 'month', 'dayofmonth', 'hourofday', 'minutes', 'seconds', 'microseconds', 'DateTime', 'latitude(deg)', 'longitude(deg)', 'GGAstring','RMCstring']]
    out_df = out_df.rename(columns={'Q': 'fix_quality', 'ns': 'nmbr_satellites','ggahdop':'HDOP','height(m)':'altitudeB', 'GGAstring':'GGAstringOriginal','RMCstring':'RMCstringOriginal'})
    
    ##Write out new file
    print("Writing out: %s" % out_csv_fn)
    out_df.to_csv(out_csv_fn, index=False)
    #out_df.to_feather(out_csv_fn.split(".")[0]+".feather")
    return plotTitle, outDF


def HeightPlot(plotTitle, outDF):
    ## Remove extraneous points
    filtered_df = outDF.iloc[:]
    filtered_df = filtered_df.reset_index(drop=True)
    

    ## Plot the GPS height data and SnowDepth vs Time
    fig, ax1 = plt.subplots(figsize=(10, 6))

    # Plot height vs. time on the first y-axis
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Height (m)', color='blue')
    ax1.plot(outDF['DateTime'], outDF['height(m)'], color='blue', label='GPS Height')
    ax1.tick_params(axis='y', labelcolor='blue')
    ax1.legend(loc='upper left')
    # Plot Depth vs. time on the second y-axis
    ax2 = ax1.twinx()
    ax2.set_ylabel('Depth (cm)', color='red')
    ax2.plot(filtered_df['DateTime'], filtered_df['DepthCm'], color='red', label='MagnaProbe Depth (cm)')
    ax2.tick_params(axis='y', labelcolor='red')
    ax2.legend(loc='upper right')
    
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
    
    
if __name__ == "__main__":
    plotTitle, outDF = process()
    # plot(plotTitle, outDF)
    HeightPlot(plotTitle, outDF)


