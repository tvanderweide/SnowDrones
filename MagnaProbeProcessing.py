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
#import numpy as np
from datetime import timedelta
import math
import matplotlib.pyplot as plt
import seaborn as sns


def process():
    #Define CSV, POS, and output file locations
    survey_pts_csv_fn =  "P:/SnowDrones/Surveys/2024/2024-04-16_Utq/MagnaProbe/MagnaProbe.csv"
    ppk_pos_fn = "P:/SnowDrones/Surveys/2024/2024-04-16_Utq/GPSDATA/MagnaProbe/SEPT1050_Emlid_OPUS_Forward_35_GloON_15Deg.pos"
    out_csv_fn = 'P:/SnowDrones/Surveys/2024/2024-04-16_Utq/MagnaProbe/MagnaProbe_corrected.csv'
    plotTitle = '2024-04-16_Utq'
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
    #print("Writing out: %s" % out_csv_fn)
    #out_df.to_csv(out_csv_fn, index=False)
    #out_df.to_feather(out_csv_fn.split(".")[0]+".feather")
    return plotTitle, ppk_pos, survey_pts, outDF

def plot(plotTitle, outDF):
    print("Entering plot func")
    # create boxplot with a different y scale for different rows
    #bins = [0,0.25,1,max(outDF['sde(m)'])]
    labels = ["Fix","Float","Single"]
    #outDF['value_group'] = pd.cut(outDF['Lateral RMS'],bins = bins,labels=labels)
    groups = outDF.groupby('Q')
    fig, axes = plt.subplots(1, 3)
    i = 0
    for ID, group in groups:
        ax = sns.boxplot(y=3*group['sdne(m)'].abs(), ax=axes.flatten()[i])
        ax.set_xlabel(labels[i])
        ax.set_ylim(3*group['sdne(m)'].abs().min()*0.85, 3*group['sdne(m)'].abs().max()*1.15)
        if i == 0:
            ax.set_ylabel('3*SDNE or Location Uncertainty (m)')
        else:
            ax.set_ylabel('')
        ax.text(0.1, 0.942, "N= " + str(len(group)), transform=ax.transAxes, size=10, weight='bold')
        ax.text(0.1, 0.9, "Mean= " + str(3*round(group['sdne(m)'].abs().mean() ,2)) + " m", transform=ax.transAxes, size=10, weight='bold')
        i += 1
    fig.suptitle(plotTitle + ' MagnaProbe GPS Location Precision', fontweight ='bold') 
    fig.text(0.5, 0.02, 'GPS Location Solution', ha='center', va='center')
    plt.subplots_adjust(left=0.1,
                    right=0.9,
                    wspace=0.4, 
                    hspace=0.4)
    plt.show()

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
    
    
if __name__ == "__main__":
    plotTitle, ppk_pos, survey_pts, outDF = process()
    # plot(plotTitle, outDF)
    HeightPlot(plotTitle, ppk_pos, survey_pts)


