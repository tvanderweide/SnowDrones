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
#from datetime import timedelta
import math


def process():
    #Define CSV, POS, and output file locations
    survey_pts_csv_fn =  "C:/Users/Endure1/Documents/SnowDrones/MagnaProbe/CR800Series_OperatorViewCSV.csv"
    ppk_pos_fn = "C:/Users/Endure1/Documents/SnowDrones/MagnaProbe/SEPT3550.pos"
    out_csv_fn = 'C:/Users/Endure1/Documents/SnowDrones/MagnaProbe/CR800Series_OperatorViewCSV_corrected.csv'
    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, header = 1)
    header = 'Date GPST latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio'
    print('Loading: %s' % ppk_pos_fn)
    ppk_pos = pd.read_csv(ppk_pos_fn, comment='%', delim_whitespace=True, names=header.split(), parse_dates=[[0,1]])
    # Format the dataFrames
    survey_pts = survey_pts.iloc[2:]
    survey_pts['rmcutc'] = survey_pts["rmcutc"].astype(float)
    survey_pts['rmcutc'] = survey_pts["rmcutc"].apply(lambda x: '{0:.2f}'.format(x))
    survey_pts['Combined'] = survey_pts['rmcutcdate'] + survey_pts['rmcutc']
    survey_pts['DateTime'] = pd.to_datetime(survey_pts['Combined'], format='%d%m%y%H%M%S.%f')
    ppk_pos['DateTime'] = pd.to_datetime(ppk_pos['Date_GPST'], format='%Y-%m-%d %H:%M:%S.%f') #- timedelta(seconds=18)
    outDF = pd.merge(survey_pts, ppk_pos, on='DateTime', how='left')
    
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
    
    #Write out new file
    print("Writing out: %s" % out_csv_fn)
    #out_df.to_csv(out_csv_fn, index=False)
    #out_df.to_feather(out_csv_fn.split(".")[0]+".feather")
    return outDF

def plot(outDF):
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # create boxplot with a different y scale for different rows
    #bins = [0,0.25,1,max(outDF['sde(m)'])]
    labels = ["Fix","Float","Single"]
    #outDF['value_group'] = pd.cut(outDF['Lateral RMS'],bins = bins,labels=labels)
    groups = outDF.groupby('Q')
    fig, axes = plt.subplots(1, len(labels))
    i = 0
    for ID, group in groups:
        ax = sns.boxplot(y=group['sdne(m)'].abs(), ax=axes.flatten()[i])
        ax.set_xlabel(labels[i])
        ax.set_ylim(group['sdne(m)'].abs().min()*0.85, group['sdne(m)'].abs().max()*1.15)
        if i == 0:
            ax.set_ylabel('Location Uncertainty (m)')
        else:
            ax.set_ylabel('')
        ax.text(0.1, 0.942, "N= " + str(len(group)), transform=ax.transAxes, size=10, weight='bold')
        ax.text(0.1, 0.9, "Mean= " + str(group['sdne(m)'].abs().mean().round(2)) + " m", transform=ax.transAxes, size=10, weight='bold')
        i += 1
    fig.suptitle('2023-12-20_MCS MagnaProbe GPS Location Precision', fontweight ='bold') 
    fig.text(0.5, 0.02, 'GPS Location Solution', ha='center', va='center')
    plt.subplots_adjust(left=0.1,
                    right=0.9,
                    wspace=0.4, 
                    hspace=0.4)
    plt.show()

    
if __name__ == "__main__":
    outDF = process()
    plot(outDF)