# -*- coding: utf-8 -*-
"""
Update Radar Trace CSV with PPK position output from corrected .pos file from EmlidStudio

Inputs:
    Radar Survey .csv file and .pos file

Notes:
    If it seems off it could be due to the time-offset (line 40) between UTC and GPST
    
@author: Thomas Van Der Weide
12/27/2023
"""

import pandas as pd


if __name__ == "__main__":
    
    #Define CSV, POS, and output file locations
    radar_trace_csv_fn =  "P:/SnowDrones/Surveys/2024/2023-12-21_MCS/Radar/data1.csv"
    ppk_pos_fn = "P:/SnowDrones/Surveys/2024/2023-12-21_MCS/Radar/ReachRadar_raw_20231221195150.pos"
    out_csv_fn = 'P:/SnowDrones/Surveys/2024/2023-12-21_MCS/Radar/data1_geotagged.csv'
    
    #Load the Files
    print('Loading: %s' % radar_trace_csv_fn)
    radarTrace = pd.read_csv(radar_trace_csv_fn, header = None)
    header = 'Date GPST latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio'
    print('Loading: %s' % ppk_pos_fn)
    ppk_pos = pd.read_csv(ppk_pos_fn, comment='%', delim_whitespace=True, names=header.split(), parse_dates=[[0,1]])
    
    
    ## Time sync Radar with GPS .pos file
    ## For example
    survey_pts['DateTime'] = pd.to_datetime(survey_pts['Combined'], format='%d%m%y%H%M%S.%f')
    ppk_pos['DateTime'] = pd.to_datetime(ppk_pos['Date_GPST'], format='%Y-%m-%d %H:%M:%S.%f') #- timedelta(seconds=18)
    outDF = pd.merge(survey_pts, ppk_pos, on='DateTime', how='left')
    
    ##Write out new file
    #print("Writing out: %s" % out_csv_fn)
    #out_df.to_csv(out_csv_fn, index=False)
    
    # Format the dataFrames and plot rover altitude range
    ppkpos['Heights'] = ppkpos.rolling(window=10)['height(m)'].mean()
    ppkpos['Heights'].plot(style='.')
    
    
    