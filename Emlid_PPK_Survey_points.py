#! /usr/bin/env python

"""
Update Emlid Reach CSV Survey points with PPK position accuracy
Scatter plot of points used in final position

Author: David Shean
Updated By: Thomas Van Der Weide
"""

import numpy as np
import pandas as pd
import glob
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from geopy import Point, distance
from math import radians, cos, sin, asin, sqrt
import os


#Hack to update solution status
def get_solution_status(Q):
    Q = np.round(Q)
    out = None
    if Q == 1.0: 
        out = 'FIX'
    elif Q == 2.0: 
        out = 'FLOAT'
    elif Q == 5.0:
        out = 'SINGLE'
    return out

# Plot all the points collected within the survey time for each GCP
def plotScatter(pt, distances, GCPID, pt_ppk_pos, dists, save_fold):
    x_min = min(pt_ppk_pos['longitude(deg)'])
    x_max = max(pt_ppk_pos['longitude(deg)'])
    y_min = min(pt_ppk_pos['latitude(deg)'])
    y_max = max(pt_ppk_pos['latitude(deg)'])
    
    # # Plot Histograms of LatLonHeight
    # fig, ax = plt.subplots(figsize=(16,8))
    # pt_ppk_pos['longitude(deg)'].hist(grid=False)
    # plt.title(GCPID + " Longitudes")
    
    # fig, ax = plt.subplots(figsize=(16,8))
    # pt_ppk_pos['latitude(deg)'].hist(grid=False)
    # plt.title(GCPID + " Latitudes")
    
    # fig, ax = plt.subplots(figsize=(16,8))
    # pt_ppk_pos['height(m)'].hist(grid=False)
    # plt.title(GCPID + " heights")
    
    #Plot the Errors from average
    pd.DataFrame(distances).plot(y='path distance', use_index=True, title = GCPID)
    if saveImgs:
        plt.savefig(save_fold + GCPID + '_Dist.png', bbox_inches='tight')
    
    # Plot points and their elevation as color map
    fig, ax = plt.subplots(figsize=(16,8))
    im = ax.scatter(pt_ppk_pos['longitude(deg)'], pt_ppk_pos['latitude(deg)'], c=pt_ppk_pos["height(m)"])
    ax.set_xlabel('Lat')
    ax.set_ylabel('Lon')
    ax.set_title(GCPID)
    ax.text((x_min + 0.02*(x_max-x_min)), (y_min + 0.125*(y_max-y_min)), pt['solution status'], ha='left')
    ax.text((x_min + 0.02*(x_max-x_min)), (y_min + 0.1*(y_max-y_min)), 'SD Lat: ' + str(round(pt['sde'],3)) + " (m)", ha='left')
    ax.text((x_min + 0.02*(x_max-x_min)), (y_min + 0.075*(y_max-y_min)), 'SD Lon: ' + str(round(pt['sdn'],3)) + " (m)", ha='left')
    ax.text((x_min + 0.02*(x_max-x_min)), (y_min + 0.05*(y_max-y_min)), 'SD Hgt: ' + str(round(pt['sdz'],3)) + " (m)", ha='left')
    plt.axis([x_min, x_max, y_min, y_max])
    fig.colorbar(im, ax=ax)
    if saveImgs:
        plt.savefig(save_fold + GCPID + '_Scatter.png', bbox_inches='tight')
    
    return


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    # Radius of earth in kilometers is 6371
    m = 6371 * c / 1000
    return m


def calc_distances(coords: pd.DataFrame,
                  col_lat='latitude(deg)',
                  col_lon='longitude(deg)',
                  point_obj=Point) -> pd.DataFrame:
    traces = len(coords)
    distances = [None] * (traces)
    for i in range(traces):
        # start = point_obj((coords.iloc[0][col_lat], coords.iloc[0][col_lon]))
        # Find the distance from the mean location
        start = point_obj((coords[col_lat].mean(), coords[col_lon].mean()))
        finish = point_obj((coords.iloc[i][col_lat], coords.iloc[i][col_lon]))
        distances[i] = {
            'start': start,
            'finish': finish,
            'path distance': distance.geodesic(start, finish).meters,
        }
    
    distVec = pd.DataFrame(distances)["path distance"]

    return distances, distVec

def call_RMS(x):
    return np.sqrt((x**2).sum()/len(x))


if __name__ == "__main__":
    # Save files?
    saveGPS = 1
    saveAgi = 1
    saveSDP = 0
    plotImgs = 1
    saveImgs = 1
    
    #Define CSV and POS file locations
    main_fold = "E:/SnowDrones_HDrive/2022_Data/SurveyDates/22-4-07_PilotsPeak_GNSS/RS2/"
    rover_fold = main_fold + "rinex_reach*/"
    
    save_fold = main_fold + "Images/"
    if saveImgs:
        if not os.path.exists(save_fold):
            os.makedirs(save_fold)
    
    # Find the CSV file
    for csvFN in sorted(glob.iglob(main_fold + "*.csv")):
        csvFN = csvFN.replace('\\', '/')
        print("csvFN: ", csvFN)
        survey_pts_csv_fn =  csvFN
    
    # Find the POS file
    for posfile in sorted(glob.iglob(rover_fold + "*.pos")):
        posfile = posfile.replace('\\', '/')
        print("posFile: ", posfile)
        ppk_pos_fn = posfile
    
    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, parse_dates=[5,6], index_col=0)
    # Pos header
    header = 'Date UTC latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio'
    print('Loading: %s' % ppk_pos_fn)
    ppk_pos = pd.read_csv(ppk_pos_fn, comment='%', delim_whitespace=True, names=header.split(), parse_dates=[[0,1]])
    
    out_pt = []
    Agi_out = []
    print('Processing %i input points' % survey_pts.shape[0])
    for GCP, pt in survey_pts.iterrows():
        print(GCP)
        Agi_temp_dict = {}
        #Extract start/stop times for the point
        hours_offset = pt['Averaging start'].rpartition(" UTC-0")[2][:-3]
        start_str = pt['Averaging start'].rpartition(" UTC")[0]
        start = datetime.strptime(start_str,'%Y-%m-%d %H:%M:%S.%f') + timedelta(hours=int(hours_offset)) + timedelta(seconds=18)
        end_str = pt['Averaging end'].rpartition(" UTC")[0]
        end = datetime.strptime(end_str,'%Y-%m-%d %H:%M:%S.%f') + timedelta(hours=int(hours_offset)) + timedelta(seconds=18)
        
        #Determine indices in ppk pos file for corresponding datetime
        ppk_pos_idx1 = (ppk_pos['Date_UTC'] >= start) & (ppk_pos['Date_UTC'] < end)
        ppk_pos_idx2 = (ppk_pos['Date_UTC'] >= start) & (ppk_pos['Date_UTC'] < end) & (ppk_pos['Q'] == 1)
        
        # Use 'Fixed' points if possible
        if (np.count_nonzero(ppk_pos_idx2) >= 5):
            ppk_pos_idx = ppk_pos_idx2
            print("True")
        else:
            ppk_pos_idx = ppk_pos_idx1
            print("False")
        
        #Pull out corresponding ppk positions
        pt_ppk_pos = ppk_pos[ppk_pos_idx]
        #Should check that pt['sample count'] == pt_ppk_pos.count()[0]
        # pt_ppk_pos.reset_index(drop=True, inplace=True)
        pt_ppk_pos["index"] = range(len(pt_ppk_pos))
        # pt_ppk_posT = pt_ppk_pos.iloc[300:]
        # Plot the .pos file locations
        distances, dists = calc_distances(pt_ppk_pos)
        normalized_dists=(dists-dists.min())/(dists.max()-dists.min()) * len(pt_ppk_pos)

        #Compute statistics for pos
        pt_ppk_pos_mean = pt_ppk_pos.mean()
        pt_ppk_pos_std = pt_ppk_pos.std()

        #Update fields for point
        pt['Lon'] = pt_ppk_pos_mean['longitude(deg)']
        pt['Lat'] = pt_ppk_pos_mean['latitude(deg)']
        pt['Alt'] = pt_ppk_pos_mean['height(m)']
        Agi_temp_dict['Location'] = GCP
        Agi_temp_dict['Description'] = pt["Description"]
        Agi_temp_dict['Lon'] = pt_ppk_pos_mean['longitude(deg)']
        Agi_temp_dict['Lat'] = pt_ppk_pos_mean['latitude(deg)']
        Agi_temp_dict['Alt'] = pt_ppk_pos_mean['height(m)']

        #Assume antenna height is not included in PPK solution pos output
        pt['Alt'] -= pt['Antenna height']

        #Calculate mean solution status
        #Could add a Q = 1 filter 
        pt['solution status'] = get_solution_status(pt_ppk_pos_mean['Q'])

        #Compute standard deviation of all positions at this point 
        #Should really convert to a local coord system, but this estimate works for now
        lat_m = 111000
        lon_m = np.cos(np.radians(pt['Lat']))*lat_m
        pt['sde'] = pt_ppk_pos_std['longitude(deg)'] * lon_m
        pt['sdn'] = pt_ppk_pos_std['latitude(deg)'] * lat_m
        pt['sdz'] = pt_ppk_pos_std['height(m)']
        if pt['solution status'] == "FIX":
            mult = 2
        else:
            mult = 1.2
        Agi_temp_dict['sde'] = round(pt_ppk_pos_std['longitude(deg)'] * lon_m * mult, 4)
        Agi_temp_dict['sdn'] = round(pt_ppk_pos_std['latitude(deg)'] * lat_m * mult, 4)
        Agi_temp_dict['sdz'] = round(pt_ppk_pos_std['height(m)'] * mult, 4)

        out_pt.append(pt)
        Agi_out.append(Agi_temp_dict)
        if plotImgs:
            plotScatter(pt, distances, GCP, pt_ppk_pos, normalized_dists, save_fold)
    
    
    if saveGPS: 
        out_df = pd.DataFrame(out_pt)
        #Write out new file
        out_fn = main_fold + csvFN.rpartition("/")[2].rpartition(".")[0] + '_ppk_pos_PROCESSED.csv'
        print("Writing out: %s" % out_fn)
        out_df.to_csv(out_fn)
    
    # Write out Agisoft File
    if saveAgi:
        Agi_out_df = pd.DataFrame(Agi_out)
        Agi_out_df.drop('Description', axis=1, inplace=True)
        Agi_out_fn = main_fold + 'GCPs_Agisoft_Processed.csv'
        print("Writing out: %s" % Agi_out_fn)
        Agi_out_df.to_csv(Agi_out_fn, index=False)
        
    # Write out Agisoft File
    if saveSDP:
        SDP_out_df = pd.DataFrame(Agi_out)
        SDP_out_fn = main_fold + 'SDP_Processed.csv'
        print("Writing out: %s" % SDP_out_fn)
        SDP_out_df.to_csv(SDP_out_fn, index=False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    