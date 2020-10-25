#! /usr/bin/env python

"""
Update Emlid Reach Survey points with PPK position output from RTKLIB

"""

import os
import argparse
import numpy as np
import pandas as pd

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


def process():
    #Define CSV and POS file locations
    main_fold = 'C:/Users/Endure1/Documents/UAVFlights/GPSLogs/2020-10-16/Rover/'
    csv_fn = 'LDP_10-16-2020'
    survey_pts_csv_fn =  main_fold + csv_fn + '.csv'
    ppk_pos_fn = main_fold + 'raw_202010161801_RINEX-3_03/raw_202010161801.pos'

    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, parse_dates=[5,6], index_col=0)
    header = 'Date UTC latitude(deg) longitude(deg)  height(m)   Q  ns   sdn(m)   sde(m)   sdu(m)  sdne(m)  sdeu(m)  sdun(m) age(s)  ratio'
    print('Loading: %s' % ppk_pos_fn)
    ppk_pos = pd.read_csv(ppk_pos_fn, comment='%', delim_whitespace=True, names=header.split(), parse_dates=[[0,1]])

    out_pt = []
    print('Processing %i input points' % survey_pts.shape[0])
    for index, pt in survey_pts.iterrows():
        #Extract start/stop times for the point
        start = pt['collection start'].replace(tzinfo=None)
        end = pt['collection end'].replace(tzinfo=None)
        #Determine indices in ppk pos file for corresponding datetime
        ppk_pos_idx1 = (ppk_pos['Date_UTC'] >= start) & (ppk_pos['Date_UTC'] < end)
        ppk_pos_idx2 = (ppk_pos['Date_UTC'] >= start) & (ppk_pos['Date_UTC'] < end) & (ppk_pos['Q'] == 1)
        
        # Only used 'Fixed' points unless it gets rid of too many
        if (np.count_nonzero(ppk_pos_idx2) > (0.2 * np.count_nonzero(ppk_pos_idx1))):
            ppk_pos_idx = ppk_pos_idx2
            print("True")
        else:
            ppk_pos_idx = ppk_pos_idx1
            print("False")
            
        #Pull out corresponding ppk positions
        pt_ppk_pos = ppk_pos[ppk_pos_idx]
        #Should check that pt['sample count'] == pt_ppk_pos.count()[0]

        #Compute statistics for pos
        pt_ppk_pos_mean = pt_ppk_pos.mean()
        pt_ppk_pos_std = pt_ppk_pos.std()
        #pt_ppk_pos_med = pt_ppk_pos.median()

        #Update fields for point
        pt['longitude'] = pt_ppk_pos_mean['longitude(deg)']
        pt['latitude'] = pt_ppk_pos_mean['latitude(deg)']
        pt['height'] = pt_ppk_pos_mean['height(m)']

        #Assume antenna height is not included in PPK solution pos output
        pt['height'] -= pt['antenna height']

        #Calculate mean solution status
        #Could add a Q = 1 filter 
        pt['solution status'] = get_solution_status(pt_ppk_pos_mean['Q'])

        #Compute standard deviation of all positions at this point 
        #Should really convert to a local coord system, but this estimate works for now
        lat_m = 111000
        lon_m = np.cos(np.radians(pt['latitude']))*lat_m
        pt['sde_samples'] = pt_ppk_pos_std['longitude(deg)'] * lon_m
        pt['sdn_samples'] = pt_ppk_pos_std['latitude(deg)'] * lat_m
        pt['sdu_samples'] = pt_ppk_pos_std['height(m)']

        #Compute mean of std values from input positions
        pt['sde_mean'] = pt_ppk_pos_mean['sde(m)']
        pt['sdn_mean'] = pt_ppk_pos_mean['sdn(m)']
        pt['sdu_mean'] = pt_ppk_pos_mean['sdu(m)']

        out_pt.append(pt)

    out_df = pd.DataFrame(out_pt)
    #Write out new file
    out_fn = main_fold + csv_fn + '_ppk_pos.csv'
    print("Writing out: %s" % out_fn)
    out_df.to_csv(out_fn)

if __name__ == "__main__":
    process()