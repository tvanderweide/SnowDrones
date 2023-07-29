#! /usr/bin/env python

"""
Update Emlid Reach CSV Survey points with PPK position accuracy

Author: Thomas Van Der Weide
"""

import pandas as pd
import glob
from geopy import Point, distance
import os


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



if __name__ == "__main__":   
    #Define CSV and POS file locations
    main_fold = "E:/SnowDrones_HDrive/2022_Data/SurveyDates/22-3-17_PilotsPeak/SDP/rinex_reach_raw_202203172156/"
    
    # Find the CSV file
    for csvFN in sorted(glob.iglob(main_fold + "Comp.csv")):
        csvFN = csvFN.replace('\\', '/')
        print("csvFN: ", csvFN)
        survey_pts_csv_fn =  csvFN
    
    
    #Load the Files
    print('Loading: %s' % survey_pts_csv_fn)
    survey_pts = pd.read_csv(survey_pts_csv_fn, index_col=0)
    out_pt = []
    print('Processing %i input points' % survey_pts.shape[0])
    for GCP, pt in survey_pts.iterrows():
        print(GCP)
        PythonPPK = Point((pt["PyLat"], pt["PyLon"]))
        EmlidPPK = Point((pt["EmLat"], pt["EmLon"]))
        print(distance.geodesic(PythonPPK, EmlidPPK).meters)
  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    