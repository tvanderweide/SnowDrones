# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 09:14:32 2022
Script to read in corrected Emlid Studio CSV file and output CSV formatted for Agisoft

@author: Endure1
"""
import glob
import pandas as pd

mainFold = "E:/SnowDrones_HDrive/2022_Data/SurveyDates/"

#Iterate over all corrected CSV files
for csvFile in sorted(glob.iglob(mainFold + '*/RS2/rinex*/*corrected.csv')):
    csvFile = csvFile.replace('\\', '/') #ex .../2021/
    csvPath = csvFile.rpartition("/")[0]
    csvName = csvFile.rpartition("/")[2].rpartition(".")[0]
    
    # Reformat the CSV
    csvDF = pd.read_csv(csvFile, index_col=0)
    newCSV = csvDF[['Longitude', 'Latitude','Ellipsoidal height', 'Easting RMS', 'Northing RMS', 'Elevation RMS']]
    
    # Save the new CSV
    out_fn = csvPath + "/" + csvName + '_Agisoft.csv'
    print("Writing out: %s" % out_fn)
    newCSV.to_csv(out_fn)