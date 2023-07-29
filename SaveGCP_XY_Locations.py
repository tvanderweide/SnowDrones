# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 14:43:06 2022

@author: Endure1
"""
import Metashape as PhotoScan
import pandas as pd

##Commands to run things in console
doc = PhotoScan.app.document
chunk = doc.chunk
chunk.label

# Create the DF
df = pd.DataFrame(columns=['Project','Marker_name','Img_Name','Img_X','Img_Y'])
for camera in chunk.cameras:
    for marker in chunk.markers:
        if not marker.projections[camera]:
            continue
        else:
            x0, y0 = marker.projections[camera].coord
            new_row = {'Project':chunk.label,'Marker_name':marker.label,'Img_Name':camera.label,'Img_X':x0,'Img_Y':y0}
            df = df.append(new_row, ignore_index=True)

fileName = "/home/thomasvanderweide/scratch/SnowDrones/2022/LowQual2022/GCP_Img_XY_log.csv"
df.to_csv(fileName, sep=',',index=False)