# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 17:58:05 2021

Save a CSV file from Agisoft with GCP coordinates on Images in X,Y pixels

@author: Thomas Van Der Weide
"""

import Metashape as PhotoScan
import pandas as pd
from os import path


#####----------------------------------------------------------------------------------------------------------------------#######
def findDir(chunk):
    current_doc = Metashape.app.document.path
    current_dir = path.dirname(current_doc)
    csvFN = current_dir.replace('SnowDrones_2021_Processing', 'SnowDrones').rpartition("/")[0] + "/" + str(chunk.label) + "/RGB/GCPwithImageLocations_new.csv"
    return csvFN
    
#####----------------------------------------------------------------------------------------------------------------------#######
def writeCameraMarkerCSV(mc_dict, csv_FN):
    df.to_csv(csv_FN, sep=',', index = False)
        

#####------------------------------------------------------------------------------#####
#####------------------------------- MAIN -----------------------------------------#####
#####------------------------------------------------------------------------------#####
if __name__ == '__main__':
    print("Starting Program...")
    
    # Document should already be open
    doc = PhotoScan.app.document
    chunk_list = doc.chunks
    
    
    for chunk in list(doc.chunks):
        camList = []
        csv_FN = findDir(chunk)
        for camera in chunk.cameras:
            if(camera.enabled):
                for marker in chunk.markers:
                    if not marker.projections[camera]:
                        continue
                    else:
                        x0, y0 = marker.projections[camera].coord
                        cameraDictionary = {"Date" : str(chunk.label), "Marker_Name" : str(marker.label), "Img_Name" : str(camera.label), "Img_X" : x0, "Img_Y" : y0}
                        camList.append(cameraDictionary)
                        
        df = pd.DataFrame(camList, columns=["Date", "Marker_Name","Img_Name","Img_X","Img_Y"])
        writeCameraMarkerCSV(df, csv_FN)
        
    print("Finished Program...")
        
        
        
        
        
        
        