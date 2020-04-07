# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 10:41:58 2020
Display the images with GCP's in them
Click GCP location within image
Export csv file with list of images and target locations
@author: Endure1
"""

import csv
from shapely.geometry import Point, Polygon
import exiftool
import glob
from datetime import datetime
import pyproj
import math
import pandas as pd

####--------------------------------------------------------------------------------------------------------------------####
def read_RTK(csvFN):
    # Read in the RTK (UTM) coordinates
    GCP_array = []
    with open(csvFN, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['Location'][0:3] == 'gcp':
                GCP_array.append(row)
    return GCP_array


####--------------------------------------------------------------------------------------------------------------------####
def GSD_calc(x,y,f,h):
#    x = int(metadata['EXIF:ExifImageWidth'])
#    y = int(metadata['EXIF:ExifImageHeight'])
#    f = 28 #28mm focal length camera (in relation to 35mm sensor)
#    h = 55 #height above ground [meters]
    gsd = (35 * h * 100)/(f*x)
    ax = gsd * x / 100 #[meters]
    ay = gsd * y / 100 #[meters]
    
    return ax, ay


####--------------------------------------------------------------------------------------------------------------------####
def getCorners(ax, ay, cenLat, cenLon, heading):
    #Vector math to calculate GPS coords of corners
    R = 6378.1 #Radius of the Earth
    
    cornerAng = [45,135,225,315]
    xlat = []
    ylon = []
    for ang in cornerAng:
        brng = math.radians(float(heading) + int(ang)) #Bearing is camera direction + angle to corners
        d = math.sqrt((ax/2)**2 + (ay/2)**2) / 1000 #Distance in km (function of ax/2 and ay/2)
        
        lat1 = math.radians(float(cenLat)) #Current lat point converted to radians
        lon1 = math.radians(float(cenLon)) #Current long point converted to radians
        
        lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
             math.cos(lat1)*math.sin(d/R)*math.cos(brng))
        
        lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                     math.cos(d/R)-math.sin(lat1)*math.sin(lat2))
        
        lat2 = math.degrees(lat2)
        lon2 = math.degrees(lon2)
        xlat.append(lat2)
        ylon.append(lon2)
    
    #Convert these into shapely box for testing
    coords = [list(a) for a in zip(ylon, xlat)]
    # Create a Polygon
    poly = Polygon(coords)
    
    return poly


####--------------------------------------------------------------------------------------------------------------------####
def iter_folders(mainFold, fieldSite, imgType, Date1):
    # Find the Image folder
    input_folder = mainFold + fieldSite + "/"
    
    img_list = []
    #Iterate through all the folders
    for day in sorted(glob.iglob(input_folder + '*')):
        day = day.replace('\\', '/')
        day_folder = day.rpartition("/")[2]
#        print(day)
        if day_folder[0] != "L":
            date = datetime.strptime(day_folder, '%m-%d-%Y').date()
#            print(date)
            if date == Date1:
                print(day_folder)
                ImgTypeFold = day + "/" + imgType + "/"
                for imgFolder in sorted(glob.iglob(ImgTypeFold + '*')):
#                    print(imgFolder)
                    imgFolder = imgFolder.replace('\\', '/')
                    img_folder = imgFolder.rpartition("/")[2]
                    if img_folder == "100MEDIA":
                        for imgs in sorted(glob.iglob(imgFolder + '/*')):
                            imgs = imgs.replace('\\', '/')
                            img_list.append(imgs)

                csvFN = input_folder + 'LDP_unprocessed_header.csv'
                GCP_array = read_RTK(csvFN)
                
                # Get image metadata
                metadata = get_metadata(img_list)
                # Dictionary that lists the images whose footprint covers the GCP
                imgTargetList = get_imgTargetList(GCP_array, metadata)
                            
    return imgTargetList


####--------------------------------------------------------------------------------------------------------------------####
def get_imgTargetList(GCP_array, md):
    f = int(md[0]['EXIF:FocalLengthIn35mmFormat'])
    h = 55
    x = int(md[0]['EXIF:ExifImageWidth'])
    y = int(md[0]['EXIF:ExifImageHeight'])
    
    # Test that the GCP Point is within the orthoPhoto (geom)
    imgTargetList = {}
    # Transform to project coordinates from WGS84 to UTM
    p = pyproj.Proj("+proj=utm +zone=11T +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    for GCP in GCP_array:
        if GCP['Location'][-1:] != 't':
            imgList = []
            lat_test = float(GCP['Lat'])
            lon_test = float(GCP['Lon'])
            x,y = p(lat_test, lon_test, inverse=True)
            p1 = Point(x, y)
            
            for img in md:
                cenLat = img['Composite:GPSLatitude']
                cenLon = img['Composite:GPSLongitude']
                heading = img['XMP:GimbalYawDegree']
                
                ax, ay = GSD_calc(x,y,f,h)
                poly = getCorners(ax, ay, cenLat, cenLon, heading)
                
                if p1.within(poly):
                    imgList.append(img['SourceFile'])
                    
            imgTargetList[GCP['Location']] = imgList
    
    return imgTargetList


####--------------------------------------------------------------------------------------------------------------------####
def get_metadata(imgList):
    print("Getting Metadata...")
    if type(imgList)==list:
        with exiftool.ExifTool() as et:
            metadata = et.get_metadata_batch(imgList)
    else:
        with exiftool.ExifTool() as et:
            metadata = et.get_metadata(imgList)
            
    print("Retrived Metadata...")
    return metadata


####--------------------------------------------------------------------------------------------------------------------####
def populateDF(Date1, imgTargetList, imgTargets_df):
    for key in imgTargetList:
        for entry in imgTargetList[key]:
#            print(entry)
            imgName = entry.rpartition("/")[2]
            #Append the mean to the DF
            temp_df = pd.DataFrame({"Date": Date1,
                                    "Marker_Name" : key,
                                    "Img_Name": entry,
                                    "Img_X": 0,
                                    "Img_Y": 0,}, index = [0])
            imgTargets_df = imgTargets_df.append(temp_df, ignore_index = True)
        
    return imgTargets_df


######-------------------------------MAIN-----------------------------------------######
fieldSiteList = ["LDP","BogusRidge","BullTrout","Headwall","PoleCat","TableRock","Treeline"]
fieldSite = fieldSiteList[0]
imgTypeList = ["RGB","Multispec","Thermal"]
imgType = imgTypeList[0]
mainFold = "/SNOWDATA/SnowDrones-Processing/"
Date1 = "02-04-2020"
datetime_obj = datetime.strptime(Date1, '%m-%d-%Y').date()
    
imgTargetList = iter_folders(mainFold, fieldSite, imgType, datetime_obj)

# Create the DF
imgTargets_df = pd.DataFrame({"Date": Date1,
                       "Marker_Name" : 0,
                       "Img_Name": 0,
                       "Img_X": 0,
                       "Img_Y": 0,}, index=[0])

imgTargets_df = populateDF(Date1, imgTargetList, imgTargets_df)
imgTargets_df = imgTargets_df.iloc[1:]
imgTargets_df = imgTargets_df.reset_index()
imgTargets_df = imgTargets_df.drop(['index'], axis=1)

# Save the target list as feather
imgTargets_outfile = mainFold + fieldSite + "/" + Date1 + "/" + imgType + "/imgTargets_df.feather"
print("Saving df to: " + imgTargets_outfile)
imgTargets_df.to_feather(imgTargets_outfile)


print("Finished saving feather")






