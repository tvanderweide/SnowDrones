# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 09:45:49 2020

Correct Gopro image metadata

@author: Thomas Van Der Weide
"""

# import exiftool
import glob
import os
from exif import Image


####--------------------------------------------------------------------------------------------------------------------####
def get_metadata(imgList):
    if type(imgList)==list:
        with exiftool.ExifTool("C:/Users/Endure1/Documents/SnowDrones/exiftool-11.91/exiftool(-k).exe") as et:
            metadata = et.get_metadata_batch(imgList)
    else:
        with exiftool.ExifTool("C:/Users/Endure1/Documents/SnowDrones/exiftool-11.91/exiftool(-k).exe") as et:
            metadata = et.get_metadata(imgList)

    return metadata

####--------------------------------------------------------------------------------------------------------------------####
def get_metadata_mod(imgList):
    if type(imgList)==list:
        with exiftool.ExifTool("C:/Users/Endure1/Documents/SnowDrones/exiftool-11.91/exiftool(-k).exe") as et:
            for img in imgList:
                # et.execute("-GPSAltitude=25", img)
                et.execute("EXIF:FocalLength=4", img)
                et.execute("EXIF:FocalLengthIn35mmFormat=22", img)
                et.execute("Composite:FocalLength35efl=22.6", img)
                et.execute("Composite:FOV=82",img)
                et.execute("EXIF:GPSAltitude=25", img)
                
    else:
        with exiftool.ExifTool("C:/Users/Endure1/Documents/SnowDrones/exiftool-11.91/exiftool(-k).exe") as et:
            metadata = et.get_metadata(imgList)

    return metadata




####--------------------------------------------------------------------------------------------------------------------####
input_folder = "C:/Users/Endure1/Documents/SnowDrones/Code/RGB/geotagged_minusfour/"

FocalFolder = input_folder + "focalCorr/"
if not os.path.isdir(FocalFolder):
    os.makedirs(FocalFolder)
CorrFolder = input_folder + "/bothCorr/"
if not os.path.isdir(CorrFolder):
    os.makedirs(CorrFolder)

#Iterate through all the folders
for img in sorted(glob.iglob(input_folder + '*.JPG')):
    img = img.replace('\\', '/')
    print(img)
    imgFocal = FocalFolder + img.rpartition("/")[2]
    imgCorr = CorrFolder + img.rpartition("/")[2]

    
    # Open the Image file
    with open(img, 'rb') as image_file:
        origImg = Image(image_file)
    
    #Write image with focal length corrections
    with open(imgFocal, 'wb') as new_image_file:
        origImg.focal_length = 4
        origImg.focal_length_in_35mm_film = 22
        # my_image.gps_altitude = 25
        new_image_file.write(origImg.get_file())
      
    #Write image with both focal length corrections and altitude correction
    with open(imgCorr, 'wb') as new_image_file2:
        origImg.gps_altitude = 25
        new_image_file2.write(origImg.get_file())
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

