# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:08:51 2020

Rename DJI Images so they are unique

@author: Thomas Van Der Weide
"""
import os
import glob
from datetime import datetime


####--------------------------------------------------------------------------------------------------------------------####
def rename(dir, pattern, titlePattern):
    for pathAndFilename in glob.iglob(os.path.join(dir, pattern)):
        title, ext = os.path.splitext(os.path.basename(pathAndFilename))
        os.rename(pathAndFilename, 
                  os.path.join(dir, titlePattern % title + ext))
        
####--------------------------------------------------------------------------------------------------------------------####
def iter_folders(mainFold, fieldSite, imgType, DateObj):
    # Find the Image folder
    input_folder = mainFold + fieldSite + "/"

    #Iterate through all the folders
    for day in sorted(glob.iglob(input_folder + '*/')):
        day = day.replace('\\', '/')
        day_folder = day.rpartition("/")[0]
        if os.path.isdir(day_folder):
            tempDate = day_folder.rpartition("/")[2]
            date = datetime.strptime(tempDate, '%m-%d-%Y').date()
            # print(date)
            if date == DateObj:
                ImgTypeFold = day + imgType + "/"
                if os.path.isdir(ImgTypeFold):  
                    for imgFolder in sorted(glob.iglob(ImgTypeFold + '10*')):
                        imgFolder = imgFolder.replace('\\', '/')
                        print(imgFolder)
                        if os.path.isdir(imgFolder):  
                            img_folder = imgFolder.rpartition("/")[2]
                            for imgs in sorted(glob.iglob(imgFolder + '/*.JPG')):
                                imgs = imgs.replace('\\', '/')

                                checkName = imgs.rpartition("/")[2].rpartition(".")[0]
                                # print(checkName)
                                # If the name hasn't already been changed to include MEDIA at the end
                                if checkName[-1] == "A":
                                    pass
                                else:
                                    newName = imgs.rpartition("/")[0] + "/" + imgs.rpartition("/")[2].rpartition(".")[0] + "_" + img_folder + "." + imgs.rpartition(".")[2]
                                    # print(newName)
                                    os.rename(imgs, newName)

    return


######-------------------------------MAIN-----------------------------------------######
fieldSiteList = ["LDP","BogusRidge","BullTrout","Headwall","PoleCat","TableRock","Treeline"]
fieldSite = fieldSiteList[1]
imgTypeList = ["RGB","Multispec","Thermal"]
imgType = imgTypeList[0]
# mainFold = "/SNOWDATA/SnowDrones_Processing_Tom/"
mainFold = "F:/SnowDrones/"
# mainFold = "/scratch/thomasvanderweide/SnowDrones/"

locFold = mainFold + fieldSite + "/"
AllDates = [dI for dI in sorted(os.listdir(locFold)) if os.path.isdir(os.path.join(locFold,dI))]
# AllDates = AllDates[8:9]


for ProcessDate in AllDates:
    datetime_obj = datetime.strptime(ProcessDate, '%m-%d-%Y').date()
    imgList = iter_folders(mainFold, fieldSite, imgType, datetime_obj)












