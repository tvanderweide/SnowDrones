# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 11:08:51 2020

Rename DJI Images so they are unique
/SNOWDATA/SnowDrones/
for 2021 processing

@author: Thomas Van Der Weide
"""
import os
import glob


####--------------------------------------------------------------------------------------------------------------------####
def rename(dir, pattern, titlePattern):
    for pathAndFilename in glob.iglob(os.path.join(dir, pattern)):
        title, ext = os.path.splitext(os.path.basename(pathAndFilename))
        os.rename(pathAndFilename, 
                  os.path.join(dir, titlePattern % title + ext))
        
####--------------------------------------------------------------------------------------------------------------------####
def iterFold():
    # Go through folder structure
    for foldPath in sorted(glob.iglob(locFold + '*/')):
        foldPath = foldPath.replace('\\', '/') #ex .../2021/
        
        for dirPath in sorted(glob.iglob(foldPath + '*/')):
            dirPath = dirPath.replace('\\', '/') #ex .../2021/2021-03-23/
            if os.path.isdir(dirPath):  
                for imgFolder in sorted(glob.iglob(dirPath + 'RGB/10*')):
                    imgFolder = imgFolder.replace('\\', '/')
                    print(imgFolder)
                    if os.path.isdir(imgFolder):  
                        img_folder = imgFolder.rpartition("/")[2]
                        for imgs in sorted(glob.iglob(imgFolder + '/*.DNG')):
                            imgs = imgs.replace('\\', '/')

                            checkName = imgs.rpartition("/")[2].rpartition(".")[0]
                            # If the name hasn't already been changed to include MEDIA at the end
                            if checkName[-1] == "A":
                                # print(checkName)
                                pass
                            else:
                                newName = imgs.rpartition("/")[0] + "/" + imgs.rpartition("/")[2].rpartition(".")[0] + "_" + img_folder + "." + imgs.rpartition(".")[2]
                                # print(newName)
                                os.rename(imgs, newName)
            

    return



######-------------------------------MAIN-----------------------------------------######
## Define header folder
siteLoc = "BogusRidge"
locFold  = "/SNOWDATA/SnowDrones/" + siteLoc + "/"

imgList = iterFold()












