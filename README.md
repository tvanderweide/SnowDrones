# SnowDrones
Aerial Imagery Processing

## Summary
These scripts were created to process a time series worth of data in Agisoft Metashape. The GCP selection is done using 'imgTargetList' and 'GCP_Locations' on your own PC because doing it on Agisoft over VPN was lagging hard. The main script is 'Agisoft_Metashape_Processing' which needs to be run from Metashape using Tools -> Run Script. 

## Install
The install is used to run 'imgTargetList' and 'GCP_Locations'. If it doesn't work just install the packages to whichever enviornment you want one-by-one since there are only a handful. Once you have anaconda installed and envrionment.yml downloaded  
**Use the terminal for the following steps:**  
```
conda env create -f environment.yml
```

**Verify that it installed**
```
conda env list
```

**Activate this environment**
```
conda activate SnowDrones
```

**Install Exiftool onto system**  
(One way)  
change to download/install directory
```
git clone git://github.com/smarnach/pyexiftool.git
cd pyexif..
python setup.py install
```

## Code Overview
### RenameDJI_Images.py
This should be run first. It renames all the images taken with a DJI Mavic Pro so Agisoft can recognize the imported GCP locations correctly. Files simply have the folder they are in appened to the end of the file name.  
EX) DJI_0001_100MEDIA.JPG

### imgTargetList.py
This script uses GPS coordinates from the image metadata to save a feather database file with all images that theoretically could contain at least one GCP.  
The feather file is saved to:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Date>/\<ImgType>/imgTargets_df.feather**  
The csv containing GCP locations is expected to be at:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Site>.csv**  
User defines the survey site, image type, main folder path, and which days to process.  

### GCP_locations_qt.py
This script uses the database of images from imgTargetList.py to open images that may contain GCPs and allow the user to manually define the center of the target pixel for each image to import into Agisoft Metashape.  
Images database is expected to be at:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Date>/\<ImgType>/imgTargets_df.feather**  
CSV file to read into Agisoft:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Date>/\<ImgType>/GCPwithImageLocations.csv**  
User defines the survey site, image type, main folder path, and which days to process.  

### Agisoft_Metashape_Processing.py
This script is used to semi-automate the Agisoft processing workflow. It's completely configurable and allows for consistent processing parameters between surveys. It can also be used for sensitivity testing of parameters.  
The csv containing GCP locations is expected to be at:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Site>.csv**  
CSV file with GCP center pixel defined for each Image at:  
**/SNOWDATA/SnowDrones-Processing/\<Site>/\<Date>/\<ImgType>/GCPwithImageLocations.csv**  
User defines the survey site, image type, main folder path, days to process, and Agisoft processing parameters and output.  
#### Main Agisoft Processing Parameters   
1. Accuracy (for photo alignment)  
   1. Key Points  
   1. Tie Points  
   1. Filtering  
1. Tie-Point Quality Settings  
   1. Image Quality Filter  
      1. Default: True
      1. Removes images with less than 0.5 image quality
      1. Optimizes cameras and compensates for rolling-shutter
   1. Reduce Errors
      1. Reconstruction Uncertainty
         1. Default threshold = 10
         1. At most will remove up to 50% of the original points
         1. Threshold goes up by 1 until no more than 50% of original points are removed
      1. Projection Accuracy
         1. Default threshold = 2.0
         1. Tries to remove up to 50% of the remaining points
         1. Threshold goes up by 0.1 until no more than 50% of points are removed
         1. Max Threshold is set to 5.0 at which it will remove all points above that criteria even if it's more than 50% of the points
      1. Reprojection Error
         1. Default threshold = 0.3
         1. Threshold goes up by 0.01 until 20% or less of the remaining tie-points are removed
1. Quality (for dense cloud)



