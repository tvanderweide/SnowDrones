#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:21:23 2020
@author: Yu-Hsuan Tu
@author: Thomas Van Der Weide
This Python Script is developed for Agisoft MetaShape 1.5.5
Python core is 3.8.1

This script runs through all chunks and will do the following:
    1. Align Photos if there's no tie point
    2. Do the standard process if there is tie point
When aligning photos, users can decide whether using image quality to disable bad photos

Prerequisites for standard workflow:
    1. Set CRS
    2. Photo alignment
    3. Import GCP
    4. Optimse Camera
    5. Set Region
The standard workflow includes:
    Build dense point cloud
    Point cloud classification
    Build DEM
    Build orthomosaic
    
All chunks will be applied.
The DEM will be generated in duplicated chunk: "chunk name"_DEM respectively
Therefore, please avoid "_DEM" in your chunk name. Otherwise, it will not be processed.
"""
"""
Best Practice to align all images together and then split into chunks.
To-do: 
    Merge and then split into chunks
    Merge Models to include Thermal
    (https://github.com/agisoft-llc/metashape-scripts/blob/master/src/align_model_to_model.py)
    
"""
import os
import re
import glob
try:
    import Metashape as PhotoScan
except ImportError:
    import PhotoScan


#####-----------------------------------------------------------------------------------------------------------------------------------#######
## get the photo (.tif) list in specified folder
def getPhotoList(root_path, photoList, img_type):
    pattern = '.' + img_type + '$'
    for root, dirs, files in os.walk(root_path):
        for name in files:
            if re.search(pattern,name):
                cur_path = os.path.join(root, name)
                #print (cur_path)
                photoList.append(cur_path)
    return

#####----------------------------------------------------------------------------------------------------------------------#######
def AlignPhoto(chunk, Accuracy, Key_Limit, Tie_Limit, QualityFilter, QualityCriteria):
    if QualityFilter:
        if chunk.cameras[0].meta['Image/Quality'] is None:
            chunk.estimateImageQuality()
        for band in [band for camera in chunk.cameras for band in camera.planes]:
            if float(band.meta['Image/Quality']) < QualityCriteria:
                band.enabled = False
    
    originalImgList = list()
    for camera in chunk.cameras:
      if not camera.transform:
            originalImgList.append(camera)
    originalCount = len(originalImgList)
    
    chunk.matchPhotos(accuracy=Accuracy, 
                      generic_preselection=True, 
                      reference_preselection=True, 
                      filter_mask=False, 
                      keypoint_limit=Key_Limit, 
                      tiepoint_limit=Tie_Limit)
    chunk.alignCameras()
    
    # Realign non-aligned images
    realign_list = list()
    for camera in chunk.cameras:
          if not camera.transform:
                realign_list.append(camera)
    chunk.alignCameras(cameras = realign_list)
    
    realign_list = list()
    for camera in chunk.cameras:
          if not camera.transform:
                realign_list.append(camera)
    chunk.alignCameras(cameras = realign_list)
    
    realign_list = list()
    for camera in chunk.cameras:
          if not camera.transform:
                realign_list.append(camera)
    
    # if more than 40% of images were aligned
    if len(realign_list) <= originalCount * 0.4:
        successAlignment = True
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                              fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                              adaptive_fitting=False, tiepoint_covariance=False)
    else:
        successAlignment = False
    
    #Remove unused variables
    originalImgList = None
    originalCount = None
    
    return successAlignment
    

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDenseCloud(chunk, Quality, FilterMode):
    try:
        chunk.buildDenseCloud(quality=Quality, 
                              filter= FilterMode, 
                              keep_depth=False, 
                              reuse_depth=False)
    except:
        chunk.buildDepthMaps(quality=Quality,
                             filter=FilterMode,
                             reuse_depth=False)
        chunk.buildDenseCloud(point_colors=True)

#####----------------------------------------------------------------------------------------------------------------------#######  
def ClassifyGround(chunk, Max_Angle, Cell_Size):
    DEM_resolution, Image_resolution = GetResolution(chunk)
    chunk.dense_cloud.classifyGroundPoints(max_angle=Max_Angle, 
                                           max_distance=2*Image_resolution, 
                                           cell_size=Cell_Size)

#####----------------------------------------------------------------------------------------------------------------------#######  
def BuildModel(chunk):
    try:
        chunk.buildModel(surface=Surface, 
                         interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                         face_count=PhotoScan.FaceCount.HighFaceCount, 
                         source=SurfaceSource, 
                         vertex_colors=True)
    except:
        chunk.buildModel(surface=Surface, 
                         interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                         face_count=PhotoScan.FaceCount.HighFaceCount, 
                         source=PhotoScan.DataSource.DenseCloudData, 
                         vertex_colors=True)

#####----------------------------------------------------------------------------------------------------------------------#######  
def BuildDSM(chunk):
    try:
        chunk.buildDem(source=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                       projection = chunk.crs)
    except:
        chunk.buildDem(source=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation)
        
#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDEM(chunk):
    try:
        chunk.buildDem(source=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                       projection = chunk.crs,
                       classes=[PhotoScan.PointClass.Ground])
    except:
        chunk.buildDem(source=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                       projection = chunk.crs)

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildMosaic(chunk, BlendingMode, save_ortho):
    try:
        chunk.buildOrthomosaic(surface=PhotoScan.DataSource.ElevationData, 
                               blending=BlendingMode, 
                               color_correction=Color_correction, 
                               fill_holes=True, 
                               projection= chunk.crs)
    except:
        if Color_correction:
            chunk.calibrateColors(source_data=PhotoScan.DataSource.ModelData, color_balance=Color_balance)
        chunk.buildOrthomosaic(surface=PhotoScan.DataSource.ElevationData, 
                               blending=BlendingMode,  
                               fill_holes=True, 
                               projection= chunk.crs)
    doc.save()
    

#####----------------------------------------------------------------------------------------------------------------------#######
def add_altitude(chunk, flightHeightFile):
    """
    Adds user-defined altitude for camera instances in the Reference pane
    """   
    # Get the flight height
    try:
        # flightHeightFile = "/SNOWDATA/SnowDrones-Processing/LDP/01-31-2020/RGB/100MEDIA/FlightHeight.txt"
        with open(flightHeightFile , 'r') as myfile:
            data = myfile.read()
        alt = int(data)
    except:
        alt = int(55)

    # Update flight altitudes
    for camera in chunk.cameras:
        if camera.reference.location:
            coord = camera.reference.location
            camera.reference.location = PhotoScan.Vector([coord.x, coord.y, alt])

#####----------------------------------------------------------------------------------------------------------------------#######
def StandardWorkflow(doc, chunk, saveOrtho, **kwargs):
    doc.save()
    
    # Skip the chunk if it is the DEM chunk we created
    if '_DEM' in chunk.label:
        pass
    else:
        if chunk.dense_cloud is None:
            BuildDenseCloud(chunk, kwargs['Quality'], kwargs['FilterMode'])
            doc.save()
            
            # #Classification
            # ClassifyGround(chunk, kwargs['Max_Angle'], kwargs['Cell_Size'])
            # doc.save()
            
        # #Build Mesh
        # if chunk.model is None:
        #     BuildModel(chunk)
        # doc.save()
        
        # #Build DSM
        # if chunk.elevation is None:
        #     BuildDSM(chunk)
        
        #Build DEM
        if chunk.elevation is None:
            BuildDEM(chunk)
            doc.save()
        
        #Create OrthoMosaic
        if chunk.orthomosaic is None:
            BuildMosaic(chunk, kwargs['BlendingMode'], saveOrtho)
        doc.save()
        
        # Export Orthomosaic
        if os.path.exists(saveOrtho):
            os.remove(saveOrtho)
        chunk.exportOrthomosaic(saveOrtho, image_format = PhotoScan.ImageFormatTIFF)
        
    return

#####----------------------------------------------------------------------------------------------------------------------#######
def GetResolution(chunk):
    DEM_resolution = float(chunk.dense_cloud.meta['dense_cloud/resolution']) * chunk.transform.scale
    Image_resolution = DEM_resolution / int(chunk.dense_cloud.meta['dense_cloud/depth_downscale'])
    return DEM_resolution, Image_resolution

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RU(chunk, init_threshold=10):
    minRemainingPercent = 50 #percentage of left points
    points = chunk.point_cloud.points
    minRemainingPoints = int(len(points) * (minRemainingPercent/100))
    
    # This is used to reduce error based on reconstruction uncertainty
    tie_points = chunk.point_cloud
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReconstructionUncertainty)
    threshold = init_threshold
    while fltr.max_value > 10:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
#        if nselected >= len(tie_points.points) / 2 and threshold <= 50:
        if nselected >= minRemainingPoints and threshold <= 50:
            fltr.resetSelection()
            threshold += 1
            continue
        UnselectPointMatch(chunk)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected <= 5:
            break
        print('Delete {} tie point(s)'.format(nselected))
        tie_points.removeSelectedPoints()
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                              fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                              adaptive_fitting=False, tiepoint_covariance=False)
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ReconstructionUncertainty)
        threshold = init_threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_PA(chunk, init_threshold=2.0):
    # This is used to reduce error based on projection accuracy
    tie_points = chunk.point_cloud
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ProjectionAccuracy)
    threshold = init_threshold
    while fltr.max_value > 2.0:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected >= len(tie_points.points) / 2 and threshold <= 3.0:
            fltr.resetSelection()
            threshold += 0.1
            continue
        UnselectPointMatch(chunk)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected == 0:
            break
        print('Delete {} tie point(s)'.format(nselected))
        tie_points.removeSelectedPoints()
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                              fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                              adaptive_fitting=False, tiepoint_covariance=False)
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ProjectionAccuracy)
        threshold = init_threshold
    # This is to tighten tie point accuracy value
    chunk.tiepoint_accuracy = 0.1
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=False, tiepoint_covariance=False)

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RE(chunk, init_threshold=0.3):
    # This is used to reduce error based on repeojection error
    tie_points = chunk.point_cloud
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReprojectionError)
    threshold = init_threshold
    while fltr.max_value > 0.3:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected >= len(tie_points.points) / 10:
            fltr.resetSelection()
            threshold += 0.01
            continue
        UnselectPointMatch(chunk)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected == 0:
            break
        print('Delete {} tie point(s)'.format(nselected))
        tie_points.removeSelectedPoints()
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                              fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                              adaptive_fitting=False, tiepoint_covariance=False)
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ReprojectionError)
        threshold = init_threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def UnselectPointMatch(chunk, *band):
    point_cloud = chunk.point_cloud
    points = point_cloud.points
    point_proj = point_cloud.projections
    npoints = len(points)

    n_proj = dict()
    point_ids = [-1] * len(point_cloud.tracks)


    for point_id in range(0, npoints):
        point_ids[points[point_id].track_id] = point_id

    # Find the point ID using projections' track ID
    for camera in chunk.cameras:
        if camera.type != PhotoScan.Camera.Type.Regular:
            continue
        if not camera.transform:
            continue

        for proj in point_proj[camera]:
            track_id = proj.track_id
            point_id = point_ids[track_id]
            if point_id < 0:
                continue
            if not points[point_id].valid:
                continue

            if point_id in n_proj.keys():
                n_proj[point_id] += 1
            else:
                n_proj[point_id] = 1

    # Unselect points which have less than three projections
    for i in n_proj.keys():
        if n_proj[i] < 3:
            points[i].selected = False
            


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def mark_Images(chunk, GCPpath, path):
    # Load the GCP Locations
    try:
        chunk.loadReference(GCPpath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyz",create_markers=True)
    except:
        print("GCP Coordinate File Not Found.")
    
    file = open(path, "rt")	#input file
    eof = False
    line = file.readline()
    line = file.readline()
    if not len(line):
        eof = True
        
    while not eof:	
        sp_line = line.rsplit(",", 4)   #splitting read line by five parts
        print(sp_line)
        y = float(sp_line[4])			#y- coordinate of the current projection in pixels
        x = float(sp_line[3])			#x- coordinate of the current projection in pixels
        path = sp_line[2]				#camera label
        marker_name = sp_line[1]		    #marker label
        
        flag = 0
        for i in range (len(chunk.cameras)):	
            if chunk.cameras[i].label == path:		#searching for the camera
                print("found camera")
                for j in range (len(chunk.markers)):	#searching for the marker (comparing with all the marker labels in chunk)
                    if chunk.markers[j].label == marker_name:
                        print("Found Marker")
                        chunk.markers[j].projections[chunk.cameras[i]] =  PhotoScan.Marker.Projection(PhotoScan.Vector([x,y]), True)		#setting up marker projection of the correct photo)
                        flag = 1
                        break
                
                if not flag:
                    print("Not Flag")
                    marker = chunk.addMarker()
                    marker.label = marker_name
                    marker.projections[chunk.cameras[i]] =  PhotoScan.Marker.Projection(PhotoScan.Vector([x,y]), True)
                break
            
        line = file.readline()		#reading the line from input file
        if not len(line):
            eof = True
            break # EOF
        
    file.close()

    
#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Run the main code
print("Starting Program...")
#######################################################
# User variables
#
# Variables for image quality filter
# QualityFilter: True, False
# QualityCriteria: float number range from 0 to 1 (default 0.5)
QualityFilter = False
QualityCriteria = 0.5
#
# Variables for photo alignment
# Accuracy: HighestAccuracy, HighAccuracy, MediumAccuracy, LowAccuracy, LowestAccuracy
Accuracy = PhotoScan.Accuracy.LowestAccuracy
Key_Limit = 60000
Tie_Limit = 0
#
# Variables for building dense cloud
# Quality: UltraQuality, HighQuality, MediumQuality, LowQuality, LowestQuality
# Filter: AggressiveFiltering, ModerateFiltering, MildFiltering, NoFiltering
Quality = PhotoScan.Quality.LowestQuality
FilterMode = PhotoScan.FilterMode.MildFiltering
#
# Variables for dense cloud ground point classification
# Maximum distance is usually twice of image resolution
# Which will be calculated later
Max_Angle = 13
Cell_Size = 10
#
# Variable for building 3D mesh
# Surface: Arbitrary, HeightField
# SurfaceSource: PointCloudData, DenseCloudData, DepthMapsData
Surface = PhotoScan.SurfaceType.Arbitrary
SurfaceSource = PhotoScan.DataSource.DepthMapsData
#
# Variable for building orthomosaic
# Since 1.4.0, users can choose performing color correction (vignetting) and balance separately.
# Blending: AverageBlending, MosaicBlending, MinBlending, MaxBlending, DisabledBlending
# Color_correction: True, False
# Color_balance: True, False
BlendingMode = PhotoScan.BlendingMode.MosaicBlending
Color_correction = False
Color_balance = False
#
#######################################################

if __name__ == '__main__':
    # Define Folders to Process
    mainFold = "/SNOWDATA/SnowDrones-Processing/"
    Loc = "LDP"
    locFold = mainFold + Loc + "/"
    
    # Read in the GCP RTK Target GPS coordinates
    GCP_coordFN = locFold + Loc + "_unprocessed.csv"
            
    # AllDates = [dI for dI in sorted(os.listdir(locFold)) if os.path.isdir(os.path.join(locFold,dI))]
    # AllDates = AllDates[4:5]
    AllDates = ["02-14-2020"]
    DataType = ["RGB", "Thermal"]
    
    # Iter through all folders
    for ProcessDate in AllDates:
        print(ProcessDate)
        dateFolder =  locFold + ProcessDate + "/"
        #Where to save the metashape Project file
        saveprojName = Loc + "_" +  ProcessDate +"_LowestTest.psx"
        psxfile = dateFolder + saveprojName
        
        # Clear the Console
        PhotoScan.app.console.clear()
        # construct the document class
        doc = PhotoScan.app.document
        ## Open existing project or save project
        openDoc = 0
        try:
            doc.open( psxfile, read_only=False, ignore_lock=True )
            openDoc = 1
        except:
            doc.save( psxfile )
        
        for imgType in DataType:
            typeFolder =  locFold + ProcessDate + "/" + imgType + "/"
            if os.path.isdir(typeFolder):               
                #Define the subfolder names
                if imgType == "RGB":
                    var = '10*/'
                    imgExt = ".JPG"
                elif imgType == "Thermal":
                    var = '20*/'
                    imgExt = ".tiff"
                else:
                    print("imgType is non-valid for folder iteration")
                chunk = doc.addChunk()
                chunkName = imgType
                
                # #If the document was opened and if the chunk already exists then skip the creation
                # chunk_list = doc.chunks
                # if all(openDoc == 1 and (chunkName in chunkList)):
                #     print("TRUE")
                # else:
                #     print("FALSE")
                
                chunk.label = chunkName
                for fold in sorted(glob.iglob(typeFolder + var)):
                    fold = fold.replace('\\', '/')
                    print(fold)
                    
                    ### Photo List ###
                    photoList = []
                    getPhotoList(fold, photoList, imgExt)
                    # Add photos to list
                    if photoList:
                        chunk.addPhotos(photoList)
                        doc.save()
                    else:
                        print("Photo list is empty.")
            
                 
        # Remove any empty chunks
        for chunk in list(doc.chunks):
            if not len(chunk.cameras):
                doc.remove(chunk)
                
        chunk_list = doc.chunks
        print(chunk_list)
        # Loop for all initial chunks
        for chunk in chunk_list:
            doc.chunk = chunk
            
            if chunk.label == "RGB":
                print("Processing RGB")
            
                # Align Photo only if it is not done yet
                if chunk.point_cloud is None:
                    for i in range (len(chunk.cameras)):	
                        print(chunk.cameras[i].label)
                
                    # Change Image Projection to UTM11
                    new_crs = PhotoScan.CoordinateSystem("EPSG::32611") #UTM11N
                    # wgs_84 = PhotoScan.CoordinateSystem("EPSG::4326")
                    for camera in chunk.cameras:
                        camera.reference.location = new_crs.project(chunk.crs.unproject(camera.reference.location))
                    
                    
                    # Read in the GCP locations on the Images
                    typeFolder =  locFold + ProcessDate + "/" + str(chunk.label) + "/"
                    TargetFN = typeFolder + "GCPwithImageLocations.csv"
                    
                    filelist = [GCP_coordFN, TargetFN]
                    if all([os.path.isfile(f) for f in filelist]):
                        mark_Images(chunk, GCP_coordFN, TargetFN)
                    else:
                        print("GCP_coordFN does not exist.")
                    
    #                    add_altitude(chunk, flightHeightFile)
                        
                        
                    # Set different quality parameters for alignment of images for Thermal and RGB
                    if chunk.label == "Thermal":
                        Accuracy = PhotoScan.Accuracy.HighAccuracy
                    else:
                        Accuracy = PhotoScan.Accuracy.LowestAccuracy
                    
                    # successfulAlignment = AlignPhoto(chunk, Accuracy, Key_Limit, Tie_Limit, QualityFilter, QualityCriteria)
                    successfulAlignment = False
                    # ReduceError_RU(chunk)
                    # ReduceError_PA(chunk)
                # If there already is a point cloud
                else:
                    successfulAlignment = True
                    
                # Do the rest when there's tie point
                if successfulAlignment:
                    # Define the ortho file name and save location
                    saveOrthoLoc =  locFold + ProcessDate + "/" + str(chunk.label) + "/"
                    saveOrthoName = Loc + "_" +  ProcessDate + "_" + str(chunk.label) +"_LowestOrtho.tif"
                    saveOrtho = saveOrthoLoc + saveOrthoName
        #            print(saveOrtho)
        
                    # ReduceError_RE(chunk)
                    
                    # Set different quality parameters for Thermal and RGB
                    if chunk.label == "Thermal":
                        # For building dense cloud
                        # Quality: UltraQuality, HighQuality, MediumQuality, LowQuality, LowestQuality
                        Quality = PhotoScan.Quality.HighQuality
                    else:
                        Quality = PhotoScan.Quality.LowestQuality
                    
                    # try:
                    #     StandardWorkflow(doc, chunk, saveOrtho,
                    #                    Quality=Quality, FilterMode=FilterMode, 
                    #                    Max_Angle=Max_Angle, Cell_Size=Cell_Size, 
                    #                    BlendingMode=BlendingMode)
                    # except:
                    #     print("Could not finish processing " + str(chunk.label) + "dataset.")
                
        print("Finished Processing" + psxfile)
        # doc.clear()
        
    # #Close the App
    # app = PhotoScan.Application()
    # app.quit()



# Split into chunks if needed for highest res
#    https://www.agisoft.com/forum/index.php?topic=6587.0
