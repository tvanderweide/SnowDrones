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
    Split into chunks for building point cloud
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
# import split_in_chunks_python
import csv
import time


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# def splitInChunks(doc, chunk):
#     # Call the split in chunks script
#     #######----------------------- List of Parameter options --------------------------------------------######
#     FILTERING = {"3": PhotoScan.NoFiltering,
#                  "0": PhotoScan.MildFiltering,
#                  "1": PhotoScan.ModerateFiltering,
#                  "2": PhotoScan.AggressiveFiltering}
    
#     MESH = {"Arbitrary": PhotoScan.SurfaceType.Arbitrary,
#             "Height Field": PhotoScan.SurfaceType.HeightField}
    
#     DENSE = {"Ultra":  1,
#              "High":   2,
#              "Medium": 4,
#              "Low":    8,
#              "Lowest": 16}
    
#     # Options
#     # Overlap percentage
#     overLap = 15
#     # How many chunks (grid layout)
#     gridX = 2
#     gridY = 2
#     # What processing to do
#     buildDenseCloud = 0
#     buildMesh = 0
#     mergeBack = 1
#     deleteChunks = 0
#     # Set Parameters
#     quality = DENSE["Low"]
#     mesh_mode = MESH["Height Field"]
#     filtering = FILTERING["1"]
#     split_in_chunks_python.splitChunks(doc, chunk, overLap, gridX, gridY, buildDenseCloud, buildMesh, mergeBack, deleteChunks, quality, mesh_mode, filtering)


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


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def mark_Images(chunk, GCPpath, path, gcp_crs, gcp_ref_acc):
    # Load the GCP Locations
    try:
        chunk.importReference(GCPpath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyz",create_markers=True, crs=gcp_crs)
    #, crs=PhotoScan.CoordinateSystem("EPSG::32611") 
    #, threshold = 3
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
        marker_name = sp_line[1]		#marker label
        
        flag = 0
        for i in range (len(chunk.cameras)):	
            if chunk.cameras[i].label == path:		#searching for the camera
                print("found camera")
                for j in range (len(chunk.markers)):	#searching for the marker (comparing with all the marker labels in chunk)
                    if chunk.markers[j].label == marker_name:
                        print("Found Marker")
                        #setting up marker projection of the correct photo
                        chunk.markers[j].projections[chunk.cameras[i]] =  PhotoScan.Marker.Projection(PhotoScan.Vector([x,y]), True)		
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
    
    ## Disable GCP-t and TCPs
    for marker in chunk.markers:
        # Set the Marker Accuracy
        marker.reference.accuracy = PhotoScan.Vector([gcp_ref_acc,gcp_ref_acc,gcp_ref_acc])
        
        if marker.label == 'Hhh':
            marker.reference.enabled = False
        elif marker.label == '100':
             marker.reference.enabled = False
        elif marker.label == 'Location':
             marker.reference.enabled = False
        elif marker.label == 'gcp0t':
             marker.reference.enabled = False
        elif marker.label == 'gcp0-t':
             marker.reference.enabled = False
        elif marker.label == 'gcp1t':
             marker.reference.enabled = False
        elif marker.label == 'gcp2t':
             marker.reference.enabled = False
        elif marker.label == 'gcp3t':
             marker.reference.enabled = False
        elif marker.label == 'gcp4t':
             marker.reference.enabled = False
        elif marker.label == 'tcp1':
             marker.reference.enabled = False
        elif marker.label == 'tcp2':
             marker.reference.enabled = False
        elif marker.label == 'tcp3':
             marker.reference.enabled = False
        elif marker.label == 'tcp4':
             marker.reference.enabled = False
        else:
             pass
    
    file.close()
    return
    

#####----------------------------------------------------------------------------------------------------------------------#######
def AlignPhoto(locFold, ProcessDate, typeFolder, chunk, Accuracy, Key_Limit, Tie_Limit, QualityFilter, QualityCriteria, img_crs, new_crs):   
    # Set the camera and marker projections
    if img_crs != new_crs:
        for camera in chunk.cameras:
            if camera.reference.location:
                camera.reference.location = PhotoScan.CoordinateSystem.transform(camera.reference.location, img_crs, new_crs)
    
    # Filter out poor qualtiy Images
    tempIQavg = 0
    if QualityFilter:
        i = 1
        tempIQ = 1
        if chunk.cameras[0].meta['Image/Quality'] is None:
            try:
                chunk.estimateImageQuality() # Metashape 1.5
            except:
                chunk.analyzePhotos() # Metashape 1.6
        for camera in chunk.cameras:
            if float(camera.meta["Image/Quality"]) < QualityCriteria:
                camera.enabled = False
                tempIQ = tempIQ + float(camera.meta["Image/Quality"])
                i += 1
        ## For Multispectral camera
        # for band in [band for camera in chunk.cameras for band in camera.planes]:
        #     if float(band.meta['Image/Quality']) < QualityCriteria:
        #         band.enabled = False
        
        tempIQavg = tempIQ / i
    
    # Get list of images in chunk
    originalImgList = list()
    for camera in chunk.cameras:
      if not camera.transform:
            originalImgList.append(camera)
    originalCount = len(originalImgList)
    
    # Metashape 1.5
    try:
        chunk.matchPhotos(accuracy=Accuracy, 
                          generic_preselection=True, 
                          reference_preselection=False, 
                          filter_mask=False, 
                          keypoint_limit=Key_Limit, 
                          tiepoint_limit=Tie_Limit)
    # Metashape 1.6
    except:
        chunk.matchPhotos(downscale=Accuracy, 
                      generic_preselection=True, 
                      reference_preselection=False, 
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
    
    aligned_list = list()
    for camera in chunk.cameras:
          if camera.transform:
                aligned_list.append(camera)
    align2 = len(aligned_list)
    
    # if more than 70% of images were aligned
    if align2 >= (originalCount * 0.7):
        successAlignment = True
        # enabling rolling shutter compensation
        try:
            for sensor in chunk.sensors:
                sensor.rolling_shutter = True
        except:
            print("Error applying rolling shutter correction")
        
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                              fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                              adaptive_fitting=False, tiepoint_covariance=False)
    else:
        successAlignment = False   
    
    #Remove unused variables
    realign_list = None
    aligned_list = None
    originalImgList = None
    # originalCount = None
    # successAlignment = False
    return successAlignment, originalCount, align2, tempIQavg
    

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
            
            if chunk.label == "Medium2":
                #Classification
                ClassifyGround(chunk, kwargs['Max_Angle'], kwargs['Cell_Size'])
                doc.save()
            
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
    print("reconstruction uncertainty")
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        # If too many points are selected for removal, clear the selection and up the threshold
        if nselected >= minRemainingPoints and threshold <= 50: # Max threshold of 50 is arbitrary
            fltr.resetSelection()
            threshold += 1
            continue
        elif nselected < minRemainingPoints:
            break
        elif threshold > 50:
            break
        else:
            break
        
    UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                          fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                          adaptive_fitting=False, tiepoint_covariance=False)


    return threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_PA(chunk, init_threshold=2.0):
    # This is used to reduce error based on projection accuracy
    minRemainingPercent = 50 #percentage of left points
    tie_points = chunk.point_cloud
    minRemainingPoints = int(len(tie_points.points) * (minRemainingPercent/100.0))
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ProjectionAccuracy)
    threshold = init_threshold
    print("projection accuracy")
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        # If we get to the threshold of 30 all selected points will be removed regardless
        if nselected >= minRemainingPoints and threshold <= 30:
            fltr.resetSelection()
            threshold += 0.1
            continue
        elif nselected < minRemainingPoints:
            break
        elif threshold > 30:
            break
        else:
            break
        
    UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                          fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                          adaptive_fitting=False, tiepoint_covariance=False)
    
    # # This is to tighten tie point accuracy value
    # chunk.tiepoint_accuracy = 0.1
    # chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
    #                       fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
    #                       fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
    #                       adaptive_fitting=False, tiepoint_covariance=False)

    return threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RE(chunk, init_threshold=0.3):
    print("repeojection error")
    # This is used to reduce error based on repeojection error
    minRemainingPercent = 80 #percentage of left points
    tie_points = chunk.point_cloud
    minRemainingPoints = int(len(tie_points.points) * (minRemainingPercent/100.0))
    
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReprojectionError)
    threshold = init_threshold
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected >= minRemainingPoints and threshold <= 4.0: #20% of total remaining points is arbitrary
            fltr.resetSelection()
            threshold += 0.01
            continue
        elif nselected < minRemainingPoints:
            break
        elif threshold > 4.0:
            break
        else:
            break
        
    UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=False, tiepoint_covariance=False)

    return threshold

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

#####----------------------------------------------------------------------------------------------------------------------#######
def markerError(chunk):
    listaErrores = []
    for marker in chunk.markers:
        if marker.reference.enabled:
           source = chunk.crs.unproject(marker.reference.location) #measured values in geocentric coordinates
           estim = chunk.transform.matrix.mulp(marker.position) #estimated coordinates in geocentric coordinates
           local = chunk.crs.localframe(chunk.transform.matrix.mulp(marker.position)) #local LSE coordinates
           error = local.mulv(estim - source)
           
           total = error.norm()      #error punto
           SumCuadrado = (total) ** 2    #cuadrado del error
           listaErrores += [SumCuadrado]      #lista de los errores
       
    #print(listaErrores)
    
    suma = sum(listaErrores)
    n = len(listaErrores)
    ErrorTotal = (suma / n) ** 0.5
    
    ErrTot = str(round(ErrorTotal,2))
    return ErrTot

#####------------------------------------------------------------------------------#####
#####---------------------------- MAIN --------------------------------------------#####
#####------------------------------------------------------------------------------#####
if __name__ == '__main__':
    print("Starting Program...")
    
    #####-----------------------------------Agisoft User variables---------------------------------------------------------------------#######
    # Variables for photo alignment
    # Accuracy: HighestAccuracy, HighAccuracy, MediumAccuracy, LowAccuracy, LowestAccuracy
    # Accuracy = PhotoScan.Accuracy.MediumAccuracy
    #  Metashape 1.6 Accuracy Variables
    ALIGN = {"Highest":  0,
               "High":   1,
              # "High":   8,
              "Medium": 2,
              "Low":    4,
              "Lowest": 8}
    
    # Variables for building dense cloud
    # Quality: UltraQuality, HighQuality, MediumQuality, LowQuality, LowestQuality
    # Filter: AggressiveFiltering, ModerateFiltering, MildFiltering, NoFiltering
    # Quality = PhotoScan.Quality.LowestQuality
    # FilterMode = PhotoScan.FilterMode.MildFiltering
    
    # Metashape 1.6 Quality Variables
    DENSE = {"Ultra":  1,
               "High":   2,
              # "High":   16,
              "Medium": 4,
              "Low":    8,
              "Lowest": 16}
    
    # Metashape 1.6 Filtering Variables
    FILTERING = {"None":  PhotoScan.NoFiltering,
                 "Mild":  PhotoScan.MildFiltering,
                 "Moderate": PhotoScan.ModerateFiltering,
                 "Aggressive": PhotoScan.AggressiveFiltering}
    
    
    # Variables for dense cloud ground point classification
    # Maximum distance is usually twice of image resolution
    # Which will be calculated later
    Max_Angle = 13
    Cell_Size = 10
    
    # Variable for building 3D mesh
    # Surface: Arbitrary, HeightField
    # SurfaceSource: PointCloudData, DenseCloudData, DepthMapsData
    Surface = PhotoScan.SurfaceType.Arbitrary
    SurfaceSource = PhotoScan.DataSource.DepthMapsData
    
    # Variable for building orthomosaic
    # Since 1.4.0, users can choose performing color correction (vignetting) and balance separately.
    # Blending: AverageBlending, MosaicBlending, MinBlending, MaxBlending, DisabledBlending
    # Color_correction: True, False
    # Color_balance: True, False
    BlendingMode = PhotoScan.BlendingMode.MosaicBlending
    Color_correction = False
    Color_balance = False
    
    # Set the project projection
    # PhotoScan.CoordinateSystem("EPSG::32611") #UTM11N
    # PhotoScan.CoordinateSystem("EPSG::4326")  #WGS84
    GCP_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of the GCPs
    img_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of Image geotags
    new_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Desired Output Coordinate system
    gcp_ref_acc = 0.1 #Estimated accuracy of GCP GPS coordinates in meters
    
    # Variables for image quality filter
    # QualityFilter: True, False
    # QualityCriteria: float number range from 0 to 1 (default 0.5)
    QualityFilter = True
    QualityCriteria = 0.5
    
    
    #####--------------------------------------------------Processing---------------------------------------------------------------#######
    # Define Folders to Process
    mainFold = "/SNOWDATA/SnowDrones-Processing/"
    Loc = "LDP"
    locFold = mainFold + Loc + "/"
    # Where to save the CSV File
    imgCount_outfile = locFold + "/imgCount_df.csv"
            
    AllDates = [dI for dI in sorted(os.listdir(locFold)) if os.path.isdir(os.path.join(locFold,dI))]
    AllDates = AllDates[0:7]
    # AllDates = ["02-04-2020"]
    # DataType = ["RGB", "Thermal"]
    DataType = ["RGB"]
    
    ## Parameters for Image Alignment
    # Accuracy: Highest, High, Medium, Low, Lowest
    accuracyLvl = "High"
    Accuracy = ALIGN[accuracyLvl]
    Key_Limit = 60000
    Tie_Limit = 6000
    
    ## Parameters for the Dense Cloud
    # Quality: Ultra, High, Medium, Low, Lowest
    qualityLevel = 'High'
    Quality = DENSE[qualityLevel]
    # Filtering: None, Mild, Moderate, Aggressive
    filterLevel = 'Mild'
    FilterMode = FILTERING[filterLevel]
    
    
    # Switches for Processing steps
    loadPhotos = 1
    markImages = 1 # Must have loadPhotos = 1
    processImgs = 1
    alignPhotos = 1 # Must have processImgs = 1
    reduceError = 1 # Must have alignPhotos = 1
    createDEM = 1
    classifyGround = 0
    createDEM2 = 0
    createOrtho = 1
    
    
    # Iter through all folders
    i = 0 # Used to save the header in the CSV file
    for ProcessDate in AllDates:
        start = time.time()
        # Clear dictionary items
        temp_dict = {"Date": 0,
                        "Accuracy" : 0,
                        "Key_Limit" : 0,
                        "Tie_Limit" : 0,
                        "Quality" : 0,
                        "Filter" : 0,
                        "Aligned" : 0,
                        "Total" : 0,
                        "Img_Quality_Avg" : 0,
                        "ProcessTime" : 0,
                        "RU_Thresh" : 0,
                        "PA_Thresh" : 0,
                        "RE_Thresh" : 0,
                        "Tie_ptsBefore" : 0,
                        "Marker_errBefore" : 0,
                        "Tie_ptsAfter" : 0,
                        "Marker_errAfter" : 0}
        
        print(ProcessDate)
        dateFolder =  locFold + ProcessDate + "/"
        #Where to save the metashape Project file
        saveprojName = Loc + "_" +  ProcessDate + "_High.psx"
        psxfile = dateFolder + saveprojName
        
        # Clear the Console
        # PhotoScan.app.console.clear()
        # construct the document class
        # doc = PhotoScan.Document(psxfile)
        doc = PhotoScan.app.document
        doc.save( psxfile )
        
        if loadPhotos == 1:
            # List of existing chunks for each project
            chunk_list = doc.chunks
            
            # Create a chunk for each image type
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
                    
                    photoList = []
                    for fold in sorted(glob.iglob(typeFolder + var)):
                        fold = fold.replace('\\', '/')
                        print(fold)
                        
                        ### Photo List ###
                        # photoList = []
                        getPhotoList(fold, photoList, imgExt)
                                
                                
                    # Create imgType chunk
                    if photoList:
                        chunk = doc.addChunk()
                        chunk.label = qualityLevel
                        chunk.addPhotos(photoList)
                        doc.save()
                
                        if markImages == 1:                            
                            # Read in the GCP RTK Target GPS coordinates
                            GCP_coordFN = locFold + Loc + "_Agisoft.csv"
    
                            # Read in the GCP locations on the Images
                            typeFolder =  locFold + ProcessDate + "/RGB/"
                            TargetFN = typeFolder + "GCPwithImageLocations.csv"
                            
                            filelist = [GCP_coordFN, TargetFN]
                            print(filelist)
                            if all([os.path.isfile(f) for f in filelist]):
                                mark_Images(chunk, GCP_coordFN, TargetFN, GCP_crs, gcp_ref_acc)
                            else:
                                print("One of the GCP_coord csv files does not exist.")
                        
            # Remove any empty chunks in doc
            for chunk in list(doc.chunks):
                if not len(chunk.cameras):
                    doc.remove(chunk)
            chunk_list = doc.chunks
        doc.save()

        # Process the images
        if processImgs == 1:            
            for chunk in chunk_list:
                doc.chunk = chunk
                chunk.crs = img_crs
                # Align Photo only if it is not done yet
                if alignPhotos == 1:
                    successfulAlignment, tot_Img, align_Img, ImgQual = AlignPhoto(locFold, ProcessDate, typeFolder, chunk, Accuracy, Key_Limit, Tie_Limit, QualityFilter, QualityCriteria, img_crs, new_crs)
                    print("Successful Alignment")
                    doc.save()
                    if reduceError == 1:
                        print("Enter Reduce Error Stages")
                        # Save the new marker error
                        beforeRE = markerError(chunk)
                        beforeRE_Tie = str(len(chunk.point_cloud.points))  
                        
                        RUT = ReduceError_RU(chunk)
                        PAT = ReduceError_PA(chunk)
                        RET = ReduceError_RE(chunk)
                        
                        # Save the new marker error
                        afterRE = markerError(chunk)
                        afterRE_Tie = str(len(chunk.point_cloud.points))
                        
                else:
                    successfulAlignment = True
                print("Past Alignment Stage")
                doc.save()
                # Do the rest when there's tie point
                if successfulAlignment:
                    if createDEM == 1:
                        # Define the ortho file name and save location
                        saveOrthoLoc =  locFold + ProcessDate + "/RGB/"
                        saveOrthoName = Loc + "_" + ProcessDate + "_" + str(chunk.label) + "_Ortho.tif"
                        saveOrtho = saveOrthoLoc + saveOrthoName
            
                        # if there are over 1000 RGB images:
                        #     splitInChunks(doc, chunk)

                        try:
                            StandardWorkflow(doc, chunk, saveOrtho,
                                            Quality=Quality, FilterMode=FilterMode, 
                                            Max_Angle=Max_Angle, Cell_Size=Cell_Size, 
                                            BlendingMode=BlendingMode)
                        except:
                            print("Could not finish processing " + str(chunk.label) + "dataset.")
        
        end = time.time()
        processTime = end - start
        temp_dict= {"Date": ProcessDate,
                    "Accuracy" : accuracyLvl,
                    "Key_Limit" : Key_Limit,
                    "Tie_Limit" : Tie_Limit,
                    "Quality" : qualityLevel,
                    "Filter" : filterLevel,
                    "Aligned" : align_Img,
                    "Total" : tot_Img,
                    "Img_Quality_Avg" : ImgQual,
                    "ProcessTime" : processTime,
                    "RU_Thresh" : RUT,
                    "PA_Thresh" : PAT,
                    "RE_Thresh" : RET,
                    "Tie_ptsBefore" : beforeRE_Tie,
                    "Marker_errBefore" : beforeRE,
                    "Tie_ptsAfter" : afterRE_Tie,
                    "Marker_errAfter" : afterRE}

        doc.save()
        print("Finished Processing" + psxfile)
        
    
        # Save the results to a CSV file
        with open(imgCount_outfile, 'a', newline='') as csv_file:
            fieldnames = ['Date', 'Accuracy', 'Key_Limit', 'Tie_Limit','Quality','Filter','Aligned','Total', 'Img_Quality_Avg', "ProcessTime", 'RU_Thresh', 'PA_Thresh','RE_Thresh',
                          'Tie_ptsBefore','Marker_errBefore','Tie_ptsAfter','Marker_errAfter']
            writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
            if i == 0:
                writer.writeheader()
                i += 1
            writer.writerow(temp_dict)
        
        # Clear the doc before processing next day
        doc.clear()
    
    #End ProcessDays For Loop

    # #Close the App
    # app = PhotoScan.Application()
    # app.quit()



