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
# import xml.etree.ElementTree as ET


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
def img_qualityControl(chunk, QualityCriteria):
        # Tell Agisoft to calculate each images quality
    if chunk.cameras[0].meta['Image/Quality'] is None:
        try:
            chunk.estimateImageQuality() # Metashape 1.5
        except:
            chunk.analyzePhotos() # Metashape 1.6
    
    # Disable cameras that are lower than the defined QualityCritera
    print(len(chunk.cameras))
    for camera in chunk.cameras:
        if float(camera.meta["Image/Quality"]) < QualityCriteria:
            camera.enabled = False
    x = 0

    for camera in chunk.cameras:
    	if not camera.enabled:
    		x = x + 1
    		pass
    	else:
    		pass
    		
    print ('Disabled images: %d ' % x)
    
    ## For Multispectral camera
    # for band in [band for camera in chunk.cameras for band in camera.planes]:
    #     if float(band.meta['Image/Quality']) < QualityCriteria:
    #         band.enabled = False
        
    # Calculate the average remaining image quality
    images_quality = 0
    i = 0
    if len(chunk.cameras) > 0:
        for camera in chunk.cameras:
            if camera.enabled:
                images_quality = images_quality + float(camera.meta['Image/Quality'])
                i += 1
    
    if i > 0:
        images_quality_avg = images_quality / i
    else:
        images_quality_avg = 0
    
    return images_quality_avg, x


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def mark_Images(chunk, GCPpath, path, gcp_crs, gcp_ref_acc):
    # Load the GCP Locations
    try:
        chunk.importReference(GCPpath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyz",create_markers=True, crs=gcp_crs)
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
        elif marker.label == 'BASE':
            marker.reference.enabled = False
        elif marker.label == '100':
             marker.reference.enabled = False
        elif marker.label == 'gcp0center':
             marker.reference.enabled = False
        elif marker.label == 'Location':
             marker.reference.enabled = False
        # elif marker.label == 'gcp0t':
        #      marker.reference.enabled = False
        elif marker.label == 'gcp0-t':
             marker.reference.enabled = False
        # elif marker.label == 'gcp1t':
        #      marker.reference.enabled = False
        # elif marker.label == 'gcp2t':
        #      marker.reference.enabled = False
        # elif marker.label == 'gcp3t':
        #      marker.reference.enabled = False
        # elif marker.label == 'gcp4t':
        #      marker.reference.enabled = False
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
    if QualityFilter:
        # ImgQualAvg, lowQualImg = img_qualityControl(chunk, QualityFilter)
                # Tell Agisoft to calculate each images quality
        if chunk.cameras[0].meta['Image/Quality'] is None:
            try:
                chunk.estimateImageQuality() # Metashape 1.5
            except:
                chunk.analyzePhotos() # Metashape 1.6
        
        # Disable cameras that are lower than the defined QualityCritera
        print(len(chunk.cameras))
        for camera in chunk.cameras:
            if float(camera.meta["Image/Quality"]) < QualityCriteria:
                camera.enabled = False
        
        lowQualImg = 0
        for camera in chunk.cameras:
        	if not camera.enabled:
        		lowQualImg += 1
        		pass
        	else:
        		pass
        		
        print ('Disabled images: %d ' % lowQualImg)
        
        ## For Multispectral camera
        # for band in [band for camera in chunk.cameras for band in camera.planes]:
        #     if float(band.meta['Image/Quality']) < QualityCriteria:
        #         band.enabled = False
            
        # Calculate the average remaining image quality
        images_quality = 0
        i = 0
        if len(chunk.cameras) > 0:
            for camera in chunk.cameras:
                if camera.enabled:
                    images_quality = images_quality + float(camera.meta['Image/Quality'])
                    i += 1
        
        if i > 0:
            ImgQualAvg = images_quality / i
        else:
            ImgQualAvg = 0
        
            
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
    return successAlignment, originalCount, align2, ImgQualAvg, lowQualImg
    

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDenseCloud(chunk, Quality, FilterMode, savePointCloud, savePtCloudFN):
    
    print("Build Dense Cloud")
    try:
        ### Metashape 1.5
        if chunk.dense_cloud is None:
            chunk.buildDepthMaps(quality=Quality,
                         filter=FilterMode,
                         reuse_depth=False)
            chunk.buildDenseCloud(point_colors=True)
        
        # Export point cloud
        if savePointCloud == 1:
            task = PhotoScan.Tasks.ExportPoints()
            task.source_data = PhotoScan.DataSource.DenseCloudData
            task.format = PhotoScan.PointsFormat.PointsFormatLAS
            task.save_colors = True
            task.crs = PhotoScan.CoordinateSystem("EPSG::32611")
            task.path = savePtCloudFN
            task.apply(chunk)

        
    except:
        # Metashape 1.6.2
        if chunk.dense_cloud is None:
            task = PhotoScan.Tasks.BuildDepthMaps()
            task.downscale = Quality
            task.filter_mode = FilterMode
            task.reuse_depth = False
            task.subdivide_task = True
            task.apply(chunk)
            
            task = PhotoScan.Tasks.BuildDenseCloud()
            task.max_neighbors = 100
            task.subdivide_task = True
            task.point_colors = True
            task.apply(chunk)
        
        # Export point cloud
        if savePointCloud == 1:
            task = PhotoScan.Tasks.ExportPoints()
            task.source_data = PhotoScan.DataSource.DenseCloudData
            task.format = PhotoScan.PointsFormat.PointsFormatLAS
            task.save_colors = True
            task.crs = PhotoScan.CoordinateSystem("EPSG::32611")
            task.path = savePtCloudFN
            task.apply(chunk)
            

#####----------------------------------------------------------------------------------------------------------------------#######  
def ClassifyGround(chunk, Max_Angle, Cell_Size):
    DEM_resolution, Image_resolution = GetResolution(chunk)
    chunk.dense_cloud.classifyGroundPoints(max_angle=Max_Angle, 
                                           max_distance=2*Image_resolution, 
                                           cell_size=Cell_Size)

#####----------------------------------------------------------------------------------------------------------------------#######  
def BuildModel(chunk):
    try:
        chunk.buildModel(surface_type=Surface, 
                         interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                         face_count=PhotoScan.FaceCount.HighFaceCount, 
                         source_data=SurfaceSource, 
                         vertex_colors=True)
    except:
        chunk.buildModel(surface_type=Surface, 
                         interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                         face_count=PhotoScan.FaceCount.HighFaceCount, 
                         source_data=PhotoScan.DataSource.DenseCloudData, 
                         vertex_colors=True)

#####----------------------------------------------------------------------------------------------------------------------#######  
def BuildDSM(chunk):
    try:
        chunk.buildDem(source_data=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                       projection = chunk.crs)
    except:
        chunk.buildDem(source_data=PhotoScan.DataSource.DenseCloudData, 
                       interpolation=PhotoScan.Interpolation.EnabledInterpolation)
        
#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDEM(chunk):
    try:
        # If the point cloud is classified
        chunk.buildDem(source_data=PhotoScan.DataSource.DenseCloudData, 
                        interpolation=PhotoScan.Interpolation.EnabledInterpolation,
                        classes=[PhotoScan.PointClass.Ground])
        # task = PhotoScan.Tasks.BuildDem()
        # task.source_data = PhotoScan.DataSource.DenseCloudData
        # task.interpolation = PhotoScan.Interpolation.Extrapolated
        # task.projection = chunk.crs
        # task.classes = [PhotoScan.PointClass.Ground]
        # # task.network_distribute = True
        # task.apply(chunk)

    except:
        # If the point cloud is not classified
        chunk.buildDem(source_data=PhotoScan.DataSource.DenseCloudData, 
                        interpolation=PhotoScan.Interpolation.EnabledInterpolation)
        # task = PhotoScan.Tasks.BuildDem()
        # task.source_data = PhotoScan.DataSource.DenseCloudData
        # task.interpolation = PhotoScan.Interpolation.Extrapolated
        # task.projection = chunk.crs
        # # task.network_distribute = True
        # task.apply(chunk)
        
    # Export DEM
    # compression = PhotoScan.ImageCompression()
    # compression.tiff_big = True
    # chunk.exportRaster(path = saveOrthoFN, source_data=PhotoScan.ModelData, image_compression = compression)

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildMosaic(chunk, BlendingMode, saveOrthoFN, save_ortho):
    try:
        # Metashape 1.5
        chunk.buildOrthomosaic(surface_data=PhotoScan.DataSource.ElevationData, 
                                blending_mode=BlendingMode,
                                fill_holes=True, 
                                projection= chunk.crs)
        
        if save_ortho == 1:
            # Export Orthomosaic
            if os.path.exists(saveOrthoFN):
                os.remove(saveOrthoFN)
            try:
                #Metashape 1.5
                chunk.exportOrthomosaic(saveOrthoFN, image_format = PhotoScan.ImageFormatTIFF)
            except:
                #Metashape 1.6
                compression = PhotoScan.ImageCompression()
                compression.tiff_big = True
                chunk.exportRaster(path = saveOrthoFN, source_data=PhotoScan.OrthomosaicData, image_compression = compression)
    except:
        # Metashape 1.6.2 got rid of the 'projection' tag
        chunk.buildOrthomosaic(surface_data=PhotoScan.DataSource.ElevationData, 
                                blending_mode=BlendingMode,
                                fill_holes=True)
        # task = PhotoScan.Tasks.BuildOrthomosaic()
        # task.surface_data = PhotoScan.DataSource.ElevationData
        # # task.resolution = 0.05
        # task.fill_holes = True
        # task.projection = chunk.crs
        # task.blending_mode = PhotoScan.BlendingMode.MosaicBlending
        # task.subdivide_task = True
        # task.apply(chunk)
    
        if save_ortho == 1:
            # Export Orthomosaic
            if os.path.exists(saveOrthoFN):
                os.remove(saveOrthoFN)
            try:
                #Metashape 1.5
                chunk.exportOrthomosaic(saveOrthoFN, image_format = PhotoScan.ImageFormatTIFF)
            except:
                #Metashape 1.6
                compression = PhotoScan.ImageCompression()
                compression.tiff_big = True
                chunk.exportRaster(path = saveOrthoFN, source_data=PhotoScan.OrthomosaicData, image_compression = compression)
        
    return
    

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
def GetResolution(chunk):
    DEM_resolution = float(chunk.dense_cloud.meta['dense_cloud/resolution']) * chunk.transform.scale
    Image_resolution = DEM_resolution / int(chunk.dense_cloud.meta['dense_cloud/depth_downscale'])
    return DEM_resolution, Image_resolution

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RU(chunk, init_threshold=10):
    # This is used to reduce error based on reconstruction uncertainty
    # Min percentage of points to be remaining after RU
    minRemainingPercent = 50
    tie_points = chunk.point_cloud
    temp_points = tie_points.points
    minRemainingPoints = int(len(temp_points) * (minRemainingPercent/100))
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReconstructionUncertainty)
    threshold = init_threshold
    print("reconstruction uncertainty")
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in temp_points if p.selected])
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
        
    # UnselectPointMatch(chunk)
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
        
    # UnselectPointMatch(chunk)
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
    print("reprojection error")
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
    
    # UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=False, tiepoint_covariance=False)

    return threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def UnselectPointMatch(chunk, *band):
    # Unselect points that have fewer than 3 projections
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

#####----------------------------------------------------------------------------------------------------------------------#######
def secondsToText(secs):
    secs = round(secs)
    days = secs//86400
    hours = (secs - days*86400)//3600
    minutes = (secs - days*86400 - hours*3600)//60
    seconds = secs - days*86400 - hours*3600 - minutes*60
    result = ("{0} day{1} ".format(days, "s" if days!=1 else "") if days else "") + \
    ("{0} hour{1} ".format(hours, "s" if hours!=1 else "") if hours else "") + \
    ("{0} minute{1} ".format(minutes, "s" if minutes!=1 else "") if minutes else "") + \
    ("{0} second{1}".format(seconds, "s" if seconds!=1 else "") if seconds else "")
    return result


#####------------------------------------------------------------------------------#####
#####------------------------------- MAIN -----------------------------------------#####
#####------------------------------------------------------------------------------#####
if __name__ == '__main__':
    print("Starting Program...")
    
    #####-----------------------------------Agisoft User variables---------------------------------------------------------------------#######
    ## Variables for photo alignment
    ## Accuracy: HighestAccuracy, HighAccuracy, MediumAccuracy, LowAccuracy, LowestAccuracy
    ## Metashape 1.5 Accuracy Variables
    # ALIGN = {"Highest":  PhotoScan.Accuracy.HighestAccuracy,
    #         "High":   PhotoScan.Accuracy.HighAccuracy,
    #       "Medium": PhotoScan.Accuracy.MediumAccuracy,
    #       "Low":    PhotoScan.Accuracy.LowAccuracy,
    #       "Lowest": PhotoScan.Accuracy.LowestAccuracy}
    ##  Metashape 1.6 Accuracy Variables
    ALIGN = {"Highest":  0,
                "High":   1,
               # "High":   8,
              "Medium": 2,
              "Low":    4,
              "Lowest": 8}
    
    # Variables for building dense cloud
    # Quality: UltraQuality, HighQuality, MediumQuality, LowQuality, LowestQuality
    # # Metashape 1.5 Quality Variables
    # DENSE = {"Ultra":  PhotoScan.Quality.UltraQuality,
    #        "High":   PhotoScan.Quality.HighQuality,
    #       "Medium": PhotoScan.Quality.MediumQuality,
    #       "Low":    PhotoScan.Quality.LowQuality,
    #       "Lowest": PhotoScan.Quality.LowestQuality}
    
    # Metashape 1.6 Quality Variables
    DENSE = {"Ultra":  1,
               "High":   2,
              # "High":   16,
              "Medium": 4,
              "Low":    8,
              "Lowest": 16}
    
    # Filter: AggressiveFiltering, ModerateFiltering, MildFiltering, NoFiltering
    # Metashape 1.5 and 1.6 Filtering Variables
    FILTERING = {"None":  PhotoScan.NoFiltering,
                 "Mild":  PhotoScan.MildFiltering,
                 "Moderate": PhotoScan.ModerateFiltering,
                 "Aggressive": PhotoScan.AggressiveFiltering}
    
    
    # Variables for dense cloud ground point classification
    Max_Angle = 13
    Cell_Size = 10
    
    # Variable for building 3D mesh
    # Surface: Arbitrary, HeightField
    # SurfaceSource: PointCloudData, DenseCloudData, DepthMapsData
    Surface = PhotoScan.SurfaceType.Arbitrary
    SurfaceSource = PhotoScan.DataSource.DepthMapsData
    
    # Variable for building orthomosaic
    # Blending: AverageBlending, MosaicBlending, MinBlending, MaxBlending, DisabledBlending
    BlendingMode = PhotoScan.BlendingMode.MosaicBlending
    
    # Set the project projection
    # PhotoScan.CoordinateSystem("EPSG::32611") #UTM11N
    # PhotoScan.CoordinateSystem("EPSG::4326")  #WGS84
    GCP_crs = PhotoScan.CoordinateSystem("EPSG::32611") #Coordinate System of the GCPs
    img_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of Image geotags
    new_crs = PhotoScan.CoordinateSystem("EPSG::32611") #Desired Output Coordinate system
    
    #Estimated accuracy of GCP GPS coordinates in meters
    gcp_ref_acc = 0.1
    
    # Variables for image quality filter
    # QualityFilter: True, False
    # QualityCriteria: float number range from 0 to 1
    QualityFilter = True
    QualityCriteria = 0.7
        
        
    #####--------------------------------------------------Processing---------------------------------------------------------------#######
    # Define Folders to Process
    # mainFold = "/SNOWDATA/SnowDrones-Processing/"
    mainFold = "/scratch/thomasvanderweide/SnowDrones/"
    Loc = "BogusRidge"
    locFold = mainFold + Loc + "/"
    # Where to save the CSV File with processing information
    imgCount_outfile = locFold + Loc + "_stats_top.csv"
    # Identifyer for the project  ex) LDP_10-19-2020_projID
    projID = "_HighTop.psx"
    # Point Cloud file name Id ex) LDP_10-19-2020_High_<ptCloudID>
    ptCloudID = "_PtCloud_11-24_Top.las"
    # orthoMosaic file name Id ex) LDP_10-19-2020_High_<orthoID>
    orthoID = "_Ortho_11-24_Top.tif"
    # Read in the GCP RTK Target GPS coordinates
    # GCP_coordFN = locFold + Loc + "_Agisoft.csv"
    GCP_coordFN = locFold + Loc + "_GCPs_Topcon.csv"
    
    # Define the Dates to process
    AllDates = [dI for dI in sorted(os.listdir(locFold)) if os.path.isdir(os.path.join(locFold,dI))]
    AllDates = AllDates[1:3]
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
    
    
    #####------   Switches for Processing steps   ---------#####
    loadPhotos = 1
    markImages = 1 # Must have loadPhotos = 1
    
    processImgs = 1
    alignPhotos = 1 
    reduceError = 1 # Must have alignPhotos = 1
    writeCSV = 1  # Must have alignPhotos = 1, otherwise all these values will be zero
    
    # Must have processImgs = 1 for any of the below
    processPointCloud = 1
    savePointCloud = 1 # Must have processPointCloud == 1
    buildMesh = 0
    buildDSM = 0
    createDEM = 1
    buildOrthomosaic = 1
    saveOrtho = 1
    
    # classifyGround = 0
    # createDEM2 = 0
    
    # Iter through all folders
    j = 0 # Used to save the header in the CSV file
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
                        "lowQualityImg" : 0,
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
        saveprojName = Loc + "_" +  ProcessDate + projID
        psxfile = dateFolder + saveprojName
        
        # Clear the Console
        # PhotoScan.app.console.clear()
        
        # Open new Metashape document
        doc = PhotoScan.app.document
        # Try to open an existing project file
        try: 
            doc.open( psxfile, read_only=False, ignore_lock=True )
        # Save to a new project file
        except: 
            doc.save( psxfile )
            
        
        # List of existing chunks for each project
        chunk_list = doc.chunks
            
        if loadPhotos == 1:
            
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
                        
                        ### Photo List
                        getPhotoList(fold, photoList, imgExt)
                                
                    # Create imgType chunk
                    if photoList:
                        chunk = doc.addChunk()
                        chunk.label = qualityLevel
                        chunk.addPhotos(photoList)
                        doc.save()
                
                if markImages == 1:                            
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
                chunk.crs = new_crs
                # Align Photo only if it is not done yet
                if alignPhotos == 1:
                    successfulAlignment, tot_Img, align_Img, ImgQual, lowQualImg = AlignPhoto(locFold, ProcessDate, typeFolder, chunk, Accuracy, Key_Limit, Tie_Limit, QualityFilter, QualityCriteria, img_crs, new_crs)
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
                    # if there are over 1000 RGB images:
                    #     splitInChunks(doc, chunk)
                    if processPointCloud == 1:
                        savePtCloudFN = locFold + ProcessDate + "/RGB/" + Loc + "_" + ProcessDate + "_" + str(chunk.label) + ptCloudID
                        BuildDenseCloud(chunk, Quality, FilterMode, savePointCloud, savePtCloudFN)
                        doc.save()
                    
                    if buildMesh == 1:
                        #Build Mesh
                        if chunk.model is None:
                            BuildModel(chunk)
                        doc.save()
                        
                    if buildDSM == 1:
                        # #Build DSM
                        if chunk.elevation is None:
                            BuildDSM(chunk)
                            doc.save()
                        
                    if createDEM == 1:
                        # Build DEM
                        if chunk.elevation is None:
                            BuildDEM(chunk)
                            doc.save()
                    
                    if buildOrthomosaic == 1:
                        #Create OrthoMosaic
                        if chunk.orthomosaic is None:
                            orthoFN = locFold + ProcessDate + "/RGB/" + Loc + "_" + ProcessDate + "_" + str(chunk.label) + orthoID
                            BuildMosaic(chunk, BlendingMode, orthoFN, saveOrtho)
                        doc.save()
        
        end = time.time()

        doc.save()
        print("Finished Processing" + psxfile)
        
        # Save the results to a CSV file
        if writeCSV == 1:
            procTime = end - start
            processTime = secondsToText(procTime)
            temp_dict= {"Date": ProcessDate,
                        "Accuracy" : accuracyLvl,
                        "Key_Limit" : Key_Limit,
                        "Tie_Limit" : Tie_Limit,
                        "Quality" : qualityLevel,
                        "Filter" : filterLevel,
                        "Aligned" : align_Img,
                        "Total" : tot_Img,
                        "lowQualityImg" : lowQualImg,
                        "Img_Quality_Avg" : ImgQual,
                        "ProcessTime" : processTime,
                        "RU_Thresh" : RUT,
                        "PA_Thresh" : PAT,
                        "RE_Thresh" : RET,
                        "Tie_ptsBefore" : beforeRE_Tie,
                        "Marker_errBefore" : beforeRE,
                        "Tie_ptsAfter" : afterRE_Tie,
                        "Marker_errAfter" : afterRE}
        
            with open(imgCount_outfile, 'w+', newline='') as csv_file:
                fieldnames = ['Date', 'Accuracy', 'Key_Limit', 'Tie_Limit','Quality','Filter','Aligned','Total','lowQualityImg', 'Img_Quality_Avg', "ProcessTime", 'RU_Thresh', 'PA_Thresh','RE_Thresh',
                              'Tie_ptsBefore','Marker_errBefore','Tie_ptsAfter','Marker_errAfter']
                reader = csv.DictReader(csv_file)
                writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
                rows = list(reader)
                if len(rows) == 0:
                    print('csv file is empty')
                    writer.writeheader() # write the header
                
                writer.writerow(temp_dict)
        
        # Clear the doc before processing next day
        doc.clear()
    
    #End ProcessDays For Loop

    # #Close the App
    # app = PhotoScan.Application()
    # app.quit()

    
# ## Attempt to import a polygon to define the bounding box (NOT WORKING)
## https://www.agisoft.com/forum/index.php?topic=8985.0
#     kmlFile = "/scratch/thomasvanderweide/SnowDrones/LDP/LDP.kml"
#     tree = ET.parse(kmlFile)
#     root = tree.getroot()
#     nmsp = '{http://www.opengis.net/kml/2.2}'
    
#     for pm in tree.iterfind('.//{0}Placemark'.format(nmsp)):
#         print(pm.find('{0}name'.format(nmsp)).text)

#         for ls in pm.iterfind('{0}Polygon/{0}outerBoundaryIs/{0}LinearRing/{0}coordinates'.format(nmsp)):
#             strArr = ls.text.strip().split(" ")
            
#     corners = ["c1", "c2","c3","c4"]
    

#     newList = []
#     for line in strArr:
#         line = line.split(",")
#         newList.append(line)

#     cornerPts = dict(zip(corners, newList))
    
#     for marker_name in cornerPts:
#         marker = chunk.addMarker()
#         marker.label = marker_name
#         marker.reference.location = PhotoScan.Vector((float(cornerPts[marker_name][0]),float(cornerPts[marker_name][1]),0))
    
    
