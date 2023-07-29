#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 3 2021
@author: Thomas Van Der Weide
This Python Script is developed for Agisoft MetaShape 1.6
Python core is 3.8.1
Purpose: Automated processing of SfM timeseries
"""

"""
This script will do the following:
    1. Import Photos
    2. Align Photos
    3. Find coded markers and then defines the bounding box region
    4. Gradual Selection and camera optimization
    5. Standard processing
When aligning photos, users can decide whether using image quality to disable bad photos


The standard processing includes:
    Build dense point cloud
    Point cloud classification
    Build DEM
    Build orthomosaic
"""

"""
To-do: 
    For large datasets it's best Practice to align all images together and then split into chunks.
        Split into chunks for building point cloud (https://github.com/agisoft-llc/metashape-scripts/blob/master/src/align_model_to_model.py)
    Parameter sensitivity iterations
        Optional format to re-run the script with different parameters for comparison
"""
import os
import re
import glob
try:
    import Metashape as PhotoScan
except ImportError:
    import PhotoScan
import csv
import time
import statistics
import math

##Commands to run things in console
# import Metashape as PhotoScan
# doc = PhotoScan.app.document
# chunk = doc.chunk
# chunk.label

#####-----------------------------------------------------------------------------------------------------------------------------------#######
## Find images of .ext type and load into project
def loadImgs(chunk, root_path, img_type):
    # If the Images have already been loaded
    if len(chunk.cameras) > 0:
        imgCount = len(chunk.cameras)

    # If there are no images in chunk
    else:  
        print("Loading Images...")
        photoList = []
        pattern = '.' + img_type + '$'
        
        for fold in sorted(glob.iglob(root_path + "Flight*/")):
            fold = fold.replace('\\', '/')
                        
            for root, dirs, files in os.walk(fold):
                for name in files:
                    if re.search(pattern,name):
                        cur_path = os.path.join(root, name)
                        photoList.append(cur_path)
                        
        # Load photos
        if photoList:
            chunk.addPhotos(photoList)
            imgCount = len(photoList)
        else:
            imgCount = "N/A"
        
    # Set the defined camera accuracy
    for camera in chunk.cameras:
        camera.reference.accuracy = PhotoScan.Vector([cam_ref_acc, cam_ref_acc, cam_ref_acc])
    
    return imgCount


#####-----------------------------------------------------------------------------------------------------------------------------------#######
## Remove Image lower than user defined threshold
def lowQualImgCount(chunk):
    lowQualImg = 0
    for camera in chunk.cameras:
    	if not camera.enabled:
    		lowQualImg += 1
    		pass
    	else:
    		pass
    return lowQualImg

def estImgQual(chunk):
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
    
    return images_quality_avg

def img_qualityControl(chunk):
    # Tell Agisoft to calculate each images quality
    if chunk.cameras[0].meta['Image/Quality'] is None:
        chunk.analyzePhotos() # Metashape 1.6
    
    # Disable cameras that are lower than the defined QualityCritera
    if QualityFilter:
        for camera in chunk.cameras:
            if float(camera.meta["Image/Quality"]) < QualityCriteria:
                camera.enabled = False
    
    lowQualImg = lowQualImgCount(chunk)
    images_quality_avg = estImgQual(chunk)

    
    return images_quality_avg, lowQualImg


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def auto_mark_Images(chunk, gcpPath):
    markerCount = len(chunk.markers)
    if markerCount == 0:
        print("Marking Images...")
        global crs   # should this be crs= chunk.crs ???
        chunk.detectMarkers(target_type=PhotoScan.CircularTarget12bit, tolerance = 25)
        chunk.importReference(path=gcpPath, format=PhotoScan.ReferenceFormatCSV, columns='nxyz',delimiter=',',
                              crs = PhotoScan.CoordinateSystem('LOCAL_CS["Local CS",LOCAL_DATUM["Local Datum",0],UNIT["metre",1]]'))
        chunk.updateTransform()
        
        # Find number of targets (markers) found
        markerCount = len(chunk.markers)
        # If there aren't enough targets, relax the tolerance
        if markerCount < 4:
            chunk.detectMarkers(target_type=PhotoScan.CircularTarget12bit, tolerance = 50)
            chunk.importReference(path=gcpPath, format=PhotoScan.ReferenceFormatCSV, columns='nxyz',delimiter=',',
                                  crs = PhotoScan.CoordinateSystem('LOCAL_CS["Local CS",LOCAL_DATUM["Local Datum",0],UNIT["metre",1]]'))
            chunk.updateTransform()
            
            # Find number of targets (markers) found
            markerCount = len(chunk.markers)
        
        # If there aren't enough targets, relax the tolerance
        if markerCount < 4:
            chunk.detectMarkers(target_type=PhotoScan.CircularTarget12bit, tolerance = 75)
            chunk.importReference(path=gcpPath, format=PhotoScan.ReferenceFormatCSV, columns='nxyz',delimiter=',',
                                  crs = PhotoScan.CoordinateSystem('LOCAL_CS["Local CS",LOCAL_DATUM["Local Datum",0],UNIT["metre",1]]'))
            chunk.updateTransform()
            
            # Find number of targets (markers) found
            markerCount = len(chunk.markers)
        
        
        # Set the bounding box around the markers
        if customBB:
            BBCameraLoc(chunk)        
    
    return markerCount


#####----------------------------------------------------------------------------------------------------------------------#######
# Define bounding box by camera locations
# https://www.agisoft.com/forum/index.php?topic=10102.0
def cross(a, b):
	result = PhotoScan.Vector([a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y *b.x])
	return result.normalized()

def BBCameraLoc(chunk):
    try:
        # Test if the photos have already been aligned
        BUFFER = 25 #percent
        new_region = PhotoScan.Region()
        xcoord = PhotoScan.Vector([10E10, -10E10])
        ycoord = PhotoScan.Vector([10E10, -10E10])
        avg = [[],[]]
        T = chunk.transform.matrix
        s = chunk.transform.matrix.scale()
        crs = chunk.crs
        z = PhotoScan.Vector([0,0])
        
        for marker in chunk.markers:
        	if marker.position:
        		coord = crs.project(T.mulp(marker.position))
        		xcoord[0] = min(coord.x, xcoord[0])
        		xcoord[1] = max(coord.x, xcoord[1])
        		ycoord[0] = min(coord.y, ycoord[0])
        		ycoord[1] = max(coord.y, ycoord[1])
        		z[0] += (coord.z)
        		z[1] += 1
        		avg[0].append(coord.x)
        		avg[1].append(coord.y)
        z = z[0] / z[1]
        avg = PhotoScan.Vector([statistics.median(avg[0]), statistics.median(avg[1]), z])
        
        corners = [PhotoScan.Vector([xcoord[0], ycoord[0], z]),
        			PhotoScan.Vector([xcoord[0], ycoord[1], z]),
        			PhotoScan.Vector([xcoord[1], ycoord[1], z]),
        			PhotoScan.Vector([xcoord[1], ycoord[0], z])]
        corners = [T.inv().mulp(crs.unproject(x)) for x in list(corners)]			
        
        side1 = corners[0] - corners[1]
        side2 = corners[0] - corners[-1]
        side1g = T.mulp(corners[0]) - T.mulp(corners[1])
        side2g = T.mulp(corners[0]) - T.mulp(corners[-1])
        side3g = T.mulp(corners[0]) - T.mulp(PhotoScan.Vector([corners[0].x, corners[0].y, 0]))
        new_size = ((100 + BUFFER) / 100) * PhotoScan.Vector([side2g.norm()/s, side1g.norm()/s, 3*side3g.norm() / s]) ##
        
        xcoord, ycoord, z = T.inv().mulp(crs.unproject(PhotoScan.Vector([sum(xcoord)/2., sum(ycoord)/2., z - 2 * side3g.z]))) #
        new_center = PhotoScan.Vector([xcoord, ycoord, z]) #by 4 corners
        
        horizontal = side2
        vertical = side1
        normal = cross(vertical, horizontal)
        horizontal = -cross(vertical, normal)
        vertical = vertical.normalized()
        
        R = PhotoScan.Matrix ([horizontal, vertical, -normal])
        new_region.rot = R.t()
        
        new_region.center = new_center
        new_region.size = new_size
        chunk.region = new_region
    except:
        print("##########----------- Could not redefine the bounding box -----------#############")
    
#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def loadGCP(chunk, gcpPath, dirPath): 
    # Remove any existing markers if markerRenew flag is == 1
    if markerRenew:
        chunk.remove(chunk.markers)
    
    if loadMarkers:
        # How many columns are in the CSV?
        with open(gcpPath, 'r') as csv:
            first_line = csv.readline()
        ncol = first_line.count(',') + 1
            
        print("Loading reference targets location...")
        # Load the GCP Locations
        if ncol == 7:
            chunk.importReference(gcpPath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyzXYZ", create_markers=True, crs=gcp_crs)
        else:
            chunk.importReference(gcpPath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyz", create_markers=True, crs=gcp_crs)
        
        
        
        ## Disable Header 'GCP' and validation targets
        for marker in chunk.markers:
            # Set the Marker Accuracy
            if ncol == 7:
                pass
            else:
                marker.reference.accuracy = PhotoScan.Vector([gcp_ref_acc,gcp_ref_acc,gcp_ref_acc])
                
            if marker.label == 'Location':
                 chunk.remove(marker)
            elif marker.label == 'Name':
                 chunk.remove(marker)
            # LDP GCP
            elif marker.label.startswith('VGCP'):
                marker.reference.enabled = False
            elif marker.label.startswith('rad'):
                marker.reference.enabled = False
        
        # Load markers and existing image XY target locations
        markerFN_csv = dirPath + markerFN
        if markerFN == "False":
            pass
        elif os.path.isfile(markerFN_csv):
            mark_Images(chunk, markerFN_csv)
        else:
            print("############ Could not find the Image Marker File. ##############")
            
    # Find number of targets (markers) loaded
    markerCount = len(chunk.markers)  
    
    return markerCount


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Load saved GCP target locations in images
def mark_Images(chunk, markerPath):
    ftype = markerPath.rpartition(".")[2].lower()
    if ftype == "csv":
        # Open the CSV with GCP locations saved as pixel coordinates on images
        file = open(markerPath, "rt")	   #input file
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
            cam_name = sp_line[2]			#camera label
            marker_name = sp_line[1]		#marker label
            
            flag = 0
            for i in range (len(chunk.cameras)):	
                if chunk.cameras[i].label == cam_name:		#searching for the camera
                    print("found camera")
                    for j in range (len(chunk.markers)):	    #searching for the marker (comparing with all the marker labels in chunk)
                        if chunk.markers[j].label.lower() == marker_name.lower():
                            print("Found Marker")
                            #setting up marker projection of the correct photo
                            chunk.markers[j].projections[chunk.cameras[i]] =  PhotoScan.Marker.Projection(PhotoScan.Vector([x,y]), True)		
                            flag = 1
                            break
                    
                    if not flag:
                        print(marker_name)
                        print(" not found")
                        # marker = chunk.addMarker()
                        # marker.label = marker_name
                        # marker.projections[chunk.cameras[i]] =  PhotoScan.Marker.Projection(PhotoScan.Vector([x,y]), True)
                    break
            
            line = file.readline()		#reading the line from input file
            if not len(line):
                eof = True
                break # EOF
        
        file.close()
    
    return
    

#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Reduce Error in images
def reduceErr(chunk, RU_thresh, PA_thresh, RE_thresh, RU_minRemainingPercent, PA_minRemainingPercent, RE_minRemainingPercent, markerCSV):
    print("Enter Reduce Error Stages")
    if uncheckCams:
        for camera in chunk.cameras:
            camera.reference.enabled = False
            
    nestedDictonary = {} # Create empty directory for saving marker error CSV File
    # Save the starting marker error
    beforeRE, nestedDictonary = markerError(chunk, markerCSV, nestedDictonary, 1)
    beforeRE_Tie = str(len(chunk.point_cloud.points))
    
    if gradualSelectionFlag == 1:
        chunk.updateTransform()
        RUT = ReduceError_RU(chunk, RU_thresh, RU_minRemainingPercent)
        PAT = ReduceError_PA(chunk, PA_thresh, PA_minRemainingPercent)
        RET = ReduceError_RE(chunk, RE_thresh, RE_minRemainingPercent)        
        chunk.updateTransform()
    else:
        fltr = PhotoScan.PointCloud.Filter()
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ReconstructionUncertainty)
        RUT = fltr.max_value
        
        fltr = PhotoScan.PointCloud.Filter()
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ProjectionAccuracy)
        PAT = fltr.max_value
        
        fltr = PhotoScan.PointCloud.Filter()
        fltr.init(chunk, PhotoScan.PointCloud.Filter.ReprojectionError)
        RET = fltr.max_value
    
    ## Disable GCP if they have fewer than 3 projections
    markersUsed = len(list(chunk.markers))
    if removeMarkersFlag:
        for marker in chunk.markers:        
            if len(marker.projections.keys()) <= 2:
                marker.reference.enabled = False
                markersUsed -= 1
                    
    # Save the marker error to csv
    afterRE, nestedDictonary = markerError(chunk, markerCSV, nestedDictonary, 0)
    afterRE_Tie = str(len(chunk.point_cloud.points))
    
    return RUT, PAT, RET, beforeRE_Tie, beforeRE, afterRE_Tie, afterRE, markersUsed


#####----------------------------------------------------------------------------------------------------------------------#######
def AlignPhoto(chunk):   
    # Test if the photos have already been aligned
    if chunk.point_cloud is None:
        chunk.crs = img_crs
        # Might want to redefine the camera coordinate system before alignment if the output coordinate system is different!
        # Set the camera and marker projections
        # if img_crs != out_crs:
        #     for camera in chunk.cameras:
        #         if camera.reference.location:
        #             camera.reference.location = PhotoScan.CoordinateSystem.transform(camera.reference.location, img_crs, out_crs)
        
        if len(chunk.cameras) > 0:
            ImgQualAvg, lowQualImg = img_qualityControl(chunk)
            
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
        
            # enabling rolling shutter compensation
            try:
                if camRollingShutter:
                    for sensor in chunk.sensors:
                        sensor.rolling_shutter = True
            except:
                print("Error applying rolling shutter correction")
            
            #Remove unused variables
            realign_list = None
            aligned_list = None
            
            # Set bounding box to area where camera positions are
            if customBB:
                BBCameraLoc(chunk)
        # If there are no cameras found
        else:
            align2 = 'N/A'
            ImgQualAvg = 'N/A'
            lowQualImg = 'N/A'
        
    # If the images were already aligned
    else:
        if len(chunk.cameras) > 0:
            ImgQualAvg, lowQualImg = img_qualityControl(chunk)
        aligned_list = list()
        for camera in chunk.cameras:
            if camera.enabled:
                if camera.transform:
                    aligned_list.append(camera)
            else:
                camera.enabled = True
                if camera.transform:
                    aligned_list.append(camera)
                camera.enabled = False
                
        align2 = len(aligned_list)
        
    return align2, ImgQualAvg, lowQualImg
    

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDenseCloud(chunk, Quality, FilterMode):
    print("Build Dense Cloud")
    try:
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
        
            task = PhotoScan.Tasks.FilterDenseCloud()
            task.point_spacing = 0.001
            task.apply(chunk)
    except:
        print("Error building dense cloud")
    

#####----------------------------------------------------------------------------------------------------------------------#######  
def BuildModel(chunk, modelFN):
    try:
        if chunk.model is None:
            if chunk.dense_cloud is None:
                print("No pointcloud to build model from.")
            else:
                chunk.buildModel(surface_type=Surface, 
                                 interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                                 face_count=PhotoScan.FaceCount.HighFaceCount, 
                                 source_data=SurfaceSource, 
                                 vertex_colors=True)
    except:
        if chunk.model is None:
            if chunk.dense_cloud is None:
                print("No pointcloud to build model from.")
            else:
                chunk.buildModel(surface_type=Surface, 
                                 interpolation=PhotoScan.Interpolation.EnabledInterpolation, 
                                 face_count=PhotoScan.FaceCount.HighFaceCount, 
                                 source_data=PhotoScan.DataSource.DenseCloudData, 
                                 vertex_colors=True)

        
#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDEM(chunk, demFN):
    try:
        if chunk.elevation is None:
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
        if chunk.elevation is None:
            # If the point cloud is not classified
            chunk.buildDem(source_data=PhotoScan.DataSource.DenseCloudData, 
                            interpolation=PhotoScan.Interpolation.EnabledInterpolation, resolution=0.01)
            # task = PhotoScan.Tasks.BuildDem()
            # task.source_data = PhotoScan.DataSource.DenseCloudData
            # task.interpolation = PhotoScan.Interpolation.Extrapolated
            # task.projection = chunk.crs
            # # task.network_distribute = True
            # task.apply(chunk)
        
    # Export DEM
    if saveDEM == 1:
        if not chunk.elevation:
            print("No Elevation dataset")
        else:
            # If file already exists, remove it
            if os.path.exists(demFN):
                os.remove(demFN)
            try:
                # Metashape 1.5
                chunk.exportDem(demFN, image_format=PhotoScan.ImageFormatTIFF,nodata=-99999)
            except:
                # Metashape 1.6
                compression = PhotoScan.ImageCompression()
                compression.tiff_big = True
                # chunk.exportRaster(path = demFN, source_data=PhotoScan.ModelData, image_compression = compression)
                chunk.exportRaster(demFN, source_data=PhotoScan.ElevationData,image_format=PhotoScan.ImageFormatTIFF,nodata_value=-99999)


#####----------------------------------------------------------------------------------------------------------------------#######
def BuildMosaic(chunk, BlendingMode, saveOrthoFN, save_ortho):
    try:
        if chunk.orthomosaic is None:
            # Metashape 1.6
            chunk.buildOrthomosaic(surface_data=PhotoScan.DataSource.ElevationData, 
                                    blending_mode=BlendingMode,
                                    fill_holes=True, 
                                    projection= chunk.crs)
        

    except:
        if chunk.orthomosaic is None:
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
           
    # Export Orthomosaic    
    if save_ortho == 1:
        if not chunk.orthomosaic:
            print("No orthomosaic to save.")
        else:
            # If file already exists, remove it
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
                
                # WGS_84 = Metashape.CoordinateSystem("EPSG::4326")
            
                # # # Degree using EPSG:4326
                # # X_1M_IN_DEG = 1.13747e-05
                # # Y_1M_IN_DEG = 9.0094e-06
            
                # # # Degree using EPSG:4326
                # # X_5M_IN_DEG = 5.76345e-05
                # # Y_5M_IN_DEG = 4.50396e-05
            
                # # EXPORT_DEFAULTS = dict(
                # #     image_format=Metashape.ImageFormat.ImageFormatTIFF,
                # #     projection=WGS_84,
                # #     dx=X_1M_IN_DEG,
                # #     dy=Y_1M_IN_DEG,
                # # )
    return
    

#####----------------------------------------------------------------------------------------------------------------------#######  
def ClassifyGround(chunk, Max_Angle, Cell_Size):
    DEM_resolution, Image_resolution = GetResolution(chunk)
    chunk.dense_cloud.classifyGroundPoints(max_angle=Max_Angle, 
                                           max_distance=2*Image_resolution, 
                                           cell_size=Cell_Size)
    
#####----------------------------------------------------------------------------------------------------------------------#######
def ClassifyPointCloud(chunk):
    chunk.dense_cloud.classifyPoints( source = PhotoScan.PointClass.Created, target = [PhotoScan.PointClass.Ground, PhotoScan.PointClass.HighVegetation, PhotoScan.PointClass.Building], confidence = 0.1)

#####----------------------------------------------------------------------------------------------------------------------#######
def savePtCloudFun(chunk, savePtCloudFN, new_crs):
    # Export point cloud
    if not chunk.dense_cloud:
        print("No point cloud to save")
    else:
        # If file already exists, remove it
        if os.path.exists(savePtCloudFN):
            os.remove(savePtCloudFN)
            
        task = PhotoScan.Tasks.ExportPoints()
        task.source_data = PhotoScan.DataSource.DenseCloudData
        #task.format = PhotoScan.PointsFormat.PointsFormatPLY
        task.format = PhotoScan.PointsFormat.PointsFormatLAS
        task.save_colors = True
        task.crs = new_crs
        task.path = savePtCloudFN
        task.apply(chunk)
        

#####----------------------------------------------------------------------------------------------------------------------#######
def GetResolution(chunk):
    DEM_resolution = float(chunk.dense_cloud.meta['BuildDenseCloud/resolution']) * chunk.transform.scale
    Image_resolution = DEM_resolution / int(chunk.dense_cloud.meta['dense_cloud/depth_downscale'])
    return DEM_resolution, Image_resolution

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RU(chunk, RU_thresh, minRemainingPercent):
    # This is used to reduce error based on reconstruction uncertainty
    # Min percentage of points to be remaining after RU
    tie_points = chunk.point_cloud
    temp_points = tie_points.points
    remainingPoints = int(len(temp_points) * (minRemainingPercent/100.0))
    maxRemovedPoints = int(len(temp_points) - remainingPoints)
    
    
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReconstructionUncertainty)
    threshold = RU_thresh
    print("reconstruction uncertainty")
    nselected = 0
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in temp_points if p.selected])
        # If too many points are selected for removal, clear the selection and up the threshold
        if nselected >= maxRemovedPoints and threshold <= 40:
            fltr.resetSelection()
            threshold += 1
            continue
        elif nselected < maxRemovedPoints:
            break
        elif threshold > 40:
            break
        else:
            break
        
    # UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=True, tiepoint_covariance=False)


    return threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_PA(chunk, PA_thresh, minRemainingPercent):
    # This is used to reduce error based on projection accuracy
    tie_points = chunk.point_cloud
    temp_points = tie_points.points
    remainingPoints = int(len(tie_points.points) * (minRemainingPercent/100.0))
    maxRemovedPoints = int(len(temp_points) - remainingPoints)
    
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ProjectionAccuracy)
    threshold = PA_thresh
    print("projection accuracy")
    nselected = 0
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        # If we get to the threshold of 30 all selected points will be removed regardless
        if nselected >= maxRemovedPoints and threshold <= 6:
            fltr.resetSelection()
            threshold += 0.25
            continue
        elif nselected < maxRemovedPoints:
            break
        elif threshold > 6:
            break
        else:
            break
        
    # UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=True, tiepoint_covariance=False)
    
    # # This is to tighten tie point accuracy value
    # chunk.tiepoint_accuracy = 0.1

    return threshold

#####----------------------------------------------------------------------------------------------------------------------#######
def ReduceError_RE(chunk, RE_thresh, minRemainingPercent):
    print("reprojection error")
    # This is used to reduce error based on repeojection error
    # minRemainingPercent = 80 #percentage of left points
    tie_points = chunk.point_cloud
    temp_points = tie_points.points
    remainingPoints = int(len(tie_points.points) * (minRemainingPercent/100.0))
    maxRemovedPoints = int(len(temp_points) - remainingPoints)
    
    fltr = PhotoScan.PointCloud.Filter()
    fltr.init(chunk, PhotoScan.PointCloud.Filter.ReprojectionError)
    threshold = RE_thresh
    print("Reprojection Error")
    nselected = 0
    while fltr.max_value > threshold:
        fltr.selectPoints(threshold)
        nselected = len([p for p in tie_points.points if p.selected])
        if nselected >= maxRemovedPoints  and threshold <= 2.0:
            fltr.resetSelection()
            threshold += 0.05
            continue
        elif nselected < maxRemovedPoints:
            break
        elif threshold > 2.0: # Some arbitrarily high value
            break
        else:
            break
    
    # UnselectPointMatch(chunk)
    print('Delete {} tie point(s)'.format(nselected))
    tie_points.removeSelectedPoints()
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=True, tiepoint_covariance=False)

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
def markerError(chunk, markerCSV, nestedDictonary, BeforeFlag):
    #total_error = PhotoScan.Vector([0,0,0])
    total_error = 0
    n = 0
    for marker in chunk.markers:
        if(marker.reference.enabled):
            if(marker.position):
               	est = chunk.transform.matrix.mulp(marker.position)
                ref = chunk.crs.unproject(marker.reference.location)
                m = chunk.crs.localframe(chunk.transform.matrix.mulp(marker.position))
               	error = m.mulv(est - ref)
               	print(marker.label, error.norm(), error.x, error.y, error.z)
               	#total_error += PhotoScan.Vector([error.x ** 2, error.y **2, error.z **2])
                total_error += error.norm()**2
               	n += 1
                
                # Calculate pixel error for each marker
                total =0
                num = 0
                for camera in chunk.cameras:
                    if not camera in marker.projections.keys():
                        continue
                    elif camera.enabled and camera.transform:
                        try:
                            v_proj = marker.projections[camera].coord #2 dimensional vector of the marker projection on the photo
                            v_reproj = camera.project(marker.position) #2 dimensional vector of projected 3D marker position
                        
                            diff = (v_proj - v_reproj).norm() #reprojection error for current photo
                            total += diff ** 2
                            num +=1
                        except:
                            continue
                    else:
                        continue
                try:
                    pixErr = math.sqrt(total / num)
                except:
                    pixErr = 999999
                        
                
                
                if markerErrOutput:
                    if BeforeFlag:
                        marker_dict = {"Name" : str(marker.label),
                                    "Proj_keys" : len(marker.projections.keys()),
                                    "Used_4_Recon" : str(marker.reference.enabled),
                                    "X Error before (m)" : error.x,
                                    "Y Error before (m)" : error.y,
                                    "Z Error before (m)" : error.z,
                                    "Total Error before (m)" : error.norm(),
                                    "Reprojection Error (pix)" : pixErr}
                        nestedDictonary[str(marker.label)] = marker_dict
                    else:
                        nestedDictonary[str(marker.label)]["Used_4_Recon"] = str(marker.reference.enabled)
                        nestedDictonary[str(marker.label)]["X Error after (m)"] = error.x
                        nestedDictonary[str(marker.label)]["Y Error after (m)"] = error.y
                        nestedDictonary[str(marker.label)]["Z Error after (m)"] = error.z
                        nestedDictonary[str(marker.label)]["Total Error after (m)"] = error.norm()
                        nestedDictonary[str(marker.label)]["Reprojection Error (pix)"] = pixErr
                        
                        writeMarkerCSV(nestedDictonary[str(marker.label)], markerCSV)
    # End Marker for loop
    
    # Calculate RMSE of norms (Overall marker errors)
    try:
        ErrorTotal = math.sqrt( total_error / n)
    except:
        ErrorTotal = 99999
    
    return ErrorTotal, nestedDictonary


#####----------------------------------------------------------------------------------------------------------------------#######
def cameraError(chunk, cameraErrOutput, cameraCSV):
    for camera in chunk.cameras:
        camera.reference.enabled = True
        
    #### Write camera errors
    cameraDictonary = {}
    total_error = PhotoScan.Vector([0,0,0])
    n = 0
    fieldnames = ["Name", "X Error(m)", "Y Error(m)", "Z Error(m)", "Total Error(m)"]
    
    for camera in chunk.cameras:
        if(camera.enabled):
            if camera.center:
                est = chunk.transform.matrix.mulp(camera.center)
                ref = chunk.crs.unproject(camera.reference.location)
                m = chunk.crs.localframe(chunk.transform.matrix.mulp(camera.center))
                if (est and ref):
                    error = m.mulv(est - ref)
                    camera_dict = {"Name" : str(camera.label), 
                                   "X Error(m)" : error.x, 
                                   "Y Error(m)" : error.y, 
                                   "Z Error(m)" : error.z, 
                                   "Total Error(m)" : error.norm()}
                    cameraDictonary[str(camera.label)] = camera_dict
                                            
                    total_error += PhotoScan.Vector([error.x ** 2, error.y **2, error.z **2])
                    n += 1
                    
                    with open(cameraCSV, 'a', newline='') as csv_file:
                        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
                        if os.path.getsize(cameraCSV) == 0:
                            print('csv file is empty')
                            writer.writeheader()
                        writer.writerow(camera_dict)
    
    for i in range(len(total_error)): #printing total X, Y, Z errors
        ErrorTotal = math.sqrt(total_error[i] / n)
        

    return total_error, ErrorTotal

#####----------------------------------------------------------------------------------------------------------------------#######
def writeMarkerCSV(m_dict, mcsv_FN):
    with open(mcsv_FN, 'a', newline='') as csv_file:
        fieldnames = ['Name', 'Proj_keys', "Used_4_Recon", "X Error before (m)", "Y Error before (m)", "Z Error before (m)", "Total Error before (m)",
                      "Reprojection Error (pix)", "X Error after (m)", "Y Error after (m)", "Z Error after (m)", "Total Error after (m)", "Reprojection Error (pix)"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
        # If there is no header, write the header
        if os.path.getsize(mcsv_FN) == 0:
            print('csv file is empty')
            writer.writeheader()
        
        writer.writerow(m_dict)

#####----------------------------------------------------------------------------------------------------------------------#######
def writeCameraCSV(c_dict, ccsv_FN):
    fieldnames = ["Name", "X Error(m)", "Y Error(m)", "Z Error(m)", "Total Error(m)"]
    with open(ccsv_FN, 'a', newline='') as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
        # If there is no header, write the header
        if os.path.getsize(ccsv_FN) == 0:
            print('csv file is empty')
            writer.writeheader()
        
        writer.writerow(c_dict)
    
#####----------------------------------------------------------------------------------------------------------------------#######
def writeCSV(temp_dict, csv_FN):
    with open(csv_FN, 'a', newline='') as csv_file:
        fieldnames = ['Name', 'Accuracy', 'Key_Limit', 'Tie_Limit','Quality','Filter','Aligned','Total','lowQualityImg', 'Img_Quality_Avg', "MarkerCount",
                      "MarkersUsed", "ProcessTime", 'RU_Thresh', 'PA_Thresh','RE_Thresh', 'Tie_ptsBefore','Marker_errBefore','Tie_ptsAfter','Marker_errAfter',
                      "Cam Tot err", 'Cal_f', 'Cal_x', 'Cal_y']

        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
        # If there is no header, write the header
        if os.path.getsize(csv_FN) == 0:
            print('csv file is empty')
            writer.writeheader()
        
        writer.writerow(temp_dict)
    
    
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


####--------------------------------------------------------------------------------------------------------------------####
def processImgSet(doc, chunk, dirPath, saveFold, gcpPath, csv_outfile, markerCSV, cameraCSV):
    start = time.time()
    # Define or reset dictionary items
    align_Img = tot_Img = lowQualImg = ImgQual = markerCount = markersUsed = processTime = 0
    RUT = PAT = RET = beforeRE_Tie = beforeRE = afterRE_Tie = afterRE = ErrorTotal = calf = calx = caly = 0
    total_error = [0,0,0] #marker X,Y,Z errors
    
    print("Started Processing " + str(chunk.label))
    
    # Find and load photos
    tot_Img = loadImgs(chunk, dirPath, imgExt)
    
    # Align the Photos
    align_Img, ImgQual, lowQualImg = AlignPhoto(chunk)
    doc.save()
    
    ## Automatic Target Detection
    # markerCount = auto_mark_Images(chunk, gcpPath)
           
    # Load GCPs
    markerCount = loadGCP(chunk, gcpPath, dirPath)
    doc.save()
    
    # Reduce Error through gradual selection if enabled
    # Calculates marker errors and outputs to a CSV if enabled
    RUT, PAT, RET, beforeRE_Tie, beforeRE, afterRE_Tie, afterRE, markersUsed = reduceErr(chunk, RU_thresh, PA_thresh, RE_thresh,
                                    RU_minRemainingPercent, PA_minRemainingPercent, RE_minRemainingPercent, markerCSV)
    doc.save()
    
    # Calculates and outputs a CSV file with individual camera errors (useful for statistics reporting)
    if cameraErrOutput:
        total_error, ErrorTotal = cameraError(chunk, cameraErrOutput, cameraCSV)
        doc.save()
    
    # Build Dense Cloud
    # If the dense cloud already exists it will just save
    if buildDenseFlag:
        BuildDenseCloud(chunk, Quality, FilterMode)
        doc.save()
    
    if classifyPointsFlag:
        ClassifyPointCloud(chunk)
    
    if classifyGroundFlag:
        ClassifyGround(chunk, Max_Angle, Cell_Size)
    
    if savePointCloud:
        savePtCloudFN = saveFold + str(chunk.label) + ".las" # was saveFold + saveID + ".las"
        savePtCloudFun(chunk, savePtCloudFN, out_crs)
    
    # Build Mesh
    if buildMeshFlag:
        modelFN = saveFold + str(chunk.label) + "_Model.obj"
        BuildModel(chunk, modelFN)
        doc.save()
    
    # Build DEM
    if buildDEMFlag:
        demFN = saveFold + str(chunk.label) + "_DEM.tif"
        BuildDEM(chunk, demFN)
        doc.save()
    
    # Create OrthoMosaic
    if createOrthoFlag:
        orthoFN = saveFold + str(chunk.label) + "_ortho.tif"
        BuildMosaic(chunk, BlendingMode, orthoFN, saveOrtho)
        doc.save()
    
    #Save the camera calibration to csv
    if chunk.sensors:
        if chunk.sensors[0].calibration:
            calf = chunk.sensors[0].calibration.f
            calx = chunk.sensors[0].calibration.cx
            caly = chunk.sensors[0].calibration.cy
    
    
    end = time.time()
    print("Finished Processing " + str(chunk.label))
    
    # Save the results to a CSV file
    procTime = end - start
    processTime = secondsToText(procTime)

    
    temp_dict= {"Name" : str(chunk.label),
                "Accuracy" : accuracyLvl,
                "Key_Limit" : Key_Limit,
                "Tie_Limit" : Tie_Limit,
                "Quality" : qualityLevel,
                "Filter" : filterLevel,
                "Aligned" : align_Img,
                "Total" : tot_Img,
                "lowQualityImg" : lowQualImg,
                "Img_Quality_Avg" : ImgQual,
                "MarkerCount" : markerCount,
                "MarkersUsed" : markersUsed,
                "ProcessTime" : processTime,
                "RU_Thresh" : RUT,
                "PA_Thresh" : PAT,
                "RE_Thresh" : RET,
                "Tie_ptsBefore" : beforeRE_Tie,
                "Marker_errBefore" : beforeRE,
                "Tie_ptsAfter" : afterRE_Tie,
                "Marker_errAfter" : afterRE,
                "Cam Tot err" : ErrorTotal,
                "Cal_f" : calf,
                "Cal_x" : calx,
                "Cal_y" : caly}
        
    writeCSV(temp_dict, csv_outfile)
    
    if saveReport:
        saveReportFN = saveFold + saveID + "_" + str(chunk.label) + ".pdf"
        try:
            chunk.exportReport(saveReportFN, title=str(chunk.label), description=saveID, page_numbers=True)
        except:
            print("Failed to save the Report")
    


####--------------------------------------------------------------------------------------------------------------------####
def iterFold():
    savePath = saveDir + saveID + "/"   #ex <saveDir> + /2021/ + <saveID> + /
    directory = os.path.dirname(savePath)
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    # Where to save the CSV File
    csv_outfile = savePath + saveID + ".csv"
    # Where to save the PSX File
    psxfile = savePath + saveID + ".psx"
    
    # Open new Metashape document
    doc = PhotoScan.Document()
    # Try to open an existing project file
    try: 
        doc.open( psxfile, read_only=False, ignore_lock=True )
    # Save to a new project file
    except: 
        doc.save( psxfile )
            
    # Go through folder structure
    for foldPath in sorted(glob.iglob(mainFold + '2023-06-26*/')):
        foldPath = foldPath.replace('\\', '/') #ex .../2021/
        chunkName = foldPath.rpartition("/")[0].rpartition("/")[2] #ex 2021-03-23
        cNameExists = 0
        # Check if the chunk already exists
        if len(doc.chunks):
            for item in list(doc.chunks):
                if item.label == chunkName:
                    chunk = item
                    cNameExists = 1
                    break
        if cNameExists == 1:
            pass
        else:
            chunk = doc.addChunk()
            chunk.label = chunkName
            #Define the project coordinate system
            chunk.crs = out_crs
        
        # Save marker CSV file location
        markerCSV = savePath + saveID + "_MarkerErrors_" + chunkName + ".csv"
        cameraCSV = savePath + saveID + "_CameraErrors_" + chunkName + ".csv"
        
        #GCP File.
        gcpPath = glob.glob(foldPath + "GPSDATA/Rover/*_Agisoft.csv")[0]
        gcpPath = gcpPath.replace('\\', '/')
        
        # Process the data
        # dirPath = foldPath + 'P4RTK/'
        dirPath = foldPath+ 'SfM/'
        processImgSet(doc, chunk, dirPath, savePath, gcpPath, csv_outfile, markerCSV, cameraCSV)           


    return

                  
        
""" 
#####----------------------------------------------------------------------------------------------------------------------------------------#######
#####----------------------------------------------------------------------------------------------------------------------------------------#######
#####----------------------------------------------------- Main Method ----------------------------------------------------------------------#######
#####----------------------------------------------------------------------------------------------------------------------------------------#######
#####----------------------------------------------------------------------------------------------------------------------------------------#######
"""
if __name__ == '__main__':
    print("Starting Script...")
    
    #####-----------------------------------Agisoft variables---------------------------------------------------------------------#######
    # Do not change these unless needed for a different version of Agisoft Metashape
    
    ## Variables for photo alignment
    ##  Metashape 1.6 Accuracy Variables
    ALIGN = {"Highest":  0,
              "High":   1,
              "Medium": 2,
              "Low":    4,
              "Lowest": 8}
    
    # Variables for building dense cloud
    # Metashape 1.6 Quality Variables
    DENSE = {"Ultra":  1,
              "High":   2,
              "Medium": 4,
              "Low":    8,
              "Lowest": 16}
    
    # Filter: AggressiveFiltering, ModerateFiltering, MildFiltering, NoFiltering
    # Metashape 1.5 and 1.6 Filtering Variables
    FILTERING = {"None":  PhotoScan.NoFiltering,
                 "Mild":  PhotoScan.MildFiltering,
                 "Moderate": PhotoScan.ModerateFiltering,
                 "Aggressive": PhotoScan.AggressiveFiltering}
    
    
    """
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####------------------------------------------------- User Defined Parameters --------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    """
    ## Set the Coodinate Systems
    # PhotoScan.CoordinateSystem("EPSG::32611") #UTM11N
    # PhotoScan.CoordinateSystem("EPSG::4326")  #WGS84
    # PhotoScan.CoordinateSystem("EPSG:7912")  #ITRF2014
    # PhotoScan.CoordinateSystem('LOCAL_CS["Local CS",LOCAL_DATUM["Local Datum",0],UNIT["metre",1]]') # Local coordinate system
    gcp_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of the GCPs
    gcp_crs_UTM = PhotoScan.CoordinateSystem("EPSG::32611") #Coordinate System of the GCPs
    img_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of Image geotags
    out_crs = PhotoScan.CoordinateSystem("EPSG::32611") #Desired Output Coordinate system
    
    
    ## Set GCP and camera accuracies in meters
    gcp_ref_acc = 0.2
    cam_ref_acc = 3
    
    ## Image quality filter
    # QualityFilter: True, False
    # QualityCriteria: float number range from 0 to 1
    QualityFilter = True
    QualityCriteria = 0.7
    
    ## Parameters for Image Alignment
    # Accuracy: Highest, High, Medium, Low, Lowest
    accuracyLvl = "High"
    Accuracy = ALIGN[accuracyLvl]
    Key_Limit = 40000
    Tie_Limit = 4000
    # Enable rolling shutter correction
    camRollingShutter = False
    
    ## Parameters for the Dense Cloud
    # Quality: Ultra, High, Medium, Low, Lowest
    qualityLevel = 'High'
    Quality = DENSE[qualityLevel]
    # Filtering: None, Mild, Moderate, Aggressive
    filterLevel = 'Mild'
    FilterMode = FILTERING[filterLevel]
    # Variables for dense cloud ground point classification
    Max_Angle = 13           # Metashape default is 13
    Cell_Size = 10           # Metashape default is 10
    
    
    # Variables for building 3D mesh
    # SurfaceType: Arbitrary, HeightField
    # DataSource: PointCloudData, DenseCloudData, DepthMapsData
    Surface = PhotoScan.SurfaceType.Arbitrary
    SurfaceSource = PhotoScan.DataSource.DepthMapsData
    
    # Variable for building orthomosaic
    # BlendingMode: AverageBlending, MosaicBlending, MinBlending, MaxBlending, DisabledBlending
    BlendingMode = PhotoScan.BlendingMode.MosaicBlending
    
    ## Reduce Error Parameters
    RU_thresh = 1  #Start at this threshold, go up by 1 until it reaches RU_minRemainingPercent of remaining points OR RU_thresh = 40
    PA_thresh = 0.5   #Start at this threshold, go up by 0.25 until it reaches PA_minRemainingPercent of remaining points OR  PA_thresh = 6
    RE_thresh = 0.2 #Start at this threshold, go up by 0.05 until it reaches RE_minRemainingPercent of remaining points OR RE_thresh = 2
    # Minimum percentage of points remaining after each RE step
    RU_minRemainingPercent = 70
    PA_minRemainingPercent = 70
    RE_minRemainingPercent = 90
    # Flag to uncheck cameras before gradual selection to reduce error
    uncheckCams = 0
    # After gradual selection remove markers that have fewer than 3 projections
    removeMarkersFlag = 0

    
    ## Flags for what steps to run
    saveID = "HighQual2023v2" # Used to name the project file (e.g. <saveID>.psx)
    markerRenew = 0          # Deletes existing markers, used to load updated markers
    # Processes to run after GCP Flag refinement
    FlagsRefined = 1
    if FlagsRefined:
        gradualSelectionFlag = 0 # Set to 0 if gradual selection of the point cloud was already run to skip this step
        customBB = 1             # Set a custom bounding box around camera locations
        loadMarkers = 0          # Load GCPs
        markerErrOutput = 1      # Save a CSV file with details on individual marker accuracies
        cameraErrOutput = 1      # Save a CSV file with details on individual camera location accuracies
        buildDenseFlag = 1       # Build Dense Cloud
        savePointCloud = 1       # Save a point cloud
        classifyPointsFlag = 1   # Must have build dense cloud flag
        classifyGroundFlag = 0   # Classify PointCloud
        buildMeshFlag = 1        # Build Mesh
        buildDEMFlag = 1         # Build DEM
        saveDEM = 1              # Save a DEM
        createOrthoFlag = 1      # Create OrthoMosaic
        saveOrtho = 1            # Save an orthophoto
        saveReport = 1           # Save the automatically generated Agisoft report as a PDF
    else:
        gradualSelectionFlag = 0 # Set to 0 if gradual selection of the point cloud was already run to skip this step
        customBB = 0             # Set a custom bounding box around camera locations
        loadMarkers = 1          # Load GCPs
        markerErrOutput = 0      # Save a CSV file with details on individual marker accuracies
        cameraErrOutput = 0      # Save a CSV file with details on individual camera location accuracies
        buildDenseFlag = 0       # Build Dense Cloud
        savePointCloud = 0       # Save a point cloud
        classifyPointsFlag = 0   # Must have build dense cloud flag
        classifyGroundFlag = 0   # Classify PointCloud
        buildMeshFlag = 0        # Build Mesh
        buildDEMFlag = 0         # Build DEM
        saveDEM = 0              # Save a DEM
        createOrthoFlag = 0      # Create OrthoMosaic
        saveOrtho = 0            # Save an orthophoto
        saveReport = 0           # Save the automatically generated Agisoft report as a PDF
    
    ## Folder Name and Location
    # header folder
    # mainFold = "/SNOWDATA/SnowDrones-Processing/2022_Data/SurveyDates/Process/"
    # mainFold = "/SNOWDATA/SnowDrones-Processing/GrandMesa_2023/SfM/DCIM/SURVEY/"
    mainFold = "/SNOWDATA/SnowDrones-Processing/2023_FineFuels/FieldDays/" # This is the folder where the dates surveyed resides
    # mainFold = "E:/SnowDrones_HDrive/2022_Data/SurveyDates/"
    imgExt = ".JPG"
    
    # Define where to save the photoscan project and related files (e.g. ortho)
    # saveDir = "/scratch/thomasvanderweide/SnowDrones/2022/"
    saveDir = "/scratch/thomasvanderweide/SnowDrones/2023_FineFuels/"
    # saveDir = "E:/SnowDrones_HDrive/2022_Data/SurveyDates/Output/"
    
    # Define existing marker file name if the flags have already been set for images
    # markerFN = "GCP_Img_XY_loc.csv" # Set to False if there is no file to load
    markerFN = "False" # Set to False if there is no file to load
    
    # Run the script over all folders in mainFold
    iterFold()
    
    # Closes the existing doc
    doc = PhotoScan.Document()
    
    # #Close the App
    # app = PhotoScan.Application()
    # app.quit()

