#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 3 2021
@author: Thomas Van Der Weide
This Python Script is developed for Agisoft MetaShape 1.6 (can be made be backward compatible with 1.5.5)
Python core is 3.8.1
Purpose: Automated processing of BCAL small SfM biomass project with automatic target detection

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
    Import a polygon to define the bounding box
    For large datasets it's best Practice to align all images together and then split into chunks.
        Split into chunks for building point cloud (https://github.com/agisoft-llc/metashape-scripts/blob/master/src/align_model_to_model.py)
    Parameter sensitivity iterations
        Optional format to re-run the code with different parameters for comparison
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
import math
import statistics


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
        
        for fold in sorted(glob.iglob(root_path)):
            fold = fold.replace('\\', '/')
                            
            for root, dirs, files in os.walk(fold):
                for name in files:
                    if re.search(pattern,name):
                        cur_path = os.path.join(root, name)
                        photoList.append(cur_path)
                        
        # Load photos
        if photoList:
            chunk.addPhotos(photoList)
            # Set the defined camera accuracy
            for camera in chunk.cameras:
                camera.reference.accuracy = PhotoScan.Vector([cam_ref_acc, cam_ref_acc, cam_ref_acc])
            
            imgCount = len(photoList)
        else:
            imgCount = "N/A"
    
    return imgCount


#####-----------------------------------------------------------------------------------------------------------------------------------#######
def lowQualImgCount(chunk):
    lowQualImg = 0
    if QualityFilter:
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

def img_qualityControl(chunk, QualityCriteria, QualityFilter):
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
    
    lowQualImg = lowQualImgCount(chunk)
    images_quality_avg = estImgQual(chunk)

    
    return images_quality_avg, lowQualImg


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def auto_mark_Images(chunk, gcpPath):
    markerCount = len(chunk.markers)
    if markerCount == 0:
        print("Marking Images...")
        marker_repeat = 0
        global crs
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
        
        ##Resize Bounding box region
        # T = chunk.transform.matrix
    
        # if chunk.crs:
        	# v_t = T * PhotoScan.Vector( [0,0,0,1] )
        	# v_t.size = 3
        	# m = chunk.crs.localframe(v_t)
        # else:
        	# m = PhotoScan.Matrix().diag([1,1,1,1])
        
        # m = m * T
        
        # s = math.sqrt(m[0,0]**2 + m[0,1]**2 + m[0,2]**2) #scale factor
        
        # R = PhotoScan.Matrix( [[m[0,0],m[0,1],m[0,2]], [m[1,0],m[1,1],m[1,2]], [m[2,0],m[2,1],m[2,2]]])
        # R = R * (1. / s)
        
        # reg = chunk.region
        # reg.rot = R.t()
        # chunk.region = reg
    
        R = chunk.region.rot     # Bounding box rotation matrix
        C = chunk.region.center  # Bounding box center vector
    
        if chunk.transform.matrix:
            T = chunk.transform.matrix
            s = math.sqrt(T[0, 0] ** 2 + T[0, 1] ** 2 + T[0, 2] ** 2)  # scaling # T.scale() 
            S = PhotoScan.Matrix().Diag([s, s, s, 1])                  # scale matrix
        else:
            S = PhotoScan.Matrix().Diag([1, 1, 1, 1])
    
        T = PhotoScan.Matrix([[R[0, 0], R[0, 1], R[0, 2], C[0]],
                              [R[1, 0], R[1, 1], R[1, 2], C[1]],
                              [R[2, 0], R[2, 1], R[2, 2], C[2]], #was C[2]
                              [      0,       0,       0,    1]])
    
        chunk.transform.matrix = S * T.inv()  # resulting chunk transformation matrix
    # If the markers have already been found
    else:
        markerCount = len(chunk.markers)
        marker_repeat = 1

    return markerCount, marker_repeat


#####----------------------------------------------------------------------------------------------------------------------#######
# Define bounding box by camera locations
# https://www.agisoft.com/forum/index.php?topic=10102.0
def cross(a, b):
	result = PhotoScan.Vector([a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y *b.x])
	return result.normalized()

def BBCameraLoc(chunk):   
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
    
#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Mark GCP target locations in images
def loadGCP(chunk, gcpPath):
    markerCount = len(chunk.markers)
    if markerCount == 0:
        print("Loading reference targets location...")
        # Load the GCP Locations
        try:
            chunk.importReference(gcpPath, format = PhotoScan.ReferenceFormatCSV, delimiter=",", columns="nxyz",create_markers=True, crs=gcp_crs)
            for marker in chunk.markers:
                # Set the Marker Accuracy
                marker.reference.accuracy = PhotoScan.Vector([gcp_ref_acc,gcp_ref_acc,gcp_ref_acc])
        except:
            print("GCP Coordinate File Not Found.")
            
        ## Disable Header 'GCP'
        for marker in chunk.markers:
            if marker.label == 'Location':
                 marker.reference.enabled = False
            # LDP GCP
            elif marker.label == 'VGCP0':
                 marker.reference.enabled = False
            elif marker.label == 'VGCP1':
                 marker.reference.enabled = False
            #Bogus VGCP
            elif marker.label == 'VGCP':
                 marker.reference.enabled = False
            elif marker.label == 'VGCP4':
                 marker.reference.enabled = False
    
    # Find number of targets (markers) loaded
    markerCount = len(chunk.markers)
    
    # If the targets have already been marked on the images
    markerFN = gcpPath.rpartition("/")[0] + "/RGB/GCPwithImageLocations.csv"
    if os.path.isfile(markerFN):
        mark_Images(chunk, markerFN)
    
    marker_repeat = 0
    return markerCount, marker_repeat


#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Load saved GCP target locations in images
def mark_Images(chunk, path):
    # Open the CSV with GCP locations saved as pixel coordinates on images
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
                    if chunk.markers[j].label.lower() == marker_name.lower():
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
    
    file.close()
    return
    

#####-----------------------------------------------------------------------------------------------------------------------------------#######
# Reduce Error in images
def reduceErr(chunk, RU_thresh, PA_thresh, RE_thresh, marker_repeat, RU_minRemainingPercent, PA_minRemainingPercent, RE_minRemainingPercent, markerErrOutput, markerCSV, removeMarkersFlag):
    print("Enter Reduce Error Stages")
    nestedDictonary = {} # Create empty directory for saving marker error CSV File
    # Save the new marker error
    beforeRE, nestedDictonary = markerError(chunk, markerErrOutput, markerCSV, nestedDictonary, 1)
    beforeRE_Tie = str(len(chunk.point_cloud.points))
    
    if marker_repeat == 0:
        RUT = ReduceError_RU(chunk, RU_thresh, RU_minRemainingPercent)
        PAT = ReduceError_PA(chunk, PA_thresh, PA_minRemainingPercent)
        RET = ReduceError_RE(chunk, RE_thresh, RE_minRemainingPercent)
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
            if len(marker.projections.keys()) <= 3:
                marker.reference.enabled = False
                markersUsed -= 1
    
    # Save the new marker error
    afterRE, nestedDictonary = markerError(chunk, markerErrOutput, markerCSV, nestedDictonary, 0)
    afterRE_Tie = str(len(chunk.point_cloud.points))
    
    return RUT, PAT, RET, beforeRE_Tie, beforeRE, afterRE_Tie, afterRE, markersUsed


#####----------------------------------------------------------------------------------------------------------------------#######
def AlignPhoto(chunk):   
    # Test if the photos have already been aligned
    if chunk.point_cloud is None:
        # Set the camera and marker projections
        if img_crs != out_crs:
            for camera in chunk.cameras:
                if camera.reference.location:
                    camera.reference.location = PhotoScan.CoordinateSystem.transform(camera.reference.location, img_crs, out_crs)
        
        ImgQualAvg, lowQualImg = img_qualityControl(chunk, QualityCriteria, QualityFilter)
        
        
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
        
        chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
                              fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
                              fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
                              adaptive_fitting=False, tiepoint_covariance=False)
        
        #Remove unused variables
        realign_list = None
        aligned_list = None
        
        # Set bounding box to area where camera positions are
        if customBB:
            BBCameraLoc(chunk)
        
    # If the images were already aligned
    else:
        aligned_list = list()
        for camera in chunk.cameras:
              if camera.transform:
                    aligned_list.append(camera)
        align2 = len(aligned_list)
            
        lowQualImg = lowQualImgCount(chunk)
        ImgQualAvg = estImgQual(chunk)
        
        
    return align2, ImgQualAvg, lowQualImg
    

#####----------------------------------------------------------------------------------------------------------------------#######
def BuildDenseCloud(chunk, Quality, FilterMode, savePointCloud, savePtCloudFN, new_crs):

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
        print("farts")
            
    # Export point cloud
    if savePointCloud == 1:
        if not chunk.dense_cloud:
            print("No point cloud to save")
        else:
            # If file already exists, remove it
            if os.path.exists(savePtCloudFN):
                os.remove(savePtCloudFN)
                
            task = PhotoScan.Tasks.ExportPoints()
            task.source_data = PhotoScan.DataSource.DenseCloudData
            task.format = PhotoScan.PointsFormat.PointsFormatPLY      #task.format = PhotoScan.PointsFormat.PointsFormatLAS
            task.save_colors = True
            task.crs = new_crs
            task.path = savePtCloudFN
            task.apply(chunk)



#####----------------------------------------------------------------------------------------------------------------------#######  
def ClassifyGround(chunk, Max_Angle, Cell_Size):
    DEM_resolution, Image_resolution = GetResolution(chunk)
    chunk.dense_cloud.classifyGroundPoints(max_angle=Max_Angle, 
                                           max_distance=2*Image_resolution, 
                                           cell_size=Cell_Size)

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
def GetResolution(chunk):
    DEM_resolution = float(chunk.dense_cloud.meta['dense_cloud/resolution']) * chunk.transform.scale
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
        if nselected >= maxRemovedPoints and threshold <= 50:
            fltr.resetSelection()
            threshold += 1
            continue
        elif nselected < maxRemovedPoints:
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
    # chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=False, fit_b2=False, 
    #                       fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=False, 
    #                       fit_p1=True, fit_p2=True, fit_p3=False, fit_p4=False, 
    #                       adaptive_fitting=False, tiepoint_covariance=False)
    
    # # This is to tighten tie point accuracy value
    # chunk.tiepoint_accuracy = 0.1
    chunk.optimizeCameras(fit_f=True, fit_cx=True, fit_cy=True, fit_b1=True, fit_b2=True, 
                          fit_k1=True, fit_k2=True, fit_k3=True, fit_k4=True, 
                          fit_p1=True, fit_p2=True, fit_p3=True, fit_p4=True, 
                          adaptive_fitting=False, tiepoint_covariance=False)

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
def markerError(chunk, markerErrOutput, markerCSV, nestedDictonary, BeforeFlag):
    errorList = []
    for marker in chunk.markers:
        if(marker.position):
            est = chunk.crs.project(chunk.transform.matrix.mulp(marker.position))  # Gets estimated marker coordinate
            ref = marker.reference.location

            if est and ref:
                error = (est - ref).norm()  # The .norm() method gives the total error. Removing it gives X/Y/Z error
                print(marker.label, error)
                errorList += [error**2]
                
                if markerErrOutput:
                    xyzerr = est - ref
                    if BeforeFlag:
                        marker_dict = {"Name" : str(marker.label),
                                    "Proj_keys" : len(marker.projections.keys()),
                                    "Used_4_Recon" : str(marker.reference.enabled),
                                    "X Error before (m)" : xyzerr[0],
                                    "Y Error before (m)" : xyzerr[1],
                                    "Z Error before (m)" : xyzerr[2],
                                    "Total Error before (m)" : error}
                        nestedDictonary[str(marker.label)] = marker_dict
                    else:
                        nestedDictonary[str(marker.label)]["Used_4_Recon"] = str(marker.reference.enabled)
                        nestedDictonary[str(marker.label)]["X Error after (m)"] = xyzerr[0]
                        nestedDictonary[str(marker.label)]["Y Error after (m)"] = xyzerr[1]
                        nestedDictonary[str(marker.label)]["Z Error after (m)"] = xyzerr[2]
                        nestedDictonary[str(marker.label)]["Total Error after (m)"] = error
                        
                        writeMarkerCSV(nestedDictonary[str(marker.label)], markerCSV)
    
    # Calculate RMSE of norms
    suma = sum(errorList)
    n = len(errorList)
    try:
        ErrorTotal = (suma / n) ** 0.5
    except:
        ErrorTotal = 99999
    
    return ErrorTotal, nestedDictonary


#####----------------------------------------------------------------------------------------------------------------------#######
def writeMarkerCSV(marker_dict, mcsv_FN):
    with open(mcsv_FN, 'a', newline='') as csv_file:
        fieldnames = ['Name', 'Proj_keys', "Used_4_Recon", "X Error before (m)", "Y Error before (m)", "Z Error before (m)", "Total Error before (m)"
                      , "X Error after (m)", "Y Error after (m)", "Z Error after (m)", "Total Error after (m)"]
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
        # If there is no header, write the header
        if os.path.getsize(mcsv_FN) == 0:
            print('csv file is empty')
            writer.writeheader()
        
        writer.writerow(marker_dict)


#####----------------------------------------------------------------------------------------------------------------------#######
def writeCSV(temp_dict, csv_FN):
    with open(csv_FN, 'a', newline='') as csv_file:
        fieldnames = ['Name', 'Accuracy', 'Key_Limit', 'Tie_Limit','Quality','Filter','Aligned','Total','lowQualityImg', 'Img_Quality_Avg', "MarkerCount", "MarkersUsed",
                      "ProcessTime", 'RU_Thresh', 'PA_Thresh','RE_Thresh', 'Tie_ptsBefore','Marker_errBefore','Tie_ptsAfter','Marker_errAfter', 'Cal_f', 'Cal_x', 'Cal_y']
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
def iterFold():
    # Go through folder structure
    for foldPath in sorted(glob.iglob(mainFold + '*/')):
        foldPath = foldPath.replace('\\', '/') #ex .../2021/
        saveName = foldPath.rpartition("/")[0].rpartition("/")[2] #ex 2021
        savePath = saveDir + saveName + "/" + saveID + "/"   #ex <saveDir> + /2021/ + <saveID> + /
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
        
        # Process sub-directories as chunks in project
        for dirPath in sorted(glob.iglob(foldPath + '*/')):
            dirPath = dirPath.replace('\\', '/') #ex .../2021/2021-03-23/
            chunkName = dirPath.rpartition("/")[0].rpartition("/")[2] #ex 2021-03-23
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
            markerCSV = savePath + saveID + "_MarkerError_" + chunkName + ".csv"
            
            #GCP File
            gcpPath = dirPath + "GCPs_uncorrected.csv"
            
            # Process the data
            processImgSet(doc, chunk, dirPath, savePath, gcpPath, csv_outfile, markerCSV)

    return

####--------------------------------------------------------------------------------------------------------------------####
def processImgSet(doc, chunk, dirPath, saveFold, gcpPath, csv_outfile, markerCSV):
    start = time.time()
    # Define or reset dictionary items
    align_Img = tot_Img = lowQualImg = ImgQual = markerCount = markersUsed = processTime = 0
    RUT = PAT = RET = beforeRE_Tie = beforeRE = afterRE_Tie = afterRE = marker_repeat = calf = calx = caly = 0

    # Find and load photos
    tot_Img = loadImgs(chunk, dirPath, imgExt)
    
    ## Automatic Target Detection
    ## marker_repeat = 1 if markers were already loaded. This will assume reduce error was already run and skip this step
    ## Bounding box is automatically defined around detected targets
    # markerCount, marker_repeat = auto_mark_Images(chunk, gcpPath)
    
    # Load GCPs
    markerCount, marker_repeat = loadGCP(chunk, gcpPath)
    doc.save()
    
    # Align the Photos
    align_Img, ImgQual, lowQualImg = AlignPhoto(chunk)
    doc.save()
    
    # # Reduce Alignment Error
    # # marker_repeat = 1 if markers were already loaded. This will assume reduce error was already run and skip this step
    # RUT, PAT, RET, beforeRE_Tie, beforeRE, afterRE_Tie, afterRE, markersUsed = reduceErr(chunk, RU_thresh, PA_thresh, RE_thresh, marker_repeat,
    #                                 RU_minRemainingPercent, PA_minRemainingPercent, RE_minRemainingPercent, markerErrOutput, markerCSV, removeMarkersFlag)
    # doc.save()
    
    # # Build Dense Cloud
    # savePtCloudFN = saveFold + str(chunk.label) + ".ply" # was saveFold + saveID + ".las"
    # BuildDenseCloud(chunk, Quality, FilterMode, savePointCloud, savePtCloudFN, out_crs)
    # doc.save()
    
    # # # Build Mesh
    # # modelFN = saveFold + str(chunk.label) + "_Model.obj"
    # # BuildModel(chunk, modelFN)
    # # doc.save()
    
    # # Build DEM
    # demFN = saveFold + str(chunk.label) + "_DEM.tif"
    # BuildDEM(chunk, demFN)
    # doc.save()
    
    # # Create OrthoMosaic
    # orthoFN = saveFold + str(chunk.label) + "_ortho.tif"
    # BuildMosaic(chunk, BlendingMode, orthoFN, saveOrtho)
    # doc.save()
    
    #Save the camera calibration to csv
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
                
                            
        
    
#####------------------------------------------------------------------------------#####
#####------------------------------- MAIN -----------------------------------------#####
#####------------------------------------------------------------------------------#####
if __name__ == '__main__':
    print("Starting Program...")
    
    #####-----------------------------------Agisoft variables---------------------------------------------------------------------#######
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
    
    """
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####------------------------------------------------- User Defined Parameters --------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    #####----------------------------------------------------------------------------------------------------------------------------------------#######
    """
    ## Set the projection
    # PhotoScan.CoordinateSystem("EPSG::32611") #UTM11N
    # PhotoScan.CoordinateSystem("EPSG::4326")  #WGS84
    # PhotoScan.CoordinateSystem('LOCAL_CS["Local CS",LOCAL_DATUM["Local Datum",0],UNIT["metre",1]]') # Local coordinate system
    gcp_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of the GCPs
    img_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Coordinate System of Image geotags
    out_crs = PhotoScan.CoordinateSystem("EPSG::4326") #Desired Output Coordinate system
    
    ## Set GCP and camera accuracies in meters
    gcp_ref_acc = 5.0
    cam_ref_acc = 10.0
    
    ## Variables for image quality filter
    # QualityFilter: True, False
    # QualityCriteria: float number range from 0 to 1
    QualityFilter = True
    QualityCriteria = 0.7
    
    # Enable rolling shutter correction
    camRollingShutter = True
    
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
    
    ## Reduce Error Parameters
    RU_thresh = 10  #Start at this threshold, go up by 1 until it reaches RU_minRemainingPercent of remaining points OR RU_thresh = 50
    PA_thresh = 2  #Start at this threshold, go up by 0.25 until it reaches PA_minRemainingPercent of remaining points OR  PA_thresh = 6
    RE_thresh = 0.3 #Start at this threshold, go up by 0.05 until it reaches RE_minRemainingPercent of remaining points OR RE_thresh = 2
    # Minimum percentage of points remaining after each RE step
    RU_minRemainingPercent = 50
    PA_minRemainingPercent = 50
    RE_minRemainingPercent = 90
    removeMarkersFlag = 1  # Remove markers that have fewer than 3 projections
    
    ## Define header folder
    siteLoc = "BogusRidge"
    mainFold = "/SNOWDATA/SnowDrones/" + siteLoc + "/"
    imgExt = ".JPG"
    
    # Define where to save the photoscan project and related files (e.g. ortho)
    saveDir = "/SNOWDATA/SnowDrones_2021_Processing/" + siteLoc + "/"

    # Used to name the project file (e.g. <saveID>.psx)
    saveID = "HighQuality" # Project to open
    customBB = 1        # Set a custom bounding box around camera locations
    markerErrOutput = 0  # Save a CSV file with details on individual marker accuracies
    savePointCloud = 0   # Save a point cloud
    saveDEM = 0
    saveOrtho = 0   # Save an orthophoto
    saveReport = 0
    
    # Run the script over all folders in mainFold
    iterFold()
    
    
    # #Close the App
    # app = PhotoScan.Application()
    # app.quit()



