# Script to run Metashape and create a pointcloud from given image folder (adapted from CB)
# last edits: BG 08/05/2026

# edits needed:
# ----------------------------
# non hard coded paths
# fix exporting pointcloud
# GCPs need to be added to the photo, this still isn't working.... maybe have metashape scan for GCPs? test this out
# run code in terminal:  /Applications/MetashapePro.app/Contents/MacOS/MetashapePro -r /Users/bagaenzl/MasonBEAST/aur_meta.py
import Metashape
import os
# Monday Aug 12, 2024
image_folder=r'/Volumes/Elements/RawImages/1702830601159/' # change to be not hard coded just for test run
#gcp_path=r'/Volumes/Elements/GCPs/GCPsfromNCSUcomp/GCPS_08_12_2024.txt'
pcpath="/Volumes/Elements/1702830601159/1702830601159_ptcld"
cam_a=r'/Volumes/Elements/MetashapeFiles/camA_precal.xml'
cam_b=r'/Volumes/Elements/MetashapeFiles/camB_precal.xml'

# initialize project document
doc=Metashape.app.document
psx_path="/Volumes/Elements/MetashapeFiles/1702830601159.psx" # needs to be changed to not be hard coded
doc.save(psx_path)

print('-----------------------------------------')
print('Step 1: Create Chunk and Upload Images')
print('-----------------------------------------')
# create chunk if there are none
if len(doc.chunks):
    chunk = doc.chunk
    # WIPE THE SLATE CLEAN: Remove old cached cameras/sensors from previous failed runs
    for c in list(chunk.cameras):
        chunk.remove(c)
    for s in list(chunk.sensors):
        chunk.remove(s)
    print("Cleared existing chunk data from previous test run.")
else:
    chunk = doc.addChunk()

# split cams
imgs_a=sorted([
     os.path.join(image_folder, f) 
    for f in os.listdir(image_folder) 
    if f.lower().endswith('.tiff') and not f.startswith('.') and "cameraa" in f.lower()
])
imgs_b=sorted([
     os.path.join(image_folder, f) 
    for f in os.listdir(image_folder) 
    if f.lower().endswith('.tiff') and not f.startswith('.') and "camerab" in f.lower()
])

print(f" found {len(imgs_a)} images for Cam A and {len(imgs_b)} for Cam B")

chunk.addPhotos(imgs_a, layout=Metashape.MultiframeLayout)
chunk.addPhotos(imgs_b,layout=Metashape.MultiframeLayout)

# pair images
#paired_imgs=[list(pair) for pair in zip(imgs_a,imgs_b)]
#flat_alt_imgs=[]
#file_groups_map=[]
#for frame_idx,(a,b) in enumerate(zip(imgs_a,imgs_b)):
    #flat_alt_imgs.extend([a,b])
    #file_groups_map.extend([frame_idx,frame_idx])

# import images
#chunk.addPhotos(imgs_a,layout=Metashape.MultiframeLayout,strip_extensions=True)
#chunk.addPhotos(imgs_b,layout=Metashape.MultiframeLayout,strip_extensions=True)
#chunk.addPhotos(filenames=flat_alt_imgs,filegroups=file_groups_map,layout=Metashape.MultiframeLayout)

print('--------------------------------------------')
print('Step 2: Separate Cameras and Import Intrinsics/Extrinsics')
print('--------------------------------------------')

doc.save()

sensor_cam_a=chunk.sensors[0] if len(chunk.sensors)>0 else chunk.addSensor()
sensor_cam_b=chunk.addSensor()

sensor_cam_b.width=sensor_cam_a.width
sensor_cam_b.height=sensor_cam_a.height

# import cam intrinsics
sensor_cam_a.label="Camera A"
sensor_cam_a.type=Metashape.Sensor.Type.Frame
sensor_cam_a.calibration.load(cam_a,format=Metashape.CalibrationFormatXML)
sensor_cam_a.fixed_calibration=True #fixed imported values

sensor_cam_b.label="Camera B"
sensor_cam_b.type=Metashape.Sensor.Type.Frame
sensor_cam_b.calibration.load(cam_b,format=Metashape.CalibrationFormatXML)
sensor_cam_b.fixed_calibration=True #fixed imported values

gcp_crs=Metashape.CoordinateSystem("EPSG::6347")
chunk.crs=gcp_crs
chunk.camera_crs=gcp_crs

pos_camA=Metashape.Vector([239766.4,3784761.5,10.324904]) # measured cam extrinsics from drive
pos_camB=Metashape.Vector([239759.3,3784754.9,10.392904])

for frame in chunk.frames:
    camera_a=frame.cameras[0]
    camera_a.sensor=sensor_cam_a
    camera_a.label="Camera A"
    camera_a.reference.location=pos_camA
    camera_a.reference.enabled=True

    camera_b=frame.cameras[1]
    camera_b.sensor=sensor_cam_b
    camera_b.label="Camera B"
    camera_b.reference.location=pos_camB
    camera_b.reference.enabled=True

doc.save()

print("loaded images and calibrated cameras")

print('----------------------------------')
print('Step 3: Align Photos')
print('-------------------------------------')

# allign photos
print("begining allignment")
count=0 
chunk=Metashape.app.document.chunk
#chunk.matchPhotos(downscale=1, #accuracy
                      #generic_preselection=True,
                      #reference_preselection=False,
                      #filter_mask=False,
                      #mask_tiepoints=False,
                      #keypoint_limit=60000,
                      #tiepoint_limit=0,
                      #guided_matching=True,
                      #keep_keypoints=True
                      #)
for frame in chunk.frames:
    count+=1
    frame.matchPhotos(downscale=1, #accuracy
                      generic_preselection=True,
                      reference_preselection=False,
                      filter_mask=False,
                      mask_tiepoints=False,
                      keypoint_limit=60000,
                      tiepoint_limit=0,
                      keep_keypoints=True
                      )
    print('matching frame number:',str(count))

# Accuracy options 0 = ultra high (upscaled 2x), 1= highest accuracy (original size of image), 2=medium (downscaled 2x), 4 = lowest (downscaled 4x)

print('photo matching complete... aligning cams')
chunk.alignCameras(adaptive_fitting=False)
aligned_count = sum(1 for c in chunk.cameras if c.transform)
print(f"Cameras successfully aligned: {aligned_count} out of {len(chunk.cameras)}")

if aligned_count < 3:
    print("CRITICAL ERROR: Not enough cameras aligned to build a 3D dense cloud.")
    print("Please check if your images have enough overlap or if they are blurry/dark.")
    #sys.exit(1)

print("Cameras aligned successfully.")
doc.save()

print('---------------------------------------')
print('Step 4: Bring in GCPs and Set Coordinate System')
print('----------------------------------------')
doc.save()
# This is not working yet so needs to be manually done
# scan images to detect GCPs (refine later)?
#print("detecting markers.............")
#chunk.detectMarkers(target_type=Metashape.CrossTarget)
#import csv
#gcp_reference_coords={} # dictionary of gcp locations (measured in field)
#with open(gcp_path, mode='r') as f:
    #reader=csv.reader(f,delimmiter=" ")
    #for row in reader: 
        #if not row or row[0].startswith("#"):
            #continue # skip any empty lines or comments (usually there aren't any)
        #label=row[0].strip()
        #x=float(row[1]) # easting
        #y=float(row[2]) # northing
        #z=float(row[3]) # altitude

        #gcp_reference_coords[label]=[x,y,z]
# apply coords to ID'ed markers
#for marker in chunk.markers:
    #if marker.label in gcp_reference_coords:
        #coords=gcp_reference_coords[marker.label]
        #marker.reference.location=Metashape.Vector(coords)
        #marker.reference.enabled=True
        #print(f"assigned coords to {marker.label}")
    #else:
        #print(f"no coordinates found for marker {marker.label}")

# update chunk and optimize
#chunk.updateTransform()
#chunk.optimizeCameras()
#print("gcps imported and alignment optimized")
# bring in GCPs
#gcp_crs=Metashape.CoordinateSystem("EPSG::6347")
#chunk.crs=gcp_crs
#chunk.marker_crs=gcp_crs

#chunk.importReference(path=gcp_path,
                      #format=Metashape.ReferenceFormatCSV,
                      #columns="nxyz",
                      #delimiter=" ",
                      #group_delimiters=True,
                      #crs=gcp_crs,
                      #items=Metashape.ReferenceItemsMarkers,
                      #create_markers=True)
#chunk.updateTransform() # update based on GCPs
#doc.save()

#T=chunk.transform.matrix

#for marker in chunk.markers:
#     pinned_ct=0
#     print(f"total markers found {len(chunk.markers)}")
#     if not marker.reference.location:
#         print(f"skipping {marker.label}: empty ref loc")
#         continue
#     #transform gcps from world to metashape geocentric coors
#     gcp_geo=chunk.crs.unproject(marker.reference.location)
#     # transform from geocentric to internal 3D chunk local coords
#     gcp_internal=T.inv().mulp(gcp_geo)
#     #marker.position=gcp_internal #3D anchor position

#     for frame in chunk.frames:
#         for camera in chunk.cameras:
#             if not camera.transform:
#                 continue
#         # convert 3d chunk to local cam coord
#             #cam_pt=camera.transform.inv().mulp(gcp_internal)
#             #if cam_pt.z<=0: # in case behind cameras
#                 #continue
#         #project 3d coord onto 2d image sensor aka pixel space
#             #pixel=camera.sensor.calibration.project(cam_pt)
#             pixel=camera.project(gcp_internal)
#             print(f"Projected Pixel for {marker.label} on {camera.label}: X={pixel.x}, Y={pixel.y}") 
#         # put marker on gcp
#             if 0<= pixel.x < camera.sensor.width and 0<= pixel.y<camera.sensor.height:
#                 pixel_vector=Metashape.Vector([pixel.x,pixel.y])
#                 marker.projections[camera]=Metashape.Marker.Projection(pixel_vector,pinned=True)
#                 pinned_ct+=1
#                 print(f"pinned gcp {marker.label} onto image frame {camera.label} at pixel x: {round(pixel.x,1)}, y: {round(pixel.y,1)}")

# print(f'pinned {pinned_ct} gcps onto images')
#chunk.optimizeCameras()
#doc.save()
#print('gcps added and optimized...building pointclouds')

#chunk.buildDepthMaps(downscale=1,
                     #filter_mode=Metashape.ModerateFiltering,
                     #) # 2=high quality, moderate depth filtering)
#doc.save()

#chunk.buildPointCloud(source_data=Metashape.DepthMapsData,
                      #point_colors=True,
                      #point_confidence=True,
                      #keep_depth=True)
##doc.save()
#print('pointclouds complete...exporting pointclouds')
#chunk.exportPointCloud(path=pcpath,
                       #source_data=Metashape.PointCloudData,
                       #format=Metashape.PointCloudFormatXYZ,
                       #crs=chunk.crs,
                       #save_point_color=True,
                       #precision=6)

#print("pointclouds exported")