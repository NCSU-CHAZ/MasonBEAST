import Metashape
import os

# initialize project document
doc=Metashape.app.document
psx_path="/Volumes/Elements/MetashapeFiles/1723489201189.psx" # needs to be changed to not be hard coded
doc.save(psx_path)
chunk=doc.addChunk()

# import cam intrinsics
cam_a=r'/Volumes/Elements/MetashapeFiles/camA_precal.xml'
cam_b=r'/Volumes/Elements/MetashapeFiles/camB_precal.xml'

# split cams
imgs_a=[]
imgs_b=[]
# Monday Aug 12, 2024
image_folder=r'/Volumes/Elements/RawImages/1723489201189/' # change to be not hard coded just for test run
for filename in os.listdir(image_folder):
    if filename.lower().endswith(('.tiff')):
        full_path=os.path.join(image_folder,filename)

        if "CameraA" in filename.lower():
            imgs_a.append(full_path)
        elif "CameraB" in filename.lower():
            imgs_b.append(full_path)
        else:
            print(f" file {filename} isn't an image")

print(f" found {len(imgs_a)} images for Cam A and {len(imgs_b)} for Cam B")

# import images
chunk.addPhotos(imgs_a)
cams_group_a=[c for c in chunk.cameras if c.photo.path in imgs_a]

chunk.addPhotos(imgs_b)
cams_group_b=[c for c in chunk.cameras if c.photo.path in imgs_b]

# import cam intrinsics
sensor_cam_a=chunk.addSensor()
sensor_cam_a.label="Camera A"
sensor_cam_a.type=Metashape.Sensor.Type.Frame
sensor_cam_a.calibration.load(cam_a,format=Metashape.CalibrationFormat.XML)
sensor_cam_a.fixed_calibration=True #fixed imported values

for camera in cams_group_a:
    camera.sensor=sensor_cam_a

sensor_cam_b=chunk.addSensor()
sensor_cam_b.label="Camera B"
sensor_cam_b.type=Metashape.Sensor.Type.Frame
sensor_cam_b.calibration.load(cam_b,format=Metashape.CalibrationFormat.XML)
sensor_cam_b.fixed_calibration=True #fixed imported values

for camera in cams_group_b:
    camera.sensor=sensor_cam_b

# clear any default sensors (recomended)
for s in list(chunk.sensors):
    assigned_cams=[c for c in chunk.camera if c.sensor==s]
    if len(assigned_cams)==0:
        chunk.remove(s)

print("loaded images and calibrated cameras")

# allign photos
print("begining allignment")
chunk.matchPhotos(generic_preselection=True,reference_preselection=False,tiepoint_limit=60000,downscale=1) # downscale =1 is high accuracy
print('photo matching complete... aligning cams')
chunk.alignCameras()
print("cams alligned")

# bring in GCPs
gcp_crs=Metashape.CoordinateSystem("ESPG::6347")
chunk.crs=gcp_crs
chunk.marker_crs=gcp_crs
gcp_path=r'/Volumes/Elements/GCPs/GCPsfromNCSUcomp/GCPS_08_12_2024.txt'
chunk.importReference(path=gcp_path,format=Metashape.ReferenceFormatCSV,columns="nxyz",delimiter=" ",crs=gcp_crs,items=Metashape.ReferenceItemsMarkers,create_markers=True)

chunk.updateTransform()
chunk.optimizeCameras()
print('gcps added and optimized...building pointclouds')

chunk.buildDepthMaps(downscale=2,filter_mode=Metashape.ModerateFiltering) # 2=high quality, moderate depth filtering)
doc.save()

chunk.buildPointCloud(source_data=Metashape.DepthMapsData,point_colors=True,point_confidence=True,keep_depth=True)
doc.save()
print('pointclouds complete...exporting pointclouds')
pcpath="/Volumes/Elements/1723489201189/1723489201189_ptcld"
chunk.exportPointCloud(path=pcpath,source_data=Metashape.PointCloudData,format=Metashape.PointCloudFormatTXT,crs=chunk.crs,save_colors=True,precision=6)

print("pointclouds exported")