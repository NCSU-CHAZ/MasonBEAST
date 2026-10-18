#export_ptcld.py

# Script to execute pointcloud export
# path to run in terminal:
# /Applications/MetashapePro.app/Contents/MacOS/MetashapePro -r /Users/bagaenzl/MasonBEAST/export_ptcld.py
import os
import Metashape
# paths
epochnum='1702841401078'
genpath=r'/Volumes/Elements/'
psxpath=r'/Volumes/Elements/MetashapeFiles/1702841401078.psx'
savepath=r'/Volumes/rsstu/users/k/kanarde/MasonBEAST/data/PointClouds/1702841401078/'#'/Volumes/Elements/PointClouds/1702841401078/'

# function to export
def export_ptcld_func(chunk,savepath):
    if not chunk.point_cloud:
        print(f"skip chunk: '{chunk.label}. Does not have any frames")
        return
    
    print(f"Chunk: {chunk.label} has {len(chunk.frames)} frames. Starting export...")

    os.makedirs(savepath,exist_ok=True)

    for frame_idx,frame in enumerate(chunk.frames):
        cloud=frame.point_cloud

        if not cloud:
            print(f"Frame {frame_idx}: Null, skipping...")
            continue
        if cloud.point_count==0:
            print(f" Frame {frame_idx}: Pointcloud is empty (0) points. Skipping.....")
            continue
        print(f"Frame {frame_idx}: Found {cloud.point_count} points. Exporting....")

        filename=os.path.join(savepath,f'{epochnum}_ptcld{frame_idx:04d}.txt')

        try:
            chunk.exportPointCloud(path=filename,
                               source_data=Metashape.PointCloudData,
                               format=Metashape.PointCloudFormatXYZ,
                               crs=chunk.crs,
                               save_point_normal=True,
                               save_point_color=True
                               )
            print(f" Frame {frame_idx} successfully exported")
        except Exception as e:
            print(f"Frame {frame_idx} failed to export: {e}")


doc=Metashape.Document()
doc.open(psxpath,ignore_lock=True)

for chunk in doc.chunks:
    export_ptcld_func(chunk,savepath)
    print(f"Pointclouds exported to {savepath}")


