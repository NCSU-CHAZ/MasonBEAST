%% Script to process monthly GEMs (camera location comparison)

% paths 
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage 
ptcld_path=append(genpath,'PointClouds'); % path to pointclouds
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(genpath,'/Measured'); % save measured GEMs
metashape_path=append(genpath,'/Metashape'); % save metashape GEMs

% GEM specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
dxy = 0.5; % meter
numframes=2;

% Loop through Metashape GEMs
spec='*meta_ptcld';
for i=1:length(ptcld_path)
    figpath=append(metashape_path,'/Figures');
    [meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,dxy,numframes,ptcld_path,measured_path,figpath,spec);
end