%% Script to process monthly GEMs (camera location comparison)

% paths 
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage 
pcpath=append(genpath,'/PointClouds/'); % path to pointclouds
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % save measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % save metashape GEMs

% GEM specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
dxy = 0.5; % meter
numframes=2;

% Loop through Metashape pointclouds
spec='*meta_ptcld*';
figpath=append(metashape_path,'Figures');
[meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,dxy,numframes,pcpath,metashape_path,figpath,spec);

% compare to hand surveys (need to figure out how to grab surveys through
% loop) 
HS_savepath=append(genpath,'/Surveys/rotated_surveys');
surveypath=append(CAM_analysispath,'/Surveys');
listofsurveys=dir(surveypath);
for i=2:length(listofsurveys)
    listofsurvs(i)=append(listofsurveys(i).folder,listofsurveys(i).name);
end

listofFiles=dir(metashape_path);
for i=1:length(listofFiles)
    GEMmatpath=append(listofFiles(i).folder,'/',listofFiles(i).name);
    [handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
end

% Loop through Measured pointclouds
spec='*meas_ptcld*';
figpath=append(measured_path,'Figures');
GEMsavepath=measured_path;
[meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,dxy,numframes,pcpath,measured_path,figpath,spec);

