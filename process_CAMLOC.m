%% Script to process monthly GEMs (camera location comparison)

% paths 
genpath='/Volumes/kanarde-1/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
pcpath=append(genpath,'/PointClouds/'); % path to pointclouds
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % save measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % save metashape GEMs

% GEM specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.2; % meter
numframes=2;

% Loop through Metashape pointclouds
spec='*meta_ptcld*';
figpath=append(metashape_path,'Figures');
[meanGEMz,medGEMz,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,numframes,pcpath,metashape_path,figpath,spec);


% compare to hand surveys (need to figure out how to grab surveys through
% loop) 
%HS_savepath=append(genpath,'/Surveys/rotated_surveys');
%surveypath=append(CAM_analysispath,'/Surveys');
%listofsurveys=dir(surveypath);
%for i=2:length(listofsurveys)
    %listofsurvs(i)=append(listofsurveys(i).folder,listofsurveys(i).name);
%end

%listofFiles=dir(metashape_path);
%for i=1:length(listofFiles)
    %GEMmatpath=append(listofFiles(i).folder,'/',listofFiles(i).name);
    %[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
%end

% Loop through Measured pointclouds
spec='*meas_ptcld*';
figpath=append(measured_path,'Figures');
[meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,numframes,pcpath,measured_path,figpath,spec);


% GEM comparison manually
HS_savepath=append(genpath,'/Surveys/rotated_surveys');
surveypath=append(CAM_analysispath,'/Surveys');
% epoch nums: 1708030441768, 1711479601066,1713452401874, 1717156801902,
% 1719428401397,1721934001780, 1723489201189, 1726772401770, 1728561601419,
% 1730736001873, 1734181201371, 1735848001721, 1738526401587, 1741978801668
% survey names: '2024_03_26_Transects_UTM.xlsx'
% '2024_05_31_Transects_UTM.xlsx' '2024_06_26_Transects_UTM.xlsx'
% '2024_07_25_transects_UTM.xlsx' '2024_08_12_Transects_UTM.xlsx'
% '2024_09_18_Transects_UTM.xlsx' '2024_10_01_Transects_UTM.xlsx'
% '2024_11_04_Transects_UTM.xlsx' '2025_01_29_Transects_UTM.xlsx'
% '2025_03_14_Transects_UTM.xlsx'

% survey path
hand_survey_path=append(surveypath,'/2025_03_14_Transects_UTM.xlsx');
% metashape
GEMmat_path=append(metashape_path,'1741978801668/');
figpath=append(metashape_path,'Figures');
[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
% measured
GEMmat_path=append(measured_path,'1741978801668/');
figpath=append(measured_path,'Figures');
[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);

%% Pull a transect from GEMs
% specs
ypick = [7.3 -20]; % (location on grid)
yavg=1.2; % width (m)
savepath=append(metashape_path,'/Transects');
figpath=append(savepath,'/Figures');
[ztran]=transect_pull(GEMpath,numframes,yavg,dxy,ypick,savepath,figpath);

