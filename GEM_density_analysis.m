% GEM region density analysis
% BG 6/25/25

% TESTING CODE ------------
% Paths
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % metashape GEMs
GEMpath_meta=append(metashape_path,'/1723489201189/');
GEMpath_meas=append(measured_path,'/1723489201189/');

listofFiles=dir(measured_path);
filenames = cell(length(listofFiles), 1); % Preallocate cell arra

for i=1:length(listofFiles)
    baseFilename=listofFiles(i).name;
    fullFilename=fullfile(listofFiles(i).folder,baseFilename);
    filenames{i}=fullFilename;
end

% grab only GEM files (have at least one number in name)
epochfile=regexp(filenames,'\d','once'); 
epochfile=~cellfun('isempty',epochfile);
epochstring=filenames(epochfile);

for k=1:length(epochstring)
    % load GEM mat file
    GEMpath=string(append(epochstring(k),'/meanGEMz.mat'));
    GEMz=load(GEMpath);
    GEMz=GEMz.meanGEMz;
    
end



figpath=append(measured_path,'Figures');