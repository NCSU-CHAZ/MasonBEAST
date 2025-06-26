% GEM region density analysis
% BG 6/26/25

% Paths
% ---------------
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % metashape GEMs

% Constants
% ---------------
regions={'dune','upperbeachface','lowerbeachface','RBR','shoreline'}; % regions to calculate density

% Process GEMs
% ----------------
path=measured_path; % or metashape_path
% list of files in folder
listofFiles=dir(path);
filenames = cell(length(listofFiles), 1); % Preallocate cell array
% grab file names 
for i=1:length(listofFiles)
    baseFilename=listofFiles(i).name;
    fullFilename=fullfile(listofFiles(i).folder,baseFilename);
    filenames{i}=fullFilename;
end

% grab only GEM files (have at least one number in name)
epochfile=regexp(filenames,'\d','once'); 
epochfile=~cellfun('isempty',epochfile);
epochstring=filenames(epochfile);

density=zeros(length(epochstring),length(regions)); % preallocate density matrix

for k=1:length(epochstring)
    % load GEM mat file
    GEMpath=string(epochstring(k));
    figpath=append(path,'Figures'); 
    for j=1:length(regions)
        region=regions(j);
        % calculate density and save for each GEM
        density(i,j)=GEM_region_density_calc(GEMpath,region,figpath);
    end
end

