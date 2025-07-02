% GEM region density analysis
% BG 6/26/25

% EDITS NEEDED
% ---------------
% GEM names are plotted in reverse order, need to fix
% geomorphic region will change for each GEM - need to figure out how to
% implement in function

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
path=metashape_path; % or measured_path
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
    % save GEM name
    GEMname=split(GEMpath,'/');
    GEMname=GEMname(9,1);
    GEMnames(k)=string(GEMname);
    for j=1:length(regions)
        region=regions(j);
        % calculate density and save for each GEM
        density(k,j)=GEM_region_density_calc(GEMpath,region,figpath);
    end
end

% plot density histograms
% --------------------------
% add  row and column of zeros for plotting
density_rbr=density(:,4); % grab RBR density
density_norbr=density(:,[1:3,5]); % just densities of geomorphic locations

zc=zeros(size(density_rbr,1),1);
density_wzeros=[density_rbr,zc];
zr=zeros(1,size(density_wzeros,2));
densityrbr_wzeros=[zr;density_wzeros];
regions_rbr={'RBR'}; % regions w/ rbr for plotting

% plot rbr density
[x,y]=meshgrid(1:2,1:length(epochstring)+1);
z=densityrbr_wzeros;
fig=figure(1);
clf; pcolor(x,y,z); hold on; 
a=colorbar; clim([0 1]);a.Label.String='Normalized Density (-)';
xticks([0.5 1.5]); yticks([0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5]); % need to generalize this
yticklabels([0,flip(GEMnames)]);xticklabels([0,regions_rbr]); title('GEM RBR Density');fontsize(gcf,16,"points");


% plot geomorphic region density
zc=zeros(size(density_norbr,1),1);
density_wzeros=[density_norbr,zc];
zr=zeros(1,size(density_wzeros,2));
density_norbr_wzeros=[zr;density_wzeros];

regions_no_rbr={'dune','upperbeachface','lowerbeachface','shoreline'}; % regions w/o rbr for plotting

[x,y]=meshgrid(1:length(regions_no_rbr)+1,1:length(epochstring)+1);
z=density_norbr_wzeros;
fig=figure(2);
clf; pcolor(x,y,z); hold on; 
a=colorbar; clim([0 1]);a.Label.String='Normalized Density (-)';
xticks([0.5 1.5 2.5 3.5 4.5]); yticks([0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5]); % need to generalize this
yticklabels([0,flip(GEMnames)]);xticklabels([0,regions_no_rbr]); title('GEM Region Density');fontsize(gcf,16,"points");

figpath=fullfile(figpath,'GEMdensityplot');
saveas(fig,figpath,'png');
close(fig);

