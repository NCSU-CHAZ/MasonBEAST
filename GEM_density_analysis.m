% GEM region density analysis
% BG 07/28/2025

% EDITS NEEDED
% ---------------
% fix big in GEM name labeling in density plot
% plot RBR density
% fix bug in rbr density calculation

% Paths
% ---------------
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % metashape GEMs

% Constants
% ---------------
regions={'dune','upperbeachface','lowerbeachface','RBR','shoreline'}; % regions to calculate density
% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;
% create grid
dxy=0.5;
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
rbrloc = [239779.25, 3784738.325]; % location of RBR in swash
% rbr rotation and transformation
[rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);
% camera rotation and transformation
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
[CamAx, CamAy] = rotateCoordinates(camlocA(1), camlocA(2), Xloc, Yloc, rotang);
[CamBx, CamBy] = rotateCoordinates(camlocB(1), camlocB(2), Xloc, Yloc, rotang);


% Process GEMs
% ----------------
path=metashape_path; % or measured_path metashape_path
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

% Preallocate density array for density
density=zeros(length(regions),1); % preallocate density matrix

%for k=1:length(regions)
k=14; % change for each GEM wanted
    % load GEM mat file
    GEMpath=string(epochstring(k));
    figpath=append(path,'Figures');
    % save GEM name
    GEMname=split(GEMpath,'/');
    GEMname=GEMname(9,1);
    GEMnames=string(GEMname);
    GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
    GEMdate(k)=string(GEMdate(1,1));

    meanGEMz=load(fullfile(GEMpath,'meanGEMz_1.mat'));
    meanGEMz=[meanGEMz.meanGEMz];

    % Plot GEM to ID regions
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanGEMz(:,:,1)); grid off;%box on;hold on
    %scatter(GCPx,GCPy,60,'fill','r','MarkerEdgeColor','k') (need to figure
    %out how to specifiy GCPs/do we need?)
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    scatter(rbrx,rbry,60,'fill','sq','g','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'g');
    
    % calculate density for regions
    for j=1:length(regions)
        region=regions(j);
        if  strcmp(region, 'dune')
            cols=1:21;
            rows=1:10;
            rect_x=0; % lower left corner
            rect_y=0; % lower left corner
            width=9; % 4.5 m
            height=20; % 10 m
            region_specs={cols,rows,rect_x,rect_y,width,height};
        elseif strcmp(region, 'upperbeachface')
            cols=10:45;
            rows=7:70;
            rect_x=7; % lower left corner
            rect_y=-10; % lower left corner
            width=70-7; % (divide by 2) m
            height=35; 
            region_specs={cols,rows,rect_x,rect_y,width,height};
        elseif strcmp(region, 'lowerbeachface')
            cols=-80:10;
            rows=7:70;
            rect_x=7; % lower left corner
            rect_y=-80; % lower left corner
            width=70-7; % 
            height=70; % 
            region_specs={cols,rows,rect_x,rect_y,width,height};
        elseif strcmp(region, 'shoreline')
            cols=-80:24;
            rows=70:75;
            rect_x=70; % lower left corner
            rect_y=-80; % lower left corner
            width=5; % 
            height=104; % 
            region_specs={cols,rows,rect_x,rect_y,width,height};
        elseif strcmp(region,'RBR')
            row=rbrx;
            col=rbry;
            rect_x=rbrx-1;
            rect_y=rbry-1;
            height=2;
            width=2;
            region_specs={cols,rows,rect_x,rect_y,width,height};
        end 
        density(j)=GEM_region_density_calc(GEMpath,region,region_specs,figpath); % dune, upper, lower, RBR, shoreline
    end
%end

% save density to pull later
matname=fullfile(GEMpath,append('Density','.mat')); % append(GEMpath,'/Density.mat') to load later
save(matname,'density')

% calculatae rbr density seperate 
% need to add in

% plot density histograms
% ------------------------------------
path=metashape_path; % or measured_path metashape_path
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

% preallocate density arrays
dune_density=zeros(length(epochstring),1);
upper_density=zeros(length(epochstring),1);
lower_density=zeros(length(epochstring),1);
shore_density=zeros(length(epochstring),1);

for j=1:length(epochstring)
    GEMpath=string(epochstring(j));
    figpath=append(path,'Figures');
    % save GEM name
    GEMname=split(GEMpath,'/');
    GEMname=GEMname(9,1);
    GEMnames(j)=string(GEMname);
    GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
    GEMdate(j)=string(GEMdate(1,1));

    % grab densities
    densitypath=append(GEMpath,'/Density.mat');
    density=load(densitypath);
    density=density.density;

    dune_density(j)=density(1,1);
    upper_density(j)=density(2,1);
    lower_density(j)=density(3,1);
    shore_density(j)=density(5,1);
end

% array of densities with zeros for plotting
zc=zeros(size(shore_density,1),1);
density_wzeros=[dune_density,upper_density,lower_density,shore_density,zc];
zr=zeros(1,5);
density_wzeros=[density_wzeros;zr];
% bug here in naming
for i=1:length(GEMnames)
GEMdates(i)=datetime(str2num(GEMnames(i)),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMchar=char(GEMdates(i));
GEMmonth=split(GEMchar,'-');
GEMmonth(i)=GEMmonth(2,1);
end 

% plot 
[x,y]=meshgrid(1:5,1:15);
z=density_wzeros;
regions={'Dune','UpperBeachFace','LowerBeachFace','Shoreline'};
fig=figure(3);
clf; pcolor(x,y,z); hold on; colormap('summer');
a=colorbar; clim([0 0.5]);a.Label.String='Normalized Density (-)';
xticks([0.5 1.5 2.5 3.5 4.5]); yticks([0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5]); % need to generalize this
xticklabels([0,regions]); title('GEM Region Density');fontsize(gcf,16,"points"); %yticklabels([0,GEMdates])

figpath=fullfile(figpath,'GEMdensityplot');
saveas(fig,figpath,'png');
close(fig);

%% Plot Histograms of density
fig=figure(8);
subplot(4,1,1); hist(dune_density); title("Dune"); ylabel( "# of GEMs"); xlim([0,1]);ylim([0,6]);
subplot(4,1,2); hist(upper_density); title("Upper Beach Face"); ylabel( "# of GEMs"); xlim([0,1]);ylim([0,6]);
subplot(4,1,3); hist(lower_density); title("Lower Beach Face"); ylabel( "# of GEMs"); xlim([0,1]); ylim([0,6]);
subplot(4,1,4); hist(shore_density); title("Shore");ylabel( "# of GEMs");ylim([0,6]);
xlabel("Density"); xlim([0,1]);
% OLD _____________________________________
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

