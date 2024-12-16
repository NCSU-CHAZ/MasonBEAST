%% Code to grab just a transect from GEM 
% BG and CB
% reads in point cloud and grabs transect ( a shortened version of
% proc_ptcloud.m)
clc
clear all
close all
format long g

addpath(genpath('/Users/bagaenzl/Desktop/CIRN-Quantitative-Coastal-Imaging-Toolbox'))
%addpath(genpath('/Users/cmbaker9/Documents/Programs/MTOOLS'))
%addpath('/Users/cmbaker9/Documents/Research/MasonBEAST/code/imagery/support_routines')

%% STEP 1: Create paths, files and naming

% general path and names
datapath    = '/Users/bagaenzl/Desktop/MasonBEAST Data/';
trialname   = '1730736001873';%'1702827001820';%'1699459201959';%'1702827001820';%'1698951602001';
tempname   = '1730736001873';
%list of epoch nums using
% 1717167601959 = May 31 11 AM 2024
% 1719428401397 = June 26 3PM 2024
% 1721934001780 = July 25 3PM 2024
% '1723489201189' = Aug 12 2024 x
% 1726056001556 = Sept 11 8 AM 2024
% '1726412401052' = Sept 15 11AM 2024 x
% '1726772401770' = Sept 19 3PM 2024 x
% 1726858801716 = Sept 20 3PM 2024
% '1728327601158' = Oct 7 3PM 2024 x
% 1728388801457 = Oct 8 8AM 2024
% '1728561601419' = Oct 10 8AM 2024 x
% '1730736001873' = Nov 4 11PM 2024 x
% 1732464001262 = Nov 24 11 AM 2024

%gcpname = [datapath,'GCP Locations txts/gcps/GCPs_09_18_2024.txt']; 
% (don't need for this)

%orthopath = [datapath,'stereo/'];
pcpath = [datapath,'PointClouds/',trialname]; % need to make this a bit more consise and intuitive
%numframes = length(dir(fullfile(pcpath, '*.txt')));
numframes=2;
%% STEP 2: Define variables

% define local coordinate system origin and rotation angle
% Xloc = min(ptcl(:,1));
% Yloc = max(ptcl(:,2));
Xloc = 239737;
Yloc = 3784751;
rotang = 35; % i'm a little bit confused on where these numbers came from. Is the rotation angle the rough angle of the alongshore?

% create grid
dxy = 0.2;
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[x,y] = meshgrid(gridX,gridY);

% create matrix for storing z values
z = NaN(size(x,1),size(x,2),numframes);
numpts = z;

% Creating transects
ypick = [7.3 -20]; % in m, define locations to pick transects
% Also best place to grab this?

yavg = 1.2; % in m, defining how wide in the alongshore to average over
% This is 1.2 m. Do we want to do that? How should we pick the width value?
% Should I redo it with varying width and see results?

iyavg = round(yavg/dxy/2); % indices to pull out before and after transect to average over

%% STEP 3: Read and rotate GCPs

% GCP Input
%GCPs = readmatrix(gcpname);
%GCPz = GCPs(:,4); % I seperate the x, y and z but is not required

% gcp rotation and transformation
%[GCPx, GCPy] = rotateCoordinates(GCPs(:,2), GCPs(:,3), Xloc, Yloc, rotang); 

% Cam x,y
%CamA = [239766.1, 3784761.9];%, 10.37
%CamB = [239759.4, 3784755.0];%, 10.26];

% camera rotation and transformation
%[CamAx, CamAy] = rotateCoordinates(CamA(1), CamA(2), Xloc, Yloc, rotang);
%[CamBx, CamBy] = rotateCoordinates(CamB(1), CamB(2), Xloc, Yloc, rotang);

%% STEP 4: Read data and grid data

for i = 1:numframes

    % % read dem
    % [A,R] = readgeoraster(fname); % read in tif file (flips file over and you have to flip it again later)
    % A(A<-7) = NaN; % Create reasonable z limit
    % Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
    % Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
    % [X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
    % A = flipud(A); % When tif files are read in this way they need to be flipped

    % read point cloud
    if i == 25
        ptcl = [0 0 0];
    else
        ptcl = readmatrix([pcpath,'/',tempname,'_ptcld',num2str(i),'.txt']); % columns x,y,z (This file path is not working FIGURE OUT WHY!!) Possibly rewrite code with own file naming on seperate file.
    end

    % % read point cloud
    % [ortho, orthoinfo] = imread([orthopath,trialname,'_orthomosaic',num2str(i),'.tif']); % columns x,y,z

    % % Dem coordinate rotation and transformation
    % [Xrot, Yrot] = rotateCoordinates(X, Y, Xloc, Yloc, rotang);

    % pt cloud coordinate rotation and transformation x,y to cross- and alongshore
    [Xrot, Yrot] = rotateCoordinates(ptcl(:,1), ptcl(:,2), Xloc, Yloc, rotang);

    % grid point cloud
    [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),x,y,@median); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
    ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
    ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan

    % store output into a matrix
    z(:,:,i) = ztemp;
    numpts(:,:,i) = ntemp;
    clear ztemp ntemp

end
%% STEP 6: Extract a cross-shore transect

% extract transect at location ypick(1)
[~,iy] = min(abs(y(:,1)-ypick(1)));
ztran = median(z(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran);
ztran = movmean(ztran,2,1,'omitnan'); % save ztran


%% Try creating transects by looping through Pointcloud folder
folderpath='/Users/bagaenzl/Desktop/MasonBEAST Data/PointClouds/';
filePattern=fullfile(folderpath,'*.txt');
listofFiles=dir(filePattern);
% Need to define variables still
% define local coordinate system origin and rotation angle
% Xloc = min(ptcl(:,1));
% Yloc = max(ptcl(:,2));
Xloc = 239737;
Yloc = 3784751;
rotang = 35; 

% create grid
dxy = 0.2;
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[x,y] = meshgrid(gridX,gridY);

% create matrix for storing z values
numframes=2;
z = NaN(size(x,1),size(x,2),numframes);
numpts = z;

% Creating transects
ypick = [7.3 -20]; % in m, define locations to pick transects
% Also best place to grab this?

yavg = 1.2; % in m, defining how wide in the alongshore to average over
% This is 1.2 m. Do we want to do that? How should we pick the width value?
% Should I redo it with varying width and see results?

iyavg = round(yavg/dxy/2); % indices to pull out before and after transect to average over


% loop through files in PointCloud folder
for k=1:length(listofFiles)
    baseFilename=listofFiles(k).name;
    fullFilename=fullfile(listofFiles(k).folder,baseFilename);
    epochname=split(baseFilename,'.');
    epochname=string(epochname(1));
    transectFilename=append(epochname,'tran');
    %transectname=[baseFilename,'transect'];
    for i = 1:numframes
    % read point cloud
        if i == 25
            ptcl = [0 0 0];
        else
            ptcl = readmatrix(fullFilename); % columns x,y,z
        end
    % pt cloud coordinate rotation and transformation x,y to cross- and alongshore
        [Xrot, Yrot] = rotateCoordinates(ptcl(:,1), ptcl(:,2), Xloc, Yloc, rotang);

    % grid point cloud
        [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),x,y,@mean); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
        ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
        ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan

        % store output into a matrix
        z(:,:,i) = ztemp;
        numpts(:,:,i) = ntemp;
        clear ztemp ntemp
    end
    % extract transect at location ypick(1)
    [~,iy] = min(abs(y(:,1)-ypick(1)));
    ztran = median(z(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
    ztran = squeeze(ztran);
    ztran = movmean(ztran,2,1,'omitnan'); 
    save(transectFilename,"ztran",'-mat') % save ztran
end

%% Loop through to create graph of transects 
% not quite sure how to do this yet tbh so here is the harder way to do it

% load transects and put in array
tranfolderpath='/Users/bagaenzl/Desktop/MasonBEAST Data/MasonBEAST';
tranfilePattern=fullfile(tranfolderpath,'*tran');
tranlistofFiles=dir(tranfilePattern);


for k=1:length(tranlistofFiles)
    tranbaseFilename=tranlistofFiles(k).name;
    tranfullFilename=fullfile(tranlistofFiles(k).folder,tranbaseFilename);
    Ztran_struct=load(tranfullFilename);
    for i=1:length(Ztran_struct)
        Ztran=Ztran_struct.ztran(:,1);
    end
    ztran_array(:,k)=Ztran;
end

% plotting specs
x=0:0.15:(550*0.15); % x meters across the beach

ftsz = [24 20];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
% iylim = [-13,13];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 1;
xlab = 'Cross-Shore (m)';
ylab = 'Elev. (m)';

ax1pos = [0.1 0.2 0.6 0.6];
% I want to plot these
% '1723489201189' = Aug 12 2024 x k 5
% 1726056001556 = Sept 11 8 AM 2024 k 6
% '1726412401052' = Sept 15 11AM 2024 x 7
% '1726772401770' = Sept 19 3PM 2024 x k 8
% 1726858801716 = Sept 20 3PM 2024 9
% '1728327601158' = Oct 7 3PM 2024 x 10
% 1728388801457 = Oct 8 8AM 2024  11
% '1728561601419' = Oct 10 8AM 2024 x 12 k
% '1730736001873' = Nov 4 11PM 2024 x 13 k
% 1732464001262 = Nov 24 11 AM 2024 14 k 
% read image and process it
margins = [7 7 15 11]; % N S E W
A = imread('https://www.mathworks.com/matlabcentral/answers/uploaded_files/853000/image.png');
A = A(1+margins(1):end-margins(2),1+margins(4):end-margins(3),:);
CT0 = permute(mean(im2double(A),1),[2 3 1]);
CT0 = CT0([true; ~all(diff(CT0,1,1)==0,2)],:); % remove duplicate rows
N = 10; % specify the number of colors in table
na = size(CT0,1);
CT = interp1(linspace(0,1,na),CT0,linspace(0,1,N));
figure('units','inches','position',[1 1 12 6],'color','w');
ax1 = axes('Position',ax1pos);
hold on
for i=1:69
    if i==37||i==40 
        plot(x,ztran_array(:,i),"Linewidth",3.5)
        hold on
    elseif i==47||i==50
        plot(x,ztran_array(:,i),"Linewidth",3.5)
        hold on
    elseif i==56
        plot(x,ztran_array(:,i),"Linewidth",3.5)
        hold on
    end
end
colororder(CT);
box on
ylim([0 4.5]);
ylabel(ylab);
xlabel(xlab);
xlim([0 40]);
clim([0 11]);
legend('Sept 11','Sept 19','Sept 20');
% set(ax2,'YDir','normal','Color','k')
colormap("spring");

% plot TC8
figure('units','inches','position',[1 1 12 6],'color','w');
ax1 = axes('Position',ax1pos);
hold on
plot(x,ztran_array(:,37),"Linewidth",3.5); hold on; plot(x,ztran_array(:,47),"Linewidth",3.5); hold on;
plot(x,ztran_array(:,50),"Linewidth",3.5); hold on; plot(x,ztran_array(:,56),"Linewidth",3.5); hold on; 
colororder(CT);
box on
ylim([0 4.5]);
ylabel(ylab);
xlabel(xlab);
xlim([0 40]);
clim([0 11]);
legend('Sept 11','Sept 19','Sept 20','Oct 8');
title('TC8 Transects');