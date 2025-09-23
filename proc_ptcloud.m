% reads point cloud and creates a plots
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
trialname   = '1702827001820';%'1702827001820';%'1699459201959';%'1702827001820';%'1698951602001';
tempname   = '1702827001820';
gcpname = [datapath,'GCP Locations txts/gcps/GCPs_08_01_2024.txt'];

%orthopath = [datapath,'stereo/'];
pcpath = [datapath,'PointClouds/',trialname]; % need to make this a bit more consise and intuitive
%numframes = length(dir(fullfile(pcpath, '*.txt')));
numframes=2;

figfolder = [genpath,'/Figures'];
%% BG
genpath='/Volumes/kanarde-1/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
stormCHAZerspath=append(genpath,'/StormCHAZerz Data');
decNoreastherpath=append(stormCHAZerspath,'/Dec2023Noreaster_Processed');

epochnum='1702827001820'; % epoch number 

ptcldpath=append(decNoreastherpath,'/',epochnum); % path to pointclouds

spec='*_ptcld*';
GEMsavepath=append(ptcldpath,'/GEMs/');
figpath=append(ptcldpath,'/Figures');

numframes = 120;


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
GCPs = readmatrix(gcpname);
GCPz = GCPs(:,4); % I seperate the x, y and z but is not required

% gcp rotation and transformation
[GCPx, GCPy] = rotateCoordinates(GCPs(:,2), GCPs(:,3), Xloc, Yloc, rotang); 

% Cam x,y
CamA = [239766.1, 3784761.9];%, 10.37
CamB = [239759.4, 3784755.0];%, 10.26];

% camera rotation and transformation
[CamAx, CamAy] = rotateCoordinates(CamA(1), CamA(2), Xloc, Yloc, rotang);
[CamBx, CamBy] = rotateCoordinates(CamB(1), CamB(2), Xloc, Yloc, rotang);

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
        ptcl = readmatrix([pcpath,'/',epochnum,'_ptcld',num2str(i),'.txt']); % columns x,y,z (This file path is not working FIGURE OUT WHY!!) Possibly rewrite code with own file naming on seperate file.
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

%% STEP 5: Create a figure

scsz = 20;
ftsz = [22 18];
tklen = 0.5;
lw = 1.2;
ixlim = [-120,0];
iylim = [-100,0];
% iylim = [-13,13];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 1.5;
xlab = 'Cross-Shore (m)';
ylab = 'Alongshore (m)';
crange = [-1 4.5];
cmapplt = demcmap(crange,60);%cmocean('topo',20);

figure('units','inches','position',[0 0 10 6],'color','w')
pcolor(x,y,z(:,:,1))
box on
hold on
scatter(GCPx,GCPy,60,'fill','r','MarkerEdgeColor','k')
scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k')
scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k')
shading interp
axis equal
ylim([-60 30]);
xlim([-10 90]);
clim([0.1 3.8])
colormap(cmapplt)
clim(crange+[2 0])
h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
set(h1,'ytick',[-60:20:20]);
hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
text(98,-6,'Elev. (m)','interpreter','latex','fontsize',ftsz(1));

sname = ['GEM_',trialname];
print([figfolder,sname],'-dpng')
exportgraphics(gcf, [figfolder,sname,'.png']);
exportgraphics(gcf, [figfolder,sname,'.pdf']);

%% STEP 6: Extract a cross-shore transect and plot transects

% extract transect at location ypick(1)
[~,iy] = min(abs(y(:,1)-ypick(1)));
ztran = median(z(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran);
ztran = movmean(ztran,2,1,'omitnan'); 

ftsz = [20 16];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 0.5;
xlab = 'Cross-shore (m)';
ylab = 'Elev. (m)';
cmap = cmocean('thermal',numframes);

figure('units','inches','position',[1 1 10 3],'color','w');
hold on
for i = 1:10 %numframes
    plot(x(1,:),ztran(:,i),'LineWidth',1.5) %Color',cmap(i,:)
end
box on
clim([0 (10-1)/2]) %clim([0 (numframes-1)/2])
colormap(cmap)
h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
hc = colorbar('Location','eastoutside','Position', [0.93 0.225 0.03 0.4],'orientation','vertical','YAxisLocation','right');
set(hc,'fontsize',ftsz(2),'linewidth',lw);
text(164,3.3,'$t$ (s)','interpreter','latex','fontsize',ftsz(1));

sname = 'storm_transects';
print([figfolder,sname],'-dpng')
exportgraphics(gcf, [figfolder,sname,'.png']);
exportgraphics(gcf, [figfolder,sname,'.pdf']);

%% Extract transects for free surface and beach

% extranct transect 1
[~,iy] = min(abs(y(:,1)-ypick(1)));
ztran = median(z(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran);
ztran = movmean(ztran,2,1,'omitnan');

% compute beach as minimum (THIS CAN BE IMPROVED)
zbeach = min(ztran,[],2,'omitnan'); % beach elevation
[~,ix] = min(abs(x(1,:)-44));
zbeachmean = mean(ztran,2,'omitnan');
[~,ixmean] = min(abs(x(1,:)-12));
ix2 = find(x(1,:) > 34 & x(1,:) < 36);

[~,ixtran] = min(abs(x(1,:)-14));
[~,ixon] = min(abs(x(1,:)-2));

ztran(ztran-zbeach < 0.04)  = NaN; % NaN out z elevations if they're within 4 cm of beach elevation

% remove values that are less than 5 indeces long (removes features that
% are intermittant -- likely remove this in future. This is just to improve
% the video visualization
numstring = 5;
ztran_filt = ztran;
for i = 1:size(ztran,2)
    ztran_filt(:,i) = replaceShortNonNanSequences(squeeze(ztran(:,i)), numstring); % this is water elevation
end
 
%% STEP 6: Create a video of the transects

ftsz = [20 16];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 0.5;
xlab = 'Cross-shore (m)';
ylab = 'Elev. (m)';
cmap = cmocean('thermal',numframes);

sname =[trialname,'_timeseries'];
% I need to get Videowriter folder
v = VideoWriter([figfolder,'\',sname], 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)

figure('units','inches','position',[1 1 10 3],'color','w');

for i = 1:numframes+1
    clf
    plot(x(1,ixon:ix),zbeach(ixon:ix),'LineWidth',3,'Color',[148, 116, 27]/256)
    hold on
    plot(x(1,ix:end),ztran(ix:end,i),'LineWidth',3,'Color','b')
    %plot(x(1,ixtran:ix),ztran_filt(ixtran:ix,i),'LineWidth',3,'Color','b')

    box on
    ylim([0 4.5]);
    xlim([0 100]);
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    title(['$t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)

%% Extract transect 2 and create a plot
% What is transect 2 for?

% extract transect 2
[~,iy] = min(abs(y(:,1)-ypick(2)));
ztran2 = median(z(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran2 = squeeze(ztran2);
ztran2 = movmean(ztran2,2,1,'omitnan');

% compute beach as minimum (THIS CAN BE IMPROVED)
zbeach = min(ztran,[],2,'omitnan');
[~,ix] = min(abs(x(1,:)-44));
zbeachmean = mean(ztran,2,'omitnan');
[~,ixmean] = min(abs(x(1,:)-12));
ix2 = find(x(1,:) > 34 & x(1,:) < 36);

[~,ixtran] = min(abs(x(1,:)-14));
[~,ixon] = min(abs(x(1,:)-2));

ztran(ztran-zbeach < 0.04)  = NaN; % NaN out z elevations if they're within 4 cm of beach elevation

% remove values that are less than 5 indeces long (removes features that
% are intermittant -- likely remove this in future. This is just to improve
% the video visualization
numstring = 5;
ztran_filt = ztran;
for i = 1:size(ztran,2)
    ztran_filt(:,i) = replaceShortNonNanSequences(squeeze(ztran(:,i)), numstring);
end


%% STEP 6: Plot transect 2

ftsz = [20 16];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
% iylim = [-13,13];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 0.5;
xlab = '$x$ (m)';
ylab = '$z$ (m)';
cmap = cmocean('thermal',numframes);

sname =[trialname,'_timeseries_2'];
v = VideoWriter([figfolder,'\',sname], 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)

figure('units','inches','position',[1 1 10 3],'color','w');

for i = 1:numframes+1
    clf
    plot(x(1,ix2:ix),zbeach(ix2:ix),'LineWidth',3,'Color',[0.4 0.4 0.4])
hold on
plot(x(1,1:ix2),zbeachmean(1:ix2),'LineWidth',3,'Color',[0.4 0.4 0.4])
    plot(x(1,ixmean:end),ztran2(ixmean:end,i),'LineWidth',3,'Color','b')

box on
ylim([0 4.5]);
xlim([0 100]);
clim([0 (numframes-1)/2])
colormap(cmap)
h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
title(['Transect 2: $t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
pause(0.1)
writeVideo(v,getframe(gcf))
end
close(v)

%% 
t = [0:0.5:(numframes-1)/2];
figure; 
[tgrid,xgrid] = meshgrid(t,x(1,:));
pcolor(xgrid,tgrid,ztran)
shading flat


%% Functions

function [RX, RY] = rotateCoordinates(X, Y, LocalX, LocalY, theta)
    % ROTATECOORDINATES Rotates the coordinates (X, Y) around the origin
    % defined by (LocalX, LocalY) by an angle theta.
    %
    % Inputs:
    %   - X      : nx1 or 1xn vector of X coordinates
    %   - Y      : nx1 or 1xn vector of Y coordinates
    %   - LocalX : scalar value for the X coordinate of the local origin
    %   - LocalY : scalar value for the Y coordinate of the local origin
    %   - theta  : scalar value of the rotation angle in degrees
    %
    % Outputs:
    %   - Xout   : nx1 or 1xn vector of rotated X coordinates
    %   - Yout   : nx1 or 1xn vector of rotated Y coordinates
    
    % Translate coordinates to the origin
    Xo = X - LocalX;
    Yo = Y - LocalY;
    
    % Perform the rotation
    RY = Yo .* cosd(theta) + Xo .* sind(theta);
    RX = Xo .* cosd(theta) - Yo .* sind(theta);

    % Translate to local coordinates
    RX = RX-21;
end

function dataOut = replaceShortNonNanSequences(data, numstring)
    % Initialize output as a copy of the input data
    dataOut = data;
    
    % Find indices where data is not NaN
    isNotNan = ~isnan(data);
    
    % Find start and end indices of consecutive non-NaN sequences
    d = diff([0; isNotNan; 0]); % Add 0 at start and end to capture edge sequences
    starts = find(d == 1);    % Indices where sequences start
    ends = find(d == -1) - 1; % Indices where sequences end
    
    % Loop through each sequence
    for i = 1:length(starts)
        seqLength = ends(i) - starts(i) + 1;
        if seqLength < numstring
            % Replace short sequences with NaN
            dataOut(starts(i):ends(i)) = NaN;
        end
    end
end


function h1 = plotstyleCMB(gca,xlab,ylab,ftsize,ticklen,lnwidth,tickminor,tickdir)
    % function to style plot in CMB (and Nirni) fashion
    % INPUT:
    % gca:      figure
    % xlab,ylab: labels for the x- and y-axis (strings) (default: x/y)
    % ftsize:   fontsize for labels and numbers [lab, num] (default: [20,15])
    % ticklen:  tick length (default: 1)
    % lnwid:    linewidth for figure (default: 1.2)
    %
    % Created by Baker, C.M. on 2022-Nov-30
    
    if ~exist('ftsize','var'); ftsize=[20 15]; end
    if ~exist('ticklen','var'); ticklen = 1; end
    if ~exist('lnwidth','var'); lnwidth = 1.2; end
    if ~exist('tickminor','var'); tickit = 0; else; tickit = 1; end
    if ~exist('tickdir','var'); tickdir = 'in'; end
    
    h1=gca;
    set(h1,'ticklength',ticklen*get(h1,'ticklength'),'tickdir',tickdir);
    set(h1,'fontsize',ftsize(2),'LineWidth',lnwidth);
    ylabel(ylab,'interpreter','latex','fontsize',ftsize(1))
    xlabel(xlab,'interpreter','latex','fontsize',ftsize(1))
    if tickit == 1%isstring(tickminor)
        set(h1,'xminortick',tickminor,'yminortick',tickminor);
    end
    set(h1,'layer','top')
end
