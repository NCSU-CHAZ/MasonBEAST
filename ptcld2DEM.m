function[meanDEMz,medDEMz]=ptcld2DEM(camlocA,camlocB,rbrloc,dxy,pcpath,savepath,spec)
% function[meanDEMz_matrix,medDEMz_matrix]=ptcld2DEM(camlocA,camlocB,rbrloc,dxy,pcpath,savepath,spec)
% % --------------------------------------------------------------------------
% This function takes in a pointcloud txt file from Metashape and creates a
% Digital Elevation Model (DEM). The DEM is saved as a matrix .mat file
% within the savepath. Static DEM plots of each point in time are saved
% within a figure folder in the savepath, along with a video of the DEM
% over time.
% --------------------------------------------------------------------------
% INPUTS:
% -------
% camloc = camera A location coordinates [Ax,Ay]
% camloc = camera B location coordinates [Bx,By]
% rbrloc = rbr location coordinates [x,y]
% dxy = grid bin size (m)
% pcpath = path to Metashape pointclouds (.txt file) 
% savepath = path to save DEM matrix mat file, figures are stored under
% here as well
% spec = string of specified pointcloud to look at (i.e, '*ptcld1',
% '*meta_ptcld',or a specific epoch number)
%
% OUTPUTS:
% --------
% meanDEMz_matrix = matrix of mean elevation values from given
% pointcloud as a mat file
% medDEMz_matrix = matrix of median elevation values from given
% pointcloud as a mat file
%
% Last Edits: BG 08/17/2026
% 
% EDITS NEEDED: 
% ---------------------------
% leave space in matrix for null pointcloud
% add video of DEM timesteps
% 
% POTENTIALLY ADD IF NEEDED
% % GCPpath = path to GCP location (.txt file) THIS IS NOT ADDED IN YET BUT
% EVENTUALLY
% -GCP file must have elevation values in 4th column
% ----------------------------------------------------------------------------
format long g

% CONSTANTS
% -------------------------
% define figure paths
DEMsavepath=savepath; % for matrix of DEMs
DEMfigpath=append(savepath,'Figures'); % for static plots and video
% create save folders
if ~isfolder(DEMfigpath)
    mkdir(DEMfigpath);
end
if ~isfolder(DEMsavepath)
    mkdir(DEMsavepath);
end
% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;
% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
% camera rotation and transformation
[CamAx, CamAy] = rotateCoordinates(camlocA(1), camlocA(2), Xloc, Yloc, rotang);
[CamBx, CamBy] = rotateCoordinates(camlocB(1), camlocB(2), Xloc, Yloc, rotang);
% rbr rotation and transformation
[rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);

% LOAD POINTCLOUDS
% ------------------------------
% Loop through given pointclouds
filePattern=fullfile(pcpath,spec);
listofFiles=dir(filePattern);
listofFiles=listofFiles(~startsWith({listofFiles.name},'.')); % exclude gen files
sortfilenames=natsortfiles({listofFiles.name});
numframes=length(sortfilenames);

% CREATE DEMs (STATIC PLOTS AND .MAT FILE)
%----------------------------------------------
% create matrix for storing z values
meanDEMz = NaN(size(Xgrid,1),size(Xgrid,2),numframes);
numpts = meanDEMz;
medDEMz = NaN(size(Xgrid,1),size(Xgrid,2),numframes);

% Create median and mean GEM using the Pointcloud (Modified from CM Baker)
% point cloud is in NAVD83 (2011) UTM Zone 18 N EPSG 6347
% save z values in a matrix for each pointcloud
for i = 1:length(listofFiles) 
    % read point cloud
    baseFilename=string(sortfilenames(i));
    fullFilename=fullfile(listofFiles(i).folder,baseFilename);
    fprintf(1, 'Now reading %s\n', fullFilename);
    ptcl = readmatrix(append(listofFiles(i).folder,'/',baseFilename)); % columns x,y,z
    % pt cloud coordinate rotation and transformation x,y to cross- and alongshore
    [Xrot, Yrot] = rotateCoordinates(ptcl(:,1), ptcl(:,2), Xloc, Yloc, rotang);
    % grid point cloud
    % Mean 
    [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),Xgrid,Ygrid,@mean); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
     ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
     ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan
     % store output into a matrix
     meanDEMz(:,:,i) = ztemp; % median (or mean) values of all of the point cloud frames
     numpts(:,:,i) = ntemp;
     clear ztemp ntemp
    % Median
    [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),Xgrid,Ygrid,@median); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
    ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
    ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan
     % store output into a matrix
     medDEMz(:,:,i) = ztemp; % median (or mean) values of all of the point cloud frames
     numpts(:,:,i) = ntemp;
     clear ztemp ntemp 

    DEMname=split(fullFilename,'/');
    DEMname=DEMname(6,1);
    DEMname=split(DEMname,'_');
    DEMname=DEMname(1,1);
    DEMname=string(DEMname);
    DEMdate=datetime(str2num(DEMname),'ConvertFrom','epochtime','TicksPerSecond',1000); % epochs in milliseconds
    DEMdate=string(DEMdate);
    DEMtitle=append(DEMname,',',DEMdate);
    DEMfilepath=append(DEMsavepath,'DEMz_matrix'); 
    
% Plotting DEM with mean elevation
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanDEMz(:,:,i)); grid off;box on;hold on;
    %yvals=-60:2:6; %this is for showing transect lines for averaging
    %yline(yvals,'--','Color',[0.7 0.7 0.7],'Linewidth',1.5); hold on;
    %yvals=8.3:2:22.3;
    %yline(yvals,'--','Color',[0.7 0.7 0.7],'Linewidth',1.5); hold on;
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    scatter(rbrx,rbry,60,'fill','sq','k','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'k');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw); hc.Label.String = 'Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);
    title(DEMtitle);
    meanfigpath=append(DEMfigpath,"/meanDEM",DEMname,'_',num2str(i));
    saveas(fig,meanfigpath,'png');
    fprintf('saved meanDEM fig: %s to %s\n', DEMtitle, meanfigpath);
    close(fig);

    % Plotting GEM with median elevation
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,medDEMz(:,:,i)); grid off;box on;hold on
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm'); 
    scatter(rbrx,rbry,60,'fill','sq','k','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'k');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw); hc.Label.String = 'Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);
    title(DEMtitle);
    medfigpath=append(DEMfigpath,"/medDEM",DEMname,'_',num2str(i));
    saveas(fig,medfigpath,'png');
    fprintf('saved medDEM fig: %s to %s\n', DEMtitle, medfigpath);
    close(fig);
end

 % save mean and median elevation values
matname=fullfile(savepath,append('/meanDEMz_matrix','.mat'));
save(matname,'meanDEMz')
fprintf('saved meanDEMz matrix to: %s\n', matname);
matname=fullfile(savepath,append('/medDEMz_matrix','.mat'));
save(matname,'medDEMz')
fprintf('saved medDEMz matrix to: %s\n', matname);


