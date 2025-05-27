function[meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,dxy,numframes,pcpath,GEMsavepath,figpath,spec)
% function[meanGEMz,medGEMz]=ptcld_to_GEM(datapath,savepath)
% --------------------------------------------------------------------------
% This function takes in a pointcloud txt file from Metashape and creates a
% Gridded Elevation Model and saves it as a matrix in a mat file to a 
% desired location. It also creates a figure and saves it to desired
% location. 
% --------------------------------------------------------------------------
% INPUTS:
% -------
% 
% camloc = camera A location coordinates [Ax,Ay]
% camloc = camera B location coordinates [Bx,By]
% dxy = grid bin size (m)
% numframes = number of frames to process per poinmtcloud
%
% pcpath = path to Metashape pointcloud (.txt file) 
% GCPpath = path to GCP location (.txt file) THIS IS NOT ADDED IN YET BUT
% EVENTUALLY
% -GCP file must have elevation values in 4th column
% GEMsavepath = path to save GEM mat file
% figpath = path to save figure
%
% spec = string of specified pointcloud to look at (i.e, '*ptcld1',
% '*meta_ptcld',or a specific epoch number)
% 
% OUTPUTS:
% --------
% meanGEMz = matrix of mean elevation values from given
% pointcloud as a mat file
% medGEMz = matrix of median elevation values from given
% pointcloud as a mat file
%
% --------------------------------------------------------------------------
format long g

% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;

% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);

% % create matrix for storing z values
% z = NaN(size(x,1),size(x,2),numframes);
% numpts = z;

% GCP Input
%GCPs = readmatrix(GCPpath);
%GCPz = GCPs(:,4); % I seperate the x, y and z but is not required
% gcp rotation and transformation
%[GCPx, GCPy] = rotateCoordinates(GCPs(:,2), GCPs(:,3), Xloc, Yloc, rotang); 

% camera rotation and transformation
[CamAx, CamAy] = rotateCoordinates(camlocA(1), camlocA(2), Xloc, Yloc, rotang);
[CamBx, CamBy] = rotateCoordinates(camlocB(1), camlocB(2), Xloc, Yloc, rotang);
% 
% Loop through given pointclouds
filePattern=fullfile(pcpath,spec);
listofFiles=dir(filePattern);
ptcld_struct=struct(); % structure to store wanted ptclds in

for k=1:length(listofFiles)
    baseFilename=listofFiles(k).name;
    fullFilename=fullfile(listofFiles(k).folder,baseFilename);
    fprintf(1, 'Now reading %s\n', fullFilename);
    % for i=1:length(dates_wanted)
    %     transectFilename=[tranfolderpath,transects_wanted(i)];
    %     TF=contains(baseFilename,dates_wanted(i));
    %     if TF==1
    %         ptcld_struct(k).ptclds=fullFilename; % This adds ptclouds we want into a structure
    %         listofFiles(k).transects=transectFilename; % Adds the path to the hand survey transects (transect
    %     end
    % end
end

%
% Create median and mean GEM using the Pointcloud (Modified from CM Baker)
% point cloud is in NAVD83 (2011) UTM Zone 18 N EPSG 6347

for j=1:length(listofFiles)
    for i = 1:numframes
    % read point cloud
        if i == 25
            ptcl = [0 0 0];
        else
            %ptcld_toread=load(ptcld_struct(j).ptclds);
            ptcl = readmatrix(append(listofFiles(j).folder,'/',listofFiles(j).name)); % columns x,y,z
        end
        % pt cloud coordinate rotation and transformation x,y to cross- and alongshore
        [Xrot, Yrot] = rotateCoordinates(ptcl(:,1), ptcl(:,2), Xloc, Yloc, rotang);

    % grid point cloud
    % Mean 
        [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),Xgrid,Ygrid,@mean); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
        ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
        ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan
        % store output into a matrix
        meanGEMz(:,:,i) = ztemp; % median (or mean) values of all of the point cloud frames
        numpts(:,:,i) = ntemp;
        clear ztemp ntemp
    % Median
        [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),Xgrid,Ygrid,@median); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
        ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
        ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan
        % store output into a matrix
        medGEMz(:,:,i) = ztemp; % median (or mean) values of all of the point cloud frames
        numpts(:,:,i) = ntemp;
        clear ztemp ntemp

    end

    % save mean and median elevation values
    GEMname=split(listofFiles(j).name,'_');
    GEMname=GEMname(1,1);
    GEMname=string(GEMname);
    GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
    GEMdate=string(GEMdate);
    GEMtitle=append(GEMname,',',GEMdate);
    GEMfilepath=append(GEMsavepath,GEMname,'/');

    % check if folder exists and save
    if isfolder(GEMfilepath) == true
        save GEMfilepath meanGEMz medGEMz;
    else
        mkdir(GEMfilepath);
    end 

    matname=fullfile(GEMfilepath,'meanGEMz.mat');
    save(matname,'meanGEMz')
    matname=fullfile(GEMfilepath,'medGEMz.mat');
    save(matname,'medGEMz')

    % Plotting GEM with mean elevation
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanGEMz(:,:,1)); grid off;box on;hold on
    %scatter(GCPx,GCPy,60,'fill','r','MarkerEdgeColor','k') (need to figure
    %out how to specifiy GCPs/do we need?)
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    title(GEMtitle);
    filename=append('mean',GEMname);
    meanfigpath=fullfile(figpath, filename);
    %meanfigpath=append(figpath,"/mean",GEMname);
    saveas(fig,meanfigpath,'png');

    % Plotting GEM with median elevation
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,medGEMz(:,:,1)); grid off;box on;hold on
    %scatter(GCPx,GCPy,60,'fill','r','MarkerEdgeColor','k') (need to figure
    %out how to specifiy GCPs/do we need?)
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm'); shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    title(GEMtitle);
    filename=append('med',GEMname);
    medfigpath=fullfile(figpath, filename);
    %meanfigpath=append(figpath,"/mean",GEMname);
    saveas(fig,medfigpath,'png');

end


%% TEST
% camlocA = [239766.1, 3784761.9];%, 10.37
% camlocB = [239759.4, 3784755.0];%, 10.26];
% dxy = 0.5; % meter
% numframes=2;
% 
% pcpath='/Users/bagaenzl/Desktop/MasonBEAST Data/PointClouds/';
% /Volumes/kanarde-1/MasonBEAST/data/PointClouds
% GEMsavepath='/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_mat_files/';
%/Volumes/kanarde-1/MasonBEAST/data/GEMs/Camera_Location_Analysis
% figpath='/Users/bagaenzl/Desktop/MasonBEAST Data/Figures';
% 
% spec='1702827001820_ptcld1.txt';
% 
% [meanGEMz,medGEMz]=ptcld_to_GEM(camlocA,camlocB,dxy,numframes,pcpath,GEMsavepath,figpath,spec);


%% Edits NEEDED
% increase font size 
% sve a figure per frame