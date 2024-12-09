%% Compare GEM Transect with Hand Surveys
% Compare using Roundgridfun with 1, 0.5, and 0.2 m steps (bins)
% general path and names
datapath    = '/Users/bagaenzl/Desktop/MasonBEAST Data/';
trialname   = '1723489201189';%1722524460162';%'1702827001820';%'1699459201959';%'1702827001820';%'1698951602001';
tempname   = '1723489201189';
gcpname = [datapath,'GCP Locations txts/gcps/GCPs_08_01_2024.txt'];

%orthopath = [datapath,'stereo/'];
pcpath = [datapath,'PointClouds/',trialname]; % need to make this a bit more consise and intuitive
%numframes = length(dir(fullfile(pcpath, '*.txt')));
numframes=2;
format long
surveyname = [datapath,'surveys/2024_08_12_TransectsUTM.xls'];

% create grid 
grid_step=1; % bin size change to 1 m, 0.5 m, and 0.2 m
Xgrid=-7:grid_step:60; % m from GEM and transect plot bounds
Ygrid=-30:grid_step:27; % m from GEM and transect bounds

% local x and y coords and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35; % rotate to cross-shore

% read in hand surveys
handtran=readmatrix(surveyname);
Xtran=handtran(:,5); % Northing
Ytran=handtran(:,4); % Easting
Ztran=handtran(:,2); % meters

% Rotate Hand Survey 
% coordinate rotation and transformation x,y to cross- and alongshore
[Xrottran, Yrottran] = rotateCoordinates(Xtran,Ytran, Xloc, Yloc, rotang);

% GEM Data using the Pointcloud ( By CM Baker)
% point cloud is in NAVD83 (2011) UTM Zone 18 N EPSG 6347
for i = 1:numframes
    % read point cloud
    if i == 25
        ptcl = [0 0 0];
    else
        ptcl = readmatrix([pcpath,'/',tempname,'_ptcld',num2str(i),'.txt']); % columns x,y,z = 1,2,3
    end

    % pt cloud coordinate rotation and transformation x,y to cross- and alongshore
    [Xrot, Yrot] = rotateCoordinates(ptcl(:,1), ptcl(:,2), Xloc, Yloc, rotang);

    % grid point cloud
    [ztemp,ntemp]  = roundgridfun(Xrot,Yrot,ptcl(:,3),Xgrid,Ygrid,@mean); % computes median or mean of binned point cloud with xpt, ypt, zpt values at resolution of xgrid, ygrid
    ztemp(ztemp == 0) = NaN; % z is the gridded elevations, rounding grid function sets locations without points equal to zero, switching to nan
    ntemp(ntemp == 0) = NaN; % n is the number of points per bin, rounding grid function sets locations without points equal to zero, switching to nan

    % store output into a matrix
    MeanGEMz(:,:,i) = ztemp; % median (or mean) values of all of the point cloud frames
    numpts(:,:,i) = ntemp;
    clear ztemp ntemp

end

[ZtranMean,TranNumpts]=roundgridfun(Xrottran,Yrottran,Ztran,Xgrid,Ygrid,@mean);
ZtranMean(ZtranMean==0)=NaN;

% Plot
fig=figure; ax1=subplot(1,3,1);pcolor(Xgrid,Ygrid,MeanGEMz(:,:,1));
hold on; title("Averaged 08/12/24 GEM Elevation Values with 1 m Bins"); colorbar;c1=caxis; 
hold on; ax2=subplot(1,3,2);pcolor(Xgrid,Ygrid,ZtranMean); title("Averaged 08/12/24 Hand Transect Elevation Values with 1 m Bins");
colorbar;c2 = caxis; 
%c2.Label.String = 'Elevation (m NAD83 (2011))'
c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
subplot(1,3,1); colorbar off; subplot(1,3,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
hold on; ax3=subplot(1,3,3); pcolor(Xgrid,Ygrid,ZtranMean-MeanGEMz(:,:,1)); title("Difference with 1 m Bins")
c3 = colorbar; colormap(ax3,"hot");
c3.Label.String = 'Difference in Elevation (m NAD83 (2011))';

% Roundgridfun example
% x = rand(1000,1)*10; y = rand(1000,1)*10; z = x.*y;
%     xgi = 0:1:10; ygi = 0:1:10;
%     [val,numpts]=roundgridfun(x,y,z,xgi,ygi,@min); % Minimum Z Surface
%     figure; pcolor(xgi,ygi,val);









