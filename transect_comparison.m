%% Compare GEM Transect with Hand Surveys
% Compare using Roundgridfun with 1, 0.5, and 0.2 m steps (bins)
clc
clear all
close all
format long
% general path and names
datapath    = '/Users/bagaenzl/Desktop/MasonBEAST Data/';
trialname   = '1721934001780';%1722524460162';%'1702827001820';%'1699459201959';%'1702827001820';%'1698951602001';
tempname   = '1721934001780';
%list of epoch nums using
% % 1717167601959 = May 31 11 AM 2024 x
% 1719428401397 = June 26 3PM 2024 x
% 1721934001780 = July 25 3PM 2024 x looked super weird
% '1723489201189' = Aug 12 2024
% 1726056001556 = Sept 11 8 AM 2024
% '1726412401052' = Sept 15 11AM 2024
% '1726772401770' = Sept 19 3PM 2024
% 1726858801716 = Sept 20 3PM 2024
% '1728327601158' = Oct 7 3PM 2024
% 1728388801457 = Oct 8 8AM 2024
% '1728561601419' = Oct 10 8AM 2024
% '1730736001873' = Nov 4 11PM 2024
% 1732464001262 = Nov 24 11 AM 2024
%gcpname = [datapath,'GCP Locations txts/gcps/GCPs_09_18_2024.txt'];

%orthopath = [datapath,'stereo/'];
pcpath = [datapath,'PointClouds/',trialname];
%numframes = length(dir(fullfile(pcpath, '*.txt')));
numframes=2;
format long
surveyname = [datapath,'surveys/2024_07_25_Transects_UTM.xlsx'];

% create grid
grid_step=0.5; % bin size change to 1 m, 0.5 m, and 0.2 m
Xgrid=-7:grid_step:60; % m from GEM and transect plot bounds
Ygrid=-30:grid_step:27; % m from GEM and transect bounds

% local x and y coords and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35; % rotate to cross-shore

% read in hand surveys
handtran=readmatrix(surveyname);
Xtran=handtran(:,3); % Easting
Ytran=handtran(:,2); % Northing
Ztran=handtran(:,1); % meters

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
        % (with point cloud folders)ptcl = readmatrix([pcpath,'/',tempname,'_ptcld',num2str(i),'.txt']); % columns x,y,z = 1,2,3
        ptcl = readmatrix([pcpath,'_ptcld',num2str(i),'.txt']); % columns x,y,z = 1,2,3
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
fig=figure; ax1=subplot(1,4,1);pcolor(Xgrid,Ygrid,MeanGEMz(:,:,1)); grid off; shading flat;
hold on; title("Averaged GEM Elevation Values"); colorbar;c1=clim;
hold on; ax2=subplot(1,4,2);pcolor(Xgrid,Ygrid,ZtranMean); grid off; shading flat; title("Averaged Hand Transect Elevation Values");
colorbar;c2 = clim;
%c2.Label.String = 'Elevation (m NAD83 (2011))'
c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
subplot(1,4,1); colorbar off; subplot(1,4,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
hold on; ax3=subplot(1,4,3); pcolor(Xgrid,Ygrid,ZtranMean-MeanGEMz(:,:,1)); grid off; shading flat; title("Difference")
c3 = colorbar; colormap(ax3,"hot");
c3.Label.String = 'Difference in Elevation (m NAD83 (2011))'; hold on;
ax4=subplot(1,4,4); pcolor(Xgrid,Ygrid,rmse(MeanGEMz(:,:,1),ZtranMean,[58 68],"omitnan")); grid off; shading flat; title("RMSE per Bin");
c4=colorbar; colormap(ax4,"summer");c4.Label.String = 'RMSE (m NAD83 (2011))';
linkaxes([ax1 ax2 ax3 ax4]);
sgtitle("07/25/24 with 0.5 m Bins");set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen

% Roundgridfun example
% x = rand(1000,1)*10; y = rand(1000,1)*10; z = x.*y;
%     xgi = 0:1:10; ygi = 0:1:10;
%     [val,numpts]=roundgridfun(x,y,z,xgi,ygi,@min); % Minimum Z Surface
%     figure; pcolor(xgi,ygi,val);

A=rmse(MeanGEMz(:,:,1),ZtranMean,[58 68],"omitnan");
save('1721934001780_RMSE.mat','A');




%% RMSE comparison over time
% load files
rm_1717167601959=load('1717167601959_RMSE.mat');
rm_1719428401397=load('1719428401397_RMSE.mat');
%rm_1723489201189=load('1723489201189_RMSE.mat');

% plot
fig=figure; ax1=subplot(1,2,1);pcolor(Xgrid,Ygrid,rm_1717167601959.A); grid off; shading flat;
hold on; title("May 31"); colorbar;c1=clim;
hold on; ax2=subplot(1,2,2);pcolor(Xgrid,Ygrid,rm_1719428401397.A); grid off; shading flat; title("June 26");
colorbar;c2 = clim; hold on;
%ax3=subplot(1,3,3); pcolor(Xgrid,Ygrid,rm_1730736001873.A); grid off, shading flat; title("November 4th"); colorbar;c3=clim;
%hold on;
%c2.Label.String = 'Elevation (m NAD83 (2011))'
c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
subplot(1,2,1); colorbar off; subplot(1,2,2); colorbar off;% subplot(1,3,3);colorbar off;
cbar=colorbar('EastOutside'); cbar.Label.String = 'RMSE (m)'; colormap("spring"); caxis([0 , 0.5])


%% Comparing GEMs with and without manual camera locations
% Paths
folderpath='/Users/bagaenzl/Desktop/MasonBEAST Data/PointClouds/';
filePattern=fullfile(folderpath,'*.txt');
listofFiles=dir(filePattern);
tranfolderpath='/Users/bagaenzl/Desktop/MasonBEAST Data/surveys/';
% What GEMs you are looking at in epoch time
dates_wanted=["1723489201189" "1726772401770" "1730736001873"];
% Transects correlated to each GEM epoch time
transects_wanted=["2024_08_12" "2024_09_18" "2024_11_04" ];

% Define variables
Xloc = 239737; % define local coordinate system origin and rotation angle
Yloc = 3784751;
rotang = 35; % rotation angle to alongshore

% create grid
grid_step=0.5; % bin size change to 1 m, 0.5 m, and 0.2 m
Xgrid=-7:grid_step:60; % m from GEM and transect plot bounds
Ygrid=-30:grid_step:27; % m from GEM and transect bounds


ptcld_struct=struct; % structure to store wanted ptclds in
% Loop through given epoch numbers
for k=1:length(listofFiles)
    baseFilename=listofFiles(k).name;
    fullFilename=fullfile(listofFiles(k).folder,baseFilename);
    for i=1:length(dates_wanted)
        transectFilename=[tranfolderpath,transects_wanted(i)];
        TF=contains(baseFilename,dates_wanted(i));
        if TF==1
            ptcld_struct(k).ptclds=fullFilename; % This adds ptclouds we want into a structure
            ptcld_struct(k).transects=transectFilename; % Adds the path to the hand survey transects (transect
        end
    end
    fun = @(s) all(structfun(@isempty,s)); % found off of matlab forum (https://www.mathworks.com/matlabcentral/answers/613956-how-to-check-if-a-element-of-a-struct-is-empty)
    idx = arrayfun(fun,ptcld_struct);
    ptcld_struct(idx)=[]; % remove the empty elements
end
% (works fine up until here)

% % create matrix for storing z values
numframes=2;
% z = NaN(size(x,1),size(x,2),numframes);
% numpts = z;
format long g

% loop through ptcld structure and create GEMs
for j=1:length(ptcld_struct)
    for i = 1:numframes
    % read point cloud
        if i == 25
            ptcl = [0 0 0];
        else
            %ptcld_toread=load(ptcld_struct(j).ptclds);
            ptcl = readmatrix([ptcld_struct(j).ptclds]); % columns x,y,z
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

    GEMname=split(ptcld_struct(j).ptclds,'/');
    GEMname=split(GEMname(7,1),'.');
    GEMname=GEMname(1,1);
    GEMname=string(GEMname);
    GEMfilepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/',GEMname);
    save(GEMfilepath,"MeanGEMz",'-mat') % save z values from GEM as mat file

    % read in hand surveys
    transectfilepath=append(ptcld_struct(j).transects(1),ptcld_struct(j).transects(2),'_Transects_UTM.xlsx');
    handtran=readmatrix(transectfilepath);
    Xtran=handtran(:,3); % Easting
    Ytran=handtran(:,2); % Northing
    Ztran=handtran(:,1); % meters

    % Rotate Hand Survey
    % coordinate rotation and transformation x,y to cross- and alongshore
    [Xrottran, Yrottran] = rotateCoordinates(Xtran,Ytran, Xloc, Yloc, rotang);
    % Calculate average Z value for hand surveys
    [ZtranMean,TranNumpts]=roundgridfun(Xrottran,Yrottran,Ztran,Xgrid,Ygrid,@mean);
    ZtranMean(ZtranMean==0)=NaN;
    ZtranMean_filepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files',ptcld_struct(j).transects(2));
    save(ZtranMean_filepath,"ZtranMean",'-mat'); % save ztran values as mat file

    % Plot GEM with Hand Survey, Difference, and RMSE
    comp_figfilepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/',GEMname,'_comparison');
    % Need to specify size, and place to store it
    fig=figure; ax1=subplot(1,4,1);pcolor(Xgrid,Ygrid,MeanGEMz(:,:,1)); grid off; shading flat;
    hold on; title("Averaged GEM Elevation Values"); colorbar;c1=clim;
    hold on; ax2=subplot(1,4,2);pcolor(Xgrid,Ygrid,ZtranMean); grid off; shading flat; title("Averaged Hand Transect Elevation Values");
    colorbar;c2 = clim;
    %c2.Label.String = 'Elevation (m NAD83 (2011))'
    c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
    subplot(1,4,1); colorbar off; subplot(1,4,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
    hold on; ax3=subplot(1,4,3); pcolor(Xgrid,Ygrid,ZtranMean-MeanGEMz(:,:,1)); grid off; shading flat; title("Difference")
    c3 = colorbar; colormap(ax3,"hot");
    c3.Label.String = 'Difference in Elevation (m NAD83 (2011))'; hold on;
    ax4=subplot(1,4,4); pcolor(Xgrid,Ygrid,rmse(MeanGEMz(:,:,1),ZtranMean,[58 68],"omitnan")); grid off; shading flat; title("RMSE per Bin");
    c4=colorbar; colormap(ax4,"summer");c4.Label.String = 'RMSE (m NAD83 (2011))';
    linkaxes([ax1 ax2 ax3 ax4]);
    sgtitle(GEMname(1,1));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
    saveas(fig,comp_figfilepath,'png');

    % Compare cam man and norm
    TF1=contains(ptcld_struct(j).ptclds,dates_wanted(1));
    TF2=contains(ptcld_struct(j).ptclds,dates_wanted(2));
    TF3=contains(ptcld_struct(j).ptclds,dates_wanted(3));
    % loop through and create plots to compare 1st frame of man cam and
    % norm GEM for each date
        if TF1==1
        cam_man=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1723489201189_cam_man_ptcld1']); % load manually camer loc GEM
        norm=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1723489201189_ptcld1']); % load GEM
        figfilepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/',dates_wanted(1),'_camera_loc_comparison');
        %plot
        fig=figure; ax1=subplot(1,3,1);pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,2)); grid off; shading flat;
        hold on; title("Manual Camera Location"); colorbar;c1=clim;
        hold on; ax2=subplot(1,3,2);pcolor(Xgrid,Ygrid,norm.MeanGEMz(:,:,1)); grid off; shading flat; title("No Camera Location Input");
        colorbar;c2 = clim;
        c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
        subplot(1,3,1); colorbar off; subplot(1,3,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
        hold on; ax3=subplot(1,3,3); pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,1)-norm.MeanGEMz(:,:,1)); grid off; shading flat; title("Difference")
        c3 = colorbar; colormap(ax3,"hot");
        c3.Label.String = 'Difference in Elevation (m NAD83 (2011))'; hold on;
        linkaxes([ax1 ax2 ax3]);
        sgtitle(dates_wanted(1));
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
        saveas(fig,figfilepath,'png');

        elseif TF2==1
            cam_man=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1726772401770_cam_man_ptcld1']); % load manually camer loc GEM
            norm=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1726772401770_redo_ptcld1']); % load GEM
            figfilepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/',dates_wanted(2),'_camera_loc_comparison');
            %plot
            fig=figure; ax1=subplot(1,3,1);pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,1)); grid off; shading flat;
            hold on; title("Manual Camera Location"); colorbar;c1=clim;
            hold on; ax2=subplot(1,3,2);pcolor(Xgrid,Ygrid,norm.MeanGEMz(:,:,1)); grid off; shading flat; title("No Camera Location Input");
            colorbar;c2 = clim;
            c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
            subplot(1,3,1); colorbar off; subplot(1,3,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
            hold on; ax3=subplot(1,3,3); pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,1)-norm.MeanGEMz(:,:,1)); grid off; shading flat; title("Difference")
            c3 = colorbar; colormap(ax3,"hot");
            c3.Label.String = 'Difference in Elevation (m NAD83 (2011))'; hold on;
            linkaxes([ax1 ax2 ax3]);
            sgtitle(dates_wanted(2));
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
            saveas(fig,figfilepath,'png');

        elseif TF3==1
            cam_man=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1730736001873_cam_man_ptcld1']); % load manually camer loc GEM
            norm=load(['/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/','1730736001873_ptcld1']); % load GEM
            figfilepath=append('/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_comp_files/',dates_wanted(3),'_camera_loc_comparison');
            %plot
            fig=figure; ax1=subplot(1,3,1);pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,1)); grid off; shading flat;
            hold on; title("Manual Camera Location"); colorbar;c1=clim;
            hold on; ax2=subplot(1,3,2);pcolor(Xgrid,Ygrid,norm.MeanGEMz(:,:,1)); grid off; shading flat; title("No Camera Location Input");
            colorbar;c2 = clim;
            c=[min([c1(1) c2(1)]),max([c1(2) c2(2)])];
            subplot(1,3,1); colorbar off; subplot(1,3,2); colorbar off; cbar=colorbar('EastOutside'); cbar.Label.String = 'Elevation (m NAD83 (2011))';
            hold on; ax3=subplot(1,3,3); pcolor(Xgrid,Ygrid,cam_man.MeanGEMz(:,:,1)-norm.MeanGEMz(:,:,1)); grid off; shading flat; title("Difference")
            c3 = colorbar; colormap(ax3,"hot");
            c3.Label.String = 'Difference in Elevation (m NAD83 (2011))'; hold on;
            linkaxes([ax1 ax2 ax3]);
            sgtitle(dates_wanted(3));
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
            saveas(fig,figfilepath,'png');

        end
end









