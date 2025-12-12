%% Script to process monthly GEMs (camera location comparison)

% paths 
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
pcpath=append(genpath,'/PointClouds/'); % path to pointclouds
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % save measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % save metashape GEMs

% GEM specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.5; % meter change to 0.2v normally
numframes=2;

% Loop through Metashape pointclouds
spec='*meta_ptcld*';
figpath=append(metashape_path,'Figures');
[meanMAPzmatrix,medMAPzmatrix,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,metashape_path,figpath,spec);


% compare to hand surveys (need to figure out how to grab surveys through
% loop) 
%HS_savepath=append(genpath,'/Surveys/rotated_surveys');
%surveypath=append(CAM_analysispath,'/Surveys');
%listofsurveys=dir(surveypath);
%for i=2:length(listofsurveys)
    %listofsurvs(i)=append(listofsurveys(i).folder,listofsurveys(i).name);
%end

%listofFiles=dir(metashape_path);
%for i=1:length(listofFiles)
    %GEMmatpath=append(listofFiles(i).folder,'/',listofFiles(i).name);
    %[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
%end

% Loop through Measured pointclouds
spec='*meas_ptcld*';
figpath=append(measured_path,'Figures');
[meanMAPzmatrix,medMAPzmatrix]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,measured_path,figpath,spec);


% GEM comparison manually
HS_savepath=append(genpath,'/Surveys/rotated_surveys');
surveypath=append(CAM_analysispath,'/Surveys');
% epoch nums: 1708030441768, 1711479601066,1713452401874, 1717156801902,
% 1719428401397,1721934001780, 1723489201189, 1726772401770, 1728561601419,
% 1730736001873, 1734181201371, 1735848001721, 1738526401587, 1741978801668
% survey names: '2024_03_26_Transects_UTM.xlsx'
% '2024_05_31_Transects_UTM.xlsx' '2024_06_26_Transects_UTM.xlsx'
% '2024_07_25_transects_UTM.xlsx' '2024_08_12_Transects_UTM.xlsx'
% '2024_09_18_Transects_UTM.xlsx' '2024_10_01_Transects_UTM.xlsx'
% '2024_11_04_Transects_UTM.xlsx' '2025_01_29_Transects_UTM.xlsx'
% '2025_03_14_Transects_UTM.xlsx'

% epochs
epochs={'1711479601066','1713452401874', '1717156801902','1719428401397','1721934001780', '1723489201189', '1726772401770', '1728561601419','1730736001873','1734181201371', '1735848001721', '1738526401587', '1741978801668'};
surveynames={'2024_03_26_Transects_UTM.xlsx', '2024_03_26_Transects_UTM.xlsx','2024_05_31_Transects_UTM.xlsx' ,'2024_06_26_Transects_UTM.xlsx','2024_07_25_transects_UTM.xlsx', '2024_08_12_Transects_UTM.xlsx','2024_09_18_Transects_UTM.xlsx', '2024_10_01_Transects_UTM.xlsx','2024_11_04_Transects_UTM.xlsx', '2025_01_29_Transects_UTM.xlsx','2025_01_29_Transects_UTM.xlsx','2025_01_29_Transects_UTM.xlsx','2025_03_14_Transects_UTM.xlsx'};
% diff matrix and rmse
diff_matrix = NaN(211,221,length(epochs)); % this is hard coded
hs_matrix = NaN(211,221,length(epochs)); % this is hard coded
rmse_vals=NaN(length(epochs),1);

for i=1:length(epochs)
    hand_survey_path=append(surveypath,'/',surveynames(i)); % assign hand survey
    epoch=epochs(i); % epoch num
    for j=3:2:size(meanMAPzmatrix,3)
        MAPz=meanMAPzmatrix(:,:,j); % assign map z
    end
    % compare
    [hs_gridmean,hs_gridmed,diff,rmse_val]=GEM_compare(MAPz,epoch,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
    diff_matrix(:,:,i)=diff;
    rmse_vals(i,1)=rmse_val;
    hs_matrix(:,:,i)=hs_gridmean;
end

%% Create Video of Differences 
% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
epochsnum=str2double(epochs);

vname=append('measured_','diff_video');
v = VideoWriter(append(figpath,'/',vname),'MPEG-4');
v.FrameRate = 0.5; 
v.Quality = 100;
open(v)

fig=figure('units','inches','position',[0 0 14 10],'color','w');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
for i=1:length(epochs)
    GEMdate=datetime(epochsnum(i),'ConvertFrom','epochtime','TicksPerSecond',1000); 
    clf
    t=tiledlayout('horizontal');ax1=nexttile;
    %ax1=subplot(1,2,1);
    pcolor(Xgrid,Ygrid,meanMAPzmatrix(:,:,i)); grid off; shading flat;
    hold on; title("Averaged Stereo Elevation Values"); cb1=colorbar(ax1);hold on;cb1.Label.String = 'Elevation (m NAVD83 (2011))';
    hold on; set(gca,'fontsize',14); xlim([0 50]); ylim([-45 20]);
    %ax2=subplot(1,2,2);
    ax2=nexttile;
    pcolor(Xgrid,Ygrid,hs_matrix(:,:,i)); grid off; shading flat; title("Gridded Hand Transect Elevation Values");
    cb2=colorbar(ax2); cb2.Label.String = 'Elevation (m NAVD83 (2011))';
    %linkaxes([ax1 ax2]);sgtitle(append(GEMname,',',GEMdate)); 
    title(t,append(string(epochs(i)),string(GEMdate)),'FontSize',16);ylabel(t,'Alongshore (m)');xlabel(t,'Cross-shore (m)');
    hold on; set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
    ax3=nexttile;
    pcolor(Xgrid,Ygrid,diff_matrix(:,:,i)); grid off; shading flat; title('Stereo minus Transects');
    colormap(ax3,'cool'); cb3=colorbar(ax3); cb3.Label.String = 'Difference in Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
    % Enlarge figure to full screen
    %rmse_txt=num2str(rmse_val); rmse_txt=append('RMSE = ', rmse_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',rmse_txt,'EdgeColor','none','FontSize',28);
    clim(ax1,[0 5]);clim(ax2,[0 5]), clim(ax3,[-0.15 0.15]);
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)


%% OLD____
%{
% survey path
hand_survey_path=append(surveypath,'/2025_01_29_Transects_UTM.xlsx');
% metashape
GEMmat_path=append(metashape_path,'1735848001721/');
figpath=append(metashape_path,'Figures');
[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
% measured
GEMmat_path=append(measured_path,'1735848001721/');
figpath=append(measured_path,'Figures');
[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(GEMmat_path,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
%}
%% Pull a transect from GEMs
% specs
ypick = [7.3 -20]; % (location on grid)
yavg=1.2; % width (m)
savepath=append(metashape_path,'/Transects');
figpath=append(savepath,'/Figures');
[ztran]=transect_pull(GEMpath,numframes,yavg,dxy,ypick,savepath,figpath);

