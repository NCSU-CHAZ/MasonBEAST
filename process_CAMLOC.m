%% Script to process monthly GEMs (camera location comparison)

% paths 
genpath='/Volumes/kanarde-1/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
pcpath=append(genpath,'/PointClouds/'); % path to pointclouds
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % save measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % save metashape GEMs (using an input reference point)
forced_ext_path=append(CAM_analysispath,'/Forced/'); % saved forced extrinsics

% GEM specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.5; % meter change to 0.2v normally
numframes=2;

% Loop through Metashape pointclouds
spec='*meta_ptcld*';
figpath=append(metashape_path,'Figures');
[meta_meanMAPzmatrix,meta_medMAPzmatrix,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,metashape_path,figpath,spec);


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
[ref_meanMAPzmatrix,ref_medMAPzmatrix]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,measured_path,figpath,spec);

% Loop through Forced pointclouds
spec='*precal*';
figpath=append(forced_ext_path,'Figures');
[precal_meanMAPzmatrix,precal_medMAPzmatrix]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,forced_ext_path,figpath,spec);

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

% for quick comp of March 24 and 25
%epochs={'1711479601066','1741978801668'};
%surveynames={'2024_03_26_Transects_UTM.xlsx','2025_03_14_Transects_UTM.xlsx'};

% diff matrix and rmse
precal_diff_matrix = NaN(211,221,length(epochs)); % this is hard coded
precal_hs_matrix = NaN(211,221,length(epochs)); % this is hard coded
precal_rmse_vals=NaN(length(epochs),1);

meanMAPzmatrix=precal_meanMAPzmatrix;
% this is hard coded and needs to be fixed
precal_meanMAPzmatrix_oneframe=cat(3,meanMAPzmatrix(:,:,1),meanMAPzmatrix(:,:,3),meanMAPzmatrix(:,:,5),meanMAPzmatrix(:,:,7),meanMAPzmatrix(:,:,9),meanMAPzmatrix(:,:,11),meanMAPzmatrix(:,:,13),meanMAPzmatrix(:,:,15),meanMAPzmatrix(:,:,17),meanMAPzmatrix(:,:,19),meanMAPzmatrix(:,:,21),meanMAPzmatrix(:,:,23),meanMAPzmatrix(:,:,25));
%meanMAPzmatrix_oneframe=cat(3,meanMAPzmatrix(:,:,1),meanMAPzmatrix(:,:,3));

for i=1:length(epochs)
    hand_survey_path=append(surveypath,'/',surveynames(i)); % assign hand survey
    epoch=epochs(i); % epoch num
    MAPz=precal_meanMAPzmatrix_oneframe(:,:,i);
    disp(string(hand_survey_path));
    disp(string(epoch));
    %for j=3:2:size(meanMAPzmatrix,3)
        %MAPz=meanMAPzmatrix(:,:,j); % assign map z
    %end
    % compare
    [hs_gridmean,hs_gridmed,diff,rmse_val]=GEM_compare(MAPz,epoch,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
    precal_diff_matrix(:,:,i)=diff;
    precal_rmse_vals(i,1)=rmse_val;
    precal_hs_matrix(:,:,i)=hs_gridmean;
end

% plot RMSE vals

dates_label ={'Mar 24', 'April 24', 'May 24','June 24', 'July 24', 'Aug 24', 'Sept 24', 'Oct 24' ,'Nov 24', 'Dec 24', 'Jan 25', 'Feb 25', 'Mar 25'};
fig=figure(1);
clf
plot(1:13,precal_rmse_vals,'b-o','LineWidth',4); ylim([0 0.35]);
ylabel('RMSE (m)');xlabel('Date'); xticks(1:13);xticklabels(dates_label); set(gca, 'FontSize', 16);

%% Create Video of Differences 
% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
epochsnum=str2double(epochs);

vname=append('forced_ext_','diff_video');
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
    pcolor(Xgrid,Ygrid,precal_meanMAPzmatrix_oneframe(:,:,i)); grid off; shading flat;
    hold on; title("Averaged Stereo Elevation Values"); cb1=colorbar(ax1);hold on;cb1.Label.String = 'Elevation (m NAVD83 (2011))';
    hold on; set(gca,'fontsize',14); xlim([0 50]); ylim([-45 20]);
    %ax2=subplot(1,2,2);
    ax2=nexttile;
    pcolor(Xgrid,Ygrid,precal_hs_matrix(:,:,i)); grid off; shading flat; title("Gridded Hand Transect Elevation Values");
    cb2=colorbar(ax2); cb2.Label.String = 'Elevation (m NAVD83 (2011))';
    %linkaxes([ax1 ax2]);sgtitle(append(GEMname,',',GEMdate)); 
    title(t,append(string(epochs(i)),string(GEMdate)),'FontSize',16);ylabel(t,'Alongshore (m)');xlabel(t,'Cross-shore (m)');
    hold on; set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
    ax3=nexttile;
    pcolor(Xgrid,Ygrid,precal_diff_matrix(:,:,i)); grid off; shading flat; title('Stereo minus Transects');
    colormap(ax3,'cool'); cb3=colorbar(ax3); cb3.Label.String = 'Difference in Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
    % Enlarge figure to full screen
    %rmse_txt=num2str(rmse_val); rmse_txt=append('RMSE = ', rmse_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',rmse_txt,'EdgeColor','none','FontSize',28);
    clim(ax1,[0 5]);clim(ax2,[0 5]), clim(ax3,[-0.15 0.15]);
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)

%% Forced Extrinsics (OLD)
% Loop through Measured pointclouds
% spec='*1726772401770_meas_ptcld*';
% figpath=append(measured_path,'Figures');
% [sept_meanMAPzmatrix,medMAPzmatrix]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,measured_path,figpath,spec);
% 
% spec='*1711479601066_meas_ptcld*';
% figpath=append(measured_path,'Figures');
% [march24_meanMAPzmatrix,medMAPzmatrix]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,pcpath,measured_path,figpath,spec);
% 
% % combine matrices 
% meanMAPzmatrix=cat(3,sept_meanMAPzmatrix(:,:,1),sept_meanMAPzmatrix(:,:,3),march24_meanMAPzmatrix(:,:,1));
% 
% % GEM comparison manually
% HS_savepath=append(genpath,'/Surveys/rotated_surveys');
% surveypath=append(CAM_analysispath,'/Surveys');
% 
% % epochs
% epochs={'1726772401770','1726772401770','1711479601066'};
% surveynames={'2024_09_18_Transects_UTM.xlsx','2024_09_18_Transects_UTM.xlsx','2024_03_26_Transects_UTM.xlsx'};
% % diff matrix and rmse
% diff_matrix = NaN(211,221,length(epochs)); % this is hard coded
% hs_matrix = NaN(211,221,length(epochs)); % this is hard coded
% rmse_vals=NaN(length(epochs),1);
% 
% for i=1:length(epochs)
%     hand_survey_path=append(surveypath,'/',surveynames(i)); % assign hand survey
%     epoch=epochs(i); % epoch num
%     disp(string(epoch));
%     disp(string(hand_survey_path));
%     %j=0:2:size(meanMAPzmatrix,3)-1;
%     %k=i+j;
%     MAPz=meanMAPzmatrix(:,:,i); % assign map z
%     %disp(string(j));
%     %disp(string(k));
% 
%     % compare
%     [hs_gridmean,hs_gridmed,diff,rmse_val]=GEM_compare(MAPz,epoch,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath);
%     diff_matrix(:,:,i)=diff;
%     rmse_vals(i,1)=rmse_val;
%     hs_matrix(:,:,i)=hs_gridmean;
% end

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


%% Plot Camera Locations per epoch to show movement

%% Yaw, Pitch, and Roll over time
yaw_A=[167.498,167.986,167.668,167.753,167.635,167.724,167.683,167.88,167.466,167.501,167.99,167.541,168.157];
yaw_B=[145.95,145.09,145.055,146.117,146.103,146.265,145.712,145.817,145.501,145.299,145.449,145.51,145.341];

pitch_A=[61.522,61.278,61.535,61.486,61.327,61.126,61.429,61.45,61.212,61.446,61.567,61.541,61.49];
pitch_B=[62.552,62.496,62.625,62.585,62.46,62.248,62.628,62.626,62.338,62.584,62.684,62.706,62.679];

roll_A=[-0.134,0.294,0.064,0.065,0.055,-0.025,0.037,0.249,-0.147,-0.061,0.267,0.515,0.465];
roll_B=[-3.566,-4.526,-3.463,-3.516,-3.394,-3.353,-3.891,-3.689,-4.076,-4.207,-4.191,-3.995,-4.274];

fig=figure;
clf
subplot(3,2,1); plot(yaw_A,'-o','Linewidth',2,'Color','m'); hold on; ylabel('Yaw (degrees)'); title ('Cam A');xticks(1:13);
set(gca,"XTickLabel",[]);set(gca,'FontSize',16);ylim([167.4 169.4]);
subplot(3,2,2);plot(yaw_B,'-o','Linewidth',2,'Color','c'); hold on; title('Cam B'); xticks(1:13); set(gca,"XTickLabel",[]);set(gca,'FontSize',16);ylim([145 147]);
subplot(3,2,3);plot(pitch_A,'-o','Linewidth',2,'Color','m'); hold on; ylabel('Pitch (degrees)');xticks(1:13); set(gca,"XTickLabel",[]);set(gca,'FontSize',16);ylim([61.1 63.1]);
subplot(3,2,4);plot(pitch_B,'-o','Linewidth',2,'Color','c'); hold on; xticks(1:13); set(gca,"XTickLabel",[]);set(gca,'FontSize',16);ylim([62.2 64.2]);
subplot(3,2,5);plot(roll_A,'-o','Linewidth',2,'Color','m'); hold on; ylabel('Roll (degrees)');xticks(1:13);xticklabels(dates_label); set(gca,'FontSize',16);ylim([-0.2 1.8]);
subplot(3,2,6);plot(roll_B,'-o','Linewidth',2,'Color','c'); hold on; xticks(1:13); xticklabels(dates_label);ylim([-4.4 -2.4]);
set(gca,'FontSize',16);





