function[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(MAPz,epoch,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath)
% function[handsurvey_grid_mean,handsurvey_grid_med]=GEM_compare(MAPz,epoch,camlocA,camlocB,dxy,hand_survey_path,HS_savepath,figpath)
% --------------------------------------------------------------------------
% This function takes a mat file of a griided GEM (has gone through
% ptcld_to_GEM) and compares it's values to a hand survey. This is done by
% puting the hand survey on the same grid. The survey serves as the true
% value and RMSE is calculated to show how the GEM compares to the hand
% survey. 
% --------------------------------------------------------------------------
% INPUTS:
% -------
% 
% MAPz = matrix of elevation values from stereo map
% epoch = epoch # of map
% camlocA = camera A location coordinates [Ax,Ay]
% camlocB = camera B location coordinates [Bx,By]
% dxy = grid bin size (m)
%
% hand_survey_path = path to hand survey folder
%%%%% Hand surveys need to be in the order of (Elevation, Northing,
%%%%% Easting)
% GCPpath = path to GCP location (.txt file) THIS IS NOT ADDED IN YET BUT
% EVENTUALLY
% -GCP file must have elevation values in 4th column
% HS_savepath = path to save gridded hand survey mat file
% figpath = path to save figure
% 
% OUTPUTS:
% --------
% handsurvey_grid_mean = gridded mean hand survey as a mat file
% handsurvey_grid_med = gridded median hand survey as a mat file
%
% EDITS NEEDED:
% -------------------------------------------
% epoch to time conversion is wrong
% 
% LAST UPDATED: BG 11/11/2025
%
%
% EXAMPLE:
% -------------------------------------------
% GEMmat_path = '/Users/bagaenzl/Desktop/MasonBEAST Data/GEM_mat_files/1702827001820/';
% camlocA = [239766.1, 3784761.9];%, 10.37
% camlocB = [239759.4, 3784755.0];%, 10.26];
% dxy = 0.5; % meter
% hand_survey_path = '/Users/bagaenzl/Desktop/MasonBEAST Data/surveys/2024_06_26_Transects_UTM.xlsx';
% HS_savepath = '/Users/bagaenzl/Desktop/MasonBEAST Data/surveys/rotated_surveys';
% figpath = '/Users/bagaenzl/Desktop/MasonBEAST Data/Figures';


format long g

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

%-------------------------------------------------------------------------------
%{ 
OLD FOR LOADING GEM
% Load in GEM med and mean
medGEMz=load(fullfile(GEMmat_path,'medGEMz_1.mat'));
medGEMz=[medGEMz.medGEMz];
meanGEMz=load(fullfile(GEMmat_path,'meanGEMz_1.mat'));
meanGEMz=[meanGEMz.meanGEMz];

% GEM name and date
GEMname=split(GEMmat_path,'/');
GEMname=GEMname(9,1);
GEMname=string(GEMname);

GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMdate=string(GEMdate);
GEMtitle=append(GEMname,',',GEMdate);
%}
epoch=epoch{1};
GEMdate=datetime(epoch,'ConvertFrom','epochtime','TicksPerSecond',1000); % TIME IS WRONG BY LIKE A DAY
GEMtitle=append(string(epoch),',',string(GEMdate));
%--------------------------------------------------------------------------
% read in hand surveys
    hand_survey_path=hand_survey_path{1};
    handtran=readmatrix(hand_survey_path);
    Xtran=handtran(:,3); % Easting
    Ytran=handtran(:,2); % Northing
    Ztran=handtran(:,1); % Elevation (meters)

    % Hand survey date
    HSdate=split(hand_survey_path,'/');
    HSdate=HSdate(9,1);
    HSdate=split(HSdate,'_');
    HSdate=append(HSdate(1,1),HSdate(2,1),HSdate(3,1));
    HSdate=string(HSdate);

    % Rotate Hand Survey 
    % coordinate rotation and transformation x,y to cross- and alongshore
    [Xrottran, Yrottran] = rotateCoordinates(Xtran,Ytran, Xloc, Yloc, rotang);

    % Calculate average Z value for hand surveys
    [ZtranMean,TranNumpts]=roundgridfun(Xrottran,Yrottran,Ztran,Xgrid,Ygrid,@mean);
    ZtranMean(ZtranMean==0)=NaN;
    handsurvey_grid_mean=ZtranMean;
    HSGmean_filepath=fullfile(HS_savepath,[append(HSdate,'mean','.mat')]);
    save(HSGmean_filepath,"handsurvey_grid_mean") % save mean z values from hand survey as mat file

    % Calculate the median Z value for hand surveys
    [ZtranMed,TranNumpts]=roundgridfun(Xrottran,Yrottran,Ztran,Xgrid,Ygrid,@median);
    ZtranMed(ZtranMed==0)=NaN;
    handsurvey_grid_med=ZtranMed;
    HSGmed_filepath=fullfile(HS_savepath,[append(HSdate,'med','.mat')]);
    save(HSGmed_filepath,"handsurvey_grid_med") % save mean z values from hand survey as mat file


% --------------------------------------------------------------------------
% Calculate RMSE for GEM vs Hand Survey
rmse_array=rmse(MAPz,handsurvey_grid_mean,'omitnan');
rmse_val=mean(rmse_array,'omitnan');

% PLOT MAP, SURVEY, AND DIFF
% --------------------------------------------------------------------------
compfigpath=append(figpath,"/comp_",GEMname,'_',string(i));
fig=figure;
t=tiledlayout('horizontal');ax1=nexttile;
    %ax1=subplot(1,2,1);
pcolor(Xgrid,Ygrid,MAPz); grid off; shading flat;
hold on; title("Averaged Stereo Elevation Values"); cb1=colorbar(ax1);hold on;cb1.Label.String = 'Elevation (m NAVD83 (2011))';
hold on; set(gca,'fontsize',14); xlim([0 50]); ylim([-45 20]);
    %ax2=subplot(1,2,2);
ax2=nexttile;
pcolor(Xgrid,Ygrid,ZtranMean); grid off; shading flat; title("Gridded Hand Transect Elevation Values");
cb2=colorbar(ax2); cb2.Label.String = 'Elevation (m NAVD83 (2011))';
    %linkaxes([ax1 ax2]);sgtitle(append(GEMname,',',GEMdate)); 
title(t,append(GEMname,',',GEMdate));ylabel(t,'Alongshore (m)');xlabel(t,'Cross-shore (m)');
hold on; set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
ax3=nexttile;
pcolor(Xgrid,Ygrid,MAPz-ZtranMean); grid off; shading flat; title('Stereo minus Transects');
colormap(ax3,'hot'); cb3=colorbar(ax3); cb3.Label.String = 'Difference in Elevation (m NAVD83 (2011))';
set(gca,'fontsize',14);xlim([0 50]); ylim([-45 20]);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
    %rmse_txt=num2str(rmse_mean(i)); rmse_txt=append('RMSE = ', rmse_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',rmse_txt,'EdgeColor','none','FontSize',28);
clim(ax1,[0 5]);clim(ax2,[0 5]), clim(ax3,[-0.15 0.15]);
saveas(t,compfigpath,'png');
close(fig);

% ---------------------------------------------------------------------------
%{
OLD FOR MED OR MEAN
numframes=2;
for i=1:numframes
    rmse_array=rmse(MAPz,handsurvey_grid_med,'omitnan'); % forecasted, observed, omitnans
    rmse_med=mean(rmse_array,'omitnan');
end

% mean
for i=1:numframes
    meanGEMz_r=meanGEMz;
    rmse_array=rmse(meanGEMz_r,handsurvey_grid_mean,'omitnan'); % forecasted, observed, omitnans
    rmse_mean(i)=mean(rmse_array,'omitnan');
end

% Plot and save figure
for i=1:numframes
    % Mean
    meanfigpath=append(figpath,"/mean_comp_",GEMname,'_',string(i));
    fig=figure;
    t=tiledlayout('horizontal');nexttile;
    %ax1=subplot(1,2,1);
    pcolor(Xgrid,Ygrid,meanGEMz); grid off; shading flat;
    hold on; title("Averaged GEM Elevation Values"); colorbar;caxis([0 5]);c1=clim; hold on;
    %ax2=subplot(1,2,2);
    nexttile;
    pcolor(Xgrid,Ygrid,ZtranMean); grid off; shading flat; title("Gridded Hand Transect Elevation Values");
    a=colorbar();caxis([0 5]);c2 = clim; a.Label.String = 'Elevation (m NAD83 (2011))';
    %linkaxes([ax1 ax2]);sgtitle(append(GEMname,',',GEMdate)); 
    title(t,append(GEMname,',',GEMdate));ylabel(t,'Alongshore (m)');xlabel(t,'Cross-shore (m)');
    set(gca,'fontsize',14);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
    %rmse_txt=num2str(rmse_mean(i)); rmse_txt=append('RMSE = ', rmse_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',rmse_txt,'EdgeColor','none','FontSize',28);
    saveas(t,meanfigpath,'png');
    close(fig);

    % Median
    medfigpath=append(figpath,"/med_comp_",GEMname,'_',string(i));
    fig=figure;
    t=tiledlayout('horizontal');nexttile;
    % ax1=subplot(1,2,1);
    pcolor(Xgrid,Ygrid,medGEMz(:,:,i)); grid off; shading flat;
    hold on; title("Median GEM Elevation Values"); colorbar;caxis([0 5]);c1=clim; hold on;
    %ax2=subplot(1,2,2);
    nexttile;
    pcolor(Xgrid,Ygrid,ZtranMean); grid off; shading flat; title("Gridded Hand Transect Elevation Values");
    a=colorbar();caxis([0 5]);c2 = clim; a.Label.String = 'Elevation (m NAD83 (2011))';
    %linkaxes([ax1 ax2]);sgtitle(append(GEMname,',',GEMdate));
    title(t,append(GEMname,',',GEMdate));ylabel(t,'Alongshore (m)');xlabel(t,'Cross-shore (m)');
    set(gca,'fontsize',14);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);% Enlarge figure to full screen
    %rmse_txt=num2str(rmse_med(i)); rmse_txt=append('RMSE = ', rmse_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',rmse_txt,'EdgeColor','none','FontSize',28);
    saveas(t,medfigpath,'png');
    close(fig);
end

%}