% Read in Metashape DEMs, Interpolate to LiDAR Mesh using min,mean,and max
% values within a cell, compare Metashape values to LiDAR

% generalpath = 'C:/MASONbeast/';
% datapath = [generalpath,'data/MetashapeTIFs/DEMs/']; % setting paths import/export
% figfolder = [generalpath,'data/Figures/DEMfigs/'];
 fname = '1699992001273v20_DEM.tif';  % 11/16/2023
% fname = '1708099201978v2_DEM.tif';  % 02/15/2024
% fname = '1711479601066_DEM.tif';  % 03/26/2024

% [A,R] = readgeoraster([datapath,fname]); % read in tif file
[A,R] = readgeoraster([fname]);
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A);
%% Now Read in LiDAR
% LidarPath = 'C:/MASONbeast/LIDAR/';
% fname2 = 'mase_lidar_nov152023_15cm.tif';
fname2 = 'mase_021524_15cm_FOV.tif';
% [AL,RL] = readgeoraster([fname2]); % read in tif file  02/15/2024
[AL,RL] = readgeoraster([fname2]);   % 11/16/2023
AL(AL<-7) = NaN; % Create reasonable z limit
XLin = RL.XWorldLimits(1):RL.CellExtentInWorldX:RL.XWorldLimits(2);
YLin = RL.YWorldLimits(1):RL.CellExtentInWorldY:RL.YWorldLimits(2);
[XL,YL] = meshgrid(XLin(1:end-1),YLin(1:end-1)); % develop mesh
AL = flipud(AL);
AL1 = AL;
AL2 = AL;
AL3 = AL;
figname = {'Mean Interp Error', 'Min Interp Error', 'Max Interp Error', 'STD Interp'};
GCPpath = 'C:/Temp/crl1884/WORK/';
%%
 % GCPs = readtable([GCPpath,'20231116.txt']);  % 11/16/2023 gcps 
GCPs = readtable(['GCPs_02_15_2024.txt']);    % 02/15/2024 gcps
GCPs = table2array(GCPs);
GCPeast = GCPs(:,2);
GCPnorth = GCPs(:,3);
GCPZ = GCPs(:,4);
% GCPeast = [239775.56, 239760.33, 239761.18, 239763.20, 239758.85, 239768.07];
% GCPnorth = [3784747.94, 3784729.84, 3784739.58, 3784745.31, 3784749.42, 3784752.25]; % GCPs used 

%% Mean surface
[MeanZ,Meannumpts] = roundgridfun(X,Y,A,XL,YL,@mean);% takes mean of all values in cell 
meanfunctCorrect = find(MeanZ == 0);
MeanZ(meanfunctCorrect) = NaN; % rounding grid function creates false zeros, making them nan
MeanNAN = isnan(MeanZ); % resolving nan values from Metshape DEM
findMeanNAN = find(MeanNAN > 0);
AL1(findMeanNAN)= NaN; 
MeanCOMP = AL1 - MeanZ; % Compare DEMs, take difference in elevation
MeanRMSE = rmse(AL1,MeanZ);
%% Interpolating Minimum values within each cell
[MinZ,Minnumpts] = roundgridfun(X,Y,A,XL,YL,@min);% takes mean of all values in cell 
MinNAN = isnan(MinZ);
findMinNAN = find(MinNAN > 0);
AL2(findMinNAN)= NaN; 
MinCOMP = AL2 - MinZ;
MinRMSE = rmse(AL2,MinZ);
%% Interpolating Maximum values within each cell
[MaxZ,Maxnumpts] = roundgridfun(X,Y,A,XL,YL,@max);% takes mean of all values in cell 
MaxNAN = isnan(MaxZ);
findMaxNAN = find(MaxNAN > 0);
AL3(findMaxNAN)= NaN; 
MaxCOMP = AL3 - MaxZ;
MaxRMSE = rmse(AL3,MaxZ);
%% Standerd deviation of Interpolation
[stdZ,stdnumpts] = roundgridfun(X,Y,A,XL,YL,@std);% takes mean of all values in cell 


%% Plot all figures
figure(2); % Minimum value within each cell 
pcolor(X,Y,A); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('5/7/2024');
ylim([3784680,3784770]);
% plot(XL(2360,:),YL(2360,:),'r-') % Berm Window 
% plot(XL(2250,:),YL(2250,:),'r-')
% plot(XL(:,150),YL(:,150),'r-')
% plot(XL(:,50),YL(:,50),'r-')
hold off;
%%
figure(3); % Mean value within each cell 
pcolor(XL,YL,MeanZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('Mean Z Interp ');
ylim([3784680,3784770]);
hold off;
%%
figure(4); % Maximum value within each cell 
pcolor(XL,YL,MaxZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('Max Z Interp');
ylim([3784680,3784770]);
hold off;
%% 

figure(5); % Comparing Mean Metashape DEM to Lidar DEM
pcolor(XL,YL,MeanCOMP); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-1 1]);
colormap hsv; 
h = colorbar; h.YLabel.String = 'Difference [m]';
ylim([3784680,3784770]);
title(figname(1),'FEB 15');
hold off;
%%

figure(6); % Comparing Minimum Metashape DEM to Lidar DEM
pcolor(XL,YL,MinCOMP); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-1 1]);
colormap hsv;
h = colorbar; h.YLabel.String = 'Difference [m]';
title(figname(2));
ylim([3784680,3784770]);
hold off;


figure(7); % Comparing Maximum Metashape DEM to Lidar DEM
pcolor(XL,YL,MaxCOMP); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-1 1]);
colormap hsv;
h = colorbar; h.YLabel.String = 'Difference [m]';
title(figname(3));
ylim([3784680,3784770]);
hold off;
%% Display standerd deviation 
figure(8);
pcolor(XL,YL,stdZ); hold on;
shading flat;
axis image;
demcmap([-1 1]);
h = colorbar; h.YLabel.String = 'STD';
title(figname(4));
ylim([3784680,3784770]);
hold off;

%% Rotate DEM 
% lims E : 239700-239825       N :   3784600-3784775
% LocalX = XL(1,1); % Local Origin - 
% LocalY = YL(1,1);
LocalX = 239750;  % 11/16/2023
LocalY = 3784710;
Xo = XL - LocalX; % subtract origin
Yo = YL - LocalY; 
theta = 35;
RY = Yo.*cosd(theta) + Xo.*sind(theta);
RX = Xo.*cosd(theta) - Yo.*sind(theta);
XLout = RX + LocalX;
YLout = RY + LocalY;
%%
figure(9); % Display rotated DEM 
pcolor(XLout,YLout,MeanCOMP); hold on;
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation NAVD88 [m]';
title('Rotated DEM');
xlim auto;
% ylim([3784650,3784710]);
% xlim([239550,239610])
% yline(3784680,'r');
% yline(3784730,'r');
xline(239780,'r')
hold off;
%% Seperate Beach from waves 
waves = find(XLout > 239780);
Beach = MeanCOMP;
Beach(waves) = NaN;

figure(40);
pcolor(XLout,YLout,Beach); hold on;
shading flat;
axis image;
demcmap([-0.25 0.25]);
h = colorbar; h.YLabel.String = 'Difference [m]';
title('Lidar - Stereo');
xlim auto;
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)');
xlim([239730 239790])
ylim([3784680 3784760])
hold off;


%%
figure(10); % first display figure of DEM comparison to show where transects are taken from 
pcolor(XL,YL,MeanZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
ylim([3784680,3784770]);
title('Transect Locations');
hold off;


CrossX = 239750:0.15:239820; % Develop X matrix of transects from domain range
CrossY(1,:) = ((-20/24)*CrossX) + 3984570; % Create Line oriented perpendicular to coastline
CrossZ(1,:) = interp2(X,Y,A,CrossX,CrossY(1,:)); % Interp first Z from Metashape
CrossZL(1,:) = interp2(XL,YL,AL,CrossX,CrossY(1,:)); % Interp frist Z from Lidar

% Build off first transect in south direction
for i = 2:10;
CrossY(i,:) = CrossY(i-1,:) - (5);
CrossZ(i,:) = interp2(X,Y,A,CrossX,CrossY(i,:));
CrossZL(i,:) = interp2(XL,YL,AL,CrossX,CrossY(i,:));
end

figure(10); hold on; % Show Transects taken from DEM 
plot(CrossX,CrossY(2,:),'b-');
plot(CrossX,CrossY(3,:),'y-');
plot(CrossX,CrossY(4,:),'g-');
plot(CrossX,CrossY(5,:),'c-');
plot(CrossX,CrossY(6,:),'m-');
plot(CrossX,CrossY(7,:),'k-');
hold off;

figure(12);
plot(CrossZ(2,:),'b','linewidth',2); hold on;
plot(CrossZL(2,:),'r');
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 1'); hold off;

figure(13);
plot(CrossZ(3,:),'y','linewidth',2); hold on;
plot(CrossZL(3,:),'r');
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 2');hold off;

figure(14);
plot(CrossZ(4,:),'g','linewidth',2); hold on;
plot(CrossZL(4,:),'r'); 
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 3'); hold off;

figure(15);
plot(CrossZ(5,:),'c','linewidth',2); hold on;
plot(CrossZL(5,:),'r');
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 4'); hold off;

figure(16);
plot(CrossZ(6,:),'m','linewidth',2); hold on;
plot(CrossZL(6,:),'r'); 
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 5'); hold off;


figure(17);
plot(CrossZ(7,:),'k','linewidth',2); hold on;
plot(CrossZL(7,:),'r'); 
xlabel('Grid Points(15 cm)')
ylabel('Elevation m')
title('Transect 6'); hold off;

%% RMSE / error of transects 
% T1 
Transectrmse(1) = rmse(CrossZL(2,:),CrossZ(2,:));
Transectrmse(2) = rmse(CrossZL(3,:),CrossZ(3,:));
Transectrmse(3) = rmse(CrossZL(4,:),CrossZ(4,:));
Transectrmse(4) = rmse(CrossZL(5,:),CrossZ(5,:));
Transectrmse(5) = rmse(CrossZL(6,:),CrossZ(6,:));
Transectrmse(6) = rmse(CrossZL(7,:),CrossZ(7,:));


%% Compare GCP Locations
for i =1:length(GCPeast)
    MetaGCPz(i) = interp2(X,Y,A,GCPeast(i),GCPnorth(i));
end 
% 
gcpnan = isnan(MetaGCPz);
GCPnan =find(gcpnan==0); 
GCPdiff = GCPZ(GCPnan)'-MetaGCPz(GCPnan);
avediff = mean(GCPdiff);
%%  Compare to previous transects
load('mase_camera_profiles_20231016.mat');
% Convert all transects to UTM for comparison 
[T1(:,1),T1(:,2)] = ll2utm(tran1(:,2),tran1(:,1),18);
T1(:,3) = tran1(:,3);
[T2(:,1),T2(:,2)] = ll2utm(tran2(:,2),tran2(:,1),18);
T2(:,3) = tran2(:,3);
[T3(:,1),T3(:,2)] = ll2utm(tran3(:,2),tran3(:,1),18);
T3(:,3) = tran3(:,3);
[T4(:,1),T4(:,2)] = ll2utm(tran4(:,2),tran4(:,1),18);
T4(:,3) = tran4(:,3);
[T5(:,1),T5(:,2)] = ll2utm(tran5(:,2),tran5(:,1),18);
T5(:,3) = tran5(:,3);
% Problem Point
T2(60,:) = NaN;
%% Interpolate metashape DEM to transect locations 
% MetaTran1Z = interp2(X,Y,A,T1(:,1),T1(:,2));
% MetaTran2Z = interp2(X,Y,A,T2(:,1),T2(:,2));
% MetaTran3Z = interp2(X,Y,A,T3(:,1),T3(:,2));
% MetaTran4Z = interp2(X,Y,A,T4(:,1),T4(:,2));
% MetaTran5Z = interp2(X,Y,A,T5(:,1),T5(:,2));
%%  Instead of interpolating I am taking the nearest point from Meta shape DEM
% create 3 column vectors of XYZ 
sizeAL1 = size(AL1);
sizeAL1 = sizeAL1(1)*sizeAL1(2);
XYZlidar = zeros(sizeAL1,3);
% XYZlidar(:,1) = XL;
% XYZlidar(:,2) = YL;
% XYZlidar(:,3) = AL1;
XYZMeta = zeros(sizeAL1,3);
XYZMeta(:,1) = reshape(XL',[],1);
XYZMeta(:,2) = reshape(YL',[],1);
XYZMeta(:,3) = reshape(MeanZ',[],1);;

k = 1; % Number of nearest neighbors to take 
% [MetaT1Ind,MetaT1dis] = findNearestNeighbors(XYZMeta,T1,k);

% Not working yet
% ************


%% Taking differences of transects, Most only pull from indicies that are not Nan
T1nan = isnan(MetaTran1Z);
T1ind = find(T1nan==0);
T1zcomp = T1(T1ind,:) - MetaTran1Z(T1ind);

T2nan = isnan(MetaTran2Z);
T2ind = find(T2nan==0);
T2zcomp = T2(T2ind,:) - MetaTran2Z(T2ind);

T3nan = isnan(MetaTran3Z);
T3ind = find(T3nan==0);
T3zcomp = T3(T3ind,:) - MetaTran3Z(T3ind);

T4nan = isnan(MetaTran4Z);
T4ind = find(T4nan==0);
T4zcomp = T4(T4ind,:) - MetaTran4Z(T4ind);

T5nan = isnan(MetaTran5Z);
T5ind = find(T5nan==0);
T5zcomp = T5(T5ind,:) - MetaTran5Z(T5ind);

MeanTranEr(1) = mean(T1zcomp(:,3));
MeanTranEr(2) = mean(T2zcomp(:,3));
MeanTranEr(3) = mean(T3zcomp(:,3));
MeanTranEr(4) = mean(T4zcomp(:,3));
MeanTranEr(5) = mean(T5zcomp(:,3));





%% Plot hand transects vs Metashape
figure(18);
plot(MetaTran1Z,'b');
hold on
plot(T1(:,3),'r');
xlabel('Cross-shore')
ylabel('Elevation NAVD88(m)')
title('Transect 1 Comp')
legend('MetaShape','Hand Survey')
hold off

figure(19);
plot(MetaTran2Z,'b');
hold on
plot(T2(:,3),'r');
xlabel('Cross-shore')
ylabel('Elevation NAVD88(m)')
title('Transect 2 Comp')
legend('MetaShape','Hand Survey')
hold off

figure(20);
plot(MetaTran3Z,'b');
hold on
plot(T3(:,3),'r');
xlabel('Cross-shore')
ylabel('Elevation NAVD88(m)')
title('Transect 3 Comp')
legend('MetaShape','Hand Survey')
hold off

figure(21);
plot(MetaTran4Z,'b');
hold on
plot(T4(:,3),'r');
xlabel('Cross-shore')
ylabel('Elevation NAVD88(m)')
title('Transect 4 Comp')
legend('MetaShape','Hand Survey')
hold off

figure(22);
plot(MetaTran5Z,'b');
hold on
plot(T5(:,3),'r');
xlabel('Cross-shore')
ylabel('Elevation NAVD88(m)')
title('Transect 5 Comp')
legend('MetaShape','Hand Survey')
hold off

%% Where are transects taken from? 
figure(23); % first display figure of DEM comparison to show where transects are taken from 
pcolor(XL,YL,MeanZ); hold on;
plot(GCPeast,GCPnorth,'ro')
plot(T1(:,1),T1(:,2),'r')
plot(T2(:,1),T2(:,2),'r')
plot(T3(:,1),T3(:,2),'r')
plot(T4(:,1),T4(:,2),'r')
plot(T5(:,1),T5(:,2),'r')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
ylim([3784680,3784770]);
ylabel('UTM Norhting(m)')
xlabel('UTM Easting(m)')
title('Transect Locations');
hold off;

%% Plot levels of error by indicy (to show that the error is much lower the
% closer to the dune (cameras)
figure(24);
plot(T1ind,T1zcomp(:,3)); hold on;
plot(T2ind,T2zcomp(:,3))
plot(T3ind,T3zcomp(:,3))
plot(T4ind,T4zcomp(:,3))
plot(T5ind,T5zcomp(:,3))
legend('T1','T2','T3','T4','T5')
xlabel('Cross-shore Index of Transect')
ylabel('Difference (m)')
title('Difference: Transect - Meta')
hold off;

%% Plot number of points within each cell 
Meannumpts(findMeanNAN)= NaN; 
figure(25);
pcolor(XL,YL,Meannumpts);
shading flat;
axis image;
demcmap([25 36]);
ylim([3784680,3784770]);
% colormap hsv;
h = colorbar; h.YLabel.String = 'Number of points';
title('Number of Points Within each Cell, Mean Interp');
hold off;

Minnumpts(findMinNAN)= NaN; 
figure(26);
pcolor(XL,YL,Minnumpts);
shading flat;
axis image;
demcmap([25 36]);
ylim([3784680,3784770]);
% colormap hsv;
h = colorbar; h.YLabel.String = 'Number of points';
title('Number of Points Within each Cell, Min Interp');
hold off;

Maxnumpts(findMaxNAN)= NaN; 
figure(27);
pcolor(XL,YL,Maxnumpts);
shading flat;
axis image;
demcmap([25 36]);
ylim([3784680,3784770]);
% colormap hsv;
h = colorbar; h.YLabel.String = 'Number of points';
title('Number of Points Within each Cell, Max Interp');
hold off;

%% 
% Extract Dune aqnd berm and compare rmse
MinDuneInd = find(MinZ>3);
MinBermInd = find(MinZ<3 & MinZ>1.2);
% MinDuneComp = AL2(MinDuneInd) - MinZ(MinDuneInd); % if needed
MinDrmse = rmse(AL2(MinDuneInd),MinZ(MinDuneInd)); % RMSE of Dune
MinBrmse = rmse(AL2(MinBermInd),MinZ(MinBermInd));% RMSE of Berm 
%% 
MeanDuneInd = find(MeanZ>3);
MeanBermInd = find(MeanZ<3 & MeanZ>1.2);
% MeanDuneComp = AL1(MeanDuneInd) - MeanZ(MeanDuneInd); % if needed
MeanDrmse = rmse(AL1(MeanDuneInd),MeanZ(MeanDuneInd)); % RMSE of Dune (using AL1 for mean,AL2-min,AL3-max because I seperated them earlier) 
MeanBrmse = rmse(AL1(MeanBermInd),MeanZ(MeanBermInd));% RMSE of Berm 
%% 
MaxDuneInd = find(MaxZ>3);
MaxBermInd = find(MaxZ<3 & MaxZ>1.2);
% MaxDuneComp = AL3(MaxDuneInd) - MaxZ(MaxDuneInd); % if needed
MaxDrmse = rmse(AL3(MaxDuneInd),MaxZ(MaxDuneInd)); % RMSE of Dune
MaxBrmse = rmse(AL3(MaxBermInd),MaxZ(MaxBermInd));% RMSE of Berm 
%%
figure(28); % Display Berm and Dune Extraction (MINIMUM INTERP)
pcolor(XL,YL,MinZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('Min Dune and Berm Extraction');
ylim([3784680,3784770]);
plot(XL(MinDuneInd),YL(MinDuneInd),'g.');
plot(XL(MinBermInd),YL(MinBermInd),'b.');
legend('DEM','GCPs','Dune','Berm')
hold off;


figure(29); % Display Berm and Dune Extraction (MEAN INTERP)
pcolor(XL,YL,MeanZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('Mean Dune and Berm Extraction');
ylim([3784680,3784770]);
plot(XL(MeanDuneInd),YL(MeanDuneInd),'g.');
plot(XL(MeanBermInd),YL(MeanBermInd),'b.');
legend('DEM','GCPs','Dune','Berm')
hold off;

figure(30); % Display Berm and Dune Extraction (MAX INTERP)
pcolor(XL,YL,MaxZ); hold on;
plot(GCPeast,GCPnorth,'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('Max Dune and Berm Extraction');
ylim([3784680,3784770]);
plot(XL(MaxDuneInd),YL(MaxDuneInd),'g.');
plot(XL(MaxBermInd),YL(MaxBermInd),'b.');
legend('DEM','GCPs','Dune','Berm')
hold off;


%%  seperate north and south portions 
% NorthPfind = find(Yout>3784680); 02/15/2024
% SouthPfind = find(Yout<3784680);
% AL1Berm = AL1(MeanBermInd);
% MeanZBerm = MeanZ(MeanBermInd);
% AL1Dune = AL1(MeanDuneInd);
% MeanZDune = MeanZ(MeanDuneInd);

NorthPfind = find(Yout>3784730); % 11/16/2023
SouthPfind = find(Yout<3784730);
NorthPcomp = AL1(NorthPfind)- MeanZ(NorthPfind);
SouthPcomp = AL1(SouthPfind) - MeanZ(SouthPfind);
% Full Portions RMSE
NorthPrmse = rmse(AL1(NorthPfind),MeanZ(NorthPfind));
SouthPrmse = rmse(AL1(SouthPfind),MeanZ(SouthPfind));


% % North and South of BERM vs Dune
% NorthPrmseBERM = rmse(AL1Berm(NorthPfind),MeanZBerm(NorthPfind));
% SouthPrmseBERM = rmse(AL1Berm(SouthPfind),MeanZBerm(SouthPfind));
% NorthPrmseDune = rmse(AL1Dune(NorthPfind),MeanZDune(NorthPfind));
% SouthPrmseDune = rmse(AL1Dune(SouthPfind),MeanZDune(SouthPfind));


% % Show portions on figure 
% figure(31);
% plot(XL(NorthPfind),YL(NorthPfind),'r');
% plot(XL(SouthPfind),YL(SouthPfind),'b');
% legend('North Portion','South Portion');
% title('Display North vs South used in Comp');
% xlabel('UTM Easting(m)');
% ylabel('UTM Northing(m)');

%% 03/26/2024 Transect Comp 

GCPs = readtable('GCPs_03_26_2024.txt');    % 02/15/2024 gcps
GCPs = table2array(GCPs);
GCPeast = GCPs(:,2);
GCPnorth = GCPs(:,3);
GCPZ = GCPs(:,4);

Tran1 = readtable('C:\Temp\crl1884\WORK\03262024Transects.xlsx');
Tran = table2array(Tran1(2:end,:));
XT = Tran(:,1);
YT = Tran(:,2);
ZT = Tran(:,3);

figure(32);
pcolor(X,Y,A); hold on;
plot(GCPeast(1:6),GCPnorth(1:6),'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
title('3/26/24 Metashape DEM');
ylim([3784650,3784770]);
plot(XT,YT,'r')
hold off;




%% Rotate MetaShape DEM
LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
Xo = X - LocalX; % subtract origin
Yo = Y - LocalY; 
theta = 35;
RY = Yo.*cosd(theta) + Xo.*sind(theta);
RX = Xo.*cosd(theta) - Yo.*sind(theta);
Xout = RX + LocalX;
Yout = RY + LocalY;
% Rotate Transects too!
XoT = XT - LocalX; % subtract origin
YoT = YT - LocalY; 
theta = 35;
RYT = YoT.*cosd(theta) + XoT.*sind(theta);
RXT = XoT.*cosd(theta) - YoT.*sind(theta);
XTout = RXT + LocalX;
YTout = RYT + LocalY;
%%
figure(33); % Display rotated DEM 
pcolor(Xout,Yout,A); hold on;
plot(XTout,YTout,'r')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation NAVD88 [m]';
title('11/16/2023','DEM');
xlim auto;
% ylim([3784650,3784710]);
% xlim([239550,239610])
hold off;
%% Now Compare DEM to Transects
META(:,1) = X(:);
META(:,2) = Y(:);
META(:,3) = A(:);

% Initialize variables to store results
numPoints = size(Tran, 1);
averagedZValues = NaN(numPoints, 1);

% Define threshold for proximity in Y direction (15 cm)
thresholdY = 0.15; % in meters

% Iterate through each point in Transect
for i = 1:numPoints
    % Extract Y-coordinate of the current Transect point
    yCoord = Tran(i, 2);
    
    % Find DEM points within threshold distance of the current Y-coordinate
    withinThresholdIndices = find(abs(META(:, 2) - yCoord) <= thresholdY);
    
    % Remove NaN values from the DEM data
    demPoints = META(withinThresholdIndices, :);
    demPoints = demPoints(~isnan(demPoints(:, 3)), :);
    
    % Calculate average Z value if there are points within threshold distance
    if ~isempty(demPoints)
        averagedZValues(i) = mean(demPoints(:, 3));
    end
end
%%
TranRMSE = rmse(averagedZValues,ZT);
figure(34);
plot(Tran(:,3),'r')
hold on
plot(averagedZValues,'b')
title('Compare Walking Transects to Averaged Z values within 15cm')
legend('Walking Transects', 'MetaShape')
ylabel('Elevation NAVD88 (m)')
hold off
%%
Zdiff = Tran(:,3) - averagedZValues;
figure(35);
plot(Zdiff,'r')
title('Hand Surveyed Z values - Avg. Meta Z val within 15 cm')
ylabel('Error (m)')
xlabel('Idicies')
yline(0)
%%
figure(36)
subplot(2,1,1)
plot(Tran(:,3),'r')
hold on
plot(averagedZValues,'b')
title('Compare Walking Transects to Averaged Z values within 15cm')
legend('Walking Transects', 'MetaShape')
ylabel('Elevation NAVD88 (m)')
xlim([1,550])
hold off
subplot(2,1,2)
plot(Zdiff,'r')
title('Hand Surveyed Z values - Avg. Meta Z val within 15 cm')
ylabel('Error (m)')
xlabel('Index')
yline(0)
xlim([1,550])
%% 


%% Extract each transect 
T1 = 1:80;
T2 = 92:145;
T3 = 158:208;
T4 = 220:266;
T5 = 273:325;
T6 = 347:380;
T7 = 403:440;
T8 = 456:494;
T9 = 504:545;
TransectRMSE(1) = rmse(ZT(T1),averagedZValues(T1));
TransectRMSE(2) = rmse(ZT(T2),averagedZValues(T2));
TransectRMSE(3) = rmse(ZT(T3),averagedZValues(T3));
TransectRMSE(4) = rmse(ZT(T4),averagedZValues(T4));
TransectRMSE(5) = rmse(ZT(T5),averagedZValues(T5));
TransectRMSE(6) = rmse(ZT(T6),averagedZValues(T6));
TransectRMSE(7) = rmse(ZT(T7),averagedZValues(T7));
TransectRMSE(8) = rmse(ZT(T8),averagedZValues(T8));
TransectRMSE(9) = rmse(ZT(T9),averagedZValues(T9));


%% 
figure(37);
plot(ZT(T1),'r'); hold on
plot(averagedZValues(T1),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 1')
xlabel('Index')
hold off

figure(38);
plot(ZT(T2),'r'); hold on
plot(averagedZValues(T2),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 2')
xlabel('Index')
set ( gca, 'xdir', 'reverse' )
hold off

figure(39);
plot(ZT(T3),'r'); hold on
plot(averagedZValues(T3),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 3')
xlabel('Index')
hold off

figure(40);
plot(ZT(T4),'r'); hold on
plot(averagedZValues(T4),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 4')
xlabel('Index')
set ( gca, 'xdir', 'reverse' )
hold off

figure(41);
plot(ZT(T5),'r'); hold on
plot(averagedZValues(T5),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 5')
xlabel('Index')
hold off

figure(42);
plot(ZT(T6),'r'); hold on
plot(averagedZValues(T6),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 6')
xlabel('Index')
set ( gca, 'xdir', 'reverse' )
hold off

figure(43);
plot(ZT(T7),'r'); hold on
plot(averagedZValues(T7),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 7')
xlabel('Index')
hold off

figure(44);
plot(ZT(T8),'r'); hold on
plot(averagedZValues(T8),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 8')
xlabel('Index')
set ( gca, 'xdir', 'reverse' )
hold off

figure(45);
plot(ZT(T9),'r'); hold on
plot(averagedZValues(T9),'b');
legend('Walking Transects', 'MetaShape');
ylabel('Elevation NAVD88 (m)');
title('Transect 9')
xlabel('Index')
hold off







% make DEM into column vectors
% Xmeta = X(:);
% Ymeta = Y(:);
% Zmeta = A(:);
% PwithinThreshind = [];
% for i = 1:545;
%     Thresh = 0.25;
%      dis = abs(Yout - YTout(i));
%     PwithinThreshind = find(dis <= deltaY);
% if ~isempty(PwithinThreshind)
%     ZvalDEM = Zmeta(PwithinThreshind);
%     averageDEMZ = mean(ZvalDEM);
%     MetaTranZ(i) = averageDEMZ;
% end
% clear dis 
% clear PwithinThreshind
% clear ZvalDEM
% clear averageDEMZ
%    
%  end
    
%   MetaZ(i) = mean(A(findmindis));















%%
Xmeta = X(:);
Ymeta = Y(:);
TrancompXY(:,1) = Tran(:,1);
TrancompXY(:,2) = Tran(:,2);
T = delaunayn([Xmeta,Ymeta]);
k = dsearchn([Xmeta,Ymeta],T,TrancompXY);

MetaTran(:,1) = X(k);
MetaTran(:,2) = Y(k);
MetaTran(:,3) = A(k);
% dist  =  sqrt( (X(k)-).^2 + (yy(k)-y_data).^2 );



% [xx,yy] = meshgrid([0:10],[20:30]');
% data = [5 , 28; 8, 29.5];
% T = delaunayn([xx(:), yy(:)]);
% k = dsearchn([xx(:) yy(:)], T, data);
% [xx(k), yy(k)]
% 
% error = z_data - z_dem(nearest);
% also look at 
% dist  =  sqrt( (xx(k)-x_data).^2 + (yy(k)-y_data).^2 );
% 
% if you want to evaluate the slope contribution. First, focus on cross-shore slope...

% error_slope = delta_x * local_beach_slope
%                          = (x_data(2:end-1)-xx(k(2:end-1))) * (z_data(3:end)-z_data(1:end-2)) / ( x_data(3:end) - x_data(1:end-2)  )

%% Curious if I try interp instead of swath averaged

MetaTranInterp =  interp2(X,Y,A,XT,YT); % Interp first Z from Metashape
Tnan = isnan(MetaTranInterp);
Tind = find(Tnan==0);
Tzcomp = ZT(Tind) - MetaTranInterp(Tind);
InterpTranRMSE = rmse(ZT(Tind),MetaTranInterp(Tind));

[MetaTranInterp2, numpts] = roundgridfun(X,Y,A,XT,YT,@mean);


%% May Pres figs

figure(45); % Comparing Mean Metashape DEM to Lidar DEM
pcolor(X,Y,A); hold on;
plot(GCPeast,GCPnorth,'ro')
plot(cameast,camnorth,'bo')
shading flat;
axis image;
demcmap([-2 5]);
legend('','Control Points','Cameras')
% colormap hsv; 
h = colorbar; h.YLabel.String = 'Elevation (m), NAVD88';
ylim([3784680,3784770]);
title('November 16th, 2023','DEM');
hold off;

