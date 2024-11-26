clear;
clc;
%% Metashape DEM Input   (STEP 1) 
fname = '1717167601959_DEM.tif';  % define tif file used
[A,R] = readgeoraster([fname]); % read in tif file 
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A); % When tif files are read in this way they need to be flipped

% Transect Input (STEP 1)
Tran1 = readtable('C:\Temp\crl1884\WORK\05312024_trans_UTM.xlsx'); % my current file structure for transects
Tran = table2array(Tran1(2:end,:));
XT = Tran(:,3);
YT = Tran(:,2);
ZT = Tran(:,4); % make sure these are the correct x,y,z values (the excel sheets can vary for transect data)

% GCP Input (STEP 1)
GCPs = readtable(['GCPs_02_15_2024.txt']);    % 02/15/2024 : Define correct GCPs
GCPs = table2array(GCPs);
GCPeast = GCPs(:,2);
GCPnorth = GCPs(:,3);
GCPZ = GCPs(:,4); % I seperate the x, y and z but is not required

%% Rotate X and Y so that the grid can be oriented in the cross and long shore 
% (STEP 2)
LocalX = 239740;  % 3/26/2024
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

%% Display rotated DEM (STEP 3)
Q = 1; % query point I use to plot over the DEM and transects (when hardcoding 
% to find start/end of each transect 
figure(1); 
pcolor(Xout,Yout,A); hold on; % MetaShape DEM
plot(XTout,YTout,'r')  % Transects
plot(XTout(Q),YTout(Q),'go') % Query Point
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation NAVD88 [m]';
title('3/26/2024','DEM');
xlim auto;
hold off;
%% 5/31 Transect extraction  (STEP 3)
% When finding unique transects these indicies will change 
T1 = 1:135;
T2 = 145:270;
T3 = 280:390;
T4 = 400:505;
T5 = 525:615;
T6 = 630:700;
T7 = 725:782;
T8 = 800:863;

%% Generate grid oriented in cross/longshore to compare transects to Metashape (STEP 4)
gridX = 239695:1:239750; % 1 m resolution
gridY = 3784715:1:3784770;
[Gx,Gy] = meshgrid(gridX,gridY);

%% Rotate it back ** OPTIONAL**
% LocalX = 239740;  % 11/16/2023
% LocalY = 3784650;
% XoGr = Gx - LocalX; % subtract origin
% YoGr = Gy - LocalY; 
% RYGr = YoGr.*cosd(-theta) + XoGr.*sind(-theta);
% RXGr = XoGr.*cosd(-theta) - YoGr.*sind(-theta);
% TgridX = RXGr + LocalX;
% TgridY = RYGr + LocalY;

%% Regrid Stereo to comp grid (STEP 5)
[MeanZ,Meannumpts] = roundgridfun(Xout,Yout,A,Gx,Gy,@mean);% takes mean of all values in cell 
MeanZ(MeanZ == 0) = NaN; % rounding grid function creates false zeros, making them nan
% regrid transects to the same grid above (STEP 5)    
[TZ,Tnumpts] = roundgridfun(XTout,YTout,ZT,Gx,Gy,@mean);% takes mean of all values in cell 
TZ(TZ == 0) = NaN;

%%
Qi = 35;
figure(3);
pcolor(Gx,Gy,MeanZ);
hold on
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('3/26/2024',' DEM');
plot(XTout(T1),YTout(T1),'r')  % T1
plot(Gx(Qi,:),Gy(Qi,:),'k')
% plot(Gx(228,:),Gy(228,:),'k')

%% plot all transects vs Stereo 3/26
figx = 1:1:56; % based on grid size 
figure(4); % T1
plot(figx,TZ(34,:),'b');
hold on
plot(figx,MeanZ(34,:),'r');
legend('Transect 1', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T1 Comp')
hold off

% here is a Katherine edit




