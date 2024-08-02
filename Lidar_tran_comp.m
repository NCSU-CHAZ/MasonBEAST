fname = '1708099201978v2_DEM.tif';  % 02/15/2024
[A,R] = readgeoraster([fname]);
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A);
%% read in Lidar
fname2 = 'mase_021524_15cm_FOV.tif';
[AL,RL] = readgeoraster([fname2]);   % 02/15/24
AL(AL<-7) = NaN; % Create reasonable z limit
XLin = RL.XWorldLimits(1):RL.CellExtentInWorldX:RL.XWorldLimits(2);
YLin = RL.YWorldLimits(1):RL.CellExtentInWorldY:RL.YWorldLimits(2);
[XL,YL] = meshgrid(XLin(1:end-1),YLin(1:end-1)); % develop mesh
AL = flipud(AL);
%% GCPs
GCPs = readtable(['GCPs_02_15_2024.txt']);    % 02/15/2024 gcps
GCPs = table2array(GCPs);
GCPeast = GCPs(:,2);
GCPnorth = GCPs(:,3);
GCPZ = GCPs(:,4);

%% Put to new grid
load('15cmgrid2.mat');
[MeanZ,Meannumpts] = roundgridfun(X,Y,A,TgridX,TgridY,@mean);

[LidarZ,Meannumpts2] = roundgridfun(XL,YL,AL,TgridX,TgridY,@mean);
ZeroInd = find(LidarZ == 0);
LidarZ(ZeroInd) = NaN;

% takes mean of all values in cell 
% meanfunctCorrect = find(MeanZ == 0);
% MeanZ(meanfunctCorrect) = NaN; % rounding grid function creates false zeros, making them nan
% MeanNAN = isnan(MeanZ); % resolving nan values from Metshape DEM
% findMeanNAN = find(MeanNAN > 0);
% AL1(findMeanNAN)= NaN; 
% MeanCOMP = AL1 - MeanZ; % Compare DEMs, take difference in elevation
%% MeanRMSE = rmse(AL1,MeanZ);

figure(40);
pcolor(X,Y,A); hold on;
plot(TgridX(190,:),TgridY(190,:),'r','Linewidth',1)
plot(TgridX(end,:),TgridY(end,:),'b','Linewidth',1)
% plot(GCPeast(1:6),GCPnorth(1:6),'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('02/15/2024',' DEM');
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)');
hold off;





%% Lidar Profile vs Meta Profile
row = 20; %which transect 
medmeta = zeros(1,367);
medlidar = zeros(1,367);
MetaComp = zeros(3,367);
MetaComp(1,:) = MeanZ(row,:);
MetaComp(2,:) = MeanZ(row+1,:);
MetaComp(3,:) = MeanZ(row+2,:);
LidarComp =  zeros(3,367);
LidarComp(1,:) = LidarZ(row,:);
LidarComp(2,:) = LidarZ(row+1,:);
LidarComp(3,:) = LidarZ(row+2,:);

for i = 1:367;
    
    if isnan(MetaComp(1,i)) && isnan(MetaComp(2,i)) && isnan(MetaComp(3,i))
        medmeta(i) = nan;
    elseif isnan(MetaComp(1,i)) && isnan(MetaComp(2,i))
        medmeta(i) = MetaComp(3,i);
    elseif isnan(MetaComp(1,i)) && isnan(MetaComp(3,i))
        medmeta(i) = MetaComp(2,i);
    elseif isnan(MetaComp(2,i)) && isnan(MetaComp(3,i))
        medmeta(i) = MetaComp(1,i);
    else
        medmeta(i) = median(MetaComp(:,i), 'omitnan');
    end
end
   
for i = 1:367;
    
    if isnan(LidarComp(1,i)) && isnan(LidarComp(2,i)) && isnan(LidarComp(3,i))
        medlidar(i) = nan;
    elseif isnan(LidarComp(1,i)) && isnan(LidarComp(2,i))
        medlidar(i) = LidarComp(3,i);
    elseif isnan(LidarComp(1,i)) && isnan(LidarComp(3,i))
        medlidar(i) = LidarComp(2,i);
    elseif isnan(LidarComp(2,i)) && isnan(LidarComp(3,i))
        medlidar(i) = LidarComp(1,i);
    else
        medlidar(i) = median(LidarComp(:,i), 'omitnan');
    end
end

 %%   row 20 on recent
Px = 0:0.15:55;
figure(11)
plot(Px,LidarZ(20,:),'r','Linewidth',2)
hold on
plot(Px,medmeta,'b','Linewidth',2)
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m), NAVD88');
title('Regrided Lidar Comparison','2/15/2024')
legend('Lidar','Stereo')
hold off





%% Transects 3/26/24
fname = '1711479601066_DEM.tif';  % 03/26/2024
[A,R] = readgeoraster([fname]); 
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A);
% % 
Tran1 = readtable('C:\Temp\crl1884\WORK\03262024Transects.xlsx');
Tran = table2array(Tran1(2:end,:));
XT = Tran(:,1);
YT = Tran(:,2);
ZT = Tran(:,3);

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

T1 = 1:80;
T2 = 92:145;
T3 = 158:208;
T4 = 220:266;
T5 = 273:325;
T6 = 347:380;
T7 = 403:440;
T8 = 456:494;
T9 = 504:545;

T1gridX = 239

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


%%
gridX = 239695:0.15:239750;
gridY = 3784715:0.15:3784770;
[Gx,Gy] = meshgrid(gridX,gridY);
% Rotate it back 
LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
XoGr = Gx - LocalX; % subtract origin
YoGr = Gy - LocalY; 
RYGr = YoGr.*cosd(-theta) + XoGr.*sind(-theta);
RXGr = XoGr.*cosd(-theta) - YoGr.*sind(-theta);
TgridX = RXGr + LocalX;
TgridY = RYGr + LocalY;





