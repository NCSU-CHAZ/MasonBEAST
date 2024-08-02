% for i = 1:6
% filename =  '1711479601066v2_DEM.tif';
%     filename = {'1699459201959_DEM.tif','1702483201914_DEM.tif','1704225601966_DEM.tif','1706904001127_DEM.tif','1709409601439_DEM.tif','1712242801208_DEM.tif'};
%     filename = '1699459201959_DEM.tif';
%     filename = '1702483201914_DEM.tif';
%     filename = '1704225601966_DEM.tif';
% filename = '1706904001127_DEM.tif';
% filename = '1709409601439_DEM.tif';
%     filename = '1712242801208_DEM.tif';
% filename = '1702670401385_DEM.tif' ; 
% filename = '1703088001887_DEM.tif' ;
% filename = '1696014001553_DEM.tif' ; 
% filename = '1698159601837_DEM.tif' ; % 10/24/2023 11am
% filename = '1708099201978v2_DEM.tif';  % 02/15/2024
% filename = '1697482801712_DEM.tif' ; % 10/16/2023
% filename = '1713898801714_DEM.tif';  % 04/23/2023
filename = '1715108401562_DEM.tif';  % 05/7/2024

% figname = {'11/08/2023','12/13/2023','1/2/2024','2/2/2024','3/2/2024','4/4/2024'};
    [A,R] = readgeoraster([filename]);
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A);

%%
 % GCPs = readtable([GCPpath,'20231116.txt']);  % 11/16/2023 gcps 
GCPs = readtable(['GCPs_02_15_2024.txt']);    % 02/15/2024 gcps
GCPs = table2array(GCPs);
GCPeast = GCPs(:,2);
GCPnorth = GCPs(:,3);
GCPZ = GCPs(:,4);
%% Rotate

LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
Xo = X - LocalX; % subtract origin
Yo = Y - LocalY; 
theta = 35;
RY = Yo.*cosd(theta) + Xo.*sind(theta);
RX = Xo.*cosd(theta) - Yo.*sind(theta);
Xout = RX + LocalX;
Yout = RY + LocalY;
%% Rotate GCPS
LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
XoG = GCPeast - LocalX; % subtract origin
YoG = GCPnorth - LocalY; 
RYG = YoG.*cosd(theta) + XoG.*sind(theta);
RXG = XoG.*cosd(theta) - YoG.*sind(theta);
GCPeastout = RXG + LocalX;
GCPnorthout = RYG + LocalY;
% cams 
CamX= [239766.5, 239759.38];
CamY = [3784761.5,3784754.9];
XoC = CamX - LocalX; % subtract origin
YoC = CamY - LocalY; 
RYC = YoC.*cosd(theta) + XoC.*sind(theta);
RXC = XoC.*cosd(theta) - YoC.*sind(theta);
CamXout = RXC + LocalX;
CamYout = RYC + LocalY;
%% Rotate Lidar 
theta = 35;
XoL = XL - LocalX; % subtract origin
YoL = YL - LocalY; 
RYL = YoL.*cosd(theta) + XoL.*sind(theta);
RXL = XoL.*cosd(theta) - YoL.*sind(theta);
XLout = RXL + LocalX;
YLout = RYL + LocalY;


%%
CompX(1,:) = [239700:0.15:239750];
CompX(2,:) = [239700:0.15:239750];
CompY(1,:)= zeros(1,length(CompX));
CompY(2,:)= zeros(1,length(CompX));
% CompY(1,:) = CompY(1,:) + 3784745.15;
% CompY(2,:) = CompY(2,:) + 3784745;
CompY(1,:) = CompY(1,:) + 3784750.15;
CompY(2,:) = CompY(2,:) + 3784750;

%% Develop grid 
% gridX = 239685:0.15:239780;
% gridY = 3784660:0.15:3784760;
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

%% Display DEM
figure(30);
pcolor(Xout,Yout,A); hold on;
plot(GCPeastout(1:6),GCPnorthout(1:6),'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('November 14th, 2023',' DEM');
ylim([3784650,3784770]);
xlim([239680,239800]);
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)');
% legend('','Transect Taken')
plot(CamXout,CamYout,'bs')
legend('','Ground Control Points','Cameras')
% plot(TgridX(235,:),TgridY(235,:),'r','Linewidth',1)
% plot(TgridX(end,:),TgridY(end,:),'r','Linewidth',0.2)
hold off;


%%
N = scatteredInterpolant(Xout(:),Yout(:),double(A(:)),'nearest');
IntN = N(CompX,CompY);
Z = zeros(1,length(CompX));
for i = 1:334;
%     disregard(1,:) = isnan(IntN(1,:));
%     disregard(2,:) = isnan(IntN(2,:));
    if isnan(IntN(1,i)) && isnan(IntN(2,i))
        Z(i) = nan;
    elseif isnan(IntN(1,i))
        Z(i) = IntN(2,i);
    elseif isnan(IntN(2,i))
        Z(i) = IntN(1,i);
    else
        Z(i) = mean(IntN(:,i), 'omitnan');
    end
end
%% Lidar 
N = scatteredInterpolant(XLout(:),YLout(:),double(AL(:)),'nearest');
IntN = N(CompX,CompY);
Z = zeros(1,length(CompX));
for i = 1:334;
%     disregard(1,:) = isnan(IntN(1,:));
%     disregard(2,:) = isnan(IntN(2,:));
    if isnan(IntN(1,i)) && isnan(IntN(2,i))
        Z(i) = nan;
    elseif isnan(IntN(1,i))
        Z(i) = IntN(2,i);
    elseif isnan(IntN(2,i))
        Z(i) = IntN(1,i);
    else
        Z(i) = mean(IntN(:,i), 'omitnan');
    end
end
%% 
figx = 0:0.15:50;
figure(8);
plot(figx,Z,'r','Linewidth',3)
hold on;
plot(figx,ZMeta0215,'b','Linewidth',3)
pbaspect([3 2 1])
% title('02/15/2024','Lidar Comparison')
xlabel('cross-shore distance (m)')
ylabel('elevation (m NAVD88)');
legend('Lidar','Stereo');
xlim([2,30])
ylim([0.5,4])
hold off;





%% Display ALL 
load ZApr.mat
load ZDec.mat
load ZFeb.mat
load Zjan.mat
load Zmar.mat
load ZNov.mat
load ZSept.mat
load ZOct.mat
load ZDec2.mat
load ZDec3.mat
load Colors.mat


figx = 0:0.15:50;
%%
figure(75);
plot(figx,ZSept,'ko'); 
% hold on;
% plot(figx,ZNov,'c');
% plot(figx,ZDec,'r');
% plot(figx,ZJan,'k');
% plot(figx,ZFeb,'g');
% plot(figx,ZMar,'m');
% plot(figx,ZApr,'b');
% plot(figx,ZDec2,'k','Linewidth',3);
% plot(figx,ZDec3,'k-','Linewidth',3);
legend('September','November','December','January','February','March','April');
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m), NAVD88');
title('Transects over time');
hold off

%% Define shades of red
red_shades = linspace(1, 0, 4); % Adjust the number to change the intensity

% Up to December
% figure(2);
% plot(figx,ZSept,'Color', [1, red_shades(1), red_shades(1)]);
% hold on
% plot(figx,ZOct,'Color', [1, red_shades(2), red_shades(2)]);
% plot(figx,ZNov,'Color', [1, red_shades(3), red_shades(3)]);
% plot(figx,ZDec,'Color', [1, red_shades(4), red_shades(4)]);
figure(3);
plot(figx,ZSept,'r','Linewidth',2)
hold on
plot(figx,ZOct,'b','Linewidth',2)
plot(figx,ZNov,'g','Linewidth',2)
plot(figx,ZDec,'m','Linewidth',2)
legend('9/28/2023','10/16/2023','11/8/2023','12/13/2023');
title('September - December Transects');
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m), NAVD88');
xlim([0,27.5])
hold off
%%
TranZ = zeros(4,334);
TranZ(1,:) = ZSept;
TranZ(2,:) = ZOct;
TranZ(3,:) = ZNov;
TranZ(4,:) = ZDec;
n = 4; % number of plots
cmap = jet(6);
figure('Color','w');
hold on
for i = 1:4
    plot(figx,TranZ(i,:),'Color',cmap(i,:),'LineWidth',2)
end
colormap(cmap)
caxis([1 4])
c = colorbar;
c.Label.String = 'Month';
ticks = linspace(1, 4, 4); % Assuming 1 corresponds to September, 4 corresponds to December
ticklabels = {'Sept', 'Oct', 'Nov', 'Dec'};
c.Ticks = ticks;
c.TickLabels = ticklabels;
title('September - December Transects');
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m), NAVD88');
xlim([0,27.5])
hold off;
%% Same fig with Pre and Post Storm
n = 6; % number of plots
TranZ(5,:) = ZDec2;
TranZ(6,:) = ZDec3;
cmap = jet(6);
figure('Color','w');
hold on
for i = 1:6;
    plot(figx,TranZ(i,:),'Color',cmap(i,:),'LineWidth',2)
end
colormap(cmap)
caxis([1 6])
c = colorbar;
c.Label.String = 'Month';
ticks = linspace(1, 6, 6); % Assuming 1 corresponds to September, 4 corresponds to December
ticklabels = {'Sept', 'Oct', 'Nov', 'Dec','Pre-Storm','Post-Storm'};
c.Ticks = ticks;
c.TickLabels = ticklabels;
title('12/17 Pre and Post Storm Comparison');
xlabel('Cross Shore Distance (m)')
ylabel('Elevation (m), NAVD88');
xlim([0,27.5])
hold off;
%% Same plots but one by one 
figure(20);
% plot(figx,ZSept,'Color',d,'Linewidth',2);
% hold on
% plot(figx,ZOct,'Color',b,'Linewidth',2)
% plot(figx,ZNov,'Color',c,'Linewidth',2)
% plot(figx,ZDec,'Color',g,'Linewidth',2)
plot(figx,ZDec2,'k','Linewidth',2)
hold on
plot(figx,ZDec3,'Color',f,'Linewidth',2)
pbaspect([3 2 1])
plot(figx,ZJan,'Color',h,'Linewidth',2)
plot(figx,ZFeb,'Color',i,'Linewidth',2)
plot(figx,ZMar,'Color',j,'Linewidth',2)
plot(figx,ZApr,'Color',k,'Linewidth',2)
% legend('9/28/2023','10/16/2023','11/8/2023','12/13/2023','Pre-Storm','Post-Storm','1/2/2024','2/02/2024','3/02/2024','4/04/2024','location','southwest');
legend('Pre-Storm','Post-Storm','1/02/2024','2/02/2024','3/02/2024','4/04/2024','location','southwest');
% legend('Pre-Storm','Post-Storm')
title('2024 Winter - Spring');
xlabel('cross-shore distance (m)')
ylabel('elevation (m NAVD88)');
xlim([2,30])
ylim([0.5 4])
hold off
