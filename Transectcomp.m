%% Transects 3/26/24
% fname = '1711479601066_DEM.tif';  % 03/26/2024
fname = '1717167601959_DEM.tif';  % 05/31/2024
[A,R] = readgeoraster([fname]); 
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh
A = flipud(A);
%% 
% Tran1 = readtable('C:\Temp\crl1884\WORK\03262024Transects.xlsx');
Tran1 = readtable('C:\Temp\crl1884\WORK\05312024_trans_UTM.xlsx');
Tran = table2array(Tran1(2:end,:));
XT = Tran(:,3);
YT = Tran(:,2);
ZT = Tran(:,4);
%% Rotate X and Y
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
%% 3/26/24
T1 = 1:80;
T2 = 92:145;
T3 = 158:208;
T4 = 220:266;
T5 = 273:325;
T6 = 347:380;
T7 = 403:440;
T8 = 456:494;
T9 = 504:545;
%% 5/31/24 and 3/26/24 Transects are together so extracting 5/31
% T1 : most north 
XTout = XTout(550:end);
YTout = YTout(550:end);
ZT = ZT(550:end);
%% 5/31
T1 = 1:135;
T2 = 145:270;
T3 = 280:390;
T4 = 400:505;
T5 = 525:615;
T6 = 630:700;
T7 = 725:782;
T8 = 800:863;








%%
figure(33); % Display rotated DEM 
pcolor(Xout,Yout,A); hold on;
plot(XTout,YTout,'r')
plot(XTout(780),YTout(780),'go')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation NAVD88 [m]';
title('3/26/2024','DEM');
xlim auto;
% ylim([3784650,3784710]);
% xlim([239550,239610])
hold off;


%%
% gridX = 239695:0.15:239750; % 15 cm res
% gridY = 3784715:0.15:3784770;

gridX = 239695:1:239750; % 1 m res
gridY = 3784715:1:3784770;

[Gx,Gy] = meshgrid(gridX,gridY);
theta = 35;
%% Rotate it back 
LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
XoGr = Gx - LocalX; % subtract origin
YoGr = Gy - LocalY; 
RYGr = YoGr.*cosd(-theta) + XoGr.*sind(-theta);
RXGr = XoGr.*cosd(-theta) - YoGr.*sind(-theta);
TgridX = RXGr + LocalX;
TgridY = RYGr + LocalY;


%% 
figure(1);
pcolor(Xout,Yout,A)
shading flat
axis image; hold on;
plot(Gx,Gy,'r')
% plot(TgridX(225,:),TgridY(225,:),'r-');
% plot(XT(1:80),YT(1:80),'g')
hold off
%%
figure(2);
pcolor(X,Y,A); hold on;
% plot(GCPeast(1:6),GCPnorth(1:6),'ro')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('3/26/2024',' DEM');
% ylim([3784700,3784760]);
% xlim([239750,239820]);
xlabel('UTM Easting (m)')
ylabel('UTM Northing (m)');
plot(XT(T1),YT(T1),'r')
% plot(TgridX(218,:),TgridY(218,:),'k')
% plot(TgridX(228,:),TgridY(228,:),'k')

plot(XT(T2),YT(T2),'r')
plot(TgridX(195,:),TgridY(195,:),'k')
plot(TgridX(205,:),TgridY(205,:),'k')

plot(XT(T3),YT(T3),'r')
plot(TgridX(160,:),TgridY(160,:),'k')
plot(TgridX(170,:),TgridY(170,:),'k')

plot(XT(T5),YT(T5),'r')
plot(TgridX(95,:),TgridY(95,:),'k')
plot(TgridX(105,:),TgridY(105,:),'k')

hold off

%% Regrid Stereo to comp grid
[MeanZ,Meannumpts] = roundgridfun(Xout,Yout,A,Gx,Gy,@mean);% takes mean of all values in cell 
% meanfunctCorrect = find(MeanZ == 0);
% MeanZ(meanfunctCorrect) = NaN; % rounding grid function creates false zeros, making them nan
% MeanNAN = isnan(MeanZ); % resolving nan values from Metshape DEM
% findMeanNAN = find(MeanNAN > 0);
%% regrid transects 
[TZ,Tnumpts] = roundgridfun(XTout,YTout,ZT,Gx,Gy,@mean);% takes mean of all values in cell 



%%
[Vq,N] = gridbin(x,y,z,Gx,Gy) ;
%%

figure(3);
pcolor(Gx,Gy,MeanZ);
hold on
% plot(TgridX,TgridY,'r')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('3/26/2024',' DEM');
% ylim([3784700,3784760]);
% xlim([239750,239820]);
% pcolor(Gx,Gy,TZ)


plot(XTout(T1),YTout(T1),'r')  % T1
plot(Gx(35,:),Gy(35,:),'k')
% plot(Gx(228,:),Gy(228,:),'k')
%%
plot(Gx(34,:),Gy(34,:),'k')
% plot(Gx(35,:),Gy(35,:),'k')



plot(XTout(T2),YTout(T2),'r') % T2
% plot(Gx(30,:),Gy(30,:),'k')
plot(Gx(31,:),Gy(31,:),'k')


% % plot(XT(T3),YT(T3),'r')
% % plot(TgridX(160,:),TgridY(160,:),'k')
% % plot(TgridX(170,:),TgridY(170,:),'k')
plot(XTout(T3),YTout(T3),'r')  %T3
% plot(Gx(25,:),Gy(25,:),'k')
plot(Gx(26,:),Gy(26,:),'k')

% plot(Gx(160,:),Gy(160,:),'k') % old T3
% plot(Gx(170,:),Gy(170,:),'k')
% 
% % plot(XT(T5),YT(T5),'r')
% % plot(TgridX(95,:),TgridY(95,:),'k')
% % plot(TgridX(105,:),TgridY(105,:),'k')
plot(XTout(T5),YTout(T5),'r')
% plot(Gx(15,:),Gy(15,:),'k')
plot(Gx(16,:),Gy(16,:),'k')

% plot(Gx(95,:),Gy(95,:),'k')
% plot(Gx(105,:),Gy(105,:),'k')
hold off;

%% now make transect from taking median of regrided data
% Transect1 = zeros(1,56);
% for i = 1:56;
% %     if isnan(TZ(34,i))
% %         Transect1(i) = TZ(35,i);
% %     elseif isnan(TZ(35,i))
% %         Transect1(i) = TZ(34,i);
% %     elseif isnan(TZ(34,i)) && isnan(TZ(35,i))
% %          Transect1(i) = nan;
% %     end
%     
% Transect1(i) = nanmean(TZ(34,i),TZ(35,i));
% %     Transect2(i) = mean(TZ(30,i),TZ(31,i),'omitnan');
% %     Transect3(i) = mean(TZ(25,i),TZ(26,i),'omitnan');
% %     Transect4(i) = mean(TZ(15,i),TZ(16,i),'omitnan');
% end 
TZ(TZ == 0) = nan;


figx = 1:1:56;
figure(6); % T1
plot(figx,TZ(34,:),'b');
hold on
plot(figx,MeanZ(34,:),'r');
legend('Transect 1', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T1 Comp')
hold off

figure(7); % T2
plot(figx,TZ(31,:),'b');
hold on
plot(figx,MeanZ(31,:),'r');
legend('Transect 2', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T2 Comp')
hold off

figure(8); % T3
plot(figx,TZ(26,:),'b');
hold on
plot(figx,MeanZ(26,:),'r');
legend('Transect 3', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T3 Comp')
hold off

figure(9); % T4
plot(figx,TZ(16,:),'b');
hold on
plot(figx,MeanZ(16,:),'r');
legend('Transect 4', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T4 Comp')
hold off


%% Trying all Transects
gridX = 239700:1:239765; % 1 m res
gridY = 3784685:1:3784750;

[Gx,Gy] = meshgrid(gridX,gridY);
    
%% Regrid Stereo to comp grid
[MeanZ,Meannumpts] = roundgridfun(Xout,Yout,A,Gx,Gy,@mean);% takes mean of all values in cell 
% meanfunctCorrect = find(MeanZ == 0);
% MeanZ(meanfunctCorrect) = NaN; % rounding grid function creates false zeros, making them nan
% MeanNAN = isnan(MeanZ); % resolving nan values from Metshape DEM
% findMeanNAN = find(MeanNAN > 0);
%% regrid transects 
[TZ,Tnumpts] = roundgridfun(XTout,YTout,ZT,Gx,Gy,@mean);% takes mean of all values in cell 
%%
figure(10);
%pcolor(Gx,Gy,MeanZ);
hold on
% plot(TgridX,TgridY,'r')
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
title('3/26/2024',' DEM');
% ylim([3784700,3784760]);
% xlim([239750,239820]);
pcolor(Gx,Gy,TZ)
%plot(XTout,YTout,'r')

shading flat;
axis image;

plot(Gx(64,:),Gy(64,:),'k') % most north 
plot(Gx(61,:),Gy(61,:),'k')
plot(Gx(56,:),Gy(56,:),'k')
plot(Gx(50,:),Gy(50,:),'k')
plot(Gx(40,:),Gy(40,:),'k')
plot(Gx(36,:),Gy(36,:),'k')
plot(Gx(28,:),Gy(28,:),'k')
plot(Gx(13,:),Gy(13,:),'k')
plot(Gx(7,:),Gy(7,:),'k') % most south 
hold off
%%
% plot all transects vs Stereo 3/26
TZ(TZ == 0) = nan;
figx = 1:66;

figure(11); % T1
plot(figx,TZ(64,:),'b');
hold on
plot(figx,MeanZ(64,:),'r');
legend('Transect 1 ', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T1 Comp 5/31')
hold off
%%
figure(12); % T2
plot(figx,TZ(61,:),'b');
hold on
plot(figx,MeanZ(61,:),'r');
legend('Transect 2', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T2 Comp')
hold off
%%
figure(13); % T3
plot(figx,TZ(56,:),'b');
hold on
plot(figx,MeanZ(56,:),'r');
legend('Transect 3', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T3 Comp')
hold off

figure(14); % T4
plot(figx,TZ(50,:),'b');
hold on
plot(figx,MeanZ(50,:),'r');
legend('Transect 4', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T4 Comp')
hold off

figure(15); % T5
plot(figx,TZ(40,:),'b');
hold on
plot(figx,MeanZ(40,:),'r');
legend('Transect 5', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T5 Comp')
hold off

figure(16); % T6
plot(figx,TZ(36,:),'b');
hold on
plot(figx,MeanZ(36,:),'r');
legend('Transect 6', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T6 Comp')
hold off

figure(17); % T7
plot(figx,TZ(28,:),'b');
hold on
plot(figx,MeanZ(28,:),'r');
legend('Transect 7', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T7 Comp')
hold off

figure(18); % T8
plot(figx,TZ(13,:),'b');
hold on
plot(figx,MeanZ(13,:),'r');
legend('Transect 8', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T8 Comp')
hold off

figure(19); % T9
plot(figx,TZ(7,:),'b');
hold on
plot(figx,MeanZ(7,:),'r');
legend('Transect 9', 'Stereo');
xlabel('cross-shore distance (m)');
ylabel('elevation (m,NAVD88)');
title('T9 Comp')
hold off
%%
R1 = rmse(TZ(64,:),MeanZ(64,:));
R2 = rmse(TZ(61,:),MeanZ(61,:));
R3 = rmse(TZ(56,:),MeanZ(56,:));
R4 = rmse(TZ(50,:),MeanZ(50,:));
R5 = rmse(TZ(40,:),MeanZ(40,:));
R6 = rmse(TZ(36,:),MeanZ(36,:));
R7 = rmse(TZ(28,:),MeanZ(28,:));
R8 = rmse(TZ(13,:),MeanZ(13,:));
R9 = rmse(TZ(7,:),MeanZ(7,:));


%% 05/31 work 
% MetaZ = zeros(66,66);
% find = TZ ~= 0;
% Meta = (Gx,Gy,MeanZ);
% Transects = (Gx,Gy,TZ);

extractedMeta = MeanZ(find);
extractedTransects = TZ(find);
fx = 1:488;
figure(66);
plot(fx,extractedMeta,'r');
hold on
plot(fx,extractedTransects,'b');
hold off
rmse(extractedTransects,extractedMeta)














%% look at rotation
LocalX = 239740;  % 11/16/2023
LocalY = 3784650;
Xo = X - LocalX; % subtract origin
Yo = Y - LocalY; 
theta = 35;
RY = Yo.*cosd(theta) + Xo.*sind(theta);
RX = Xo.*cosd(theta) - Yo.*sind(theta);
Xout = RX + LocalX;
Yout = RY + LocalY;

figure(5);
pcolor(Gx,Gy,Vq);
hold on;
shading flat;
axis image;
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m] NAVD88';
% plot(XTout,YTout,'b')
hold off


