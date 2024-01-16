% Genereating DEM figures from Metashape
generalpath = 'C:/MASONbeast/';
datapath = [generalpath,'data/MetashapeTIFs/DEMs/']; % setting paths import/export
figfolder = [generalpath,'data/Figures/DEMfigs/'];
fname = '1699992001273v7_DEM.tif';
figname = '1699992001273v7_DEMplot';
[A,R] = readgeoraster([datapath,fname]); % read in tif file
A(A<-7) = NaN; % Create reasonable z limit
Xin = R.XWorldLimits(1):R.CellExtentInWorldX:R.XWorldLimits(2);
Yin = R.YWorldLimits(1):R.CellExtentInWorldY:R.YWorldLimits(2);
[X,Y] = meshgrid(Xin(1:end-1),Yin(1:end-1)); % develop mesh

LocalX =  239750; % Local Origin - will vary before everything is automated
LocalY = 3784295;
% Rotation 
Xo = X - LocalX; % subtract origin
Yo = Y - LocalY; 
theta = -30; % angle of coastline(roughly)
RY = Yo.*cosd(theta) + Xo.*sind(theta);
RX = Xo.*cosd(theta) - Yo.*sind(theta);
Xout = RX + LocalX;
Yout = RY + LocalY;

figure(1); % Generate figure
pcolor(RX,RY,A); hold on;
shading flat;
axis image;
xlim([0,80]);
ylim([0, 80]);
demcmap([-2 5]);
h = colorbar; h.YLabel.String = 'Elevation [m]';
xlabel('Cross-shore (m)');
ylabel('Along-shore (m)');
title(figname);
hold off;
figure(1);
Sname = [figfolder,figname]; % saving figure
print(Sname,'-dpng');
