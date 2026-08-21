function[DEMvid]=DEMvideo(DEMzmatrixpath,epochnum,camlocA,camlocB,rbrloc,dxy,savepath)
%function[DEMvid]=DEMvideo(DEMzmatrixpath,camlocA,camlocB,rbrloc,dxy,savepath)
% -------------------------------------------------------------------------
% This function generates a video of DEM timeseries and saves it as a MP4
% file.
%
% INPUTS:
% ------------
% DEMzmatrixpath = path to DEM elevation mat file matrix
% epochnum = epoch number of DEM
% camloc = camera A location coordinates [Ax,Ay]
% camloc = camera B location coordinates [Bx,By]
% rbrloc = rbr location coordinates [x,y]
% dxy = matrix bin size
% savepath = path to save the video to
%
% OUTPUTS:
% -------------
% DEMvid = video of DEM over time
% 
% Last Updated: 08/17/26 BG
% ------------------------------------------------------------------------
% Load elevation map
DEMz=load(DEMzmatrixpath);
DEMz=DEMz.meanDEMz;
% DEM name
DEMdate=datetime(str2num(epochnum),'ConvertFrom','epochtime','TicksPerSecond',1000);
DEMdate=string(DEMdate);
DEMtitle=append(epochnum,',',DEMdate);
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
% rbr rotation and transformation
[rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);

sname =append(epochnum,'_DEM');
% this is epochnum_yloc#_width#
numframes=size(DEMz,3);

% create timestep video
v = VideoWriter(append(savepath,sname), 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)

fig=figure('units','inches','position',[0 0 10 6],'color','w');
xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';

for i=1:numframes
    clf(fig);
    Z=DEMz(:,:,i);
    h=pcolor(Xgrid,Ygrid,Z); grid off;box on;hold on;
    set(h,'EdgeColor','none');
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    scatter(rbrx,rbry,60,'fill','sq','k','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'k');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw); hc.Label.String = 'Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);
    title(DEMtitle);
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)



