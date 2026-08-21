function[v,points_array,ztran_matrix]=tran_video(DEMzmatrixpath,epochnum,yavg,dxy,ypick,savepath)
%function[v,points_array,ztran_matrix]=tran_video(DEMzmatrixpath,epochnum, yavg,dxy,ypick,savepath)
% 
% This function takes in a GEM and pulls a cross-shore transect of a 
% specified alongshore average, creates a video over time, and saves 
% the transect elevation values as a matrix over time. The lowest transect
% elevation values are used as the bed surface of the beach, whereas the 
% rest are depicted as incoming waves. 
% Each frame of the wave video is saved as an individual
% photo (.png) and evaluated for "fullness".
% 
% INPUTS:
% ---------
% DEMzmatrixpath = path to DEM elevation mat file matrix
% epochnum = epoch number of DEM
% yavg = width of transect in alongshore direction (m)
% dxy = matrix bin size
% ypick = transect location (m in the alongshore)
% savepath = path to save transect matrix, video, and individual timesteps
% to 
% 
% OUTPUTS: 
% ----------
% v = video of wave by wave interactions with the bed
% points_array = an array of values corresponding to each timesteps 
% number of points resolved (0 to 1, = (number of  points/number of (points+NaNs))
% ztran_matrix = a matrix of transect elevation values over time

% EDITS NEEDED
% -------------
% change to moving average of lowest points for beach face every 5 min? 
% remove points within beachface range from water values

% Last Updated: 08/17/26 BG
% -------------------------------------------------------------------------
% CONSTANTS
% -------------------------------------------------------------------------
xloc=239737;
yloc=3784751;
rotang=35;
iyavg=round(yavg/dxy/2); % indices to pull out before and after transect to average over
ztran_savepath=savepath; % save transect elevation matrix
figpath=append(savepath,'Figures/','alongshore_loc',num2str(ypick),'/'); % subfolder for figures
videopath=append(figpath,'video/'); % save video
tsteppath=append(figpath,'timesteps/'); % individual transect pngs and points_array
if ~isfolder(figpath)
    mkdir(figpath);
    mkdir(videopath);
    mkdir(tsteppath);
end
ftsz = [20 16];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 0.5;
xlab = 'Cross-shore (m)';
ylab = 'Elev. (m)';
%--------------------------------------------------------------------------
% STEP 1: Load DEM and Pull Transect at Given Location/Width
% -------------------------------------------------------------------------
% Load elevation map
DEMz=load(DEMzmatrixpath);
DEMz=DEMz.meanDEMz;
% DEM name
DEMdate=datetime(str2num(epochnum),'ConvertFrom','epochtime','TicksPerSecond',1000);
DEMdate=string(DEMdate);
DEMtitle=append(epochnum,',',DEMdate);
% create and rotate grid 
gridX=0:dxy:110;
gridY=-80:dxy:25;
[x,y]=meshgrid(gridX,gridY);
% extranct transect in all frames
[~,iy] = min(abs(y(:,1)-ypick));
ztran = median(DEMz(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran); % one matrix
ztran = movmean(ztran,2,1,'omitnan');
% compute beach as minimum (THIS CAN BE IMPROVED)
zbeach = min(ztran,[],2,'omitnan'); % beach elevation
[~,ix] = min(abs(x(1,:)-44));
zbeachmean = mean(ztran,2,'omitnan');
[~,ixmean] = min(abs(x(1,:)-12));
ix2 = find(x(1,:) > 34 & x(1,:) < 36);
[~,ixtran] = min(abs(x(1,:)-14));
[~,ixon] = min(abs(x(1,:)-2));
ztran_matrix=ztran+zbeach;
%ztran(ztran-zbeach < 0.02)  = NaN; % NaN out z elevations if they're within 2 cm of beach elevation
%--------------------------------------------------------------------------
% STEP 2: Create Video of the Transect Over Time
% -------------------------------------------------------------------------
sname =append(epochnum,'_ypick',string(ypick(1)),'_width',string(yavg));
% this is epochnum_yloc#_width#
numframes=size(ztran,2);

% create timestep video
v = VideoWriter(append(videopath,sname), 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)
figure('units','inches','position',[1 1 10 3],'color','w');
for i = 1:numframes %+1
    clf
    plot(x(1,ixon:ix)*dxy,zbeach(ixon:ix),'LineWidth',3,'Color',[148, 116, 27]/256) % plots beach (minimum transect)
    hold on
    plot(x(1,ix:end)*dxy,ztran(ix:end,i),'LineWidth',3,'Color','b') % plots water in swash
    plot(x(1,ixtran:ix)*dxy,ztran(ixtran:ix,i),'LineWidth',3,'Color','b') % plots water in uprush
    box on
    ylim([0 4.5]);
    xlim([0 23]);
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    title(['$yloc$ = ',num2str(ypick),'m',' $t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
    xlabel("Cross-Shore (m)",'Fontsize',18,'Fontname','Times New Roman'); ylabel("Elevation (m)",'Fontsize',18,'Fontname','Times New Roman')
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)
close(figure);
fprintf('Saved video of transect to %s\n', append(videopath,sname));

% -------------------------------------------------------------------------
% STEP 3: Plot Timesteps Individually & Calculate Fullness
% -------------------------------------------------------------------------
points_array=zeros(numframes,1);

fig=figure('units','inches','position',[1 1 10 3],'color','w');
for i = 1:numframes%+1
    clf(fig);
    plot(x(1,ixon:ix),zbeach(ixon:ix),'LineWidth',3,'Color',[148, 116, 27]/256) % plots beach (minimum transect)
    hold on
    plot(x(1,ix:end),ztran(ix:end,i),'LineWidth',3,'Color','b') % plots water in swash
    plot(x(1,ixtran:ix),ztran(ixtran:ix,i),'LineWidth',3,'Color','b') % plots water in uprush
    box on
    ylim([0 4.5]);
    xlim([0 100]);
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    title(['$yloc$ = ',num2str(ypick),'m',' $t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
    xlabel("Cross-Shore (m)",'Fontsize',18,'Fontname','Times New Roman'); ylabel("Elevation (m)",'Fontsize',18,'Fontname','Times New Roman')
    % save fig as png
    filename=append('tstep_',num2str(i),'.png');
    timefigpath=fullfile(tsteppath, filename);
    saveas(fig,timefigpath,'png');
    % calclate number of points resolved
    nannum=sum(isnan(ztran(:,i)));
    numpoints=length(ztran(:,i));
    points_array(i)=(numpoints-nannum)/numpoints; % (later use to separate usable transects for MAE prediction)
   
end
close(fig);
matname=fullfile(ztran_savepath,append('ztranmatrix','_ypick',string(ypick(1)),'_width',string(yavg),'.mat'));
save(matname,'ztran_matrix');
matname=fullfile(tsteppath,append('pts_resolved','_ypick',string(ypick(1)),'_width',string(yavg),'.mat'));
save(matname,'ztran_matrix');
fprintf('saved ztran matrix and points resolved');
% plot # points resolved over time
fig=figure('units','inches','position',[1 1 7 5],'color','w');clf;
plot(points_array,'o','Color','k','MarkerFaceColor','m','MarkerSize',4); hold on;
hold on; legend('Resolved data points per timestep'); xlabel('Timestep (time/2)');
ylim([0 1]);xlim([0 numframes+1]); ax=gca; ax.XTick=unique(round(ax.XTick));
ylabel('valid points/total points in transect');
ptsname=append(tsteppath,'points_resolved');
saveas(fig,ptsname,'png');
close(fig);
fprintf('saved plot of points resolved');
