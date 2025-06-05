function[wave_video]=tran_video(GEMpath,yavg,dxy,ypick,figpath)
%function[wave_video]=tran_video(ztranpath)
% 
% This function takes in a GEM and creates a video showing timesteps of 1 
% second a specified transect. The lowest transect values are used as 
% the bed surface of the beach, whereas the rest are depicted as 
% incoming waves.
% 
% INPUTS:
% ---------
% GEMpath = path to GEM mat files 
% yavg = width of transect in alongshore direction (m)
% dxy = matrix bin size
% ypick = transect location [x y] format
% figpath = path to save video to
% 
% OUTPUTS: 
% ----------
% wave_video = video of wave by wave interactions with the bed

% Constants
xloc=239737;
yloc=3784751;
rotang=35;
iyavg=round(yavg/dxy/2); % indices to pull out before and after transect to average over
figfolder=figpath;

% Load GEM
GEMz=load(GEMpath);
GEMz=GEMz.meanGEMz;

% GEM name
GEMname=split(GEMpath,'/');
GEMname=GEMname(9,1);
GEMname=string(GEMname);
GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMdate=string(GEMdate);
GEMtitle=append(GEMname,',',GEMdate);


% create and rotate grid 
gridX=0:dxy:110;
gridY=-80:dxy:25;
[x,y]=meshgrid(gridX,gridY);

% extranct transect 1
[~,iy] = min(abs(y(:,1)-ypick(1)));
ztran = median(GEMz(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran);
ztran = movmean(ztran,2,1,'omitnan');

% compute beach as minimum (THIS CAN BE IMPROVED)
zbeach = min(ztran,[],2,'omitnan'); % beach elevation
[~,ix] = min(abs(x(1,:)-44));
zbeachmean = mean(ztran,2,'omitnan');
[~,ixmean] = min(abs(x(1,:)-12));
ix2 = find(x(1,:) > 34 & x(1,:) < 36);

[~,ixtran] = min(abs(x(1,:)-14));
[~,ixon] = min(abs(x(1,:)-2));

ztran(ztran-zbeach < 0.04)  = NaN; % NaN out z elevations if they're within 4 cm of beach elevation

% remove values that are less than 5 indeces long (removes features that
% are intermittant -- likely remove this in future. This is just to improve
% the video visualization
%numstring = 5;
%ztran_filt = ztran;
%for i = 1:size(ztran,2)
    %ztran_filt(:,i) = replaceShortNonNanSequences(squeeze(ztran(:,i)), numstring); % this is water elevation
%end

% constants
numframes=size(ztran,2);
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

sname =append(GEMname,'_timeseries');

% need to test this
v = VideoWriter(append(figfolder,'\',sname), 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)

figure('units','inches','position',[1 1 10 3],'color','w');

for i = 1:numframes+1
    clf
    plot(x(1,ixon:ix),zbeach(ixon:ix),'LineWidth',3,'Color',[148, 116, 27]/256)
    hold on
    plot(x(1,ix:end),ztran(ix:end,i),'LineWidth',3,'Color','b')
    plot(x(1,ixtran:ix),ztran_filt(ixtran:ix,i),'LineWidth',3,'Color','b')

    box on
    ylim([0 4.5]);
    xlim([0 100]);
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    title(['$t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)
