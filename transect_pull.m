function[ztran]=transect_pull(GEMpath,numframes,yavg,dxy,ypick,savepath,figpath)
% function[ztran]=transect_pull(GEMpath,dxy,transect_loc)
% 
% This function pulls a transect for each GEM given and saves it as a mat
% file. It also creates a plot and saves it as a figure.
% 
% INPUTS:
% --------
% GEMpath = path to GEM mat files 
% numframes = number of frames in GEM
% yavg = width of transect in alongshore direction (m)
% dxy = matrix bin size
% ypick = transect location [x y] format
% savepath = path to save elevation values of the transect
% figpath = path to save transect plot
%
% OUTPUTS:
% ---------
% ztran = transect elevation values at given location and width

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

[~,iy]=min(abs(y(:,1)-ypick(1)));
ztran=median(GEMz(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran=squeeze(ztran);
ztran=movmean(ztran,2,1,'omitnan');

save(savepath,'ztran','-mat');

% plotting specs 
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

% plot and save for each frame
for j=1:numframes
    fig=figure('units','inches','position',[1 1 10 3],'color','w');
    hold on
    plot(x(1,:),ztran(:,j),'LineWidth',1.5);
    box on
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    %hc = colorbar('Location','eastoutside','Position', [0.93 0.225 0.03 0.4],'orientation','vertical','YAxisLocation','right');
    %set(hc,'fontsize',ftsz(2),'linewidth',lw);
    text(164,3.3,'$t$ (s)','interpreter','latex','fontsize',ftsz(1));

    sname=append(GEMname,'_tran_',string(j));
    figname=append(figfolder,'/',sname,'.png');
    exportgraphics(fig,figname);
    close(fig)
end 





