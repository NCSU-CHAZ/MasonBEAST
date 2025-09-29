function[wavetimeseries]=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,savepath)
% function[timeseries,plot]=GEMmatrix_to_timeseries(GEMmatrixpath,pt,transect)
% -------------------------------------------------------------------------
% This function takes in a matrix of GEM elevation values (in sequential
% order) and pulls out a time series of water surface elevation at a single
% point. 
% --------------------------------------------------------------------------
% INPUTS:
% -------
% GEMmatrixpath = path to GEM matrix
% pt = point where to take time series
% dxy = resolution of GEM
% ypick = transect locations
% yavg = alongshore avg of transect
% savepath = path to save timeseries .mat file and plots
% 
% OUTPUTS:
% --------
% wavetimeseries = .mat file of water elevation values at a single point over
% time (saved as savepath/timeseries_ypick(1)_ptval_numframes.mat, pt_val is rounded
% down for no decimals)
%

% Load GEMz matrix
GEMz=load(GEMmatrixpath);
GEMz=GEMz.meanGEMz;
numframes=size(GEMz,3);

% GEM name
GEMname=split(GEMmatrixpath,'/');
GEMname=GEMname(8,1);
GEMname=string(GEMname);
GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMdate=string(GEMdate);
GEMtitle=append(GEMname,',',GEMdate);

% create and rotate grid 
gridX=0:dxy:110;
gridY=-80:dxy:25;
[x,y]=meshgrid(gridX,gridY);
iyavg=round(yavg/dxy/2); % indices to pull out before and after transect to average over

% extranct transect in all frames
[~,iy] = min(abs(y(:,1)-ypick(1)));
ztran = median(GEMz(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran = squeeze(ztran); % one matrix
ztran = movmean(ztran,2,1,'omitnan');
zbeach = min(ztran,[],2,'omitnan'); % beach elevation

% extract one point over time
timeseries=ztran(pt,:);
zbeachts=ones(numframes);
zbeachts=zbeachts(:,1);
zbeachts=zbeachts*zbeach(pt,:);

wavetimeseries=timeseries-zbeachts;
wavetimeseries=wavetimeseries(1,:);

pt_val=pt*dxy; 
pt_val=string(pt_val);
pt_valsave=string(floor(pt*dxy));
matname=fullfile(savepath,append('/Matfiles/timeseries_',string(ypick(1)),'_',pt_valsave,'_',string(numframes),'.mat'));
save(matname,'wavetimeseries')

time=0.5:0.5:(numframes/2);

% plot water elevation time series (subtract beach)
fig=figure('units','inches','position',[1 1 10 3],'color','w');
hold on
plot(time,wavetimeseries,'LineWidth',3); hold on; ylabel("Elevation (m)");
xlabel("Time (s)"); title(append("GEM Derived Wave Time Series at ",pt_val,"m From Camera System"));
ylim([0 0.5]);
filename=append('timeseries_plot_',string(ypick(1)),'_',pt_valsave,'_',string(numframes));
figpath=append(savepath,'/Plots/',filename);
saveas(fig,figpath,'png');
close(fig);

% Plot location of point with GEM
xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
fig=figure('units','inches','position',[0 0 10 6],'color','w');
pcolor(x,y,GEMz(:,:,1)); grid off;box on;hold on
shading interp;
axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
set(hc,'fontsize',ftsz(2),'linewidth',lw);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
title(GEMtitle);
hold on; yline(ypick(1),'k--','LineWidth',3);hold on,
scatter(pt*dxy,ypick(1),140,'fill','sq','m','MarkerEdgeColor','k');
filename=append('GEM_pt_plot_',pt_valsave,'_',string(numframes));
figpath=append(savepath,'/Plots/',filename);
saveas(fig,figpath,'png');
close(fig);

