function[timeseries]=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,savepath)
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
% timeseries = .mat file of water elevation values at a single point over
% time (saved as savepath/timeseries_ptval_numframes.mat, pt_val is rounded
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

% extract one point over time
timeseries=ztran(pt,:);
pt_val=pt*dxy; 
pt_val=string(pt_val);
pt_valsave=string(floor(pt*dxy));
matname=fullfile(savepath,append('/Matfiles/timeseries_',pt_valsave,'_',string(numframes),'.mat'));
save(matname,'timeseries')

% plot time series
fig=figure('units','inches','position',[1 1 10 3],'color','w');
hold on
plot(timeseries,'LineWidth',3); hold on; ylabel("Elevation (m)");
xlabel("Time (s)"); title(append("GEM Derived Wave Time Series at ",pt_val,"m From Camera System"));
filename=append('timeseries_plot_',pt_valsave,'_',string(numframes));
figpath=append(savepath,'/Plots/',filename);
saveas(fig,figpath,'png');
close(fig);
