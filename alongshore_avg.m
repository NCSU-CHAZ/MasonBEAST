% Alongshore averaging variation workflow
% B Gaenzle
% 11/4/25 Updated

%___________________________________________________________________________
% This code takes in pointclouds and creates a wave transect video
% for a chosen transect of varying alongshore widths (0.5m, 1.2m,1.5m,2m).
% This gappiness is quantified through quality of transects. 
%___________________________________________________________________________
% Running on Dec 17th, 2023 Nor'Easter 3:30 PM (epoch = 1702827001820)

% Paths
% ---------------
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
stormCHAZerspath=append(genpath,'/StormCHAZerz Data');
decNoreastherpath=append(stormCHAZerspath,'/Dec2023Noreaster_Processed');

epochnum='1702827001820'; % epoch number 
ptcldpath=append(decNoreastherpath,'/ptclds/',epochnum,'_allframes');
figpath=append(ptcldpath,'/Figures');
savepath=append(ptcldpath,'/MAPs/');
% GEM specs
% --------------
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.2; % meter
% Process Pointclouds
% ------------------------
spec='*_ptcld*';
[meanMAPz_matrix,medMAPz_matrix,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,ptcldpath,savepath,figpath,spec);

% Create Transect Video 
% ------------------------
ypick = [-30 -20]; % in m, define locations to pick transects
MAPzmatrixpath=append(savepath,'/meanMAPz.mat');
% 0.5 m alongshore avg
yavg = 0.5;
newfigpath=append(figpath,'/Transect_yneg30','/05yavg');
[v_05,quality_array05]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

yavg = 1.2; % in m, defining how wide in the alongshore to average over
newfigpath=append(figpath,'/Transect_yneg30','/12yavg');
[v_12,quality_array12]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

yavg=1.5;
newfigpath=append(figpath,'/Transect_yneg30','/15yavg');
[v_15,quality_array15]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

yavg=2.0;
newfigpath=append(figpath,'/Transect_yneg30','/20yavg');
[v_20,quality_array20]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

yavg=2.5;
newfigpath=append(figpath,'/Transect_yneg30','/25yavg');
[v_25,quality_array25]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

yavg=3.0;
newfigpath=append(figpath,'/Transect_yneg30','/30yavg');
[v_30,quality_array30]=tran_video(MAPzmatrixpath,yavg,dxy,ypick,newfigpath);

% plot all quality on 1 plot
fig=figure('units','inches','position',[1 1 10 7],'color','w');clf;
plot(quality_array05,'o','Color','r','MarkerFaceColor','r','MarkerSize',6); hold on;
plot(quality_array12,'o','Color','m','MarkerFaceColor','m','MarkerSize',6); hold on;
plot(quality_array15,'o','Color','y','MarkerFaceColor','y','MarkerSize',6); hold on;
plot(quality_array20,'o','Color','g','MarkerFaceColor','g','MarkerSize',6); hold on;
plot(quality_array25,'o','Color','c','MarkerFaceColor','c','MarkerSize',6); hold on;
plot(quality_array30,'o','Color','b','MarkerFaceColor','b','MarkerSize',6); hold on;
hold on; ylabel('Quality Value (# NaNs/# points in transect)');
legend('yavg = 0.5m','yavg = 1.2m','yavg = 1.5m','yavg = 2.0m','yavg = 2.5m','yavg = 3.0m'); xlabel('Image # (time/2)');
ylim([0 1]);xlim([0 numframes+1]); ax=gca; ax.XTick=unique(round(ax.XTick));title('Quality Value for Transects');
% Plot Time Series of Various Transects at a Point
% -------------------------------------------------
TS_path='/Volumes/kanarde/MasonBEAST/data/StormCHAZerz Data/Dec2023Noreaster_Processed/1702827001820_allframes/WaveTimeSeries';
pt=73;
timeseries73=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
pt=76;
timeseries76=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
pt=170;
timeseries170=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
% plot GEM and points on it (NEED TO ADD TO FUNCTION)
