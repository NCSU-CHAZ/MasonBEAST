% Code to process pointclouds and create GEMs and wave video
% BG 08/04/2025

% Running on Dec 17th, 2023 Nor'Easter 3:30 PM (epoch = 1702827001820)

% Paths
% ---------------
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
stormCHAZerspath=append(genpath,'/StormCHAZerz Data');
decNoreastherpath=append(stormCHAZerspath,'/Dec2023Noreaster_Processed');

epochnum='1702827001820'; % epoch number 

ptcldpath=append(decNoreastherpath,'/',epochnum);
% GEM specs
% --------------
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.2; % meter
numframes=120;

% Process Pointclouds
% ------------------------
spec='*meas_ptcld*';
GEMsavepath=append(ptcldpath,'/GEMs/meas');
figpath=append(ptcldpath,'/Figures');
[meanGEMz_matrix,medGEMz_matrix,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,numframes,ptcldpath,GEMsavepath,figpath,spec);

% Create Transect Video
% ------------------------
ypick = [7.3 -20]; % in m, define locations to pick transects
yavg = 1.2; % in m, defining how wide in the alongshore to average over
%GEMpath=append(GEMsavepath,'meanGEMz_',num2str(1),'.mat');
GEMmatrixpath=append(GEMsavepath,'/meanGEMz.mat');
[v,quality_array]=tran_video(GEMmatrixpath,yavg,dxy,ypick,figpath);

% Plot Time Series of Various Transects at a Point
% -------------------------------------------------
TS_path='/Volumes/rsstu/users/k/kanarde/MasonBEAST/data/StormCHAZerz Data/Dec2023Noreaster_Processed/1702827001820/WaveTimeSeries';
pt=73;
timeseries73=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
pt=76;
timeseries76=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
pt=170;
timeseries170=GEMmatrix_to_timeseries(GEMmatrixpath,pt,dxy,ypick,yavg,TS_path);
% plot GEM and points on it (NEED TO ADD TO FUNCTION)



%% LATER 
% See Wave Propogate on One Plot (add in later)
% ------------------------------------------------
figure('units','inches','position',[1 1 10 3],'color','w');
hold on
plot(x(1,ixon:ix),zbeach(ixon:ix),'LineWidth',3,'Color',[148, 116, 27]/256);
hold on
for i = 1:4
    plot(x(1,ix:end),ztran(ix:end,i),'LineWidth',3) % plots water in swash
    plot(x(1,ixtran:ix),ztran(ixtran:ix,i),'LineWidth',3) % plots water in uprush
    %plot(x(1,:),ztran(:,i),'LineWidth',1.5) %Color',cmap(i,:)
end
legend("","t=1s","t=2s","t=3s","t=4s");
box on
clim([0 (numframes-1)/2])

% Total GEM Video (later)
% ------------------------------
format long g

% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;

% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);

GEMz=load(GEMmatrixpath);
GEMz=GEMz.meanGEMz;

% need to test this
v = VideoWriter(append(figfolder,'\GEMvidTEST'), 'MPEG-4');
v.FrameRate=2;%12;
v.Quality = 100;
open(v)

xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
fig=figure('units','inches','position',[0 0 10 6],'color','w');

for i = 1:numframes %+1
    clf
    pcolor(Xgrid,Ygrid,GEMz(:,:,i)); grid off;box on;
    box on
    ylim([0 4.5]);
    xlim([0 100]);
    clim([0 (numframes-1)/2])
    %colormap(cmap)
    %h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
    title(['$t$ = ',num2str(round(i/2)),' s'],'interpreter','latex','fontsize',ftsz(1));
    %scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    %text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    %scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    %text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    %scatter(rbrx,rbry,60,'fill','sq','g','MarkerEdgeColor','k');
    %text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'g');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    title(GEMtitle);
    pause(0.1)
    writeVideo(v,getframe(gcf))
end
close(v)