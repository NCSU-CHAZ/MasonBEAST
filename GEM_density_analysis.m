% GEM region density analysis
% BG 6/25/25

% TESTING CODE ------------
% Paths
genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data;
CAM_analysispath=append(genpath,'/GEMs/Camera_Location_Analysis'); % path to camera analysis files  
measured_path=append(CAM_analysispath,'/Measured/'); % measured GEMs
metashape_path=append(CAM_analysispath,'/Metashape/'); % metashape GEMs
GEMpath_meta=append(metashape_path,'/1723489201189/');
GEMpath_meas=append(measured_path,'/1723489201189/');


% plot GEM to ID regions
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanGEMz(:,:,1)); grid off;box on;hold on
    %scatter(GCPx,GCPy,60,'fill','r','MarkerEdgeColor','k') (need to figure
    %out how to specifiy GCPs/do we need?)
    %scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    %text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    %scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    %text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    scatter(rbrx,rbry,60,'fill','sq','g','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'g');
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    %title(GEMtitle);
    %filename=append('mean',GEMname,'_',string(i));
    %meanfigpath=fullfile(figpath, filename);
    %meanfigpath=append(figpath,"/mean",GEMname);
    %saveas(fig,meanfigpath,'png');
    %close(fig);