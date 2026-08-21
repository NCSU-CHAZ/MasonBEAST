function[transect_matrix]=DEM2transects(DEMzmatrixpath,camlocA,camlocB,rbrloc,dxy,ypick,yavg,savepath)
% function[transect_matrix]=DEM2transects(DEMzmatrixpath,dxy,ypick,yavg,savepath)
% % --------------------------------------------------------------------------
% This function takes in a DEM matrix mat file and pulls transects at
% given alongshore locations and stores them as a matlab structure. 
% The transects are averaged at X m resolution 
% per user specification. 
% Videos of each alongshore transect (in the swash zone) are gnereated and saved in subfolders,
% along with an assessment of points resolved per timestep.
% ----------------------------------------------------------------
% INPUTS:
% -------
% DEMzmatrixpath = path to DEM elevation mat file matrix
% camloc = camera A location coordinates [Ax,Ay]
% camloc = camera B location coordinates [Bx,By]
% rbrloc = rbr location coordinates [x,y]
% dxy = matrix bin size
% ypick = transect location(s) (m in the alongshore) 
% yavg = width of transect in alongshore direction (m)
% savepath = path to save figures and transect structure to
%
% OUTPUTS:
% --------
% transect_matrix = matrix of alongshore averaged transects save as .m file
%
% 
% EDITS NEEDED: 
% --------------------
% - 

% Last Edits: BG 08/17/2026
% -------------------------------------------------------------------------
%  CONSTANTS
% -------------------------------------------------------------------------
transavepath=append(savepath,'Transects/');

% -------------------------------------------------------------------------
% Load Pull Transects at Specified Alongshore Locations
% -------------------------------------------------------------------------
tran_struct=struct();

for i=1:length(ypick)
    current_y=ypick(i);
    %figpath=append(transavepath,'transect_y',num2str(current_y));
    fprintf('evaluating transect at ypick: %s ', num2str(current_y));

    [v,pts,ztran_matrix]=tran_video(DEMzmatrixpath,epochnum,yavg,dxy,current_y,transavepath); 
    tran_struct(i).name=append('transect_y',num2str(current_y));
    tran_struct(i).data=ztran_matrix;
end

save(append(transavepath,'alongshore_transects.mat'),'tran_struct');% save transect structure
fprintf('saved transect tructure');


% -------------------------------------------------------------------------
% Plot Transects in New DEM
% -------------------------------------------------------------------------
% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;
% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
%[Xgrid,Ygrid] = meshgrid(gridX,gridY);
% camera rotation and transformation
[CamAx, CamAy] = rotateCoordinates(camlocA(1), camlocA(2), Xloc, Yloc, rotang);
[CamBx, CamBy] = rotateCoordinates(camlocB(1), camlocB(2), Xloc, Yloc, rotang);
% rbr rotation and transformation
[rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);

% matrix for plotting
cat_trans=cat(3,tran_struct.data); % (pts, tsteps, transect name)
combined_tran_matrix=permute(cat_trans,[3,1,2]); % (transect name, pts, tsteps)
[Xgrid,Ygrid] = meshgrid(gridX,ypick);
xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';

for i=1:size(combined_tran_matrix,3) % loop over timesteps
    Z=combined_tran_matrix(:,:,i); % all transects at time t
    if size(Z,1) ~=length(ypick)
        Z=Z';
    end

    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    h=pcolor(Xgrid,Ygrid,Z); 
    grid off;box on;hold on;
    set(h,'EdgeColor','none');
    scatter(CamAx,CamAy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamAx(1)+0.5, CamAy(1), 'Cam A', 'FontSize', 12, 'Color', 'm');
    scatter(CamBx,CamBy,60,'fill','sq','m','MarkerEdgeColor','k');
    text(CamBx(1)+0.5, CamBy(1), 'Cam B', 'FontSize', 12, 'Color', 'm');
    scatter(rbrx,rbry,60,'fill','sq','k','MarkerEdgeColor','k');
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 12, 'Color', 'k');
    
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw); hc.Label.String = 'Elevation (m NAVD83 (2011))';
    set(gca,'fontsize',14);

    DEMtitle=sprintf('2m Alongshore Avg DEM (%s) - Timestep %d', num2str(epochnum),i);
    title(DEMtitle)
    figpath=append(transavepath,'Alongshore_avg_DEM_',num2str(i));
    saveas(fig,figpath,'png');
    fprintf('saved alongshore avg transect DEM: %s to %s\n',DEMtitle,figpath);

    close(fig);

end

