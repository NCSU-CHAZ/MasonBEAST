function[density,region_fig]=GEM_region_density_calc(GEMpath,region,savepath);
% function[density,region_fig]=GEM_region_density_calc(GEMpath,region,savepath);
% 
% This function takes in a GEM .mat file of gridded elevation values and
% calculates the normalized density for the requested geographic region of
% the GEM. The normalized density is defined as (# of points - # of NaNs)/#
% of points. The GEMs are on a 0.2 m grid. 
% 
% 
% INPUTS:
% --------------------
% GEMpath = path to GEM .mat file
% region = string of region of interest, the options are:
%               - 'dune'
%               - 'RBR'
%               - 'upperbeachface'
%               - 'lowerbeachface'
%               - 'shoreline'
% savepath = path to save figures and files to
% 
% OUTPUTS: 
% ---------------------
% density = array of density values per geomorphic region requested
% region_fig = GEM with regions outlined 
% 
% ------------------------------------------
format long g

% Load GEM
meanGEMz=load(fullfile(GEMmat_path,'meanGEMz.mat'));
meanGEMz=[meanGEMz.meanGEMz];

% GEM name and date
GEMname=split(GEMmat_path,'/');
GEMname=GEMname(10,1);
GEMname=string(GEMname);

GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMdate=string(GEMdate);
GEMtitle=append(GEMname,',',GEMdate);
% ------------------------------------------
% Constants
dxy=0.2; % m
numframes=length(meanGEMz,3);

% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;

% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
% --------------------------------------------
% Set region indices
if region == 'dune'
        rows=1:9;
        cols=1:9;
        rect_x=3;
        rect_y=3;
        width=2;
        height=2;
    elseif region == 'upperbeachface'

    elseif region == 'lowerbeachface'

    elseif region == 'shoreline'

    elseif region == 'RBR'
        rbrloc = [239779.25, 3784738.325]; % location of RBR in swash
        % rbr rotation and transformation
        [rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);

end 

% Calculate Normalized Density per Region and plot GEM with region outline
for i=1:numframes
    % calculate normalized density
    indexed_GEM=meanGEMz(rows,cols,i); % grab portion of GEM
    numNans=nansum(indexed_GEM); % number of NaNs in section
    pts=sum(indexed_GEM)-numNans; % number of points resolved
    density=pts/sum(indexed_GEM); % normalized density of region

    % plot
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanGEMz(:,:,i)); grid off;box on;hold on
    scatter(rbrx,rbry,70,'fill','sq','g','MarkerEdgeColor','k'); hold on;
    rectangle('Position',[rect_x,rect_y,width,height],'EdgeColor','k','Linewidth',2); % plot region outline
    hold on; text(rect_x+1,rect_y+1,region,'FontSize',16,'Color','k'); % label region
    text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 16, 'Color', 'g'); % label rbr
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    title(GEMtitle);
    density_txt=num2str(density); density_txt=append('Normalized Density = ', density_txt); annotation('textbox',[0.531589801274837,0.07001239157373,0.100000000000001,0.2],'String',density_txt,'EdgeColor','none','FontSize',28);
    filename=append('density_',region,'_',GEMname,'_',string(i));
    figpath=fullfile(figpath, filename);
    saveas(fig,figpath,'png');
    close(fig);

end




