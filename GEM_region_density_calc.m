function[density]=GEM_region_density_calc(GEMpath,region,region_specs,savepath)
% function[density,region_fig]=GEM_region_density_calc(GEMpath,region,savepath);
% 
% This function takes in a GEM .mat file of gridded elevation values and
% calculates the normalized density for the requested geographic region of
% the GEM. The normalized density is defined as (# of points - # of NaNs)/#
% of points. The GEMs are on a 0.5 m grid. 
% 
% 
% INPUTS:
% --------------------
% GEMpath = path to GEM .mat file
% region= string of region of interest, the options are:
%               - 'dune'
%               - 'RBR'
%               - 'upperbeachface'
%               - 'lowerbeachface'
%               - 'shoreline'
% region_specs = array of region locations:
%               [cols, rows, rect_x, rect_y, width, height]
%               - cols and rows are grabbed from the original GEM
%               - rect_x and rect_y are the lower left corner of plotted
%               rectangle boundary
%               - width and height are the width and height of the plotted
%               rectangle boundary
% savepath = path to save figures and files to
% 
% OUTPUTS: 
% ---------------------
% density = density values of geomorphic region requested
% 
% ------------------------------------------
format long g

% Load GEM
meanGEMz=load(fullfile(GEMpath,'meanGEMz.mat'));
meanGEMz=[meanGEMz.meanGEMz];

% GEM name and date
GEMname=split(GEMpath,'/');
GEMname=GEMname(9,1);
GEMname=string(GEMname);

GEMdate=datetime(str2num(GEMname),'ConvertFrom','epochtime','TicksPerSecond',1000);
GEMdate=string(GEMdate);
GEMtitle=append(GEMname,',',GEMdate);
% ------------------------------------------
% Constants
dxy=0.5; % m
numframes=1;%size(meanGEMz,3);
% define local coordinate system origin and rotation angle
Xloc = 239737;
Yloc = 3784751;
rotang = 35;

rbrloc = [239779.25, 3784738.325]; % location of RBR in swash
% rbr rotation and transformation
[rbrx, rbry] = rotateCoordinates(rbrloc(1), rbrloc(2), Xloc, Yloc, rotang);


% create grid
gridX = 0:dxy:110;
gridY = -80:dxy:25;
[Xgrid,Ygrid] = meshgrid(gridX,gridY);
% --------------------------------------------
% Set region indices
GEMY=1:size(meanGEMz,2);
GEMX=1:size(meanGEMz,1);
x=Xgrid(1,:);
y=Xgrid(:,1);
y=1:size(Ygrid,2);
if  strcmp(region, 'dune')
        cols=1:21;
        rows=1:10;
        [row,col]=meshgrid(rows,cols);
        rect_x=0; % lower left corner
        rect_y=0; % lower left corner
        width=9; % 4.5 m
        height=20; % 10 m
    elseif strcmp(region, 'upperbeachface')
        cols=10:45;
        rows=9:33;
        [row,col]=meshgrid(rows,cols);
        rect_x=9; % lower left corner
        rect_y=-10; % lower left corner
        width=24; % 12 m
        height=35; % 17.5 m
    elseif strcmp(region, 'lowerbeachface')
        cols=-80:10;
        rows=9:33;
        [row,col]=meshgrid(rows,cols);
        rect_x=9; % lower left corner
        rect_y=-80; % lower left corner
        width=24; % 12 m
        height=70; % 35 m
    elseif strcmp(region, 'shoreline')
        cols=-80:24;
        rows=35:45;
        [row,col]=meshgrid(rows,cols);
        rect_x=35; % lower left corner
        rect_y=-80; % lower left corner
        width=10; % 5 m
        height=104; % 52 m
    elseif strcmp(region,'RBR')
        row=rbrx;
        col=rbry;
        rect_x=rbrx-1;
        rect_y=rbry-1;
        height=2;
        width=2;

end 

% Calculate Normalized Density per Region and plot GEM with region outline
for i=1:numframes
    % calculate normalized density
    indexed_GEM=interp2(Xgrid,Ygrid,meanGEMz(:,:,i),row,col);% meanGEMz(rows,cols,i); % grab portion of GEM
    numNans=nnz(isnan(indexed_GEM)); % number of NaNs in section
    pts=numel(indexed_GEM)-numNans; % number of points resolved
    density=pts/numel(indexed_GEM); % normalized density of region
    density_txt=num2str(density); density_txt=append('Normalized Density = ', density_txt);

    % plot
    xlab = 'Cross-Shore (m)';ylab = 'Alongshore (m)';
    fig=figure('units','inches','position',[0 0 10 6],'color','w');
    pcolor(Xgrid,Ygrid,meanGEMz(:,:,i)); grid off;box on;hold on
    scatter(rbrx,rbry,70,'fill','sq','g','MarkerEdgeColor','k'); hold on;
    rectangle('Position',[rect_x,rect_y,width,height],'EdgeColor','k','Linewidth',2); % plot region outline
    hold on; 
    if strcmp(region, 'dune')
        text(rect_x+1,rect_y+1,region,'FontSize',20,'Color','k'); % label region
        text(rect_x-4,rect_y-2,density_txt,'FontSize',20,'Color','k'); hold on;
        text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 20, 'Color', 'g'); % label rbr
    elseif strcmp(region,'upperbeachface')
        text(rect_x+1,rect_y+1,region,'FontSize',20,'Color','k'); % label region
        text(rect_x-4,rect_y-2,density_txt,'FontSize',20,'Color','k'); hold on;
        text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 20, 'Color', 'g'); % label rbr
    elseif strcmp(region,'RBR')
        text(rect_x+3,rect_y+1,region,'FontSize',20,'Color','k'); % label region
        text(rect_x-4,rect_y-2,density_txt,'FontSize',20,'Color','k'); hold on;
    else
        text(rect_x+1,rect_y+50,region,'FontSize',20,'Color','k'); % label region
        text(rect_x,rect_y+60,density_txt,'FontSize',20,'Color','k'); hold on;
        text(rbrx(1)+0.5, rbry(1), 'RBR', 'FontSize', 20, 'Color', 'g'); % label rbr
    end
    
    shading interp;
    axis equal;ylim([-60 30]); ylabel(ylab);xlabel(xlab);xlim([-10 90]);clim([0 3.8]);
    ftsz = [22 18]; lw = 1.2; hc = colorbar('Location','eastoutside','Position', [0.83 0.14 0.035 0.4],'orientation','vertical','YAxisLocation','right');
    set(hc,'fontsize',ftsz(2),'linewidth',lw);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]); % Enlarge figure to full screen
    title(GEMtitle);
    filename=append('density_',region,'_',GEMname,'_',string(i));
    figpath=fullfile(savepath, filename);
    saveas(fig,figpath,'png');
    close(fig);

end




