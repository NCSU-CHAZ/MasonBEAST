function[ztran]=transect_pull(GEMpath,yavg,dxy,ypick,savepath,figpath)
% function[ztran]=transect_pull(GEMpath,dxy,transect_loc)
% 
% This function pulls a transect for each GEM given and saves it as a mat
% file. It also creates a plot and saves it as a figure.
% 
% INPUTS:
% --------
% GEMpath = path to GEM mat file
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

% load GEM elevation values
GEM=load(GEMpath);

% create and rotate grid
gridX=0:dxy:110;
gridY=-80:dxy:25;
[x,y]=meshgrid(gridX,gridY);
[xrot,yrot]=rotateCoordinates(GEM(:,1),GEM(:,2),Xloc,Yloc,rotang);

% load GEM elevation values
GEM=load(GEMpath);

[~,iy]=min(abs(y(:,1)-ypick(1)));
ztran=median(GEM(iy-iyavg:iy+iyavg,:,:),1,'omitnan');
ztran=squeeze(ztran);
ztran=movemean(ztran,2,1,'omitnan');
save(savepath,'ztran','-mat');
