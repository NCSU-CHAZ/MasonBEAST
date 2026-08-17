function[transect_matrix]=DEM2transects(DEMzmatrixpath,dxy,ypick,yavg,savepath)
% function[transect_matrix]=DEM2transects(DEMzmatrixpath,dxy,ypick,yavg,savepath)
% % --------------------------------------------------------------------------
% This function takes in a DEM matrix mat file and pulls transects at
% given alongshore locations and stores them as a matlab structure. 
% The transects are averaged at X m resolution 
% per user specification. 
% Videos of each alongshore transect are gnereated and saved in subfolders,
% along with an assessment of points resolved per timestep.
% ----------------------------------------------------------------
% INPUTS:
% -------
% DEMzmatrixpath = path to DEM elevation mat file matrix
% yavg = width of transect in alongshore direction (m)
% dxy = matrix bin size
% ypick = transect location(s) (m in the alongshore) 
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

