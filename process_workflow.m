% Workflow for processing full epoch folders
% Last edits: 05/27/2026 BG


% paths,epochnum, specs
% -------------------------------
% Remember to change campath AND spec before running
epochnum='1702827001820'; % epoch number
% paths 
%genpath='/Volumes/kanarde/MasonBEAST/data';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
%pcpath=append(genpath,'/PointClouds/'); % path to pointclouds
% paths if using external drive
genpath='/Volumes/Elements';
stormCHAZerspath=append(genpath,'/StormCHAZerz Data');
decNoreastherpath=append(stormCHAZerspath,'/Dec2023Noreaster_processed');
epochpath=append(decNoreastherpath,'/',epochnum);
ptcldpath=append(epochpath,'/ptclds');
figpath=append(epochpath,'/Figures/');
GEMsavepath=append(epochpath,'/GEMs/');
vidsavepath=append(epochpath,'/videos/');

% specs
camlocA = [239766.1, 3784761.9];%, 10.37
camlocB = [239759.4, 3784755.0];%, 10.26];
rbrloc = [239779.25, 3784738.325];
dxy = 0.2; % meter grid resolution size

% Process Pointclouds
% ------------------------
spec='1702827001820*_ptcld*';
[meanZ_matrix,medZ_matrix,Xrot,Yrot]=ptcld_to_GEM(camlocA,camlocB,rbrloc,dxy,ptcldpath,GEMsavepath,figpath,spec);

% Create Transect Video
% -------------------------
ypick = [7.3 -20]; % in m, define locations to pick transects
Zmatrixpath=append(GEMsavepath,'/meanMAPz.mat');
% 0.5 m alongshore avg
yavg = 1.2;
newfigpath=append(vidsavepath,'Transect_y7_3','/1_2yavg');
[v_12,quality_array12]=tran_video(Zmatrixpath,yavg,dxy,ypick,newfigpath);





