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
%               - 'all'
% savepath = path to save figures and files to
% 
% OUTPUTS: 
% ---------------------
% density = array of density values per geomorphic region requested
% region_fig = GEM with regions outlined 
% 



