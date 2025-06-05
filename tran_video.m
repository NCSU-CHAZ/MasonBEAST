function[wave_video]=tran_video(ztranpath,figpath)
%function[wave_video]=tran_video(ztranpath)
% 
% This function takes in a matrix of transect elevation values per frame at the same
% location and creates a video showing timesteps of 1 second. The lowest
% transect values are used as the bed surface of the beach, whereas the
% rest are depicted as incoming waves.
% 
% INPUTS:
% ---------
% ztran = matrix of transect elevation values (m) [elevation values (z),frame number]
% figpath = path to save video to
% 
% OUTPUTS: 
% ----------
% wave_video = video of wave by wave interactions with the bed


