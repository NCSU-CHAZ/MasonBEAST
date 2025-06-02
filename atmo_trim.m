function[trimmed_atmos]=atmo_trim(atmospress, atmostime,RBR_structure,start_time,end_time,savepath)
%function[trimmed_atmos]=atmo_trim(atmospress,atmostime,RBR_structure,start_time,end_time,savepath)
% 
% This function takes in atmospheric pressure data and converts it to dbar,
% along with linearly interpolating the data to be the same length as the
% RBR pressure data. 
% 
% INPUTS:
% ---------
% atmospress = array of atmospheirc pressure data in mbar
% atmostime = array of time points correlating to atmospress
% RBR_structure = structure from RSKreaddata
% start_time = collection start time in string format 'yyyy-mm-dd-HH-MM'
% end_time = collection end time in string format 'yyyy-mm-dd-HH-MM'
% savepath = path to save trimmed atmos data
% 
% OUTPUTS:
% ---------
% trimmed_atmos = vector of trimmed, linearly interpolated atmospheric
% pressure data in dbar


