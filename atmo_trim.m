function[Patmos]=atmo_trim(atmospress, atmostime,RBR_structure,start_time,end_time)
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
% 
% OUTPUTS:
% ---------
% Patmos = vector of trimmed, linearly interpolated atmospheric
% pressure data in dbar

% Constants
mbar2dbar=0.01; 

% Convert start and end time to dateime
format_in='yyyy-mm-dd-HH-MM';
start_time=datenum(start_time,format_in);
end_time=datenum(end_time,format_in);


% Trim to start and end times
time_trim=find(atmostime>=start_time & atmostime<= end_time);
atmostime=atmostime(time_trim);
atmospress=atmospress(time_trim);

% convert from mbar to dbar
atmospress=atmospress.*mbar2dbar;

% linearly interpolate vector to be the same length as RBR pressure data
Patmos = interp1(atmospress,RBR_structure.data.values);
