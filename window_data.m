function[win_data,nbursts, NFFT]=window_data(fs,window_length,data)
% function[win_data]=window_data(fs,window_length,data)
% 
% This function takes in a data array of one column of a creates a matrix
% of data, where each column represents a window of data.
% 
% INPUTS:
% --------
% fs = sampling rate (Hz)
% window_length = length of window (s)
% data = array of data (one column)
% 
% OUTPUTS:
% ---------
% win_data = windowed data matrix
% nbursts = number of bursts per window
% NFFT = number of FFT points
% 

window=window_length*fs; 
nbursts = floor(length(data)/window); % number of bursts in datasets (rounded to nearest whole number)
N=floor(log2(length(data)));
NFFT= 2^N; % number of FFT points
win_data= deal(NaN(window,nbursts)); % preallocate array


windct=0;
for i=1:nbursts
    win_data(:,i)=data((1+windct):(windct+window));
    windct=windct+window; % increase counter
end
