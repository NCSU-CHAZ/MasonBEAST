function[pxx,f]=calc_spectra(data,nw,nbursts,taper_type)
% function[pxx,f]=calc_spectra(data,window_length,fs,taper_type)
% 
% This function takes a windowed pressure time series and calculates the
% power spectra using the multitaper method. 
% 
% INPUTS:
% --------
% data = data matrix 
% nw = bandwidth product
% nbursts = number of bursts per window
% taper_type = string of sepcified taper type (ex: 'Slepian')
% 
% OUTPUTS:
% ---------
% pxx = power spectra
% f = frequency

% constants
