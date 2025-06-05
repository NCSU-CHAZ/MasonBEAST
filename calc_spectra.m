function[pxx,f]=calc_spectra(win_data,nw,nbursts,NFFT,fs,taper_type)
% function[pxx,f]=calc_spectra(data,window_length,fs,taper_type)
% 
% This function takes a windowed pressure time series and calculates the
% power spectra using the multitaper method. 
% 
% INPUTS:
% --------
% win_data = windowed data matrix 
% nw = bandwidth product
% nbursts = number of bursts per window
% NFFT = number of FFT points (2^n = # of sampling points)
% fs = sampling frequency
% taper_type = string of sepcified taper type (ex: 'Slepian')
% 
% OUTPUTS:
% ---------
% pxx = power spectra
% f = frequency

% preallocate arrays
pxx=deal(NaN((NFFT/2)+1,nbursts));
f=deal(NaN((NFFT/2)+1,nbursts));

% calculate power spectra using multitaper method
for i=1:nbursts
    [pxx(:,i),f(:,i)]=pmtm(win_data(:,i),nw,NFFT,fs,'Tapers',taper_type);
end

