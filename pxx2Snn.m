function[Snn]=pxx2Snn(pxx,rho,nbursts)
%function[Snn]=pxx2Snn(pxx,rho,g)
% 
% This function converts pressure head power spectra to elevation power
% spectra.

% INPUTS:
% -------
% pxx = matrix of windowed power head pressure spectra data (Pa^2/Hz)
% rho = density (kg/m^3)
% nbursts = number of bursts
% 
% OUTPUTS:
% --------
% Snn = elevation power spectra (m^2/Hz)

g=9.81; % gravitational constant (m/s^2)

% calculate power spectra using multitaper method
for i=1:nbursts
    Snn(:,i)=pxx(:,i)/((rho^2)*g^2);
end
