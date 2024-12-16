%% Analyze data from buried pressure transducer
% BG 12/12/2024
clc
clear all
close all
format long g
%% Load Tools
addpath(genpath('/Users/bagaenzl@CCEE-LT-229/Desktop/MasonBEAST Data/MatlabTools/'));
% folder contains RSKtools and CIRN-Quantitative-Coastal-Imaging-Toolbox
%% Load in data
rsk_filename='209244_20240918_1723.rsk';
rsk_filepath=['/Users/bagaenzl/Desktop/MasonBEAST Data/RBR Data/',rsk_filename];
rsk=RSKopen(rsk_filepath); % RSK tools not working
Pseries=rsk;
% load atmos data (need to find)
% start and end dates

%% Define Variables
srho=; % Sea water density
ptH=; % height of PT?
sgrav=9.81; %[m/s^2] gravity constant
T=; % wave period
%% Calculate Mean Water Depth
 % dynamic pressure
dynP=Pseries-atmP;
meanPseries=mean(Pseries);
% mean water level
meanH=(meanPseries/(srho*sgrav))+ptH; 

%% Account for Attenuation
% Linearly interpolate burial depth
ptZ=; % interpolated matrix
% Calculate deep water wave length
deepL=(sgrav*T^2)/(2*pi);
% Initial wave number guess
k0=(2*pi)/deepL;
% Solve for wave number, k, with initial guess and mean water depth
k=((k0*meanH)*(1+(k0*meanH)^(1.3)-exp(-(1.1+2.0*k0*meanH))))/(meanH*sqrt(tanh(k0*meanH)));
% Calculate pressure response factor, Kp
Kp=cosh(k*ptH)./cosh(K.*meanH);
% Calculate depth attenuation factor (e^kz)
ekz=exp(k.*ptZ);
