%% TC8 Processing Script 
% BG 06/02/2025
% This script reads atmospheric and RBR data and saves RBR guage data.

% Paths
addpath('/Users/bagaenzl/Desktop/rbr-rsktools-7a76410a599a');
genpath='/Volumes/kanarde/MasonBEAST/data/';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
stormpath=append(genpath,'StormCHAZerz Data/'); % path to storm Chazers data
pcpath=append(genpath,'/PointClouds/','TC8/'); % path to pointclouds
atmospath=append(stormpath,'AtmosPressure/','TC8_MaseN_atmospressure.csv'); % path to atmospheric pressure data from CORMP - MSNB_N
RBRpath=append(stormpath,'RBRs/','209244_20240918_1723.rsk'); % path to RBR data
savepath=append(stormpath,'TC8_processed');

% Constants
dbar2Pa = 1000; % conversion from dbar to Pa
fs=16; % sampling frequency
window_length=3600; % one hour windows (s)
nw = 6; % bandwidth product
taper_type='slepian'; % slepian tapers
h_inst=0; % height of instrument (m) (0 because it's bottom mounted)

% Read in data
RBR_structure=RSKopen(RBRpath);
RBR_structure=RSKreaddata(RBR_structure); 

atmos=readtable(atmospath,'NumHeaderLines',6);
atmospress=table2array(atmos(:,2)); % atmos pressure in dbar every 6 minutes
atmostime=table2array(atmos(:,1));

% Trim and interpolate atmospheric pressure data, save
start_time='2024-09-05-16-48-15'; % sept 5 2024, 4:36:15 PM (start of bucket test)
end_time='2024-09-18-17-22-00'; % sept 18 2024, 5:22:00 PM (end of final bucket test)
[Patmos]=atmo_trim(atmospress,atmostime,RBR_structure,start_time,end_time);
clear atmos atmospress atmostime;

% Calculate guage pressure and save
RBR_structure=RSKderiveseapressure(RBR_structure,'patm',Patmos); % (remove atmospheric pressure)
RBR_pressure = RBR_structure.data.values(:,2).*dbar2Pa; % guage pressure (Pa)
RBR_gpress_filepath=fullfile(savepath,[append('RBR_guage_09052024_09182024','.mat')]);
save(RBR_gpress_filepath,"RBR_pressure") ;

% Calculate dynamic pressure and save
[dyn_press,MWL]=dynpress(RBR_pressure,h_inst);
RBR_dpress_filepath=fullfile(savepath,[append('RBR_dyn_09052024_09182024','.mat')]);
save(RBR_dpress_filepath,"dyn_press") ;
clear RBR_pressure;

% Calculate pressure power spectra
[P_window,nburst,NFFT]=window_data(fs,window_length,dyn_press); % window guage pressure data
[pxx,f]=calc_spectra(P_window,nw,nburst,NFFT,fs,taper_type); % calculate power spectra
