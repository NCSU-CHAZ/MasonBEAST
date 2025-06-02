%% TC8 Processing Script 
% BG 06/02/2025

% Paths
genpath='/Volumes/kanarde/MasonBEAST/data/';% path to Research storage /Volumes/kanarde-1/MasonBEAST/data /Volumes/rsstu/users/k/kanarde/MasonBEAST/data
stormpath=append(genpath,'StormCHAZerz Data/'); % path to storm Chazers data
pcpath=append(genpath,'/PointClouds/','TC8/'); % path to pointclouds
atmospath=append(stormpath,'AtmosPressure/','TC8_MaseN_atmospressure.csv'); % path to atmospheric pressure data from CORMP - MSNB_N
RBRpath=append(stormpath,'RBRs/','209244_20240918_1723'); % path to RBR data

% Constants
dbar2Pa = 1000; % conversion from dbar to Pa
fs=16; % sampling frequency
window_length=3600; % one hour windows (s)
nw = 6; % bandwidth product
taper_type='Slepian'; % slepian tapers

% Read in data
RBR_structure=RSKopen(RBRpath);
RBR_structure=RSKreadata(RBR_structure); 

atmos=readtable(atmospath);
atmospress=table2array(atmos(:,2)); % atmos pressure in dbar every 6 minutes
atmostime=table2array(atmos(:,1));


% Trim and interpolate atmospheric pressure data
start_time='2024-09-05-16-48-15'; % sept 5 2024, 4:36:15 PM (start of bucket test)
end_time='2024-09-18-17-22-00'; % sept 18 2024, 5:22:00 PM (end of final bucket test)
[Patmos]=atmo_trim(atmospress,atmostime,RBR_structure,start_time,end_time);

% Calculate guage pressure
RBR_structure=RSKderiveseapressure(RBR_structure,'patm',Patmos); % (remove atmospheric pressure)
RBR_pressure = RBR_structure.data.values(:,2).*dbar2Pa; % guage pressure (Pa)

% Calculate pressure power spectra
[P_window,nburst,NFFT]=window_data(fs,window_length,RBR_pressure); % window guage pressure data
[pxx,f]=calc_spectra(P_window,nw,nbursts,NFFT,fs,taper_type); % calculate power spectra
