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
zbeg=0; % using starting RBR point as bed surface
zend=1.639667-1.66567; % initial top point - final top point (m)
rho =1024; % density of seawater (kg/m^3)
g=9.81;

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
[dyn_press,MWL]=dynpress(RBR_pressure,h_inst); % Pa
RBR_dpress_filepath=fullfile(savepath,[append('RBR_dyn_09052024_09182024','.mat')]);
save(RBR_dpress_filepath,"dyn_press") ;
clear RBR_pressure;

% Plot Pressure Head (m)
time = datetime(RBR_structure.data.tstamp, 'ConvertFrom', 'datenum'); 
format_in='yyyy-mm-dd-HH-MM-ss';
start_time='2024-09-16-00-00-00';
end_time='2024-09-17-00-00-00';
start_time=datenum(start_time,format_in);
end_time=datenum(end_time,format_in);

xlabs=[datenum('2024-09-16-01-00-00',format_in),datenum('2024-09-16-03-00-00',format_in),datenum('2024-09-16-05-00-00',format_in),datenum('2024-09-16-07-00-00',format_in),datenum('2024-09-16-09-00-00',format_in),datenum('2024-09-16-11-00-00',format_in)datenum('2024-09-16-13-00-00',format_in),datenum('2024-09-16-15-00-00',format_in),datenum('2024-09-16-17-00-00',format_in),datenum('2024-09-16-19-00-00',format_in),datenum('2024-09-16-21-00-00',format_in),datenum('2024-09-16-23-00-00',format_in)];
xlabs=datetime(xlabs,'ConvertFrom','datenum');

fig=figure(1); plot(RBR_structure.data.tstamp,dyn_press./(rho*g),'Linewidth',2);datetick('x'); xlim([start_time end_time]); xlabel('Time');ylabel('Pressure Head (m)');
hold on;hold on; yline(zbeg,'-','Before storm bed level','Linewidth',1.5); hold on;
yline(zend,'-','Post storm bed level','Linewidth',1.5); set(gca,'Fontsize',16);%yline(MWL,':','Mean Water Level','Linewidth',1.5);
title('Pressure Head Time Series');legend('Pressure Head', 'Before storm bed level', 'Post storm bed level');


% Calculate pressure power spectra
[P_window,nburst,NFFT]=window_data(fs,window_length,dyn_press); % window guage pressure data
[pxx,f]=calc_spectra(P_window,nw,nburst,NFFT,fs,taper_type); % calculate power spectra (Pa^2/Hz)

% Elevation spectra and depth attenuation
Snn=pxx2Snn(pxx,rho,nburst); % (m^2/Hz)
[Snn_d,kp,ekz]=depth_att(Snn,f,zbeg,zend,h_inst,MWL,nburst,NFFT,savepath); % this is wrong NaNs after row 6

% Plotting function (need to write)
fig=figure(1); subplot(2,1,1);plot(Snn_d(2,:),'Linewidth',2);ylabel('Elevation spectra (m^2)'); 
hold on; subplot(2,1,2); plot(dyn_press./(rho*g),'Linewidth',2); xlabel('Sample Number');ylabel('Dynamic Pressure (m)');

fig=figure(2); plot(Snn_d(2,:),'Linewidth',2);ylabel('Elevation spectra (m^2)'); xlabel('Hours');

fig=figure(3); plot(dyn_press./(rho*g),'Linewidth',2); xlabel('Sample Number');ylabel('Dynamic Pressure (m)');