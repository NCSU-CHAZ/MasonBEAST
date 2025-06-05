function[Snn_d,kp,ekz]=depth_att(Snn,f,zbeg,zend,h_inst,MWL,savepath)
% function[Snn_d]=depth_att(Snn)
% 
% This function corrects the elevation spectra for depth attenuation.
% 
% INPUTS:
% --------
% Snn = windowed eleavtion spectra (m^2/Hz)
% f = windowed frequency values
% zbeg = depth below the bed surface at the begining of the given time
% series (positive in downward direction) (m)
% zend = depth below the bed surface at the end of the time series (m)
% h_inst = height of instrument (PT = 0 since bottom mounted)
% MWL = mean water level (m)
% savepath = path to save correct Snn to
% 
% OUTPUTS:
% ---------
% Snn_d = windowed depth attenuated elevation spectra
% kp = depth attenuation factor
% ekz = depth attenuation factor by burial (e^kz)
