function[Snn_d,kp,ekz]=depth_att(Snn,f,zbeg,zend,h_inst,MWL,nbursts,NFFT,savepath)
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
% nbursts = number of bursts per window
% NFFT = number of points
% savepath = path to save correct Snn to
% 
% OUTPUTS:
% ---------
% Snn_d = windowed depth attenuated elevation spectra
% kp = depth attenuation factor
% ekz = depth attenuation factor by burial (e^kz)
% --------------------------------------------------------------
% specs and constants
format long 
g=9.81; % gravitational constant m/s^2 

%--------------------------------------------------------
% linearly interpolate burial depth
zpts=[zbeg, zend]; % burial depth points
dxy=1:1/(nbursts-1):2; % bin size
z=interp1(zpts,dxy); % linearly interpolated burial depth
z=repmat(z,(NFFT/2)+1,1); % window matrix compatible

% ---------------------------------------------------------------
% Solve for wave number using dispersion relationship (not accounting for
% dopple effects)

% initiate k matrix
k=deal(NaN((NFFT/2)+1,nbursts));
kp=deal(NaN((NFFT/2)+1,nbursts));
ekz=deal(NaN((NFFT/2)+1,nbursts));
Snn_d=deal(NaN((NFFT/2)+1,nbursts));

% loop through windows and solve iteratively
for i=1:nbursts
    T(:,i)=1./f(:,i); % wave period (s)
    rad_f(:,i)=(2*pi)./T(:,i); % radial frequency 
    for j=1:10%length(T(:,i))
        % solve for wavelength iteratively
        Tcalc=T(:,i); % simplify calculations
        L0(j)=(g*Tcalc(j)^2)/(2*pi); % deep water wavelength (m)

        errortol=0.001; % error tolerance
        err=10; % error check
        Lguess=L0(j); % use deep water wavelength as initial guess
        ct=0; % initiate counter
        while err>errortol && ct<10000
            ct=ct+1; % count up
            L(j)=((g*Tcalc(j)^2)/(2*pi))*tanh((2*pi*MWL)/Lguess); % calculate new L
            err=abs(L(j)-Lguess); % check for error
            Lguess=L(j); % update guess
        end 
        L(j)=Lguess;
        k(j,i)=(2*pi)/L(j); % wave number windowed matrix (s^-1)

        % Calculate depth attenuation factor
        kp(j,i)=cosh(k(j,i)*h_inst)/cosh(k(j,i)*MWL);
        % Calculate depth attenuation by burial
        ekz(j,i)=exp(k(j,i)*z(j,i));
        % calculate corrected surface elevation spectrum
        Snn_d(j,i)=Snn(j,i)*((ekz(j,i)^2)/(kp(j,i)^2));

    end
    
end

% ---------------------------------------------------
% save corrected elevation spectrum
Snn_d_filepath=fullfile(savepath,[append('Snn_d_09052024_09182024','.mat')]);
save(Snn_d_filepath,"Snn_d") ;

