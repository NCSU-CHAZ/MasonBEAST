%% Plotting Monthly DEMS Over Time
% read and create plots transects
clc
clear all
close all
format long g

%addpath(genpath('/Users/cmbaker9/Documents/Research/Lab_Experiments/codes/CIRN-Quantitative-Coastal-Imaging-Toolbox'))
%addpath(genpath('/Users/cmbaker9/Documents/Programs/MTOOLS'))
%addpath('/Users/cmbaker9/Documents/Research/MasonBEAST/code/imagery/support_routines')

%% STEP 1: Create paths, files and naming

% general path and names
datapath    = '/Users/bagaenzl/Desktop/MasonBEAST Data/transects/';
% trialname   = '1699459201959';%'1702827001820';%'1698951602001';
% fname = [datapath,'stereo/dem/',trialname,'/',trialname,'_dem1.tif'];
% gcpname = [datapath,'gcps/GCPs_02_15_2024.txt'];
% figfolder = '/Users/cmbaker9/Documents/Research/MasonBEAST/figures/';

%% Call previous months
% 2023
Zsept=load([datapath,'ZSept.mat']);
Zoct=load([datapath,'ZOct.mat']);
Znov=load([datapath,'ZNov.mat']);
Zdec=load([datapath,'ZDec.mat']);
% 2024
Zjan=load([datapath,'Zjan.mat']);
Zfeb=load([datapath,'ZFeb.mat']);
Zmar=load([datapath,'Zmar.mat']);
Zapr=load([datapath,'ZApr.mat']);
Zmay=load([datapath,'Zmay.mat']);
Zjun=load([datapath,'Zjune.mat']);
Zjul=load([datapath,'Zjuly.mat']);
Zaug=load([datapath,'Zaug1.mat']);
Zsept_preTC8=load([datapath,'Zsep_predurTC8.mat']);
Zsept_postTC8=load([datapath,'Zsept19_postTC8.mat']);
% there is another storm sept27
% no image collection sept 27-oct 7
Zoct_7=load([datapath,'Zoct7.mat']);
Zoct_prestorm=load([datapath, 'Zoct10_prestorm.mat']);
Znov_24=load([datapath, 'Znov_24.mat']);
% storm
Zprestorm=load([datapath,'ZDec2.mat']);
Zpoststorm=load([datapath,'ZDec3.mat']);

%% Combine into a matrix

% combine
Z(:,1) = Zsept.ZSept;  % September
Z(:,2) = Zoct.ZOct;   % October
Z(:,3) = Znov.ZNov;   % November
Z(:,4) = Zdec.ZDec;   % December
Z(:,5) = Zjan.ZJan;   % January
Z(:,6) = Zfeb.ZFeb;   % February
Z(:,7) = Zmar.ZMar;   % March
Z(:,8) = Zapr.ZApr;   % April
Z(:,9) = Zmay.Z;   % May
Z(:,10) = Zjun.Z;  % June
Z(:,11) = Zjul.Z;  % July
% GEM method transects
Zgem(:,1) = Zaug.ztran(:,1);  % August
Zgem(:,2) = Zsept_preTC8.ztran(:,1); % sept pre storm
Zgem(:,3) = Zsept_postTC8.ztran(:,1); % sept post storm
Zgem(:,4) = Zoct_prestorm.ztran(:,1); % October
Zgem(:,5) = Znov_24.ztran(:,1); % November

Zstorm(:,1) = Zprestorm.ZDec2;
Zstorm(:,2) = Zpoststorm.ZDec3;

%% X meters across beach
x=0:0.15:(550*0.15); % need to check this out

%%

ftsz = [24 20];
lw = 1.5;
ixlim = [-120,0];
iylim = [-100,0];
% iylim = [-13,13];
iclim = [-0.201 0.201];
tickminor = 'on';
tickdir ='in';
ticklen = 1;
xlab = 'Cross-Shore (m)';
ylab = 'Elev. (m)';
%cmap = slanCM('thermal-2',12);%cmocean('phase',13);

ax1pos = [0.1 0.2 0.6 0.6];
%cmapplt = flipud(cmocean('deep',20));
% cmapplt = cmapplt(3:end,:);

figure('units','inches','position',[1 1 12 6],'color','w');
ax1 = axes('Position',ax1pos);
hold on
for i = 5:size(Z,2)
    if i == 2
        [~,ix] = min(abs(x-27));
        plot(x(1:ix),Z(1:ix,i),'LineWidth',3)
    else
        plot(x,Z(:,i),'LineWidth',3)
    end
end
plot(x,Zstorm(:,1),'LineWidth',3,'Color','k')
plot(x,Zstorm(:,2),'LineWidth',3,'Color',[0.3 0.6 0])

box on
ylim([0 4.5]);
xlim([0 40]);
clim([0 11])
% set(ax2,'YDir','normal','Color','k')
colormap(cmap)
h1 = plotstyleCMB(gca,xlab,ylab,ftsz,ticklen,lw,tickminor,tickdir);
% set(h1,'yticklabel','');
% rectangle(ax2,'Position',[27 11.8 1.4 1.5],'EdgeColor','none','FaceColor',[1 1 1 0.8])
% text(-118,-7,'(c)','interpreter','latex','fontsize',ftsz(1),'Color','k');
% hc = colorbar('Location','eastoutside','Position', [0.93 0.225 0.03 0.4],'orientation','vertical','YAxisLocation','right');
%     set(hc,'fontsize',ftsz(2),'linewidth',lw);
% text(164,3.3,'$t$ (s)','interpreter','latex','fontsize',ftsz(1));
% leglab = {'Sept','Oct','Nov','Dec','Jan','Feb','Mar','May','Apr','Jun','Jul','Pre-Storm','Post-Storm'};
% legend(leglab, 'Position', [0.67 0.45 0.2 0.1],...
%     'EdgeColor','none','Color','none','interpreter','latex','fontsize',ftsz(2));

    sname = 'monthly_transects_7';
    print([figfolder,sname],'-dpng')
    exportgraphics(gcf, [figfolder,sname,'.png']);
    exportgraphics(gcf, [figfolder,sname,'.pdf']);

%% Plot Data
figure(1)
clf
plot(cross_dist,Znov.ZNov(1,7:end),'r','Linewidth',1);
hold on
plot(cross_dist,Zdec.ZDec(1,7:end),'Color',[225/255 165/225 0],'Linewidth',1);
hold on
plot(cross_dist,Zjan.ZJan(1,7:end),'y','Linewidth',1);
hold on
plot(cross_dist,Zfeb.ZFeb(1,7:end),'g','Linewidth',1);
hold on
plot(cross_dist,Zmar.ZMar(1,7:end),'c','Linewidth',1);
hold on
plot(cross_dist,Zapr.ZApr(1,7:end),'b','Linewidth',1);
hold on
plot(cross_dist,Zmay.Z(1,7:end),'Color',[144/255 103/255 167/255],'Linewidth',1);
hold on
plot(cross_dist,Zjune.Z(1,7:end),'Color',[1 0.5 0.8],'Linewidth',1);
hold on
plot(cross_dist,Zjuly.Z(1,7:end),'k','Linewidth',1);
hold on
legend("November","December","January","February","March","May","April","June","July");

ylabel("Elevation (m), NAVD88",'Fontsize',16);
xlabel("Cross Shore Distance (m)",'Fontsize',16);
hold on
%ylim([0 4.5]);
%yticks([0.5 1 1.5 2 2.5 3 3.5 4]);
%xlim([0 50]);
%xticks([0 5 10 15 20 25 30 35 40 45 50]);
hold on
title("Transects Over Time",'Fontsize',16);

%% Plot Individual Months
figure(2)
clf
plot(cross_dist,Zjuly.Z(1,7:end),'Color',[144/255 103/255 167/255],'Linewidth',2);
hold on
ylabel("Elevation (m), NAVD88",'Fontsize',16);
xlabel("Cross Shore Distance (m)",'Fontsize',16);
hold on
title("7/04/2024","Transect",'Fontsize',16);

% Plot August to Nov 2024
figure(2)
clf
plot(x,Zgem(:,1),'r','Linewidth',3.5);
hold on
plot(x,Zgem(:,2),'color',[0.85 0.3 0.09],'Linewidth',3.5);
hold on
plot(x,Zgem(:,3),'y','Linewidth',3.5);
hold on
plot(x,Zgem(:,4),'g','Linewidth',3.5);
hold on
plot(x,Zgem(:,5),'b','Linewidth',3.5);
legend("August","September Pre TC8","September Post TC8","October","November");
ylabel("Elevation (m), NAVD88",'Fontsize',16);
xlabel("Cross Shore Distance (m)",'Fontsize',16);
hold on
title("Transects Aug to Nov",'Fontsize',16);