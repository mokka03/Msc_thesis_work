%% Magnonic Crystals: From Simple Models toward Applications Jaros?awW. K?os and Maciej Krawczyk (pg. 288)
% out-of-plane magnetization
clear;

%% Parameters
YIG331 = struct;
YIG331.A_exch = 3.65E-12; % exchange coefficient (J/m)
YIG331.d = 100e-9;  % film thickness
YIG331.theta = 0;  % H0 angle, 0 if oop
YIG331.Bext = 215e-3;  % applied field in T

%% Measured data

dose_levels = (0:8)*1e12;
dose_levels = dose_levels/1e12;
load('data/Wavelength.mat')
YIG331.lambda = lambda;
YIG331.f0 = ones(size(YIG331.lambda))*2.1e9;

figure(1)
plot(dose_levels,YIG331.lambda*1e6,'o-')

xlabel('Ion Dose (10^{12} ions/cm^2)');
ylabel('Wavelength (\mum)');
% title({['Dependence of wavelength'] ['upon Ga+ ion dose']});

set(gca,'FontSize',15);

% SaveFig('FIBregions2100MHz_10dbm_H92/figure/','Wavelength',gcf);

%% Calculate Ms levels
% Get Ms for irradiated areas
Ms_levels = getMsLevels(YIG331);

%% Dose - Ms plot
[dose_Ms_fit,~,mu] = polyfit(dose_levels,Ms_levels,4);
dosefit = (0:0.1:9);
Ms_from_fit = polyval(dose_Ms_fit,dosefit,[],mu);

figure(3)
plot(dose_levels,Ms_levels*1e-3,'.')
hold on;
plot(dosefit,Ms_from_fit*1e-3)
hold off;
xlabel('Ion Dose (10^{12} ions/cm^2)');
ylabel('Saturation Magnetization (kA/m)');
% title({['Dependence of Ga+ ion dose upon'] ['saturation magnetization']});
xlim([min(dosefit) max(dosefit)]);
set(gca,'FontSize',15);

% SaveFig('FIBregions2100MHz_10dbm_H92/figure/','dose_Ms_fit',gcf);

%%
min(min(Ms_levels))
max(max(Ms_levels))
dose_Ms_data = [dose_levels*1e12; Ms_levels];
% save('FIBregions2100MHz_10dbm_H92/data/dose_Ms_data.mat','dose_Ms_data');