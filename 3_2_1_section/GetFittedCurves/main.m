clear;

%% Fit curve to dMs-dose data
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\YIG_samples\YIG328\data\dose_Ms_data_2050MHz.mat')

%%
dose_levels = dose_Ms_data(1,1:end);
Ms_levels = dose_Ms_data(2,1:end);

[Ms_dose_fit,mu_Ms,dose_Ms_fit,mu_dose] = getFittedCurves(dose_levels,Ms_levels,true);

% save("data/Ms_dose_fit.mat","Ms_dose_fit");
% save("data/mu_Ms.mat","mu_Ms");
% save("data/dose_Ms_fit.mat","dose_Ms_fit");
% save("data/mu_dose.mat","mu_dose");

%%
% Ms_fit = (116:0.01:119)*1e3;
% doseFromFit = polyval(Ms_dose_fit,Ms_fit,[],mu_Ms);
% 
% idx = 7;
% figure(1)
% plot(Ms_fit,doseFromFit)
% hold on
% plot(Ms_levels(1:idx),dose_levels(1:idx),'.')
% hold off