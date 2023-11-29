clear;

%% Fit curve to dMs-dose data
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\YIG_samples\YIG331\myStuff\FIBregions2100MHz_10dbm_H92\data\dose_Ms_data.mat')

%%
dose_levels = dose_Ms_data(1,1:end);
Ms_levels = dose_Ms_data(2,1:end);

[Ms_dose_fit,mu_Ms,dose_Ms_fit,mu_dose] = getFittedCurves(dose_levels,Ms_levels,true);

% save("data/Ms_dose_fit.mat","Ms_dose_fit");
% save("data/mu_Ms.mat","mu_Ms");
% save("data/dose_Ms_fit.mat","dose_Ms_fit");
% save("data/mu_dose.mat","mu_dose");