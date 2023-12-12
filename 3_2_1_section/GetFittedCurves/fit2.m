clear;
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\YIG_samples\YIG328\data\dose_Ms_data_2050MHz.mat')
dose_levels = dose_Ms_data(1,1:end);
Ms_levels = dose_Ms_data(2,1:end);
%%
% Ion dose dMs curve
[dose_Ms_fit,~,mu_dose] = polyfit(dose_levels,Ms_levels,2);
dose_fit = min(dose_levels)-1e12:1e10:max(dose_levels)*2;
MsFromFit = polyval(dose_Ms_fit,dose_fit,[],mu_dose);

figure(21)
plot(dose_levels/1e12,Ms_levels/1e3,'.')
hold on;
% plot(doseFromFit,dMsFromFit','o');
plot(dose_fit/1e12,MsFromFit/1e3);
hold off;
ylabel('Ms (kA/m)');
xlabel('Ion Dose (10^{12} ions/cm^2)');
% xlim([min(dose_levels) max(dose_levels)]/1e12);
legend({'Measured data', 'Polynomial fit'},'Location','northwest')
set(gca,'FontSize',15);
%%
wall_Ms = polyval(dose_Ms_fit,20e12,[],mu_dose)