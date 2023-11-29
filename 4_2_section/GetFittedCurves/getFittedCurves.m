function [Ms_dose_fit,mu_Ms,dose_Ms_fit,mu_dose] = getFittedCurves(dose_levels,Ms_levels,showPlot)
%GETFITTEDCURVES Summary of this function goes here
%   Detailed explanation goes here
    
    % Ion dose dMs curve
    [dose_Ms_fit,~,mu_dose] = polyfit(dose_levels,Ms_levels,4);
    dose_fit = min(dose_levels)-1e12:1e10:max(dose_levels);
    MsFromFit = polyval(dose_Ms_fit,dose_fit,[],mu_dose);

    % Ms ion dose curve
    [~, idx] = max(MsFromFit);
    Ms_fit = MsFromFit(1:idx);

    [Ms_dose_fit,~,mu_Ms] = polyfit(Ms_fit,dose_fit(1:size(Ms_fit,2)),12);
    doseFromFit = polyval(Ms_dose_fit,Ms_fit,[],mu_Ms);
    
    if showPlot
        figure(21)
        plot(dose_levels/1e12,Ms_levels/1e3,'.')
        hold on;
        % plot(doseFromFit,dMsFromFit','o');
        plot(dose_fit/1e12,MsFromFit/1e3);
        hold off;
        ylabel('Ms (kA/m)');
        xlabel('Ion Dose (10^{12} ions/cm^2)');
        xlim([min(dose_levels) max(dose_levels)]/1e12);
        legend({'Measured data', 'Polynomial fit'},'Location','northwest')
        set(gca,'FontSize',15);

        figure(22)
        plot(Ms_levels/1e3,dose_levels/1e12,'.')
        hold on;
        plot(Ms_fit/1e3,doseFromFit/1e12);  % expfit
%         plot(Ms_fit/1e3,dose_fit(1:size(Ms_fit,2))/1e12);
        hold off;
        xlabel('Ms (kA/m)');
        ylabel('Ion Dose (10^{12} ions/cm^2)');
        ylim([0 6]);
        legend({'Measured data', 'Polynomial fit'},'Location','northwest')
        set(gca,'FontSize',15);
        max(MsFromFit)
        

    end
end

