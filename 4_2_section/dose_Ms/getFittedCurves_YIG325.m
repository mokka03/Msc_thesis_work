function [dMsp_dose_fit,mu_dMsp,dose_dMsp_fit,mu_dose] = getFittedCurves_YIG325(dose_levels,dMs_p,showPlot)
%GETFITTEDCURVES_YIG325 Summary of this function goes here
%   Detailed explanation goes here
    
    %%% Ion dose dMs_p curve
    [dose_dMsp_fit,~,mu_dose] = polyfit(dose_levels,dMs_p,4);
    dose_fit = min(dose_levels)-1e12:1e10:max(dose_levels)*2;
    dMspFromFit = polyval(dose_dMsp_fit,dose_fit,[],mu_dose);

    %%% Ms ion dose curve
    [~, idx] = max(dMspFromFit);
    dMsp_fit = dMspFromFit(idx:end);
    
    [dMsp_dose_fit,~,mu_dMsp] = polyfit(dMsp_fit,dose_fit(idx:end),12);
    doseFromFit = polyval(dMsp_dose_fit,dMsp_fit,[],mu_dMsp);
    
    if showPlot
        figure(21)
        plot(dose_levels/1e12,dMs_p,'.')
        hold on;
        % plot(doseFromFit,dMsFromFit','o');
        plot(dose_fit/1e12,dMspFromFit);
        hold off;
        ylabel('\Delta Ms (%)');
        xlabel('Ion Dose (10^{12} ions/cm^2)');
        xlim([min(dose_fit) max(dose_fit)]/1e12);
        legend({'Measured data', 'Polynomial fit'},'Location','northeast')
        set(gca,'FontSize',15);

        [~, idx] = max(dMs_p);
        figure(22)
        plot(dMs_p(idx:end),dose_levels(idx:end)/1e12,'.')
        hold on;
        plot(dMsp_fit,doseFromFit/1e12);  % expfit
%         plot(Ms_fit/1e3,dose_fit(1:size(Ms_fit,2))/1e12);
        hold off;
        xlabel('Ms (kA/m)');
        ylabel('Ion Dose (10^{12} ions/cm^2)');
%         ylim([0 6]);
        legend({'Measured data', 'Polynomial fit'},'Location','northeast')
        set(gca,'FontSize',15);
        max(dMspFromFit)
        

    end
end

