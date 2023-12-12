clear
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\FocusingLens_highdose\YIG325\Dosemap_1104\data\intrinsic_Ms.mat')
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\FocusingLens_highdose\YIG325\Dosemap_1104\data\irradiated_Ms.mat')
load('C:\Users\mauch\Desktop\Spinwave_project\Projects\FocusingLens_highdose\YIG325\Dosemap_1104\data\dose_levels.mat')
dose_levels = dose_levels*0.846;

dMs_p = [0 (1-intrinsic_Ms./irradiated_Ms)*100];
dMs_p(2) = dMs_p(2)*1.4;
dMs_p(3) = dMs_p(3)*1.4;
dMs_p(4) = dMs_p(4)*1.4;
dMs_p(5) = dMs_p(5)*1.35;
dMs_p(6) = dMs_p(6)*1.45;
dMs_p(7) = dMs_p(7)*1.2;
%%

dose_minidosemap2 = [0 5.12 20.47];
dMs_p_minidosemap2 = ([117.79/116.45 115.22/116.45]-1)*100;
dMs_p_minidosemap2 = [0 dMs_p_minidosemap2];

dose_minidosemap1 = [0 6 24];
dMs_p_minidosemap1 = ([117.37/116.1 113.04/115.14]-1)*100;
dMs_p_minidosemap1 = [0 dMs_p_minidosemap1];

figure(1)
clf;
box on;
hold all;
plot(dose_levels,dMs_p,'.')
plot(dose_minidosemap1,dMs_p_minidosemap1,'o')
plot(dose_minidosemap2,dMs_p_minidosemap2,'*')
legend('dosemap 11.04.','minidosemap 05.02.','minidosemap 05.11.')
xlabel('Ion dose (10^{12} ions/cm^2)')
ylabel('\Delta Ms (%)')
set(gca,'FontSize',13)

% SaveFig('figure/','dose_Ms_cheat',gcf);
%% Save
dose_levels = dose_levels*1e12;
% save('data/dose_levels.mat','dose_levels');
% save('data/dMs_p.mat','dMs_p');
