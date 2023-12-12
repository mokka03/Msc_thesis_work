clear;
load('data/my_finescan40_198p8mT.mat')
meas_output = my(180:680,355:525)+1e-6;
meas_slice = sum(abs(my(180:680,425:455)+1e-6),2);
meas_slice = meas_slice/max(meas_slice);

my_path = ...
'models\MagnSim_400\my.mat';
load(my_path)
sim_output = my(120:620,385:555);
sim_slice = sum(abs(my(120:620,445:475)),2);
sim_slice = sim_slice/max(sim_slice);

x = (1:501)/10;
y = (1:171)/10;
probe = nsidedpoly(1000, 'Center', [x(end)/2 y(end)/2], 'Radius', 1.45);
rectangle_pos = [0,7,50,3];

figure(43)
subplot(3,1,1)
pcolor(x,y,sim_output'); shading interp;
hold on;
plot(probe,'EdgeColor','r','FaceColor','none','LineWidth',1);
rectangle('Position',rectangle_pos,'EdgeColor','k','LineWidth',0.8)
hold off;
clim([-1 1]*2e-2);
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
xlabel('y (\mum)')
ylabel('x (\mum)')
title('Simulation')
set(gca,'FontSize',10)

subplot(3,1,2)
pcolor(x,y,meas_output'); shading interp;
hold on;
plot(probe,'EdgeColor','r','FaceColor','none','LineWidth',1);
rectangle('Position',rectangle_pos,'EdgeColor','k','LineWidth',0.8)
hold off;
clim([-1 1]*8e-6)
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
xlabel('y (\mum)')
ylabel('x (\mum)')
title('Measurement')
set(gca,'FontSize',10)

subplot(3,1,3)
plot(x,sim_slice)
hold on;
plot(x,meas_slice)
hold off;
legend('Simulation', 'Measurement')
xlim([x(1) x(end)]);
xlabel('y (\mum)')
ylabel('I_{norm}')
title('Normalized spin-wave intensity')
set(gca,'FontSize',10)

% SaveFig('figure/focus_YIG325_v12_400/198p8mT/','output', gcf);