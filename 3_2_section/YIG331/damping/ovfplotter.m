%% Get data
clear;

path = 'mumax\';
name = 'a11e-4.out';

% my
file_path = strcat(path,name);
data = oommf2matlab(fullfile(file_path,sprintf('m%6.6i.ovf',2)));
absval = data.datax.^2+data.datay.^2+data.dataz.^2;
data.datax(absval ==0.0) = NaN;
data.datay(absval ==0.0) = NaN;
data.dataz(absval ==0.0) = NaN;
my = data.datay';

%% Parameters
nx = size(my,2);
ny = size(my,1);
x_limit = [0 nx+1];
y_limit = [0 ny+1];
x_ticks = 0:200:nx;
y_ticks = 0:200:ny;
x_ticklabels = {'0','20', '40', '60','80','100'};
y_ticklabels = {'0','20', '40', '60','80','100'};


%% my
f2 = plotGeom(12,my,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
c = colorbar;
clim([-0.04 0.04]/2);
c.Label.String = "my";

% SaveFig('figure_v9/','my_mx3_o_large', gcf);
% save('m_100umx100um_binary.mat','m')

%% Slice
my_slice = sum(my(400:600,120:end-120));
my_slice = my_slice/max(my_slice);



[~,locs] = findpeaks(my_slice,'MinPeakHeight',0);
my_slice = my_slice(locs(1):end);
clear locs;

[pks,locs] = findpeaks(my_slice,'MinPeakHeight',0);
locs = locs*1e-7;

x = (0:size(my_slice,2)-1)*1e-7;
f_sim = fit(locs',pks','exp1')

figure(645)
plot(x*1e6,my_slice)
hold on
plot(locs*1e6,pks,'.')
plot(x*1e6,f_sim(x))
yline(0)
hold off

%% Measurement data
load('data/my_meas.mat')
load('data/f_meas.mat')

x_meas = (0:size(my_slice_meas,2)-1)*1e-7;
% [pks_meas,locs_meas] = findpeaks(my_slice_meas,'MinPeakHeight',0);
% locs_meas = locs_meas*1e-7;
% f_meas = fit(locs_meas',pks_meas','exp1')

figure(646)
plot(x,my_slice)
hold on
% plot(locs,pks,'.')
plot(x,f_sim(x))
plot(x_meas,my_slice_meas)
plot(x_meas,f_meas(x_meas)+0.2)
yline(0)
hold off
xlabel('x (\mum)')
ylabel('Normalized amplitude')
legend('Simulated wave', 'Fit on simulation', 'Measured wave', 'Fit on measurement')
set(gca,'fontsize', 15)

% SaveFig('figure/Oneline2100and2150MHz_12dbm_H92/','a11e-4', gcf);
