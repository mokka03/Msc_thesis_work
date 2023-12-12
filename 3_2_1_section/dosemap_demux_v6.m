%% Using YIG331 sample
clear;

%% Fit curve to dMs-dose data
load('GetFittedCurves\data\Ms_dose_fit.mat')
load('GetFittedCurves\data\mu_Ms.mat')
load('GetFittedCurves\data\dose_Ms_fit.mat')
load('GetFittedCurves\data\mu_dose.mat')


%% Load trained Msat
load('python\models\demux_v6\Msat.mat')
trainedMsat = Msat(121:620,61:460);

% Plot parameters
nx = size(trainedMsat,2);
ny = size(trainedMsat,1);
x_ticks = linspace(0,500,6);
y_ticks = linspace(0,500,6);
x_ticklabels = {'0','10', '20', '30', '40', '50'};
y_ticklabels = {'0','10', '20', '30', '40', '50'};

% Plot
f1 = plotGeom(23,trainedMsat/1e3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap summer;
c = colorbar;
c.Label.String = "Saturation magnetization (kA/m)";
% clim([114.5 116.72]*1e3);
% SaveFig('figure/demux_v6/dosemap/','continuous_Msat', gcf);

%% Claculate dosemap
dosemap = polyval(Ms_dose_fit,trainedMsat,[],mu_Ms);

% Plot
f2 = plotGeom(24,dosemap/1e12,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap autumn;
c = colorbar;
c.Label.String = 'Ion Dose (10^{12} ions/cm^2)';
% SaveFig('figure/demux_v6/dosemap/','continuous_dose_distribution', gcf);

%% Discretise dosemap
discrete_dosemap = zeros(size(dosemap));
dose_max = max(max(dosemap));

discrete_dosemap(dosemap > dose_max/2) = dose_max;

% Plot
f3 = plotGeom(25,discrete_dosemap/1e12,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap autumn;
c = colorbar;
c.Label.String = 'Ion Dose (10^{12} ions/cm^2)';
% SaveFig('figure/demux_v6/dosemap/','discrete_dose_distribution', gcf);

%% Just for plot dose distribution
discrete_dosemap_full = zeros(size(Msat));
discrete_dosemap_full(121:620,61:460) = discrete_dosemap;
discrete_dosemap_full(100:120,61:end) = 20e12;
discrete_dosemap_full(621:640,61:end) = 20e12;


x_ticks_ = linspace(0,700,8);
y_ticks_ = 0:100:700;
x_ticklabels_ = {'0','10', '20', '30', '40', '50', '60', '70'};
y_ticklabels_ = {'0','10', '20', '30', '40', '50', '60', '70'};
f30 = plotGeom(250,discrete_dosemap_full/1e12,x_ticks_,y_ticks_,x_ticklabels_,y_ticklabels_);
colormap autumn;
c = colorbar;
clim([0 8])
c.Label.String = 'Ion Dose (10^{12} ions/cm^2)';
% SaveFig('figure/demux_v6/dosemap/','discrete_dose_distribution_full', gcf);

%% Calculate discrete Msat

discrete_trainedMsat = polyval(dose_Ms_fit,discrete_dosemap,[],mu_dose);

% Plot
f4 = plotGeom(26,discrete_trainedMsat/1e3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap summer;
c = colorbar;
c.Label.String = "Saturation magnetization (kA/m)";
% clim([114 116.72]*1e3);
% SaveFig('figure/demux_v6/dosemap/','discrete_Msat', gcf);



%% Thresholding for mumax
Msat_binary_thresholding = zeros(size(Msat));
Msat_binary_thresholding(121:620,61:460) = discrete_trainedMsat;
Msat_binary_thresholding(Msat_binary_thresholding < 117e3) = 1;
Msat_binary_thresholding(Msat_binary_thresholding > 117e3) = 0;

figure(30)
imshow(Msat_binary_thresholding);


% imwrite(Msat_binary_thresholding, 'figure/demux_v6/mumax/binary_pictures_thresholding/Msat_level_0.png');

%% Binary pictures for FIB thresholding
distr_thresholding = 1 - Msat_binary_thresholding(121:620,61:460);
FIB_thresholding = zeros(540,640);
FIB_thresholding(21:520,1:400) = distr_thresholding;

figure(31)
imshow(FIB_thresholding);

% imwrite(FIB_thresholding, 'figure/demux_v6/files_for_FIB/distribution_thresholding.png');
% write_xbm(FIB_thresholding,"figure/demux_v6/files_for_FIB/distribution_thresholding.xbm");
% save('figure/demux_v6/files_for_FIB/distribution_thresholding.mat','FIB_thresholding')

%% Walls for FIB

distr_wall = zeros(540,640);
distr_wall(1:20,:) = 1;
distr_wall(521:end,:) = 1;

figure(32)
imshow(distr_wall);
%%
% wall_Ms = polyval(dose_Ms_fit,20e12,[],mu_dose)

% imwrite(distr_wall, 'figure/demux_v43_lr1/files_for_FIB/wall.png');
% write_xbm(distr_wall,"figure/demux_v43_lr1/files_for_FIB/wall.xbm");
% save('figure/demux_v43_lr1/files_for_FIB/wall.mat','distr_wall')
