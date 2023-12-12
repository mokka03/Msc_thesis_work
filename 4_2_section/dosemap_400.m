clear;

%% Fit curve to dMs-dose data
load('dose_Ms\data\dMsp_dose_fit.mat')
load('dose_Ms\data\mu_dMsp.mat')
load('dose_Ms\data\dose_dMsp_fit.mat')
load('dose_Ms\data\mu_dose.mat')

%% Load trained Msat
load('models\focus_YIG325_v12_400\Msat.mat')
trainedMsat = Msat(91:650,61:end);

% Plot parameters
nx = size(trainedMsat,2);
ny = size(trainedMsat,1);
x_ticks = 0:200:nx;
y_ticks = 0:200:ny;
x_ticklabels = {'0','20', '40', '60', '80', '100'};
y_ticklabels = {'0','20', '40', '60', '80', '100'};

% Plot
f1 = plotGeom(23,trainedMsat,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap summer;
c = colorbar;
c.Label.String = "Saturation magnetization [kA/m]";
% clim([114.5 116.72]*1e3);
% SaveFig('figure/focus_YIG325_v12_400/dosemap/','trainedMsat', gcf);

%% Claculate dosemap
Ms0 = 116.45e3;
Msat0 = ones(size(trainedMsat))*Ms0;

dMsat_percent = (1-Msat0./trainedMsat)*100;

dosemap = polyval(dMsp_dose_fit,dMsat_percent,[],mu_dMsp);
dosemap(31:530,401:end) = 0;

% Plot
f2 = plotGeom(24,dosemap,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap autumn;
c = colorbar;
c.Label.String = 'Ion Dose [10^{12} ions/cm^2]';

% SaveFig('figure/focus_YIG325_v12_400/dosemap/','continuous_dosemap', gcf);

%% RGB images for FIB
grayscaled = uint8(dosemap/max(max(dosemap))*255);
RGB = zeros([size(grayscaled),3]);
RGB(:,:,1) = grayscaled;
RGB(:,:,2) = grayscaled;
RGB(:,:,3) = grayscaled;

figure(25)
imshow(grayscaled);

% imwrite(uint8(RGB),'figure\focus_YIG325_v12_400\for_FIB\dosemap.bmp');

%% Discrete dosemap
discrete_dosemap = double(grayscaled)/255*max(max(dosemap));

f3 = plotGeom(26,discrete_dosemap,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap autumn;
c = colorbar;
c.Label.String = 'Ion Dose [10^{12} ions/cm^2]';

% SaveFig('figure/focus_YIG325_v12_400/dosemap/','discrete_dosemap', gcf);
% save('figure/focus_YIG325_v12_400/dosemap/discrete_dosemap.mat','discrete_dosemap');

%% Discrete Msat

discrete_Msat = polyval(dose_dMsp_fit,discrete_dosemap,[],mu_dose);
discrete_Msat = (1 + discrete_Msat/100).*Msat0;

f3 = plotGeom(27,discrete_Msat,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap summer;
c = colorbar;
c.Label.String = "M_s (kA/m)";
% clim([114.5 116.72]*1e3);

% SaveFig('figure/focus_YIG325_v12_400/dosemap/','discrete_Msat', gcf);

%%
Msat(91:650,61:end) = discrete_Msat;

f3 = plotGeom(28,Msat*1e-3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
colormap summer;
c = colorbar;
c.Label.String = "M_s (kA/m)";
% clim([114.5 116.72]*1e3);

% save('figure/focus_YIG325_v12_400/dosemap/Msat_trained.mat','Msat');