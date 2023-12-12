%% Get data
clear;

path = 'mumax\';
% name = 'demux_v6_2050MHz.out';
% name = 'demux_v6_2100MHz.out';
name = 'demux_v6_doubleFreq.out';

% Msat
file_path = strcat(path,name);
data = oommf2matlab(fullfile(file_path,sprintf('Msat%6.6i.ovf',0)));
absval = data.datax.^2+data.datay.^2+data.dataz.^2;
data.datax(absval ==0.0) = NaN;
data.datay(absval ==0.0) = NaN;
data.dataz(absval ==0.0) = NaN;
Msat = flipud(data.datax');

% my
file_path = strcat(path,name);
data = oommf2matlab(fullfile(file_path,sprintf('m%6.6i.ovf',2)));
absval = data.datax.^2+data.datay.^2+data.dataz.^2;
data.datax(absval ==0.0) = NaN;
data.datay(absval ==0.0) = NaN;
data.dataz(absval ==0.0) = NaN;
my = flipud(data.datay');

%% Parameters
nx = size(Msat,2);
ny = size(Msat,1);
x_limit = [0 nx+1];
y_limit = [0 ny+1];
x_ticks = 0:400:nx;
y_ticks = 0:400:ny;
x_ticklabels = {'0','20', '40', '60', '80'};
y_ticklabels = {'0','20', '40', '60', '80'};

probe1 = nsidedpoly(1000, 'Center', [nx-240*2 290+4*60], 'Radius', 30);
probe2 = nsidedpoly(1000, 'Center', [nx-240*2 290+11*60], 'Radius', 30);
rectangle_pos = [60*2,120*2,801,1001];

%% Msat
f1 = plotGeom(11,Msat*1e-3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
hold on
rectangle('Position',rectangle_pos,'EdgeColor','r')
rectangle('Position',[120,200,1281,41],'EdgeColor','k')
rectangle('Position',[120,1240,1281,41],'EdgeColor','k')
plot(probe1,'EdgeColor','r','FaceColor','none','LineWidth',1);
plot(probe2,'EdgeColor','r','FaceColor','none','LineWidth',1);
hold off;
c = colorbar;
clim([114 117.4])
% clim([0 1]*1e-6)
colormap summer;
c.Label.String = "Msat (kA/m)";

% SaveFig('figure/demux_v6/mumax/thresholding/','Msat', gcf);
% save('figure/demux_v43_lr1/mumax/Msat_mx3.mat','Msat')

%% my
f2 = plotGeom(12,flipud(my),x_ticks,y_ticks,x_ticklabels,y_ticklabels);
hold on
rectangle('Position',rectangle_pos,'EdgeColor','r')
rectangle('Position',[120,200,1281,41],'EdgeColor','k')
rectangle('Position',[120,1240,1281,41],'EdgeColor','k')
plot(probe1,'EdgeColor','r','FaceColor','none','LineWidth',1);
plot(probe2,'EdgeColor','r','FaceColor','none','LineWidth',1);
hold off;
c = colorbar;
% clim([-0.01 0.01]*1.2);
clim([-0.01 0.01]*.6);
c.Label.String = "my";

% SaveFig('figure/demux_v6/mumax/thresholding/','my_2050MHz', gcf);
% SaveFig('figure/demux_v6/mumax/thresholding/','my_doubleFreq', gcf);
% save('m_100umx100um_binary.mat','m')

%%
lambda = WavelengthFourier(my(1:10,120:end),5e-8)
