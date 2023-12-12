%% Get data
clear;

path = 'mumax\';
name = 'minidosemap_final.out';
d = 2;

% Msat
file_path = strcat(path,name);
data = oommf2matlab(fullfile(file_path,sprintf('Msat%6.6i.ovf',0)));
absval = data.datax.^2+data.datay.^2+data.dataz.^2;
data.datax(absval ==0.0) = NaN;
data.datay(absval ==0.0) = NaN;
data.dataz(absval ==0.0) = NaN;
Msat = data.datax';

% my
file_path = strcat(path,name);
data = oommf2matlab(fullfile(file_path,sprintf('m%6.6i.ovf',2)));
absval = data.datax.^2+data.datay.^2+data.dataz.^2;
data.datax(absval ==0.0) = NaN;
data.datay(absval ==0.0) = NaN;
data.dataz(absval ==0.0) = NaN;
my = data.datay';

%% Parameters
nx = size(Msat,2);
ny = size(Msat,1);
x_limit = [0 nx+1];
y_limit = [0 ny+1];
x_ticks = 0:400:nx;
y_ticks = 0:400:ny;
x_ticklabels = {'0','20', '40', '60', '80', '100', '120'};
y_ticklabels = {'0','20', '40', '60', '80', '100'};

%% Msat
f1 = plotGeom(11,Msat*1e-3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
c = colorbar;
colormap summer;
c.Label.String = "Msat [kA/m]";

% SaveFig('figure/','Msat_20', gcf);
% save('Msat_mx3.mat','Msat')

%% my
% my = my(1:200,:);
f2 = plotGeom(12,my,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
hold on
rectangle_pos = [200,440,1001,601];
rectangle('Position',rectangle_pos,'EdgeColor','r')
rectangle_pos = [200,1080,1001,601];
rectangle('Position',rectangle_pos,'EdgeColor','r')
rectangle_pos = [200,400,2001,41];
rectangle('Position',rectangle_pos,'EdgeColor','k')
rectangle_pos = [200,1040,2001,41];
rectangle('Position',rectangle_pos,'EdgeColor','k')
rectangle_pos = [200,1680,2001,41];
rectangle('Position',rectangle_pos,'EdgeColor','k')
hold off;
c = colorbar;
clim([-0.04 0.04]*0.08);
c.Label.String = "my";

% SaveFig('figure/','my_final', gcf);
% save('m_100umx100um_binary.mat','m')

%%
my1 = my(1100:1600,120:1120);
figure(8)
image(my1,'CDataMapping','scaled')

lambda = WavelengthFourier(my1,5e-8)
