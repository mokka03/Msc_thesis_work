clear; %close all;

basedir = ...
'python\models\demux_v5_B215mT\';

load([basedir 'Msat.mat'])

nx = size(Msat,2);
ny = size(Msat,1);
x_ticks = 0:200:nx;
y_ticks = 0:200:ny;
x_ticklabels = {'0','20', '40', '60', '80'};
y_ticklabels = {'0','20', '40', '60', '80', '100'};

probe1 = nsidedpoly(1000, 'Center', [nx-240 145+4*30], 'Radius', 15);
probe2 = nsidedpoly(1000, 'Center', [nx-240 145+11*30], 'Radius', 15);
rectangle_pos = [60,120,401,501];

%% Msat
f1 = plotGeom(1,Msat*1e-3,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
hold on;
rectangle('Position',rectangle_pos,'EdgeColor','r')
rectangle('Position',[60,100,641,21],'EdgeColor','k')
rectangle('Position',[60,620,641,21],'EdgeColor','k')
plot(probe1,'EdgeColor','r','FaceColor','none','LineWidth',1);
plot(probe2,'EdgeColor','r','FaceColor','none','LineWidth',1);
hold off;
colormap summer;
c = colorbar;
clim([114 117.39]);
c.Label.String = "Saturation magnetization (kA/m)";

% SaveFig('figure/demux_v5_B215mT/spintorch/','Msat', gcf);

%% m
for i = 0:1
    load([basedir 'my_input' num2str(i) '.mat'])
    % m = m./Msat;
    
    f2 = plotGeom(i+2,my,x_ticks,y_ticks,x_ticklabels,y_ticklabels);
    hold on;
    rectangle('Position',rectangle_pos,'EdgeColor','r')
    rectangle('Position',[60,100,641,21],'EdgeColor','k')
    rectangle('Position',[60,620,641,21],'EdgeColor','k')
    plot(probe1,'EdgeColor','r','FaceColor','none','LineWidth',1);
    plot(probe2,'EdgeColor','r','FaceColor','none','LineWidth',1);
    hold off;
    c = colorbar;
    c.Label.String = "my";
    clim([-0.03 0.03]*0.25);
    
%     SaveFig('figure/demux_v5_B215mT/spintorch/','my_2150MHz', gcf);
end

%% Source
% load([basedir 'src.mat'])
% timesteps = 8000;
% dt = 20e-12;
% t = 0:dt:timesteps*dt;
% 
% figure(3)
% plot(t*1e9,src*1e3)
% xlim([t(1) t(end)]*1e9)
% xlabel('Time (ns)')
% ylabel('Excitation field (mT)')
% set(gca,'FontSize',14)
% % SaveFig('figure/v1/o6/doubleFreq/','src', gcf);
