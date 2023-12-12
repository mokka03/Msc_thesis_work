clear;
load("data\fine_scan_noplot.mat")

my = rot90(scan_real(:,1:y_size),3);
my = imresize(my,[900 1250]);
my = fliplr(my(:,1:860));

x = (1:860)/10;
y = (1:900)/10;

% Filter
my = imgaussfilt(my,5);

figure(1)
pcolor(x,y,my); axis equal; shading interp;
hold on
rectangle('Position',[7,8,50,30],'EdgeColor','r', 'LineWidth',1)
rectangle('Position',[7,40,50,30],'EdgeColor','r', 'LineWidth',1)
rectangle('Position',[7,6,100,2],'EdgeColor','k')
rectangle('Position',[7,38,100,2],'EdgeColor','k')
rectangle('Position',[7,70,100,2],'EdgeColor','k')
hold off
c = colorbar;
clim([-1.5 1.5]*3e-6)
c.Label.String = "my (a.u.)";
xlabel('x (\mum)')
ylabel('y (\mum)')
xlim([x(1) x(end)])
set(gca,'FontSize',15);

% SaveFig('figure/','fine_scan_noplot', gcf);
% save('data/my_fine_scan_noplot.mat','my');