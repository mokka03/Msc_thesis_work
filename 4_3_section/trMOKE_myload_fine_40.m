clear;
load("data\finescan40_198p8mT.mat")

my = rot90(scan_real(:,1:y_size),3);
my = imresize(my,[800 1250]);

for i = 1:size(my,1)
    my(i,:) = circshift(my(i,:),floor(i/15));
end

xx = 1000;
my = fliplr(my(:,100:xx));

x = (1:xx-99)/10;
y = (1:800)/10;

rectangle_pos = [4,18,40,50];
probe = nsidedpoly(1000, 'Center', [44 43], 'Radius', 1.5);

% Filter 1
my = imgaussfilt(my,5);

figure(1)
pcolor(x(1:700),y(50:790),my(50:790,1:700)); axis equal; shading interp;
hold on;
rectangle('Position',rectangle_pos,'EdgeColor','r','LineWidth',0.8)
plot(probe,'EdgeColor','r','FaceColor','none','LineWidth',1);
hold off;
c = colorbar;
c.Label.String = "my (a.u.)";
clim([-1 1]*5e-6)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca,'FontSize',15);
xlim([x(1) 70]);

% SaveFig('figure/focus_YIG325_v12_400/198p8mT/','finescan40_198p8mT', gcf);
% save('data/my_finescan40_198p8mT.mat','my');