clear;

load('data\my_fine_scan_noplot_1s_2.mat')
my = imresize(my,[1000 1250]);

x = (1:1000)*1e-1;
y = (1:1000)*1e-1;
figure(83)
pcolor(x,y,fliplr(my(:,1:1000))); axis equal; shading interp;
hold on
rectangle('Position',[1,18,50,30],'EdgeColor','r', 'LineWidth',1)
rectangle('Position',[1,50,50,30],'EdgeColor','r', 'LineWidth',1)
rectangle('Position',[1,16,100,2],'EdgeColor','k')
rectangle('Position',[1,48,100,2],'EdgeColor','k')
rectangle('Position',[1,80,100,2],'EdgeColor','k')
hold off
c = colorbar;
clim([-1.5 1.5]*1e-5)
c.Label.String = "my (a.u.)";
xlabel('x (\mum)')
ylabel('y (\mum)')
xlim([x(1) x(end)])
set(gca,'FontSize',15);

% SaveFig('figure\','my_fine_scan_noplot_1s_2',gcf);

%%
intr = my(850:1000,1:1000);
short = my(500:800,500:1000);
long = my(200:450,500:1000);

figure(84)
image(short,'CDataMapping','scaled')
colorbar;
clim([-1.5 1.5]*1e-5)

%% Average
intr = intr(:,25:end);
for_slice = short;

figure(85)
clf;
hold all;
for i = 1:size(for_slice,1)
    plot(for_slice(i,:))
end

slice_intr = sum(intr);
slice_short = sum(short);

figure(86)
plot(slice_short)

%%
dx = 100e-9;
lambda_intr = WavelengthFourier(slice_intr,dx);
lambda_short = WavelengthFourier(slice_short,dx);

% save('data/sliceIntr_fine_scan_noplot_1s_2.mat','slice_intr');
% save('data/sliceShort_fine_scan_noplot_1s_2.mat','slice_short');

