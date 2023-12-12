clear;

%% Create Msat
d = 20;

Msat = ones(1060,1200)*115.65e3;
Ms_ = [116.72 114.25 112]*1e3;
y_length = 500;

Msat(221:520,101:y_length+100) = ones(300,y_length)*Ms_(1);
Msat(521+d:520+d+300,101:y_length+100) = ones(300,y_length)*Ms_(2);

Msat(201:220,101:y_length*2+100)  = ones(20,y_length*2)*Ms_(3);
Msat(521:540,101:y_length*2+100)  = ones(20,y_length*2)*Ms_(3);
Msat(841:860,101:y_length*2+100)  = ones(20,y_length*2)*Ms_(3);


figure(1)
pcolor(Msat*1e-3); axis equal; shading interp;
xlim([1 size(Msat,2)])
ylim([1 size(Msat,1)])
c = colorbar;

% save('data/Msat_minidosemap1.mat','Msat');

%% Binaryfigures for mumax
for i = 1:3
    Msat_binary = Msat == Ms_(i);
    Msat_binary = 1 - Msat_binary;
%     imwrite(Msat_binary,['mumax\binary_pictures\20\Msat_level_' num2str(i) '.png']);
end

%% dosemap
Msat_mod = Msat(201:860,101:1100);
doses = [6 24 34]*1e12;
dosemap = ones(660,1000);
for i = 1:3
    dosemap(Msat_mod == Ms_(i)) = doses(i);
end

figure(2)
pcolor((0.1:0.1:100),(0.1:0.1:66),flipud(dosemap*1e-12)); axis equal; shading interp;
c = colorbar;
colormap summer
c.Label.String = "Ion dose (10^{12} ions/cm^2)";
xlim([1 size(Msat_mod,2)]/10)
ylim([1 size(Msat_mod,1)]/10)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca,'FontSize',14)

% SaveFig('for_FIB/','dosemap', gcf);

%% Binary figures for FIB
for i = 1:3
    binary_distr = dosemap == doses(i);
%     binary_distr = 1 - binary_distr;
%     imwrite(binary_distr,['for_FIB\binary_pictures\dose_level_' num2str(i) '.bmp']);
end

%%
im = imread('for_FIB\binary_pictures\dose_level_1.bmp');

%%
grayscaled = uint8(dosemap/max(max(dosemap))*255);
RGB = zeros([size(grayscaled),3]);
RGB(:,:,1) = grayscaled;
RGB(:,:,2) = grayscaled;
RGB(:,:,3) = grayscaled;

% imwrite(uint8(RGB),'for_FIB\binary_pictures_grayscaled\dosemap.bmp');

%%
sum(sum(grayscaled))/255*54.4*1e-6
