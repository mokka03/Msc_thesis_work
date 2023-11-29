addpath('C:\femm42\mfiles');

openfemm;
newdocument(0);

%geometrical parameters in micrometers
thickness = 0.3;
width_ground = 2.0;
width_signal = 4.0;
gap = 2.0;

%enclosure dimensions in micrometers
encH = 1000;
encV = 1000;

%y-position of the cross section line of flux density
Ypos = -0.270;

%dimensions in mumax
dx = 1e-6;
dy = 1e-6;
dz = 540e-9;

Nx = 500;
Ny = 400;

Xsize = dx*Nx;
Ysize = dy*Ny;

%frequency
f = 3e9;

%define the striplines
p_left = [0,0
          width_ground,0
          width_ground,thickness
          0,thickness];

p_middle = [width_ground + gap,0
            width_ground + width_signal + gap,0
            width_ground + width_signal + gap,thickness
            width_ground + gap,thickness];

p_right = [width_ground + width_signal + 2.0*gap,0
            2*width_ground + width_signal + 2.0*gap,0
            2*width_ground + width_signal + 2.0*gap,thickness
            width_ground + width_signal + 2.0*gap,thickness];

DrawClosedPolygon(p_left);
DrawClosedPolygon(p_middle);
DrawClosedPolygon(p_right);

%add materials
mi_getmaterial('Copper')
mi_modifymaterial( 'Copper', 5, 41); % Gold conductivity
mi_getmaterial('Air');
    
%define circuits
mi_addcircprop('Ground',5.0e-4,1);
mi_addcircprop('Source',1.0e-3,1);

%add material labels
mi_addblocklabel(width_ground/2.0,thickness/2);
mi_selectlabel(width_ground/2.0,thickness/2);
mi_setblockprop('Copper', 0, 0, 'Ground', 0, 0, -1);
mi_clearselected();

mi_addblocklabel(width_ground + gap + 0.5*width_signal,thickness/2);
mi_selectlabel(width_ground + gap + 0.5*width_signal,thickness/2);
mi_setblockprop('Copper', 0, 0, 'Source', 0, 0, 1);
mi_clearselected();

mi_addblocklabel(2.0*gap + 1.5*width_ground + width_signal,thickness/2);
mi_selectlabel(2.0*gap + 1.5*width_ground + width_signal,thickness/2);
mi_setblockprop('Copper', 0, 0, 'Ground', 0, 0, -1);
mi_clearselected();


mi_addblocklabel(gap + 0.5*width_signal + width_ground,6.0*thickness);
mi_selectlabel(gap + 0.5*width_signal + width_ground,6.0*thickness);
mi_setblockprop('Air', 0, 0, 'None', 0, 0, 0);
mi_clearselected();

%add boundary condition
mi_addboundprop('Zero_A',0,0,0,0,0,0,0,0,0,0,0);

%draw enclosure
enc = [-encH/2.0, -encV/2.0
       encH/2.0, -encV/2.0
       encH/2.0, encV/2.0
       -encH/2.0,encV/2.0];
DrawClosedPolygon(enc);

enc = [enc;enc(1,:)];
for i = 1:4
    enc_mid = (enc(i,:) + enc(i+1,:))/2;
    mi_selectsegment(enc_mid(1),enc_mid(2));
end

mi_setsegmentprop('Zero_A',0,0,0,0);
mi_clearselected();

mi_probdef(f,'micrometers','planar',1e-8,120,30,0);

mi_saveas(char(strrep(cd + "\CPW_design_03.fem", '\', '\\')));

mi_analyze(0);
mi_loadsolution();

%getting the cross section of the flux density along a line
%columns of the resulting array:
% x-coordinate, y coordinate, Re(Bx), Re(By), Im(Bx), Im(By), abs(Bx),
% abs(By)
flux_density = [];
for j = 1:Nx
    flux_density(j, 1) = j; %x-coordinate
    flux_density(j, 2) = Ypos; %y-coordinate
    tmp = mo_getb(-74.5 + j, Ypos);
    flux_density(j, 3) = real(tmp(1,1)); %Re(Bx)
    flux_density(j, 4) = real(tmp(1,2)); %Re(By)
    flux_density(j, 5) = imag(tmp(1,1)); %Im(Bx)
    flux_density(j, 6) = imag(tmp(1,2)); %Im(By)
%    flux_density(j, 7) = sqrt((real(tmp(1,1)))^2+(imag(tmp(1,1)))^2);
%    flux_density(j, 8) = sqrt((real(tmp(1,2)))^2+(imag(tmp(1,2)))^2);
end

%generate the mumax input file for the real part of Bx and Bz
%y-->z in mumax coordinate system

fname = cat(2,'CPW_DL_Real_3GHz','.ohf');
ohffile = fopen(fname,'w');

fprintf(ohffile,'# OOMMF: rectangular mesh v1.0\n');
fprintf(ohffile,'# Segment count: 1\n');
fprintf(ohffile,'# Begin: Segment\n');
fprintf(ohffile,'# Begin: Header\n');
fprintf(ohffile,'# Title: %s\n',fname);
fprintf(ohffile,'# meshtype: rectangular\n');
fprintf(ohffile,'# meshunit: m\n');
fprintf(ohffile,'# xbase: %E\n',dx/2.0);
fprintf(ohffile,'# ybase: %E\n',dz/2.0);
fprintf(ohffile,'# zbase: %E\n',dy/2.0);
fprintf(ohffile,'# xstepsize: %E\n',dx);
fprintf(ohffile,'# ystepsize: %E\n',dy);
fprintf(ohffile,'# zstepsize: %E\n',dz);
fprintf(ohffile,'# xnodes: %i\n',Nx);
fprintf(ohffile,'# ynodes: %i\n',Ny);
fprintf(ohffile,'# znodes: 1\n');
fprintf(ohffile,'# xmin: 0\n');
fprintf(ohffile,'# ymin: 0\n');
fprintf(ohffile,'# zmin: 0\n');
fprintf(ohffile,'# xmax: %E\n',Xsize);
fprintf(ohffile,'# ymax: %E\n',Ysize);
fprintf(ohffile,'# zmax: %E\n',dz);
fprintf(ohffile,'# valueunit: T\n');
fprintf(ohffile,'# valuemultiplier: 1.0\n');
fprintf(ohffile,'# ValueRangeMinMag: %E\n',min((min([flux_density(:,3), flux_density(:,4)]))));
fprintf(ohffile,'# ValueRangeMaxMag: %E\n',max((max([flux_density(:,3), flux_density(:,4)]))));
fprintf(ohffile,'# End: Header\n');
fprintf(ohffile,'# Begin: Data Text\n');

for ii = 1 : Ny
    
    for jj = 1 : Nx
        
        
        fprintf(ohffile,'%E %E %E\n', [flux_density(jj,3) 0.0 flux_density(jj,4)]);
        
    end
    
end




fprintf(ohffile,'# End: Data Text\n');
fprintf(ohffile,'# End: Segment\n');

fclose all;

%generate the mumax input file for the imaginary part of Bx and Bz
%y-->z in mumax coordinate system

fname = cat(2,'CPW_DL_Imag_3GHz','.ohf');
ohffile = fopen(fname,'w');

fprintf(ohffile,'# OOMMF: rectangular mesh v1.0\n');
fprintf(ohffile,'# Segment count: 1\n');
fprintf(ohffile,'# Begin: Segment\n');
fprintf(ohffile,'# Begin: Header\n');
fprintf(ohffile,'# Title: %s\n',fname);
fprintf(ohffile,'# meshtype: rectangular\n');
fprintf(ohffile,'# meshunit: m\n');
fprintf(ohffile,'# xbase: %E\n',dx/2.0);
fprintf(ohffile,'# ybase: %E\n',dz/2.0);
fprintf(ohffile,'# zbase: %E\n',dy/2.0);
fprintf(ohffile,'# xstepsize: %E\n',dx);
fprintf(ohffile,'# ystepsize: %E\n',dy);
fprintf(ohffile,'# zstepsize: %E\n',dz);
fprintf(ohffile,'# xnodes: %i\n',Nx);
fprintf(ohffile,'# ynodes: %i\n',Ny);
fprintf(ohffile,'# znodes: 1\n');
fprintf(ohffile,'# xmin: 0\n');
fprintf(ohffile,'# ymin: 0\n');
fprintf(ohffile,'# zmin: 0\n');
fprintf(ohffile,'# xmax: %E\n',Xsize);
fprintf(ohffile,'# ymax: %E\n',Ysize);
fprintf(ohffile,'# zmax: %E\n',dz);
fprintf(ohffile,'# valueunit: T\n');
fprintf(ohffile,'# valuemultiplier: 1.0\n');
fprintf(ohffile,'# ValueRangeMinMag: %E\n',min((min([flux_density(:, 5), flux_density(:, 6)]))));
fprintf(ohffile,'# ValueRangeMaxMag: %E\n',max((max([flux_density(:, 5), flux_density(:, 6)]))));
fprintf(ohffile,'# End: Header\n');
fprintf(ohffile,'# Begin: Data Text\n');

for ii = 1 : Ny
    
    for jj = 1 : Nx
        
        
        fprintf(ohffile,'%E %E %E\n', [flux_density(jj, 5) 0.0 flux_density(jj, 6)]);
        
    end
    
end




fprintf(ohffile,'# End: Data Text\n');
fprintf(ohffile,'# End: Segment\n');

fclose all;

%this is the function to draw the polygons
function DrawClosedPolygon(p)
for i = 1:size(p,1)
    mi_addnode(p(i,1),p(i,2));
end

for i = 1:size(p,1)-1
    mi_addsegment(p(i,1),p(i,2),p(i+1,1),p(i+1,2));
end
mi_addsegment(p(size(p,1),1),p(size(p,1),2),p(1,1),p(1,2));
end