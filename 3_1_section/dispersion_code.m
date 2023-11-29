%% Magnonic Crystals: From Simple Models toward Applications Jaros?awW. K?os and Maciej Krawczyk (pg. 288)

% set magnetization directin here:
oop = 1;    % (1: out-of-plane, 0: in-plane)


%% constants
mu0 = 4*pi*1e-7;
if oop
    theta = 0;  % H0 angle, 0 if oop
else
    theta = pi/2;
end

d = 20e-9;      % film thickness
if oop
    Bext = 300e-3;   % applied field in T
else
    Bext = 30e-3;   % applied field in T
end
%% material parameters
Ms_CoFe = 1.4e6;
A_CoFe = 30E-12;
Ms_YIG = 1.4e5;
A_YIG = 3.65E-12;
Ms_Py = 8.6e5;
A_Py = 13E-12;

Ms = Ms_YIG;
A = A_YIG;
gamma = 2*pi*28e9; % Hz/T gamma*mu0=2.21e5 used in OOMMF

%% 3D plot and contour
if oop
    H0 = Bext/mu0-Ms; %% effective field
else
    H0 = Bext/mu0; %% effective field
end
kx = 2*pi*[-4e6:1e4:4e6];
ky = 2*pi*[-4e6:1e4:4e6]';
KX = repmat(kx,length(ky),1);
KY = repmat(ky,1,length(kx));
KK = abs(KX+1i*KY);
Phi = angle(KX+1i*KY);
omegaM = gamma*mu0*Ms;
omegaH = gamma*mu0*(H0+2*A/(mu0*Ms).*KK.^2);
P = 1-(1-exp(-KK*d))./(KK*d);
omega = sqrt(omegaH.*(omegaH+omegaM*(P+sin(theta)^2*(1-P.*(1+cos(Phi).^2)+omegaM./omegaH.*(P.*(1-P).*sin(Phi).^2)))));

figure(7)
surfc(kx/2/pi*1e-6,ky/2/pi*1e-6,(omega)*1e-9/2/pi,'EdgeColor','none')
xlabel('k_x (2\pi/\mum)')
ylabel('k_y (2\pi/\mum)')
zlabel('f (GHz)')



%% 1D plot
for i = 0:5
    Bext = Ms*mu0+i*20e-3; %% applied field in T
    H0 = Bext/mu0-Ms; %% effective field
    kxx = 2*pi*[0:1e3:10e6];
    omegaHx = gamma*mu0*(H0+2*A/(mu0*Ms).*kxx.^2);
    Px = 1-(1-exp(-abs(kxx)*d))./(abs(kxx)*d);
    Phi2 = pi/2*0;
    omegax = sqrt(omegaHx.*(omegaHx+omegaM*(Px+sin(theta)^2*(1-Px.*(1+cos(Phi2).^2)+omegaM./omegaHx.*(Px.*(1-Px).*sin(Phi2).^2)))));

    figure(8)
    plot(kxx/2/pi*1e-6,omegax*1e-9/2/pi)
    hold on
    xlabel('k_x (2\pi/\mum)')
    ylabel('f (GHz)')
end
hold off
legend('0 mT','20 mT','40 mT','60 mT','80 mT','100 mT','location','northwest')

set(gca,'FontSize',14)


