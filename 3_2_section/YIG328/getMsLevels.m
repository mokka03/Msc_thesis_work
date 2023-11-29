function MsLevels = getMsLevels(params)
%GETMSLEVELS Summary of this function goes here
%   Detailed explanation goes here
    % Parameters
    mu0 = 4*pi*1e-7;
    gamma = 2*pi*28e9; % Hz/T gamma*mu0=2.21e5 used in OOMMF
    A_exch = params.A_exch;
    d = params.d;
    theta = params.theta;
    Bext = params.Bext;
    lambda = params.lambda;
    f0= params.f0;
    
    
    % Calculate Ms levels
    kxx_ = 2*pi./lambda;
    Ms_test = (20e3):100:(Bext/mu0);

    f = zeros(size(Ms_test));
    MsLevels = zeros(size(lambda));

    j = 1;
    for kxx = kxx_ % Iterate through the measured wavelength
        % Calculate f-Ms data
        jj = 1;
        for Ms = Ms_test
            H0 = Bext/mu0-Ms; % effective field
            omegaM = gamma*mu0*Ms;
            omegaHx = gamma*mu0*(H0+2*A_exch/(mu0*Ms).*kxx.^2);
            Px = 1-(1-exp(-abs(kxx)*d))./(abs(kxx)*d);
            Phi2 = pi/2*0;
            omegax = sqrt(omegaHx.*(omegaHx+omegaM*(Px+sin(theta)^2*(1-Px.*(1+cos(Phi2).^2)+omegaM./omegaHx.*(Px.*(1-Px).*sin(Phi2).^2)))));
            f(jj) = omegax/2/pi;
            jj=jj+1;
        end

        [f_Ms_fit,~,mu] = polyfit(f,Ms_test,23); % Fit curve to f Ms data
        Ms_calculated = polyval(f_Ms_fit,f0(j),[],mu); % Get the Ms at f0 frequncy
        MsLevels(j) = Ms_calculated;
        j=j+1;
    end

%     figure(164)
%     plot(f,Ms_test,'.')
%     hold on;
%     plot(f,polyval(f_Ms_fit,f,[],mu))
%     plot(f0,polyval(f_Ms_fit,f0,[],mu),'o')
%     hold off;
end