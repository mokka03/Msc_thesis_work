function [lambda] = WavelengthFourier(m, dx)
%WAVELENGTHFOURIER Summary of this function goes here
%   Detailed explanation goes here

Fs = 1/dx;  % Sampling frequence
P = zeros(1,round(size(m,2)/2)+1);

for i = 1:size(m,1)
    slice = m(i,:);
    % Fourier transform
    L = length(slice);
    Y = fft(slice);

    P2 = abs(Y/L);
    Pi = P2(1:L/2+1);
    Pi(2:end-1) = 2*Pi(2:end-1);

    P = P + Pi/max(Pi);
end

f = Fs*(0:(L/2))/L; % Frequencies

% Plot
f1 = figure('visible','off');
plot(1./f*1e6,P/max(P));
title("Wavelength in x direction");
xlabel('Wavelength (\mum)');
set(gca,'FontSize',15);
xlim([0 20]);


 [~, argmax] = max(P);
 lambda = 1/f(argmax);
end

