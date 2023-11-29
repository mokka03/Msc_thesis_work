clear;
load('data/m_intr_slice_meas.mat')
load('models/MagnSim/my_2141_0.mat')

figure(73)
image(my,'CDataMapping','scaled')

%%
m_slice_sim = sum(my(400:600,35:end-100));
m_slice_sim = m_slice_sim/max(m_slice_sim);

l = length(m_slice_meas);
x = linspace(0,l/10,l);

figure(74)
clf;
plot(x(1:566),m_slice_sim)
hold on;
plot(x,m_slice_meas)
hold off;
xlabel('x (\mum')
ylabel('m_{norm}')
set(gca,'FontSize',15)

% SaveFig('figure/FIBregions_2050MHz_m5dbm_Hpos95p7/','damping', gcf);

lambda1 = WavelengthFourier(m_slice_sim,100e-9)
lambda2 = WavelengthFourier(m_slice_meas,100e-9)