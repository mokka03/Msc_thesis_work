clear;

basedir = 'python\models\demux_v6_frequency_sweep\';
f_MHz = 1950:5:2200;
desired_outs = zeros(2,length(f_MHz));

for i = 1:length(f_MHz)
    %%% Load data
    load([basedir 'outputs_' num2str(f_MHz(i)) 'MHz.mat']);
    desired_outs(1,i) = outputs(5);
    desired_outs(2,i) = outputs(12);
    %%% Normalize
    desired_outs(:,i) = normalize(desired_outs(:,i),"norm",1);
end

figure(735)
clf;
hold on;
plot(f_MHz*1e-3,desired_outs(1,:),'.-');
plot(f_MHz*1e-3,desired_outs(2,:),'.-');
hold off;
xlabel('Frequency (GHz)')
ylabel('Ratio of power on the outputs')
xlim([2 2.15])
legend('Output 1','Output 2')
set(gca,'FontSize',13)

% SaveFig('figure/demux_v6/spintorch/','fsweep',gcf);