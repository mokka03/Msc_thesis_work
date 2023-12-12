clear;

basedir = 'python\models\demux_v5_B215mT_fsweep\';
f_MHz = 2000:5:2250;
desired_outs = zeros(2,length(f_MHz));

for i = 1:length(f_MHz)
    %%% Load data
    load([basedir 'outputs_' num2str(f_MHz(i)) 'MHz.mat']);
    desired_outs(1,i) = outputs(5);
    desired_outs(2,i) = outputs(12);
    %%% Normalize
    desired_outs(:,i) = normalize(desired_outs(:,i),"norm",1);
end

load('data\lambda_sweep.mat')

figure(735)
clf;
box on;
hold on;
plot(f_MHz*1e-3,desired_outs(1,:),'.-');
plot(f_MHz*1e-3,desired_outs(2,:),'.-');
hold off;
xlabel('Frequency (GHz)')
ylabel('Ratio of power on the outputs')
xlim([2.05 2.2])
legend('Output 1','Output 2')
set(gca,'FontSize',13)

% figure(735)
% clf;
% box on;
% hold on;
% plot(lambda_sweep*1e6,desired_outs(1,:),'.-');
% plot(lambda_sweep*1e6,desired_outs(2,:),'.-');
% hold off;
% xlabel('Frequency (GHz)')
% ylabel('Ratio of power on the outputs')
% xlim([1.8 3.5])
% % xlim([lambda_sweep(41) lambda_sweep(11)]*1e6)
% legend('Output 1','Output 2')
% set(gca,'FontSize',13)

% SaveFig('figure/demux_v5_B215mT/spintorch/','fsweep_YIG331',gcf);