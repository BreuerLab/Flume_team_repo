%% LDV Profiles Data Export 20220901
% masked trigger profile

clear;

% From the SPEED file:
%
% 1     Device Time     (msec)
% 2     Device Time     (usec)
% 3     Speed           (m/s)
% 4     SNR             

%% load data

% LDV data directory:
FOLDERNAME = ('R:\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220915_LDV_trigger_testing\');
% Speed data file name:
FILE = ('trigger_long_daq_test_2022-09-15-16-23-01.SPEED.MSEBP.txt');
FILENAME = fullfile(FOLDERNAME,FILE); % full file name for loading data
LDV = importdata(FILENAME); % importing data into a double var array

U   = LDV(:,3);
SNR = LDV(:,4);

n = 1;

for i = 1:length(LDV(:,3)) % filtering the low SNR measurements
    if LDV(i,4) > 6
        U_filtered(n,2) = LDV(i,3);
        U_filtered(n,1) = LDV(i,1);
        n = n + 1;
    end
end

U_mean = mean(U);

% [b,a] = butter(6, experiment_frequency*cutoff/(sampling_freq/2), 'low');
% U_ultra_filtered = filtfilt(b, a, U_filtered(:,2));

%% Plots

figure(1)

subplot(1,2,1) % flow velocity overview
plot(U_filtered(:,1), U_filtered(:,2), 'k'); hold on;
yline(U_mean, 'r', 'LineWidth', 3); hold off;
% plot(U_filtered(:,1), U_ultra_filtered, 'xb'); hold off;
xlabel('time (ms)', 'interpreter', 'latex')
ylabel('$U_{\infty} (m/s)$', 'interpreter', 'latex')
legend('$U_{raw}$', '$U_{mean}$', 'interpreter', 'latex', 'location', 'southeast');
set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

subplot(1,2,2) % measurements histogram
histogram(U_filtered(:,2), 'FaceAlpha', 0.4, 'FaceColor', 'b', 'EdgeColor', 'none');
xlabel('$U_{\infty}$ (m/s)', 'interpreter', 'latex');
% xlim([0.15, 0.45])
legend('LDV','Interpreter','latex','location','northwest');
set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

% subplot(2,2,3) % flow velocity zoom-in
% plot(U_filtered(:,1), U_filtered(:,2), 'g'); % hold on;
% % plot(U_filtered(:,1), U_ultra_filtered, 'b', 'linewidth', 2); hold off;
% xlabel('time (ms)', 'interpreter', 'latex')
% ylabel('$U_{\infty}$ (m/s)', 'interpreter', 'latex')
% xlim([2000000, 2300000])
% set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

% subplot(2,2,4) % welch transform
% plot(f,10*log10(pxxLDV), 'k', 'linewidth', 2);
% xlabel('Hz', 'interpreter', 'latex')
% ylabel('PSD', 'interpreter', 'latex')
% xlim([0, experiment_frequency*60])
% set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

