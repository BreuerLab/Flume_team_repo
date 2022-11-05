%% LDV Stations Data Export 20220920
% synchronous profile

clear;

% From the SPEED file:
%
% 1     Device Time     (msec)
% 2     Device Time     (usec)
% 3     Speed           (m/s)
% 4     SNR             

%% load data

% LDV data directory:
FOLDERNAME = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220915_LDV_profile_long_test\Stations');
cd(FOLDERNAME);
% Speed data file name:
FILE = dir('*.SPEED.MSEBP.txt');

N = 1:1:length(FILE); % total number of files

figure(1); hold on;

for j = 1:1:length(N)
    filename = ['300_experiments_20220915_', num2str(N(j),'%03.0f'), '.SPEED.MSEBP.txt'];
    LDV = importdata(filename);
    
    devT = LDV(:,1);
    U   = LDV(:,3);
    SNR = LDV(:,4);
    
    n = 1;
    
    for i = 1:length(LDV(:,3)) % filtering the low SNR measurements
        if LDV(i,4) > 7
            U_filtered(n,2) = LDV(i,3);
            U_filtered(n,1) = LDV(i,1);
            n = n + 1;
        end
    end
    
    U_mean(j) = mean(U_filtered(:,2));
    U_std(j) = std(U_filtered(:,2));
    
    subplot(1,2,1); hold on;
    histogram(U_filtered(:,2), 'FaceAlpha', 0.4, 'EdgeColor', 'none');
    
end

xlabel('$U_{\infty}$ (m/s)', 'interpreter', 'latex');
ylabel('Measurements', 'interpreter', 'latex');
legend('S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13');
set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

subplot(1,2,2);
errorbar(U_mean,U_std,'-o','MarkerSize',10,'MarkerFaceColor',[0,0,1],'MarkerEdgeColor','none');
xlabel('Stops', 'interpreter', 'latex')
ylabel('$U_{\infty} (m/s)$', 'interpreter', 'latex')
% legend('$U_{\infty}$', 'interpreter', 'latex', 'location', 'southeast');
set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

hold off;
