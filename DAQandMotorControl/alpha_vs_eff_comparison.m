% Using the calculated wake velocity from the actuator disk model to
% calculate the efficiency of the trailing foil.

% 2022 07 26

% NOTE: the vectrino flow measurement is compared with the corrected
% free-stream flow from the blockage correction of the leading foil, and
% with the calculated wake velocity from the leading foil corrected to
% account for the trailing foil.

% NOTE: The physical efficiency of the trailing foil (normalized by the
% calculated wake velocity from the leading foil) is compared with the
% corrected trailing efficiency (also normalized by the wake calculated
% from the leading foil), and with the corrected trailing efficiency
% (normalized by the corrected free-stream flow of the leading foil).

% NOTE: true comparison should be made between uncorrected calculations
% using the physical flow measurements (vectrino and LDV), and the
% corrected calculations.

clear;

cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\');
addpath(genpath("Libraries"));

% load('20220726_TandemFoil_calcWake_efficiency_A2E_a15_PHPh.mat'); index = 7;
load('20220726_TandemFoil_calcWake_efficiency_A2E_a33_PHPh.mat'); index = 7;
% load('20220726_TandemFoil_calcWake_efficiency_A2E_a68_PHPh.mat'); index = 7;

alphaT4_w = reshape(alphaT4_3(:,index,:),[],1);
U_vectrino = reshape(U_flow(:,index,:),[],1);
U_corrected2 = reshape(U_2prime(:,index,:),[],1);
U_corrected3 = reshape(U_3prime(:,index,:),[],1); 
eff_phy_w = reshape(Eff_phys_3(:,index,:),[],1);
eff_cor_w = reshape(Eff_corr_3(:,index,:),[],1);

% load('20220706_BarnsleyWellicome_TandemFoil_efficiency_A2E_a33_PHPh.mat'); index = 7;
% load('20220720_Houslby_TandemFoil_efficiency_A2E_a33_PHPh.mat'); index = 7;

% alphaT4_b = reshape(alphaT4_3(:,index,:),[],1);
% U_vectrino = reshape(U_flow(:,index,:),[],1);
% U_corrected = reshape(U_3prime(:,index,:),[],1);

Eff_corr_3 = Eff_corr_3.*(U_3prime./U_2prime).^3;
% eff_phy = reshape(Eff_phys_3(:,index,:),[],1);
eff_cor = reshape(Eff_corr_3(:,index,:),[],1);

%% Plotting

figure(1)


subplot(1,2,1) % flow velocities

plot(alphaT4_w,U_vectrino, 'bo', 'MarkerSize', 16, 'linewidth', 1.5); hold on;
plot(alphaT4_w,U_corrected2, 'kx', 'MarkerSize', 16, 'linewidth', 1.5);
plot(alphaT4_w,U_corrected3, 'r+', 'MarkerSize', 16, 'linewidth', 1.5); hold off;

% ylim([0.37,0.5])
xlabel('$\alpha_{T/4}$', 'Interpreter', 'latex');
ylabel('$U$ [m/s]','Interpreter','latex');

legend("MEASURED $U_{\infty}$", "CORRECTED $U_{\infty}$", "CORRECTED $U_{wake}$",...
    'Interpreter', 'latex', 'Location', 'West')

set(gca,'FontSize', 22, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(1,2,2) % efficiencies

plot(alphaT4_w,eff_phy_w, 'm^', 'MarkerSize', 16, 'linewidth', 1.5); hold on;
plot(alphaT4_w,eff_cor_w, 'r+', 'MarkerSize', 16, 'linewidth', 1.5);
plot(alphaT4_w,eff_cor, 'kx', 'MarkerSize', 16, 'linewidth', 1.5); hold off;

% ylim([0.37,0.5])
xlabel('$\alpha_{T/4}$', 'Interpreter', 'latex');
ylabel('$\eta_{tr}$','Interpreter','latex');

legend("MEASURED $\eta_{wake,tr}$", "CORRECTED $\eta_{wake,tr}$", "CORRECTED $\eta_{\infty,tr}$",...
    'Interpreter', 'latex', 'Location', 'Southeast')

set(gca,'FontSize', 22, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

