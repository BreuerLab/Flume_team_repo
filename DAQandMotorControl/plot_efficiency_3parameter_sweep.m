%% Plotting uncorrected efficiency

clear;

cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\');
addpath(genpath("Libraries"));

% load('20220619_TandemFoil_efficiency_A2E_a15_PHPh.mat');
load('20220619_TandemFoil_efficiency_A2E_a33_PHPh.mat');
% load('20220619_TandemFoil_efficiency_A2E_a68_PHPh.mat');

global_phase = rad2deg((2*pi*6*0.0762)/(0.39*(1/0.63)) + deg2rad(ph));

lvlstp = 0.005;

[X, Y] = meshgrid(p3,h3);
% [X, Y] = meshgrid(global_phase,h3);

ylbl = ('$h_{tr}$');
xlbl = ('$\theta_{tr}$');
% xlbl = ('$\Phi_{1-2}$');
% xlbl = ('$\psi_{1-2}$');

%%
figure(1)

mt = ['Corrected $\eta_{tr}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_corr_2,2)-1
    
    Eff3 = squeeze(Eff_corr_3(:,i,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff_phys_3,[],'all'), max(Eff_phys_3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\psi_{1-2}$ = ',num2str(ph(i)), '$^o$, $\eta_{tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(2)

mt = ['$\eta_{sys}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,2)-1
    
    Eff3 = squeeze(Eff_sys(:,i,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff_sys,[],'all'), max(Eff_sys,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\psi_{1-2}$ = ',num2str(ph(i)), '$^o$, $\eta_{sys}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(3)

mt = ['$\eta_{sys}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,2)-1
    
    Eff3 = squeeze(Eff_sys_corr(:,i,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff_sys_corr,[],'all'), max(Eff_sys_corr,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\psi_{1-2}$ = ',num2str(ph(i)), '$^o$, $\eta_{sys}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

