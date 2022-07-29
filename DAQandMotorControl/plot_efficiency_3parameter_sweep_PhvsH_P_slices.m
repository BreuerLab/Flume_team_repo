%% Plotting uncorrected efficiency (phase vs heave in pitch slices)

clear;

cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\');
addpath(genpath("Libraries"));

% load('20220720_Houslby_TandemFoil_efficiency_A2E_a15_PHPh.mat');
load('20220720_Houslby_TandemFoil_efficiency_A2E_a33_PHPh.mat');
% load('20220720_Houslby_TandemFoil_efficiency_A2E_a68_PHPh.mat');

global_phase = rad2deg((2*pi*6*0.0762)/(0.39*(1/0.63)) + deg2rad(ph));

lvlstp = 0.005;

[X, Y] = meshgrid(ph,h3);
% [X, Y] = meshgrid(global_phase,h3);

ylbl = ('$h_{tr}$');
xlbl = ('$\psi_{1-2}$');

%%
figure(8)

mt = ['\textbf{Uncorrected} $\eta_{sys}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_corr_2,1)
    
    Eff3 = squeeze(Eff_sys(i,:,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff3,[],'all'), max(Eff3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\theta_{tr}$ = ',num2str(p3(i)), '$^o$, $\eta_{sys}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(9)

mt = ['\textbf{Corrected} $\eta_{sys}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_corr_2,1)
    
    Eff3 = squeeze(Eff_sys_corr(i,:,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff3,[],'all'), max(Eff3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\theta_{tr}$ = ',num2str(p3(i)), '$^o$, $\eta_{sys}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(3)

mt = ['\textbf{Uncorrected} $\eta_{tr}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,1)
    
    Eff3 = squeeze(Eff_phys_3(i,:,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff3,[],'all'), max(Eff3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\theta_{tr}$ = ',num2str(p3(i)), '$^o$, $\eta_{tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(4)

mt = ['\textbf{Corrected} $\eta_{tr}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,1)
    
    Eff3 = squeeze(Eff_corr_3(i,:,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff3,[],'all'), max(Eff3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\theta_{tr}$ = ',num2str(p3(i)), '$^o$, $\eta_{tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

%%
figure(5)

mt = ['$C_{P,tr}$, $\alpha_{T/4}$ = ', num2str(alphaT4,2)];
sgtitle(mt, 'FontSize', 20, 'Interpreter', 'latex');
colormap('magma')

for i = 1:size(Eff_sys,1)
    
    Eff3 = squeeze(CPP3(i,:,:)+CPH3(i,:,:))';

    subplot(2,3,i);
    contourf(X, Y, Eff3, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
    
    % Find maximum
    
    [M, I] = max(Eff3,[],'all','linear');
%     xline(X(I),'--y');
%     yline(Y(I),'--y');
    plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
    
    caxis([min(Eff3,[],'all'), max(Eff3,[],'all')]);
    set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
    
    xlabel(xlbl,'Interpreter','latex');
    ylabel(ylbl,'Interpreter','latex');
    
    t = ['$\theta_{tr}$ = ',num2str(p3(i)), '$^o$, $C_{P,tr}$ = ', num2str(M,3)];
    title(t, 'Interpreter', 'latex', 'FontSize', 18);
    
end

colorbar();

