%% Plotting efficiency p2 vs h2
% 20220927 - errik 'andi

clear;

cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\');
addpath(genpath("Libraries"));

p2 = [65  70  75]; lp2 = length(p2);
h2 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh2 = length(h2);

%% Title

% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220922_TandemThursday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20220922_SingleFoil_PHSweep_A3E_p2=75_h2=1.2.mat');
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_SingleFoil_PHSweep_A3E_p2=80_h2=1.2.mat')

% out(:,5) = Prof_out_angle(:,5); % for wallace pitch that came undone
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, fs, 3, 6*chord);

maintitle = ['$f^* =$ ', num2str(par.fred), ', $c =$ ', num2str(foil.chord),...
    ' m, Re = ', num2str(round(par.Re,-4)/1000), 'k, $U_{\infty} =$ ', num2str(round(par.U,4)), ' m/s'];

%% Main data

load('20221011_SingleFoil_efficiency_A3E_PH.mat')

p2 = [65  70  75 80]; lp2 = length(p2);
h2 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh2 = length(h2);

lvlstp = 0.005;

[X, Y] = meshgrid(p2,h2);
% [X, Y] = meshgrid(global_phase,h3);

ylbl = ('$h$');
xlbl = ('$\theta$');

%% Plotting

figure(1)
colormap('viridis')

mt = ['\textbf{Raw} $\eta$, ', maintitle];
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

% parameter = CPP2 + CPH2; % variable to plot, total power coefficient
parameter = Eff_2;

contourf(X, Y, parameter, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;

[M, I] = max(parameter,[],'all','linear');
plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
hold off;

caxis([min(parameter,[],'all'), max(parameter,[],'all')]);
% caxis([0.3, 0.41]);

set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

xlabel(xlbl,'Interpreter','latex');
ylabel(ylbl,'Interpreter','latex');

colorbar();

%%

figure(2)
colormap('magma')

mt = ['\textbf{Corr} $\eta$, ', maintitle];
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

parameter = Eff_2prime; % variable to plot

contourf(X, Y, parameter, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;

[M, I] = max(parameter,[],'all','linear');
plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
hold off;

caxis([min(parameter,[],'all'), max(parameter,[],'all')]);
% caxis([0.3, 0.41]);

set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

xlabel(xlbl,'Interpreter','latex');
ylabel(ylbl,'Interpreter','latex');

colorbar();
