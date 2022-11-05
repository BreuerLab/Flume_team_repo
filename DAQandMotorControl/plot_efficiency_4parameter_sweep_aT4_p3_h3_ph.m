%% Plotting efficiency aT4 vs p3 vs h3 vs ph
% 20220815 - errik 'andi

% THIS IS THE BEST ONE TILL NOW 20220929

clear;

cd('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\');
addpath(genpath("Libraries"));

%% Title

% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220929_TandemThursday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20220929_TandemFoil_APHPhaseSweep_A3E_alpha=0.155_p3=65_h3=1.1c_phase=-180.mat');
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221003_TandemMonday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221003_TandemFoil_APHPhaseSweep_A3E_alpha=0.155_p3=65_h3=0.7c_phase=-180.mat')
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_TandemFoil_APHPhaseSweep_A3E_alpha=0.679_p3=75_h3=1.2c_phase=120.mat')

% out(:,5) = Prof_out_angle(:,5); % for wallace pitch that came undone
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, fs, 3, 6*chord);

maintitle = ['$f^* =$ ', num2str(par.fred), ', $c =$ ', num2str(foil.chord),...
    ' m, Re = ', num2str(round(par.Re,-4)/1000), 'k, $U_{\infty} =$ ', num2str(round(par.U,4)), ' m/s'];

%% Main data

% load('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpWake_EffWake_SysEffBoth.mat');
load('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm_SysEffFrstrm_wBaseline.mat');
% load('20220929_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm.mat')

aT4 = [0.155, 0.33, 0.679]; laT4 = length(aT4);
p3 = [65  70  75]; lp3 = length(p3);
h3 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh3 = length(h3);
ph = [-180  -120   -60     0    60   120]; lph = length(ph);

% ph_rad = [-pi -(2/3)*pi -(1/3)*pi 0 (1/3)*pi (2/3)*pi];
% global_phase = 2*pi*(6/0.86)*0.12 + ph_rad; % based on Yunxing's thesis

lvlstp = 0.005;

% xticks_range = [-pi -(2/3)*pi -(1/3)*pi 0 (1/3)*pi (2/3)*pi];
% xticks_range = [(3/4)*pi, pi, (5/4)*pi, (3/2)*pi, (7/4)*pi, 2*pi, (9/4)*pi, (5/2)*pi];

[X, Y] = meshgrid(ph,h3);
% [X, Y] = meshgrid(ph_rad,h3);
% [X, Y] = meshgrid(global_phase,h3);

ylbl = ('$h^*_{tr} = H_{0,tr}/c$');
xlbl = ('$\psi_{1-2}$ (deg)');

%%

figure(1)
colormap('inferno')

mt = ['\textbf{Measured} $C_{P,tr}$, ', maintitle];
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

for n = 1:size(Yp,1) % every alphaT4
    
    variable = CPP3 + CPH3;
    par = squeeze(variable(n,:,:,:)); % change this variable to the desired one
    m = n;
    for i = 1:size(par,1) % every pitch angle

        par_loop = squeeze(par(i,:,:));
        
        m = n + 3*(i-1); % for the subplot index
        subplot(3,3,m); 
        contourf(X, Y, par_loop, 'LevelStep', 0.02, 'LineStyle', 'none'); hold on;

        [M, I] = max(par_loop,[],'all','linear');
        plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
        hold off;
        
        caxis([min(variable,[],'all'), max(variable,[],'all')]);
        
%         xticks(xticks_range)
%         xticklabels({'$-\pi$','$-2/3\pi$','$-1/3\pi$','$0$','$1/3\pi$','$2/3\pi$','$\pi$'})
%         xticklabels({'$3/4\pi$', '$\pi$', '$5/4\pi$', '$3/2\pi$', '$7/4\pi$', '$2\pi$', '$9/4\pi$', '$5/2\pi$'})

        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        xlabel(xlbl,'Interpreter','latex');
        ylabel(ylbl,'Interpreter','latex');
        
        subtitle = ['$\alpha_{T/4} =$ ', num2str(aT4(n),2), ', $\theta_{tr} =$ ', num2str(p3(i)), '$^\circ$'];
        title(subtitle, 'FontSize', 26, 'Interpreter', 'latex');
    end
end

c = colorbar();
c.TickLabelInterpreter = 'latex';

%%
figure(2)
colormap('magma')

mt = ['\textbf{Measured} $\eta_{tr}$, ', maintitle]; % change this to appropriate
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

for n = 1:size(Yp,1) % every alphaT4
    
    variable = Eff_3;
    par = squeeze(variable(n,:,:,:)); % change this variable to the desired one
    m = n;
    for i = 1:size(par,1) % every pitch angle

        par_loop = squeeze(par(i,:,:));
        
        m = n + 3*(i-1); % for the subplot index
        subplot(3,3,m);
        contourf(X, Y, par_loop, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;

        [M, I] = max(par_loop,[],'all','linear');
        plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
        hold off;
        
        caxis([min(variable,[],'all'), max(variable,[],'all')]);
        
%         xticks(xticks_range)
%         xticklabels({'$-\pi$','$-2/3\pi$','$-1/3\pi$','$0$','$1/3\pi$','$2/3\pi$','$\pi$'})
%         xticklabels({'$3/4\pi$', '$\pi$', '$5/4\pi$', '$3/2\pi$', '$7/4\pi$', '$2\pi$', '$9/4\pi$', '$5/2\pi$'})
        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        xlabel(xlbl,'Interpreter','latex');
        ylabel(ylbl,'Interpreter','latex');
        
        subtitle = ['$\alpha_{T/4} =$ ', num2str(aT4(n),2), ', $\theta_{tr} =$ ', num2str(p3(i)), '$^\circ$'];
        title(subtitle, 'FontSize', 26, 'Interpreter', 'latex');
    end
end

c = colorbar();
c.TickLabelInterpreter = 'latex';

%%
figure(3)
colormap('magma')

mt = ['\textbf{Measured} $\eta_{sys}$, ', maintitle]; % change this to appropriate
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

for n = 1:size(Yp,1) % every alphaT4
    
    variable = Eff_sys;
    par = squeeze(variable(n,:,:,:)); % change this variable to the desired one
    m = n;
    for i = 1:size(par,1) % every pitch angle

        par_loop = squeeze(par(i,:,:));
        
        m = n + 3*(i-1); % for the subplot index
        subplot(3,3,m);
        contourf(X, Y, par_loop, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;
        
        [M, I] = max(par_loop,[],'all','linear');
        plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
        hold off;
        
        caxis([min(variable,[],'all'), max(variable,[],'all')]); % change this to desired variable for correct color scaling
        
%         xticks(xticks_range)
%         xticklabels({'$-\pi$','$-2/3\pi$','$-1/3\pi$','$0$','$1/3\pi$','$2/3\pi$','$\pi$'})
%         xticklabels({'$3/4\pi$', '$\pi$', '$5/4\pi$', '$3/2\pi$', '$7/4\pi$', '$2\pi$', '$9/4\pi$', '$5/2\pi$'})
        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        xlabel(xlbl,'Interpreter','latex');
        ylabel(ylbl,'Interpreter','latex');
        
        subtitle = ['$\alpha_{T/4} =$ ', num2str(aT4(n),2), ', $\theta_{tr} =$ ', num2str(p3(i)), '$^\circ$, $\eta_{sys} =$ ', num2str(M,3)];
        title(subtitle, 'FontSize', 22, 'Interpreter', 'latex');
    end
end

c = colorbar();
c.TickLabelInterpreter = 'latex';

%%
figure(4)
colormap('magma')

mt = ['\textbf{Measured} $\eta_{le}$, ', maintitle]; % change this to appropriate
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

for n = 1:size(Yp,1) % every alphaT4
    
    variable = Eff_2;
    par = squeeze(variable(n,:,:,:)); % change this variable to the desired one
    m = n;
    for i = 1:size(par,1) % every pitch angle

        par_loop = squeeze(par(i,:,:));
        
        m = n + 3*(i-1); % for the subplot index
        subplot(3,3,m);
        contourf(X, Y, par_loop, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;

        [M, I] = max(par_loop,[],'all','linear');
        plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
        hold off;
        
        caxis([min(variable,[],'all'), max(variable,[],'all')]); % change this to desired variable for correct color scaling
        
%         xticks(xticks_range)
%         xticklabels({'$-\pi$','$-2/3\pi$','$-1/3\pi$','$0$','$1/3\pi$','$2/3\pi$','$\pi$'})
%         xticklabels({'$3/4\pi$', '$\pi$', '$5/4\pi$', '$3/2\pi$', '$7/4\pi$', '$2\pi$', '$9/4\pi$', '$5/2\pi$'})
        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        xlabel(xlbl,'Interpreter','latex');
        ylabel(ylbl,'Interpreter','latex');
        
        subtitle = ['$\alpha_{T/4} =$ ', num2str(aT4(n),2), ', $\theta_{tr} =$ ', num2str(p3(i)), '$^\circ$'];
        title(subtitle, 'FontSize', 26, 'Interpreter', 'latex');
    end
end

c = colorbar();
c.TickLabelInterpreter = 'latex';

%%
figure(5)
colormap('viridis')

mt = ['$\bar{U}_{wake}$, ', maintitle]; % change this to appropriate
sgtitle(mt, 'FontSize', 28, 'Interpreter', 'latex');

for n = 1:size(Yp,1) % every alphaT4
    
    variable = Uwake; % CHANGE THIS TO THE DESIRED VARIABLE
    lvlstp = 0.0005; % to account for very small variations in the mean wake flow
    par = squeeze(variable(n,:,:,:));
    m = n;
    mean(mean(mean(par)))
    for i = 1:size(par,1) % every pitch angle

        par_loop = squeeze(par(i,:,:));
        
        m = n + 3*(i-1); % for the subplot index
        subplot(3,3,m);
        contourf(X, Y, par_loop, 'LevelStep', lvlstp, 'LineStyle', 'none'); hold on;

        [M, I] = max(par_loop,[],'all','linear');
        plot(X(I),Y(I),'p','MarkerSize',18,'MarkerEdgeColor','k','MarkerFaceColor','y');
        hold off;
        
        caxis([min(variable,[],'all'), max(variable,[],'all')]); % change this to desired variable for correct color scaling
        
%         xticks(xticks_range)
%         xticklabels({'$-\pi$','$-2/3\pi$','$-1/3\pi$','$0$','$1/3\pi$','$2/3\pi$','$\pi$'})
%         xticklabels({'$3/4\pi$', '$\pi$', '$5/4\pi$', '$3/2\pi$', '$7/4\pi$', '$2\pi$', '$9/4\pi$', '$5/2\pi$'})
        set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
        
        xlabel(xlbl,'Interpreter','latex');
        ylabel(ylbl,'Interpreter','latex');
        
        subtitle = ['$\alpha_{T/4} =$ ', num2str(aT4(n),2), ', $\theta_{tr} =$ ', num2str(p3(i)), '$^\circ$'];
        title(subtitle, 'FontSize', 26, 'Interpreter', 'latex');
    end
end

c = colorbar();
c.TickLabelInterpreter = 'latex';
