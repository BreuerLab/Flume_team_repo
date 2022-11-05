%% Quick 2 Foil Analysis

addpath(genpath("Libraries"));

foils_analyzed = input('\nSelect analysis: \nLeading foil [1], Trailing foil [2], Both foils [3] : ');

% NOTE: got to fix the blockage correction

% clear;
% 
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220619_TandemSunday_AlphaSweep_APHPhase_A2E_a33_a67\20220619_TandemFoil_APHPhaseSweep_A2E_alpha=0.33_p3=60_h3=0.7c_phase=0.mat')
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=60_h3=1.15c_phase=-180.mat');
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=80_h3=0.55c_phase=180.mat');

% out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220812
foiltype = 'A3E';
[foil, ~, ~] = foils_database(foiltype);

% [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, 3, 6*0.061, flume_height, U_wake); % uncomment when using wake velocity calculation or measurement
[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, transientcycs, 6*foil.chord, flume_height);
% [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, srate, transientcycs);

% out = filter_motor_noise_wallace(out, par.freq, par.srate, 30); % to show nice data, doesn't affect the efficiency calculation
% out = filter_motor_noise_gromit(out, par.freq, par.srate, 30); % to show nice data, doesn't affect the efficiency calculation
res = calculate_forces(par, kin, out);

%% Plotting

if foils_analyzed == 1 || foils_analyzed == 3
    [toverT1, pitch_cyc2, CL2_cyc] = cycle_avg_data(kin.p2_comm, res.CL2); % LiftC2
    [toverT2, pitch_cyc2, CM2_cyc] = cycle_avg_data(kin.p2_comm, res.CM2); % TorqueC2
    [toverT3, pitch_cyc2, CD2_cyc] = cycle_avg_data(kin.p2_comm, (res.CD2) ); % DragC2
    [toverT4, pitch_cyc2, CP2_cyc] = cycle_avg_data(kin.p2_comm, (res.CPH2+res.CPP2) ); % PowerC2
end

if foils_analyzed == 2 || foils_analyzed == 3
    [toverT5, pitch_cyc3, CL3_cyc] = cycle_avg_data(kin.p3_comm, res.CL3); % LiftC3
    [toverT6, pitch_cyc3, CM3_cyc] = cycle_avg_data(kin.p3_comm, res.CM3); % TorqueC3
    [toverT7, pitch_cyc3, CD3_cyc] = cycle_avg_data(kin.p3_comm, (res.CD3) ); % DragC3
    [toverT8, pitch_cyc3, CP3_cyc] = cycle_avg_data(kin.p3_comm, (res.CPH3+res.CPP3)); % PowerC3
end

%% General parameters and results

figure();

sgtitle('\textbf{Experiment Properties}', 'FontSize', 22, 'interpreter', 'latex');

exp_properties = {['$\eta_{le}$ = ', num2str(res.Eff_2,3)], ['$\eta_{tr}$ = ', num2str(res.Eff_3,3)],...
    ['$\eta_{sys}$ = ', num2str((res.Eff_2+res.Eff_3),3)],...
    ['$\eta_{le,corr}$ = ', num2str(res.Eff_2prime,3)], ['$\eta_{tr,corr}$ = ', num2str(res.Eff_3prime,3)],...
    ['$\Phi_{1-2}=$ ', num2str(par.global_phase,3), '$^o$'], ['$\psi_{1-2} =$ ', num2str(par.phase13), '$^o$'], ['$f^* =$ ', num2str(par.fred)]...
    ['$\alpha_{T/4,le}$ = ', num2str(par.alphaT4_2,3)], ['$\alpha_{T/4,tr}$ = ', num2str(par.alphaT4_3,3)],...
    ['$\theta_{le}|\theta_{tr} =$ ', num2str(par.P2,2), '$^o|$', num2str(par.P3, '%.f'), '$^\circ$'],...
    ['$h_{le}|h_{tr} =$ ', num2str(par.H2/foil.chord), '$c|$', num2str(par.H3/foil.chord), '$c$'],...
    ['$c =$ ', num2str(foil.chord), ' m'], ['Re = ', num2str(round(par.Re,-4)/1000),'k'],...
    ['$U$ = ', num2str(round(par.U,3)), ' m/s'], ['$U_{wake}$ = ', num2str(round(par.U_wake,3)), ' m/s']};

t_box = annotation('textbox', [0.3, 0.01, 0.5, 0.9], 'String', exp_properties, 'FitBoxToText', 'off', 'FontSize', 18, 'interpreter', 'latex');
t_box.LineStyle = 'none';

%% Cycle-averaged Force Measurements

figure();

maintitle = ['$\alpha_{T/4}$ = ', num2str(par.alphaT4_2,3),...
    ', $\Phi_{1-2}=$ ', num2str(par.global_phase,3), '$^o, \psi_{1-2} =$ ', num2str(par.phase13), '$^o, f^* =$ ', num2str(par.fred)...
    ', $\theta_{le}|\theta_{tr} =$ ', num2str(par.P2), '$^o|$', num2str(par.P3,2), '$^\circ$',...
    ', $h_{le}|h_{tr} =$ ', num2str(par.H2/foil.chord), '$c|$', num2str(par.H3/foil.chord), '$c$',...
    ', $c =$ ', num2str(foil.chord), ' m, Re = ', num2str(round(par.Re,-4)/1000), 'k, $U_{\infty}$ = ', num2str(round(par.U,4))];

sgtitle(maintitle, 'Interpreter', 'latex', 'FontSize', 24);

if foils_analyzed == 1 || foils_analyzed == 3
subplot(2,4,1) % leading lift
yyaxis left
ylabel('$C_{L,le}$', 'Interpreter', 'latex')
shadedErrorBar(toverT1, CL2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT1, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,2) % leading torque
yyaxis left
ylabel('$C_{M,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT2, CM2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT2, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,3) % leading drag
yyaxis left
ylabel('$C_{D,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT3, CD2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT3, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,4) % leading power
yyaxis left
ylabel('$C_{P,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT4, CP2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT4, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
end

if foils_analyzed == 2 || foils_analyzed == 3
subplot(2,4,5) % trailing lift
yyaxis left
ylabel('$C_{L,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT5, CL3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT5, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,6) % trailing torque
yyaxis left
ylabel('$C_{M,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT6, CM3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT6, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,7) % trailing drag
yyaxis left
ylabel('$C_{D,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT7, CD3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT7, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,4,8) % trailing power
yyaxis left
ylabel('$C_{P,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT8, CP3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT8, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
end

%% PSD

% par.freq = 0.1848;
par.srate = 1000;

figure()

[pxxFn_tr,f1] = pwelch(out(:,7),[],[],[],par.srate); % trailing normal force
[pxxFt_tr,f2] = pwelch(out(:,8),[],[],[],par.srate);
[pxxFn_le,f3] = pwelch(out(:,17),[],[],[],par.srate);
[pxxFt_le,f4] = pwelch(out(:,18),[],[],[],par.srate);

plot(f1,10*log10(pxxFn_le), 'linewidth', 1.5); hold on;
plot(f2,10*log10(pxxFt_le), 'linewidth', 1.5);
plot(f3,10*log10(pxxFn_tr), 'linewidth', 1.5);
plot(f4,10*log10(pxxFt_tr), 'linewidth', 1.5); hold off;

legend('$F_{N,tr}$', '$F_{T,tr}$', '$F_{N,le}$', '$F_{T,le}$', 'interpreter', 'latex', 'fontsize', 18);
xlabel('Hz', 'interpreter', 'latex')
ylabel('PSD', 'interpreter', 'latex')
% xlim([0, par.freq*30])
set(gca,'FontSize',16, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');
