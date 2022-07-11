%% Quick Trailing Foil Analysis

addpath(genpath("Libraries"));

% NOTE: got to fix the blockage correction

% clear;
% 
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220619_TandemSunday_AlphaSweep_APHPhase_A2E_a33_a67\20220619_TandemFoil_APHPhaseSweep_A2E_alpha=0.33_p3=60_h3=0.7c_phase=0.mat')
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=60_h3=1.15c_phase=-180.mat');
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=80_h3=0.55c_phase=180.mat');

out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220617 - 20220622
foil_separation = 6;

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, 3, foil_separation);
% [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, srate, transientcycs);
% out = filter_motor_noise_gromit(out, par.freq, par.srate, 30); % to show nice data, doesn't affect the efficiency calculation
res = calculate_forces(par, kin, out);

%% Plotting

[toverT5, pitch_cyc3, CL3_cyc] = cycle_avg_data(kin.p3_comm, res.CL3, 1000, 1); % LiftC3
[toverT6, pitch_cyc3, CM3_cyc] = cycle_avg_data(kin.p3_comm, res.CM3, 1000, 1); % TorqueC3
[toverT7, pitch_cyc3, CD3_cyc] = cycle_avg_data(kin.p3_comm, res.CD3, 1000, 1); % DragC3
[toverT8, pitch_cyc3, CP3_cyc] = cycle_avg_data(kin.p3_comm, (res.CPH3+res.CPP3), 1000, 1); % PowerC3

%% Cycle-averaged Force Measurements

figure();

maintitle = ['$\eta_{le}$ = ', num2str(res.Eff_2,3), ', $\eta_{tr}$ = ', num2str(res.Eff_3,3),...
    ', $\Phi_{1-2}=$ ', num2str(par.global_phase,3), '$^o, \psi_{1-2} =$ ', num2str(par.phase13), '$^o, f^* =$ ', num2str(par.fred)...
    ', $\theta_{le}|\theta_{tr} =$ ', num2str(par.P2), '$^o|$', num2str(par.P3), '$^o$',...
    ', $h_{le}|h_{tr} =$ ', num2str(par.H2/foil.chord), '$c|$', num2str(par.H3/foil.chord), '$c$',...
    ', $c =$ ', num2str(foil.chord), ' m, Re = ', num2str(round(par.Re,-4)/1000), 'k'];

sgtitle(maintitle, 'Interpreter', 'latex', 'FontSize', 24);


subplot(2,2,1)
yyaxis left
ylabel('$C_{L,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT5, CL3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT5, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,2)
yyaxis left
ylabel('$C_{M,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT6, CM3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT6, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,3)
yyaxis left
ylabel('$C_{D,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT7, CD3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT7, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,4)
yyaxis left
ylabel('$C_{P,tr}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT8, CP3_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT8, mean(pitch_cyc3));
ylabel('$\theta(t/T)_{tr}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


