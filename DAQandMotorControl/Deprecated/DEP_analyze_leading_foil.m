%% Quick Leading Foil Analysis
% DEPRECATED

addpath(genpath("Libraries"));

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, 1000, 3, 6);
% out = filter_motor_noise_gromit(out, par.freq, par.srate, 40); % to show nice data, doesn't affect the efficiency calculation
res = calculate_forces(par, kin, out);

%% Plotting

[toverT1, pitch_cyc2, CL2_cyc] = cycle_avg_data(kin.p2_comm, res.CL2); % LiftC2
[toverT2, pitch_cyc2, CM2_cyc] = cycle_avg_data(kin.p2_comm, res.CM2); % TorqueC2
[toverT3, pitch_cyc2, CD2_cyc] = cycle_avg_data(kin.p2_comm, (res.CD2) ); % DragC2
[toverT4, pitch_cyc2, CP2_cyc] = cycle_avg_data(kin.p2_comm, (res.CPH2+res.CPP2) ); % PowerC2

%% Cycle-averaged Force Measurements

figure();

maintitle = ['$\eta_{le}$ = ', num2str(res.Eff_2,3), ', $\eta_{corr}$ = ', num2str(res.Eff_2prime,3),...
    ', $\Phi_{1-2}=$ ', num2str(par.global_phase,3), '$^o, \psi_{1-2} =$ ', num2str(par.phase13), '$^o, f^* =$ ', num2str(par.fred)...
    ', $\theta_{le}|\theta_{tr} =$ ', num2str(par.P2), '$^o|$', num2str(par.P3), '$^o$',...
    ', $h_{le}|h_{tr} =$ ', num2str(par.H2/foil.chord), '$c|$', num2str(par.H3/foil.chord), '$c$',...
    ', $c =$ ', num2str(foil.chord), ' m, Re = ', num2str(round(par.Re,-3)/1000), 'k'];

sgtitle(maintitle, 'Interpreter', 'latex', 'FontSize', 24);


subplot(2,2,1)
yyaxis left
ylabel('$C_{L,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT1, CL2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT1, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,2)
yyaxis left
ylabel('$C_{M,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT2, CM2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT2, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,3)
yyaxis left
ylabel('$C_{D,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT3, CD2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT3, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


subplot(2,2,4)
yyaxis left
ylabel('$C_{P,le}$', 'FontSize', 18, 'Interpreter', 'latex')
shadedErrorBar(toverT4, CP2_cyc,{@mean, @std},'lineprops',{'-b','LineWidth',1.5},'transparent',true,'patchSaturation',0.2);

yyaxis right
plot(toverT4, mean(pitch_cyc2));
ylabel('$\theta(t/T)_{le}$', 'Interpreter', 'latex')

set(gca,'FontSize',18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');


