%% Quick Analysis

% NOTE: got to fix the blockage correction

% clear;
% 
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo\20220619_TandemSunday_AlphaSweep_APHPhase_A2E_a33_a67\20220619_TandemFoil_APHPhaseSweep_A2E_alpha=0.33_p3=60_h3=0.7c_phase=0.mat')
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=60_h3=1.15c_phase=-180.mat');
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\20220617_TandemFoil_PHPhaseSweep_A2E_p3=80_h3=0.55c_phase=180.mat');

out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220617 - 20220622

addpath(genpath("Libraries"));

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out);
res = calculate_forces(par, kin, out);

%% Plotting

[toverT1, pitch_cyc2, CL2_cyc] = cycle_avg_data(kin.p2_comm, res.CL2); % LiftC2
[toverT2, pitch_cyc2, CM2_cyc] = cycle_avg_data(kin.p2_comm, res.CM2); % TorqueC2
[toverT3, pitch_cyc2, CP2_cyc] = cycle_avg_data(kin.p2_comm, (res.CPH2+res.CPP2) ); % PowerC2

[toverT4, pitch_cyc3, CL3_cyc] = cycle_avg_data(kin.p3_comm, res.CL3); % LiftC3
[toverT5, pitch_cyc3, CM3_cyc] = cycle_avg_data(kin.p3_comm, res.CM3); % TorqueC3
[toverT6, pitch_cyc3, CP3_cyc] = cycle_avg_data(kin.p3_comm, (res.CPH3+res.CPP3)); % PowerC3

%% Force Measurements

figure();

subplot(2,3,1)
yyaxis left
shadedErrorBar(toverT1, CL2_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT1, mean(pitch_cyc2));
title('Leading C_L');

subplot(2,3,2)
yyaxis left
shadedErrorBar(toverT2, CM2_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT2, mean(pitch_cyc2));
title('Leading C_M');

subplot(2,3,3)
yyaxis left
shadedErrorBar(toverT3, CP2_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT3, mean(pitch_cyc2));
title('Leading C_P');

subplot(2,3,4)
yyaxis left
shadedErrorBar(toverT4, CL3_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT4, mean(pitch_cyc3));
title('Trailing C_L');

subplot(2,3,5)
yyaxis left
shadedErrorBar(toverT5, CM3_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT5, mean(pitch_cyc3));
title('Trailing C_M');

subplot(2,3,6)
yyaxis left
shadedErrorBar(toverT6, CP3_cyc,{@mean, @std},'lineprops','-b','transparent',true,'patchSaturation',0.2);
yyaxis right
plot(toverT6, mean(pitch_cyc3));
title('Trailing C_P');

maintitle = ['Leading eff = ', num2str(res.Eff_2), ', Trailing eff = ', num2str(res.Eff_3)];

sgtitle(maintitle);

