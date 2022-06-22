%% Quick Analysis

% NOTE: got to fix the blockage correction

% clear;
% 
% load(['R:\ENG_Breuer_Shared\ehandyca\Data_main_repo\20220617_TandemFriday_AlphaSweep_PHPhase_A2E_a15\' ...
%     '20220617_TandemFoil_PHPhaseSweep_A2E_p3=60_h3=0.85c_phase=120.mat'])
% 
% % For data taken on 20220617:
% EP.srate = 1000;
% EP.P2 = EP.pitch2;
% EP.H2 = EP.heave2;
% EP.P3 = EP.pitch3;
% EP.H3 = EP.heave3;
% % continue

addpath(genpath("Libraries"));

[kin, par, foil] = extract_measurements_2rigsV2(foiltype, Prof_out_angle, out, srate, transientcycs);
res = calculate_forces(par, kin, out);

%% Plotting

[toverT1, pitch_cyc2, CL2_cyc] = phase_avg_data(kin.p2_meas, eff.CL2); % LiftC2
[toverT2, pitch_cyc2, CM2_cyc] = phase_avg_data(kin.p2_meas, eff.CM2); % TorqueC2
[toverT3, pitch_cyc2, CP2_cyc] = phase_avg_data(kin.p2_meas, (eff.PwrH2+eff.PwrP2) ); % PowerC2

[toverT4, pitch_cyc3, CL3_cyc] = phase_avg_data(kin.p3_meas, eff.CL3); % LiftC3
[toverT5, pitch_cyc3, CM3_cyc] = phase_avg_data(kin.p3_meas, eff.CM3); % TorqueC3
[toverT6, pitch_cyc3, CP3_cyc] = phase_avg_data(kin.p3_meas, (eff.PwrH3+eff.PwrP3)); % PowerC3

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

maintitle = ['Leading eff = ', num2str(eff.Eff_2_prime), ', Trailing eff = ', num2str(eff.Eff_3_prime)];

sgtitle(maintitle);

