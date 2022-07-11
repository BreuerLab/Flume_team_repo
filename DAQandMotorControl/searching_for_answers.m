%% Searching for answers

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, EP.transientcycs, foil_separation);
res = calculate_forces(par, kin, out);

alpha_e = atan(-kin.h2_vel./par.U) + kin.p2_meas;

figure(1)

plot(kin.p2_meas); hold on;
plot(alpha_e); hold off;

[toverT, pitch_cycle, data_cycle] = cycle_avg_data(kin.p2_comm, alpha_e, 1000, 0);

figure(2)
shadedErrorBar(toverT, data_cycle,{@mean, @std},'lineprops',{'-b','LineWidth',1},'transparent',true,'patchSaturation',0.2);