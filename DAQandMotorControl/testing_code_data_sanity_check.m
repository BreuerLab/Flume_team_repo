%% Testing code - Data sanity check
% 2022/10/28
% a simple code thqat runs a simple experiment to check basic calculations

foiltype = 'A3E';

[foil, ~, ~] = foils_database(foiltype);
chord = foil.chord;
span = foil.span;

U = 0.33;

fred = [0.08,0.09,0.1,0.11,0.12,0.13,0.14];
freq = fred*U/chord;
% freq = 1;

H2 = 0;
P2 = 0;
H3 = 0.8;
P3 = 70;

num_cyc = 30;

for trial = 1:length(fred)

    % Take another force sensor tare measurement right before the trial starts
    [~,bias_newloaded,~] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil);
    % bias_trial -> bias_new - bias_loaded + bias
    bias_trial.Wallace = bias_newloaded.Wallace - bias_loaded.Wallace + bias.Wallace;
    bias_trial.Gromit = bias_newloaded.Gromit - bias_loaded.Gromit + bias.Gromit;
    bias_trial.accmeter = bias_newloaded.accmeter - bias_loaded.accmeter + bias.accmeter;
    bias_trial.pitch = bias.pitch;
    
    [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq_trial, pitch2, heave2, pitch3, heave3,phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,bias_trial,foiltype, freq(trial), P2, H2, P3, H3, 0, -90,...
        num_cyc, 3, 0, 0, 0, 5);
    
    current_time = clock;
    
%     fx_w = mean(out(:,17)); fx_w_std = std(out(:,17));
%     fy_w = mean(out(:,18)); fy_w_std = std(out(:,18));
%     fz_w = mean(out(:,19)); fz_w_std = std(out(:,19));
%     Mx_w = mean(out(:,20)); Mx_w_std = std(out(:,20));
%     My_w = mean(out(:,21)); My_w_std = std(out(:,21));
%     Mz_w = mean(out(:,22)); Mz_w_std = std(out(:,22));
%     
%     fx_g = mean(out(:,7)); fx_g_std = std(out(:,7));
%     fy_g = mean(out(:,8)); fy_g_std = std(out(:,8));
%     fz_g = mean(out(:,9)); fz_g_std = std(out(:,9));
%     Mx_g = mean(out(:,10)); Mx_g_std = std(out(:,10));
%     My_g = mean(out(:,11)); My_g_std = std(out(:,11));
%     Mz_g = mean(out(:,12)); Mz_g_std = std(out(:,12));
    
    filename = ['20221102_BiasTesting_', foiltype, 'fred=', num2str(fred(trial),3), '_p3=', num2str(pitch3,2), '_h3=', num2str(H3,2), 'c_U=0.33.mat'];
    save(filename);

end

[pxx, f] = pwelch(out(:,17),[],[],1000);
plot(f,10*log10(pxx));

%%


% load(filename);
% 
% raaange = 1:length(out);
% % raaange = 6364:125713;
% 
% [foil, ~, ~] = foils_database(foiltype);
% chord = foil.chord;
% span = foil.span;
% 
% %% leading foil
% 
% heave_commanded = Prof_out_angle(raaange,4);
% heave_measured = Prof_out_angle(raaange,4);
% % heave_star_measured = heave_measured/chord;
% pitch_measured = deg2rad(Prof_out_angle(raaange,3));
% moment_z0 = -out(raaange,22); % negative due to the change in frame of reference
% force_x0 = out(raaange,17);
% force_y0 = out(raaange,18);
% force_D = force_y0.*cos(pitch_measured) - force_x0.*sin(pitch_measured);
% force_L = force_x0.*cos(pitch_measured) + force_y0.*sin(pitch_measured);
% inertialload_y = out(raaange,23);
% flowspeed_measured = mean(abs(out(raaange,13)));
% 
% T = 1/1000;
% 
% heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),100);
% heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),100);
% 
% pitch_velo = movmean((1/T)*gradient(squeeze(pitch_measured)),100);
% pitch_accel = movmean((1/T)*gradient(squeeze(pitch_velo)),100);
% 
% P_p = force_L.*heave_velo;
% P_h = moment_z0.*pitch_velo;
% 
% yp1 = heave_measured + chord*0.5*sin(pitch_measured);
% yp2 = heave_measured - chord*0.5*sin(pitch_measured);
% Yp = 2*max(max(yp1),max(yp2));
% 
% P_flow = 0.5*1000*flowspeed_measured^3*Yp*span;
% 
% Eff_le = mean(P_p + P_h)/P_flow
% L_max_le = max(abs(force_L))
% 
%% trailing foil

heave_commanded = Prof_out_angle(raaange,6);
heave_measured = out(raaange,6);
% heave_star_measured = heave_measured/chord;
pitch_measured = out(raaange,5);
moment_z0 = -out(raaange,12); % negative due to the change in frame of reference
force_x0 = out(raaange,7);
force_y0 = out(raaange,8);
force_D = force_y0.*cos(pitch_measured) - force_x0.*sin(pitch_measured);
force_L = force_x0.*cos(pitch_measured) + force_y0.*sin(pitch_measured);
inertialload_y = out(raaange,23);
flowspeed_measured = mean(abs(out(raaange,13)));

T = 1/1000;

heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),100);
heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),100);

pitch_velo = movmean((1/T)*gradient(squeeze(pitch_measured)),100);
pitch_accel = movmean((1/T)*gradient(squeeze(pitch_velo)),100);

P_p = force_L.*heave_velo;
P_h = moment_z0.*pitch_velo;

yp1 = heave_measured + chord*0.5*sin(pitch_measured);
yp2 = heave_measured - chord*0.5*sin(pitch_measured);
Yp = 2*max(max(yp1),max(yp2));

P_flow = 0.5*1000*flowspeed_measured^3*Yp*span;

Eff_tr = mean(P_p + P_h)/P_flow
L_max_tr = max(abs(force_L))
% 
% 
