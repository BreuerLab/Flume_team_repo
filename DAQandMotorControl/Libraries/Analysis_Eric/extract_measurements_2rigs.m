%% Extract Measurements 2 Rigs V2

% EP --> Experiment Parameters (saved after the experiment)

function    [kin, par, foil] = extract_measurements_2rigsV2(foiltype, Prof_out_angle, out, srate, transientcycs)

    [foil, ~, ~] = foils_database(foiltype);

    if ~exist('srate','var')
        % parameter does not exist, default it to the following:
        srate = 1000; % sampling frequency [Hz]
    end
    
    if ~exist('flume_height','var')
        % parameter does not exist, default it to the following:
        flume_height = 0.55; % assumed flume water depth [m]
    end
    
    if ~exist('transientcycs','var')
        % parameter does not exist, default it to the following:
        transientcycs = 3; % assumed transient cycles (ramp up and down)
    end
    
    kin.T = 1/srate; % duration of each timestep
    Tsamp = 1/srate; % duration of each timestep (sampling time)
    g = 9.80665; % acceleration of gravity
    nu = 1.00e-6; % water kinematic viscosity @22 degC (20220531)
    
    %% Flow speed
    
    U = round(mean(abs(out(:,13))),4)*1.026; % flow speed (with LDV correction)
    
    %% Frequency
    
    if mean(abs(Prof_out_angle(:,4))) == 0
        heave_profile = Prof_out_angle(:,6); % profile used to determine the frequency, named generally to avoid confusion in subsequent calculations
    else
        heave_profile = Prof_out_angle(:,4); % if trailing pitch is zero, use leading pitch to determine frequency
    end
    
    min_peak = max(abs(heave_profile))*0.99; % determines the minimum value that the peaks should have
    [pks,locs] = findpeaks(heave_profile,'MinPeakHeight',min_peak);
    
    if ~exist('num_cyc','var')
%         num_cyc = length(T) - transientcycs; % total number of cycles
        num_cyc = length(locs) - 1; % total number of cycles
    end
    
%     NOTE: This section is for debugging the number of cycles calculation.
%     [pks2,locs2]=findpeaks(pitch_profile);
%     figure;
%     plot(pitch_profile); hold on;
%     plot(locs, pks, 'ro');
%     plot(locs2,pks2,'bx');
%     
%     total_cyc = length(pks2) - 1;
%     exp_cyc = total_cyc - 2*transientcycs;
%     
%     if exp_cyc ~= num_cyc
%         error('Number of cycles incorrect. Verify cycles profile.');
%     end

%     total_cyc = num_cyc + 2*transientcycs; % total number of cycles including transients
    T = NaN(num_cyc,1); % period of each cycle
    time = 0:Tsamp:length(heave_profile)/srate-Tsamp; time = time'; % time of the experiment
    
    for i = 1:num_cyc%start_cyc:num_cyc+transientcycs-1
        T(i) = time(locs(i+1)) - time(locs(i)); % time period for each cycle
    end
    
    period = mean(T);  % mean time-period considering only the cycles that will be phase-averaged
    freq = round(1/period,3); % real frequency from the period
    fred = round(freq*foil.chord/U,2); % reduced frequency
    
    %% Kinematics
    
    % Heave and Pitch amplitudes
    
    P2 = max(abs(Prof_out_angle(:,3))); % [deg]
    H2 = max(abs(Prof_out_angle(:,4))); % [m]
    P3 = max(abs(Prof_out_angle(:,5))); % [deg]
    H3 = max(abs(Prof_out_angle(:,6))); % [m]
    
    % Commanded motion
    
    % Define timesteps for the main experimental data in the run, excluding ramp up/down
    tsteps = length(out); % total timesteps in the experiment
    tstep_start = round(transientcycs*(period/Tsamp)+1); % first step after ramp up transient
    tstep_end = round(tsteps-transientcycs*(period/Tsamp)); % first step after ramp down transient
    tsteps_exp = tstep_end - tstep_start; % useful data timesteps
    
    times = Tsamp*(1:tstep_end-tstep_start+1); % time of the usable experimental data
    time_star = times/period; % non-dim time
    
    p2_comm = Prof_out_angle(tstep_start:tstep_end,3); % leading commanded heave [deg]
    h2_comm = Prof_out_angle(tstep_start:tstep_end,4); % leading commanded pitch [m]
    p3_comm = Prof_out_angle(tstep_start:tstep_end,5); % trailing commanded heave [deg]
    h3_comm = Prof_out_angle(tstep_start:tstep_end,6); % trailing commanded pitch [m]
    
    % Measured motion
    
    p2_meas = out(tstep_start:tstep_end,3); % leading measured pitch [rad]
    h2_meas = out(tstep_start:tstep_end,4); % leading measured heave [m]
    p3_meas = out(tstep_start:tstep_end,5); % trailing measured pitch [rad]
    h3_meas = out(tstep_start:tstep_end,6); % trailing measured heave [m]
    
    p2_vel = movmean((1/Tsamp)*gradient(squeeze(p2_meas)),100); % leading pitch velocity (smoothed with moving mean)
    h2_vel = movmean((1/Tsamp)*gradient(squeeze(h2_meas)),100); % leading heave velocity (smoothed with moving mean)
    h2_acc = movmean((1/Tsamp)*gradient(squeeze(h2_vel)),100); % leading heave acceleration (smoothed with moving mean)
    
    p3_vel = movmean((1/kin.T)*gradient(squeeze(p3_meas)),100); % trailing pitch velocity (smoothed with moving mean)
    h3_vel = movmean((1/kin.T)*gradient(squeeze(h3_meas)),100); % trailing heave velocity (smoothed with moving mean)
    h3_acc = movmean((1/kin.T)*gradient(squeeze(h3_vel)),100); % trailing heave acceleration (smoothed with moving mean)
    
    % Other calculated parameters
    
    Re = U*foil.chord/nu; % Reynolds number
    Fr = U/sqrt(g*flume_height); % Froude number
    St = freq*foil.chord/U; % Strouhal number
    % St = EP.fred*H2/(pi*foil.chord); % Another definition of St for pluging flight, considering the reduced frequency
    %                                % and the heaving amplitude as the characteristic length.
    alphaT4 = atan(-2*pi*(H2/foil.chord)*fred) + deg2rad(P2); % relative angle of attack at T/4 (only relevant to the front foil)
    
    %% Storing variables
    
    par.foiltype = foiltype; % from foils_database [string]
    par.flume_height = flume_height; % depth of water in the flume [m]
    par.num_cyc = num_cyc; % number of full amplitude experimental cycles
    par.transientcycs = transientcycs; % transient cycles (ramp up/down)
    par.U = U; % flow speed [m/s]
    par.srate = srate; % sampling frequency [Hz]
    par.freq = freq; % experiment frequency [Hz]
    par.fred = fred; % reduced frequency
    par.Re = Re; % Reynolds
    par.Fr = Fr; % Froude
    par.St = St; % Strouhal
    par.alphaT4 = alphaT4; % angle of attack at 1/4 stroke
    par.P2 = P2; par.H2 = H2; % leading pitch [deg] and heave [m]
    par.P3 = P3; par.H3 = H3; % trailing pitch [deg] and heave [m]
    
    kin.p2_comm = p2_comm; % commanded leading pitch [deg]
    kin.h2_comm = h2_comm; % commanded trailing heave [m]
    kin.p3_comm = p3_comm; % commanded leading pitch [deg]
    kin.h3_comm = h3_comm; % commanded trailing heave [m]
    
    kin.p2_meas = p2_meas; % measured leading pitch [rad]
    kin.h2_meas = h2_meas; % measured trailing heave [m]
    kin.p3_meas = p3_meas; % measured leading pitch [rad]
    kin.h3_meas = h3_meas; % measured trailing heave [m]
    
    kin.p2_vel = p2_vel; % leading pitch velocity [rad/s]
    kin.h2_vel = h2_vel; % leading heave velocity [m/s]
    kin.h2_acc = h2_acc; % leading heave acceleration [m/s^2]
    
    kin.p3_vel = p3_vel; % trailing pitch velocity [rad/s]
    kin.h3_vel = h3_vel; % trailing heave velocity [m/s]
    kin.h3_acc = h3_acc; % trailing heave acceleration [m/s^2]
    
    kin.time_star = time_star; % non-dimensional time of the usable experimental data
    kin.tsteps = tsteps;
    kin.tstep_start = tstep_start;
    kin.tstep_end = tstep_end;
    kin.tsteps_exp = tsteps_exp;
    
end