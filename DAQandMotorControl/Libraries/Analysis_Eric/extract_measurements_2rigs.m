%% Extract Measurements 2 Rigs

% EP --> Experiment Parameters (saved after the experiment)

function    [kin, res, par] = extract_measurements_2rigs(dq, foiltype, Prof_out_angle, out, EP)

    kin.T = 1/dq.Rate; % cycle period
    g = 9.80665; % acceleration of gravity
    nu = 1.00e-6; % water kinematic viscosity @22 degC (20220531)
    [foil, ~, ~] = foils_database(foiltype); % foil-specific parameters
    
    if ~exist('EP.flume_height','var')
        % parameter does not exist, default it to the following:
        EP.flume_height = 0.55; % assumed flume water depth [m]
    end
    
    % Define timesteps for each subtrial, excluding ramp up/down
    
    tsteps = length(out); % total timesteps in the experiment
    tstep_start = round(transientcycs/(freq*T))+1; % first step after ramp up transient
    tstep_end = round(tsteps-transientcycs/(freq*T)); % first step after ramp down transient
    times = T*(1:tstep_end-tstep_start+1);
    kin.time_star = times*EP.freq; % time per cycle
    
    % Kinematics ("kin" variable)
    
    kin.p2_comm = Prof_out_angle(tstep_start:tstep_end,3); % leading commanded heave [deg]
    kin.h2_comm = Prof_out_angle(tstep_start:tstep_end,4); % leading commanded pitch [m]
    kin.p3_comm = Prof_out_angle(tstep_start:tstep_end,5); % trailing commanded heave [deg]
    kin.h3_comm = Prof_out_angle(tstep_start:tstep_end,6); % trailing commanded pitch [m]
    
    kin.p2_meas = out(tstep_start:tstep_end,3); % leading measured pitch [rad]
    kin.h2_meas = out(tstep_start:tstep_end,4); % leading measured heave [m]
    kin.p3_meas = out(tstep_start:tstep_end,5); % trailing measured pitch [rad]
    kin.h3_meas = out(tstep_start:tstep_end,6); % trailing measured heave [m]
    
    kin.p2_vel = movmean((1/T)*gradient(squeeze(kin.p2_meas)),100); % leading pitch velocity (smoothed with moving mean)
    kin.h2_vel = movmean((1/T)*gradient(squeeze(kin.h2_meas)),100); % leading heave velocity (smoothed with moving mean)
    kin.h2_acc = movmean((1/T)*gradient(squeeze(kin.h2_acc)),100); % leading heave acceleration (smoothed with moving mean)
    
    kin.p3_vel = movmean((1/T)*gradient(squeeze(kin.p3_meas)),100); % trailing pitch velocity (smoothed with moving mean)
    kin.h3_vel = movmean((1/T)*gradient(squeeze(kin.h3_meas)),100); % trailing heave velocity (smoothed with moving mean)
    kin.h3_acc = movmean((1/T)*gradient(squeeze(kin.h3_acc)),100); % trailing heave acceleration (smoothed with moving mean)
    
    % Resulting forces ("res" variable)
    
    res.force2_x0 = out(tstep_start:tstep_end,17); % leading normal force [N]
    res.force2_y0 = out(tstep_start:tstep_end,18); % leading tangential force [N]
    res.torque2_z0 = out(tstep_start:tstep_end,22); % leading pitch axis torque [N*m]
    res.force3_x0 = out(tstep_start:tstep_end,7); % trailing normal force [N]
    res.force3_y0 = out(tstep_start:tstep_end,8); % trailing tangential force [N]
    res.torque3_z0 = out(tstep_start:tstep_end,12); % trailing pitch axis torque [N*m]
    
    res.inertialload_y = out(tstep_start:tstep_end,23); % trailing inertial load [N?]
    
    res.Drag2 = res.force2_y0.*cos(kin.p2_meas) - res.force2_x0.*sin(kin.p2_meas); % leading Drag force [N]
    res.Lift2 = res.force2_x0.*cos(kin.p2_meas) + res.force2_y0.*sin(kin.p2_meas) - (foil.mass1+0.6)*kin.h2_acc; % leading Lift force [N]
    res.Drag3 = res.force3_y0.*cos(kin.p3_meas) - res.force3_x0.*sin(kin.p3_meas); % trailing Drag force [N]
    res.Lift3 = res.force3_x0.*cos(kin.p3_meas) + res.force3_y0.*sin(kin.p3_meas) - (foil.mass2+0.6)*kin.h3_acc; % trailing Lift force [N]
    
    % Parameters calculated
    
    par.U = abs(out(tstep_start:tstep_end,13)); % flow speed from the vectrino [m/s]
    EP.fred = (EP.freq*foil.chord)/par.U; % reduced frequency
    
    par.Re = par.U*foil.chord/nu; % Reynolds number
    par.Fr = par.U/sqrt(g*EP.flume_height); % Froude number
    par.St = EP.freq*foil.chord/par.U; % Strouhal number
    % St = EP.fred*H2/(pi*foil.chord); % Another definition of St for pluging flight, considering the reduced frequency
    %                                % and the heaving amplitude as the characteristic length.
    par.alphaT4 = atan(-2*pi*(EP.H2/foil.chord)*EP.fred) + deg2rad(EP.P2); % relative angle of attack at T/4 (only relevant to the front foil)
    
end
