function    [time_star,heave_commanded,heave_measured,heave_star_measured,pitch1_measured,pitch2_measured,force_D,force_L,inertialload_y,...
    torque_x0,flowspeed_measured,heave_velo,heave_accel] = extract_measurements(transientcycs,freq,T,Prof_out_angle,out,thcknss)

    % Define timesteps for each subtrial, excluding ramp up/down
    timesteps = length(out);
    timestep_start = round(transientcycs/(freq*T))+1;
    timestep_end = round(timesteps-transientcycs/(freq*T));
    times = T*(1:timestep_end-timestep_start+1);
    time_star = times*freq;

    heave_commanded = Prof_out_angle(timestep_start:timestep_end,6);
    heave_measured = out(timestep_start:timestep_end,6);
    heave_star_measured = heave_measured/thcknss;
    pitch1_measured = out(timestep_start:timestep_end,3);
    pitch2_measured = out(timestep_start:timestep_end,5);
    force_x0 = out(timestep_start:timestep_end,7);
    force_y0 = out(timestep_start:timestep_end,8);
    force_z0 = out(timestep_start:timestep_end,9);
    torque_x0 = out(timestep_start:timestep_end,10);
    torque_y0 = out(timestep_start:timestep_end,11);
    torque_z0 = out(timestep_start:timestep_end,12);
    force_D = force_y0.*cos(pitch2_measured) - force_x0.*sin(pitch2_measured);
    force_L = force_x0.*cos(pitch2_measured) + force_y0.*sin(pitch2_measured);
    inertialload_y = out(timestep_start:timestep_end,23);
    flowspeed_measured = abs(out(timestep_start:timestep_end,13));

    movmean_points = 100;
    heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),movmean_points);
    heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),movmean_points);

end