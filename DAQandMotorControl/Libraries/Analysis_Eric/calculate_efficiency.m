%% Calculate efficiency

% eff --> efficiency calculations

function   [eff, par] = calculate_efficiency(kin, res, par, foil, EP)
    
    rho = 1000; % density of water
    
    % Power calculation
    
    eff.PwrH2 = res.Lift2.*kin.h2_vel; % leading power due to heaving
    eff.PwrP2 = res.torque2_z0.*kin.p2_vel; % leading power due to pitching

    eff.PwrH3 = res.Lift3.*kin.h3_vel; % trailing power due to heaving
    eff.PwrP3 = res.torque3_z0.*kin.p3_vel; % trailing power due to pitching
    
    % Lift, Drag, Torque and Power coefficients
    
    eff.CL2 = res.Lift2/(0.5*rho*par.U^2*foil.span*foil.chord); % leading lift coeff
    eff.CD2 = res.Drag2/(0.5*rho*par.U^2*foil.span*foil.chord); % leading drag coeff
    eff.CM2 = res.torque2_z0/(0.5*rho*par.U^2*foil.span^2*foil.chord); % leading moment coeff
    eff.CPH2 = eff.PwrH2/(0.5*rho*par.U^3*foil.span*foil.chord); % leading heaving power coeff
    eff.CPP2 = eff.PwrP2/(0.5*rho*par.U^3*foil.span*foil.chord); % leading pitching power coeff
    
    eff.CL3 = res.Lift3/(0.5*rho*par.U^2*foil.span*foil.chord); % leading lift coeff
    eff.CD3 = res.Drag3/(0.5*rho*par.U^2*foil.span*foil.chord); % leading drag coeff
    eff.CM3 = res.torque3_z0/(0.5*rho*par.U^2*foil.span^2*foil.chord); % leading moment coeff
    eff.CPH3 = eff.PwrH3/(0.5*rho*par.U^3*foil.span*foil.chord); % leading heaving power coeff
    eff.CPP3 = eff.PwrP3/(0.5*rho*par.U^3*foil.span*foil.chord); % leading pitching power coeff
    
    % Efficiency calculation
    
    [~, ~, PwrP2_cyc] = phase_avg_data(kin.p2_meas, eff.PwrP2); % averaged data per cycle respect to pitching
    [~, ~, PwrH2_cyc] = phase_avg_data(kin.h2_meas, eff.PwrH2);

    [~, ~, PwrP3_cyc] = phase_avg_data(kin.p3_meas, eff.PwrP3); % averaged data per cycle respect to pitching
    [~, ~, PwrH3_cyc] = phase_avg_data(kin.h3_meas, eff.PwrH3);

    % Yp = max(2*(Prof_out_angle(:,6) + 0.5*foil.chord*sin(deg2rad(Prof_out_angle(:,5))))); % maximum swept area

    yp1_2 = kin.h2_meas + foil.chord*0.5*sin(kin.p2_meas);
    yp2_2 = kin.h2_meas - foil.chord*0.5*sin(kin.p2_meas);
    eff.Yp_2 = 2*max(max(yp1_2),max(yp2_2)); % maximum distance travelled by the leading edge of the leading foil
    
    yp1_3 = kin.h3_meas + foil.chord*0.5*sin(kin.p3_meas);
    yp2_3 = kin.h3_meas - foil.chord*0.5*sin(kin.p3_meas);
    eff.Yp_3 = 2*max(max(yp1_3),max(yp2_3)); % maximum distance travelled by the leading edge of the trailing foil

    eff.Yp = max(eff.Yp_2, eff.Yp_3); % assuming both foils have the same chord

    eff.EffP_2 = mean(mean(PwrP2_cyc))/(0.5*rho*par.U^3*eff.Yp*foil.span); % single foil uncorrected pitching efficiency
    eff.EffH_2 = mean(mean(PwrH2_cyc))/(0.5*rho*par.U^3*eff.Yp*foil.span); % single foil uncorrected heaving efficiency
    
    eff.EffP_3 = mean(mean(PwrP3_cyc))/(0.5*rho*par.U^3*eff.Yp*foil.span); % single foil uncorrected pitching efficiency    
    eff.EffH_3 = mean(mean(PwrH3_cyc))/(0.5*rho*par.U^3*eff.Yp*foil.span); % single foil uncorrected heaving efficiency

    eff.Eff_2 = eff.EffH_2 + eff.EffP_2; % total leading foil efficiency
    eff.Eff_3 = eff.EffH_3 + eff.EffP_3; % total trailing foil efficiency
    
    %% Blockage correction (Houlsby et al, Ross et al, Ribeiro et al)
    
    [eff.beta_2, eff.U_prime_2, eff.Eff_2_prime] = ...
        blockage_houlsby(kin.p2_meas, eff.CD2, EP.H2, eff.Eff_2, par.U, par.Fr, foil, EP.flume_depth);
    
    [eff.beta_3, eff.U_prime_3, eff.Eff_3_prime] = ...
        blockage_houlsby(kin.p3_meas, eff.CD3, EP.H3, eff.Eff_3, par.U, par.Fr, foil, EP.flume_depth);
    
end
