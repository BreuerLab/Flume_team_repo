%% Calculate Forces

function res = calculate_forces(par, kin, out)
    
    % Extracting relevant variables from 'kin' and 'par'

    tsteps_exp = kin.tsteps_exp; tstep_start = kin.tstep_start; tstep_end = kin.tstep_end; num_cyc = par.num_cyc;
    h2_comm = kin.h2_comm; h3_comm = kin.h3_comm; % commanded leading (or trailing) heave (for the cycle-avg power)
    p2_comm = kin.p2_comm; p3_comm = kin.p3_comm; % commanded leading (or trailing) pitch (for the blockage correction)
    p2_meas = kin.p2_meas; p3_meas = kin.p3_meas; % measured pitch and heave
    h2_meas = kin.h2_meas; h3_meas = kin.h3_meas;
    p2_vel = kin.p2_vel; p3_vel = kin.p3_vel; % pitch and heave velocities
    h2_vel = kin.h2_vel; h3_vel = kin.h3_vel;
    h2_acc = kin.h2_acc; h3_acc = kin.h3_acc; % pitch and heave accelerations
    U = par.U; % flow velocity
    
    [foil, ~, ~] = foils_database(par.foiltype); % foil characteristics
    rho = 1000; % density of water
    
    % Measured forces
    
    f2n = out(tstep_start:tstep_end,17); % leading normal force [N]
    f2t = out(tstep_start:tstep_end,18); % leading tangential force [N]
    tq2 = -out(tstep_start:tstep_end,22); % leading pitch axis torque [N*m] NOTE: negative due to reference transformation
    
    f3n = out(tstep_start:tstep_end,7); % trailing normal force [N]
    f3t = out(tstep_start:tstep_end,8); % trailing tangential force [N]
    tq3 = -out(tstep_start:tstep_end,12); % trailing pitch axis torque [N*m] NOTE: negative due to reference transformation

    inertialload_y = out(tstep_start:tstep_end,23); % inertial load in the y-direction (flume coordinates, NOT sensor coordinates)

    % Calcualted forces
    
    Drag2 = f2t.*cos(p2_meas) - f2n.*sin(p2_meas); % leading Drag force [N]
    Lift2 = f2n.*cos(p2_meas) + f2t.*sin(p2_meas) - (foil.mass1+0.6)*h2_acc; % leading Lift force [N]
    
    Drag3 = f3t.*cos(p3_meas) - f3n.*sin(p3_meas); % trailing Drag force [N]
    Lift3 = f3n.*cos(p3_meas) + f3t.*sin(p3_meas) - (foil.mass2+0.6)*h3_acc; % trailing Lift force [N]

    % Power calculation
    
    PwrH2 = Lift2.*h2_vel; % leading power due to heaving
    PwrP2 = tq2.*p2_vel; % leading power due to pitching

    PwrH3 = Lift3.*h3_vel; % trailing power due to heaving
    PwrP3 = tq3.*p3_vel; % trailing power due to pitching
    
    Pwr2_tot = PwrH2 + PwrP2; % total leading power extracted
    Pwr3_tot = PwrH3 + PwrP3; % total trailing power extracted
    
    if mean(abs(h2_comm)) == 0
        heave_profile = h3_comm; % profile used to determine the frequency, named generally to avoid confusion in subsequent calculations
    else
        heave_profile = h2_comm; % if trailing pitch is zero, use leading pitch to determine frequency
    end
    
    min_peak = max(abs(heave_profile))*0.99; % determines the minimum value that the peaks should have
    [~,locs] = findpeaks(heave_profile,'MinPeakHeight',min_peak);
    locs2 = locs - locs(1) + 1; % location of the starting timestep of each cycle (starting from a heave amplitude of 0)
    
    cycle_steps = NaN(1,num_cyc-1); % steps in each individual cycle
    for i = 1:length(locs2)-1
        cycle_steps(i) = locs(i+1) - locs(i);
    end
    tstep_cyc = round(mean(cycle_steps)); % this is almost exactly the same variable as 'tsteps_cyc' from "kin", but obtained differently
    
    PwrH2_cyc = NaN(tstep_cyc, num_cyc);
    PwrP2_cyc = NaN(tstep_cyc, num_cyc);
    PwrH3_cyc = NaN(tstep_cyc, num_cyc);
    PwrP3_cyc = NaN(tstep_cyc, num_cyc);
    
    for i = 1:num_cyc % cycle-averaging the power respect to the heaving
        PwrH2_cyc(:,i) = PwrH2(locs2(i):locs2(i)+tstep_cyc-1);
        PwrP2_cyc(:,i) = PwrP2(locs2(i):locs2(i)+tstep_cyc-1);
        PwrH3_cyc(:,i) = PwrH3(locs2(i):locs2(i)+tstep_cyc-1);
        PwrP3_cyc(:,i) = PwrP3(locs2(i):locs2(i)+tstep_cyc-1);
    end
    
%     figure
%     plot(mean(PwrH2_cyc,2)); hold on;
%     plot(mean(PwrH3_cyc,2));
    
    % Swept area
    
    yp1_2 = h2_meas + foil.chord*0.5*sin(p2_meas);
    yp2_2 = h2_meas - foil.chord*0.5*sin(p2_meas);
    Yp_2 = 2*max(max(yp1_2),max(yp2_2)); % maximum distance travelled by the leading edge of the leading foil
    
    yp1_3 = kin.h3_meas + foil.chord*0.5*sin(p3_meas);
    yp2_3 = kin.h3_meas - foil.chord*0.5*sin(p3_meas);
    Yp_3 = 2*max(max(yp1_3),max(yp2_3)); % maximum distance travelled by the leading edge of the trailing foil

    Yp = max(Yp_2, Yp_3); % assuming both foils have the same chord
    
    % Efficiency calculation
    
    EffP_2 = mean(mean(PwrP2_cyc,2))/(0.5*rho*U^3*Yp_2*foil.span); % single foil uncorrected pitching efficiency
    EffH_2 = mean(mean(PwrH2_cyc,2))/(0.5*rho*U^3*Yp_2*foil.span); % single foil uncorrected heaving efficiency
    
    EffP_3 = mean(mean(PwrP3_cyc,2))/(0.5*rho*U^3*Yp_3*foil.span); % single foil uncorrected pitching efficiency    
    EffH_3 = mean(mean(PwrH3_cyc,2))/(0.5*rho*U^3*Yp_3*foil.span); % single foil uncorrected heaving efficiency

    Eff_2 = EffH_2 + EffP_2; % total leading foil efficiency
    Eff_3 = EffH_3 + EffP_3; % total trailing foil efficiency
    
    Eff_sys = (mean(mean(PwrP2_cyc+PwrH2_cyc,2))+mean(mean(PwrP3_cyc+PwrH3_cyc,2)))/(0.5*rho*U^3*Yp*foil.span); % system efficiency
    
    % Non-dimensional forces
    
    CL2 = Lift2/(0.5*rho*U^2*foil.span*foil.chord); % leading lift coeff
    CD2 = Drag2/(0.5*rho*U^2*foil.span*foil.chord); % leading drag coeff
    CM2 = tq2/(0.5*rho*U^2*foil.span^2*foil.chord); % leading moment coeff
    CPH2 = PwrH2/(0.5*rho*U^3*foil.span*foil.chord); % leading heaving power coeff
    CPP2 = PwrP2/(0.5*rho*U^3*foil.span*foil.chord); % leading pitching power coeff
    
    CL3 = Lift3/(0.5*rho*U^2*foil.span*foil.chord); % leading lift coeff
    CD3 = Drag3/(0.5*rho*U^2*foil.span*foil.chord); % leading drag coeff
    CM3 = tq3/(0.5*rho*U^2*foil.span^2*foil.chord); % leading moment coeff
    CPH3 = PwrH3/(0.5*rho*U^3*foil.span*foil.chord); % leading heaving power coeff
    CPP3 = PwrP3/(0.5*rho*U^3*foil.span*foil.chord); % leading pitching power coeff
    
    % Blockage correction
    
    [beta2, U_2prime, Eff_2prime] = blockage_houlsby(p2_comm, CD2, par.H2, Eff_2, U, par.Fr, foil, par.flume_height);
    [beta3, U_3prime, Eff_3prime] = blockage_houlsby(p3_comm, CD3, par.H3, Eff_3, U, par.Fr, foil, par.flume_height);
    
    %% Store variables
    
    res.f2n = f2n; % normal force
    res.f2t = f2t; % tangential force
    res.tq2 = tq2; % pitching torque
    res.f3n = f3n;
    res.f3t = f3t;
    res.tq3 = tq3;
    res.yinertia = inertialload_y; % normal inertial force
    
    res.Lift2 = Lift2; % leading lift
    res.Drag2 = Drag2; % leading drag
    res.Lift3 = Lift3; % trailing lift
    res.Drag3 = Drag3; % trailing drag
    
    res.PwrH2 = PwrH2; % leading heaving power
    res.PwrH3 = PwrH3; % trailing heaving power
    res.PwrP2 = PwrP2; % leading pitching power
    res.PwrP3 = PwrP3; % trailing pitching power
    
    res.Pwr2_tot = Pwr2_tot; % total leading power
    res.Pwr3_tot = Pwr3_tot; % total trailing power
    
    res.Yp_2 = Yp_2; % leading swept distance
    res.Yp_3 = Yp_3; % trailing swept distance
    res.Yp = Yp; % maximum swept distance (from the two foils)
    
    res.EffP_2 = EffP_2; % leading pitching efficiency
    res.EffH_2 = EffH_2; % leading heaving efficiency
    res.EffP_3 = EffP_3; % trailing pitching efficiency
    res.EffH_3 = EffH_3; % trailing heaving efficiency
    
    res.Eff_2 = Eff_2; % total leading efficiency
    res.Eff_3 = Eff_3; % total trailing efficiency
    
    res.Eff_sys = Eff_sys; % system efficiency
    
    res.Eff_2prime = Eff_2prime; % corrected leading efficiency
    res.Eff_3prime = Eff_3prime; % corrected trailing efficiency
    res.beta2 = beta2; % leading blockage ratio
    res.beta3 = beta3; % trailing blockage ratio
    
    res.CL2 = CL2; % leading lift coeff
    res.CD2 = CD2; % leading drag coeff
    res.CM2 = CM2; % leading moment coeff
    res.CPH2 = CPH2; % leading heave power coeff
    res.CPP2 = CPP2; % leadinf pitch power coeff
    
    res.CL3 = CL3; % trailing lift coeff
    res.CD3 = CD3; % trailing drag coeff
    res.CM3 = CM3; % trailing moment coeff
    res.CPH3 = CPH3; % trailing heave power coeff
    res.CPP3 = CPP3; % trailing pitch power coeff
    
    res.U_2prime = U_2prime; % corrected flowspeed?
    res.U_3prime = U_3prime; % corrected flowspeed
    
end