%% Blockage correction (Houlsby et al, Ross et al, Ribeiro et al)
    
% pitch -- [rad]
% drag_coeff - non-dim
% heave_amp -- [m]
% U -- [m/s]
% foil -- foil characteristics
% depth -- flume depth [m]
% Eff -- total efficiency
% Fr -- Froude number

function [beta, U_prime, Eff_prime, CD_norm, CD_norm_p] = blockage_houlsby(pitch_prof, drag_coeff, heave_amp, Eff, U, Fr, foil, depth)

    width = 0.8; % flume width [m]
    g = 9.80665; % acceleration of gravity [m/s^2]
    H = heave_amp; % heaving amplitude [m]
    CD = drag_coeff;
    
    % Normalization of average CD

%     [~, ~, CD_cyc] = cycle_avg_data(pitch_prof,CD); % phase-averaged drag coeff
%     CD_avg = mean(mean(CD_cyc)); % mean drag coefficient per cycle
    CD_avg = mean(CD);

    CD_norm = CD_avg*(foil.chord/(2*H)); % normalized drag coeff to be used in calculations

    beta = (2*H*foil.span)/(width*depth); % blockage ratio

    % Iteration to find u1

    u2 = U; % initial guess for u2
    err = 1; % initializing the error
    i = 0;

    while err > 0.0001

        u1_a = (Fr^2*u2^4 - (4 + 2*Fr^2)*U^2*u2^2 + 8*U^3*u2 - 4*U^4 + 4*beta*CD_norm*U^4 + Fr^2*U^4)/...
            (-4*Fr^2*u2^3 + (4*Fr^2 + 8)*U^2*u2 - 8*U^3);

        u1_b = sqrt(u2^2 - CD_norm*U^2);

        err = abs(u1_a - u1_b);

        u2 = u2 + 0.000001; % decreased the step size
        i = i + 1;

    end

    % Parameter correction

    u1 = u1_a;

    ud = (u1*(u2-U)*(2*g*depth - u2^2 - u2*U))/(2*beta*g*depth*(u2 - u1));

    U_prime = U*((ud/U)^2 + CD_norm*0.25)/(ud/U); % corrected free-stream velocity

    CD_norm_p = CD_norm*((U/U_prime)^2); % corrected normalized drag coefficient

    Eff_prime = Eff*((U/U_prime)^3); % corrected efficiency
    
end