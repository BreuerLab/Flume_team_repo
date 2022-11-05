%% Blockage correction - Barnsley and Wellicome
    
% pitch -- [rad]
% drag_coeff - non-dim
% heave_amp -- [m]
% U -- [m/s]
% foil -- foil characteristics
% depth -- flume depth [m]
% Eff -- total efficiency
% Fr -- Froude number

function [beta, U_prime, U_wake, U_wake_prime, Eff_prime, CD_norm, CD_norm_p] = blockage_barn_well(pitch_prof, drag_coeff, heave_amp, Eff, U, Fr, foil, depth)
    
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
    
    beta4_alpha4 = 1.5;
    
    err = 1;
%     tic;
    while err > 1e-7 % decreased this
        
        alpha2_alpha4 = (-1+sqrt(1+beta*(beta4_alpha4^2 - 1)))/(beta*(beta4_alpha4 - 1));
        alpha4_1 = beta4_alpha4 - beta*(alpha2_alpha4)*(beta4_alpha4 - 1);
        beta4_alpha4_2 = sqrt(CD_norm*alpha4_1^2 + 1);

        err = abs(beta4_alpha4 - beta4_alpha4_2);

        beta4_alpha4 = beta4_alpha4_2;
    
    end
%     toc;
    alpha2 = alpha2_alpha4*(beta4_alpha4 - beta*(alpha2_alpha4)*(beta4_alpha4 - 1))^-1;
    
    U_prime = U*(alpha2^2 + CD_norm/4)/alpha2;

    U_wake = (alpha2/alpha2_alpha4)*U; % velocity behind the leading foil
    
    U_wake_prime = (2*alpha2*(U/U_prime) - 1)*U_prime; % technically a correct calculation, but probably not useful for real-world applications
    
    CD_norm_p = CD_norm*((U/U_prime)^2); % corrected normalized drag coefficient

    Eff_prime = Eff*((U/U_prime)^3); % corrected efficiency
    
end