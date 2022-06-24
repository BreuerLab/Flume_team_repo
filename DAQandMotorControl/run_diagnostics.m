function [diagnostics] = run_diagnostics(Prof_out_angle,dat,out,foiltype,fs)
    %% Check Encoder Data
    
    function [wave_mean] = verify_sinusoidal(motion_data,name)
        wave_mean = mean(motion_data);
        
        if (abs(wave_mean) > 0.1)
            warning(name + " encoder is outside margin.");
        end
    end

    diagnostics.p2_mean = verify_sinusoidal(out(:,3),"Gromit pitch");
    diagnostics.h2_mean = verify_sinusoidal(out(:,4),"Gromit heave");
    diagnostics.p3_mean = verify_sinusoidal(out(:,5),"Wallace pitch");
    diagnostics.h3_mean = verify_sinusoidal(out(:,6),"Wallace heave");

    %% Commanded Motion v. Measured

    % cross correlation
    function [lag] = calculate_lag(comm,measured)
        [c,lags] = xcorr(comm,measured);
        [max_c,index] = max(c);
        lag = lags(index) / fs; % fs = sample rate
    end
    
    diagnostics.p2_lag = calculate_lag(Prof_out_angle(:,3),out(:,3));
    diagnostics.h2_lag = calculate_lag(Prof_out_angle(:,4),out(:,4));
    diagnostics.p3_lag = calculate_lag(Prof_out_angle(:,5),out(:,5));
    diagnostics.h3_lag = calculate_lag(Prof_out_angle(:,6),out(:,6));
    
    % check amplitudes
    % standard deviation of paths

    %% Flow Rate

    %% Force Drift/Misalignment
    
    %% Force PSD
    
    [pxx,w] = pwelch(out(:,4));
    plot(w,10*log10(pxx))

    % cross correlation to check for time delay
    % flow rate — check that its constant
    % check for any force drift/misalignment
    % verify that amplitudes are expected (not too high or low)
    % compare commanded motion with measured in general ^^
    % 
end