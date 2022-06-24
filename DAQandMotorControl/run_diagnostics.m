function [diagnostics] = run_diagnostics(Prof_out_angle,dat,out,foiltype,fs)
    %% Locally Label Variables
    % NOTE: converting comm to radians
    p2n = "Gromit pitch";
    p2c = deg2rad(Prof_out_angle(:,3)); % gromit pitch comm
    p2m = out(:,3); % gromit pitch measured
    h2n = "Gromit heave";
    h2c = Prof_out_angle(:,4); % gromit heave comm
    h2m = out(:,4); % gromit heave measured
    p3n = "Wallace pitch";
    p3c = deg2rad(Prof_out_angle(:,5)); % wallace pitch comm
    p3m = out(:,5); % wallace pitch measured
    h3n = "Wallace heave";
    h3c = Prof_out_angle(:,6); % wallace heave comm
    h3m = out(:,6); % wallace heave measured

    vx = out(:,13); % vectrino x
    vy = out(:,14); % vectrino y
    vz1 = out(:,15); % vectrino z1
    vz2 = out(:,16); % vectrino z2

    %% Check Encoder Data
    
    function [wave_mean] = verify_sinusoidal(motion_data,name)
        wave_mean = mean(motion_data);
        
        if (abs(wave_mean) > 0.1)
            warning(name + " encoder is offset or drifting.");
        end
    end

    diagnostics.p2_mean = verify_sinusoidal(p2m,p2n);
    diagnostics.h2_mean = verify_sinusoidal(h2m,h2n);
    diagnostics.p3_mean = verify_sinusoidal(p3m,p3n);
    diagnostics.h3_mean = verify_sinusoidal(h3m,h3n);

    %% Commanded Motion v. Measured

    % CHECK PHASE PLOTS

    % cross correlation
    function [lag] = calculate_lag(comm,measured)
        [c,lags] = xcorr(comm,measured);
        [~,index] = max(c);
        lag = lags(index) / fs; % fs = sample rate
    end
    
    diagnostics.p2_lag = calculate_lag(p2c,p2m);
    diagnostics.h2_lag = calculate_lag(h2c,h2m);
    diagnostics.p3_lag = calculate_lag(p3c,p3m);
    diagnostics.h3_lag = calculate_lag(h3c,h3m);
    
    % check amplitudes
    function [amplitude_diff] = calculate_amplitude_diff(comm,measured,name)
        c_local_maxima = findpeaks(comm);
        [b,a] = butter(5,0.01);
        smooth_measured = filter(b,a,measured);
        [~,pos] = findpeaks(smooth_measured);
        m_local_maxima = measured(pos);

        if (length(c_local_maxima) == length(m_local_maxima))
            diffs = c_local_maxima - m_local_maxima;
            amplitude_diff = mean(diffs);
        else
            warning(name + " encoder is producing unexpected wave peaks.");
            amplitude_diff = -1;
        end
    end

    diagnostics.p2_ampl_delta = calculate_amplitude_diff(p2c,p2m,p2n);
    diagnostics.h2_ampl_delta = calculate_amplitude_diff(h2c,h2m,h2n);
    diagnostics.p3_ampl_delta = calculate_amplitude_diff(p3c,p3m,p3n);
    diagnostics.h3_ampl_delta = calculate_amplitude_diff(h3c,h3m,h3n);
    
    % difference between paths
    function [diff] = normalized_path_diff(comm,measured,name)
        diff = normalize(comm,"range") - normalize(measured,"range");
        if (max(abs(diff)) > 0.05)
            warning(name + " encoder value is over 5% off from commanded.")
        end
    end

    diagnostics.p2_diff = normalized_path_diff(p2c,p2m,p2n);
    diagnostics.h2_diff = normalized_path_diff(h2c,h2m,h2n);
    diagnostics.p3_diff = normalized_path_diff(p3c,p3m,p3n);
    diagnostics.h3_diff = normalized_path_diff(h3c,h3m,h3n);

    %% Flow Rate
    function [vectrino_std] = calculate_vectrino_std(data)
        vectrino_std = std(data);
    end

    diagnostics.vx_std = calculate_vectrino_std(vx);
    diagnostics.vy_std = calculate_vectrino_std(vy);
    diagnostics.vz1_std = calculate_vectrino_std(vz1);
    diagnostics.vz2_std = calculate_vectrino_std(vz2);

    %% Force Drift/Misalignment

    % standard deviation margin for cycle averaging
    
    %% Force PSD
    
    [pxx,w] = pwelch(out(:,4));
    %plot(w,10*log10(pxx))

    % cross correlation to check for time delay
    % flow rate — check that its constant
    % check for any force drift/misalignment
    % verify that amplitudes are expected (not too high or low)
    % compare commanded motion with measured in general ^^
    % 
end