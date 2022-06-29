function [diagnostics] = run_diagnostics(Prof_out_angle,out,fs,EP)
    addpath(genpath("Libraries"))

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

    % leading foil
    lfy = out(:,17); % leading f_n
    lfy_n = "Leading normal force";
    lfx = out(:,18); % leading f_t
    lfx_n = "Leading tangential force";
    ltz = out(:,22); % leading t_z
    ltz_n = "Leading spanwise torque";

    % trailing foil
    tfy = out(:,07); % trailing f_n
    tfy_n = "Trailing normal force";
    tfx = out(:,08); % trailing f_t
    tfx_n = "Trailing tangential force";
    ttz = out(:,12); % trailing t_z
    ttz_n = "Trailing spanwise torque";

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

    % psd
    function [] = plot_tiled_psd(data,w,psd)
        nexttile
        plot(data);
        nexttile
        plot(w/EP.freq,psd);
        xlim([0,25])
        nexttile
        plot(w/EP.freq,psd);
        xlim([25,120/EP.freq])
    end

    function [w,psd] = calculate_psd(data)
        [psd,w] = pwelch(data,[],[],[],fs);
        psd = 10*log10(psd);

        % uncomment for visualization
        plot_tiled_psd(data,w,psd);
    end

    function [idx] = find_nearest(w,value)
        [~,idx] = min(abs(w-value));
    end

    function [peaks,indices] = scan_frequencies(max_expected,psd_seg)
        minValue = min(psd_seg);
        psd_seg(psd_seg < max_expected - 20) = minValue;
        [peaks,indices] = findpeaks(psd_seg);
    end

    function [pois] = evaluate_force(force,name)
        [w,psd] = calculate_psd(force);
        w_norm = w/EP.freq;
        psd_seg = psd(find_nearest(w_norm,25):find_nearest(w_norm,120/EP.freq));
        [pks,indcs] = scan_frequencies(psd(find_nearest(w_norm,1)),psd_seg);

        if (~isempty(pks))
            warning(name + " measurement includes " + length(pks) + " potentially unexpected vibration frequencies.");
        end

        freqs = w(indcs + find_nearest(w_norm,25) - 1);
        pois = [freqs,pks]; % points of interest
    end

    % uncomment for visualization
    tiledlayout(6,3);
    diagnostics.lfy_pois = evaluate_force(lfy,lfy_n);
    diagnostics.lfx_pois = evaluate_force(lfx,lfx_n);
    diagnostics.ltz_pois = evaluate_force(ltz,ltz_n);
    diagnostics.tfx_pois = evaluate_force(tfy,tfy_n);
    diagnostics.tfy_pois = evaluate_force(tfx,tfx_n);
    diagnostics.ttz_pois = evaluate_force(ttz,ttz_n);

    % standard deviation margin for cycle averaging
    function [norm_dev] = norm_cycle_std(cycle,data,name)
        [toverT,pitch_cycle,data_cycle] = cycle_avg_data(cycle,data);

        mean_line = mean(data_cycle(4:30,:));
        dev = std(data_cycle(4:30,:));
        
        [b,a] = butter(5,20*EP.freq/(fs/2));
        f_mean_line = filtfilt(b,a,mean_line);
        f_dev = filtfilt(b,a,dev);

        [norm_mean,c,s] = normalize(f_mean_line,'range',[-1 1]);
        norm_dev = normalize(f_dev,'center',c,'scale',s);

        if (max(norm_dev) > 0.5)
            amt_over = length(find(norm_dev > 0.5));
            warning(name + "'s standard deviation is " + amt_over/length(norm_dev)*100 ... ...
                + "% over 1/2 the amplitude of the average cycle.")
        end

        nexttile
        hold on

        plot(toverT,mean(pitch_cycle(4:30,:)))

        above = norm_mean + norm_dev;
        below = norm_mean - norm_dev;
        t2 = [toverT,fliplr(toverT)];
        between = [below,fliplr(above)];
        fill(t2,between,[0 0.2235 0.3705]);

        plot(toverT,norm_mean,"LineWidth",2);
        %plot(pitch_cycle(4:30,:),data_cycle(4:30,:));
    end

    figure
    tiledlayout(2,3);
    diagnostics.lfy_normdev = norm_cycle_std(p2c,lfy,lfy_n);
    diagnostics.lfx_normdev = norm_cycle_std(p2c,lfx,lfx_n);
    diagnostics.ltz_normdev = norm_cycle_std(p2c,ltz,ltz_n);
    diagnostics.tfy_normdev = norm_cycle_std(p3c,tfy,tfy_n);
    diagnostics.tfx_normdev = norm_cycle_std(p3c,tfx,tfx_n);
    diagnostics.ttz_normdev = norm_cycle_std(p3c,ttz,ttz_n);
end