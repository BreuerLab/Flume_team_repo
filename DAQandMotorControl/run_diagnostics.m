function [diagnostics] = run_diagnostics(Prof_out_angle,dat,out,foiltype)
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
    
    %% Force PSD
    
    [pxx,w] = pwelch(kin.h2_meas);
    plot(w,10*log10(pxx))

    % cross correlation to check for time delay
    % flow rate — check that its constant
    % check for any force drift/misalignment
    % verify that amplitudes are expected (not too high or low)
    % compare commanded motion with measured in general ^^
    % 
end