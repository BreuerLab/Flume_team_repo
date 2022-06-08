% This script will run analysis on the data that is the folder/ specified by the variable "filename"

% Load data
filename = 'Data\20220526_foilandvib2\foilandvib_pitch=0deg,f2='; name2='Hz,A2=';
% filename = ['foilandvib_pitch=0deg,f2='];
    load([filename,num2str(0.3),name2,num2str(0),'cm.mat']) %,constantpitch=-90deg

% singletrial_analysis = 1;
% fstarvector = 0.14;
% Astarvector = 1.1;
% ftrials = 1; Atrials = 1;

manytrial_analysis = 1;
fstarvector = (0.06:0.02:0.24);
Astarvector = (0:0.1:1.1);
ftrials = length(fstarvector); Atrials = length(Astarvector);

% Parameters
mass_disp = pi*(chord/2)^2*span*1000;
% transientcycs = 0; % Number of cycles for wind-up and wind-down

    fvector = fstarvector*U/chord;
    Avector = Astarvector*chord*100;

% Define size of each variable
f_star_commanded = nan(ftrials,Atrials);
A_star_commanded = nan(ftrials,Atrials);
A_star_measured = nan(ftrials,Atrials);
flowspeed_measured_mean = nan(ftrials,Atrials);
power_mean = nan(ftrials,Atrials);
powercoef_mean = nan(ftrials,Atrials);
power_scale = nan(ftrials,Atrials);
flowspeed_measured_p2p = nan(ftrials,Atrials);
f_force_filtered_dom = nan(ftrials,Atrials);
force_scale = nan(ftrials,Atrials);
delay = nan(ftrials,Atrials);
num_forcespecpeaks = nan(ftrials,Atrials);

% Loop through trials with different flow speed U
for ftrial = 1:ftrials

        timesteps_persubtrial = round((num_cyc/freq)/T);
        timesteps_subtrialcropped = round((num_cyc/freq)/T);
        transient_timesteps = round((transientcycs/freq)/T); % duration in seconds of transient heaving to crop off at beginning and end
    
    % Loop through subtrials with different heave amplitude A
for Atrial = 1:Atrials

    trialname = [filename,num2str(fvector(ftrial),3), ... 
        name2,num2str(Avector(Atrial),3),'cm.mat'];
    try
        load(trialname)
    catch
        disp(['Failed to load ',trialname])
    end


    % Extract measured quantities
    [time_star,heave_commanded,heave_measured,heave_star_measured,pitch_measured,force_D,force_L,inertialload_y,...
    flowspeed_measured,heave_velo,heave_accel] = extract_measurements(transientcycs,freq,T,Prof_out_angle,out);

    flowspeed_measured_mean(ftrial,Atrial) = mean(flowspeed_measured);
    Ustar_measured_mean(ftrial,Atrial) = flowspeed_measured_mean(ftrial,Atrial)/(chord*freq);

    f_star_commanded(ftrial,Atrial) = chord*freq/flowspeed_measured_mean(ftrial,Atrial);
    A_star_commanded(ftrial,Atrial) = (max(heave_commanded)-min(heave_commanded))/(2*chord);
    A_star_measured(ftrial,Atrial) = (max(heave_measured)-min(heave_measured))/(2*chord);

% Filter force data
    force_y_corrected = force_L+inertialload_y;
    [b,a] = butter(6,10*freq*(2*T),'low'); % butterworth filter 6th order with cut-off frequency at 10*freq
    force_y_corrected_filtered = filtfilt(b,a,squeeze(force_y_corrected)); 
    force_scale(ftrial,Atrial) = 0.5*1000*chord*span*flowspeed_measured_mean(ftrial,Atrial)^2;
    liftcoef = force_y_corrected_filtered/force_scale(ftrial,Atrial);

% Calculate power
    
    power_scale(ftrial,Atrial) = 0.5*1000*chord*span*flowspeed_measured_mean(ftrial,Atrial)^3;
    power_fluid = force_y_corrected_filtered .*heave_velo;
    power_damping = 0*heave_velo.^2;
    power_total = power_fluid + power_damping;
    power_mean(ftrial,Atrial) = mean(power_fluid);
    powercoef = power_fluid/power_scale(ftrial,Atrial);
    powercoef_mean(ftrial,Atrial) = power_mean(ftrial,Atrial)/power_scale(ftrial,Atrial);

% % heave spectrum
    duration = max(time_star/freq)/T;
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(0*window_duration/2);
    heave_star_hilbert =  hilbert(heave_star_measured);
    [heave_powerspec,f_heave] = pwelch(heave_star_hilbert,window_duration,overlap,[],1/T);
    [max_power,max_index] = max(10*log10(heave_powerspec));
    f_heave_dom = f_heave(max_index);


 % force corrected+filtered spectrum using Welch's method
    duration = round(max(time_star/freq)/T);
    window_duration = round(duration/6); % size of hilbert windows measured in samples
    overlap = round(0*window_duration);
    force_hilbert =  hilbert(liftcoef);
    [force_powerspec,f_force] = pwelch(force_hilbert,window_duration,overlap,[],1/T);
    [max_force_power,max_force_index] = max(10*log10(force_powerspec));
    f_force_dom = f_force(max_force_index);  
%     spacing = (f_force(2)-f_force(1))/f;
%     findpeaks(10*log10(force_powerspec),1/spacing,'MinPeakHeight',max_power-20)

% %     % Find phase delay between heave and force
% %     maxlag = round(1/(2*freq*T));  % Maximum lag to calculate xcorr (in timesteps)
% %     [corrs,lags] = xcorr(heave_measured,force_y_corrected_filtered,maxlag);
% % %     delay = freq*T*360*finddelay(heave_measured,force_y_corrected_filtered,maxlag)
% %     hold on
% %     plot(lags*freq*T*360,corrs)
% %     [forcedelay_peak,forcedelay_peaklocs] = findpeaks(corrs,1);
% %     forcedelay_deg = forcedelay_peaklocs*T*freq*360-180; % Delay in degrees
% %     stem(forcedelay_deg,forcedelay_peak);
% %     xlabel('Lag (degrees)'); ylabel('Correlation')
% %     xlim([-180 180])
% 
% %     % Hilbert transform for instantaneous freq
% %     window = floor(duration/10);
% %     noverlap = floor(window/4);
% %     fftpoints = window;
% %     spectrogram(liftcoef,window,noverlap,fftpoints,1/T,'yaxis');
% %     caxis([-30 10])
% %     ylim([0 10*freq])
% %     xlim([0 max(times)])
% 
%    % Find peaks
% spacing = (f_force(2)-f_force(1))/freq;
%     [forcespec_peakpowers,forcespec_peaklocs] = findpeaks(10*log10(force_powerspec),1/spacing, ...
%     'MinPeakHeight',max_force_power-20,'MinPeakDistance',0.5);
%     num_forcespecpeaks = size(forcespec_peaklocs,1);
%     
% % Plot power spectrum
% St = 0.21; % estimated Strouhal number
% f_vortex = St*(flowspeed_measured_mean(ftrial,Atrial)/chord);
% plot_PrescribedMotionPowerSpectrum_MATLABin
% %(freq,f_vortex,f_heave_dom,f_heave,heave_powerspec,f_force,force_powerspec)

% Plot force and motion
    if singletrial_analysis==1
        close all
    plot_PrescribedMotionForceAndVelocity_MATLABin(time_star,heave_measured,heave_velo,force_y_corrected_filtered,power_fluid);
    drawnow
    end
end

end

