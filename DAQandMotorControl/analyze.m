% This script will run analysis on the data that is the folder/ specified by the variable "filename"

% Load data
datadir = 'R:\ENG_Breuer_Shared\jnewbolt\DAQandMotorControl\Data\';
% trialdir = 'CircCylHighRes_03-Nov-2022_19_35_11\data\'; namepart1 = 'CircCyl_pitch=0deg,f='; namepart2='Hz,A=';
% Other data locations
% datadir = 'C:\Users\Joel\Documents\Brown Data\';
% datadir = 'D:\Experiments\1foil\'; 
% trialdir = 'Enter descriptive name_03-Nov-2022_17_27_4\data\'; namepart1 = 'Vib_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'FoilAndEllipticCyl_13-Nov-2022_22_31_38\data\';namepart1 = 'FoilAndEllipticCyl_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'EllipticCylHigherFreq_06-Nov-2022_12_49_19\data\';namepart1 = 'EllipticCyl_pitch=0deg,f='; namepart2='Hz,A=';
% trialdir = 'VibHigherFreq_04-Nov-2022_19_16_8\data\';namepart1 = 'Vib_pitch=0deg,f='; namepart2='Hz,A=';
trialdir = 'FoilAndVib_13-Nov-2022_14_11_29\data\';namepart1 = 'FoilAndVib_pitch=0deg,f='; namepart2='Hz,A=';

% Combine strings to form filename and load last trial to get some necessary variable values from the trial
filename = [datadir,trialdir,namepart1];
trialfiles = dir([datadir,trialdir]);
load([datadir,trialdir,trialfiles(4).name]); freq = freq2;

singletrial_analysis = 1;
manytrial_analysis = 1;
varyphase = 0;

if manytrial_analysis==1
fstarvector = (0.1:0.01:0.4);%(0.05:0.01:0.14);
Astarvector = (0.0:0.05:1.1);
ftrials = length(fstarvector); Atrials = length(Astarvector);
    if varyphase==1
    fstarvector = 0.09;
    phase12vector = (-180:36:144);
    ftrials = length(phase12vector); 
    end
elseif singletrial_analysis==1
fstarvector = 0.12; %0.13;
Astarvector = 0.8; %1.1;
ftrials = 1; Atrials = 1;
end

addpath(genpath("Libraries"));

% Parameters
%     U=0.3;
    fvector = fstarvector*U/thcknss;
    Avector = Astarvector*thcknss*100;

% %     % Test for Eric's foil
% %     fvector = (0.4328:0.0541:0.7574); %[0.5 0.6 0.7];
% %     Avector = 0; %[0 1 2 3];
% %     thcknss = 0.061; span = 6*thcknss;
% %     ftrials = length(fvector); Atrials = length(Avector);
%     fvector = fstarvector*U/chord;
%     Avector = Astarvector*chord*100;


%% Define size of each variable
f_star_commanded = nan(ftrials,Atrials);
A_star_commanded = nan(ftrials,Atrials);
A_star_measured = nan(ftrials,Atrials);
phase12 = nan(ftrials,Atrials);
flowspeed_measured_mean = nan(ftrials,Atrials);
flowspeed_measured_stdev = nan(ftrials,Atrials);
power_mean = nan(ftrials,Atrials);
powercoef_mean = nan(ftrials,Atrials);
power_scale = nan(ftrials,Atrials);
flowspeed_measured_p2p = nan(ftrials,Atrials);
f_force_filtered_dom = nan(ftrials,Atrials);
force_scale = nan(ftrials,Atrials);
delay = nan(ftrials,Atrials);
num_forcespecpeaks = nan(ftrials,Atrials);
liftcoef_alltrials = nan(ftrials,Atrials);
liftcoef_max = nan(ftrials,Atrials);
lifcoef_min = nan(ftrials,Atrials);
dragcoef_alltrials = nan(ftrials,Atrials);
dragcoef_max = nan(ftrials,Atrials);
dragcoef_min = nan(ftrials,Atrials);
inertialload_alltrials = nan(ftrials,Atrials);
force_inphasewvelo  = nan(ftrials,Atrials);
force_inphasewaccel  = nan(ftrials,Atrials);
force_inphasewvelo_star  = nan(ftrials,Atrials);
force_inphasewaccel_star  = nan(ftrials,Atrials);
energybalance = nan(ftrials,Atrials);
energyin_star =  nan(ftrials,Atrials);
energyout_star =  nan(ftrials,Atrials);
m_star_added_eff =  nan(ftrials,Atrials);

%%

% Loop through trials with different flow speed U
for ftrial = 1:ftrials
        T = 1/samplerate;
        timesteps_persubtrial = round((num_cyc/freq)/T);
        timesteps_subtrialcropped = round((num_cyc/freq)/T);
        transient_timesteps = round((transientcycs/freq)/T); % duration in seconds of transient heaving to crop off at beginning and end
    
    % Loop through subtrials with different heave amplitude A
for Atrial = 1:Atrials

    if varyphase==0
        trialname = [filename,num2str(fvector(ftrial),3),namepart2,num2str(Avector(Atrial),3),'cm.mat'];
    elseif varyphase==1
        trialname = [filename,num2str(fvector,3), ... 
            namepart2,num2str(Avector(Atrial),3), name3,num2str(phase12vector(ftrial)),'deg.mat'];
    end

    try
        load(trialname,'transientcycs','out','Prof_out_angle','freq2','phase2'); freq = freq2;
    catch
        disp(['Failed to load ',trialname])
        break
    end
%     freq2 = freq; % Added as workaround for cases with 2 different frequencies
    % Extract measured quantities
    [time_star,heave_commanded,heave_measured,heave_star_measured,pitch_measured,force_D,force_L,inertialload_y,...
    torque_x0,flowspeed_measured,heave_velo,heave_accel] = extract_measurements(transientcycs,freq,1/samplerate,Prof_out_angle,out,thcknss);

    flowspeed_measured_mean(ftrial,Atrial) = mean(flowspeed_measured);
    flowspeed_measured_stdev(ftrial,Atrial) = std(flowspeed_measured);
    f_star_commanded(ftrial,Atrial) = thcknss*freq/flowspeed_measured_mean(ftrial,Atrial);
    A_star_commanded(ftrial,Atrial) = (max(heave_commanded)-min(heave_commanded))/(2*thcknss);
    A_star_measured(ftrial,Atrial) = (max(heave_measured)-min(heave_measured))/(2*thcknss);

% Calculate heave acceleration from heave position
%     accel_heave = del2(heave_measured,T);

% Filter force data
    force_L_corrected = force_L+(foil.mass1+0.6)*heave_accel; %+inertialload_y;
    [b,a] = butter(6,10*freq*(2*T),'low'); % butterworth filter 6th order with cut-off frequency at 10*freq
    force_L_corrected_filtered = force_L_corrected; %filtfilt(b,a,squeeze(force_L_corrected)); %
    force_D_filtered =  force_D; %filtfilt(b,a,squeeze(force_D)); %
    torque_x0_filtered = filtfilt(b,a,squeeze(torque_x0));
    force_scale(ftrial,Atrial) = 0.5*1000*thcknss*span*flowspeed_measured_mean(ftrial,Atrial)^2;
    liftcoef = force_L_corrected_filtered/force_scale(ftrial,Atrial);
    dragcoef = force_D_filtered/force_scale(ftrial,Atrial);
    dragtorquecoef = -torque_x0_filtered/(force_scale(ftrial,Atrial)*span/2);

% Calculate power
    power_scale(ftrial,Atrial) = 0.5*1000*thcknss*span*flowspeed_measured_mean(ftrial,Atrial)^3;
    power_fluid = force_L_corrected_filtered .*heave_velo;

    m_star = 10; damping_ratio = 0;
    damping_coef = 0; %4*pi*freq*m_star*(pi*(1.5*0.0254/2)^2*10*1.5*0.0254)*damping_ratio; % dimensions of kg/s
    power_damping = -damping_coef*heave_velo.^2;
    power_total = power_fluid + power_damping;
    power_mean(ftrial,Atrial) = mean(power_total);
    powercoef = power_total/power_scale(ftrial,Atrial);
    powercoef_mean(ftrial,Atrial) = power_mean(ftrial,Atrial)/power_scale(ftrial,Atrial);
%     powercoef_mean(ftrial,Atrial) = CL1*sin(phi_CL1);

% % heave spectrum
    duration = max(time_star/freq)/T;
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(0*window_duration);
    heave_star_hilbert =  hilbert(heave_star_measured);
    [heave_powerspec,f_heave] = pwelch(heave_star_hilbert,window_duration,overlap,[],1/T);
    [max_power,max_index] = max(10*log10(heave_powerspec));
    f_heave_dom = f_heave(max_index);


 % force corrected+filtered spectrum using Welch's method
    duration = round(max(time_star/freq)/T);
    window_duration = round(duration/2); % size of hilbert windows measured in samples
    overlap = round(window_duration*7/8);
    force_hilbert =  hilbert(liftcoef);
    [force_powerspec,f_force] = pwelch(force_hilbert,window_duration,overlap,[],1/T);
    [max_force_power,max_force_index] = max(10*log10(force_powerspec));
    f_force_dom = f_force(max_force_index);  
%     spacing = (f_force(2)-f_force(1))/f;
%     findpeaks(10*log10(force_powerspec),1/spacing,'MinPeakHeight',max_power-20)

    % Find phase delay between heave and force
%     [heave1,phi_heave1] = fourierSineAmpAndPhase(time_star,heave_star_measured); % Find the phase of the heave position
    [CL1,phi_CL1] = fourierSineAmpAndPhase(time_star,liftcoef);

%     % Find phase delay between heave and force
%     maxlag = round(1/(2*freq*T));  % Maximum lag to calculate xcorr (in timesteps)
%     corr_wave = heave_measured/max(heave_measured);
%     [corrs,lags] = xcorr(force_L_corrected_filtered,corr_wave,maxlag);
% %     delay = freq*T*360*finddelay(heave_measured,force_y_corrected_filtered,maxlag)
%     figure
%     hold on
%     plot(lags*freq*T*360,corrs)
%     [forcedelay_peak,forcedelay_peaklocs] = findpeaks(corrs,1);
%     forcedelay_deg = forcedelay_peaklocs*T*freq*360-180; % Delay in degrees
%     stem(forcedelay_deg,forcedelay_peak);
%     xlabel('Lag (degrees)'); ylabel('Correlation')
%     xlim([-180 180])
%     hold off
%     drawnow


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
   
    plot_PrescribedMotionForceAndVelocity(time_star,heave_star_measured,heave_velo,liftcoef,...
        dragcoef,powercoef,num_cyc,dragtorquecoef);
%     drawnow
%     plot_ForceFourierAnalysis(phi_CL1,CL1,f_star_commanded(ftrial,Atrial),A_star_measured(ftrial,Atrial))

    end

    % Energy method analysis (a la Morse and Williamson 2009)
    m_star = 10;
    m_star_added = 1;
    damping_ratio = 0;
    force_inphasewvelo(ftrial,Atrial) = mean(CL1*sin(mean(phi_CL1)));
    force_inphasewaccel(ftrial,Atrial) = mean(CL1*cos(mean(phi_CL1)));
    m_star_added_eff(ftrial,Atrial) = force_inphasewaccel(ftrial,Atrial)/(2*pi^3*A_star_measured(ftrial,Atrial))*(1/f_star_commanded(ftrial,Atrial))^2;
    energyin_star(ftrial,Atrial) = force_inphasewvelo(ftrial,Atrial);
    %/(4*pi^3*A_star_measured(ftrial,Atrial))*(1/f_star_commanded(ftrial,Atrial))^2*sqrt((m_star+m_star_added)/(m_star+m_star_added_eff));
    energyout_star(ftrial,Atrial) = (m_star+m_star_added)*damping_ratio;

% %     mass = 10; % simulated mass in kg
% %     stiffness = 10; % simulated stiffness in kg/s^2
% %     damping_coef = 0; % simulated damping in kg/s
% %     force_scale_accel = (mass*A_star_measured(ftrial,Atrial)*thcknss*(2*pi*freq)^2);
% %     energyloss = damping_coef/sqrt(stiffness*mass);
% %     force_inphasewvelo(ftrial,Atrial) = mean(CL1*sin(mean(phi_CL1)));
% %     
% % %     force_inphasewvelo_star(ftrial,Atrial) = force_inphasewvelo(ftrial,Atrial)/force_scale_accel;
% %     force_inphasewaccel(ftrial,Atrial) = mean(CL1*cos(mean(phi_CL1)));
% % %     force_inphasewaccel_star(ftrial,Atrial) = force_inphasewaccel(ftrial,Atrial)/force_scale_accel;
% % %     energybalance(ftrial,Atrial) = force_inphasewvelo_star(ftrial,Atrial)*(1+force_inphasewaccel_star(ftrial,Atrial))^(-1/2)-energyloss; % -c/sqrt(k*m)

    liftcoef_alltrials(ftrial,Atrial) = mean(liftcoef);
    liftcoef_max(ftrial,Atrial) = max(liftcoef);
    liftcoef_min(ftrial,Atrial) = min(liftcoef);
    dragcoef_alltrials(ftrial,Atrial) = mean(dragcoef);
    dragcoef_max(ftrial,Atrial) = max(dragcoef);
    dragcoef_min(ftrial,Atrial) = min(dragcoef);    
    inertialload_alltrials(ftrial,Atrial) = mean(inertialload_y);
end

 
end
   if manytrial_analysis==1
    figure
    plot_PowerCoefContour;
    drawnow
    end
