%% Calculate Tandem Foil Experiments Results
% Eric Handy, Oct 2022

clear;

% main_dir = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl');
% cd(main_dir);

% Folderpath where files are stored

% folderpath = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\');
folderpath = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221114_TandemMonday_redux_LeadingAlphaSweep\data\');

% Sweeping parameters

p2 = [40 50 70];

aT4 = [0.155, 0.33, 0.679]; laT4 = length(aT4);
p3 = [65    70    75]; lp3 = length(p3);
h3 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh3 = length(h3);
ph = [-180  -120   -60     0    60   120]; lph = length(ph);

% Initialize variables

Yp = NaN(laT4,lp3,lh3,lph); % swept area

Eff_2 = Yp;
Eff_3 = Yp;
Eff_2_std = Yp;
Eff_3_std = Yp;

Eff_2prime = Yp;
Eff_3prime = Yp;
Eff_2prime_std = Yp;
Eff_3prime_std = Yp;

Eff_sys = Yp; % uses the freestream (leading) and wake (trailing) flows
Eff_sys_frstrm = Yp; % only uses the freestream
Eff_sys_std = Yp; % from using only the fresstream

CPH2 = Yp;
CPP2 = Yp;
CPH3 = Yp;
CPP3 = Yp;
beta2 = Yp;
beta3 = Yp;
U_2prime = Yp;
U_3prime = Yp;
Uwake = Yp;

%% Load FIRST experiment in the whole set

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_TandemFoil_APHPhaseSweep_A3E_alpha=0.155_p3=65_h3=0.7c_phase=-180.mat');

transientcycs = EP.transientcycs;

start_date = currentTime(1:3); % [yyyy,mm,dd] from the first dataset in the experiment

time0 = 60*60*1000*currentTime(4) + 60*1000*currentTime(5) + 1000*currentTime(6); % starting time in [ms]
LDV_skip_cycs = offset*(1/freq)*1000; % time skipped by the ldv in [ms]

%% Using LDV data?

using_ldv = 1; % 1/0 = yes/no

%% Load LDV data

if using_ldv == 1
    
    % LDV data directory:
    FOLDERNAME = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_Eric_large_sweep_inter-foil_wake_4\');
    % Speed data file name:
    FILE = ('inter-foil_wake_A3E_aT4_p_h_ph_2022-10-11-14-16-46.SPEED.MSEBP.txt');
    FILENAME = fullfile(FOLDERNAME,FILE); % full file name for loading data
    LDV = importdata(FILENAME); % importing data into a double var array
    
    figure(1); % shows the first couple of datapoints so that the first one can be manually chosen
    xlimit_plt = round(0.05*length(LDV(:,1)));
    plot(LDV(1:xlimit_plt,1),LDV(1:xlimit_plt,3),'ko'); hold on;
    plot(LDV(1:xlimit_plt,1),LDV(1:xlimit_plt,3),'g-');
    title('LDV Experimental Data Start');
    
    % % time at which experiment data begins in ms:
    LDV_data_start = input('Review plot and determine starting time. \nTime at which experiment data begins in LDV data (ms): ');
    % LDV_data_start = 61625; % in the event that the starting datapoint has been determined
    
    [LDV_start_time, LDV_index_start] = min(abs(LDV(:,1)-LDV_data_start)); % index at which the LDV acquisition starts
    % % NOTE: basically LDV_start_time will be our "time zero", to which all experiment times will be referenced to
    
    xline(LDV(LDV_index_start,1),'r'); % draws a line on the first datapoint of the parameter sweep
    
end

%% Calculation loop

% Baseline maeasurement correction
baseline_correction = 1; % correct for drift in measurement? (only for Wallace)
if baseline_correction == 1
    % load baseline measurement data
    load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_Baseline_TandemFoils\20221011_Baseline_TandemFoil_f=0.12_p2=70_h2=0.8c_p3=75_h3=0.9c_A3E.mat');
    baseline_forceN_displacement = mean(out(:,7));
%     baseline_forceT_displacement = mean(out(:,8));
    baseline_forceM_displacement = mean(out(:,12));
end

if using_ldv == 1
    figure(2); % shows all the LDV measurements
    plot(LDV(:,1),LDV(:,3)); hold on;
end

tic;
for ii = 1:laT4
    for jj = 1:lp3
        for kk = 1:lh3
            for qq = 1:lph
                %% Load force data
                
                filename = ['20221011_TandemFoil_APHPhaseSweep_A3E_',...
                    'alpha=', num2str(aT4(ii),4), '_p3=', num2str(p3(jj),2), '_h3=', num2str(h3(kk),2),'c_phase=', num2str(ph(qq)), '.mat'];

                filepath = fullfile(folderpath,filename);
                load(filepath); % make sure that this is not a reused variable

                % Drift correction for Wallace based on initial baseline measurements
                if baseline_correction == 1
                    out(:,7) = out(:,7) - mean(out(:,7)) + baseline_forceN_displacement;
%                     out(:,8) = out(:,8) - mean(out(:,8)) + baseline_forceT_displacement; % <-- this won't work because the measurement is not symmetric
                    out(:,12) = out(:,12) - mean(out(:,12)) + baseline_forceM_displacement;
                end
                
                % For drifting encoder data:
                % out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 2022/09/22-29 (pitch encoder damaged, change cable driver)

                %% Extract corresponding LDV data
                
                if using_ldv == 1
                    
                    time_exp_start = 60*60*1000*currentTime(4) + 60*1000*currentTime(5) + 1000*currentTime(6); % starting time of the currently loaded experiment
                    if currentTime(3) > start_date(3)
                        time_exp_start = time_exp_start + 24*60*60*1000; % adds 24 hours to the time count due to a change of day in case of overnight experiments
                    end
                    
                    relative_start = time_exp_start - time0; % time start relative to the first experiment
                    relative_strt_LDV = relative_start + LDV_skip_cycs + LDV(LDV_index_start,1) ; % relative time in ms at which the experimental data begins in the LDV data
                    relative_stop_LDV = relative_strt_LDV + (num_cyc+2*transientcycs)*(1/freq)*1000 - 2*LDV_skip_cycs; % start of the experiment + total cycles in ms - offset cycles
                    
                    [~, tstrt] = min(abs(LDV(:,1)-relative_strt_LDV)); % finds the LDV data index corresponding to the start of LDV daq for the current experiment
                    [~, tstop] = min(abs(LDV(:,1)-relative_stop_LDV)); % finds the LDV data index corresponding to the stop of LDV daq for the current experiment
                    
                    U_wake = mean(LDV(tstrt:tstop,3)); % takes the mean of the wake flow velocity found for the current experiment
                    
                    xline(LDV(tstrt,1),'k'); % overlays the start, stop and mean values of the current wake flow velocity
                    xline(LDV(tstop,1),'k');
                    plot([LDV(tstrt,1),LDV(tstop,1)],[U_wake,U_wake],'r');
                    
                end
                
                %% Extract parameters and perform force calculations
                
                if using_ldv == 1
                    [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, fs, transientcycs, foil_separation, flume_height, U_wake);
                else
                    [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, fs, transientcycs, foil_separation, flume_height);
                    U_wake = par.U;
                end
                
                res = calculate_forces(par, kin, out);

                %% Store calculated values

                Yp(ii,jj,kk,qq) = res.Yp; % maximum swept distance (from the two foils)
                
                Eff_2(ii,jj,kk,qq) = res.Eff_2; % total leading efficiency
                Eff_3(ii,jj,kk,qq) = res.Eff_3; % total trailing efficiency
                
                Eff_2_std(ii,jj,kk,qq) = res.Eff_2_std; % self-explanatory
                Eff_3_std(ii,jj,kk,qq) = res.Eff_3_std;
                
                Eff_2prime(ii,jj,kk,qq) = res.Eff_2prime; % corrected leading efficiency
                Eff_3prime(ii,jj,kk,qq) = res.Eff_3prime; % corrected trailing efficiency
                
                Eff_2prime_std(ii,jj,kk,qq) = res.Eff_2prime_std; % self-explanatory
                Eff_3prime_std(ii,jj,kk,qq) = res.Eff_3prime_std;
                
                Eff_sys(ii,jj,kk,qq) = res.Eff_sys; % system efficiency
                Eff_sys_std(ii,jj,kk,qq) = res.Eff_sys_std; % obvious, duh
                
                CPH2(ii,jj,kk,qq) = mean(res.CPH2); % leading heave power coeff
                CPP2(ii,jj,kk,qq) = mean(res.CPP2); % leading pitch power coeff
                CPH3(ii,jj,kk,qq) = mean(res.CPH3); % trailing heave power coeff
                CPP3(ii,jj,kk,qq) = mean(res.CPP3); % trailing pitch power coeff
                
                beta2(ii,jj,kk,qq) = res.beta2; % leading blockage ratio
                beta3(ii,jj,kk,qq) = res.beta3; % trailing blockage ratio
                
                U_2prime(ii,jj,kk,qq) = res.U_2prime; % corrected flowspeed
                U_3prime(ii,jj,kk,qq) = res.U_3prime; % corrected flowspeed in front of the trailing foil
                
                Uwake(ii,jj,kk,qq) = U_wake; % mean wake velocity from the LDV measurements
                
            end
        end
    end
end
hold off;

toc;

%% Save data

save('20221011_TandemFoil_efficiency_A3E_a155_330_679_PHPh_CpFrstrm_EffFrstrm_SysEffFrstrm_wBaseline.mat',...
    'Yp', 'Eff_2', 'Eff_3', 'Eff_2_std', 'Eff_3_std', 'CPP2', 'CPH2', 'CPP3', 'CPH3',...
    'Eff_sys', 'Eff_sys_frstrm', 'Eff_sys_std', 'Eff_2prime', 'Eff_3prime', 'Eff_2prime_std', 'Eff_3prime_std',...
    'beta2', 'beta3', 'U_2prime', 'U_3prime', 'Uwake');

