% Calculate Single Foil Experiments Results
% Eric Handy, Sept 2022

clear;

% main_dir = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl');
% cd(main_dir);

% Folderpath where files are stored

folderpath = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\');

% Sweeping parameters

p2 = [65    70    75    80]; lp2 = length(p2);
h2 = [0.7000    0.8000    0.9000    1.0000    1.1000    1.2000]; lh2 = length(h2);

% Initialize variables

Yp = NaN(lh2,lp2); % swept area

Eff_2 = Yp;
Eff_2_std = Yp;

Eff_2prime = Yp;
Eff_2prime_std = Yp;

CPH2 = Yp;
CPP2 = Yp;
beta2 = Yp;
U_2prime = Yp;
Uwake = Yp;

%% Load FIRST experiment in the whole set

load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_SingleFoil_PHSweep_A3E_p2=65_h2=0.7.mat');

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
    
    figure(1); % shows the whole dataset so that the first experiment datapoint can be manually chosen
%     xlimit_plt = round(0.05*length(LDV(:,1)));
%     plot(LDV(1:xlimit_plt,1),LDV(1:xlimit_plt,3),'ko'); hold on;
%     plot(LDV(1:xlimit_plt,1),LDV(1:xlimit_plt,3),'g-');
    plot(LDV(:,1),LDV(:,3),'ko'); hold on;
    plot(LDV(:,1),LDV(:,3),'g-');
    title('LDV Experimental Data Start');
    
    % % time at which experiment data begins in ms:
    LDV_data_start = input('Review plot and determine starting time. \nTime at which experiment data begins in LDV data (ms): ');
    % LDV_data_start = 61625; % in the event that the starting datapoint has been determined
    
    [LDV_start_time, LDV_index_start] = min(abs(LDV(:,1)-LDV_data_start)); % index at which the LDV acquisition starts
    % % NOTE: basically LDV_start_time will be our "time zero", to which all experiment times will be referenced to
    
    xline(LDV(LDV_index_start,1),'r'); % draws a line on the first datapoint of the parameter sweep
    
end

%% Calculation loop

tic;
for ii = 1:lp2
    for jj = 1:lh2

    filename = ['20221011_SingleFoil_PHSweep_A3E_',...
        'p2=', num2str(p2(ii),2), '_h2=', num2str(h2(jj),2), '.mat'];

    filepath = fullfile(folderpath,filename);
    load(filepath);

%     out(:,5) = deg2rad(Prof_out_angle(:,5)); % for data taken on 20220812 (pitch encoder damaged, change cable driver)

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
        [kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, fs, EP.transientcycs, foil_separation, flume_height);
    end
    
    res = calculate_forces(par, kin, out);

    %% Store calculated values

    Yp(jj,ii) = res.Yp; % maximum swept distance (from the two foils)
    
    CPP2(jj,ii) = mean(res.CPP2);
    CPH2(jj,ii) = mean(res.CPH2);
    
    Eff_2(jj,ii) = res.Eff_2; % total leading efficiency
    Eff_2_std(jj,ii) = res.Eff_2_std;
    
    Eff_2prime(jj,ii) = res.Eff_2prime; % corrected leading efficiency
    Eff_2prime_std(jj,ii) = res.Eff_2prime_std;
    
    beta2(jj,ii) = res.beta2; % leading blockage ratio
    
    U_2prime(jj,ii) = res.U_2prime; % corrected flowspeed
    Uwake(jj,ii) = U_wake;
    
    end
end
toc;

save('20221011_SingleFoil_efficiency_A3E_PH.mat',...
    'Yp', 'CPP2', 'CPH2', 'Eff_2', 'Eff_2_std', 'Eff_2prime', 'Eff_2prime_std', 'beta2', 'U_2prime', 'Uwake');
