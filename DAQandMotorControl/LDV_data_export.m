%% LDV Data Export 20220901

% This one works pretty well

clear;

% From the SPEED file:
%
% 1     Device Time     (msec)
% 2     Device Time     (usec)
% 3     Speed           (m/s)
% 4     SNR             
% 
% From the STATISTICS file:
% 
% 1     Speed Mean              (m/s)
% 2     Speed RMS               (m/s)
% 3     SNR Mean                
% 4     SNR RMS                 
% 5     Datarate                (Hz)
% 6     Number of samples       
% 7     Duration                (sec)
% 8     Internal Limit Exceeded (1 = yes)

%% load data

% LDV data directory:
FOLDERNAME = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221011_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68\20221011_Eric_large_sweep_inter-foil_wake_4');
% Speed data file name:
FILE = ('inter-foil_wake_A3E_aT4_p_h_ph_2022-10-11-14-16-46.SPEED.MSEBP.txt');
FILENAME = fullfile(FOLDERNAME,FILE); % full file name for loading data
LDV = importdata(FILENAME); % importing data into a double var array
% Statistics file name:
% FILE = ('inter-foil_wake_A3E_aT4_p_h_ph_2022-09-29-19-55-03..STATISTICS.MSEBP.txt');
% FILENAME = fullfile(FOLDERNAME,FILE); % full file name for loading data
% LDVst = importdata(FILENAME); % improting data into a double var array

% n = 1;
% for i = 1:length(LDV(:,3)) % filtering measurements below a certain SNR value
%     if LDV(i,4) > 6
%         U_filt(n,2) = LDV(i,3);
%         U_filt(n,1) = LDV(i,1);
%         n = n + 1;
%     end
% end

%% Plots

figure(2)

start_time_limit = [14,00,00]; % approximate starting time of data acquisition for the axes
end_time_limit = [24,00,00]; % approximate ending time of data acquisition for the axes

start_time = [14,16,46]; % EXACT starting time of data acquisition for the conversion
LDV(:,5) = LDV(:,1)/60/60/1000 + (start_time(1) + start_time(2)/60 + start_time(3)/60/60);

plot(LDV(:,5),LDV(:,3),'k'); hold on;

% optional (change with every dataset)

xline(14.2667,'r','LineWidth',4); % time at which first set started (hrs)
xline(17.2167,'b','LineWidth',4); % time at which second set started (hrs)
xline(20.1167,'m','LineWidth',4); % time at which third set started (hrs)
xline(23.0333,'g','LineWidth',4); % time at which fourth set started (hrs)

legend('','$\alpha_{T/4} = 0.16$','$\alpha_{T/4} = 0.33$','$\alpha_{T/4} = 0.68$','Single foil','interpreter','latex','fontsize',26);

% end of optional

if end_time_limit(1) < start_time_limit(1)
    time_ticks = [start_time_limit(1):1:24,1:1:end_time_limit(1)];
    time_value = start_time_limit(1) : 1 : start_time_limit(1) + length(time_ticks);
    xticks(time_value);
else
    time_ticks = start_time_limit(1):1:end_time_limit(1);
    time_value = time_ticks;
    xticks(time_value);
end

xticklabels(string(time_ticks));
set(gca, 'FontSize', 18, 'LineWidth', 1.5, 'TickLabelInterpreter', 'latex');

xlabel('time of the day (24hr)', 'interpreter', 'latex');
ylabel('flowspeed (m/s)', 'interpreter', 'latex');

xlim([start_time_limit(1),(start_time_limit(1) + length(time_ticks))]);

