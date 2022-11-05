%% LDV Data Export 20220901

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
FOLDERNAME = ('C:\Users\ehandyca\Documents\MSE Working Directory\20220930_testing_ldv_cylinder');
% Speed data file name:
FILE = ('testing_ldv_cylinder_400_loops_2022-09-30-13-33-28.SPEED.MSEBP.txt');
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

figure()

start_time_limit = [12,00,00]; % approximate starting time of data acquisition for the axes
end_time_limit = [18,00,00]; % approximate ending time of data acquisition for the axes

start_time = [13,33,28]; % EXACT starting time of data acquisition for the conversion
LDV(:,5) = LDV(:,1)/60/60/1000 + (start_time(1) + start_time(2)/60 + start_time(3)/60/60);

plot(LDV(:,5),LDV(:,3),'k');


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

xlabel('time of the day', 'interpreter', 'latex');
ylabel('flowspeed', 'interpreter', 'latex');

xlim([start_time_limit(1),(start_time_limit(1) + length(time_ticks))]);