%% Phase Averaging

% Cycle-averages a data set over a given reference data set (raw_pitch),
% and optionally over a specific number of cycles.

% NOTE: For a sampling frequency of fs = 1000 Hz

% Eric Handy, Dec 2021
% Eric Handy, May 2022
% Eric Handy, Jun 2022 - last edit


function [toverT, pitch_cycle, data_cycle] = phase_avg_data(raw_pitch, data, cycle_number) % cycle_number refers to the number of cycles you want averaged

% fs = EP.srate; % sampling freq [Hz]
fs = 1000;
% fs = 100000; % yuanhang's data

t = 0:1/fs:length(raw_pitch)/fs-1/fs;
t = t';

% finding peaks and their locations in the pitch data
min_peak = mean(abs(raw_pitch))*0.9; % determines the minimum value that the peaks should have
[pks,locs] = findpeaks(raw_pitch,'MinPeakProminence',min_peak);
T = NaN(length(locs)-1,1);

% time period for each cycle
for i=1:1:length(locs)-1
    T(i) = t(locs(i+1)) - t(locs(i));
end

start = 1; % first cycle to be averaged, skipping the first three
tot_cycles = length(T); % total number of cycles identified by the code

if ~exist('cycle_number','var')
    % parameter does not exist, default it to the following:
    cycle_number = tot_cycles - (start-1) - 3;  % Number of cycles to be averaged: the total cycles minus the
                                                % first 3 cycles (if starting cycle is 4) and last two cycles.
                                                % To use the var "cycle_number", disable line 30(?).
end

% figure() % shows which cycles will be averaged
% plot(t,raw_pitch)
% hold on
% plot(t(locs),pks,'o') % plots all peaks
% plot(t(locs(start:cycle_number+start)), pks(start:cycle_number+start),'m*') % cycles to be averaged
% ylabel('\theta');
% xlabel('t');
% title(['Cycles averaged: ', num2str(cycle_number)]);

T_avg = mean(T(start:cycle_number));  % mean time-period considering only the cycles that will be phase-averaged
freq = 1/T_avg; % real frequency from the period (just to verify)
cycle_length = round(T_avg*fs);     % averaged number of data points in one cycle

pitch_cycle = NaN(cycle_number+1,cycle_length); % cycle-averaged pitching angle
data_cycle = NaN(cycle_number+1,cycle_length);  % cycle-averaged torque

% Data averaging loop
for i=1:1:cycle_length
    pitch_cycle(:,i) = raw_pitch( locs(start)+i : cycle_length : locs(start)+i+cycle_length*cycle_number );
    data_cycle(:,i) = data( locs(start)+i : cycle_length : locs(start)+i+cycle_length*cycle_number );
end

t_avg = 0:1/fs:length(mean(pitch_cycle))/fs-1/fs; % actual averaged time per cycle [s]
toverT = t_avg/T_avg; % normalized time over the period

end