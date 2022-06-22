%% Cycle Averaging

% Cycle-averages a data set over a given reference data set (raw_pitch),
% and optionally over a specific number of cycles.

% raw_pitch : measured/commanded pitch (must not have noise)
% data : data to be averaged
% srate : sampling rate
% cycle_number : desired number of cycles to be averaged. Default is the total number of cycles available.

% Eric Handy, Jun 2022 - last edit

function [toverT, pitch_cycle, data_cycle] = cycle_avg_data(raw_pitch, data, srate, cycle_number)

if ~exist('srate','var')
    % parameter does not exist, default it to the following:
    srate = 1000; % sampling frequency [Hz]
end

% srate = 100000; % yuanhang's data

t = 0:1/srate:length(raw_pitch)/srate-1/srate;
t = t';

% finding peaks and their locations in the pitch data
min_peak = mean(abs(raw_pitch))*0.9; % determines the minimum value that the peaks should have
[~,loc] = findpeaks(raw_pitch,'MinPeakProminence',min_peak);
locs = loc - loc(1) + 1;
T = NaN(length(locs),1);

% time period for each cycle
for i=1:length(locs)
    if i == length(locs)
        T(i) = t(end) - t(locs(i));
    else
        T(i) = t(locs(i+1)) - t(locs(i));
    end
end

if ~exist('cycle_number','var')
    % parameter does not exist, default it to the following:
    cycle_number = length(T);
else
    if cycle_number > length(T) % ensure than cycle_number is not larger than available cycles
        error('Desired number of cycles to be averages is larger than available experimental cycles.');
    end
end

% % For debugging:
% figure() % shows which cycles will be averaged
% plot(t,raw_pitch); hold on;
% xline(t(locs), '--r'); % identified cycles
% xline(t(locs(1:cycle_number)), 'g'); % cycles to be averaged
% ylabel('\theta');
% xlabel('t');
% title(['Cycles averaged: ', num2str(cycle_number)]);

T_avg = mean(T(1:cycle_number)); % mean time-period considering only the cycles that will be phase-averaged
freq = 1/T_avg; % real frequency from the period (for debugging)
cycle_length = round(T_avg*srate); % averaged number of data points in one cycle

pitch_cycle = NaN(cycle_number,cycle_length); % cycle-averaged pitching angle
data_cycle = NaN(cycle_number,cycle_length);  % cycle-averaged torque

% Data averaging loop
for i=1:cycle_length-1
    pitch_cycle(:,i) = raw_pitch( i : cycle_length : i+cycle_length*cycle_number-1 );
    data_cycle(:,i)  =      data( i : cycle_length : i+cycle_length*cycle_number-1 );
end

t_avg = 0:1/srate:length(mean(pitch_cycle))/srate-1/srate; % actual averaged time per cycle [s]
toverT = t_avg/T_avg; % normalized time over the period

end