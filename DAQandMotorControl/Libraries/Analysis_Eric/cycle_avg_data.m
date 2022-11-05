%% Cycle Averaging

% Cycle-averages a data set over a given reference data set (raw_pitch),
% and optionally over a specific number of cycles.

% raw_pitch : measured/commanded pitch (must not have noise)
% data : data to be averaged
% srate : sampling rate
% cycle_number : desired number of cycles to be averaged. Default is the total number of cycles available.

% Eric Handy, Jun 2022 - last edit

function [toverT, pitch_cycle, data_cycle] = cycle_avg_data(raw_pitch, data, srate, zero_start, cycle_number)

% zero_start == 1 to enforce starting the cycle at max pitch

if ~exist('zero_start','var')
    zero_start = 0;
end

if ~exist('srate','var')
    % parameter does not exist, default it to the following:
    srate = 1000; % sampling frequency [Hz]
end

% srate = 100000; % yuanhang's data

t = 0:1/srate:length(raw_pitch)*(1/srate)-(1/srate);
t = t';

% finding peaks and their locations in the pitch data
min_peak = mean(abs(raw_pitch))*0.9; % determines the minimum value that the peaks should have
[~,loc] = findpeaks(raw_pitch,'MinPeakHeight',min_peak);
locs = loc - loc(1) + 1;
T = NaN(length(locs)-1,1);

% time period for each cycle
for i=1:length(locs)-1
%     if i == length(locs)
% %         T(i) = t(end) - t(locs(i)); % in an ideal world this would be suficient
%         T(i) = T(i-1); % But since the last cycle may or may not be counted, better to just
%                        % make it equal to the previous, that way we avoid calculation errors.
%     else
        T(i) = t(locs(i+1)) - t(locs(i));
%     end
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
% plot(t,data);
% xline(t(locs), '--r'); % identified cycles
% xline(t(locs(1:cycle_number)), 'g'); % cycles to be averaged
% ylabel('\theta');
% xlabel('t');
% title(['Cycles averaged: ', num2str(cycle_number)]);

T_avg = mean(T(1:cycle_number)); % mean time-period considering only the cycles that will be phase-averaged
freq = 1/T_avg; % real frequency from the period (for debugging)
cycle_length = round(T_avg*srate); % averaged number of data points in one cycle

pitch_cycle = zeros(cycle_number,cycle_length); % cycle-averaged pitching angle
data_cycle = zeros(cycle_number,cycle_length);  % cycle-averaged torque

% Data averaging loop
for i=1:cycle_length-1
    pitch_cycle(:,i) = raw_pitch( i : cycle_length : i+cycle_length*cycle_number-1 );
    data_cycle(:,i)  =      data( i : cycle_length : i+cycle_length*cycle_number-1 );
end

pitch_cycle(:,end) = (pitch_cycle(:,end-1)+pitch_cycle(:,1))*0.5; % this is more of a quick-fix, might want to redo some of this code
data_cycle(:,end) = (data_cycle(:,end-1)+data_cycle(:,1))*0.5;

t_avg = 0:1/srate:length(mean(pitch_cycle))/srate-1/srate; % actual averaged time per cycle [s]
toverT = t_avg/T_avg; % normalized time over the period

% If data should start minimum value (to make plotting symmetric)
if zero_start == 1
    [~, I] = min(pitch_cycle,[],2); % find location of minimum values in every cycle
    beginning = pitch_cycle(:,1:I-1); % everything to the left of the minimum value is the beginnning
    beginning_dat = data_cycle(:,1:I-1);
    ending = pitch_cycle(:,I:end); % everything to the right of the minimum value is the ending
    ending_dat = data_cycle(:,I:end);
    pitch_cycle = [ending,beginning]; % reassign the pitch_cycle variable while reorganizing the cycle data
    data_cycle = [ending_dat,beginning_dat];
end

% t_avg = 0:1/srate:length(mean(pitch_cycle))/srate-1/srate; % actual averaged time per cycle [s]
% toverT = t_avg/T_avg; % normalized time over the period

end