function [U, P1, P2, H1, H2, f_real, f_red, phi, phase13, num_cyc, foil, rho, fs] = give_me_params(out, Prof_out_angle, foiltype)

% Exctracts all relevant parameters from a data file.
% Eric Handy, Mar 2022

% ------------------------------------------------------------------------------------------

% U  -------> average flow velocity (from vectrino) [m/s]
% P1 -------> leading pitch amp [deg]
% P2 -------> trailing pitch amp [deg]
% H1 -------> leading heave amp [meters]
% H2 -------> trailing heave amp [meters]
% f_real ---> real experiment oscillation frequency [Hz]
% f_red ----> reduced frequency [Hz]
% fs -------> sampling frequency [Hz]
% phi ------> phase between heave and pitch [deg]
% phase13 --> phase between leading and trailing foils [deg]
% num_cyc --> number of cycles
% chord ----> [m]
% AR -------> span/chord
% mass1 ----> leading foil mass [kg]
% mass2 ----> trailing foil mass [kg]
% rho ------> density of the fluid (water) [kg/m3]
% foiltype -> foil used (from the foils database) [string]

% ------------------------------------------------------------------------------------------

% >> From the "out" variable:
% - Cols 1 to 6 = pitching and heaving positions for the three rigs (starting with the foremost upstream)
% - Cols 13 to 16 = vectrino data.
%               13 --> flow velocity
% - Cols 17 to 22 = Leading foil data:
%               17 --> normal force
%               18 --> tangential force
%               22 --> spanwise torque
% - Cols 7 to 12 = Trailing foil data:
%                7 --> normal force
%                8 --> tangential force
%               12 --> spanwise torque

% ------------------------------------------------------------------------------------------

% >> From the "Prof_out_angle" variable:
% - Cols 1 to 6 = pitching and heaving positions for the three rigs (starting with the foremost upstream)
%               1 --> Shawn   pitching (defunct)
%               2 --> Shawn   heaving  (defunct)
%               3 --> Gromit  pitching (leading)
%               4 --> Gromit  heaving  (leading)
%               5 --> Wallace pitching (trailing)
%               6 --> Wallace heaving  (trailing)

% ------------------------------------------------------------------------------------------

%% Parameter extraction

% Foil properties

[foil, rho, fs] = foils_database(foiltype);

% Flow velocity

U = -mean(out(:,13)); % flow speed from the vectrino [m/s] (negative because of direction of the flow)

% Pitch and heave

P1 = round(max(abs(Prof_out_angle(:,3)))); % pitch from taking the signal to the motor and rounding the maximum absolute value
P2 = round(max(abs(Prof_out_angle(:,5))));
H1 = max(abs(Prof_out_angle(:,4))); % heave does the same thing as pitch but rounding the value differently
H2 = max(abs(Prof_out_angle(:,6))); % and multiplying times the chord length.

% % Void data determination
% 
% if P1 == 0
%     pitch = 
% 
% pitch_lead = Prof_out_angle(:,3); % leading pitch
% heave_lead = Prof_out_angle(:,4); % leading heave
% 
% pitch_trail = Prof_out_angle(:,5); % trailing pitch
% heave_trail = Prof_out_angle(:,6); % trailing heave

% Oscillation frequency and number of cycles

t = 0:1/fs:length(Prof_out_angle(:,5))/fs-1/fs; % time variable
t = t';

min_peak = mean(abs(Prof_out_angle(:,5)))*0.9; % determines the minimum value that the peaks should have
[pks,locs] = findpeaks(Prof_out_angle(:,5),'MinPeakProminence',min_peak);
T = NaN(length(locs)-1,1);

for i=1:1:length(locs)-1
    T(i) = t(locs(i+1)) - t(locs(i)); % time period for each cycle
end

start = 4; % first cycle to be averaged, skipping the first three
num_cyc = length(T); % total number of cycles identified by the code and run by the control code
good_cycles = num_cyc - (start-1) - 2; % only the cycles worth averaging

period = mean(T(start:good_cycles));  % mean time-period considering only the cycles that will be phase-averaged
f_real = round(1/period,2); % real frequency from the period

% Reduced frequency

f_red = (f_real/U)*foil.chord; % reduced frequency

% Phase diffference between trailing foil heave and pitch using cross-correlation

[xfc, lags] = xcorr(Prof_out_angle(:,5), Prof_out_angle(:,6)); % get the cross-correlation from both signals (pitching and heaving)
[~, I] = max(xfc); % obtain the index of the maximum correlation value, as we expect the signals to be most correlated at a lag
sig_lag = lags(I); % obtain the phase difference in sampled points
phi = round(sig_lag*360/((1/f_real)*1000)); % calculate the phase angle from the phase difference

% Phase between leading and trailing foil motions using cross-correlation

[xfc, lags] = xcorr(Prof_out_angle(:,6), Prof_out_angle(:,4)); % get the cross-correlation from both signals (heaving from both rigs)
[~, I] = max(xfc); % obtain the index of the maximum corrleation value, as we expect the signals to be most correlated at a lag
sig_lag = lags(I); % obtain the phase difference in sampled points
phase13 = round(sig_lag*360/((1/f_real)*1000)); % calculate the phase angle from the phase difference

end