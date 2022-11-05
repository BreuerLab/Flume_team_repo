%% Run Short

% 2022 06 13
% Quick code that runs Wallace and Gromit in tandem for a single experiment

%% Experiment Parameter Definition

foiltype = 'A3E'; % from foils_database

% freq = 0.5;  % real frequency
fred = 0.12; % reduced frequency
U = 0.33; % flow velocity
freq = (fred*U)/chord; % real frequency
% freq = 0.6492;
fred = freq*0.061/U;

P2 = 0; % Gromit, in degrees
H2 = 0;  % Gromit, in chords (non-dimensional)
P3 = 0; % Wallace, in degrees
H3 = 0;  % Wallace, in chords (non-dimensional)

phase13 = -180; % phase between both rigs
phi     = -90; % phase between heaving and pitching

num_cyc = 30;  % number of full-amplitude experimental cycles

transientcycs     = 3; % number of transient cycles
constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position

betah = 0; % trapezoidal profile parameter <--- for traingular shape
betap = 0; % <-- for trapezoid shape

%% Run experiment

time_start = clock;

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
    foiltype, freq, P2, H2, P3, H3, phase13, phi, num_cyc, transientcycs, constantamplitude, betah, betap, 4);

time_stop = clock;

%% Run diagnostics
% 
% fs = dq.Rate; % sampling frequency
% 
% [diagnostics, continue_exp] = run_diagnostics(Prof_out_angle, out, fs, freq);

%% Save

% EP.srate = dq.Rate; % extract the sampling rate from the daq setup

experiment_name = [datestr(now,'yyyymmdd'), '_EndBaseline_TandemFoil_noMotion_', foiltype, '.mat'];
% save(experiment_name, 'foiltype', 'dat', 'out', 'Prof_out_angle', 'diagnostics', 'flume_height');
save(experiment_name);

