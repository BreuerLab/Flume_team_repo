%% Run Short

% 2022 06 13
% Quick code that runs Wallace and Gromit in tandem for a single experiment

% EP --> contains curent Experiment Parameters

%% Experiment Parameter Definition

EP.foiltype = 'A2'; % from foils_database

EP.freq     = 0.5;  % real frequency
% EP.fred = 0.12; % reduced frequency
% EP.U = 0.4; % flow velocity
% EP.freq = (fred*U)/chord; % real frequency

EP.P2 = 0; % Gromit, in degrees
EP.H2 = 0;  % Gromit, in chords (non-dimensional)
EP.P3 = 30; % Wallace, in degrees
EP.H3 = 0;  % Wallace, in chords (non-dimensional)

EP.phase13 = 180; % phase between both rigs
EP.phi     = -90; % phase between heaving and pitching

EP.num_cyc  = 10;  % number of full-amplitude experimental cycles

EP.transientcycs     = 3; % number of transient cycles
EP.constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position

EP.flume_height = flume_height; % from the initial setup_DAQ

EP.betah = 0; % trapezoidal profile parameter <--- for traingular shape
EP.betap = 0; % <-- for trapezoid shape

%% Run experiment

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
    EP.foiltype, EP.freq, EP.P2, EP.H2, EP.P3, EP.H3, EP.phase13, EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude, EP.betah, EP.betap);

%% Save

EP.srate = dq.Rate; % extract the sampling rate from the daq setup

experiment_name = [datestr(now,'yyyymmdd'), '_ShortRun_Wallace_Pitch_teknic_test2', foiltype, '.mat'];
save(experiment_name, 'foiltype', 'dat', 'out', 'Prof_out_angle', 'EP');

