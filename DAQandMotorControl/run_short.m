%% Run Short

% 2022 06 13
% Quick code that runs Wallace and Gromit in tandem for a single experiment

% EP --> contains curent Experiment Parameters

%% Experiment Parameter Definition

EP.foiltype = 'A2'; % from foils_database

EP.freq     = 0.63;  % real frequency
% EP.fred = 0.12; % reduced frequency
% EP.U = 0.4; % flow velocity
% EP.freq = (fred*U)/chord; % real frequency

EP.P2 = 65; % Gromit, in degrees
EP.H2 = 1;  % Gromit, in chords (non-dimensional)
EP.P3 = 65; % Wallace, in degrees
EP.H3 = 1;  % Wallace, in chords (non-dimensional)

EP.phase13 = 180; % phase between both rigs
EP.phi     = -90; % phase between heaving and pitching

EP.num_cyc  = 30;  % number of full-amplitude experimental cycles

EP.transientcycs     = 3; % number of transient cycles
EP.constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position

%% Run experiment

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
    EP.foiltype, EP.freq, EP.P2, EP.H2, EP.P3, EP.H3, EP.phase13, EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude);

%% Save

experiment_name = [datestr(now,formatOut), '_ShortRun_', foiltype, '.mat'];
save(experiment_name, 'dq', 'foiltype', 'dat', 'out', 'Prof_out_angle', 'EP');

