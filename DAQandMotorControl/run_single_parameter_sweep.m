%% Run Single Parameter Sweep

% June 2022, Eric Handy

% Runs experiments while sweeping over a single parameter

addpath(genpath("Libraries")); % add path to the analysis code

%% Where to save

path = ('R:\ENG_Breuer_Shared\ehandyca\Main\Data\'); % path where the new directory will be saved

foldername = sprintf('%s_SingleWednesday_SingleFoil_E1_FvH', datetime('now','Format','yyyyMMdd')); % name of the new directory

mkdir('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Main\Data', foldername); % makes the new directory
folderpath = fullfile(path, foldername); % assembles the path of the new directory

%% Trial Run

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
        'A2', 0.4, 65, 1, 70, 0.8, 180, 90, 10, 3, 0); % check if phi is positive or negative
% foiltype, freq, pitch2, heave2, pitch3, heave3, phase13, phi, num_cyc, transientcycs, constantamplitude);

% save(sprintf('%s_TrialRun_%s_f=%1.2f_0.mat', datetime('now','Format','yyyyMMdd'), foiltype, freq), 'dat', 'out', 'Prof_out_angle', 'Prof_out', 'freq', 'foiltype', 'num_cyc');
trialname = [datestr(now,formatOut), '_TrialRun_', foiltype, '_f=', num2str(freq,3), '_0.mat'];
save(trialname);

fprintf('Checklist:\n  - Make sure everything is in order\n')
fprintf('Press any key to continue (Crtl+C to Cancel)\n\n')
pause

fprintf('Running main experiment...\n')

%% General parameters

foiltype = 'A2'; % from foils_database
[foil, rho, fs] = foils_database(foiltype);

nu = 1.0035e-6; % kinematic viscosity of water @ 20C (~20.85 read from the vectrino - 20220401)
U = 0.4; % flow velocity
Re = U*foil.chord/nu; % Reynolds number


%% Running parameters


% freq     = 0.5;  % real frequency
% fred = 0.12; % reduced frequency
% freq = (fred*U)/foil.chord; % real frequency

pitch2   = 65;   % Gromit, in degrees
heave2   = 1;    % Gromit, in chords
pitch3   = 75;   % Gromit, in degrees
heave3   = 0.8;  % Gromit, in chords

phase13  = 180;  % phase between both rigs
phi      = 90;   % phase between heaving and pitching

num_cyc  = 15;   % number of full-amplitude experimental cycles

transientcycs     = 3; % number of transient cycles
constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position (to enable, modify "generate_profile" from the "run_Motors" code)

%% Sweeping Parameter

fred = 0.08:0.01:0.18;
freq = (fred*U)/foil.chord;

%% Run experiment

startexp = tic;
exp = 0;

for j = 1:length(freq)
    
    alphaT4 = atan(-2*pi*heave2*fred(j)) + deg2rad(pitch2); % self-explanatory
    
    [flume, out, dat, Prof_out_angle, Prof_out, last_out, freqj, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
        foiltype, freq(j), pitch2, heave2, pitch3, heave3, phase13, phi, num_cyc, transientcycs, constantamplitude);
    
    file = [datestr(now,formatOut), '_SingleFoil_Fsweep_', foiltype, '_f=', num2str(fred(j),3), '.mat'];
    filename = fullfile(foldername, file);
    
    save(trialname, 'dat', 'Prof_out_angle', 'out', 'foiltype', 'fred', 'freqj', 'Re', 'alphaT4', 'num_cyc');
    
    pause(20)
    exp = exp + 1;
    fprintf(sprintf('\nFinished experiment %i \n\n Running next experiment... \n\n',exp));
    
end

endexp = toc;

disp(['Done \n Full run took ',num2str(endexp,3),' seconds to complete.'])
