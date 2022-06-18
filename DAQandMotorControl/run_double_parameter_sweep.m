%% Run Triple Parameter Sweep

% June 2022, Eric Handy

% Runs experiments while sweeping over two parameters

addpath(genpath("Libraries")); % add path to the analysis code

%% Where to save

path = ('R:\ENG_Breuer_Shared\ehandyca\Data_main_repo\'); % path where the new directory will be saved

foldername = sprintf('%s_TandemFriday_AlphaSweep_PHPhase_A2E_a15', datetime('now','Format','yyyyMMdd')); % name of the new directory

mkdir('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo', foldername); % makes the new directory
folderpath = fullfile(path, foldername); % assembles the path of the new directory

%% Trial Run

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
        'A2E', 0.63, 45, 0.8, 60, 0.55, -180, -90, 11, 3, 0); % check if phi is positive or negative
% foiltype, freq, pitch2, heave2, pitch3, heave3, phase13, phi, num_cyc, transientcycs, constantamplitude);

trialname = [datestr(now,'yyyymmdd'), '_TrialRun_', foiltype, '_f=', num2str(freq,3), '_2.mat'];
save(trialname);

fprintf('Checklist:\n  - Make sure everything is in order\n')
fprintf('Press any key to continue (Ctrl+C to Cancel)\n\n')
pause

fprintf('Running main experiment...\n')

%% General parameters

EP.foiltype = 'A2E'; % from foils_database
[foil, rho, fs] = foils_database(foiltype);

nu = 1.0035e-6; % kinematic viscosity of water @ 20C (~20.85 read from the vectrino - 20220401)
U = 0.4; % flow velocity
Re = U*foil.chord/nu; % Reynolds number


%% Running parameters


% freq     = 0.5;  % real frequency
fred = 0.12; % reduced frequency
EP.freq = (fred*U)/foil.chord; % real frequency

Pitch2   = 40;   % Gromit, in degrees
Heave2   = 0.8;    % Gromit, in chords
% pitch3   = 75;   % Wallace, in degrees
% heave3   = 0.8;  % Wallace, in chords

% phase13  = 180;  % phase between both rigs
EP.phi      = -90;   % phase between heaving and pitching

EP.num_cyc  = 30;   % number of full-amplitude experimental cycles

EP.transientcycs     = 3; % number of transient cycles
EP.constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position (to enable, modify "generate_profile" from the "run_Motors" code)

alphaT4 = atan(-2*pi*Heave2*fred) + deg2rad(Pitch2); % self-explanatory
EP.alphaT4 = alphaT4;

EP.flume_depth = flume_height; % from the initial setup_DAQ

%% Sweeping Parameter

% fred = 0.08:0.01:0.18;
% freq = (fred*U)/foil.chord;
Pitch3   = 60:5:80;   % Wallace, in degrees
Heave3   = 0.55:0.15:1.3;  % Wallace, in chords
phase13 = -180:60:180; % phase between Wallace and Gromit

%% Run experiment

startexp = tic;
exp = 0;

for j = 1:length(phase13)
    for k = 1:length(Pitch3)
        for i = 1:length(Heave3)
%     alphaT4 = atan(-2*pi*heave2*fred) + deg2rad(pitch2); % self-explanatory
            
            EP.pitch3 = Pitch3(k); EP.pitch2 = Pitch2;
            EP.heave3 = Heave3(i); EP.heave2 = Heave2;
            EP.phase13 = phase13(j);

            [flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3k, heave3i, phase13j, num_cyc, phi, foiltype]...
                = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
                EP.foiltype, EP.freq, Pitch2, Heave2, Pitch3(k), Heave3(i), phase13(j), EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude);
            file = ['\', datestr(now,'yyyymmdd'), '_TandemFoil_PHPhaseSweep_', EP.foiltype,...
                '_p3=', num2str(Pitch3(k)), '_h3=', num2str(Heave3(i),3), 'c_phase=', num2str(phase13(j)), '.mat'];
            filename = fullfile(folderpath, file);
            
            save(filename, 'EP', 'dat', 'Prof_out_angle', 'out', 'foiltype', 'fred', 'freq', 'pitch3k', 'heave3i', 'phase13j', 'Re', 'alphaT4', 'num_cyc');
            
            pause(20)
            exp = exp + 1;
            fprintf(sprintf('\nFinished experiment %i \n\n Running next experiment... \n\n',exp));
        end
    end    
end

endexp = toc;

disp(['Done \n Full run took ',num2str(endexp,3),' seconds to complete.'])
