%% Run single foil aT4 sweep

% Nov 2022, Eric Handy

addpath(genpath("Libraries")); % add path to the analysis code

%% Where to save

path = ('R:\ENG_Breuer_Shared\ehandyca\Data_main_repo\'); % path where the new directory will be saved

foldername = sprintf('%s_SingleSaturday_SingleFoilAlphaSweep', datetime('now','Format','yyyyMMdd')); % name of the new directory

mkdir('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo', foldername); % makes the new directory
folderpath = fullfile(path, foldername); % assembles the path of the new directory

%% Trial Run

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
        'A3E', 0.6, 75, 0.8, 0, 0, 0, -90, 10, 3, 0, 0, 0, 4); % check if phi is positive or negative
% foiltype, freq, pitch2, heave2, pitch3, heave3, phase13, phi, num_cyc, transientcycs, constantamplitude);

trialname = [datestr(now,'yyyymmdd'), '_TrialRun_', foiltype, '_f=', num2str(freq,3), '_0.mat'];
save(trialname);

fprintf('Checklist:\n  - Make sure everything is in order\n')
fprintf('Press any key to continue (Ctrl+C to Cancel)\n\n')
pause

fprintf('Preparing main experiment...\n')

%% General parameters

EP.foiltype = 'A3E'; % from foils_database
[foil, rho, fs] = foils_database(foiltype);

EP.nu = 1.0035e-6; % kinematic viscosity of water @ 20C (~21.4 read from the vectrino - 20221105)
EP.U = 0.33; % expected flow velocity
EP.Re = U*foil.chord/nu; % Reynolds number


%% Running parameters


% freq     = 0.5;  % real frequency
% fred = 0.12; % reduced frequency
% EP.freq = (fred*U)/foil.chord; % real frequency

% Pitch2   = 40;   % Gromit, in degrees
% Heave2   = 0.8;    % Gromit, in chords
Pitch3   = 0;   % Wallace, in degrees
Heave3   = 0;  % Wallace, in chords

EP.phase13  = 0;  % phase between both rigs
EP.phi      = -90;   % phase between heaving and pitching

EP.num_cyc  = 40;   % number of full-amplitude experimental cycles

EP.transientcycs     = 5; % number of transient cycles
EP.constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position (to enable, modify "generate_profile" from the "run_Motors" code)

% alphaT4 = atan(-2*pi*Heave2*fred) + deg2rad(Pitch2); % self-explanatory
% EP.alphaT4 = alphaT4;

EP.flume_height = flume_height; % from the initial setup_DAQ

EP.offset = 4; % cycles after/before which LDV begins/ends taking data

%% Sweeping Parameters

fred = [0.1,0.12,0.15];
Pitch2 = [40,45,50,55,60,65,70,75,80];
Heave2 = [0.5,0.75,1,1.25,1.5];

freq = (fred*U)/foil.chord;
total_predicted_time = 0;
for num = length(fred)
    fsweep_time = length(Pitch2)*length(Heave2)*((1/freq(num))*(EP.num_cyc + 2*EP.transientcycs) + 20);
    total_predicted_time = total_predicted_time + (fsweep_time)/60/60;
end

msg = ['Total predicted experiment time:\n', num2str(total_predicted_time),' hours\n'];
fprintf(msg);
fprintf('Press any key to continue (Ctrl+C to Cancel)\n\n')
pause

fprintf('Running main experiment...\n')

%% Run experiment

start_date = datestr(now,'yyyymmdd');
startexp = tic;
exp = 0;

for i = 1:length(fred)
    for j = 1:length(Pitch2)
        for k = 1:length(Heave2)
%     alphaT4 = atan(-2*pi*heave2*fred) + deg2rad(pitch2); % self-explanatory
            EP.fred = fred(i);
            EP.P2 = Pitch2(j);
            EP.H2 = Heave2(k);
            EP.srate = dq.Rate;
            currentTime = clock;

            alphaT4 = atan(-2*pi*Heave2(k)*fred(i)) + deg2rad(Pitch2(j)); % self-explanatory
            EP.alphaT4 = alphaT4;
            
            if rad2deg(alphaT4) < 5
                continue
            end

            [flume, out, dat, Prof_out_angle, Prof_out,last_out, freqi, pitch2j, heave2k, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
                = run_Motors(dq,last_out,bias,...
                EP.foiltype, freq(i), Pitch2(j), Heave2(k), Pitch3, Heave3, EP.phase13, EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude, 0, 0, EP.offset);

            file = ['\', start_date, '_SingleFoil_aT4_Sweep_', EP.foiltype,...
                '_fred=', num2str(fred(i)), '_p2=', num2str(Pitch(j),2), '_Heave2=', num2str(Heave(k)), 'c_aT4=', num2str(alphaT4,3), '.mat'];

            [diagnostics, continue_exp] = run_diagnostics(Prof_out_angle,out,dq.Rate,freq); % run diagnostics code

            filename = fullfile(folderpath, file);
            save(filename);%, 'EP', 'dat', 'Prof_out_angle', 'out', 'foiltype', 'fred', 'freq', 'pitch2q', 'pitch3k', 'heave3i', 'phase13j', 'Re', 'alphaT4', 'num_cyc');

            pause(20)
            exp = exp + 1;
            fprintf(sprintf('\nFinished experiment %i \n\n Running next experiment... \n\n',exp));
        end
    end    
end

endexp = toc;

disp(['Done \n Full run took ',num2str(endexp,3),' seconds to complete.\n'])
