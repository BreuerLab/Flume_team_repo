%% Run leading alpha sweep 
% Triple Parameter Sweep

% Sept 2022, Eric Handy

% Runs experiments while sweeping over three parameters

addpath(genpath("Libraries")); % add path to the analysis code

%% Where to save

path = ('R:\ENG_Breuer_Shared\ehandyca\Data_main_repo\'); % path where the new directory will be saved

foldername = sprintf('%s_TandemTuesday_AlphaSweep_APHPhase_A3E_a16_a33_a68', datetime('now','Format','yyyyMMdd')); % name of the new directory

mkdir('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Data_main_repo', foldername); % makes the new directory
folderpath = fullfile(path, foldername); % assembles the path of the new directory

%% Trial Run

[flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
        'A3E', 0.6492, 75, 0.8, 75, 1.2, -180, -90, 10, 3, 0, 0, 0, 4); % check if phi is positive or negative
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

nu = 1.0035e-6; % kinematic viscosity of water @ 20C (~20.85 read from the vectrino - 20220401)
U = 0.33; % flow velocity
Re = U*foil.chord/nu; % Reynolds number


%% Running parameters


% freq     = 0.5;  % real frequency
fred = 0.12; % reduced frequency
EP.freq = (fred*U)/foil.chord; % real frequency

% Pitch2   = 40;   % Gromit, in degrees
Heave2   = 0.8;    % Gromit, in chords
% pitch3   = 75;   % Wallace, in degrees
% heave3   = 0.8;  % Wallace, in chords

% phase13  = -180;  % phase between both rigs
EP.phi      = -90;   % phase between heaving and pitching

EP.num_cyc  = 40;   % number of full-amplitude experimental cycles

EP.transientcycs     = 3; % number of transient cycles
EP.constantamplitude = 0; % replace oscillating motion with moving the foil to a desired position (to enable, modify "generate_profile" from the "run_Motors" code)

% alphaT4 = atan(-2*pi*Heave2*fred) + deg2rad(Pitch2); % self-explanatory
% EP.alphaT4 = alphaT4;

EP.flume_height = flume_height; % from the initial setup_DAQ

offset = 4; % cycles after/before which LDV begins/ends taking data

%% Sweeping Parameter

Pitch2 = [40, 50, 70];

% fred = 0.08:0.01:0.18;
% freq = (fred*U)/foil.chord;
Pitch3   = [65,70,75];   % Wallace, in degrees
Heave3   = 0.7:0.1:1.2;  % Wallace, in chords
phase13 = -180:60:120; % phase between Wallace and Gromit

total_predicted_time = (length(Pitch2)*length(Pitch3)*length(Heave3)*length(phase13))*((1/EP.freq)*(EP.num_cyc + 2*EP.transientcycs) + 20)/60/60;

msg = ['Total predicted experiment time:\n', num2str(total_predicted_time),' hours\n'];
fprintf(msg);
fprintf('Press any key to continue (Ctrl+C to Cancel)\n\n')
pause

fprintf('Running main experiment...\n')

%% Run experiment

startexp = tic;
exp = 0;

for q = 1:length(Pitch2)

for j = 1:length(phase13)
    for k = 1:length(Pitch3)
        for i = 1:length(Heave3)
%     alphaT4 = atan(-2*pi*heave2*fred) + deg2rad(pitch2); % self-explanatory
            
            EP.P3 = Pitch3(k); EP.P2 = Pitch2(q);
            EP.H3 = Heave3(i); EP.H2 = Heave2;
            EP.phase13 = phase13(j);
            EP.srate = dq.Rate;
            currentTime = clock;
            
            alphaT4 = atan(-2*pi*Heave2*fred) + deg2rad(Pitch2(q)); % self-explanatory
            EP.alphaT4 = alphaT4;

            [flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2q, heave2, pitch3k, heave3i, phase13j, num_cyc, phi, foiltype]...
                = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
                EP.foiltype, EP.freq, Pitch2(q), Heave2, Pitch3(k), Heave3(i), phase13(j), EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude, 0, 0, offset);
            
%             file = ['\', datestr(now,'yyyymmdd'), '_TandemFoil_APHPhaseSweep_', EP.foiltype,...
            file = ['\', '20221011', '_TandemFoil_APHPhaseSweep_', EP.foiltype,...
                '_alpha=', num2str(alphaT4,3), '_p3=', num2str(Pitch3(k)), '_h3=', num2str(Heave3(i),3), 'c_phase=', num2str(phase13(j)), '.mat'];
            
            [diagnostics,continue_exp] = run_diagnostics(Prof_out_angle,out,dq.Rate,freq);

            filename = fullfile(folderpath, file);
            save(filename);%, 'EP', 'dat', 'Prof_out_angle', 'out', 'foiltype', 'fred', 'freq', 'pitch2q', 'pitch3k', 'heave3i', 'phase13j', 'Re', 'alphaT4', 'num_cyc');

            pause(20)
            exp = exp + 1;
            fprintf(sprintf('\nFinished experiment %i \n\n Running next experiment... \n\n',exp));
        end
    end    
end

end

endexp = toc;

disp(['Done \n Full run took ',num2str(endexp,3),' seconds to complete.'])

%% second experiment

Heave2 = Heave3;
Pitch2 = [65,70,75,80]; %Pitch3;

phase13 = -180;
Heave3 = 0;
Pitch3 = 0;

startexp = tic;
exp = 0;

for i = 1:length(Heave2)
    for j = 1:length(Pitch2)

        EP.P2 = Pitch2; EP.P2 = Pitch2(j);
        EP.H2 = Heave2; EP.H2 = Heave2(i);
        EP.phase13 = phase13;
        EP.srate = dq.Rate;
        currentTime = clock;
        
        alphaT4 = atan(-2*pi*Heave2*fred) + deg2rad(Pitch2(j)); % self-explanatory
        EP.alphaT4 = alphaT4;

        [flume, out, dat, Prof_out_angle, Prof_out, last_out, freq, pitch2j, heave2i, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
            = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,...
            EP.foiltype, EP.freq, Pitch2(j), Heave2(i), Pitch3, Heave3, phase13, EP.phi, EP.num_cyc, EP.transientcycs, EP.constantamplitude, 0, 0, offset);
        
%             file = ['\', datestr(now,'yyyymmdd'), '_TandemFoil_APHPhaseSweep_', EP.foiltype,...
        file = ['\', '20221011', '_SingleFoil_PHSweep_', EP.foiltype,...
            '_p2=', num2str(Pitch2(j)), '_h2=', num2str(Heave2(i),3), '.mat'];
        
        [diagnostics,continue_exp] = run_diagnostics(Prof_out_angle,out,dq.Rate,freq);

        filename = fullfile(folderpath, file);
        save(filename);%, 'EP', 'dat', 'Prof_out_angle', 'out', 'foiltype', 'fred', 'freq', 'pitch2q', 'pitch3k', 'heave3i', 'phase13j', 'Re', 'alphaT4', 'num_cyc');

        pause(20)
        exp = exp + 1;
        fprintf(sprintf('\nFinished experiment %i \n\n Running next experiment... \n\n',exp));

    end
end

endexp = toc;

disp(['Second done \n Full run took ',num2str(endexp,3),' seconds to complete.'])
