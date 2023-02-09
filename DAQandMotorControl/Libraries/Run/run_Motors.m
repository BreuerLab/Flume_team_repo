function [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq, pitch2, heave2, pitch3, heave3, pitch4, heave4, phase13, num_cyc, phi, foiltype]...
    = run_Motors(dq,last_out,bias,foiltype, freq, pitch2, heave2, pitch3, heave3, pitch4, heave4, phase13, phi,...
    num_cyc, transientcycs, constantamplitude, offset)
%%
% Given frequency [Hz], Pitch amplitude [deg] and heave amplitude [chords], this function will run 3 rigs for a set number of cycles
% NOTE: Considering the frontmost rig was removed (Shawn), it's values are defaulted to 0
%
% > freq --- oscillation frequency of the foils (assuming both oscillate at the same frequency)
% > foiltype --- selected from the function foils_database (A1, A2, E1, etc.) [char]
%
% > pitch1 and heave1 --- Shawn
% > pitch2 and heave2 --- Wallace
% > pitch3 and heave3 --- Gromit
%
% > phase12 --- phase difference between rig 1 and rig 2 [deg]
% > phase13 --- phase difference between rig 1 and rig 3 [deg]
% > num_cyc --- number of cycles for the motion profile (must be an even number, annoyingly so)
% > phi --- phase between heave and pitch [deg]

%% Initial setup
freq1 = freq; % all foils will oscillate at the same frequency
freq2 = freq; % freq
freq3 = freq;
freq4 = freq;

pitch1 = 0; % these three are 0 because Shawn (frontmost rig) is no more
heave1 = 0; % changed this temporarily to test out the measurements
phase12 = 0;

phase12r = deg2rad(phase12); % this is technically unnecessary, as it doesn't really go anywhere
phase13r = deg2rad(phase13);

if num_cyc < 10
    num_cyc = 10;
    disp('Minimum of 10 cycles for accurate data collection')
end

[foil, ~, ~] = foils_database(foiltype); % basically overwrites the values initially defined in the daq_setup
chord = foil.chord;

if isempty(dq)
    error('Run setup_DAQ')
end

write(dq,last_out)

% heave1 = heave1*thcknss; % multiply the heave amplitude times the chord
% heave3 = heave3*thcknss; % the value for the chord comes from the "daq_setup_3rigs", where it's taken from "foils_database"
% heave2 = heave2*thcknss;

% I CAN DELETE THIS
% params = [freq1, pitch1, heave1, phase12r, phase13r, 90; %Shawn (first)
%           freq1, pitch2, heave2, phase12r, phase13r, phi; %Wallace (last) % the order of this might be wrong
%           freq1, pitch3, heave3, phase12r, phase13r, phi]; %Gromit (mid)
% I CAN DELETE THIS
% last_pos = conv_last_out(last_out,bias.pitch); % dunno what this is for

%% Profile generation
% Pitch and heave values have to be taken as negative in order to have the
% motion produced in the correct frame of reference relative to the flow
% direction.

% all profiles must have same number of timesteps, even if freq are
% different
num_cyc2 = num_cyc*freq2/freq3;
transientcycs2 = transientcycs*freq2/freq3;

[t1p, Prof1p] = generate_profile(num_cyc, freq1, dq.Rate, transientcycs, transientcycs, pitch1, phi,0);          % pitch Shawn
[t1h, Prof1h] = generate_profile(num_cyc, freq1, dq.Rate, transientcycs, transientcycs, heave1, 0,0);            % heave Shawn
[t2p, Prof2p] = generate_profile(num_cyc2, freq2, dq.Rate, transientcycs2, transientcycs2, pitch2, phase12 + phi,0);% pitch Gromit
[t2h, Prof2h] = generate_profile(num_cyc2, freq2, dq.Rate, transientcycs2, transientcycs2, heave2, phase12,0);      % heave Gromit
[t3p, Prof3p] = generate_profile(num_cyc, freq3, dq.Rate, transientcycs, transientcycs, pitch3, phase12 + phi,0);% pitch Wallace
[t3h, Prof3h] = generate_profile(num_cyc, freq3, dq.Rate, transientcycs, transientcycs, heave3, phase12,0);      % heave Wallace

[t4p, Prof4p] = generate_profile(num_cyc, freq4, dq.Rate, transientcycs, transientcycs, pitch4, phase13 + phi,0);% pitch traverse
[t4h, Prof4h] = generate_profile(num_cyc, freq4, dq.Rate, transientcycs, transientcycs, heave4*chord, phase13,0);      % heave traverse


Prof_out_angle = [Prof1p, Prof1h, Prof2p, Prof2h, Prof3p, Prof3h];
Prof_out_temp = input_conv_3rigs(Prof_out_angle, freq, heave1, heave2, heave3,bias.pitch); % let's hope this works <-- It does! But something is weird with the signs.

Prof_out_angle = [Prof_out_angle, Prof4h, Prof4p]; % add the traverse conversion

% For PIV trigger

trig_signal=ones(length(Prof_out_temp(:,1)),1);
trig_signal(1:round(1/freq*dq.Rate*offset))=0;
trig_signal(end-round(1/freq*dq.Rate)*1:end)=0;

% For the traverse signal

% Traverse heave: (motion profile and lock)
[voltageCmd_heave, cmdLock_heave] = traversecmd('y', Prof4h, 0.25); % position rig starting at the center of the flume
% Traverse heave: (motion profile and lock)
[voltageCmd_pitch, cmdLock_pitch] = traversecmd('theta', Prof4p, 180); % position rig starting facing the upstream

traverse_signal = [voltageCmd_heave, cmdLock_heave, voltageCmd_pitch, cmdLock_pitch]; % assemble a matrix to append to the output variable

% Append trigger and traverse signals to the outputted signal
Prof_out = [Prof_out_temp trig_signal traverse_signal];

%% Run section
tic
dat = readwrite(dq,Prof_out,"OutputFormat","Matrix");
last_out = Prof_out(end,:);
write(dq,last_out);
scantime = toc;
disp(['Time taken to run DAQ foreground scan: ',num2str(scantime),' seconds.'])


% Convert raw voltages to useful data
[out,t]=output_conv_3rigs(dat,bias,foil); % output angle is in radians
% Convert new traverse

Length = numel(out(:,1))/1000;
flume = 0;

end
