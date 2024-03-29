% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = 'EllipticCyl';
foiltype='EC1';
chord=0.0594; % meters
thcknss = 0.0238;
U = 0.3; % m/s
num_cyc = 60; % must be even?
transientcycs = 5;
constantpitch = 0; % 1 for constant pitch during trial, only last foil
A2pitch = 0; % Pitch amplitude in degrees
A1pitch = 30; % pitch amplitude of upstream foil in degrees
% A1star = 0; % heave amplitude of upstream foil in meters
A1 = 0.024;
freq1 = 0.15*U/thcknss; 
phase2 = -148.45;
phi = -90;
offset = 0; % Time (in cycles) from start of run to start PIV
 
for fstar = 0.1:0.01:0.3 %0.1:0.01:0.3
        freq2 = fstar*U/thcknss;

    for A2star = 0:0.05:1.1 %0.0:0.1:1.1
        A2 = A2star*thcknss;
        A1 = A1star*thcknss;

        % Used to specify heave velocity and acceleration limits
        heavevelocommandmax1 = A1*2*pi*freq;
        heavevelocommandmax2 = A2*2*pi*freq;
        heaveaccelcommandmax1 = A1*(2*pi*freq)^2;
        heaveaccelcommandmax2 = A2*(2*pi*freq)^2;
        if heavevelocommandmax1 > 0.5 || heavevelocommandmax2 > 0.50 % m/s
            disp('Commanded velocity limit exceeded, skipping this trial')
            break
        elseif heaveaccelcommandmax1 > 3.5 || heaveaccelcommandmax2 > 3.5 % m/s^2
            disp('Commanded acceleration limit exceeded, skipping this trial')
            break
        elseif A1 > 0.12 || A2 > 0.12 % meters
            disp('Commanded heave limit exceeded, skipping this trial')
            break
        end
        
        % Take another force sensor tare measurement right before the trial starts
        [~,bias_newloaded,~] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil);
        % bias_trial -> bias_new - bias_loaded + bias
        bias_trial.Wallace = bias_newloaded.Wallace - bias_loaded.Wallace + bias.Wallace;
        bias_trial.Gromit = bias_newloaded.Gromit - bias_loaded.Gromit + bias.Gromit;
        bias_trial.accmeter = bias_newloaded.accmeter - bias_loaded.accmeter + bias.accmeter;
        bias_trial.pitch = bias.pitch;

        disp(['Running trial at f=',num2str(freq,3),'Hz and A=',num2str(A2*100,3),'cm'])
        
        % Runs the function that moves the motors ("run_Motors")
        [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq,pitch2, heave2, pitch3, heave3,phase13, num_cyc, phi,...
            foiltype]...
        = run_Motors(dq,last_out,bias_trial,foiltype, freq, A1pitch, A1, A2pitch, A2, phase2,...
        phi, num_cyc, transientcycs, constantpitch, offset);
        
        trialname = [fname,'\data\',experimentnamestr,'_pitch=',num2str(A2pitch,3),'deg,f=',num2str(freq,3),'Hz,A=',num2str(100*A2,3),'cm.mat'];
%         disp('Trial complete, saving data.')
        save(trialname)
    end
end
endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])