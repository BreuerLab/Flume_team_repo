% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = 'Vib';
foiltype='V1';
chord=0.0535; % meters
U = 0.2; % m/s
num_cyc = 60; % must be even?
transientcycs = 5;
constantpitch = 0; % 1 for constant pitch during trial, only last foil
A2pitch = 0; % Pitch amplitude in degrees
A1pitch = 0; % pitch amplitude of upstream foil in degrees
A1star = 0; % heave amplitude of upstream foil in meters
phase2 = 0;
phi = 0;
offset = 0; % Time (in cycles) from start of run to start PIV

for fstar = 0.3:0.02:0.3 %0.3:0.02:0.3  
        % fstar = 0.3 % bias drift test constant freq
        freq = fstar*U/thcknss;

    for A2star = 0:0.05:0 %0.0:0.1:1.1
        A2 = A2star*thcknss;
        A1 = A1star*thcknss;

        % Used to specify heave velocity and acceleration limits
        heavevelocommandmax = A2*2*pi*freq;
        heaveaccelcommandmax = A2*(2*pi*freq)^2;
        if heavevelocommandmax > 0.50 % m/s
            disp('Commanded velocity limit exceeded, skipping this trial')
            break
        elseif heaveaccelcommandmax > 3.5 % m/s^2
            disp('Commanded acceleration limit exceeded, skipping this trial')
            break
        elseif A1 > 0.12 || A2 > 0.12 % meters
            disp('Commanded heave limit exceeded, skipping this trial')
            break
        end
        
        disp(['Running trial at f=',num2str(freq,3),'Hz and A=',num2str(A2*100,3),'cm'])
        
        % Runs the function that moves the motors ("run_Motors")
        [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq,pitch2, heave2, pitch3, heave3,phase13, num_cyc, phi,...
            foiltype]...
        = run_Motors(dq,last_out,bias,foiltype, freq, A1pitch, A1, A2pitch, A2, phase2,...
        phi, num_cyc, transientcycs, constantpitch, offset);
        
        trialname = [fname,'\data\',experimentnamestr,'_pitch=',num2str(A2pitch,3),'deg,f=',num2str(freq,3),'Hz,A=',num2str(100*A2,3),'cm.mat'];
%         disp('Trial complete, saving data.')
        save(trialname)
    end
end
endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])