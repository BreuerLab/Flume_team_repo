% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = 'vibPIV';
foiltype='V1';
chord=0.024; % meters
U = 0.2; % m/s
num_cyc = 50; % must be even?
transientcycs = 5;
constantpitch = 0; % 1 for constant pitch during trial, only last foil
A2pitch = 0; % Pitch amplitude in degrees
A1pitch = 0; % pitch amplitude of upstream foil in degrees
A1star = 0; % heave amplitude of upstream foil in meters
phase2 = 0;
phi = 0;
offset = 10; % Time (in cycles) from start of run to start PIV

for fstar = 0.24:0.02:0.24 %0.06:0.02:0.24  
        freq = fstar*U/chord;

    for A2star = 0:0.1:0 %0.0:0.1:1.1
        A2 = A2star*chord*100;

        % Used to specify heave velocity and acceleration limits
        heavevelocommandmax = A2*2*pi*freq/100;
        heaveaccelcommandmax = A2*(2*pi*freq)^2/100;
        if heavevelocommandmax > 0.50 || heaveaccelcommandmax > 4.9 
            break
        end
        
        disp(['Running trial at f=',num2str(freq,3),'Hz and A=',num2str(A2,3),'cm'])
        
        % Runs the function that moves the motors ("run_Motors")
        [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq,pitch2, heave2, pitch3, heave3,phase13, num_cyc, phi,...
            foiltype]...
        = run_Motors(dq,last_out,pitch_bias,Wbias,Gbias,accbias,foiltype, freq, A1pitch, A1star, A2pitch, A2star, phase2,...
        phi, num_cyc, transientcycs, constantpitch, offset);
        
        trialname = [fname,'\data\',experimentnamestr,'_pitch=',num2str(A2pitch,3),'deg,f=',num2str(freq,3),'Hz,A=',num2str(A2,3),'cm.mat'];
%         disp('Trial complete, saving data.')
        save(trialname)
    end
end
endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])