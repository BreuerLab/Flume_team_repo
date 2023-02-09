% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = '20221116_SingleWednesday_singleAlphaSweep_A2_Extra_';
foiltype='A2';
chord=0.0762; % meters
thcknss = 0.00635;
% U = 0.33; % m/s
U = 0.26;
num_cyc = 20; % must be even?
transientcycs = 5;
constantpitch = 0; % 1 for constant pitch during trial, only last foil

% A2pitch = 0; % Pitch amplitude in degrees
A2pitch = 0; % pitch amplitude of upstream foil in degrees
A2star = 0; % heave amplitude of upstream foil in meters
A2 = A2star*chord;

phase2 = 0;
phi = -90;
offset = 0; % Time (in cycles) from start of run to start PIV
% fstar = 0.12;

fstar_vec = [0.1,0.12,0.15];
pitch1_vec = [40,45,50,55,60,65,70,75,80]; % wallace
heave1_vec = [0.5,0.75,1,1.25,1.5]; % wallace

freq_vec = fstar_vec*U/chord;

exp = 0;

bias_realigned = bias;

for A1pitch = pitch1_vec
    for A1star = heave1_vec
        A1 = A1star*chord;
        for freq = freq_vec
            fstar = freq*chord/U;

            % Used to specify heave velocity and acceleration limits
            heavevelocommandmax = A1*2*pi*freq;
            heaveaccelcommandmax = A1*(2*pi*freq)^2;
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

            % Take another force sensor tare measurement right before the trial starts
            [~,bias_newloaded,~] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil);
            % bias_trial -> bias_new - bias_loaded + bias
            bias_trial.Wallace = bias_newloaded.Wallace - bias_loaded.Wallace + bias.Wallace;
            bias_trial.Gromit = bias_newloaded.Gromit - bias_loaded.Gromit + bias.Gromit;
            bias_trial.accmeter = bias_newloaded.accmeter - bias_loaded.accmeter + bias.accmeter;
            bias_trial.pitch = bias_realigned.pitch;
            
            exp = exp + 1;
            disp(['Running trial at f=',num2str(freq,3),'Hz and A=',num2str(A1*100,3),'cm, p2 = ',num2str(A1pitch,2),'deg, exp = ', num2str(exp)])
        
            % Runs the function that moves the motors ("run_Motors")
            [flume, out, dat, Prof_out_angle, Prof_out,last_out, freqi, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi,...
                foiltype]...
            = run_Motors(dq,last_out,bias_trial,foiltype, freq, A1pitch, A1, A2pitch, A2, phase2,...
            phi, num_cyc, transientcycs, constantpitch, offset);
            
            alphaT4 = atan(-2*pi*A1star*fstar) + deg2rad((A1pitch));
            trialname = [fname,'\data\',experimentnamestr,'_aT4=',num2str(alphaT4,3),'rad,p3=',num2str(A1pitch,2),'deg,h3=',num2str(A1star,3),'c,fstar=',num2str(fstar,3),'.mat']
    %         disp('Trial complete, saving data.')
            save(trialname)

%             % realign the stupid wallace pitch motor wallace
%             out_of_the_way = [0,0,0,0.15,0,0]; % position of gromit to move out of the way
%             [~,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,out_of_the_way,5,bias_trial,foil); % move gromit out of the way
%             [last_out,bias_realigned] = find_zero_pitch_wallace(dq,last_out,bias_trial,foil,out_of_the_way); % find zero pitch for wallace
%             [~,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,[0,0,0,0,0,0],5,bias_realigned,foil); % move gromit to original position
%             % continue
            
        end
    end
end

endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])
