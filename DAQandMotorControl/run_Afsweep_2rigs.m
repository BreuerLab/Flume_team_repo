% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = '20221113_TandemSunday_leadingAlphaSweep_aT4=0.68';
foiltype='A3E';
chord=0.061; % meters
thcknss = 0.00635;
U = 0.33; % m/s
num_cyc = 20; % must be even?
transientcycs = 3;
constantpitch = 0; % 1 for constant pitch during trial, only last foil
% A2pitch = 0; % Pitch amplitude in degrees
% A1pitch = 0; % pitch amplitude of upstream foil in degrees
A1star = 0.8; % heave amplitude of upstream foil in meters
% phase2 = 0;
phi = -90;
offset = 0; % Time (in cycles) from start of run to start PIV
fstar = 0.12;
freq = fstar*U/chord;

pitch2_vec = [40,50,70];
pitch3_vec = [65,70,75];
heave3_vec = 0.7:0.1:1.2;
phase12_vec = [-180,-120,-60,0,60,120];
A1 = A1star*chord;

exp = 0;

bias_realigned = bias;

for A1pitch = pitch2_vec %0.3:0.02:0.3
    for A2pitch = pitch3_vec
        for A2star = heave3_vec
            A2 = A2star*chord;
            for phase2 = phase12_vec

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

                % Take another force sensor tare measurement right before the trial starts
                [~,bias_newloaded,~] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil);
                % bias_trial -> bias_new - bias_loaded + bias
                bias_trial.Wallace = bias_newloaded.Wallace - bias_loaded.Wallace + bias.Wallace;
                bias_trial.Gromit = bias_newloaded.Gromit - bias_loaded.Gromit + bias.Gromit;
                bias_trial.accmeter = bias_newloaded.accmeter - bias_loaded.accmeter + bias.accmeter;
                bias_trial.pitch = bias_realigned.pitch;
                
                exp = exp + 1;
                disp(['Running trial at f=',num2str(freq,3),'Hz and A=',num2str(A2*100,3),'cm, exp = ', num2str(exp)])
            
                % Runs the function that moves the motors ("run_Motors")
                [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq,pitch2, heave2, pitch3, heave3,phase13, num_cyc, phi,...
                    foiltype]...
                = run_Motors(dq,last_out,bias_trial,foiltype, freq, A1pitch, A1, A2pitch, A2, phase2,...
                phi, num_cyc, transientcycs, constantpitch, offset);
                
                alphaT4 = atan(-2*pi*A1star*fstar) + deg2rad((A1pitch));
                trialname = [fname,'\data\',experimentnamestr,'_aT4=',num2str(alphaT4,3),'rad,p3=',num2str(A2pitch,2),'deg,h3=',num2str(A2star,3),'c,ph=',num2str(phase2,3),'.mat']
        %         disp('Trial complete, saving data.')
                save(trialname)

                % realign the stupid wallace pitch motor wallace
                out_of_the_way = [0,0,0,0.15,0,0]; % position of gromit to move out of the way
                [~,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,out_of_the_way,5,bias_trial,foil); % move gromit out of the way
                [last_out,bias_realigned] = find_zero_pitch_wallace(dq,last_out,bias_trial,foil,out_of_the_way); % find zero pitch for wallace
                [~,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,[0,0,0,0,0,0],5,bias_realigned,foil); % move gromit to original position
                % continue
                
            end
        end
    end
end
endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])
