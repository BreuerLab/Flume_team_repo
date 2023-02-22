% This script will run a series of trials with specified range of dimensionless frequency and amplitude 
% and automatically save the output data

startexp = tic;
experimentnamestr = ['20230221_TandemTuesday_threeAlphaSweep_A3E_'];
foiltype='A3E';
chord=0.061; % meters
thcknss = 0.00635;
U = 0.33; % m/s
%U = 0.26;
num_cyc = 30; % must be even?
transientcycs = 3;
constantpitch = 0; % 1 for constant pitch during trial, only last foil

%A1pitch = 0; % pitch amplitude of upstream foil in degrees
A1star = 0.8; % heave amplitude of updtream foil in meters
A1 = A1star*chord;

% phase2 = 0;
phi = -90; % phase between heave and pitch
offset = transientcycs; % Time (in cycles) from start of run to start PIV (or LDV trigger)
fstar = 0.12; % reduced frequency

phase2_vec = [-120, -60, 0, 60, 120, 180]; % inter-foil phase [deg]
pitch1_vec = [40, 50, 70]; % leading foil pitch [deg]
pitch2_vec = [65, 70, 75]; % trailing foil pitch [deg]
A2star_vec = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2];%, 1.4, 1.6]; % trailing foil heave [chords]
% heave2_vec = A2star_vec*chord*1.678688;

freq = fstar*U/chord;

exp = 0;

bias_realigned = bias; % for whenever wallace has to realign
heaveGain = 1.678688; % temporary for the new traverse

for A1pitch = pitch1_vec
    for A2pitch = pitch2_vec
        for A2star = A2star_vec
            A2 = A2star*chord*heaveGain; % trailing heave accounting for the empirical gain calculated in the new traverse [m]
            for phase2 = phase2_vec
    
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
                
%                 %% testing
%                 [out,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,[0,0,45,0,45,0],15,bias_trial,foil);
%                 [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq, pitch2, heave2, pitch3, heave3, phase13, num_cyc, phi, foiltype]...
%                 = run_Motors(dq,last_out,bias_trial,foiltype, 1, 0, 0, 0, 0, 0, 0,...
%                 60, 0, 0, 1);
%                 [out,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,[0,0,0,0,0,0],15,bias_trial,foil);
%                 trialname = [fname, '\data\20230221_ReferenceTest_1.mat'];
%                 save(trialname)
%                 %% end of testing

                exp = exp + 1;
                disp(['Running trial at p1=',num2str(A1pitch,2),'deg, p2=',num2str(A2pitch,2),'deg, h2 = ',num2str(A2star,3),'c, phase2 = ',num2str(phase2,4),'c, exp = ', num2str(exp)])
            
                % Runs the function that moves the motors ("run_Motors")
                [flume, out, dat, Prof_out_angle, Prof_out,last_out, freq, pitchNewT, heaveNewT, pitchLead, heaveLead, phase13, num_cyc, phi,...
                    foiltype]...
                = run_Motors(dq,last_out,bias_trial,foiltype, freq, A2pitch, A2, A1pitch, A1, phase2,... % remember Wallace(lead) is 1 and NewTraverse(trail) is 2
                phi, num_cyc, transientcycs, constantpitch, offset);
                
                alphaT4 = atan(-2*pi*A1star*fstar) + deg2rad((A1pitch));
                trialname = [fname,'\data\',experimentnamestr,'_aT4=',num2str(alphaT4,3),'rad,p3=',num2str(A2pitch,2),'deg,h3=',num2str(A2star,3),'c,phase=',num2str(phase2),'.mat']
        %         disp('Trial complete, saving data.')
                exp_timestamp = clock;
                save(trialname)

                % realign the stupid wallace pitch motor wallace
%                 [last_out, bias_realigned] = find_zero_pitch(dq, last_out, bias_trial, foil, 1); % find zero pitch for stupid wallace
                % switched the motors out (got rid of the stupid hudson in favor of the stepper)
                bias_realigned = bias_trial; % added this just so I didn't have to change all the code

%                 % For debugging
%                 singexp = toc(startexp);
%                 disp(['Single experiment took ',num2str(singexp,3),' seconds to complete.']) % temporal

            end
        end
    end
end

endexp = toc(startexp);
disp(['Full run took ',num2str(endexp,3),' seconds to complete.'])
