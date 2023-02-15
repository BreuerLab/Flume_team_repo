function [dq] = setup_DAQ_channels(samplerate)

disp('Inititializing NI DAQs')

% Channels in s are as follows: THIS NEEDS TO BE UPDATED
% 
% Data acquisition session using National Instruments hardware:
%    No data queued.  Will run at 1000 scans/second.
%    Number of channels: 28
%       index Type Device Channel   MeasurementType        Range         Name   
%       ----- ---- ------ ------- ------------------- ---------------- ---------
%       1     ci   Dev1   ctr3    Position            n/a              Pitch 1
%       2     ci   Dev1   ctr2    Position            n/a              Heave 1
%       3     ci   Dev1   ctr1    Position            n/a              Pitch 2
%       4     ci   Dev1   ctr0    Position            n/a              Heave 2
%       5     ci   Dev2   ctr0    Position            n/a              Pitch 3
%       6     ci   Dev2   ctr1    Position            n/a              Heave 3
%       7     ai   Dev1   ai0     Voltage (Diff)      -10 to +10 Volts Wallace 1
%       8     ai   Dev1   ai1     Voltage (Diff)      -10 to +10 Volts Wallace 2
%       9     ai   Dev1   ai2     Voltage (Diff)      -10 to +10 Volts Wallace 3
%       10    ai   Dev1   ai3     Voltage (Diff)      -10 to +10 Volts Wallace 4
%       11    ai   Dev1   ai4     Voltage (Diff)      -10 to +10 Volts Wallace 5
%       12    ai   Dev1   ai5     Voltage (Diff)      -10 to +10 Volts Wallace 6
%       13    ai   Dev1   ai6     Voltage (SingleEnd) -10 to +10 Volts Vel_x
%       14    ai   Dev1   ai14    Voltage (SingleEnd) -10 to +10 Volts Vel_y
%       15    ai   Dev1   ai7     Voltage (SingleEnd) -10 to +10 Volts Vel_z1
%       16    ai   Dev1   ai15    Voltage (SingleEnd) -10 to +10 Volts Vel_z2
%       17    ai   Dev1   ai16    Voltage (Diff)      -10 to +10 Volts Gromit 1
%       18    ai   Dev1   ai17    Voltage (Diff)      -10 to +10 Volts Gromit 2
%       19    ai   Dev1   ai18    Voltage (Diff)      -10 to +10 Volts Gromit 3
%       20    ai   Dev1   ai19    Voltage (Diff)      -10 to +10 Volts Gromit 4
%       21    ai   Dev1   ai20    Voltage (Diff)      -10 to +10 Volts Gromit 5
%       22    ai   Dev1   ai21    Voltage (Diff)      -10 to +10 Volts Gromit 6
%       23    ao   Dev1   ao2     Voltage (SingleEnd) -10 to +10 Volts Pitch 1
%       24    ao   Dev1   ao0     Voltage (SingleEnd) -10 to +10 Volts Heave 1
%       25    ao   Dev1   ao3     Voltage (SingleEnd) -10 to +10 Volts Pitch 2
%       26    ao   Dev1   ao1     Voltage (SingleEnd) -10 to +10 Volts Heave 2
%       27    ao   Dev2   ao0     Voltage (SingleEnd) -10 to +10 Volts Pitch 3
%       28    ao   Dev2   ao1     Voltage (SingleEnd) -10 to +10 Volts Heave 3


%       ##    ai   Dev1   ai0     Voltage (SingleEnd) 0 to +10 Volts   cmd recording Heave New Traverse
%       ##    ai   Dev1   ai1     Voltage (SingleEnd) 0 to +10 Volts   cmd recording Pitch New Traverse
%       ##    ci   Dev1   ctr0    Position            n/a              Heave New Traverse
%       ##    ci   Dev1   crt1    Position            n/a              Pitch New Traverse

%       ##    ao   Dev1   ao0     Voltage             0 to +10 Volts   command Heave New Traverse
%       ##    di/o Dev1   P0.0    Digital             0 or +5 Volts    cmd Lock Heave New Traverse
%       ##    ao   Dev1   ao1     Voltage             0 to +10 Volts   command Pitch New Traverse
%       ##    di/o Dev1   P0.1    Digital             0 or +5 Volts    cmd Lock Pitch New Traverse

dq = daq("ni");
% dq=daq.createSession('ni');
dq.Rate = samplerate; % sampling rate (same as fs in foils_database)
T = 1/dq.Rate;


%% Counter channels for encoder inputs
% global ch1
ch1=addinput(dq,"Dev2","ctr2","Position");
% ch1=dq.addCounterInputChannel('dev3','ctr0','Position');
ch1.EncoderType='X4';
ch1.ZResetEnable=0;
ch1.Name = 'Pitch Shawn';
ch1.ZResetCondition = 'BothLow';
% disp(['Channel 1 Terminal A: ', ch1.TerminalA])
% disp(['Channel 1 Terminal B: ', ch1.TerminalB])
% disp(['Channel 1 Terminal Z: ', ch1.TerminalZ])
ch1.ZResetValue = -65;
% global ch2
ch2=addinput(dq,"Dev2","ctr3","Position");
% ch2=dq.addCounterInputChannel('dev3','ctr1','Position');
ch2.EncoderType='X4';
ch2.ZResetEnable=0;
ch2.ZResetCondition = 'BothLow';
ch2.ZResetValue = 0;
ch2.Name = 'Heave Shawn';
% disp(['Channel 2 Terminal A: ', ch2.TerminalA])
% disp(['Channel 2 Terminal B: ', ch2.TerminalB]) 
% disp(['Channel 2 Terminal Z: ', ch2.TerminalZ])

% REPLACE GROMIT ENCODER INPUTS WITH NEW TRAVERSE
% % global ch3
% ch3=addinput(dq,"Dev2","ctr3","Position");
% % ch3=s.addCounterInputChannel('Dev2','ctr3','Position');
% ch3.EncoderType='X4';
% ch3.ZResetEnable=0;
% ch3.ZResetCondition = 'BothLow';
% ch3.Name = 'Pitch Gromit';
% % disp(['Channel 3 Terminal A: ', ch3.TerminalA])
% % disp(['Channel 3 Terminal B: ', ch3.TerminalB])
% % disp(['Channel 3 Terminal Z: ', ch3.TerminalZ])
% ch3.ZResetValue = -795;
% % global ch4
% ch4=addinput(dq,"Dev2","ctr2","Position");
% % ch4=dq.addCounterInputChannel('dev2','ctr2','Position');
% ch4.EncoderType='X4';
% ch4.ZResetEnable=0;
% ch4.ZResetCondition = 'BothLow';
% ch4.ZResetValue = 0;
% ch4.Name = 'Heave Gromit';
% % disp(['Channel 4 Terminal A: ', ch4.TerminalA])
% % disp(['Channel 4 Terminal B: ', ch4.TerminalB])
% % disp(['Channel 4 Terminal Z: ', ch4.TerminalZ])

% new traverse:
ctr_theta = addinput(dq, 'Dev3', 'ctr1', 'Position');
ctr_theta.Name = 'encoder_theta';
ctr_theta.EncoderType = 'X4';

ctr_y = addinput(dq, 'Dev3', 'ctr0', 'Position');
ctr_y.Name = 'encoder_y';
ctr_y.EncoderType = 'X4';

% global ch5
ch5=addinput(dq,"Dev2","ctr1","Position");
% ch5=dq.addCounterInputChannel('Dev2','ctr1','Position');
ch5.EncoderType='X4';
ch5.ZResetEnable=0;
ch5.ZResetCondition = 'BothLow';
ch5.Name = 'Pitch Wallace';
% ch5.ZResetValue = -104;
ch5.ZResetValue = 338;
% global ch6
ch6=addinput(dq,"Dev2","ctr0","Position");
% ch6=dq.addCounterInputChannel('Dev2','ctr0','Position');
ch6.EncoderType='X4';
ch6.ZResetEnable=0;
ch6.ZResetCondition = 'BothLow';
ch6.ZResetValue = 0;
ch6.Name = 'Heave Wallace';


% s.addTriggerConnection('Dev1/PFI12','Dev4/PFI0','StartTrigger');
% s.addClockConnection('Dev1/PFI14','Dev4/PFI14', 'ScanClock');
disp('Counters done.')

%% Analog input channels
global chins1 chins2 chins3 chins4 chins6
chins1=addinput(dq,'dev2',[0 1 2 3 4 5],'Voltage'); % Wallace force sensor
chins2=addinput(dq,'dev3',[6 7 8 9],'Voltage'); % Vectrino
chins3=addinput(dq,'dev3',[0 1 2 3 4 5],'Voltage'); % Gromit force sensor
chins6=addinput(dq,'dev2',[24],'Voltage'); % accelerometer (16 dev2.1 + 8 dev2.2 = 24th channel)


% chins1=dq.addAnalogInputChannel('dev2',[0 1 2 3 4 5],'Voltage'); % Wallace force sensor
% chins2=dq.addAnalogInputChannel('dev3',[6 7 8 9],'Voltage'); % Vectrino
% % chins3=s.addAnalogInputChannel('dev2',[16 17 18 19 20 21],'Voltage');
% chins3=dq.addAnalogInputChannel('dev3',[0 1 2 3 4 5],'Voltage'); % Gromit force sensor
% chins6=dq.addAnalogInputChannel('dev2',[24],'Voltage'); % accelerometer (16 dev2.1 + 8 dev2.2 = 24th channel)
% % chins4 = s.addAnalogInputChannel('dev3',[2],'Voltage'); 

    
chins1(1).Name = 'Wallace 1';
chins1(2).Name = 'Wallace 2';
chins1(3).Name = 'Wallace 3';
chins1(4).Name = 'Wallace 4';
chins1(5).Name = 'Wallace 5';
chins1(6).Name = 'Wallace 6';


chins2(1).TerminalConfig='SingleEnded';
chins2(1).Name = 'Vel_x';
chins2(2).TerminalConfig='SingleEnded';
chins2(2).Name = 'Vel_y';
chins2(3).TerminalConfig='SingleEnded';
chins2(3).Name = 'Vel_z1';
chins2(4).TerminalConfig='SingleEnded';
chins2(4).Name = 'Vel_z2';


chins3(1).Name = 'Gromit 1';
chins3(2).Name = 'Gromit 2';
chins3(3).Name = 'Gromit 3';
chins3(4).Name = 'Gromit 4';
chins3(5).Name = 'Gromit 5';
chins3(6).Name = 'Gromit 6';

chins4(1).Name = 'Trigger';

chins6(1).Name = 'Accelerometer';
chins6(1).TerminalConfig='SingleEnded';

disp('Analog inputs done.')


%% Motor control analog output channels
chout1 = addoutput(dq,'dev3','ao0','Voltage');  % Pitch 1   Shawn
chout2 = addoutput(dq,'dev3','ao1','Voltage');  % Heave 1   Shawn
chout3 = addoutput(dq,'dev2','ao2','Voltage');  % Pitch    Gromit
chout4 = addoutput(dq,'dev2','ao3','Voltage');  % Heave   Gromit
chout5 = addoutput(dq,'dev2','ao0','Voltage');  % Pitch  Wallace
chout6 = addoutput(dq,'dev2','ao1','Voltage');  % Heave  Wallace

% global chout1 chout2 chout3 chout4 chout5 chout6 %chout7 chout8
% chout1 = dq.addAnalogOutputChannel('dev3','ao0','Voltage');  % Pitch 1   Shawn
% chout2 = dq.addAnalogOutputChannel('dev3','ao1','Voltage');  % Heave 1   Shawn
% chout3 = dq.addAnalogOutputChannel('dev2','ao2','Voltage');  % Pitch    Gromit
% chout4 = dq.addAnalogOutputChannel('dev2','ao3','Voltage');  % Heave   Gromit
% chout5 = dq.addAnalogOutputChannel('dev2','ao0','Voltage');  % Pitch  Wallace
% chout6 = dq.addAnalogOutputChannel('dev2','ao1','Voltage');  % Heave  Wallace
% % chout7 = s.addAnalogOutputChannel('dev2','ao0','Voltage');  % 
% % chout8 = s.addAnalogOutputChannel('dev3','ao0','Voltage');  % 

chout1.Name = 'Pitch Shawn';
chout2.Name = 'Heave Shawn';
chout3.Name = 'Pitch Gromit';
chout4.Name = 'Heave Gromit'; 
chout5.Name = 'Pitch Wallace';
chout6.Name = 'Heave Wallace';
disp('Analog outputs done. Syncing and Zeroing output...')

%  PIV trigger and pulse train channels
    addoutput(dq,'dev2','Port0/line14','Digital');% output a trigger signal to the PTU, Dev2.2, P0.5
    addinput(dq,'dev2','Port0/line10','Digital'); % record the trigger signal, Dev2.2, P0.1
    addinput(dq,'dev2','Port0/line9','Digital'); % record the pulse signal from the PTU, Dev2.2, P0.0

% New traverse input channels 20230206 - courtesy of Xiaowei He

ai_y_cmd = addinput(dq, 'Dev2', 20, 'Voltage');
ai_y_cmd.Name = 'y_cmd_m';
ai_y_cmd.TerminalConfig = 'SingleEnded';

ai_theta_cmd = addinput(dq, 'Dev2', 21, 'Voltage');
ai_theta_cmd.Name = 'theta_cmd_m';
ai_theta_cmd.TerminalConfig = 'SingleEnded';

% New traverse output channels 20230206 - courtesy of Xiaowei He
% This was removed in place of using the original Gromit channels

% ao_y_cmd = addoutput(dq, 'Dev1', 'ao0', 'Voltage');
% ao_y_cmd.Name = 'y_cmd';
% 
% dio_y_lock = addoutput(dq, 'Dev1', 'Port0/Line0', 'Digital');
% dio_y_lock.Name = 'y_lock';
% 
% ao_theta_cmd = addoutput(dq, 'Dev1', 'ao1', 'Voltage');
% ao_theta_cmd.Name = 'theta_cmd';
% 
% dio_theta_lock = addoutput(dq, 'Dev1', 'Port0/Line1', 'Digital');
% dio_theta_lock.Name = 'theta_lock';


% Don't know if this needs to be here

% addTriggerConnection(s,'Dev1/PFI4','Dev4/PFI0','StartTrigger');
% 
% s.addClockConnection('Dev1/PFI5','Dev4/PFI1','ScanClock');

% addTriggerConnection(s,'/Dev2/RTSI0','/Dev1/RTSI0','StartTrigger');
% addTriggerConnection(s,'/Dev2/RTSI0','/Dev3/RTSI0','StartTrigger');
% addClockConnection(s,'/Dev2/RTSI1','/Dev1/RTSI1','ScanClock');
% addClockConnection(s,'/Dev2/RTSI1','/Dev3/RTSI1','ScanClock');
% % 
% addTriggerConnection(s,'/Dev3/20MHzTimebase','/Dev1/PFI6','StartTrigger');
% addTriggerConnection(s,'/Dev3/PFI14','/Dev2/PFI14','StartTrigger');
% addClockConnection(s,'/Dev3/PFI15','/Dev1/PFI5','ScanClock');
% addClockConnection(s,'/Dev3/PFI15','/Dev2/PFI15','ScanClock');
% % 
% t1 = addTriggerConnection(s,'/Dev3/PFI14','/Dev1/PFI6','StartTrigger');
% t2 = addTriggerConnection(s,'/Dev3/PFI14','/Dev2/PFI14','StartTrigger');
% c1 = addClockConnection(s,'/Dev1/PFI5','/Dev2/PFI15','ScanClock');
% c2 = addClockConnection(s,'/Dev1/PFI5','/Dev3/PFI15','ScanClock');

end
