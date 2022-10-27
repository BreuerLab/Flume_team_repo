% This code will connect the NI DAQ devices, establish data channels,
% tare the force sensors and accelerometer, and find zero pitch angle
% relative to flume flow

samplerate = 1000; % DAQ sample rate in measurements/second

%% Experimental Setup
% default values:
chord = 0.0535;
thcknss = 0.0265;
span = 0.401;
foil_shape = 'V1';
Wall_distance_left = 0.4;
Wall_distance_right = 0.4;
flume_height = 0.53;
flume_hertz = 10.7;
Number_of_foils = 1;
foil_separation = 0; 
foil_offset = 0;
Temperature = 21.37;
pitch_axis = 0.5;
piv_var = 0;
filt_var = 0;
exp_name = 'Enter descriptive name';
defaultanswers = {num2str(chord),num2str(span),foil_shape,num2str(Wall_distance_left),num2str(Wall_distance_right),...
    num2str(flume_height),num2str(flume_hertz),num2str(Number_of_foils),num2str(foil_separation),num2str(foil_offset),...
    num2str(Temperature),num2str(pitch_axis),num2str(piv_var),num2str(filt_var),exp_name};

prompt = {'Enter chord size (in meters): ','Enter span (in meters): ','Enter foil shapes (as string): ', ...
    'Enter Mean wall distance (left, in meters): ','Enter Mean wall distance (right, in meters): ', ...
    'Enter Flume water height(in meters): ', 'Enter anticipated flume frequency (Hz): ', ... 
    'Enter number of foils in experiment: ','Enter foil separation distance (m): ','Enter foil offset distance (m): ', ...
    'Enter flume water temperature (from vectrino, in c):   ','Enter foil Pitch Axis','Using PIV? (enter 1 for yes)', ... 
    'Filter data at 60Hz? (0 or 1)','Enter Experiment name (folder Name):'};
name = 'Experiment Configuration';
num_lines = 1; 
answer = inputdlg(prompt,name,num_lines,defaultanswers);

chord = str2double(answer{1});
span = str2double(answer{2});
foil_shape = char(answer(3));
Wall_distance_left = str2double(answer{4});
Wall_distance_right = str2double(answer{5});
flume_height = str2double(answer{6});
flume_hertz = str2double(answer{7});
Number_of_foils = str2double(answer{8});
foil_separation = str2double(answer{9});
foil_offset = str2double(answer{10});
Temperature = str2double(answer{11});
pitch_axis = str2double(answer{12});
piv_var = str2double(answer{13});
filt_var = str2double(answer{14});
fname = ['D:\Experiments\',num2str(Number_of_foils),'foil\',answer{15}];
exp_name = answer{15};
Date  = date;
Time = clock;


%  Establish DAQ channels to use
addpath(genpath("Libraries"))
[dq] = setup_DAQ_channels(samplerate);

[foil, rho, fs] = foils_database(foil_shape);

if isfolder(fname)
    disp('Warning: experiment name already exists.  Apppending date and time')
    fname = [fname,'_',Date,'_',num2str(Time(4)),'_',num2str(Time(5)),'_',num2str(round(Time(6)))];
    if isfolder(fname)
        disp('Warning: folder name still taken. Appending time')
        fname = [fname,'_',num2str(Time(4)),'_',num2str(Time(5)),'_',num2str(Time(6))];
    end
end
mkdir([fname,'\data'])

% disp('Checking velocimeters')
% V(1) = system('tasklist /FI "IMAGENAME eq vectrino.exe" 2>NUL | find /I /N "vectrino.exe">NUL','-echo');
% V(2) = system('tasklist /FI "IMAGENAME eq vector.exe" 2>NUL | find /I /N "vector.exe">NUL','-echo');
% 
% if V(1)
% %     system('C:\Nortek\Vectrino\Vectrino.exe')
%     disp('Open Vectrino.exe. Start Vectrino in software.  Start data recording as well.')
% end
% if V(2)
% %     system('C:\Nortek\Vector\Vector.exe')
%     disp('Open Vector.exe. Start Vector in software.  Start data recording as well.')
% else
%     disp('ensure Vectrino and Vector are collecting data and recording to file')
% end

last_out = [0 0 0 0 0 0 0];
  write(dq,last_out)


% Find bias voltages for force and acceleration sensors
disp('All set. Turn on motor power. Press any key to run find_bias_3rigs');
[out,bias,dat] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil);

bias.pitch = [0 0 0]; 
% Find zero bias pitch angle by finding pitch with zero lift
disp('Run flume.  Click <a href="matlab: [last_out,bias] = find_zero_pitch(dq,last_out,bias,foil);">find_zero_pitch</a> when at full speed.')

