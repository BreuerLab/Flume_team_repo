function [positionAbs] = traversemove(daqObj, axisName, positionCurrent, positionTarget)
%[positionAbs] = traversemove(daqObj, axisName, positionCurrent, positionTarget)
%
% Belt traverse manual movement function.
% Commands a selected axis to move from current position to a target
% positon and hold.
% Returns the new position as position offset.
%
% v1.0.0, xh, 01/26/2023
%
% Motors must be enabled and homed.
%
% All position variables are in the absolute coordinates referenced to home
% positions, i.e., zeros are at the home positions.
%   - linear axis unit = [m], positive direction pointing toward
%                             the starboard/computer side of the flume.
%   - rotational axis unit = [deg], positive direction pointing downward
%                                   (clockwise when looking from above).
% -------------------------------------------------------------------------
% Inputs:
% daqObj = daq object created by >> daq('ni').
%
% axisName = char variable to specify the target axis
%            - 'y' = y/heave axis
%            - 'theta' = theta/pitch axis.
%
% positionCurrent = current absolut position on the selected axis
%                    referenced to its home zero.
%
% positionTarget = target absolut position of the selected axis referenced
%                   to its home zero.
% -------------------------------------------------------------------------
% Output:
% positionAbs = last absolute postion commanded by the daq, equivalent
%                to Last_out
% -------------------------------------------------------------------------
% Note:
% Change NchYcmd, NchYhold, NchTcmd, and NchThold accordingly if the output
% channels are changed in DAQ setup.
% -------------------------------------------------------------------------

%% column numbers of the output signals in the output matrix to daq
NchTotal = 11; % total number of columns of the output matrix
NchYcmd = 8; % y-axis command signal is output via 8th column of the output matrix
NchYlock = 9; % y-axis command lock signal is output via 9th column of the output matrix
NchTcmd = 10; % theta-axis command signal is output via 10th column of the output matrix
NchTlock = 11; % theta-axis command lock signal is output via 11th column of the output matrix

%% check inputs
if ~strncmp(class(daqObj), 'daq.interfaces.DataAcquisition', 30)
    error('errorType:daq', ...
          ['[positionOffset] = traversemove(daqObj, axisName, positionCurrent, positionTarget)' ...
           '\n\tdaqObj must be a DataAcquisition object,'...
           '\n\tdaqObj = daq(''ni'').'])
end

if ~isnumeric(positionCurrent) || ~isscalar(positionCurrent) || ~isnumeric(positionTarget) || ~isscalar(positionTarget)
    error('errorType:inputType', ...
          ['[positionOffset] = traversemove(daqObj, axisName, positionCurrent, positionTarget)' ...
           '\n\tpositionCurrent and positionTarget must be numeric scalar values.'])
end

if axisName == 'y' % y/heave axis
    calib = 10/0.5; % y: [m] to [V]
    K = 0.08; % maximum linear moving speed 0.08 [m/s]
    columnNum = [NchYcmd, NchYlock]; % locate column numbers in the output matrix to daq
    Ndigit = 4; % round to the 4th decimal digit for offset return (precision: 0.1 mm)
    if positionCurrent < 0 || positionCurrent > 0.5 || positionTarget < 0 || positionTarget > 0.5
        error('errorType:exceedLimit', ...
              ['[positionOffset] = traversemove(daqObj, axisName, positionCurrent, positionTarget)' ...
               '\n\tpositionCurrent and positionTarget must be within the y-axis limit,'...
               '\n\t0 <= position <= 0.5 [m] (absolut distance from home position).'])
    end
elseif axisName == 'theta' % theta/pitch axis
    calib = 10/360; % theta: [deg] to [V]
    K = 30; % maximum rotational moving speed 30 [deg/s]
    columnNum = [NchTcmd, NchTlock]; % locate column numbers in the output matrix to daq
    Ndigit = 2; % round to the 2nd decimal digit for offset return (precision: 0.01 deg)
    if positionCurrent < 0 || positionCurrent > 360 || positionTarget < 0 || positionTarget > 360
        error('errorType:exceedLimit', ...
              ['[positionOffset] = traversemove(daqObj, axisName, positionCurrent, positionTarget)' ...
               '\n\tpositionCurrent and positionTarget must be within the theta-axis limit,'...
               '\n\t0 <= position <= 360 [deg] (absolut rotation from home position).'])
    end
else
    error('errorType:axisName', ...
          ['[positionOffset] = traversemove(daqObj, axisName, positionCurrent, positionTarget)' ...
           '\n\taxisName must be either ''y'' or ''theta''.'])
end
%% data acquisition parameters
fs = daqObj.Rate; % sampling rate [samples/s]
Ts = 1/fs;
%% position profile
amp = positionTarget - positionCurrent;
as = 5; % acceleration smoothing parameter
t1 = 1;
t2 = t1 + abs(amp)/K;
t2 = t2 - (mod(t2, Ts) - Ts) - Ts*((mod(t2, Ts) - Ts/2) < 0); % round t2 to a multiple of Ts
t = (0 : Ts : t2+1)';
positionCmd = sign(amp)*(K/as*log(cosh(as*(t-t1))./cosh(as*(t-t2))))/2; % Eldredge's function
positionCmd = positionCmd -positionCmd(1) + positionCurrent; % position command in [m]

voltageCmd = positionCmd * calib; % voltage to motor [V]
cmdLock = zeros(length(voltageCmd), 1); % unlock "hold position"
cmdLock(end) = 1; % lock motor to target position at the last time step in the output

%% output to daq
positionAbs = round(positionCmd(end), Ndigit); % record the absolute position offset in [m]
                                                 % equivalent to Last_out

outScanData = zeros(length(positionCmd), NchTotal); % initiate output matrix

outScanData([NchYlock, NchTlock]) = 1; % set the motors to "command lock", motors hold current positions as defult
                                       % high-level input to motor I/O "Input A"

outScanData(:, columnNum) = [voltageCmd, cmdLock]; % store voltage command and command lock
                                                    % signals of the target axis to the output matrix

readwrite(daqObj, outScanData); % write to daq: voltage is sent to the motor now