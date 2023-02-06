function varargout = traversecmd(axisName, positionCmd, positionOffset)
%[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)
%[voltageCmd, cmdLock] = belt_traverse(axisName, positionCmd, positionOffset)
%
% Belt traverse analog control signal convertor.
% Converts position input [m] (linear axis) or [deg] (rotational axis) to
% voltage output [V] for daq output matrix for motor control.
%
% v1.0.1, added cmdLock output, xh, 01/31/2023.
% v1.0.0, xh, 01/26/2023.
% -------------------------------------------------------------------------
% Inputs:
% axisName = char variable to specify the target axis
%            - 'y' = y/heave axis
%            - 'theta' = theta/pitch axis.
%
% positionCmd = position command referenced to the relative zero
%               - y-axis in [m], positive direction pointing toward
%                 the starboard/computer side of the flume.
%               - theta-axis in [deg], positive direction pointing downward
%                 (clockwise when looking from above).
%
% positionOffset = absolute position offset referenced to the home position
%                  (absolute zero)
%                  - y-axis in [m], setting up zero heave,
%                    e.g., positionOffset = 0.25[m] at flume center.
%                  - theta-axis in [deg], setting up zero pitch,
%                    e.g., positionOffset = 180[deg] when the x-axis
%                    pointing toward the downstream direction.
% -------------------------------------------------------------------------                          
% Outputs:
% voltageCmd = voltage signal to control motor position (0 - 10 [V])
%              to analog output channels.
%
% cmdLock = command lock signal (0 = unlock, 1 = lock)
%           to digital output channels.
% -------------------------------------------------------------------------

%% arguments check
nargoutchk(1, 2)
narginchk(3, 3)

%% check inputs
if axisName == 'y' % y/heave axis
    calib = 10/0.5; % y: [m] to [V]
    if positionOffset < 0 || positionOffset > 0.5
        error('errorType:exceedLimit', ...
              ['[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)' ...
               '\n\tpositionOffset must be within the y-axis limit,'...
               '\n\t0 <= positionOffset <= 0.5 [m] (absolut distance from home position).'])
    end
    if min(positionCmd)+positionOffset < 0 || max(positionCmd)+positionOffset > 0.5
        error('errorType:exceedLimit', ...
              ['[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)' ...
               '\n\tpositon_cmd + positionOffset must be within the y-axis limit,'...
               '\n\t0 <= positionCmd + positionOffset <= 0.5 [m].'])
    end
elseif axisName == 'theta' % theta/pitch axis
    calib = 10/360; % theta: [deg] to [V]
    if positionOffset < 0 || positionOffset > 360
        error('errorType:exceedLimit', ...
              ['[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)' ...
               '\n\tpositionOffset must be within the theta-axis limit,'...
               '\n\t0 <= positionOffset <= 360 [deg] (absolute rotation from home position).'])
    end
    if min(positionCmd)+positionOffset < 0 || max(positionCmd)+positionOffset > 360
        error('errorType:exceedLimit', ...
              ['[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)' ...
               '\n\tpositon_cmd + positionOffset must be within the theta-axis limit,'...
               '\n\t0 <= positionCmd + positionOffset <= 360 [deg].'])
    end
else
    error('errorType:axisName', ...
          ['[voltageCmd] = belt_traverse(axisName, positionCmd, positionOffset)' ...
           '\n\taxisName must be either ''y'' or ''theta''.'])
end
%% signal conversion
voltageCmd = (positionCmd + positionOffset) * calib; % voltage to motor [V]
cmdLock = zeros(length(voltageCmd), 1); % unlock "hold position"
cmdLock(end) = 1; % lock motor at the last time step in the output

%% output
if nargout == 1
    % generate voltage command signal only
    varargout{1} = voltageCmd;
else
    % generate voltage command and command lock signals
    varargout{1} = voltageCmd;
    varargout{2} = cmdLock;
end

end