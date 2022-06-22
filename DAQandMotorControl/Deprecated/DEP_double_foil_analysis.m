%% Single Foil Analysis

% DEPRACATED

% Performs a full force and efficiency analysis for a single foil
% experiment, considering that the single foil's data is recorded using the
% Wallace rig.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT: FOR DATA TAKEN BEFORE 2022 05 27, MEASURED HEAVE [out(:,6)] MUST BE
% TAKEN NEGATIVE FOR CALCULATION. SUBSEQUENT DATA ACQUSITION HAS BEEN CORRECTED.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Eric Handy - May 2022

function [kin, res, par] = double_foil_analysis(out, Prof_out_angle, foiltype, rig, depth)

%%% INPUT
%
% out -------------> data variable
% Prof_out_angle --> pitching and heaving profiles in [deg] and [m]
% foiltype --------> from "foils_database.m"
% rig -------------> [2] for Gromit, [3] for Wallace
% depth -----------> (OPTIONAL) depth of flume during experiment [m]
%
%%% OUTPUT
%
% kin --> kinematics:
%
%           - heave -> heave position [m] (+ away from computer)
%           - pitch -> pitch position [rad] (+ counter clockwise)
%           - pdot --> pitch velocity [rad/s]
%           - hdot --> heave velocity [m/s]
%           - hacc --> heave acceleration [m/s2]
%
% res ---> resulting force calculations
%
%           - Tq ---> Z torque [Nm] (+ counter clockwise)
%           - Lift -> Lift force [N] (+ away from computer)
%           - Drag -> Drag force [N] (+ streamwise direction)
%           - PwrH -> Heave power [Nm/s]
%           - PwrP -> Pitch power [Nm/s]
%
%           - CM --> Moment coeff
%           - CL --> Lift coeff
%           - CD --> Drag coeff
%           - CPH -> Heave power coeff
%           - CPP -> Pitch power coeff
% 
%           - Eff -> Measured total efficiency
%           - Eff_prime -> Corrected efficiency
%
% par ---> relevant parameters
%
%           - U ---> Flow speed
%           - Re ---> Reynold's number
%           - St ---> Strouhal number
%           - Fr ---> Froude number
%           - beta -> Blockage ratio
%           - Yp ---> Maximum swept area by the leading edge of the foil (assuming a symmetrically mounted flat-plate foil)
%           - alphaT4 -> relative angle of attack at T/4 of the cycle (mid-stroke)

%           - P1 ----> leading pitching amplitude [deg]
%           - P2 ----> trailing pitching amplitude [deg]
%           - H1 ----> leading heaving amplitude [m]
%           - H2 ----> leading heaving amplitude [m]
%           - phi ---> phase between heaving and pitching [deg]
%           - phase13 -> phase between leading and trailing kinematics [deg]
%           - f_real --> experiment frequency [Hz]
%           - f_red ---> reduced frequency
%           - num_cyc -> number of cycles
%
%%%

% General static parameters

g = 9.80665; % acceleration of gravity
width = 0.8; % flume width [m]
nu = 1.00e-6; % water kinematic viscosity @22 degC (20220531)

if ~exist('depth','var')
    % parameter does not exist, default it to the following:
    depth = 0.55; % assumed flume water depth [m]
end

% Extract general experimental parameters

[U, P2, P3, H2, H3, f_real, f_red, phi, phase13, num_cyc, foil, rho, fs] = give_me_params(out, Prof_out_angle, foiltype);
U = U*1.026; % correction for the vectrino measurement (based on LDV and Vectrino measurements)

% Other general parameter calculations

Re = U*foil.chord/nu; % Reynolds number
Fr = U/sqrt(g*depth); % Froude number
St = f_real*foil.chord/U; % Strouhal number
% St = f_red*H2/(pi*foil.chord); % Another definition of St for pluging flight, considering the reduced frequency
%                                % and the heaving amplitude as the characteristic length.
alphaT4 = atan(-2*pi*(H2/foil.chord)*f_red) + deg2rad(P2); % relative angle of attack at T/4 (only relevant to the front foil)

% Extracting data from "out" variable according to the selected rig

if rig == 2 % for Gromit
    
    Fn = out(:,17); % normal force
    Ft = out(:,18); % tangential force
    Tq = -out(:,22); % z-torque (NEGATIVE BECAUSE OF THE COORDINATE SYSTEM TRANSFORMATION)
    
    heave = out(:,4); % heave position [m]
    pitch = out(:,3); % pitch angle [rad]
    
    P = P2; H = H2; mass = foil.mass1; % pitch [deg] and heave [m] amplitudes
    
elseif rig == 3 % for Wallace
    
    Fn = out(:,7); % normal force
    Ft = out(:,8); % tangential force
    Tq = -out(:,12); % z-torque (NEGATIVE BECAUSE OF THE COORDINATE SYSTEM TRANSFORMATION)
    
    heave = out(:,6); % heave position [m]
    pitch = out(:,5); % pitch angle [rad]
    
    P = P3; H = H3; mass = foil.mass2; % pitch [deg] and heave [m] amplitudes
    
else
    
    error("Rig must have value 2 (Gromit), or 3 (Wallace)");
    
end

% Fn = out(:,7); % normal force to the sensor [N]
% Ft = out(:,8); % tangential force to the sensor [N]
% Tq = out(:,12); % torque in the z direction (axis of rotation) [Nm] (negative respect to the decided positive pitching angle)

[b,a] = butter(6, f_real*10/(fs/2), 'low'); % filtering and smoothing the acceleration

% heave = out(:,6); % measured heave, this has to be flipped (pending correction) [m]
% pitch = -out(:,5); % measured pitch [rad] (positive counterclockwise) NOTE: pitch angle is opposite to forc e transducer's reference.

% Velocity and acceleration calculations

hdot = gradient(heave)*fs; % heave velocity
hdot = filtfilt(b, a, hdot); % filter heave velocity
pdot = gradient(pitch)*fs; % pitch velocity
pdot = filtfilt(b, a, pdot); % filter pitch velocity

hacc = gradient(hdot)*fs; % heave acceleration
% hacc = filtfilt(b, a, hacc);
clear a b; % clear unnecessary variables

% Lift and Drag calculations

% NOTE: Lift calculation accounts for the inertial load of the foil and
% the lower section of the force transducer

Lift = Fn.*cos(pitch) + Ft.*sin(pitch) - (mass+0.6)*hacc; % 600 gr from the lower part of the force sensor
Drag = - Fn.*sin(pitch) + Ft.*cos(pitch);

% Power calculation

PwrH = Lift.*hdot; % power due to heaving
PwrP = Tq.*pdot; % power due to pitching

% Lift, Drag, Torque and Power coefficients

CL = Lift/(0.5*rho*U^2*foil.span*foil.chord); % lift coeff
CD = Drag/(0.5*rho*U^2*foil.span*foil.chord); % drag coeff
CM = Tq/(0.5*rho*U^2*foil.span^2*foil.chord); % moment coeff
CPH = PwrH/(0.5*rho*U^3*foil.span*foil.chord); % heaving power coeff
CPP = PwrP/(0.5*rho*U^3*foil.span*foil.chord); % pitching power coeff

% Calculation of eficiency

[~, ~, PwrH_cyc] = phase_avg_data(pitch,PwrH); % averaged data per cycle respect to the pitching
[~, ~, PwrP_cyc] = phase_avg_data(pitch,PwrP);

% Yp = max(2*(Prof_out_angle(:,6) + 0.5*foil.chord*sin(deg2rad(Prof_out_angle(:,5))))); % maximum swept area

yp1 = heave + foil.chord*0.5*sin(pitch);
yp2 = heave - foil.chord*0.5*sin(pitch);
Yp = 2*max(max(yp1),max(yp2)); % this makes a big difference


EffH = mean(mean(PwrH_cyc))/(0.5*rho*U^3*Yp*foil.span); % single foil uncorrected heaving efficiency
EffP = mean(mean(PwrP_cyc))/(0.5*rho*U^3*Yp*foil.span); % single foil uncorrected pitching efficiency

Eff = EffH + EffP; % total efficiency

%% Blockage Correction (Houlsby et al, Ross et al, Ribeiro et al)

% Normalization of average CD

[~, ~, CD_cyc] = phase_avg_data(pitch,CD); % phase-averaged drag coeff
CD_avg = mean(mean(CD_cyc)); % mean drag coefficient per cycle
% H = H2; % heaving amplitude of foil (in meters)

CD_norm = CD_avg*(foil.chord/(2*H)); % normalized drag coeff to be used in calculations

beta = (2*H*foil.span)/(width*depth); % blockage ratio

% Iteration to find u1

u2 = U; % initial guess for u2
err = 1; % initializing the error
i = 0;

while err > 0.0001
    
    u1_a = (Fr^2*u2^4 - (4 + 2*Fr^2)*U^2*u2^2 + 8*U^3*u2 - 4*U^4 + 4*beta*CD_norm*U^4 + Fr^2*U^4)/...
        (-4*Fr^2*u2^3 + (4*Fr^2 + 8)*U^2*u2 - 8*U^3);
    
    u1_b = sqrt(u2^2 - CD_norm*U^2);
    
    err = abs(u1_a - u1_b);
    
    u2 = u2 + 0.00001;
    i = i + 1;
    
end

% Parameter correction

u1 = u1_a;

ud = (u1*(u2-U)*(2*g*depth - u2^2 - u2*U))/(2*beta*g*depth*(u2 - u1));

U_prime = U*((ud/U)^2 + CD_norm*0.25)/(ud/U); % corrected free-stream velocity

CD_norm_p = CD_norm*((U/U_prime)^2); % corrected normalized drag coefficient

Eff_prime = Eff*((U/U_prime)^3); % corrected efficiency

%% Recorded values:

kin.heave = heave; % heave position h(t)
kin.pitch = pitch; % pitch position p(t)
kin.hdot = hdot; % heave velocity
kin.pdot = pdot; % pitch velocity
kin.hacc = hacc; % heave acceleration

res.Torq = Tq; % Z torque
res.Lift = Lift; % Lift force
res.Drag = Drag; % Drag force
res.PwrH = PwrH; % Heave power
res.PwrP = PwrP; % Pitch power

res.CL = CL; % lift coeff
res.CD = CD; % drag coeff
res.CM = CM; % moment coeff
res.CPH = CPH; % heave power coeff
res.CPP = CPP; % pitch power coeff

res.Eff = Eff; % measured total efficiency
res.Eff_prime = Eff_prime; % corrected efficiency

par.U = U; % Flow speed
par.Re = Re; % Reynold's number
par.St = St; % Strouhal number
par.Fr = Fr; % Froude number
par.beta = beta; % blockage ratio
par.Yp = Yp; % maximum swept area by the leading edge of the foil (assuming a symmetrically mounted flat-plate foil)
par.alphaT4 = alphaT4; % relative angle of attack at T/4

par.P1 = P1; % leading pitching amplitude
par.P2 = P2; % trailing pitching amplitude
par.H1 = H1; % leading heaving amplitude
par.H2 = H2; % leading heaving amplitude
par.phi = phi; % phase between heaving and pitching
par.phase13 = phase13; % phase between leading and trailing kinematics
par.f_real = f_real; % experiment frequency
par.f_red = f_red; % reduced frequency
par.num_cyc = num_cyc; % number of cycles

end