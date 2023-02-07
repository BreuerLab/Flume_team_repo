%% Blockage Correction - Look at Ross' paper for more info on various methods

% copied from bernardo's paper
% This is the one in the force calculation code

% clear;clc;close all

function [beta, U_prime, u2, Eff_prime, CD_norm, CD_norm_p] = test_blockage_houlsby(pitch_prof, drag_coeff, heave_amp, Eff, U, Fr, foil, depth)

width = 0.8;

%% Using Houlsby correction - Blockage ratio equation is from Gauthier's paper / d0 = flume depth / Cd is also normalized according to Gauthier's paper

%In order to be a physical system: u2 > ut > u1, v0 > ut, u2 > v0

v0 = U; % Freestream velocity is 0.5m/s
error = 10;
% u2 is the estimated bypass velocity
u2 = U + 0.01; % 0.4m/s / 0.3m/s PIV - Guess (choose a value until you see convergence of Houlsby equations (Read Ross's paper for more info)
beta = (2*heave_amp*foil.span)/(width*depth); % Blockage: At/Ac = (2*h0*S)/(depth*width) = 2*0.1*0.35/(0.6*0.8) = 0.146 (h0 = 1c for all the kinematics we have)

CD = mean(drag_coeff);

Cd = CD * (foil.chord/(2*heave_amp));   %EXP: f12h1a55 - From experiments (not corrected yet) - Only extra info necessary to update efficiency value (normalization given by Gauthier's paper)
eta = Eff;                  %EXP: f12h1a55 - From experiments (not corrected yet)

% ONLY for highly loaded cases (corrects convergence)
if Cd > 1
    Cd = 1;
end
% Cd = 1; % overriding CD for all cases and not only high loading 20230109

g = 9.80665;
d0 = depth; %Putting same as channel depth (inspired from Ross's paper)
Fr = v0/(sqrt(g*d0));

while error > 0.0001
    u1_1 = ((Fr^2)*(u2^4) - (4+2*(Fr^2))*(v0^2)*(u2^2) + 8*(v0^3)*u2 - 4*(v0^4) + 4*beta*Cd*(v0^4) + (Fr^2)*(v0^4))/(-4*(Fr^2)*(u2^3) + (4*(Fr^2)+8)*(v0^2)*u2 - 8*(v0^3));
    u1_2 = sqrt((u2^2)-Cd*(v0^2));
    error = u1_1 - u1_2;
    u2 = u2 + 0.00000001;
end

ut = (u1_1*(u2-v0)*(2*g*d0 - (u2^2) - u2*v0))/(2*beta*g*d0*(u2-u1_1));

% Corrected freestream velocity:
primev0 = v0*((ut/v0)^2 + Cd/4)/(ut/v0);
%Updated turbine parameters with blockage correction:
primeCd = Cd*((v0/primev0)^2);
primeeta = eta*((v0/primev0)^3);

U_prime = primev0;
Eff_prime = primeeta;
CD_norm = Cd;
CD_norm_p = primeCd;

end

