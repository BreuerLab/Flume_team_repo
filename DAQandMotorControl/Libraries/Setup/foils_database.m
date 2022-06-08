% Foils Database

% This function acts as a database for the different foil pairs that may be
% used in the Flume. It requires a string of the desired foil/foils as an input.

% Eric Handy, Dec 2022 - last edit

function [foil, rho, fs] = foils_database(foiltype)
%% Mostly non-changing variables

rho = 1000; % water density [kg/m3]
fs = 1000; % sampling frequency [Hz]

%% Foil data

FoilProperties =...
     [ " Chord [m]   :";" AspectRatio :";" Mass1 [kg]  :";" Mass2 [kg]  :";" Profile     :";" Material    :"]; % COMMENTS:
E1 = [     0.1;              3.5;          0.450;           0.450;         "Elliptical";    "Carbon fiber"  ]; % Built by someone
A1 = [     0.1;              3.5;          1.026;           1.174;         "Rectangular";   "Aluminum"      ]; % Nick's foils
A2 = [   3*0.0254;            6;           0.800;           0.800;         "Rectangular";   "Aluminum"      ]; % Yunxing's foils
C1 = [     0.054;             8.28;        0.386;           0;         "Cylindrical";   "Carbon fiber"  ]; % Joel's cylinder
V1 = [      0.06;              6.68;       0.306;           0;          "VibrissaeBeem50x"; "PLA3DprintWepoxy"]; % Joel's vibrissae model


foils = table(FoilProperties, E1, A1, A2, C1 ,V1); % constructing a table out of the foil data

selected_data = foils.(foiltype); % identifies the column of the selected foils

foil.chord = str2double(selected_data(1)); % asigns first column's value to the chord variable
foil.AR    = str2double(selected_data(2)); % asigns second column's value to the aspect ratio var
foil.mass1 = str2double(selected_data(3)); % asigns third column's value to the leading foil's mass var
foil.mass2 = str2double(selected_data(4)); % asigns fourth column's value to the trailing foil's mass var
foil.span  = foil.chord*foil.AR; % calculates the span
foil.profile  = selected_data(5); % asigns fifth column's value to the geometric profile of the foils
foil.material = selected_data(6); % asigns sixth column's value to the geometric profile of the foils

% foil % prints out the data of the selected foil

end