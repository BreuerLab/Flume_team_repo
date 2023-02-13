%% PIV POST PROCESSING
% Updated 20230120

% This code does post processing on PIV data from DaVis, interpolating low
% correlation vectors using the Gappy POD algorithm and filtering the data
% set using Robust Principal Component Analysis RPCA in order to perform
% more extensive analysis such as modal analysis.

clear;

%% General parameters

freq_piv = 14.9316; % frequency of the piv snapshots
fpc = 23; % frames per cycle

main_folder = pwd;
% piv_folder = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export');
piv_folder = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\Correction_01\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export');

addpath(main_folder);

%% Load data

cd(piv_folder);

files = dir('*.dat'); % finds all filenames ending with .dat extension

%% Calibration

load('piv_calibration.mat');

%% Foil coordinates

%% Frame cropping

% Load the first file
[A,~,~] = importdata(files(1).name); 

% Get the flowfield size from the text data in the piv files
newStr = extract(string(A.textdata(3)),digitsPattern);
I = str2double(newStr(2));
J = str2double(newStr(3));

% load x, y, u, v, I, J into piv_frame_cropping
[range_x, range_y, window_index] = piv_frame_cropping(...
    reshape(A.data(:,1),[I,J]), reshape(A.data(:,2),[I,J]), reshape(A.data(:,3),[I,J]), reshape(A.data(:,4),[I,J]), I, J);

n = length(range_x);
m = length(range_y);

%% PIV data compiling
fprintf('Processing...\n\n');

N = 100;%length(files);

% initialize variables
x_raw = nan(n*m,1);
y_raw = nan(n*m,1);
u_raw = nan(n*m,N);
v_raw = x_raw;
vort_raw = x_raw;
corr_raw = x_raw;
isVa_raw = x_raw;
uncU_raw = x_raw;
uncV_raw = x_raw;

tic;
for k = 1:N
    filename = ['B', num2str(k,'%04.0f'), '.dat'];
    [A,~,~] = importdata(filename); % load data from file
    
    x_raw = A.data(window_index,1);
    y_raw = A.data(window_index,2);
    u_raw(:,k) = A.data(window_index,3);
    v_raw(:,k) = A.data(window_index,4);
    vort_raw(:,k) = A.data(window_index,9);
    corr_raw(:,k) = A.data(window_index,12);
    isVa_raw(:,k) = A.data(window_index,15);
    uncU_raw(:,k) = A.data(window_index,13)/0.33; % percentage uncertainty
    uncV_raw(:,k) = A.data(window_index,14)/0.33; % percentage uncertainty
    
end
toc;

X = [u_raw; v_raw]; % compile velocity components into spatio-temporal matrix

%% Gappy Proper Orthogonal Decomposition (GPOD)
% Interpolates low correlation and missing vectors

skip = 0; % change to enable GPOD interpolation

if skip == 1
    gaps = gappy_mask(vort_raw, 'corr', corr_raw, 0.4);

    Xgappy_u = u_raw.*gaps;
    Xgappy_v = v_raw.*gaps;
    Xgappy = [Xgappy_u; Xgappy_v]; % stack both velocity components on top of each other
    Wgappy = vort_raw.*gaps;

    [Xinterp, details_velo] = GPOD(Xgappy);
    [Winterp, details_vort] = GPOD(Wgappy);
end

%% Robust Principle Component Analysis (RPCA)
% De-noises the data

skip = 0; % change to enable RPCA filtering

if skip == 1
    [L,S] = RPCA(Xinterp); % performs RPCA on the velocity components

    u_L = L(1:end/2,:);
    v_L = L(end/2+1:end,:);

    [vort_L,vort_S] = RPCA(Winterp); % performs RPCA on the vorticity components
end

%% Phase-averaging
% Phase averaging for periodic flows

[u_avg, v_avg, vort_avg] = phase_average_piv(u_raw, v_raw, vort_raw, fpc, 92);

out = track_vortex(x_raw, y_raw, u_avg, v_avg, [m,n]);

save('Averaged_Flowfield.mat','u_avg','v_avg','vort_avg','x_raw','y_raw');

%% End
cd(main_folder);
