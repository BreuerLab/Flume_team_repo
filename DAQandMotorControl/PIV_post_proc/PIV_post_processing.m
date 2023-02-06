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
tic;
for k = 1:N
    filename = ['B', num2str(k,'%04.0f'), '.dat'];
    [A,~,~] = importdata(filename); % load data from file
    
    raw_x = A.data(window_index,1);
    raw_y = A.data(window_index,2);
    raw_u(:,k) = A.data(window_index,3);
    raw_v(:,k) = A.data(window_index,4);
    raw_vort(:,k) = A.data(window_index,9);
    raw_corr(:,k) = A.data(window_index,12);
    raw_isVa(:,k) = A.data(window_index,15);
    raw_uncU(:,k) = A.data(window_index,13)/0.33; % percentage uncertainty
    raw_uncV(:,k) = A.data(window_index,14)/0.33; % percentage uncertainty
    
end

X = [raw_u; raw_v]; % compile velocity components into spatio-temporal matrix

toc;
%% Gappy Proper Orthogonal Decomposition (GPOD)
% Interpolates low correlation and missing vectors

gaps = gappy_mask(raw_vort, 'corr', raw_corr, 0.4);

Xgappy_u = raw_u.*gaps;
Xgappy_v = raw_v.*gaps;
Xgappy = [Xgappy_u; Xgappy_v]; % stack both velocity components on top of each other
Wgappy = raw_vort.*gaps;

[Xinterp, details_velo] = GPOD(Xgappy);
[Winterp, details_vort] = GPOD(Wgappy);

%% Robust Principle Component Analysis (RPCA)
% De-noises the data

[L,S] = RPCA(Xinterp); % performs RPCA on the velocity components

u_L = L(1:end/2,:);
v_L = L(end/2+1:end,:);

[vort_L,vort_S] = RPCA(Winterp); % performs RPCA on the vorticity components

%% Phase-averaging
% Phase averaging for periodic flows

[u_avg, v_avg, vort_avg] = phase_average_piv(u_L, v_L, vort_L, fpc, 92);
[uL_avg, vL_avg, vortL_avg] = phase_average_piv(u_L, v_L, raw_vort, fpc, 92);
[u_avg, v_avg, vort_avg] = phase_average_piv(raw_u, raw_v, raw_vort, fpc, 92);

save('RPCA_results.mat','X','L','S','vort_L','vort_S','raw_x','raw_y');
save('WRONG_Averaged_Flowfield.mat','u_avg','v_avg','vort_avg','raw_x','raw_y');
save('WRONG_Averaged_Flowfield_RPCA.mat','uL_avg','vL_avg','vortL_avg','raw_x','raw_y');

%% Testing (temporary)

x = reshape(raw_x,[n,m]);
y = reshape(raw_y,[n,m]);

% u = reshape(u_L,[n,m,N]);
% v = reshape(v_L,[n,m,N]);
% vort = reshape(vort_L,[n,m,N]);
% vort(abs(vort)<1) = 0;

u = reshape(uL_avg,[n,m,fpc]);
v = reshape(vL_avg,[n,m,fpc]);
% vort = reshape(vortL_avg,[n,m,fpc]);

% u = reshape(u_avg,[n,m,fpc]);
% v = reshape(v_avg,[n,m,fpc]);
% vort = reshape(vort_avg,[n,m,fpc]);

U = sqrt(u.^2 + v.^2);
% U = U - mean(U,3);

figure(3)
for frame = 1:fpc
    contourf(x,y,vort(:,:,frame),'linestyle','none','levelstep',0.005);
    colormap(brewermap([],'-RdBu'))
    caxis([-0.1,0.1])

%     quiver(x,y,u(:,:,frame),v(:,:,frame),2,'k');
    
%     contourf(x,y,U(:,:,frame),'linestyle','none','levelstep',0.005);
%     colormap(brewermap([],'-YlGnBu'))
    
    pause(0.1)
end

%% End
cd(main_folder);
