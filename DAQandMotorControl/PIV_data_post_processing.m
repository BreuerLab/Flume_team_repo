%% PIV Data Post-Processing
% Courtesy of Yuanhang 20220815
% Edited by yours truly, El Eric

% This is the OG code
% does NOT implement RPCA and GPOD

clear;

%% General parameters

frame = 23; % number of piv frames per cycle
piv_freq = 14.9316; % frequency of PIV frames in Hz

%% Use Gappy POD?

use_GPOD = 0; % 1/0 = yes/no

%% PIV folder

piv_folder = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=68_p3=75_h3=0.8_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU\Export');

%% Load Force data

force_folder = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\');
force_filename = ('20221006_alpha=68_p3=75_h3=0.8_ph=-120_A3E.mat'); % force data file

load(fullfile(force_folder,force_filename));

%% Force calculations

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);

[toverT1,          ~, CL2_cyc] = cycle_avg_data(kin.p2_comm, res.CL2); % LiftC2
[toverT2, pitch_cyc2, CP2_cyc] = cycle_avg_data(kin.p2_comm, (res.CPH2+res.CPP2)); % PowerC2

[toverT3,          ~, CL3_cyc] = cycle_avg_data(kin.p3_comm, res.CL3, samplerate, 0); % LiftC3
[toverT4, pitch_cyc3, CP3_cyc] = cycle_avg_data(kin.p3_comm, (res.CPH3+res.CPP3), samplerate, 0); % PowerC3

chord = foil.chord; % foil chord [m]

foil_separation = 6*chord; % [m] (should be given in the force data file)

%% Calibration: Coordinate translation and rotation (only necessary for ONE dataset in ONE experiment in the whole PIV session)
% Assume that all data taken in the PIV session have the SAME slight misalighnment

cd(piv_folder); % load from the PIV folder

coords_le = NaN(frame,2); % initialize variables
coords_tr = NaN(frame,2);

for ii = 1:frame
    
    filename = ['B' num2str(ii,'%04.0f') '.dat'];
    [A,~,~] = importdata(filename); % load data from file
    
    data_no_x = 551; % from "text data" in the file structure (I values)
    data_no_y = 303; % from "text data" in the file structure (J values)

    raw_x = A.data(:,1);            data_x = reshape(raw_x,[data_no_x,data_no_y]); % X dimension
    raw_y = A.data(:,2);            data_y = reshape(raw_y,[data_no_x,data_no_y]); % Y dimension
    raw_isValid = A.data(:,15);     data_isValid = reshape(raw_isValid,[data_no_x,data_no_y]); % isValid
    
    figure('Name', 'FRONT and HIND Foil Coordinates', 'WindowState', 'maximized');
    contour(data_x, data_y, data_isValid);
    colormap('bone');
    FigTit = ['Frame: ',num2str(ii)];
    title(FigTit);
    
    fprintf(['\n\nSelect approximate foil center coordinates for the LEADING and TRAILING foil in that order.\n\nCurrent frame: ',num2str(ii),'\n']);
    leading = drawpoint('Label','Leading Foil');
    trailing = drawpoint('Label','Trailing Foil');
    coords_le(ii,:) = leading.Position;
    coords_tr(ii,:) = trailing.Position;
    
    close 'FRONT and HIND Foil Coordinates';
    
end

center_le(1) = mean(coords_le(:,1)); % get the centerline of the leading foil (x)
center_le(2) = (max(coords_le(:,2)) + min(coords_le(:,2)))/2; % get the centerline of the leading foil (y)

center_tr(1) = mean(coords_tr(:,1)); % get the centerline of the trailing foil (x)
center_tr(2) = (max(coords_tr(:,2)) + min(coords_tr(:,2)))/2; % get the centerline of the trailing foil (y)

center_axis(1) = (center_le(:,1) + center_tr(:,1))/2; % get the origin of all axes (x)
center_axis(2) = (center_le(:,2) + center_tr(:,2))/2; % get the origin of all axes (y)

center_le_d = center_le - center_axis; % correct the leading centerline
center_tr_d = center_tr - center_axis; % correct the trailing centerline

disp_angle = atan(center_tr_d(2)/center_tr_d(1)); % get the displacement angle

R = [cos(disp_angle), -sin(disp_angle); sin(disp_angle), cos(disp_angle)]; % rotation matrix

center_le_p = center_le_d*R; % rotates the translated center coordinates
center_tr_p = center_tr_d*R; % rotates the translated center coordinates

figure()
plot([center_le(1),center_tr(1)],[center_le(2),center_tr(2)],'--r'); hold on;
plot([center_le_p(1),center_tr_p(1)],[center_le_p(2),center_tr_p(2)],'--k'); hold off;

save('piv_calibration.mat','center_axis','R','disp_angle'); % save calibration parameters to a file

%% If previous calibration is being used:

cd(piv_folder); % load from the PIV folder

% If calibration was previously done:
load('piv_calibration.mat');

data_no_x = 551; % from "text data" in the file structure (I values)
data_no_y = 303; % from "text data" in the file structure (J values)

%% Foil coordinates
% NOTE: ideally instead of using the isValid array to determine the foil's position, a raw image should be
% used along with the calibration determined in the previous section.

% flip measurement frame of reference from the encoders to the PIV
p2 = -out(:,3); % pitch leading foil [rad] (THESE SHOULD ALL BE NEGATIVE)
h2 = -out(:,4); % heave leading foil [m]
p3 = -out(:,5); % % pitch trailing foil [rad]
h3 = -out(:,6); % heave trailing foil [m]

% Choose the leading foil's coordinates
[A,delimiterout,headerlinesout] = importdata('B0001.dat'); % load data from file

raw_x = A.data(:,1);            data_x = reshape(raw_x,[data_no_x,data_no_y]); % X dimension w/ calibration
raw_y = A.data(:,2);            data_y = reshape(raw_y,[data_no_x,data_no_y]); % Y dimension
% raw_isValid = A.data(:,15);     isValid = reshape(raw_isValid,[data_no_x,data_no_y]); % isValid
raw_vort = A.data(:,15);        vort = reshape(raw_vort,[data_no_x,data_no_y]); % vorticity

% Apply calibration
data_xt = data_x - center_axis(1); % translate
data_yt = data_y - center_axis(2); % translate
data_x = (data_xt.*cos(disp_angle) - data_yt.*sin(disp_angle))/1000; % rotate and convert to encoder data units [m]
data_y = (data_xt.*sin(disp_angle) + data_yt.*cos(disp_angle))/1000; % rotate and convert to encoder data units [m]

figure('Name', 'Leading and Trailing EDGE coordinates', 'WindowState', 'maximized');
% contour(data_x,data_y,isValid); hold on;
contourf(data_x,data_y,vort); hold on;
colormap('jet');
FigTit = 'Front Foil Leading and Trailing Edge';
title(FigTit);

fprintf('\nSelect coordinates of the Leading Edge and Trailing Edge of the FRONT FOIL, in that order.\n\n');
% pick leading and trailing edge of the Leading foil
l_e = drawpoint('Label','Leading Edge'); l_e = l_e.Position; % pick the leading edge and only use the extracted position
t_e = drawpoint('Label','Trailing Edge'); t_e = t_e.Position; % pick the leading edge and only use the extracted position

plot([l_e(1),t_e(1)],[l_e(2),t_e(2)],'k','LineWidth',8);
pause(1);
close 'Leading and Trailing EDGE coordinates';

% Save coordinates into front foil variable

foil2_coords(1,:) = [l_e(1),  l_e(2),  NaN,  NaN,  t_e(1),  t_e(2)];
                  %  le_x,    le_y,    mc_x, mc_y, te_x,    te_y

foil2_coords(1,3) = foil2_coords(1,5) + (foil2_coords(1,1) - foil2_coords(1,5))/2; % x mid-chord [m]
foil2_coords(1,4) = foil2_coords(1,6) + (foil2_coords(1,2) - foil2_coords(1,6))/2; % y mid-chord [m]

% [BEST METHOD] using the trigger signal to match the position of the foil with the force data

% [~, I0] = min(1-out(:,24)); % trigger signal

% [ALTERNATE METHOD] calculating the leading foil position in the force data

foil2_angle0 = atan((foil2_coords(1,2) - foil2_coords(1,6))/(foil2_coords(1,1) - foil2_coords(1,5))); % starting angle of the leading foil based on the coordinates for the initial position of the leading and trailing edges
foil2_angle0 = atan((foil2_coords(1,6) - foil2_coords(1,2))/(foil2_coords(1,5) - foil2_coords(1,1)));
p2_temp = p2;
transient_time = round((1/freq)*transientcycs*samplerate); % time taken up by the transients
p2_temp(1:transient_time) = NaN; % make transient data irrelevant
p2_temp(end-transient_time:end) = NaN; % make transient data irrelevant

p2_temp(gradient(p2_temp) > 0) = NaN; % [CHANGE CONDITIONAL TO SEARCH FOR ANGLE DURING THE DESIRED STROKE]
[~, I0] = min(abs(p2_temp - foil2_angle0));  % time instance at which we match the beginning of the force and piv data (NOTE: not the actual moment at which they are matched, but correct within the cycle)

% [TESTING METHOD] adding a time delay to the trigger signal

[~, I0] = min(1-out(:,24)); % trigger signal
I0 = I0 + 4700;

figure(); hold on; % time alignment [FOR DEBUGGING]
plot(rad2deg(out(:,3)),'b'); grid on;
plot(rad2deg(out(:,5)),'r');
plot(out(:,24)*max(rad2deg(out(:,3))),'g','linewidth',2);
xline(I0,'r','LineWidth',2); hold off;
xlim([(I0-2000),(I0+2000)])

% Leading foil coordinates (based only on the encoder measurement)

f2_xpos = -foil_separation/2;% - 0.15*chord; % adding  SLIGHT displacement correction

foil2_le = [f2_xpos - cos(p2(I0))*chord/2,  h2(I0) - sin(p2(I0))*chord/2];
foil2_mc = [f2_xpos,                        h2(I0)];
foil2_te = [f2_xpos + cos(p2(I0))*chord/2,	h2(I0) + sin(p2(I0))*chord/2];

foil2_coords(1,:) = [foil2_le, foil2_mc, foil2_te];
    
% Trailing foil coordinates

f3_xpos = foil_separation/2;% + 0.1*chord; % adding  SLIGHT displacement correction

foil3_le = [f3_xpos - cos(p3(I0))*chord/2,    h3(I0) - sin(p3(I0))*chord/2]; % trailing foil leading edge
foil3_mc = [f3_xpos,                          h3(I0)]; % trailing foil mid-chord
foil3_te = [f3_xpos + cos(p3(I0))*chord/2,    h3(I0) + sin(p3(I0))*chord/2]; % trailing foil trailing edge

foil3_coords(1,:) = [foil3_le, foil3_mc, foil3_te];

piv_period = 1/piv_freq; % time between each piv frame
time_steps = piv_period/(1/samplerate); % time steps between piv frames

% figure(1); % [FOR DEBUGGING]
% plot(foil2_coords(1,[1,3,5]), foil2_coords(1,[2,4,6]), 'k', 'LineWidth', 8); hold on; % plot leading foil
% plot(foil3_coords(1,[1,3,5]), foil3_coords(1,[2,4,6]), 'k', 'LineWidth', 8); hold off; % plot trailing foil
% axis equal;
% xlim([-0.250,0.250]);
% ylim([-0.175,0.100]);

tSteps(1,1) = 0;%I0; % saving the time steps

for jj = 2:frame+1
    
    time_index = I0 + round((jj-1)*time_steps);
    
    foil2_le = [f2_xpos - cos(p2(time_index))*chord/2,  h2(time_index) - sin(p2(time_index))*chord/2];
    foil2_mc = [f2_xpos,                                h2(time_index)];
    foil2_te = [f2_xpos + cos(p2(time_index))*chord/2,	h2(time_index) + sin(p2(time_index))*chord/2];
    
    foil2_coords(jj,:) = [foil2_le, foil2_mc, foil2_te];
    
    foil3_le = [f3_xpos - cos(p3(time_index))*chord/2,  h3(time_index) - sin(p3(time_index))*chord/2];
    foil3_mc = [f3_xpos,                                h3(time_index)];
    foil3_te = [f3_xpos + cos(p3(time_index))*chord/2,  h3(time_index) + sin(p3(time_index))*chord/2];
    
    foil3_coords(jj,:) = [foil3_le, foil3_mc, foil3_te];
    
%     figure(2); % [FOR DEBUGGING]
%     plot(foil2_coords(jj,[1,3,5]), foil2_coords(jj,[2,4,6]), 'k', 'LineWidth', 8); hold on; % plot leading foil
%     plot(foil3_coords(jj,[1,3,5]), foil3_coords(jj,[2,4,6]), 'k', 'LineWidth', 8); hold off; % plot trailing foil
%     axis equal;
%     xlim([-0.25,0.25]);
%     ylim([-0.15,0.15]);
%     pause(0.5)
    
    tSteps(jj,1) = time_index - I0;
    
end

% Non-dimensionalize foil coordinates
foil2_coords = foil2_coords/chord;
foil3_coords = foil3_coords/chord;

%% PIV data

fov_x0 = -3.9; % dimensions of desired window to be plotted [chords]
fov_x1 = 4.5;  % obtained from results on Davis
fov_y0 = -2;
fov_y1 = 2;

for ii = 1:frame
    
    cd(piv_folder);
    
    file = dir('*.dat'); % finds all filenames ending with .dat extension

    N = length(file); % total number of files

    n = ii:frame:N;

    for jj = 1:1:length(n)

        filename = ['B' num2str(n(jj),'%04.0f') '.dat'];
        [A,delimiterout,headerlinesout] = importdata(filename); % load data from file
        % MAKE SURE THE COLUMN NUMBER FOR EACH VARIABLE IS CORRECT
        raw_x = A.data(:,1); % import x coords, non-dim with the foil chord
        raw_y = A.data(:,2); % import y coords, non-dim with the foil chord
        raw_u(:,jj) = A.data(:,3)/U; % import u velocity, non-dim with freestream
        raw_v(:,jj) = A.data(:,4)/U; % import v velocity, non-dim with freestream [negative due to perspective of camera]
        raw_vort(:,jj) = A.data(:,9)*(chord/U); % import w vorticity, non-dim with chord and freestream [negative due to perspective of camera]
        raw_correlation(:,jj) = A.data(:,12); % import correlation values
        raw_isValid(:,jj) = A.data(:,15); % import logical vector matrix
        
        % [START OF DEBUGGING]
%         data_x = reshape(raw_x,[data_no_x,data_no_y]);
%         data_y = reshape(raw_y,[data_no_x,data_no_y]);
%         % Apply calibration
%         data_xt = data_x - center_axis(1); % translate
%         data_yt = data_y - center_axis(2); % translate
%         data_x = (data_xt.*cos(disp_angle) - data_yt.*sin(disp_angle))/1000/chord; % rotate and non-dinemsionalize
%         data_y = (data_xt.*sin(disp_angle) + data_yt.*cos(disp_angle))/1000/chord; % rotate and non-dinemsionalize
%         data_vort = reshape(squeeze(raw_vort(:,jj)),[data_no_x,data_no_y]); data_vort(abs(data_vort) < 3.5) = 0;
%         data_isValid = reshape(squeeze(raw_isValid(:,jj)),[data_no_x,data_no_y]);
%         data_vort_temp = data_vort;
%         data_vort_temp(data_isValid == 0) = NaN; % to account for the mask for the foil shadow when plotting
%         contourf(data_x, data_y, data_vort_temp, 127, 'LineStyle', 'none'); hold on;
% %       contourf(data_x, data_y, data_isValid, 'LevelStep', 0.01, 'LineStyle', 'none'); for the valid and invalid vectors
% %       colormap(bone)
%         plot(foil2_coords(ii,[1,3,5])-0.07, foil2_coords(ii,[2,4,6])-0.05, 'k', 'LineWidth', 11); % plot leading foil
%         plot(foil3_coords(ii,[1,3,5])+0.04, foil3_coords(ii,[2,4,6])-0.07, 'k', 'LineWidth', 11); % plot trailing foil
%         colorbarpwn(-30, 30,'level', 127, 'log', 0.5, 'wrs', 0.01, 'colorP', [0.5, 0 , 0], 'colorN', [0, 0, 0.5]);
%         a = gca;
%         set(a,'color',[0.9,0.9,0.9]);
        % [END OF DEBUGGING]
        
    end

    % Cycle-averaging data
    raw_u = mean(raw_u,2);
    raw_v = mean(raw_v,2);
    raw_vort = mean(raw_vort,2);
    raw_correlation = mean(raw_correlation,2);
    raw_isValid = mean(raw_isValid,2);

    % Reshape x and y
    data_x = reshape(raw_x,[data_no_x,data_no_y]);
    data_y = reshape(raw_y,[data_no_x,data_no_y]);

    % Apply calibration
    data_xt = data_x - center_axis(1); % translate
    data_yt = data_y - center_axis(2); % translate
    data_x = (data_xt.*cos(disp_angle) - data_yt.*sin(disp_angle))/1000/chord; % rotate
    data_y = (data_xt.*sin(disp_angle) + data_yt.*cos(disp_angle))/1000/chord; % rotate
    
    % Reshape data
    data_u = reshape(raw_u,[data_no_x,data_no_y]);
    data_v = reshape(raw_v,[data_no_x,data_no_y]);
    data_vort = reshape(raw_vort,[data_no_x,data_no_y]);
    data_correlation = reshape(raw_correlation,[data_no_x,data_no_y]);
    data_isValid = reshape(raw_isValid,[data_no_x,data_no_y]);

    % % Cropping a specific window from the data
    
%     figure(); % [VISUALIZATION]
%     quiver(data_x, data_y, data_u, data_v, 'r'); hold on;
    
    [~, Ix0] = min(abs(data_x(:,1) - fov_x0)); % finds the index of the approximate desired location within the x and y matrices
    [~, Ix1] = min(abs(data_x(:,1) - fov_x1));
    [~, Iy0] = min(abs(data_y(1,:) - fov_y0));
    [~, Iy1] = min(abs(data_y(1,:) - fov_y1));

    range_x = Ix0:Ix1;
    range_y = Iy1:Iy0; % flipped because in the matrices, y0 is at a larger index than y1
    
    % Save the cropped-in data into the arrays
    data_x = data_x(range_x,range_y);
    data_y = data_y(range_x,range_y);
    data_u = data_u(range_x,range_y);
    data_v = data_v(range_x,range_y);
    data_vort = data_vort(range_x,range_y); data_vort(abs(data_vort) < 0.6) = 0; % clean-up for plotting
    data_correlation = data_correlation(range_x,range_y);
    data_isValid = data_isValid(range_x,range_y);
    
%     quiver(data_x, data_y, data_u, data_v, 'b'); % [VISUALIZATION]
%     hold off;
    
    %% saving data into compiled data
    
    foil2_coords(ii,[1,3,5]) = foil2_coords(ii,[1,3,5])-0.07; % x mild corrections
    foil2_coords(ii,[2,4,6]) = foil2_coords(ii,[2,4,6])-0.05; % y
    
    foil3_coords(ii,[1,3,5]) = foil3_coords(ii,[1,3,5])+0.04; % x
    foil3_coords(ii,[2,4,6]) = foil3_coords(ii,[2,4,6])-0.04; % y
    
    compiled_data(ii).timeStep = tSteps(ii);
    compiled_data(ii).foil2_coords = foil2_coords(ii,:);
    compiled_data(ii).foil3_coords = foil3_coords(ii,:);
    compiled_data(ii).x = data_x;
    compiled_data(ii).y = data_y;
    compiled_data(ii).u = data_u;
    compiled_data(ii).v = data_v;
    compiled_data(ii).vort = data_vort;
    compiled_data(ii).corrl = data_correlation;
    compiled_data(ii).isValid = data_isValid;
    
    %% plotting
    
    time_index = ii*(1/frame); % for the force plots
    
    h = figure;
%     set(h,'unit','centimeters','position',[3,3,45,23]); % for only the vorticity plot
    set(h,'unit','centimeters','position',[2,2,45,23]); % for vorticity and forces
%     set(h,'unit','centimeters','position',[2,2,27,23]); % for vorticity and forces
    set(h,'color','w');
    
    data_vort_temp = data_vort;
    data_vort_temp(data_isValid == 0) = NaN; % to account for the mask for the foil shadow when plotting
    
    % % Vorticity plot
%     subplot(2, 2, [1,2]);
    
    contourf(data_x, data_y, data_vort_temp, 91, 'LineStyle', 'none'); hold on; %, 127
%     colorbarpwn(-30, 30,'level', 127, 'log', 0.5, 'wrs', 0.02, 'colorP', [0.6, 0 , 0], 'colorN', [0, 0, 0.6], 'label', '$w_z (1/s)$'); % for low-med alpha
%     colorbarpwn(-40,40,'log', 'wrs', 0.066, 'colorP', [0.6, 0.05, 0.0], 'colorN', [0.0, 0.37, 0.62], 'label', '$w_z$  (1/s)'); % low aT4 (both ph)
%     colorbarpwn(-30,30,'log', 'wrs', 0.064, 'colorP', [0.6, 0.05, 0.0], 'colorN', [0.0, 0.37, 0.62], 'label', '$w_z$  (1/s)'); % medium aT4 (ph = 60)
%     colorbarpwn(-30,30,'log', 'wrs', 0.038, 'colorP', [0.6, 0.05, 0.0], 'colorN', [0.0, 0.37, 0.62], 'label', '$w_z$  (1/s)'); % medium aT4 (ph = -120)
    colorbarpwn(-30, 30, 'log', 'wrs', 0.045, 'colorP', [0.6, 0.05, 0.0], 'colorN', [0.0, 0.37, 0.62], 'label', '$w_z$  (1/s)'); % for high alpha
    a = gca;
    set(a,'color',[0.8,0.8,0.8]);

    % Foil masks
    plot(foil2_coords(ii,[1,3,5]), foil2_coords(ii,[2,4,6]), 'k', 'LineWidth', 10); % plot leading foil
    plot(foil3_coords(ii,[1,3,5]), foil3_coords(ii,[2,4,6]), 'k', 'LineWidth', 10); % plot trailing foil
    
    hold off;
    axis equal
    ylim([-2,2]); % symmetric plot window
    
    xlabel('$x/c$','fontsize',22,'interpreter','Latex');
    ylabel('$y/c$','fontsize',22,'interpreter','Latex');
    set(gca,'fontname','Times New Roman','fontsize',22,'linewidth',1.5,'TickDir','in');
    
    % % Force and power plots
%     subplot(2, 2, 3);
%     yyaxis left
%     plot(toverT3(1:end/2), mean(CL3_cyc(:,1:end/2)), 'LineWidth', 2); hold on;
%     xline(time_index,'r','linewidth',2); hold off;
%     ylabel('$C_{L,tr}$','fontsize',20,'interpreter','Latex');
%     yyaxis right
%     plot(toverT3(1:end/2), mean(pitch_cyc3(:,1:end/2)), 'LineWidth', 2);
%     xlabel('$t/T$','fontsize',20,'interpreter','Latex');
%     
%     xlim([0,0.5])
%     set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');
%     grid on; box on;
%     
%     subplot(2, 2, 4);
%     yyaxis left
%     plot(toverT3(1:end/2), mean(CP3_cyc(:,1:end/2)), 'LineWidth', 2); hold on;
%     xline(time_index,'r','linewidth',2); hold off;
%     ylabel('$C_{P,tr}$','fontsize',20,'interpreter','Latex');
%     yyaxis right
%     plot(toverT3(1:end/2), mean(pitch_cyc3(:,1:end/2)), 'LineWidth', 2);
%     xlabel('$t/T$','fontsize',20,'interpreter','Latex');
%     ylabel('$\theta_{tr}$','fontsize',20,'interpreter','Latex');
%     
%     xlim([0,0.5])
%     set(gca,'fontname','Times New Roman','fontsize',20,'linewidth',1.5,'TickDir','in');
%     grid on; box on;


    % % Saving the figure
    
    figname = ['fig', num2str(ii), '.fig'];
%     saveas(gcf,figname) % alternate way of saving
    export_fig(figname); % saves the matlab figure file
    
    figname = ['fig', num2str(ii), '.tif'];
%     saveas(gcf,figname) % alternate way of saving
    export_fig(figname); % saves the .tif file
    
end

%% Gappy POD

if use_GPOD == 1
    
    % FALTA MUCHO AQUI

    % Reshape data to tall matrices
    
    parfor n = 1:length(compiled_data)
        x(:,:,n) = compiled_data(n).x;
        y(:,:,n) = compiled_data(n).y;
        u(:,:,n) = compiled_data(n).u;
        v(:,:,n) = compiled_data(n).v;
        corrl(:,:,n) = compiled_data(n).corrl;
        isValid(:,:,n) = compiled_data(n).isValid;
    end
    
    for ii = 1:size(ug,3)
        ug_temp = reshape(ug(:,:,ii),1,[]); ug_tall(:,ii) = ug_temp; % reshape each snapshot into a column vector and store it as a column of the complete x dataset
        vg_temp = reshape(vg(:,:,ii),1,[]); vg_tall(:,ii) = vg_temp; % reshape each snapshot into a column vector and store it as a column of the complete y dataset

        uf_temp = reshape(u(:,:,ii),1,[]); uf_tall(:,ii) = uf_temp; % for the fully resolved flow field
        vf_temp = reshape(v(:,:,ii),1,[]); vf_tall(:,ii) = vf_temp;
    end
    
    % Define gappy mask and foil mask
    
    foil_mask = ones(size(isValid));
    foil_mask(isValid==0) = 0; % from locations where vectors are invalid, build a mask matrix
    
    gappy_mask = ones(size(u));
    gappy_mask(corrl<0.4) = 0; % create a gappy mask where the correlation is under the desired threshold
    
    [u_interp, v_interp, G, mse, iter] = GappyPOD_test(u, v, gappy_mask, foil_mask);

end

%% Generating the gif

gifname = '20221006_TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=-120.gif';


for ii = 1:frame % 23 is the number of frames per cycle (14... Hz/0.6.. Hz)
    figname = ['fig' num2str(ii) '.tif'];
    Z = imread(figname);
    [imind,cm] = rgb2ind(Z,256); 
    
    if ii == 1 
        imwrite(imind,cm,gifname,'gif','Loopcount',inf,'DelayTime',1/15); 
    else 
        imwrite(imind,cm,gifname,'gif','WriteMode','append','DelayTime',1/15); 
    end
end

%% save compiled data

filename = ('TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=-120.mat');
save(filename,'compiled_data');

