%% Plot flowfield

clear;

piv_freq = 14.9316; % frequency of piv data
frame = 23; % piv frames per cycle

%% Loading data

% load compiled data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=16_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU_01\Export\TandemFoil_aT4=0.16_p3=75_h3=0.7c_ph=-120.mat');
% load force data file
load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\20221006_alpha=16_p3=75_h3=0.7_ph=-120_A3E.mat');

%% Force calculations

[kin, par, foil] = extract_measurements_2rigs(foiltype, Prof_out_angle, out, samplerate, transientcycs, foil_separation, flume_height);
out = filter_motor_noise_gromit(out, freq, samplerate, 20);
out = filter_motor_noise_wallace(out, freq, samplerate, 20);
res = calculate_forces(par, kin, out);

[     ~,          ~, CL3_cyc] = cycle_avg_data(kin.p3_comm, res.CL3, samplerate, 1); % LiftC3
[toverT, pitch_cyc3, CP3_cyc] = cycle_avg_data(kin.p3_comm, (res.CPH3+res.CPP3), samplerate, 1); % PowerC3

%% Stupid alignment

for ii = 1:frame
    
    figure('Name', 'Foil Coordinates', 'WindowState', 'maximized');
    title(fprintf(['Frame: ',num2str(ii)]));
    contour(compiled_data(ii).x, compiled_data(ii).y, compiled_data(ii).isValid);
    colormap('bone');
    fprintf(['\n\nSelect approximate foil center coordinates for the LEADING and TRAILING foil in that order.\n\n When done with the current frame press [Enter]\n\nCurrent frame: ',num2str(ii),'\n']);
    leading = drawpoint('Label','Leading Foil');
    trailing = drawpoint('Label','Trailing Foil');
    coords_le(ii,:) = leading.Position;
    coords_tr(ii,:) = trailing.Position;
    close 'Foil Coordinates';
    
end

center_le(1) = mean(coords_le(:,1)); % get the centerline of the leading foil (x)
center_le(2) = (max(coords_le(:,2)) + min(coords_le(:,2)))/2; % get the centerline of the leading foil (y)

center_tr(1) = mean(coords_tr(:,1)); % get the centerline of the trailing foil (x)
center_tr(2) = (max(coords_tr(:,2)) + min(coords_tr(:,2)))/2; % get the centerline of the trailing foil (y)

center_axis(1) = (mean(center_le(:,1)) + mean(center_tr(:,1)))/2; % get the origin of all axes (x)
center_axis(2) = (mean(center_le(:,2)) + mean(center_tr(:,2)))/2; % get the origin of all axes (y)

center_le_d = center_le - center_axis; % correct the leading centerline
center_tr_d = center_tr - center_axis; % correct the trailing centerline

disp_angle = atan(center_tr_d(2)/center_tr_d(1)); % get the displacement angle

R = [cos(disp_angle), -sin(disp_angle); sin(disp_angle), cos(disp_angle)]; % rotation matrix

center_le_p = center_le_d*R; % rotates the translated center coordinates
center_tr_p = center_tr_d*R; % rotates the translated center coordinates

figure(1)
plot([center_le(1),center_tr(1)],[center_le(2),center_tr(2)],'--r'); hold on;
plot([center_le_p(1),center_tr_p(1)],[center_le_p(2),center_tr_p(2)],'--k'); hold off;

%%

x = compiled_data(1).x;
y = compiled_data(1).y;
isValid = compiled_data(1).isValid;

figure('Name', 'Foil Coordinates', 'WindowState', 'maximized');
title('Front Foil Leading and Trailing Edge');
contour(x,y,isValid); hold on;
colormap('bone');

fprintf('\nSelect coordinates of the Leading Edge and Trailing Edge of the FRONT FOIL, in that order.\n\n When done press [Enter]\n');
% pick leading and trailing edge of the leading foil
l_e = drawpoint('Label','Leading Edge'); l_e = l_e.Position; % pick the leading edge and only use the extracted position
t_e = drawpoint('Label','Trailing Edge'); t_e = t_e.Position; % pick the leading edge and only use the extracted position

plot([l_e(1),t_e(1)],[l_e(2),t_e(2)],'k','LineWidth',8);
pause(1);
close 'Foil Coordinates';

foil2_angle0 = rad2deg(atan((l_e(2)-t_e(2))/(l_e(1)-t_e(1)))); % calculate the angle of the trailing foil in order to match with force data
% [ ~, pitch3_cyc, ~] = cycle_avg_data(kin.p3_comm, res.CL3, samplerate, 1); % cycle-avg the trailing pitch

% % %

% foil2_angle0 = atan((foil2_coords(1,2) - foil2_coords(1,6))/(foil2_coords(1,1) - foil2_coords(1,5))); % starting angle of the leading foil based on the coordinates for the initial position of the leading and trailing edges
p2_temp = -out(:,3);
transient_time = round((1/freq)*transientcycs*samplerate); % time taken up by the transients
p2_temp(1:transient_time) = NaN; % make transient data irrelevant
p2_temp(end-transient_time:end) = NaN; % make transient data irrelevant

p2_temp(gradient(p2_temp) > 0) = NaN; % [CHANGE CONDITIONAL TO SEARCH FOR ANGLE DURING THE DESIRED STROKE]
[~, I0] = min(abs(p2_temp - foil2_angle0));  % time instance at which we match the beginning of the force and piv data (NOTE: not the actual moment at which they are matched, but correct within the cycle)

figure(); hold on; % time alignment [FOR DEBUGGING]
plot(rad2deg(out(:,3)),'b'); grid on;
plot(rad2deg(out(:,5)),'r');
plot(out(:,24)*max(rad2deg(out(:,3))),'g','linewidth',2);
xline(I0,'r','LineWidth',2); hold off;
xlim([(I0-2000),(I0+2000)])

% % %

% [~,I0] = min(abs(p3_0 - mean(pitch_cyc3))); % find the starting time index of the cyc-avg data of the TRAILING FOIL


% [~, I0] = min(1-out(:,24)); % find the time index


tstp_per_frm = samplerate/piv_freq; % time steps per piv frame
t = 0:(1/par.srate):(1/par.srate)*(length(out)-1);

%% Plotting

for ii = 1:frame
    
    tStep = I0 + round((ii-1)*tstp_per_frm);
%     [~, IC] = min(abs(mean(pitch_cyc3) - Prof_out_angle(tStep,5)));
    
    x = compiled_data(ii).x;
    y = compiled_data(ii).y;
    % Transform the coordinates
    xp = x.*cos(disp_angle) - y.*sin(disp_angle) - center_axis(1);
    yp = x.*sin(disp_angle) + y.*cos(disp_angle) - center_axis(2);
    
    vort = compiled_data(ii).vort;
    timeStep = compiled_data(ii).timeStep;
%     foil2_coords = compiled_data(ii).foil2_coords;
%     foil3_coords = compiled_data(ii).foil3_coords;
    
    figure(1)
    
%     subplot(2,2,[1,2])
    contourf(xp, yp, vort, 'LineStyle', 'none', 'levelstep', 5); hold on;
    colorbarpwn(-30,30,'log','wrs',0.1)
    
    plot(foil2_coords(ii,[1,3,5]), foil2_coords(ii,[2,4,6]), 'k', 'LineWidth', 10); % plot leading foil
    plot(foil3_coords(ii,[1,3,5]), foil3_coords(ii,[2,4,6]), 'k', 'LineWidth', 10); hold off; % plot trailing foil
    
    axis equal;
    
%     subplot(2,2,3)
%     yyaxis left
%     plot(toverT,mean(CL3_cyc),'k','LineWidth',2); hold on;
%     yyaxis right
%     plot(toverT,mean(pitch_cyc3),'--g','LineWidth',1.5)
%     xline(toverT(tStep),'r','LineWidth',2); hold off;
    
%     subplot(2,2,4);
    % plot the pitch profile per one cycle
%     plot(t(IC:IC+tstp_per_frm*frame), Prof_out_angle(IC:IC+tstp_per_frm*frame,5), 'b', 'LineWidth', 1.5); hold on;
    % plot the point at which the flow snapshot is taken
%     xline(t(I0),'r','LineWidth',2); hold off;
end