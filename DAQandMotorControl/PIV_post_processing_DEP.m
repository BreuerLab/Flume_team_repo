%% PIV post-processing - DEPRECATED
% 20220902

clear;

addpath(genpath('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\Libraries'));

% Variable contents:

% dat(:,1) = x-coord [m]
% dat(:,2) = y-coord [m]
% dat(:,3) = u-velocity [m/s]
% dat(:,4) = v-velocity [m/s]
% dat(:,5) = |U| magnitude [m/s]
% dat(:,6) = du/dx [1/s]
% dat(:,7) = du/dy [1/s]
% dat(:,8) = dv/dx [1/s]
% dat(:,9) = dv/dy [1/s]
% dat(:,10) = z-vorticity (dv/dx - du/dy) [1/s]
% dat(:,11) = |W| vorticity magnitude [1/s]
% dat(:,12) = Divergence 2D (du/dx - dv/dy) [1/s]
% dat(:,13) = Swirling strength 2D [1/s^2]
% dat(:,14) = Correlation value
% dat(:,15) = isValid

tic;

%% load directory

directory = '\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\tandem_test_08122022\08122022_alpha=15_p3=70_h3=0.8_ph=-120_f=0.625\SideBySide_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU_01\Export';
cd(directory); % load directory

%% loading force data

load('20220812_PIV_tandemFoils_a15_p3=70_h3=08_ph=-120_f=625A3E.mat'); % force data

p2 = out(:,3); % pitch leading foil [rad]
h2 = out(:,4)*1000; % heave leading foil [mm]
p3 = out(:,5); % pitch trailing foil [rad]
h3 = out(:,6)*1000; % heave trailing foil [mm]

[foil, ~, ~] = foils_database(foiltype);
chord = foil.chord; % [m]

%% loading PIV data

file = importdata('BM0001.dat');

data_no_x = 204; % from "text data" in the file structure (I values)
data_no_y = 152; % from "text data" in the file structure (J values)

[data_x, data_y, data_u, data_v, data_U, data_wz, data_W,...
    data_dudx, data_dudy, data_dvdx, data_dvdy,...
    data_diver, data_swirl, data_correlation, data_cor_isValid] = extract_usable_PIV_data(file, data_no_x, data_no_y, -270, 260, -175, 100, 0);

%% foil profile mask

% leading foil coordinates

foil_separation = 6*chord; % [mm]

foil2_coords = [-235.94, -39.201; % leading edge [mm]
                0, 0; % mid chord [mm] (temporary value)
                -184.28, -5.6579]; % trailing edge [mm]

foil2_coords(2,1) = foil2_coords(3,1) + (foil2_coords(1,1) - foil2_coords(3,1))/2; % x mid-chord [mm]
foil2_coords(2,2) = foil2_coords(3,2) + (foil2_coords(1,2) - foil2_coords(3,2))/2; % y mid-chord [mm]

f2_xpos = foil2_coords(2,1); % for future use



foil2_angle0 = atan((foil2_coords(1,2) - foil2_coords(3,2))/(foil2_coords(1,1) - foil2_coords(3,1))); % starting angle of the leading foil based on the coordinates for the initial position of the leading and trailing edges
[~, I0] = min(abs(p2 - foil2_angle0));  % time instance at which we match the beginning of the force and piv data (NOTE: not the actual moment at which they are matched, but correct within the cycle)

% figure(2); hold on; % time alignment visulization
% plot(out(:,3),'b');
% xline(I0,'r','LineWidth',2); hold off;

frame_displ = foil2_coords(2,2) - h2(I0); % displacement between the encoder frame of reference and the piv data frame of reference

% trailing foil coordinates

foil3_le = [(foil2_coords(2,1) + foil_separation) - cos(p3(I0))*chord/2,    frame_displ + h3(I0) - sin(p3(I0))*chord/2]; % trailing foil leading edge
foil3_mc = [(foil2_coords(2,1) + foil_separation),                          frame_displ + h3(I0)]; % trailing foil mid-chord
foil3_te = [(foil2_coords(2,1) + foil_separation) + cos(p3(I0))*chord/2,    frame_displ + h3(I0) + sin(p3(I0))*chord/2]; % trailing foil trailing edge

foil3_coords = [foil3_le; % leading edge [mm]
                foil3_mc; % mid chord [mm]
                foil3_te]; % trailing edge [mm]

f3_xpos = foil3_coords(2,1); % for future use

%% Plotting

% figure(1)
% 
% contourf(data_x, data_y, data_wz, 'LevelStep', 1, 'LineStyle', 'none'); hold on;
% colorbarpwn(-50, 50, 'wrs', 0.05);
% 
% plot(foil2_coords(:,1), foil2_coords(:,2), 'k', 'LineWidth', 8); % plot leading foil
% plot(foil3_coords(:,1), foil3_coords(:,2), 'k', 'LineWidth', 8); % plot trailing foil
% axis equal;
% hold off;

%% record first data set to compiled data

compiled_data.timeStep = I0; % index in the force data to which the current piv frame corresponds to
compiled_data.angle = [p2(I0), p3(I0)];
compiled_data.foil2_coords = foil2_coords;
compiled_data.foil3_coords = foil3_coords;
compiled_data.x = data_x;
compiled_data.y = data_y;
compiled_data.u = data_u;
compiled_data.v = data_v;
compiled_data.U = data_U;
compiled_data.vort = data_wz;
compiled_data.W = data_W;
compiled_data.dudx = data_dudx;
compiled_data.dudy = data_dudy;
compiled_data.dvdx = data_dvdx;
compiled_data.dvdy = data_dvdy;
compiled_data.diver = data_diver;
compiled_data.swirl = data_swirl;
compiled_data.corrl = data_correlation;
compiled_data.isValid = data_cor_isValid;

%% loop over subsequent frames

dt = 1/15; % time between each frame (1/piv_Hz)
srate = 1000; % daq sampling frequency

t_step = round(dt*srate)-1; % time step in the daq for each frame of piv
dataNum = 1200; % total number of piv frames taken

for j = 2:dataNum
    
    tStep = I0 + t_step*(j-1);
    
    name = ['BM', num2str(j,'%04.0f') '.dat'];
    
    file = importdata(name);
    
    [data_x, data_y, data_u, data_v, data_U, data_wz, data_W,...
    data_dudx, data_dudy, data_dvdx, data_dvdy,...
    data_diver, data_swirl, data_correlation, data_cor_isValid] = extract_usable_PIV_data(file, data_no_x, data_no_y, -270, 260, -175, 100, 0);
    
    foil2_le = [f2_xpos - cos(p2(tStep))*chord/2,    frame_displ + h2(tStep) - sin(p2(tStep))*chord/2];
    foil2_mc = [f2_xpos,                            frame_displ + h2(tStep)];
    foil2_te = [f2_xpos + cos(p2(tStep))*chord/2,    frame_displ + h2(tStep) + sin(p2(tStep))*chord/2];
    
    foil2_coords = [foil2_le; % leading edge [mm]
                    foil2_mc; % mid chord [mm] (temporary value)
                    foil2_te]; % trailing edge [mm]
    
    foil3_le = [f3_xpos - cos(p3(tStep))*chord/2,    frame_displ + h3(tStep) - sin(p3(tStep))*chord/2];
    foil3_mc = [f3_xpos,                            frame_displ + h3(tStep)];
    foil3_te = [f3_xpos + cos(p3(tStep))*chord/2,    frame_displ + h3(tStep) + sin(p3(tStep))*chord/2];
    
    foil3_coords = [foil3_le; % leading edge [mm]
                    foil3_mc; % mid chord [mm] (temporary value)
                    foil3_te]; % trailing edge [mm]
    
    %% Plotting
    
%     figure(1)
%     
%     contourf(data_x, data_y, data_wz, 'LevelStep', 1, 'LineStyle', 'none'); hold on;
%     colorbarpwn(-50, 50, 'wrs', 0.05);
%     
%     plot(foil2_coords(:,1), foil2_coords(:,2), 'k', 'LineWidth', 8); % plot leading foil
%     plot(foil3_coords(:,1), foil3_coords(:,2), 'k', 'LineWidth', 8); % plot trailing foil
%     
%     axis equal;
%     hold off; 

    %% record current data set to compiled data
    
    compiled_data(j).timeStep = tStep; % index in the force data to which the current piv frame corresponds to
    compiled_data(j).angle = [p2(tStep), p3(tStep)];
    compiled_data(j).foil2_coords = foil2_coords;
    compiled_data(j).foil3_coords = foil3_coords;
    compiled_data(j).x = data_x;
    compiled_data(j).y = data_y;
    compiled_data(j).u = data_u;
    compiled_data(j).v = data_v;
    compiled_data(j).U = data_U;
    compiled_data(j).vort = data_wz;
    compiled_data(j).W = data_W;
    compiled_data(j).dudx = data_dudx;
    compiled_data(j).dudy = data_dudy;
    compiled_data(j).dvdx = data_dvdx;
    compiled_data(j).dvdy = data_dvdy;
    compiled_data(j).diver = data_diver;
    compiled_data(j).swirl = data_swirl;
    compiled_data(j).corrl = data_correlation;
    compiled_data(j).isValid = data_cor_isValid;
    
end

[interpolated_data, G] = GappyPOD_PIVpostProc(compiled_data);

toc;

%% this section might only be temporary

for j = 1:dataNum
    
    coordinates_data(j).timeStep = compiled_data(j).timeStep;
    coordinates_data(j).angle = compiled_data(j).angle;
    coordinates_data(j).foil2_coords = compiled_data(j).foil2_coords;
    coordinates_data(j).foil3_coords = compiled_data(j).foil3_coords;
    coordinates_data(j).x = compiled_data(j).x;
    coordinates_data(j).y = compiled_data(j).y;
    
    flowfield_data(j).u = compiled_data(j).u;
    flowfield_data(j).v = compiled_data(j).v;
    flowfield_data(j).vort = compiled_data(j).vort;
    flowfield_data(j).dudx = compiled_data(j).dudx;
    flowfield_data(j).dudy = compiled_data(j).dudy;
    flowfield_data(j).dvdx = compiled_data(j).dvdx;
    flowfield_data(j).dvdy = compiled_data(j).dvdy;
    
    other_data(j).diver = compiled_data(j).diver;
    other_data(j).swirl = compiled_data(j).swirl;
    other_data(j).corrl = compiled_data(j).corrl;
    other_data(j).isValid = compiled_data(j).isValid;
    
end

save('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\tandem_test_08122022\08122022_alpha=15_p3=70_h3=0.8_ph=-120_f=0.625\SideBySide_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU_01\Export\20220812_COORDINATES_a15_p3=70_h3=08_ph=-120_f=625_A3E','coordinates_data');
save('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\tandem_test_08122022\08122022_alpha=15_p3=70_h3=0.8_ph=-120_f=0.625\SideBySide_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU_01\Export\20220812_FLOWFIELD_a15_p3=70_h3=08_ph=-120_f=625_A3E','flowfield_data');
save('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\tandem_test_08122022\08122022_alpha=15_p3=70_h3=0.8_ph=-120_f=0.625\SideBySide_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU_01\Export\20220812_OTHERDATA_a15_p3=70_h3=08_ph=-120_f=625_A3E','other_data');
    

%% plotting interpolated vorticity

for n = 1000:1050%length(interpolated_data)
    
    figure(1)
    
%     U = sqrt(interpolated_data(n).ug.^2 + interpolated_data(n).vg.^2);
    
    contourf(compiled_data(n).x, compiled_data(n).y, compiled_data(n).vort, 'LevelStep', 1, 'LineStyle', 'none'); hold on;
    colorbarpwn(-50, 50, 'wrs', 0.05);
%     contourf(interpolated_data(n).xg, interpolated_data(n).yg, U, 'LineStyle', 'none'); hold on;
%     colorbarpwn(0, 0.5, 'wrs', 0.05);
%     quiver(interpolated_data(n).xg, interpolated_data(n).yg, interpolated_data(n).ug, interpolated_data(n).vg, 2); hold on;
    
    axis equal;
    
    foil2_coords = compiled_data(n).foil2_coords;
    foil3_coords = compiled_data(n).foil3_coords;
    
    plot(foil2_coords(:,1), foil2_coords(:,2), 'k', 'LineWidth', 8); % plot leading foil
    plot(foil3_coords(:,1), foil3_coords(:,2), 'k', 'LineWidth', 8); % plot trailing foil
    
    hold off;
    
end

% save('C:\Users\ehandyca\Desktop\20220812_COMPILED_PIV_tandemFoils_a15_p3=70_h3=08_ph=-120_f=625_A3A.mat','compiled_data','interpolated_data'); % <--- so this doesn't work if the file is too large

%% cycle averaging --> this aint working yet and it really requires simplification

% fpc = 24; % piv frames-per-cycle
% tot_cycs = dataNum/fpc; % total number of cycles
% 
% angle = NaN(fpc,2,tot_cycs-1);
% vorticity = NaN(fpc,tot_cycs-1);
% vort_tall = NaN(186*97,dataNum);
% 
% for j = 1:dataNum
%     vort_tall(:,j) = reshape(flowfield_data(j).vort,[],1);
% end
% 
% i = 0;
% 
% for n = 1:fpc:dataNum-fpc % cycles
%     i = i + 1;
%     angle_temp = NaN(fpc,2);
%     vort_temp = NaN(186*96,fpc);
%     k = 0;
%     for j = n:n+fpc-1 % frames
%         k = k + 1;
%         angle_temp(k,:) = coordinates_data(j).angle;
%         vort_temp(:,k) = vort_tall(:,j);
%     end
%     angle(:,:,i) = angle_temp;
%     vorticity(:,i) = vort_temp;
% end
% 
% angle_avg = mean(angle,3);
% vort_avg_tall = mean(vorticity,3);
% 
% vort_avg = NaN(size(flowfield_data(1).vort));
% 
% for j = 1:dataNum
%     vort_avg(:,:,j) = reshape(vort_tall(:,j),data_no_x,data_no_y);
% end
