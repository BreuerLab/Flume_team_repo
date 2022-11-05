
clear;

force_folder = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20221006_TandemPIV_3alphaRegimes\');
force_filename = ('20221006_alpha=16_p3=75_h3=0.7_ph=-120_A3E.mat'); % force data file

load(fullfile(force_folder,force_filename));

[foil, ~, ~] = foils_database(foiltype);

figure(1);% hold on;

chord = foil.chord;

[~, I0] = min(1-out(:,24));

dt = 1/14.93; % time between each frame (1/piv_Hz)
srate = 1000; % daq sampling frequency

separation = 6*foil.chord;
t_step = round(1*dt*srate); % time step in the daq for each frame of piv
t_step = 40;

for i = I0:t_step:20000%length(out(:,5))
    
%     step = I0 + t_step*(j-1);
    
    % leading foil
    
    p = -out(i,3);
    h = -out(i,4);
    
    % NOTE: I think the addition or subtraction of the leading and trailing
    % edges doesn't affect the trajectory. (i.e., they can swap signs and
    % have correct motion)
    
    x0 = -(chord/2)*cos(p);
    x1 = 0;
    x2 = (chord/2)*cos(p);
    
    y0 = h - (chord/2)*sin(p);
    y1 = h;
    y2 = h + (chord/2)*sin(p);
    
    posx = [x0, x1, x2];
    posy = [y0, y1, y2];
    
    figure(3);
    plot(posx, posy, 'k', 'LineWidth', 8); hold on;

    % trailing foil
    
    p = -out(i,5);
    h = -out(i,6);
    
    x0 = separation - (chord/2)*cos(p);
    x1 = separation;
    x2 = separation + (chord/2)*cos(p);
    
    y0 = h - (chord/2)*sin(p);
    y1 = h;
    y2 = h + (chord/2)*sin(p);

    posx = [x0, x1, x2];
    posy = [y0, y1, y2];
    
    plot(posx, posy, 'k', 'LineWidth', 8); hold off;
    
    axis equal
    
    xlim([-0.1,0.45])
    ylim([-0.15,0.15])
    
end