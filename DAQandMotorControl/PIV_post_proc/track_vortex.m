%% Testing
% clear;
% load('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\DATA_main_repo\20230120_Nice_Test_PIV_flowfields\TandemFoil_aT4=0.68_p3=75_h3=0.8c_ph=-120.mat');
% 
% for i = 1:length(compiled_data)
%     x(:,:) = compiled_data(i).x;
%     y(:,:) = compiled_data(i).y;
%     u(:,:,i) = compiled_data(i).u;
%     v(:,:,i) = compiled_data(i).v;
% end
% 
% u = u(:,:,[12:23,1:11]); % reorder snapshots so that the vortex starts on the mid-heave position of the foil
% v = v(:,:,[12:23,1:11]);
% 
% [Q, vortex, strength] = track_vortex_t(x, y, u, v, 1/14.9316, 23, 200);

%% Vortex core tracking
% Upadated 02/13/2023, E Handy
%
% Tracks a vortex formed and shed by a hydrofoil based on the Q-criterion and an initial interrogation window.
% 
% ----------------------------------------------------------------------------------------------------------------------------
% % Inputs (required):
%   x - x coordinates matrix for one snapshot [n,m]
%   y - y coordinates matrix for one snapshot [n,m]
%   u - x component of flow velocity for each snapshot [n,m,N]
%   v - y component of flow velocity for each snapshot [n,m,N]
%   dt - time step between each frame in the dataset (1/piv_frequency)
% 
% ----------------------------------------------------------------------------------------------------------------------------
% % Inputs (optional):
%   frames - number of frame to be analyzed (usually frames per cycle)
%   threshold - number of maximum Q-values to be used to locate vortex center within the interrogation window (default is 150)
%
% ----------------------------------------------------------------------------------------------------------------------------
% % Outputs:
%   Q - Q-field of the complete dataset
%   vortex - tracked vortex coordinates
%   strength - mean Q-value of the vortex within the frame of the interrogation window and the threshold
%              NOTE: Vortex strength metric should be redefined so that the circulation is used instead of the Q_value
% ----------------------------------------------------------------------------------------------------------------------------
function [Q, vortex, strength] = track_vortex(x, y, u, v, dt, frames, threshold)

if ~exist('frames','var')
    frames = size(u,3); % if not given, all frames will be analyzed
end

if ~exist('threshold','var')
    threshold = 150; % if not given, only the maximum 150 Q-values of the interrogation window will be considered
end

ds = x(1,1) - x(1,2); % distance between each vector

% initialize variables
Q = nan(size(u)); % Q-field
strength = nan(frames,1);
core_pos_x = nan(frames,1); % vortex core positions x
core_pos_y = nan(frames,1); % vortex core positions y

fig1 = figure();
fig1.Position = [500, 180, 1300, 700]; % for a 1080p monitor

for i = 1:frames
    
    [dudx, dudy] = gradient(u(:,:,i), ds); % calculate gradients
    [dvdx, dvdy] = gradient(v(:,:,i), ds); % calculate gradients
    
    % for a 2D flowfield
    Qi = (dudx.*dvdy - dvdx.*dudy); % calculate the current Q-field
    Q_neg = zeros(size(Qi));
    Q_neg(Qi>1) = Qi(Qi>1);
    Q_neg(log(Qi)<13) = 0; % makes visualization easier but depends on the cutoff value
    Qi = log(Q_neg);
    Qi(Qi == -Inf) = 0;
    Q(:,:,i) = Qi;
    
    contourf(x, y, Q(:,:,i), 'linestyle', 'none'); hold on; % plot of the Q-field
    colormap(flipud(hot));
    axis equal;
    title('Vortex location');
    
    if i == 1 % only for the starting location of the vortex    
        fprintf('\n Click and drag to draw a square that encompasses the vortex core. \n IMPORTANT: Try to account for growth and dissipation by selecting a large area\n\n')
        fov = drawrectangle; % select region that encompasses the vortex core
        window = fov.Position; % [xmin, ymin, width, height]
    else
        width = window(3); % take the width of the last snapshot's selected window
        height = window(4); % take the height of the last snapshot's selected window
        window = [core_pos_x(i-1) + (mean(u(range_x,range_y,i-1),'all') + mean(sqrt(u(:,:,i-1).^2+v(:,:,i-1).^2),'all'))*dt-width/2,...
            core_pos_y(i-1)+mean(v(range_x,range_y,i-1),'all')*dt-height/2,...
            width, height]; % new interrogation window based on the previous vortex location and the mean local velocity
    end
    
    rectangle('Position', window, 'EdgeColor', 'g', 'LineStyle', '--', 'LineWidth', 2); % plot location of interrogation window
    
    [~, Ix0] = min(abs(x(:,1) - window(1))); % finds the index of the approximate desired location within the x and y matrices
    [~, Ix1] = min(abs(x(:,1) - (window(1)+window(3))));
    [~, Iy0] = min(abs(y(1,:) - window(2)));
    [~, Iy1] = min(abs(y(1,:) - (window(2)+window(4))));

    range_x = Ix0:Ix1; % index range of the selected window
    range_y = Iy1:Iy0;
    
    values = maxk(Q(range_x,range_y,i),threshold); % returns the <threshold> largest Q values in the field
    values = sort(reshape(values,[],1),'descend'); % reshape into a vector and sort from largest to smallest
    
    temp = zeros(size(Qi)); % temporary array to only plot contours within the interrogation window
    temp(range_x,range_y) = Qi(range_x,range_y); % replace only Q-values within the interrogation window
    plot_level = round(values(threshold)); % define the Q-value contour of interest
    contour(x, y, temp, [plot_level, plot_level], 'linewidth', 2, 'LineColor', [0.28 0.75 0.95]); % overlay contour of interest
    
    for j = range_x % loops only over the interrogation window
        for k = range_y
            if Qi(j,k) < values(threshold) % only perform this on the temporary in-loop calculated Q-field
                Qi(j,k) = 0; % make values smaller than the 50th largest Q-value equal to zero
            end
        end
    end
     
    core_pos_x(i) = sum(Qi(range_x,range_y).*x(range_x,range_y),'all')/sum(Qi(range_x,range_y),'all'); % center of mass - x
    core_pos_y(i) = sum(Qi(range_x,range_y).*y(range_x,range_y),'all')/sum(Qi(range_x,range_y),'all'); % center of mass - y
    
    strength(i) = mean(exp(Qi(range_x,range_y)),'all'); % calculate average vortex strength in the interrogation window
    
    plot(core_pos_x(i), core_pos_y(i), 'go', 'MarkerSize', 13, 'linewidth', 1.5); % plot located vortex core location
    plot(core_pos_x(i), core_pos_y(i), 'g+', 'MarkerSize', 8, 'linewidth', 1.5); hold off;
    
    pause(0.5);
    
end

vortex = [core_pos_x, core_pos_y]; % output vortex coordinates

% Plotting the trajectory

figure()

subplot(2,1,1)
plot(core_pos_x, core_pos_y, '^', 'markerfacecolor', [0.6350 0.0780 0.1840], 'markeredgecolor', [0.6350 0.0780 0.1840], 'markersize', 12);
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
title('Vortex trajectory','interpreter','latex');
grid on; axis equal;
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');

subplot(2,1,2)
plot(strength,'s-','color',[0,0,0],'linewidth',1.5);
xlabel('Frame','interpreter','latex');
ylabel('$\overline{Q}$','interpreter','latex');
title('Vortex strength','interpreter','latex');
set(gca,'fontname','Times New Roman','fontsize',16,'linewidth',1.5,'TickDir','in');

close(fig1);

end

