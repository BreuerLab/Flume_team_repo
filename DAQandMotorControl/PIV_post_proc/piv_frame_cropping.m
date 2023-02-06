%% Frame Cropping

function [range_x, range_y, window_index] = piv_frame_cropping(x,y,u,v,I,J)

% x - x coordinates matrix
% y - y coordinates matrix
% u - x velocity
% v - y velocity
% returns index ranges for the cropped-in data

fig_name = 'Draw desired window';
figure('name',fig_name);
quiver(x,y,u,v,'k'); hold on;
fprintf('\n Click and drag to draw a square that encompasses the desired PIV field. \n\n')
fov = drawrectangle; hold off;
window = fov.Position; % [xmin, ymin, width, height]

[~, Ix0] = min(abs(x(:,1) - window(1))); % finds the index of the approximate desired location within the x and y matrices
[~, Ix1] = min(abs(x(:,1) - (window(1)+window(3))));
[~, Iy0] = min(abs(y(1,:) - window(2)));
[~, Iy1] = min(abs(y(1,:) - (window(2)+window(4))));

range_x = Ix0:Ix1;
range_y = Iy1:Iy0;

close(fig_name)

fig_name = 'Cropped-in data';
figure('name',fig_name);
quiver(x,y,u,v,'k'); hold on;
quiver(x(range_x,range_y),y(range_x,range_y),u(range_x,range_y),v(range_x,range_y),'r'); hold off;
fprintf('Press any key to continue\n\n')
pause()
close(fig_name)

index_matrix = ones(size(x)); % clone x matrix and fill with 1s
index_matrix(range_x,range_y) = nan; % make values of selected range into NaN
temp = reshape(index_matrix,[numel(index_matrix),1]); % reshape matrix into a vector
window_index = isnan(temp); % get a logical vector to extract the values in the chosen range from any matrix

end