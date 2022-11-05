%% Extracts only desired data from the PIV file loaded - 20220906

% returns formatted data in SI units

function [raw_x, raw_y, raw_u, raw_v, raw_wz,...
    raw_dudx, raw_dudy, raw_dvdx, raw_dvdy,...
    raw_diver, raw_swirl, raw_correlation, raw_cor_isValid, matrix_shape] = extract_PIV_data(file, data_no_x, data_no_y, fov_x0, fov_x1, fov_y0, fov_y1, plot_option)

% file - structure from the ".dat" PIV file
% data_no_x - from "text data" in the file structure (I values)
% data_no_y - from "text data" in the file structure (J values)
% fov_x0 - starting x measurement of the desired data fov [mm]
% fov_x1 - stopping x measurement of the desired data fov [mm]
% fov_y0 - starting y measurement of the desired data fov [mm]
% fov_y1 - stopping y measurement of the desired data fov [mm]
% plot_option - plot extracted values [1/0]

% outputs data in SI units

    raw_data_shape = [data_no_x, data_no_y];

    raw_x  = file.data(:,1); raw_x = reshape(raw_x,raw_data_shape);
    raw_y  = file.data(:,2); raw_y = reshape(raw_y,raw_data_shape);
    raw_u  = file.data(:,3); raw_u = reshape(raw_u,raw_data_shape); 
    raw_v  = file.data(:,4); raw_v = reshape(raw_v,raw_data_shape); 
    raw_U  = file.data(:,5); raw_U = reshape(raw_U,raw_data_shape); % optional, unused

    if plot_option == 1 % only plots if plotting option is valid
        figure();
        quiver(raw_x, raw_y, raw_u, raw_v, 'r'); hold on;
    end

    [~, Ix0] = min(abs(raw_x(:,1) - fov_x0)); % finds the index of the approximate desired location within the x and y matrices
    [~, Ix1] = min(abs(raw_x(:,1) - fov_x1));
    [~, Iy0] = min(abs(raw_y(1,:) - fov_y0));
    [~, Iy1] = min(abs(raw_y(1,:) - fov_y1));

    range_x = Ix0:Ix1;
    range_y = Iy1:Iy0; % flipped because in the matrices, y0 is at a larger index than y1

    raw_x = raw_x(range_x,range_y);
    raw_y = raw_y(range_x,range_y);
    raw_u = raw_u(range_x,range_y);
    raw_v = raw_v(range_x,range_y);
    raw_U = raw_U(range_x,range_y); % optional, unused

    if plot_option == 1
        quiver(raw_x, raw_y, raw_u, raw_v, 'b'); hold off;
    end

    raw_wz = file.data(:,6); raw_wz = reshape(raw_wz,raw_data_shape); raw_wz = raw_wz(range_x,range_y);
    raw_W = file.data(:,7); raw_W = reshape(raw_W,raw_data_shape); raw_W = data_W(range_x,range_y); % optional, unused

    raw_dudx = file.data(:,8); raw_dudx = reshape(raw_dudx,raw_data_shape); raw_dudx = raw_dudx(range_x,range_y);
    raw_dudy = file.data(:,9); raw_dudy = reshape(raw_dudy,raw_data_shape); raw_dudy = raw_dudy(range_x,range_y);
    raw_dvdx = file.data(:,10); raw_dvdx = reshape(raw_dvdx,raw_data_shape); raw_dvdx = raw_dvdx(range_x,range_y);
    raw_dvdy = file.data(:,11); raw_dvdy = reshape(raw_dvdy,raw_data_shape); raw_dvdy = raw_dvdy(range_x,range_y);

    raw_diver = file.data(:,12); raw_diver = reshape(raw_diver,raw_data_shape); raw_diver = raw_diver(range_x,range_y);
    raw_swirl = file.data(:,13); raw_swirl = reshape(raw_swirl,raw_data_shape); raw_swirl = raw_swirl(range_x,range_y);

    raw_correlation = file.data(:,14); raw_correlation = reshape(raw_correlation,raw_data_shape); raw_correlation = raw_correlation(range_x,range_y);
    raw_cor_isValid = file.data(:,15); raw_cor_isValid = reshape(raw_cor_isValid,raw_data_shape); raw_cor_isValid = raw_cor_isValid(range_x,range_y);
    
    matrix_shape = shape(raw_x); % shape of the utilized data

end
