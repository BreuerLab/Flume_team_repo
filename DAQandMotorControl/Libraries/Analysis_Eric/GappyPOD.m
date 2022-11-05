%% testing section

clear;

addpath(genpath('\\lrs.brown.edu\research\ENG_Breuer_Shared\ehandyca\Flume_team_repo\DAQandMotorControl\Libraries'));

frames = 23;

for dataset = 1:frames
    foldername = ('\\lrs.brown.edu\research\ENG_Breuer_Shared\group\Flume PIV\Flume PIV\tandem foils\20221006_wake-foil_interactions_alphaT4\20221006_alpha=16_p3=75_h3=0.7_ph=-120_f=0.6492\SideBySide_PIV_MPd(4x32x32_75%ov_ImgCorr)_GPU_01\Export');
    filename = ['B', num2str(dataset,'%04.0f'), '.dat'];
    [A,delimiterout,headerlinesout] = importdata(fullfile(foldername,filename)); % load data from file
    
    x(:,dataset) = A.data(:,1);
    y(:,dataset) = A.data(:,2);
    u(:,dataset) = A.data(:,3);
    v(:,dataset) = A.data(:,4);
    corr(:,dataset) = A.data(:,12); % correlation
    isValid(:,dataset) = A.data(:,15); % valid/invalid vectors
    
end

save('testing_piv_data','x','y','u','v','corr','isValid', 'frames');

%%

clear

load('testing_piv_data.mat');

foil_mask = ones(size(isValid));
foil_mask(isValid==0) = 0; % from locations where vectors are invalid, build a mask matrix

gappy_mask = ones(size(u));
gappy_mask(corr<0.4) = 0; % create a gappy mask where the correlation is under the desired threshold

u_plot = u; u_plot(gappy_mask==0) = nan;
v_plot = v; v_plot(gappy_mask==0) = nan;

% for d = 1:frames
%     figure(1); % [VISUALIZATION]
%     quiver(x(1:166953,d),y(1:166953,d),u_plot(1:166953,d),v_plot(1:166953,d),'k');
%     title('gapped flow')
% end

[u_interp, v_interp, G, mse, iter] = GappyPOD_test(u, v, gappy_mask, foil_mask);

figure(2)
bar(mse, 'FaceColor', 'k', 'EdgeColor', 'none');

u_interp_plot = u_interp;
v_interp_plot = v_interp;
u_interp_plot(gappy_mask == 1) = nan;
v_interp_plot(gappy_mask == 1) = nan;
for d = 1:frames
    figure(3); % [VISUALIZATION]
    quiver(x(:,d),y(:,d),u_plot(:,d),v_plot(:,d),'k'); hold on;
    quiver(x(:,d),y(:,d),u_interp_plot(:,d),v_interp_plot(:,d),'r'); hold off;
    
    title('reconstructed flow');
end

%% Gappy POD
% 20221014, Ed Handy

function [u_interp, v_interp, G, mse, iter] = GappyPOD_test(u, v, gappy_mask, foil_mask)

% GAPPYPOD Interpolates vectors using the Gappy POD method.
% 
% Inputs:
% u, v : spatio-temporal matrices where rows = spatial dimension, columns = time dimension,  u[X,t]
% mask : logical gappy mask with the same size as u and v, where values 1 are valid quatities and 0 are invalid quantities (gaps)
%
% Outputs:
% u_interp, v_interp : spatio temporal matrices with interpolated vectors filled-in
% G : gappyness percentage


tic;

%% Preparing the Gappy dataset

u(gappy_mask==0) = nan;
v(gappy_mask==0) = nan;

G = (sum(gappy_mask(:) == 0)/numel(gappy_mask))*100; % gappyness

xi = [u; v]; % original gappy field, x and y components stacked

gaps = isnan(xi); % gaps in the complete set (all spatial and temporal data)

%% First iteration (temporal mean guess)

% Indeces: i : rows(space), k : columns(time)

for k = 1:size(xi,2) % loop over each column (time)
    for i = 1:size(xi,1) % loop over each row (space)
        if isnan(xi(i,k)) % if a gap is encountered
            index = ~isnan(xi(i,:)); % find all the known (not gaps) temporal data in the selected location
            if index == zeros(1,size(xi,2))
                xi(i,k) = 0; % if there are no real values in the temporal data, make value zero
            else
                xi(i,k) = mean(xi(i,index)); % replace the gap with the temporal mean at that spatial location
            end
        end
    end
end

xi_1 = zeros(size(xi)); % gappy dataset, with zeros instead of NaN values

for i = 1:size(xi,1)
    for k = 1:size(xi,2)
        if gaps(i,k) == 0
            xi_1(i,k) = xi(i,k);
        end
    end
end

%% Performing POD

iterations = size(u,2)-1; % maximum number of iterations (must be < maximum num of modes obtainable from the dataset)
if iterations >= 200
    iterations = 200; % just to not make too massive the number of iterations
end

iter = 0; % initialize the interation variable
mse = 1; % mean square error
Err = zeros(size(xi,1),size(xi,2),iterations); % initialize the error matrix
Err_rel = zeros(size(xi,1),size(xi,2),iterations); % initialize the relative error matrix
P = 1; % starting number of modes used (for the initial iteration use only two modes {P = P0 + 1})

% for iter = 1:iterations % loop for a certain number of iterations
while mse > 1e-6 % loop until the error threshold is reached
    
    iter = iter + 1;
    
    [U, ~, ~] = svd(xi, 'econ'); % SVD to get spatial modes (U) and temporal eigenmodes (S*V')

    phi = U; % spatial modes

    % %% Truncation

    % From Venturi and Karniadakis:

    P = P + 1; % increase the admitted spatial modes every iteration

    % %% Reconstuction error minimization via Least Squares
    % Since we only want to consider the error between the real values and
    % the approximated values, we'll make all interpolated values zero for
    % the error calculation.

    phi_tild = phi(:,1:P); % truncated spatial modes matrix

    % Reconstruction matrix and projected field vector

    K = NaN(P,P,size(xi,2)); % reconstruction matrix ([K]ij = (phi_i,phi_j)k --> remember k is the time index

    for k = 1:size(xi,2)
        for i = 1:P
            for j = 1:P
                K(i,j,k) = dot(phi_tild(:,i).*~gaps(:,k),phi_tild(:,j));
            end
        end
    end

    f = phi_tild\(xi_1.*~gaps); % phi * a = u_tild(masked)

    a = NaN(P,size(xi,2));

    for k = 1:size(xi,2)
        a(:,k) = K(:,:,k)\f(:,k);
    end

    % %% Reconstruct the field

    xi_hat = phi_tild*a; % reconstructed field using the new temporal coefficients

    % %% Reconstruction Error calculation

    for i = 1:size(xi,1)
        for k = 1:size(xi,2)
            if gaps(i,k) == 0 % only for ungapped datapoints
                Err(i,k,iter) = (xi_hat(i,k) - xi_1(i,k))^2; % difference between the reconstructed value and the actual value
                Err_rel(i,k,iter) = (xi_hat(i,k) - xi_1(i,k))^2/(xi_1(i,k))^2; % relative error between the reconstructed value and the actual value
            end
        end
    end

    mse(iter) = mean(Err_rel(:,:,iter),'all'); % mean of the l2-error (this might be BS)
    
    % %% Overwrite old matrix with new interpolated matrix

    for i = 1:size(xi,1)
        for k = 1:size(xi,2)
            if gaps(i,k) == 1
                xi(i,k) = xi_hat(i,k);
            end
        end
    end
    
    if iter == iterations
        break % stop algorithm if iterations exceede the maximum number of iterations
    end
    
end

%% Reconstruction of the Flow Field

% Extracting the interpolated values of each component and adding the foil mask
u_interp = xi(1:end/2,:).*foil_mask;
v_interp = xi(end/2+1:end,:).*foil_mask;

toc;

end