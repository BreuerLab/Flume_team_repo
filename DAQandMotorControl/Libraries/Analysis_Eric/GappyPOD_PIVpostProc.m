%% ENGN 2912T, Gappy Proper Orthogonal Decomposition - GPOD
% Adapted mainly from Venturi and Karniadakis, 2004

function [interpolated_data, G] = GappyPOD_PIVpostProc(compiled_data)

    parfor n = 1:length(compiled_data)
        x(:,:,n) = compiled_data(n).x;
        y(:,:,n) = compiled_data(n).y;
        u(:,:,n) = compiled_data(n).u;
        v(:,:,n) = compiled_data(n).v;
        vort(:,:,n) = compiled_data(n).vort;
        corrl(:,:,n) = compiled_data(n).corrl;
    end
    
    foil_mask = ones(size(u));
    foil_mask(u==0) = nan; % from locations where there is zero velocity (assuming this as the location of the foil), build a mask matrix
    
    %% Determine gappy mask
    
    ug = u; % initialize the gappy fields
    vg = v; % initialize the gappy fields
    
    ug(corrl<0.5) = nan; % makes any values bellow the minimum correlation value == nan
    vg(corrl<0.5) = nan;
    
    spots = isnan(ug); % gap locations mask (missingness matrix)
    real_val = ~spots; % locations of real values
    
    num = sum(spots(:) == 1);
    G = (num/numel(spots))*100; % gappyness
    
    %% Preparing the complete dataset

%     tic;

    for i = 1:size(ug,3)
        ug_temp = reshape(ug(:,:,i),1,[]); ug_tall(:,i) = ug_temp; % reshape each snapshot into a column vector and store it as a column of the complete x dataset
        vg_temp = reshape(vg(:,:,i),1,[]); vg_tall(:,i) = vg_temp; % reshape each snapshot into a column vector and store it as a column of the complete y dataset

        uf_temp = reshape(u(:,:,i),1,[]); uf_tall(:,i) = uf_temp; % for the fully resolved flow field
        vf_temp = reshape(v(:,:,i),1,[]); vf_tall(:,i) = vf_temp;
    end

    xi_full = [uf_tall; vf_tall]; % original FULLY RESOLVED field, x and y components stacked, ONLY FOR VISUALIZATION

    xi = [ug_tall; vg_tall]; % original gappy field, x and y components stacked
    xi_0 = xi;

    gaps = isnan(xi); % gaps in the complete set (all spatial and temporal data)

    %% First iteration (temporal mean guess)

    % Indeces:
    % i --> space, rows
    % k --> time, columns

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

    iterations = length(compiled_data)-1; % maximum number of iterations (must be < maximum num of modes obtainable from the dataset)
    if iterations >= 50
        iterations = 50; % just to not make too massive the number of iterations
    end

    Err = zeros(size(xi,1),size(xi,2),iterations); % initialize the error matrix
    P = 1; % starting number of modes used (for the initial iteration use only two modes {P = P0 + 1})

    for main = 1:iterations

        [U, S, V] = svd(xi, 'econ'); % SVD to get spatial modes (U) and temporal eigenmodes (S*V')

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

        % f = NaN(P,size(xi,2)); % alternate way of performing this calculation
        %
        % for k = 1:size(xi,2)
        %     for i = 1:P
        %         f(i,k) = dot(xi(:,k).*gaps(:,k),phi_tild(:,i));
        %     end
        % end

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
                if gaps(i,k) == 0
                    Err(i,k,main) = abs(xi_1(i,k) - xi_hat(i,k));
                end
            end
        end

        % %% Overwrite old matrix with new interpolated matrix

        for i = 1:size(xi,1)
            for k = 1:size(xi,2)
                if gaps(i,k) == 1
                    xi(i,k) = xi_hat(i,k);
                end
            end
        end

    end
    
    %% Reconstruction of the Flow Field

    u_talla = xi(1:end/2,:);
    v_talla = xi(end/2+1:end,:);

    for i = 1:size(ug,3)
        v_temp_a = v_talla(:,i);
        vg_a(:,:,i) = reshape(v_temp_a,size(ug,1),size(ug,2)); % reconstructed v
        u_temp_a = u_talla(:,i);
        ug_a(:,:,i) = reshape(u_temp_a,size(vg,1),size(ug,2)); % reconstructed u
    end

    % Adding the mask accounting for the foil's location to the reconstructed flow

    ug_a = ug_a.*foil_mask;
    vg_a = vg_a.*foil_mask;

%     toc;
    
    %% Visualization

%     q = 8; % flow at a given instance
% 
%     figure()
%     % contourf(x(xwindow,ywindow,q), y(xwindow,ywindow,q), UU(xwindow,ywindow,q), 'LineStyle', 'none'); hold on;
%     quiver(x(:,:,q), y(:,:,q), u(:,:,q), v(:,:,q), 2, 'color', 'k')
%     % hold off;
%     axis image
%     axis off
%     title('Fully Resolved', 'Interpreter', 'latex', 'FontSize', 20);
% 
%     % UUg = sqrt(ug.^2 + vg.^2);
% 
%     figure()
%     % contourf(x(xwindow,ywindow,q), y(xwindow,ywindow,q), UUg(xwindow,ywindow,q), 'LineStyle', 'none'); hold on;
%     quiver(x(:,:,q), y(:,:,q), ug(:,:,q), vg(:,:,q), 2, 'color', 'k')
%     % hold off;
%     axis image
%     axis off
%     title('Gapped Field', 'Interpreter', 'latex', 'FontSize', 20);
% 
%     % UUg_a = sqrt(ug_a.^2 + vg_a.^2);
% 
%     figure()
%     % contourf(x(xwindow,ywindow,q), y(xwindow,ywindow,q), UUg_a(xwindow,ywindow,q), 'LineStyle', 'none'); hold on;
%     quiver(x(:,:,q), y(:,:,q), ug_a(:,:,q), vg_a(:,:,q), 2, 'color', 'k')
%     % hold off;
%     axis image
%     axis off
%     title('Reconstructed Field', 'Interpreter', 'latex', 'FontSize', 20);
    
    %% Calculate Z vorticity
    
    vortg_a = NaN(size(x)); % there will be leftover nan values in this matrix
    
%     for n = 1:length(compiled_data)
%         xn = squeeze(x(:,:,n)); yn = squeeze(y(:,:,n)); ugn = squeeze(ug_a(:,:,n)); vgn = squeeze(vg_a(:,:,n));
%         vortg_a(:,:,n) = curl(ugn, vgn)*100;
%     end
    
    for  n = 1:length(compiled_data)
        ugn = squeeze(ug_a(:,:,n)); vgn = squeeze(vg_a(:,:,n)); % this might be unnecessary
        for i = 2:(size(ugn,1)-1)
            for j = 2:(size(ugn,2)-1)
                vortg_a(i,j,n) = ((vgn(i+1,j) - vgn(i-1,j)) - (ugn(i,j+1) - ugn(i,j-1)))*100;
            end
        end
    end
    
%     TF = isnan(vortg_a);
%     vortg_a(TF) = 0;

    
    %% Save into output variable
    
    parfor j = 1:length(compiled_data)
        interpolated_data(j).xg = x(:,:,j);
        interpolated_data(j).yg = y(:,:,j);
        interpolated_data(j).ug = ug_a(:,:,j);
        interpolated_data(j).vg = vg_a(:,:,j);
        interpolated_data(j).vortg = vortg_a(:,:,j);
    end
    
end