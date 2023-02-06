%% Gappy Proper Orthogonal Decomposition
% Adapted mainly from the following references:
% > Venturi and Karniadakis, 2004
% > Bui-Thanh, Damodaran and Willcox, 2004
% > Someone else

% Updated 20230123

% X : Gapped flowfield (gaps replaced with NaN values)
% mask : (optional) Mask for a foil or other object

function [out, details] = GPOD(xi)

    tStart = tic;
    
    gaps = isnan(xi); % gaps locations logical matrix
    num = sum(gaps(:) == 1); % total number of gaps in the data
    details.gappyness = (num/numel(gaps))*100; % ratio num_gaps/num_total
    
    xi_1 = zeros(size(xi)); % gappy dataset, with zeros instead of NaN values
    xi_1(gaps==0) = xi(gaps==0);
    
    %% Temporal mean:
    % i --> space, rows
    % k --> time, columns
    
    for k = 1:size(xi,2) % loop over each column (time)
        for i = 1:size(xi,1) % loop over each row (space)
            if isnan(xi(i,k)) % if a gap is encountered
                index = ~isnan(xi(i,:)); % find all the known (not gaps) temporal data in the selected location
                if index == zeros(1,size(xi,2))
                    temp = xi(gaps(:,k)==0);
                    xi(i,k) = mean(temp); % if there are no real values in the temporal data, make value zero CHECK THIS!!!
                else
                    xi(i,k) = mean(xi(i,index)); % replace the gap with the temporal mean at that spatial location
                end
            end
        end
    end
    
%     for i = 1:size(xi,1)
%         for k = 1:size(xi,2)
%             if gaps(i,k) == 0
%                 xi_1(i,k) = xi(i,k); % writes only real values into matrix
%             end
%         end
%     end
    
    %% Performing POD

    iterations = size(xi,2)-2; % maximum number of iterations (must be < maximum num of modes obtainable from the dataset)
    if iterations >= 100
        iterations = 100; % just to not make too massive the number of iterations
    end
    
    iteration_err = zeros(size(xi)); % in-loop error matrix
    details.Err = zeros(size(xi,1),size(xi,2),iterations); % initialize the error matrix
    P = 2; % starting number of modes used (for the initial iteration use only two modes {P = P0 + 1})

    for main = 1:iterations

        [phi, ~, ~] = svd(xi, 'econ'); % SVD to get <spatial modes U = phi> and <temporal eigenmodes S*V'> (the latter unused)

        % % Truncation
        % From Venturi and Karniadakis:
        P = P + 1; % increase the admitted spatial modes every iteration

        % % Reconstuction error minimization via Least Squares
        % Since we only want to consider the error between the real values and the approximated values, we'll make all
        % interpolated values zero for the error calculation.

        phi_tild = phi(:,1:P); % truncated spatial modes matrix

        % Reconstruction matrix (kernel?) and projected field vector:
        K = NaN(P,P,size(xi,2)); % reconstruction matrix ([K]ij = (phi_i,phi_j)k --> remember k is the time index

        for k = 1:size(xi,2)
            for i = 1:P
                for j = 1:P
                    K(i,j,k) = dot(phi_tild(:,i).*~gaps(:,k),phi_tild(:,j));
                end
            end
        end

        f = phi_tild\(xi_1.*~gaps); %projected field vector phi * a = u_tild(masked)
        
        % f = NaN(P,size(xi,2)); % alternate way of performing this calculation
        %
        % for k = 1:size(xi,2)
        %     for i = 1:P
        %         f(i,k) = dot(xi(:,k).*gaps(:,k),phi_tild(:,i));
        %     end
        % end
        
        a = NaN(P,size(xi,2));
        
        for k = 1:size(xi,2)
            a(:,k) = K(:,:,k)\f(:,k);
        end

        % Reconstruct the field
        xi_hat = phi_tild*a; % reconstructed field using the new temporal coefficients

%         % Reconstruction Error calculation (old code)
%         for i = 1:size(xi,1)
%             for k = 1:size(xi,2)
%                 if gaps(i,k) == 0
%                     Err(i,k,main) = abs(xi_1(i,k) - xi_hat(i,k));
%                 end
%             end
%         end

        % Reconstruction error
        iteration_err(gaps==0) = abs(xi_1(gaps==0) - xi_hat(gaps==0));
        details.Err(:,:,main) = iteration_err;
        
        
%         % Overwrite old matrix with new interpolated matrix (old code)
%         for i = 1:size(xi,1)
%             for k = 1:size(xi,2)
%                 if gaps(i,k) == 1
%                     xi(i,k) = xi_hat(i,k);
%                 end
%             end
%         end
        
        xi(gaps==1) = xi_hat(gaps==1); % overwrite old matrix with interpolated data


    end
    
    out = xi; % write to output
    
    details.tEnd = toc(tStart);
    disp(['GPOD calculation took ', num2str(details.tEnd/60/60), ' hours.'])
    
end