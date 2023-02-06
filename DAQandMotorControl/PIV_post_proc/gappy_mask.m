%% Gappy mask
% Creates a matrix with 1 denoting valid vectors and NaN denoting missing vectors based on a given reference
% from the PIV processing or manually placed gaps.

% Returns a logical vector matrix with the same size as input X with NaN in place of undesired vectors.

% Syntax:
%
%>> X = gappy_mask(X, 'corr', correlation_values, threshold)   Places gaps in locations where vectors have
%                                                               correlation values lower than a given threshold.
%
% >> X = gappy_mask(X, 'unc', correlation_values, threshold)    Places gaps in locations where vectors have high
%                                                               uncertainty percentages relative to the given
%                                                               threshold.
%
% >> X = gappy_mask(X, 'manual', fpc, matrix_size, component)   Here fpc is the number of frames per cycle, size
%                                                               is the dimensions of the reshaped snapshots [n,m],
%                                                               and component denotes the number of stacked
%                                                               components in X (1, 2 or 3).

function mask = gappy_mask(X, mask_type, reference, threshold, component)
    
    Xg = X; % matrix that will be masked
    mask = ones(size(Xg)); % initialize logical mask matrix
    
    %% Low-correlation values
    
    if strcmp(mask_type,'corr')
        Xg(reference<threshold) = NaN;
    end
    
    %% High uncertainty values
    
    if strcmp(mask_type,'unc')
        Xg(reference>threshold) = NaN;
    end
    
    %% Manual masking
    % Creates a manually masked matrix by placing individual masks in each frame of the cycle
    % Basically this is a cycle-averaged mask
    
    if strcmp(mask_type,'manual')
        
        fpc = reference; % number of snapshots
        n = threshold(1); % row dim of flowfield matrix
        m = threshold(2); % column dim of flowfield
        mask_range = [-3,-2,-1,0,1,2,3]'; % this determines the size of each "gap" of the mask
        x_temp = 1:n;
        y_temp = 1:m;
        [y,x] = meshgrid(y_temp,x_temp);
        mask_vert = ones(size(X(1:end/component,:))); % initialize the mask matrix, has dimensions of [length,time] for only ONE component
        
        for k = 1:fpc % loops over each snapshot
            Xtemp = X(1:end/component,:); % extracts only one component to work with
            snapshot = reshape(Xtemp(:,k),[n,m]); % reshapes the snapshot vector into a snapshot matrix
            
            figname = 'Manual masking';
            figure('Name',figname) % starts the manual process of generating a mask
            lvlstp = std(std(snapshot))*0.1; % interpolates colors between datapoints based on the data standard deviation
            contourf(x,y,snapshot,'linestyle','none','levelstep',lvlstp) % generates a contourplot of the snapshot pre-masking
            tit = ['Frame ',num2str(k),' of ',num2str(fpc)];
            title(tit);
            colormap('bone')
            fprintf('\nSelect locations to create mask\n\n')
            
            prompt = 'n'; % initialize prompt variable
            j = 0; % initialize "gap" counter (will increase with the number of gaps drawn)
            snap_temp = snapshot; % temporary snapshot matrix to plot placed gaps
            while prompt == 'n' % repeats if user chooses to NOT advance to the next snapshot
                j = j+1; % increases each time a new gap is added
                point = drawpoint; % user places a point where a gap wants to be placed
                range_x(:,j) = round(mask_range + point.Position(1)); % stores indeces of the gap (matrix cannot be initialized
                range_y(:,j) = round(mask_range + point.Position(2)); % because number of gaps per frame is unknown
                snap_temp(range_x(:,j),range_y(:,j)) = nan; % replaces placed gap location with NaN values
                contourf(x,y,snap_temp,'linestyle','none','levelstep',lvlstp) % plots result
                tit = ['Frame ',num2str(k),' of ',num2str(fpc)];
                title(tit);
                colormap('bone')
                prompt = input(['Next frame? [y/n]',newline],"s"); % asks user if they want to advance to the next frame
                                                                   % or place another gap in the current frame
            end
            close(figname) % closes figure before next snapshot loop
            
            mask_matrix = ones(size(snapshot)); % initializes a ones matrix the same size as the current snapshot matrix
            for j = 1:size(range_x,2) % loops over each gap placed and replaces the corresponding value in mask_matrix with a NaN value
                mask_matrix(range_x(:,j),range_y(:,j)) = nan;
            end
            mask_vert(:,k) = reshape(mask_matrix,size(Xtemp(:,k))); % reshapes into a vertical vector mask matrix
            
            switch component
            % replace the values in the original matrix (only the current frame) with NaN values based on the vertical mask matrix
                case 1 % for one component (vorticity for example)
                    Xg(1:end,k) = Xg(1:end,k).*mask_vert(:,k); % multiplies the real values with essentially a Valid-NaN vector
                case 2 % two components (x and y for example)
                    Xg(1:end/component,k) = Xg(1:end/component,k).*mask_vert(:,k);
                    Xg(end/component+1:end,k) = Xg(end/component+1:end,k).*mask_vert(:,k);
                case 3 % three components (x, y, and z for example)
                    Xg(1:end/component,k) = Xg(1:end/component,k).*mask_vert(:,k);
                    Xg(end/component+1:2*end/component,k) = Xg(end/component+1:2*end/component,k).*mask_vert(:,k);
                    Xg(2*end/component+1:end,k) = Xg(2*end/component+1:end,k).*mask_vert(:,k);
            end
            
        end
    end
    
    mask(isnan(Xg)) = nan; % mark locations of missing vectors as with NaN
    
end
