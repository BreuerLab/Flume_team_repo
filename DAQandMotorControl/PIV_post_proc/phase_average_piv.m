%% Flowfield phase averaging
% X --> Data matrix (rows = coordinate, columns = snapshot)

% Updated 20230207

function X_avg = phase_average_piv(X, frames_per_cycle, total_frames)

if ~exist('total_frames','var')
    total_frames = size(X,2);
end

X_avg = nan(size(X,1),frames_per_cycle); % initialize the averaged data matrix
temporary = nan(size(X,1),total_frames/frames_per_cycle); % initialize the temporary matrix

for k = 1:frames_per_cycle
    j = 0;
    for i = k:frames_per_cycle:(total_frames-frames_per_cycle+k)
        j = j + 1;
        temporary(:,j) = X(:,i);
    end
    X_avg(:,k) = mean(temporary,2);
end

end