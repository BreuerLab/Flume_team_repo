%% Phase averaging

function [u_avg, v_avg, vort_avg] = phase_average_piv(u, v, vort, frames_per_cycle, total_frames)

X = [u; v; vort]; % stack all variables on top of each other to do the phase averaging faster

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

u_avg = X_avg(1:end/3,:); % phase-averaged data
v_avg = X_avg(end/3+1:2*end/3,:); % phase-averaged data
vort_avg = X_avg(2*end/3+1:end,:); % phase-averaged data

end