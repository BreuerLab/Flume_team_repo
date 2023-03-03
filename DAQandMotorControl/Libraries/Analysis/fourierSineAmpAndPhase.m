function [fit_A1,fit_phi] = fourierSineAmpAndPhase(time_star,data)

period_timesteps = round(1/(time_star(2)-time_star(1))); % number of timesteps per period of oscillation 
window_duration =  period_timesteps; % Number of timesteps contained in one window
window_step =  floor(period_timesteps/10); % number of timesteps to move the window by (must be a fraction that evenly divides each period)
timesteps = length(time_star); % Total number of timesteps
window_steps = floor((timesteps-2*(window_duration-1))/window_step); % Number of windows that are fit

% Preallocate memory for fit parameters
fit_A1 = zeros(window_steps,1);
fit_phi = zeros(window_steps,1);
window = 1;

% Check if the window is an even number of timesteps 
if mod(window_duration,2) == 0
    half_window = window_duration/2;
else 
    half_window = (window_duration-1)/2;
end

% Fit a fourier series to the data for each data window
for timestep = half_window+1:window_step:timesteps-half_window

window_starttime = timestep-half_window;
window_endtime = timestep+half_window;
time_star_windowed = time_star(window_starttime:window_endtime);
data_windowed = data(window_starttime:window_endtime);

% % Fit a fourier series of the form a1*sin(2*pi*(f*x+p1))
% fo = fitoptions('Method','NonlinearLeastSquares',... % fitting method
%                'Lower',[0,-0.5],... % lower bounds for fit parameters a1 and b1
%                'Upper',[Inf,0.5],... 
%                'StartPoint', [0,0]);
% ft = fittype('a1*sin(2*pi*(f*x+p1))','problem','f','options',fo);
% Fit a fourier series of the form a1*sin(2*pi*f*x)+b1*cos(2*pi*f*x)
fo = fitoptions('Method','NonlinearLeastSquares',... % fitting method
               'Lower',[-Inf,-Inf],... % lower bounds for fit parameters a1 and b1
               'Upper',[Inf,Inf],... 
               'StartPoint', [0,0]);
ft = fittype('a1*sin(2*pi*f*x)+b1*cos(2*pi*f*x)','problem','f','options',fo);
fourier_fit = fit(time_star_windowed',data_windowed,ft,'problem',1);
% 
% % plot the fit function against the window of data
% close all
% hold on
% plot(fourier_fit,time_star_windowed',data_windowed); ylim([min(data)-1,max(data)+1]);
% hold off
% pause(1)

% Convert fit parameters into magnitude and phase form
fit_A1(window) = sqrt(fourier_fit.a1^2+fourier_fit.b1^2);%fourier_fit.a1;%
fit_phi(window) = atan2(fourier_fit.b1,fourier_fit.a1)/(2*pi); %fourier_fit.p1;%
window = window+1;
end


end