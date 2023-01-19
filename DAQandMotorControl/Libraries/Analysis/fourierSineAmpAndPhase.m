function [fit_A1,fit_phi] = fourierSineAmpAndPhase(time_star,data,freq,num_cyc)

% T = 1/1000;
% t = (1:T:200)';
% freq = 2;
% phi = pi/2;
% input = 2*sin(2*pi*freq*t+phi);

period = 1/freq;
period_timesteps = round(1/(time_star(2)-time_star(1))); % number of timesteps per period of oscillation 
timesteps = length(time_star);

% Perform a fit of a fourier series using a sliding window with width equal to one period of prescribed oscillation

fit_A1 = zeros(num_cyc,1);
fit_phi = zeros(num_cyc,1);
windowstep = 1;

for timestep = period_timesteps/2+1:period_timesteps:timesteps-period_timesteps/2
time_star_windowed = time_star(timestep-period_timesteps/2:timestep+period_timesteps/2)-time_star(timestep-period_timesteps/2);
data_windowed = data(timestep-period_timesteps/2:timestep+period_timesteps/2);
% Fit a fourier series of the form a0 + a1*cos(x*w) + b1*sin(x*w)
fourier_fit = fit(time_star_windowed'/freq,data_windowed,'fourier1'); %,'StartPoint', [0 1 1 2*pi*freq]
% plot(data_windowed)
% hold on
% plot(fourier_fit,time_star_windowed'/freq,data_windowed); ylim([min(data)-1,max(data)+1]); xlim([0,1]);
% hold off
% pause(0.1)

% fit_mean = fourier_fit.a0;
fit_A1(windowstep) = sqrt(fourier_fit.a1^2+fourier_fit.b1^2);
fit_phi(windowstep) = atan2(fourier_fit.a1,fourier_fit.b1);
windowstep = windowstep+1;
end


end