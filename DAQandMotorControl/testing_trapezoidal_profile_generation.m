beta = -1:0.2:2.4;

figure; hold on;

for i = 1:length(beta)

    [Time, Waveform] = generate_profile_trapezoidal(2, 0.5, 1000, 2, 2, 1, 0,0,beta(i));
    plot(Time, Waveform);

end