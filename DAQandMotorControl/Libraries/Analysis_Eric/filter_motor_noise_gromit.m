function out = filter_motor_noise_gromit(out, freq, srate, cutoff)

[b,a] = butter(6, freq*cutoff/(srate/2), 'low');

out(:,17) = filtfilt(b, a, out(:,17));
out(:,18) = filtfilt(b, a, out(:,18));
out(:,22) = filtfilt(b, a, out(:,22));

end