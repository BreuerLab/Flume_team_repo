function out = filter_motor_noise_wallace(out, freq, srate, cutoff)

[b,a] = butter(6, freq*cutoff/(srate/2), 'low');

out(:,7) = filtfilt(b, a, out(:,7));
out(:,8) = filtfilt(b, a, out(:,8));
out(:,12) = filtfilt(b, a, out(:,12));

end