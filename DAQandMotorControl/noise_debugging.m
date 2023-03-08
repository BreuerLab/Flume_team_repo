% Wallace sensor

[pxx_Fxw, f_Fxw] = pwelch(out(:,7),[],[],[],samplerate);
[pxx_Fyw, f_Fyw] = pwelch(out(:,8),[],[],[],samplerate);
[pxx_Fzw, f_Fzw] = pwelch(out(:,9),[],[],[],samplerate);
[pxx_Mxw, f_Mxw] = pwelch(out(:,10),[],[],[],samplerate);
[pxx_Myw, f_Myw] = pwelch(out(:,11),[],[],[],samplerate);
[pxx_Mzw, f_Mzw] = pwelch(out(:,12),[],[],[],samplerate);

plot(f_Fxw, 10*log10(pxx_Fxw)); hold on;
plot(f_Fyw, 10*log10(pxx_Fyw));
plot(f_Fzw, 10*log10(pxx_Fzw));
plot(f_Mxw, 10*log10(pxx_Mxw));
plot(f_Myw, 10*log10(pxx_Myw));
plot(f_Mzw, 10*log10(pxx_Mzw)); hold off;
xlim([0,100])
legend('Fx','Fy','Fz','Mx','My','Mz')

title('Wallace')

% Gromit sensor

[pxx_Fxg, f_Fxg] = pwelch(out(:,17),[],[],[],samplerate);
[pxx_Fyg, f_Fyg] = pwelch(out(:,18),[],[],[],samplerate);
[pxx_Fzg, f_Fzg] = pwelch(out(:,19),[],[],[],samplerate);
[pxx_Mxg, f_Mxg] = pwelch(out(:,20),[],[],[],samplerate);
[pxx_Myg, f_Myg] = pwelch(out(:,21),[],[],[],samplerate);
[pxx_Mzg, f_Mzg] = pwelch(out(:,22),[],[],[],samplerate);

plot(f_Fxg, 10*log10(pxx_Fxg)); hold on;
plot(f_Fyg, 10*log10(pxx_Fyg));
plot(f_Fzg, 10*log10(pxx_Fzg));
plot(f_Mxg, 10*log10(pxx_Mxg));
plot(f_Myg, 10*log10(pxx_Myg));
plot(f_Mzg, 10*log10(pxx_Mzg)); hold off;
xlim([0,100])
legend('Fx','Fy','Fz','Mx','My','Mz')

title('Gromit')
