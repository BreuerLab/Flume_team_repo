function [out,bias,dat] = find_bias_3rigs(dq,last_out,flume_hertz,fname,foil)
bias_trialduration = 10;

write(dq,last_out)
fprintf('Checklist:\n  - Zero Flume Velocity\n')
fprintf('Press any key to continue\n\n')
pause
fprintf('Finding Bias ...\n')

flume_hertz_old = flume_hertz;
flume_hertz = 0;
hprof1=linspace(last_out(2),last_out(2),bias_trialduration*dq.Rate)';

pprof1=linspace(last_out(1),last_out(1),bias_trialduration*dq.Rate)';

hprof2=linspace(last_out(4),last_out(4),bias_trialduration*dq.Rate)';

pprof2=linspace(last_out(3),last_out(3),bias_trialduration*dq.Rate)';

hprof3=linspace(last_out(6),last_out(6),bias_trialduration*dq.Rate)';

pprof3=linspace(last_out(5),last_out(5),bias_trialduration*dq.Rate)';

trigger=linspace(last_out(7),last_out(7),bias_trialduration*dq.Rate)';

output = [pprof1 hprof1 pprof2 hprof2 pprof3 hprof3 trigger];  % needs additional trigger channel for PIV

% figure(3)
% plot(output)
% s.IsNotifyWhenDataAvailableExceedsAuto=true;
% prepare(s)

% dq.queueOutputData(output); 
% [dat,t] = dq.startForeground;
dat = readwrite(dq,output,"OutputFormat","Matrix");

% for ii = [5:10 15:20]
% dat(:,ii) = medfilt1(dat(:,ii),10);
% end
bias.Wallace = mean(dat(:,7:12),1);
bias.WallaceStdev = std(dat(:,7:12),1);
bias.Gromit = mean(dat(:,17:22),1);
bias.GromitStdev = std(dat(:,17:22),1);
bias.accmeter = mean(dat(:,23),1);
bias.accStdev = std(dat(:,23),1);

% Heave_voltage = mean(dat(:,9),1);
% -.05158-.05602*sin(2*pi*.2704)
out = output_conv_3rigs(dat,bias,foil);
figure(1)
subplot(2,1,2)
% plot(dat(20:end-20,5:10) - repmat(bias.Wallace,numel(out(20:end-20,3)),1),'.')
plot(out(20:end-20,7:9)./.125,'.')
hold on 
plot(out(20:end-20,10:12)*1333/10,'.')
hold off
title('Wallace (last)')
ylabel('Forces and Torques (normalized by resolution)')
subplot(2,1,1)
% plot(dat(20:end-20,15:20) - repmat(bias.Gromit,numel(out(20:end-20,3)),1),'.')
plot(out(20:end-20,17:19)./.125,'.')
hold on 
plot(out(20:end-20,20:22)*1333/10,'.')
hold off
title('Gromit (middle)')
ylabel('Forces and Torques (normalized by resolution)')


% changed first term from "out" to "dat"
bias.RMSEW = sqrt(mean((dat(:,7:12) - repmat(bias.Wallace,numel(out(:,3)),1)).^2));
if sum(bias.RMSEW>[.15 .15 .3 .1 .1 .1])>0 
        disp(bias.RMSEW);
    disp('Warning: Wallace error signal above normal. Check wiring/ grounding.')
end
bias.RMSEG = sqrt(mean((dat(:,17:22) - repmat(bias.Gromit,numel(out(:,3)),1)).^2));
if sum(bias.RMSEG>[.15 .15 .3 .1 .1 .1])>0 
    disp('Warning: Gromit error signal above normal. Check wiring/ grounding.')
end

    disp(['Accelerometer voltage: ',num2str(bias.accmeter)]);
if bias.accmeter < 1.5 || bias.accmeter > 1.8
    disp('Accelmeter voltage far from expected ~1.68V, power cycle and try again.');
end
    

Percent_fullrange_error = bias.RMSEW./[660 660 1980 60 60 60]*100;
Percent_fullrange_errorG = bias.RMSEG./[660 660 1980 60 60 60]*100;
time = clock;
d=datevec(date);

folder_name = [fname,'\data'];

numfiles = dir([folder_name,'\bias*']);
jj = numel(numfiles)+1;
filename=[folder_name,'\bias_',num2str(jj)];

save(filename,'bias','time','dat');
flume_hertz = flume_hertz_old;
end