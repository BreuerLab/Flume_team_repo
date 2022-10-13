function [last_out] = move_to_zero(dq,last_out,pitch_bias)
% n_pos is in volts, transitions position  from current to desired

% out = input_conv2(n_pos);
if isempty(dq)
    setup_DAQ;
    disp('checking daq')
    find_bias_3rigs
end

% dq.IsContinuous = false;
% load('C:\Users\Control Systems\Documents\vert_foil\last_out')
% load('C:\Users\Control Systems\Documents\vert_foil\last_out','-ascii')


pprof1=linspace(last_out(1),pitch_bias(1),dq.Rate*2)';
hprof1=linspace(last_out(2),0,dq.Rate*2)';
pprof2=linspace(last_out(3),pitch_bias(2),dq.Rate*2)';
hprof2=linspace(last_out(4),0,dq.Rate*2)';
pprof3=linspace(last_out(5),pitch_bias(3),dq.Rate*2)';
hprof3=linspace(last_out(6),0,dq.Rate*2)';
trigger=linspace(last_out(7),0,dq.Rate*2)';

output_prof = [pprof1 hprof1 pprof2 hprof2 pprof3 hprof3 trigger];

% dq.IsNotifyWhenDataAvailableExceedsAuto=true;
% dq.queueOutputData(output_prof);
% dat = dq.startForeground;
dat = readwrite(dq,output_prof,"OutputFormat","Matrix");

last_out=[pprof1(end) hprof1(end) pprof2(end) hprof2(end) pprof3(end) hprof3(end) trigger(end)];


end