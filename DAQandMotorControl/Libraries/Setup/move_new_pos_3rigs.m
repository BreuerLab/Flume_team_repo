function [out,output_prof,last_out] = move_new_pos_3rigs(dq,last_out,position,t,Wbias,Gbias,accbias,pitch_bias)
% n_pos is in volts, transitions position  from current to desired

% out = input_conv2(n_pos);

if isempty(dq)
    setup_DAQ;
    disp('checking daq');
    find_bias_3rigs;
end
n_pos = input_conv_3rigs(position,0,0,0,0,pitch_bias);
% dq.IsContinuous = false;
% load('C:\Users\Control Systems\Documents\vert_foil\last_out')
% load('C:\Users\Control Systems\Documents\vert_foil\last_out','-ascii')
 


pprof1=linspace(last_out(1),n_pos(1),dq.Rate*t)';
hprof1=linspace(last_out(2),n_pos(2),dq.Rate*t)';
pprof2=linspace(last_out(3),n_pos(3),dq.Rate*t)';
hprof2=linspace(last_out(4),n_pos(4),dq.Rate*t)';
pprof3=linspace(last_out(5),n_pos(5),dq.Rate*t)';
hprof3=linspace(last_out(6),n_pos(6),dq.Rate*t)';

output_prof = [pprof1 hprof1 pprof2 hprof2 pprof3 hprof3];

% dq.IsNotifyWhenDataAvailableExceedsAuto=true;
% dq.queueOutputData(output_prof);
% dat = dq.startForeground;
dat = readwrite(dq,output_prof,"OutputFormat","Matrix");

last_out=[pprof1(end) hprof1(end) pprof2(end),hprof2(end) pprof3(end),hprof3(end)];

[out,t]=output_conv_3rigs(dat,Wbias,Gbias,accbias);


end