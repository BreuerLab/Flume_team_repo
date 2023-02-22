function [last_out, bias] = find_zero_pitch(dq, last_out, bias, foil, traverse)
% Finds the minimum lift on the desired hydrofoil while under flow in order
% to align its pitch angle correctly.
% traverse must be 1 (leading traverse) or 2 (trailing traverse)

switch traverse
    case 1 % leading traverse --> gromit sensor --> wallace position
        fy_index = 18; % index of the perpendicular (to the streamwise direction) transducer force channel
        mz_index = 22; % index of the pitch axis transducer moment channel
        position = [0, 0, 0, 0, 5, 0]; % position commanded to the pitch axis of the selected traverse
        out_of_the_way = [0, 0, 0, 0, 0, 0]; % for trailing alignment (zeros for the leading foil)
        prof_index = 5; % index of the motion profile corresponding to the selected traverse
        p_bias_index = 3; % pitch bias index for the corresponding traverse
    case 2 % trailing traverse --> wallace sensor --> gromit position
        fy_index = 8;
        mz_index = 12;
        position = [0, 0, 5, 0, 0, 0];
        out_of_the_way = [0, 0, 0, 0, 0, 0.15]; % to command the leading traverse to move out of the way during trailing alignment
        prof_index = 3;
        p_bias_index = 2;
end

disp('finding zero.  Traverse will now move +/- 5 degrees')
[~,~,last_out] = move_new_pos_3rigs(dq,last_out,position+out_of_the_way,10,bias,foil);
scan_time = 30;
pause(10);
b1_Vtheta = last_out(prof_index);
 [out,prof,last_out] =  move_new_pos_3rigs(dq,last_out,-position+out_of_the_way,scan_time,bias,foil);
a1_Vtheta = (last_out(prof_index)-b1_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
% Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));
Torque(:,1) = out(:,mz_index);

% f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
coefL1=polyfit(1:numel(Lift),smooth(Lift,100)',1);
aL1 = coefL1(1); bL1 = coefL1(2);
% f_tq = polyfit(smooth(Torque,100)',1:numel(Torque),1);

% plot(1:numel(prof(:,1)),prof(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
% 1:numel(lift) = Lift*a + b
%Lift = (1:numel(Lift)-b)./a)
plot(1:numel(prof(:,1)),prof(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL1*(1:numel(Lift))+bL1)

b2_Vtheta = last_out(prof_index);
 [out,prof2,last_out] =  move_new_pos_3rigs(dq,last_out,position+out_of_the_way,scan_time,bias,foil);
a2_Vtheta = (last_out(prof_index)-b2_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
% Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));
Torque(:,1) = out(:,mz_index);

coefL2=polyfit(1:numel(Lift),smooth(Lift,100)',1);
aL2 = coefL2(1); bL2 = coefL2(2);
% f_tq1 = polyfit(smooth(Torque,100)',1:numel(Torque),1);

hold on
% plot(1:numel(prof2(:,1)),prof2(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
plot(1:numel(prof2(:,1)),prof2(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL2*(1:numel(Lift))+bL2)
hold off
legend('Pitch negative','Lift smoothed','Linear fit','Pitch positive','Lift smoothed','Linear fit')



if max(aL1*(1:numel(Lift))+bL1) < 0 || min(aL1*(1:numel(Lift))+bL1) > 0
    disp('Pitch out of range. Expanding search.')
 [~,~,last_out] =   move_new_pos_3rigs(dq,last_out,3*position+out_of_the_way,5,bias,foil);
 b1_Vtheta = last_out(prof_index);
 [out,prof,last_out] =  move_new_pos_3rigs(dq,last_out,-3*position+out_of_the_way,scan_time,bias,foil);
 a1_Vtheta = (last_out(prof_index)-b1_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
%     Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));

%     f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    coefL1=polyfit(1:numel(Lift),smooth(Lift,100)',1);
    aL1 = coefL1(1); bL1 = coefL1(2);
%     plot(1:numel(prof(:,1)),prof(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
    plot(1:numel(prof(:,1)),prof(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL1*(1:numel(Lift))+bL1)
    hold on

    b2_Vtheta = last_out(prof_index);
 [out,prof2,last_out] =  move_new_pos_3rigs(dq, last_out, 3*position+out_of_the_way, scan_time, bias, foil);
 a2_Vtheta = (last_out(prof_index)-b2_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
%     Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));

%     f1=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    coefL2=polyfit(1:numel(Lift),smooth(Lift,100)',1);
    aL2 = coefL2(1); bL2 = coefL2(2);
%     bias.pitch(3) = mean([prof(round(f(2)),5) prof2(round(f1(2)),5)]);
    pitch_voltbias1 = a1_Vtheta*(-bL1/aL1)+b1_Vtheta;
    pitch_voltbias2 = a2_Vtheta*(-bL2/aL2)+b2_Vtheta;
    bias.pitch(p_bias_index) = mean([pitch_voltbias1 pitch_voltbias2]);

    %     plot(1:numel(prof2(:,1)),prof2(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
    plot(1:numel(prof(:,1)),prof(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL2*(1:numel(Lift))+bL2)
    hold off
    
    
    if max(aL1*(1:numel(Lift))+bL1) < 0 || min(aL1*(1:numel(Lift))+bL1) > 0
%         error('Pitch out of range.  Manually align pitch.')
        disp('Pitch out of range.  Manually align pitch.')
    end
    
    [~,~,last_out] = move_new_pos_3rigs(dq, last_out, position+out_of_the_way, 5, bias, foil);
     b1_Vtheta = last_out(prof_index);
    [out,prof,last_out] =  move_new_pos_3rigs(dq, last_out, -position+out_of_the_way, scan_time, bias, foil);
     a1_Vtheta = (last_out(prof_index)-b1_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
%     Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));

%     f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    coefL1=polyfit(1:numel(Lift),smooth(Lift,100)',1);
    aL1 = coefL1(1); bL1 = coefL1(2);
%     plot(1:numel(prof(:,1)),prof(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
    plot(1:numel(prof(:,1)),prof(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL1*(1:numel(Lift))+bL1)
    % 1:numel(lift) = Lift*a + b
    %Lift = (1:numel(Lift)-b)./a)
    b2_Vtheta = last_out(prof_index);
    [out,prof2,last_out] =  move_new_pos_3rigs(dq,last_out,position+out_of_the_way,scan_time,bias,foil);
    a2_Vtheta = (last_out(prof_index)-b2_Vtheta)/(scan_time*dq.Rate);

Lift(:,1) = out(:,fy_index);
%     Lift(:,1) = (out(:,7).*cos(out(:,5))+out(:,8).*sin(out(:,5)));

%     f1=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    coefL2=polyfit(1:numel(Lift),smooth(Lift,100)',1);
    aL2 = coefL2(1); bL2 = coefL2(2);
    hold on
%     plot(1:numel(prof2(:,1)),prof2(:,5),1:numel(out(:,8)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
    plot(1:numel(prof(:,1)),prof(:,prof_index),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),aL2*(1:numel(Lift))+bL2)
    hold off
    
    
end

% bias.pitch(3) = mean([prof(round(f(2)),5) prof2(round(f1(2)),5)]);          % Based on zero lift
% pitch_bias(3) = mean([prof(round(f_tq(2)),5) prof2(round(f_tq1(2)),5)]);  % Based on zero torque
    pitch_voltbias1 = a1_Vtheta*(-bL1/aL1)+b1_Vtheta;
    pitch_voltbias2 = a2_Vtheta*(-bL2/aL2)+b2_Vtheta;
    bias.pitch(p_bias_index) = mean([pitch_voltbias1 pitch_voltbias2]);

disp(['Pitch Bias (volts):  ',num2str(bias.pitch)])

pitch_check = input(['Does this look alright to you? y/n',newline],"s");
if pitch_check == 'y'

    % comment out to keep motor from moving
    [~,~,last_out] =  move_new_pos_3rigs(dq,last_out,[0 0 0 0 0 0],5,bias,foil); % also returns leading traverse to zero
    % [last_out] = move_to_zero(dq,last_out,bias);
end

% disp(['Pitch Bias (deg)',num2str(conv_last_out(last_out))])

end