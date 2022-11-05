function [last_out,bias] = find_zero_pitch_Gromit(dq,last_out,bias,foil)
% Finds the minimum lift on gromit hydrofoil while under flow.

% Updated 20221102 - Eric Handy

disp('finding zero.  Gromit will now move +/- 5 degrees')
[~,~,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 5 0 0 0],2,bias,foil); % added output to get "last_out" for the next move_new_pos_3rigs
 [out,prof,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 -5 0 0 0],10,bias,foil);



Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3))); % Lift and torque calculations are performed taking into account the conversion from the
Torque(:,1) = -out(:,22);                                        % force transducer's frame of reference to the establishde flume frame of reference

f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
f_tq = polyfit(smooth(Torque,100)',1:numel(Torque),1);

plot(1:numel(prof(:,1)),prof(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
% 1:numel(lift) = Lift*a + b
%Lift = (1:numel(Lift)-b)./a)
 [out,prof2,last_out] =  move_new_pos_3rigs(dq,last_out,[0 0 5 0 0 0],10,bias,foil);


Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3)));
Torque(:,1) = -out(:,22);

f1=polyfit(smooth(Lift,100)',1:numel(Lift),1);
f_tq1 = polyfit(smooth(Torque,100)',1:numel(Torque),1);

hold on
plot(1:numel(prof2(:,1)),prof2(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
hold off




if round(f(2)) < 0 || round(f(2)) > numel(prof(:,3))
    disp('Pitch out of range. Expanding search.')
 [~,~,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 15 0 0 0],2,bias,foil);
 [out,prof,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 -15 0 0 0],10,bias,foil);


    Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3)));

    f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    plot(1:numel(prof(:,1)),prof(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
    hold on
 [out,prof2,last_out] =  move_new_pos_3rigs(dq,last_out,[0 0 15 0 0 0],10,bias,foil);


    Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3)));

    f1=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    plot(1:numel(prof2(:,1)),prof2(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
    hold off
    
    
    
    bias.pitch(2) = mean([prof(round(f(2)),3) prof2(round(f1(2)),3)]);
    
    
    if round(f(2)) < 0 || round(f(2)) > numel(prof(:,3))
        error('Pitch out of range.  Manually align pitch.')
    end
    
    [~,~,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 5 0 0 0],2,bias,foil);
    [out,prof,last_out] = move_new_pos_3rigs(dq,last_out,[0 0 -5 0 0 0],10,bias,foil);


    Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3)));

    f=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    plot(1:numel(prof(:,1)),prof(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f(2))./f(1))
    % 1:numel(lift) = Lift*a + b
    %Lift = (1:numel(Lift)-b)./a)
    [out,prof2,last_out] =  move_new_pos_3rigs(dq,last_out,[0 0 5 0 0 0],10,bias,foil);


    Lift(:,1) = (out(:,17).*cos(out(:,3)) + out(:,18).*sin(out(:,3)));

    f1=polyfit(smooth(Lift,100)',1:numel(Lift),1);
    hold on
    plot(1:numel(prof2(:,1)),prof2(:,3),1:numel(out(:,18)),smooth(Lift,100),1:numel(Lift),((1:numel(Lift))-f1(2))./f1(1))
    hold off
    
    
end

bias.pitch(2) = mean([prof(round(f(2)),3) prof2(round(f1(2)),3)]);          % Based on zero lift
% pitch_bias(2) = mean([prof(round(f_tq(2)),3) prof2(round(f_tq1(2)),3)]);  % Based on zero torque

disp(['Pitch Bias (volts):  ',num2str(bias.pitch)])

[last_out] = move_to_zero(dq,last_out,bias);

% disp(['Pitch Bias (deg)',num2str(conv_last_out(last_out))])

end
