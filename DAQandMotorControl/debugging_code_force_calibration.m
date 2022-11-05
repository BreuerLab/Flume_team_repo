%% Testing forces calibration
% Plots mean values and std of every channel in the force sensors

clear;

fred = [0.08,0.09,0.1,0.11,0.12,0.13,0.14];
fred = 0.12;
U = 0.33;
c = 0.061;
span = 6*c;
f = U*fred/c;


figure(); hold on;

for ii = 1:length(fred)

%     foldername = 'R:\ENG_Breuer_Shared\group\Eric\Enter descriptive name_03-Nov-2022_15_58_1\data\';
%     filename = ['Vib_pitch=70deg,f=', num2str(f(ii),3), 'Hz,A=4.88cm.mat'];
%     filename = ['20221102_BiasTesting_A3Efred=',num2str(fred(ii)), '_p3=70_h3=0.8c_U=0.33'];
%     load(fullfile(foldername,filename));


    load('R:\ENG_Breuer_Shared\group\Eric\GromitTest_04-Nov-2022_18_2_30\data\CircCyl_pitch=0deg,f=0.649Hz,A=0cm.mat');
%     cycle = length(out)/36;
%     range = round(cycle*5):length(out)-round(cycle*5);
%     plot(fred(ii),mean(out(range,7)),'ko');
    c= 0.061;
    span = 6*c;
    cycle_period = round(length(out)/60);
    raaange = (5*cycle_period):1:(length(out)-5*cycle_period);
%     raaange = 1:length(out);

    heave_commanded = Prof_out_angle(raaange,4);
    heave_measured = out(raaange,4);
    % heave_star_measured = heave_measured/chord;
    pitch_measured = out(raaange,3);
    moment_z0 = -out(raaange,22); % negative due to the change in frame of reference
    force_x0 = out(raaange,17);
    force_y0 = out(raaange,18);
    force_D = force_y0.*cos(pitch_measured) - force_x0.*sin(pitch_measured);
    force_L = force_x0.*cos(pitch_measured) + force_y0.*sin(pitch_measured);
    inertialload_y = out(raaange,23);
    flowspeed_measured = mean(abs(out(raaange,13)));
    
    T = 1/1000;
    
    heave_velo = movmean((1/T)*gradient(squeeze(heave_measured)),100);
    heave_accel = movmean((1/T)*gradient(squeeze(heave_velo)),100);
    
    pitch_velo = movmean((1/T)*gradient(squeeze(pitch_measured)),100);
    pitch_accel = movmean((1/T)*gradient(squeeze(pitch_velo)),100);
    
    P_p = force_L.*heave_velo;
    P_h = moment_z0.*pitch_velo;
    
    yp1 = heave_measured + chord*0.5*sin(pitch_measured);
    yp2 = heave_measured - chord*0.5*sin(pitch_measured);
    Yp = 2*max(max(yp1),max(yp2));
    
    P_flow = 0.5*1000*flowspeed_measured^3*Yp*span;
    
    Eff = mean(P_p + P_h)/P_flow
    L_max = max(abs(force_L))

%     plot(fred(ii), Eff_tr,'ko');
end

%%

% figure()
% 
% w_forces_avg = [mean(out(:,7)),mean(out(:,8)),mean(out(:,9))];
% w_torque_avg = [mean(out(:,10)),mean(out(:,11)),mean(out(:,12))];
% 
% subplot(2,2,1)
% hold on;
% yyaxis left
% yline(w_forces_avg(1),'r','LineWidth',2)
% yline(w_forces_avg(2),'b','LineWidth',2)
% yline(w_forces_avg(3),'k','LineWidth',2)
% % ylim([1.5*abs(min(w_forces_avg)),1.5*abs(max(w_forces_avg))]);
% yyaxis right
% yline(w_torque_avg(1),'--r','LineWidth',2)
% yline(w_torque_avg(2),'--b','LineWidth',2)
% yline(w_torque_avg(3),'--k','LineWidth',2)
% % ylim([-1.5*abs(min(w_torque_avg)),1.5*abs(max(w_torque_avg))]);
% hold off;
% 
% legend('Fn','Ft','Fz','Mx','My','Mz')
% title('Wallace avg')
% 
% g_forces_avg = [mean(out(:,17)),mean(out(:,18)),mean(out(:,19))];
% g_torque_avg = [mean(out(:,20)),mean(out(:,21)),mean(out(:,22))];
% 
% subplot(2,2,2)
% hold on;
% yyaxis left
% yline(g_forces_avg(1),'r','LineWidth',2)
% yline(g_forces_avg(2),'b','LineWidth',2)
% yline(g_forces_avg(3),'k','LineWidth',2)
% % ylim([-1.5*abs(min(g_forces_avg)),1.5*abs(max(g_forces_avg))]);
% yyaxis right
% yline(g_torque_avg(1),'--r','LineWidth',2)
% yline(g_torque_avg(2),'--b','LineWidth',2)
% yline(g_torque_avg(3),'--k','LineWidth',2)
% % ylim([-1.5*abs(min(g_torque_avg)),1.5*abs(max(g_torque_avg))]);
% hold off;
% 
% legend('Fn','Ft','Fz','Mx','My','Mz')
% title('Gromit avg')
% 
% w_forces_std = [std(out(:,7)),std(out(:,8)),std(out(:,9))];
% w_torque_std = [std(out(:,10)),std(out(:,11)),std(out(:,12))];
% 
% subplot(2,2,3)
% hold on;
% yyaxis left
% yline(w_forces_std(1),'r','LineWidth',2)
% yline(w_forces_std(2),'b','LineWidth',2)
% yline(w_forces_std(3),'k','LineWidth',2)
% ylim([0,1.5*abs(max(w_forces_std))]);
% yyaxis right
% yline(w_torque_std(1),'--r','LineWidth',2)
% yline(w_torque_std(2),'--b','LineWidth',2)
% yline(w_torque_std(3),'--k','LineWidth',2)
% ylim([0,1.5*abs(max(w_torque_std))]);
% hold off;
% 
% legend('Fn','Ft','Fz','Mx','My','Mz')
% title('Wallace std')
% 
% g_forces_std = [std(out(:,17)),std(out(:,18)),std(out(:,19))];
% g_torque_std = [std(out(:,20)),std(out(:,21)),std(out(:,22))];
% 
% subplot(2,2,4)
% hold on;
% yyaxis left
% yline(g_forces_std(1),'r','LineWidth',2)
% yline(g_forces_std(2),'b','LineWidth',2)
% yline(g_forces_std(3),'k','LineWidth',2)
% ylim([0,1.5*abs(max(g_forces_std))]);
% yyaxis right
% yline(g_torque_std(1),'--r','LineWidth',2)
% yline(g_torque_std(2),'--b','LineWidth',2)
% yline(g_torque_std(3),'--k','LineWidth',2)
% ylim([0,1.5*abs(max(g_torque_std))]);
% hold off;
% 
% legend('Fn','Ft','Fz','Mx','My','Mz')
% title('Gromit std')
% 
