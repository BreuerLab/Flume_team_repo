function plot_ForceFourierAnalysis(phi_CL1,CL1,f_star_commanded,A_star_measured)
figure
hold on
grid on

%% Plot fourier fit parameters
    yyaxis left
    windowsteps = (1:1:length(phi_CL1));
    h(1) = plot(windowsteps/10,phi_CL1,'LineStyle','none','Marker','.'); 
    ylim([-0.5 0.5]); 
%     title(['f*= ' num2str(f_star_commanded,2) ' A* = ' ...
%         num2str(A_star_measured,3)],'LineStyle','none');
    xlabel('Window steps')
    ylabel('Phase lag of heave \newline behind force (cycles)')
%     pause(2)
%     figure(2)
    yyaxis right
    h(2) = plot(windowsteps/10,CL1); 
    ylim([-10 10]); 
%     title(['f*= ' num2str(f_star_commanded,2) ' A* = ' num2str(A_star_measured,3)]);
    xlabel('Time (cycles)')
    ylabel('Amplitude of lift coef mode \newline with the prescribed frequency')    
    % Replace power coefficient with force in phase with velocity, CL1*sin(phi)

% % % plot(time_ND,100*heave_velo,'DisplayName','Heave velocity','Color','green');
% % h(1) = plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','blue','LineWidth',4);
% % % plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% % % plot(time_star,dragtorquecoef,'DisplayName','{\it C}_{\tau D}','Color','green','LineWidth',4)
% % h(2) = plot(time_star,powercoef,'DisplayName','{\it C}_P','Color','red','LineWidth',4);
% % h(3) = plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4);

hold off

% legend
% ylim([-5 5])
% xlim([10 40]) %xlim([10 12]) %
% xlabel('Time (cycles)')
% % ylabel('Heave (cm), Force (N), Power (W)')
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
set(gcf, 'Position',  [100, 100, 1200, 800])
% exportgraphics(gcf,'LinePlot.png')
end