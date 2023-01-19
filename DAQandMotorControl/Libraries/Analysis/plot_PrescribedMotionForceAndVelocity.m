function plot_PrescribedMotionForceAndVelocity(time_star,heave_star_measured,heave_velo,liftcoef,dragcoef,powercoef,num_cyc,dragtorquecoef)
figure
hold on
grid on

% plot(time_ND,100*heave_velo,'DisplayName','Heave velocity','Color','green');
h(1) = plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','blue','LineWidth',4);
% plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(time_star,dragtorquecoef,'DisplayName','{\it C}_{\tau D}','Color','green','LineWidth',4)
h(2) = plot(time_star,powercoef,'DisplayName','{\it C}_P','Color','red','LineWidth',4);
h(3) = plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4);

hold off

legend(h([3 1 2]))
ylim([-5 5])
xlim([10 40]) %xlim([10 12]) %
xlabel('Time (cycles)')
% ylabel('Heave (cm), Force (N), Power (W)')
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
set(gcf, 'Position',  [100, 100, 1200, 800])
% exportgraphics(gcf,'LinePlot.png')
end


