function plot_PrescribedMotionForceAndVelocity(time_star,heave_star_measured,heave_velo,liftcoef,dragcoef,powercoef,num_cyc,dragtorquecoef)
figure
hold on
grid on
plot(time_star,heave_star_measured,'DisplayName','{\it y/d}','Color','black','LineWidth',4)
% plot(time_ND,100*heave_velo,'DisplayName','Heave velocity','Color','green');
plot(time_star,liftcoef,'DisplayName','{\it C}_L','Color','red','LineWidth',4)
% plot(time_star,dragcoef,'DisplayName','{\it C}_D','Color','blue','LineWidth',4)
% plot(time_star,dragtorquecoef,'DisplayName','{\it C}_{\tau D}','Color','green','LineWidth',4)
plot(time_star,powercoef,'DisplayName','{\it C}_P','Color','blue','LineWidth',4)
hold off

legend()
ylim([-10 10])
xlim([30 35]) %xlim([10 15]) %
xlabel('Time (cycles)')
% ylabel('Heave (cm), Force (N), Power (W)')
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 
end
