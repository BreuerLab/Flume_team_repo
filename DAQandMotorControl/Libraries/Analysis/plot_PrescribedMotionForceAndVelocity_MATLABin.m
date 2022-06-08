function plot_PrescribedMotionForceAndVelocity(time_star,heave_star_measured,heave_velo,liftcoef,powercoef)

hold on
grid on
plot(time_star,heave_star_measured,'DisplayName','h/c','Color','blue')
% plot(time_ND,100*heave_velo,'DisplayName','Heave velocity','Color','green');
plot(time_star,liftcoef,'DisplayName','C_L','Color','red')
% plot(time_star,powercoef,'DisplayName','C_P','Color','magenta')
hold off

legend()
ylim([-5 5])
xlabel('Time (cycles)')
% ylabel('Heave (cm), Force (N), Power (W)')
end
