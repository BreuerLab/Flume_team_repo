close all;
plot_w_fstar = 0;
plot_w_Ustar = 1;

% Sort the trials by frequency (needed for countorf function)
[f_star_sorted,sort_index] = sortrows(f_star_commanded);
U_star_sorted = 1./f_star_sorted;
phase12_sorted = phase12(sort_index,:);
A_star_sorted = A_star_measured(sort_index,:);
powercoef_mean_sorted = powercoef_mean(sort_index,:);
movmeanpower_points = 1;
powercoef_mean_sorted_smoothed = smooth2a(powercoef_mean_sorted,movmeanpower_points,movmeanpower_points);
% power_mean_sorted = power_mean(sort_index,:);
delay_sorted = delay(sort_index,:);
% powercoef_conv = powercoef_convtest(1,sort_index,:);

if plot_w_fstar == 1
    independent_var = f_star_sorted;
    xlimits = [0.09 0.31];
    xlabelstr = '{\it f} * = {\it f D/U }';
elseif plot_w_Ustar == 1
    independent_var = 1./f_star_sorted;
    xlimits = [0 18];
    xlabelstr = '{\it U} * = {\it U/Df }';
end

hold on

% Plot Cp vs. A* and f*
% contourf(f_star_sorted,A_star_sorted,power_mean_sorted,120,'LineStyle','none')
contourf(independent_var,A_star_sorted,powercoef_mean_sorted_smoothed,120,'LineStyle','none')%,[],'LineStyle','none'
% contourf(f_star_sorted,A_star_sorted,squeeze(powercoef_conv(1,:,:))./powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% quiver(f_star_sorted,A_star_sorted,gradA,gradf)

caxis([-1 0.2])
% colorbarpwn(-6.0,0.2,'colorN',[0 0.5 1],'log',1.5)
colormap(bluewhitered)

contour(independent_var,A_star_sorted,powercoef_mean_sorted_smoothed,[0 0],'LineWidth',4,'LineColor','k','LineStyle','--')
scatter(independent_var,A_star_sorted,60,'.','k')
grid on
xlabel(xlabelstr)
ylabel('{\it A} * = {\it A/D}')
xlim(xlimits)% xlim([0.09 0.27]) %xlim([0.047 0.145]) %
ylim([-0.04 1.12])
% xticks([0.10 0.14 0.18 0.22 0.26])
set(gca, 'Layer', 'top')
grid off
% fstarvalues_extended = (0:0.025:0.3);
heaveaccelcommandlimit = 3.5;
heavevelocommandlimit = 0.5;
a_limit_curve = heaveaccelcommandlimit./(thcknss*(2*pi*(U.*fstarvector/thcknss)).^2);
v_limit_curve = heavevelocommandlimit./(thcknss*(2*pi*(U.*fstarvector/thcknss)));
plot(independent_var(:,1),a_limit_curve,'LineWidth',4,'Color','b','LineStyle','-.')
plot(independent_var(:,1),v_limit_curve,'LineWidth',4,'Color','r','LineStyle','-.')

% % Plot Cp vs. A* and U*
% % contourf(U_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % contourf(f_star_sorted,A_star_sorted,powercoef_mean_sorted,120,'LineStyle','none') %,[],'LineStyle','none'
% % quiver(f_star_sorted,A_star_sorted,gradA,gradf)
% caxis([-1.0 0.2])
% colormap(bluewhitered)
% contour(U_star_sorted,A_star_sorted,powercoef_mean_sorted,[0,0],'LineWidth',3,'LineColor','k','LineStyle','--')
% % scatter(U_star_sorted,A_star_sorted,[],'k')
% grid on
% xlabel('{\it U}* = {\it U/Df }')
% ylabel('{\it A}* = {\it A/D}')
% xlim([3 20])
% ylim([0 1.2])
% Ustarvalues_extended = (0:1:11);
% a_limit_curve = acc_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))).^2);
% v_limit_curve = v_limit./(chord*(2*pi*(flowspeed_fixed./(chord*Ustarvalues_extended))));
% plot(Ustarvalues_extended,a_limit_curve)
% plot(Ustarvalues_extended,v_limit_curve)

c=colorbar();
c.Label.String = '{\it C}_P';
% c.Label.Interpreter = 'Latex';
set(gca,"FontName","Arial"); set(gca,"FontSize",36); set(gca,"LineWidth",2); 

hold off